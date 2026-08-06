# Sketch 3 — what the calling maths is handed

**Date:** 2026-08-06. **Branch/worktree:** throwaway worktree off `main` at `1e5ffa8`.
**Plan:** `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md` §4, Sketch 3.
**Verdict: the input shape does not matter here.**

---

## The answer, in one sentence

**The input shape does not matter at the calling boundary, because the arithmetic is
167 times the data movement** — on a 50-sample tomato cohort the expectation-maximisation
step costs **999,224 instructions per called locus** while materialising that locus as an
owned record costs **5,999**, so handing the step a record instead of slices makes the
whole calling step **0.59 % more expensive**, and once the step also has to produce the
`PosteriorRecord` the VCF writer consumes, **0.10 % more expensive**. All four figures are
measured, floor-subtracted, minimum of three alternated sweeps.

Two things sharpen that. First, **nothing in the calling step resisted a columnar input** —
the tree already contains a borrowed-slice EM entry point and it had never been driven
columnar; doing so took 21 lines. What resisted is on the *output* side, and it is the
reason 0.59 % collapses to 0.10 %. Second, **167 is not a constant**: it is 116 on
biallelic-only loci and 35 on a 10-sample cohort, because the EM grows with cohort size
faster than record materialisation does. Even at its smallest measured value the arithmetic
is 35 times the movement and the record shape costs 2.6 % of the calling step.

---

## 1. What was built, and what "the same calculation twice" means

Both arms drive **the real engine**, not a proxy:
`var_calling::posterior_engine::run_em_columnar` (`src/var_calling/posterior_engine.rs:2228`)
— the E-step / M-step loop, `summarise_posteriors`, and the exact-allele-count QUAL
convolution. Neither arm reimplements any of it.

The tree turned out to already have both shapes of entry point:

| entry point | shape | file:line |
|---|---|---|
| `run_em_for_record(record: MergedRecord, …)` | owned record in | `posterior_engine.rs:2153` |
| `run_em_columnar(inputs: EmInputs<'_>, …)` | borrowed slices in | `posterior_engine.rs:2228` |

`run_em_for_record` is a 40-line shim: it destructures the record, wraps the allele vector in
a `MergedAllelesView`, and calls `run_em_columnar`. So the two arms differ **only** in how
the inputs reach `EmInputs`, and the body they reach is byte-for-byte the same code.

- **Arm A, "record".** Per locus, materialise an owned `MergedRecord` out of a block —
  one `Vec<u8>` per allele, one `Vec<MergedAllele>`, and four more `Vec`s — hand it to the
  EM by value, and drop it. This is `run_em_for_record`'s body minus the output join.
- **Arm B, "column".** Per locus, form `EmInputs` from slices into a block the caller owns
  and reuses, plus a 21-line `ColumnAllelesView`. Nothing is materialised, nothing is copied.

Both end at `EmOutputs`, so the boundary under test is isolated: everything downstream is
identical and cancels.

Two further arms measure the **realistic** stage, in which step 6 (the VCF writer) requires a
`PosteriorRecord`:

- **`record-full`** — arm A, then join the EM outputs with the record's own (moved, free)
  passthrough fields and run `prune_unsupported_alleles`.
- **`column-full`** — arm B, then build the same `PosteriorRecord` by **copying** the
  passthrough columns out of the block, and prune.

Two ablation arms materialise the input and skip the EM call, so
`record` − `record-nomath` is arithmetic and `record-nomath` − `floor` is data movement.

Sketch code: `sketch3.diff` beside this file — 953 new lines in
`examples/sketch3_calling_input.rs`, plus 8 one-word visibility changes in
`posterior_engine.rs` (`pub(crate)` → `pub`) so an example can reach the columnar entry
point. No behavioural change to the library.

## 2. Fixture — measured, and which `.psp` files

| | |
|---|---|
| source | `benchmarks/tomato1/results/ours/cohort/psp/`, first 50 files alphabetically |
| freshness | **regenerated**, `created = 2026-07-07T16:35:18Z`, carries the `kind = "snp"` tag |
| region | `SL4.0ch01:3,400,000-16,000,000`, which covers 3 of the benchmark's 100-kb BED regions |
| path to merged evidence | `PspReader::region_records` → `PerPositionMerger` → `DustFilter` (defaults) → `VariantGrouper` (defaults) → `PerGroupMerger` (defaults) → `passes_min_alt_obs(min_obs = 2)` |

**On staleness.** The alternative fixture, `/Users/jose/devel/pop_var_caller/tmp/aligned_psp/`,
is dated `created = 2026-06-02T06:38:16Z` and predates the `kind` tag. A previous session
measured that regenerating stale tomato `.psp` cut peak memory 17–20 % and wall time 26 %
with 0.06 % of calls differing (**cited**, not re-measured here). **This sketch used the
regenerated 2026-07-07 set only.**

The skew the plan §5 insists on is present and measured:

```
14,683,868 psp records  ->  224,067 merged groups
                        ->   23,552 loci reach the EM   (200,515 dropped by min_alt_obs)
                            of which 14,913 multiallelic
work shape: mean 2.85 alleles, 5.91 genotypes, 3.93 EM iterations, 23,551/23,552 converged
```

**10.5 positions in 100 that produce a merged group survive to the EM.** Only those ever
become objects in either arm — which is exactly the property that makes an all-loci synthetic
fixture misleading.

## 3. Correctness — the two arms are bit-identical

Before any timing was quoted, `--mode verify` compared the arms locus by locus on all
23,552 loci at 50 samples, and again on the 8,639-locus biallelic subset and the
6,709-locus 10-sample fixture:

```
verify:       loci=true genotypes=true gq_bits=true qual_bits=true af_bits=true -> BIT-IDENTICAL
per-locus:    0 of 23552 diverged
verify-full:  ... -> BIT-IDENTICAL
```

Compared as **raw IEEE-754 bit patterns**, not with a tolerance: `best_genotype`,
`n_genotypes`, every `gq_phred`, `qual_phred`, every `allele_frequencies` entry, and every
`posteriors` entry. **No tolerance was needed and none was allowed.** That is not luck — the
two arms reach the same accumulation loops over the same `f64` values in the same order,
because they share `run_em_columnar` verbatim.

## 4. The numbers

`instructions retired` from `/usr/bin/time -l`, floor-subtracted with a `--mode floor` run
that builds the identical fixture and executes no EM. Modes alternated round-robin within
one script (`sketch3_measure.sh`), three sweeps, **minimum** reported. Raw TSVs beside this
file. Wall clock is not quoted anywhere: the host has 6 high-performance and 12 low-energy
cores and was under heavy load from other work throughout.

### 4.1 Headline — 50 samples, all 23,552 loci

Floor run-to-run spread was 119,535,225 instructions = **254 per locus**, so any per-locus
figure below ~250 is at the noise floor.

| mode | instructions per called locus | what it is |
|---|---:|---|
| `column-nomath` | **76** | forming slices into the block — at/below the noise floor |
| `record-nomath` | **5,999** | materialising one owned `MergedRecord`, and dropping it |
| `column` | **999,300** | slices in → EM → `EmOutputs` |
| `record` | **1,005,233** | owned record in → EM → `EmOutputs` |
| `column-full` | **1,014,806** | slices in → EM → `PosteriorRecord` + prune |
| `record-full` | **1,015,817** | owned record in → EM → `PosteriorRecord` + prune |

Derived:

- **arithmetic** = `column` − `column-nomath` = **999,224 instructions per locus**
- **data movement, record shape** = **5,999**; **column shape** = **76**
- **N = arithmetic / record-shape movement = 167**
- **record penalty, EM only: +5,933 instructions = +0.59 %**
- **record penalty, realistic stage: +1,011 instructions = +0.10 %**

### 4.2 The plan's §7 table

Per called locus. "Per read at a position" does not apply — the calling step never sees a
read; its unit is a locus × cohort. Peak RSS and allocations are discussed below the table.

| arm | instructions per called locus | per read at a position | peak RSS | allocations per locus | bytes copied at the handoff | lines of source | what it cannot do |
|---|---:|---:|---:|---:|---:|---:|---|
| **A — record in** | 1,005,233 | n/a | 2,611.5 MB (fixture-bound; see §4.4) | 21.71 | 10,379 | 85 | nothing that was asked of it |
| **B — slices in** | 999,300 | n/a | 2,537.7 MB (fixture-bound) | 13.85 | 0 | 162 | cannot survive `prune_unsupported_alleles`; cannot feed the VCF writer without materialising a record anyway (§5) |
| **A-full — record in, record out** | 1,015,817 | n/a | 2,611.9 MB | 33.44 | 10,379 in, 0 out | +59 | — |
| **B-full — slices in, record out** | 1,014,806 | n/a | 2,612.2 MB | 32.44 | 0 in, ~8,014 out | +53 | — |

Allocations and bytes are **measured** with dhat (`--features dhat-heap`,
`--target-dir target-dhat`), one sweep, floor-subtracted; the JSON files are beside this
report. The two independent counts agree exactly: dhat's `record` − `column` block delta is
**184,969**, which is the sketch's own count of `n_alleles + 5` allocations over the same
23,552 records, to the block. "Bytes copied at the handoff" is the dhat byte delta between an
arm and the arm below it: 14,021.7 − 3,642.4 = **10,379 bytes per locus** for the record
shape's input; and for `column-full`, 9,797.8 (its cost over `column`) minus the 1,783.8 that
`record-full` also pays for the join and prune = **8,014 bytes per locus** of passthrough
copied out of the block.

Source lines are non-blank, non-comment, counted over the sketch:
column machinery (block struct, offsets, `ColumnAllelesView`, `slice_inputs`, the loop)
**162 lines** against the record machinery's **85**. The extra 77 lines are entirely offsets
bookkeeping and the view adapter. The `-full` rows add 59 lines to the record arm and 53 to
the column arm, on top of a `join_posterior_record` helper (32 lines) they share.

### 4.3 N is not a constant — it grows with cohort size

| fixture | loci | arithmetic | record-shape movement | **N** | record penalty, EM only | record penalty, full stage |
|---|---:|---:|---:|---:|---:|---:|
| 50 samples, all loci | 23,552 | 999,224 | 5,999 | **167** | +0.59 % | +0.10 % |
| 50 samples, biallelic only | 8,639 | 578,930 | 5,011 | **116** | +0.85 % | +0.11 % |
| 10 samples, all loci | 6,709 | 154,740 | 4,422 | **35** | +2.58 % | +0.39 % |

All three rows measured the same way. The trend has a plain cause: going from 10 to 50
samples multiplies the EM's arithmetic by 6.4 (154,740 → 999,224) but record materialisation
by only 1.36 (4,422 → 5,999), because most of a record's allocations are per-allele and
per-record fixed cost, not per-sample bytes. **So the shape matters least exactly where the
cohort is largest** — the case this pipeline is built for. The biallelic row is the
cheapest EM available (3 genotypes, mean 3.81 iterations) and still puts the arithmetic at
116 times the movement.

### 4.4 Peak resident memory — the honest reading

Process peak RSS cannot separate the arms in this sketch, and dhat says why: `At t-gmax` is
**2,403,196,548 bytes in the floor run** and 2,403,196,549 / 2,403,196,554 in the arms — the
same live heap to nine significant figures. The high-water mark belongs to the fixture build
(the per-sample `.psp` records buffered before merging), not to either arm.

What can be said, measured:

- **The block is 242.3 MB for 23,552 loci; one live record is 10.4 kB.** The column shape
  requires the whole block resident; the record shape requires one locus. That is a factor of
  23,000 in *live* working set, and it favours records.
- **`At t-end` is 43,883 bytes in 77 blocks for both arms** — neither leaks or accumulates.
- Peak RSS did differ slightly between modes (`record` 2,611.5 MB against `floor`'s
  2,498.8 MB, +112.8 MB) even though live heap did not. With 471,040 alloc/free cycles of
  ~10 kB in the record arm and dhat reporting flat live heap, the most likely cause is
  allocator page retention rather than live data — **that is interpretation, not
  measurement**. It did not reproduce at 10 samples, where the record arm's peak RSS
  (494.6 MB) is *below* the floor's (494.5 MB, i.e. within noise).

### 4.5 The profile — corroboration only, on a contended host

**Method warning.** `sample` records a thread's stack whether or not it is on a core, so
under the load this host was carrying, a stall while descheduled lands as time in whatever
frame was on top. The split in §4.1 comes from the **ablation**, which is instruction counts
and is near-immune to scheduling. The profile below is **indicative only** and was taken
against the `[profile.profiling]` build (`lto = false, codegen-units = 16`), so its shares
locate work and do **not** transfer to release.

Call-tree shares under `run_record_arm`, 50,730 samples total:

| frame | samples | share |
|---|---:|---:|
| `run_em_row` → `run_em_columnar` and below (the arithmetic) | 49,878 | 98.3 % |
| `materialise` (building the record) | 558 | 1.10 % |
| dropping the record (self time at the close of `run_em_row`) | ~187 | 0.37 % |

Inside the arithmetic, `lgamma_r` alone is 11,025 of the 50,667 top-of-stack samples
(21.8 %), reached through
`fill_log_indep_per_g_from` — the Dirichlet-multinomial term. That is what 999,224
instructions per locus is made of.

The profile puts the input shape at **1.5 %** of the arm against the ablation's **0.59 %**.
The gap is expected and points the same way: the profiling build has LTO off, and time
samples over-weight `malloc`/`memmove`, which miss cache far more than they retire
instructions. **Quote the ablation.**

## 5. What resisted a columnar input — the useful finding

**The EM itself did not resist at all, and this is worth being precise about.** The EM reads,
per locus: four scalars (`ploidy`, `n_samples`, `n_genotypes`, and the locus coordinates);
one ragged `f64` array (`log_likelihoods`); a second ragged array (`scalars`) **only when
contamination is configured**; and exactly four queries on the allele set — `len()`,
`ref_len()`, `seq_len(i)`, `is_compound(i)`. It never reads allele bytes and never reads
compound constituents. Writing `ColumnAllelesView` over two `u32`/`bool` columns was 21 lines
and worked first time.

That boundary already existed. `EmInputs<'_>` and `AllelesView`
(`posterior_engine.rs:1145-1202`) were written for a column-native consumer — the doc comment
names a `&UnifiedAllelesColumns` impl "used by the column-native worker once Phase A.2 step 2
lands". **`UnifiedAllelesColumns` appears only in comments; it was never built, and
`run_em_for_record` is `run_em_columnar`'s only caller in the tree.** This sketch is the
first columnar consumer of a boundary designed three phases ago for one.

Three things did resist, none of them the EM's arithmetic:

1. **`prune_unsupported_alleles` is record-shaped by construction**
   (`posterior_engine.rs:872`). It drops ALT alleles no sample's argmax genotype uses, then
   re-enumerates genotypes and rebuilds `alleles`, `allele_frequencies`,
   `compound_frequencies`, `scalars`, `chain_anchor_flags`, `posteriors`, `best_genotype` and
   `gq_phred` at the new width. In record shape that is a field swap on an owned object. In a
   block it changes one locus's ragged widths, so every offsets array from that locus onward
   shifts — you would need a rewrite pass, or a per-locus side table, which is precisely the
   "second mechanism" the plan §3 records production paying for. **I did not write it
   columnar.** The `column-full` arm materialises a record and prunes that.

2. **Step 6 wants a record, so the columnar arm has to build one anyway** — and this is what
   turns a 0.59 % win into a 0.10 % one. `PosteriorRecord` implements `vcf::VcfWritable`
   (`posterior_engine.rs:1019`). The record arm *moves* `alleles` / `scalars` /
   `other_scalars` / `chain_anchor_flags` from its input into its output for free; the
   columnar arm must **copy** them out of the block. Measured: `record-full` allocates 33.44
   blocks per locus, `column-full` 32.44 — **a difference of one allocation per locus**. The
   layout moved the copy; it did not remove it. That is the same sentence production's
   redesign earned (plan §3), reproduced at a different boundary.

3. **The lifetime problem that killed the merger's borrowed view does not arise here, and
   that is a fact about this stage, not about borrowing.** `EmInputs<'a>` is consumed inside
   one call and never held across a `next()`, so nothing is self-referential. Sketch 2's
   cursor problem is real; it just is not this stage's problem. A columnar EM input is
   *safe*; it is simply not *worth* anything.

**One thing worth fixing regardless of shape.** The EM allocates 13.85 blocks per locus in
**both** arms — 1.8 times what materialising the input record costs. Two sites account for
most of it, both measured by dhat and confirmed in the profile:
`validate_record_shape` calls `genotype_order(ploidy, n_alleles)`
(`posterior_engine.rs:2549` → `per_group_merger.rs:522`), which builds a `Vec<Vec<u8>>` and
sorts it purely to compare its length against `record.n_genotypes`, then discards it — that
is `n_genotypes + 2` allocations and an *O*(n log n) sort per record, and `shape_for`
(`shape.rs:76`) already caches the identical table. The other is `EmOutputs`' five
clone-out `Vec`s (`posterior_engine.rs:2511`), which is deliberate and documented. **If
anyone wants allocations out of the calling step, that is the lever — not the input shape.**

## 6. What building each arm was actually like

**Arm A took about twenty minutes.** Destructure the record, borrow its fields, call the EM.
There was no decision to make and nothing to check afterwards.

**Arm B took about two hours**, and none of that was the EM. It was the block: nine
`off[i]..off[i+1]` derivations at three call sites, deciding whether the offsets array
carries a leading zero or a trailing sentinel (it needs both — a leading zero for the ragged
per-locus arrays and a trailing sentinel for `allele_bytes_off`, and getting that wrong is a
silent off-by-one, not a compile error), and one genuinely awkward piece of Rust: `EmInputs`
holds `&'a dyn AllelesView`, so the view must outlive the borrow and cannot be constructed
inside the struct literal. That forced `slice_inputs` to return a tuple of
`(ColumnAllelesView, EmInputsParts)` where the record arm needs no such helper at all — the
21-line adapter plus a 43-line splitter to work around a lifetime, to save 0.59 %.

Both arms then agreed on all 23,552 loci, bit for bit, first try. **That is the strongest
evidence in this report that the boundary is genuinely shape-neutral**: the EM does not care
where its `f64`s come from, and its authors had already made sure of that.

**The sketch cost one day and removed a question.** Per the plan §4, that is the outcome it
was funded for.

## 7. Recommendation

**Keep records at the calling boundary.** It is the owner's stated default, it is 85 lines
against 162 for the same stage, it holds one locus resident instead of a whole block, and on
a 50-sample cohort it costs 0.10 % of the calling step once the stage is asked to produce
what the VCF writer actually consumes. The measured penalty only rises above 1 % when the
cohort falls to about 10 samples, and even there it is 2.6 % of a step whose arithmetic is
35 times its data movement.

The one caveat to carry into the sketch 1 and 2 decisions: **this result is about the calling
step and does not generalise upward.** The generator's 4,533 instructions per covered base
(**cited**, plan §3) has no arithmetic hiding it. Records are free here precisely because
`lgamma` costs a million instructions a locus; that argument is unavailable one stage up.

## 8. Files

All in this directory, all prefixed `sketch3_` because the three sketches share it.

| file | what |
|---|---|
| `sketch3_calling_input.md` | this report |
| `sketch3.diff` | the sketch: `examples/sketch3_calling_input.rs` (953 lines) + 8 visibility changes in `posterior_engine.rs` |
| `sketch3_measure.sh`, `sketch3_measure_n.sh` | the alternating instruction-count harness |
| `sketch3_reduce.sh` | reduces a sweep TSV to floor-subtracted per-locus figures |
| `sketch3_dhat.sh` | allocation counts per arm |
| `sketch3_profile.sh` | the corroborating profile (contended host, indicative only) |
| `sketch3_main_contam_off.tsv` | 50 samples, all loci — the headline sweep |
| `sketch3_biallelic_contam_off.tsv` | 50 samples, biallelic-only |
| `sketch3_n10_contam_off.tsv` | 10 samples, all loci |
| `sketch3_pilot_contam_off.tsv` | first sweep, before the `-full` arms existed; kept for the record |
| `sketch3_dhat-*.json` | per-arm heap profiles |
| `sketch3_sample-record.txt` | the sampling profile of the record arm |

**Not measured, and why.** The contamination-on configuration — the only one in which the EM
reads the `scalars` column as well as `log_likelihoods` — was left unmeasured. It adds a
mixture pre-pass, so it moves the ratio further in the direction the report already argues
(more arithmetic against the same movement), and host contention made every extra sweep
expensive. The sketch supports it behind `--contamination` if anyone wants the number.

**How to reproduce.**

Apply `sketch3.diff` to a checkout at `1e5ffa8`, then, with the paths inside the scripts
pointed at that checkout:

```
cargo build --release --example sketch3_calling_input
bash sketch3_measure.sh out.tsv --repeats 20     # ~25 min, 21 runs
bash sketch3_reduce.sh out.tsv 23552 20
```

Correctness on its own:

```
./target/release/examples/sketch3_calling_input \
    --psp-dir benchmarks/tomato1/results/ours/cohort/psp --n-samples 50 \
    --reference $HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
    --chrom SL4.0ch01 --start 3400000 --end 16000000 --mode verify
```
