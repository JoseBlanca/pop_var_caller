# Review — the arithmetic and the values

Branch `ng-prepass-handover`, six commits `629e84ff..346abaaf`, reviewed in the worktree
`/Users/jose/devel/pop_var_caller-review-arith`. Brief: numbers only — test design and prose
belong to the other two reviewers.

Eleven findings, none of them a Blocker or a Major. **No defect in the seam's arithmetic**: every
one of the seven inputs lands where the doc says, the seed's two moments are in the right order,
and no value is dropped or overwritten on any path a caller can reach today. Five wrong numbers
(F1–F5), all in prose about correct code. The one that is about the code is F8 — an axis
`from_prepass` does not check, where the same mistake on the neighbouring axis is caught loudly.

---

## What I re-derived, and what I took on trust

**Measured myself, in the container:**

- the library suite at HEAD: `4928 passed; 0 failed; 11 ignored` (56.7 s);
- `cargo doc --no-deps --lib`: 25 `unresolved link to`, exit 101;
- `examples/ng_prepass_handover_footprint` at three shapes: default (`338 1 78`), `-- 141 1 78`,
  `-- 0 1 78`, plus one shape the report does not quote, `-- 338 1 0`, to separate the genotype
  table from the rest of a stratum's record;
- test-function counts per file at HEAD and at `629e84ff`, and `RunParameters::assemble` call
  sites at `629e84ff`, via `git show`/`git grep`;
- the test *case* counts per module from `cargo test --lib -- --list` (4,939 listed = 4,928 run +
  11 ignored).

**Derived on paper from the diff:** the branch point's 4,920. The diff adds 8 `#[test]` functions
— `generic/estimate.rs` +1 (10→11), `ssr/mod.rs` +1 (90→91), `parameter_estimation/mod.rs` +1,
`calling/run_parameters.rs` +5 (31→36) — and 4,928 − 8 = 4,920. **Consistent.** I could not
measure it: the brief forbids a checkout of `629e84ff`.

**Taken on trust:**

- that `cargo doc` also reported 25 at the branch point (commit `9ebe86bc`'s message says it went
  26 → 25; I verified only that it is 25 now, and that none of the 25 sits in a file this branch
  touched — the nearest is `src/ng/parameter_estimation/ssr/mod.rs:639`, well outside the two
  ranges the diff adds, at 743 and 6101);
- `cargo fmt` / `cargo clippy` exit 0 — not run;
- every mutation-kill claim ("emptying the map fails it", the five deliberate defects, "the
  permutation") — not re-run; that is the test-design brief;
- "one library, which is every sample of both benchmark cohorts here" — not checked;
- the two cited stratum counts, which I confirmed are *present* in their sources but did not
  re-derive: `ng_parameter_prepass_ssr_e5_2026-08-13.md:122` ("a genome's 338 strata per read
  group") and `census_tract_grain_b4_2026-08-14.md:58` ("462,701 kept tracts in 141 strata").

---

## 1. Every number in the report and in the doc comments

### Correct, re-derived

| claim | where | check |
| --- | --- | --- |
| 4,928 passing, 0 failing, 11 ignored | report §"The library suite" | measured |
| 4,920 at the branch point | same | 4,928 − 8 added tests |
| 25 unresolved links, exit 101 | same | measured |
| "8 reads a site at 7 nats against 12 at 9" | `estimate.rs:598` | `UNLIKE_LIBRARIES = [(RG0, 8, −7.0), (RG1, 12, −9.0)]` |
| 2, 5, 10 mismatched bases in a thousand | `ssr/mod.rs:6104` | 1,200 loci × (1, 3, 7) mismatching reads × 1 base ÷ compared |
| over 600,000 / 720,000 / 840,000 bases | same | 1,200 × 50 reads × (10, 12, 14) bases of `AT`-repeat |
| the whole footprint table, all four rows | report §"One sample and a thousand" | reproduced byte for byte |
| 1,141,704 at one sample; 1,141,251 at a thousand | same | measured |
| 473,484 bytes a sample at 141 strata | same | measured (`-- 141 1 78`, 1,000-sample row) |
| 3,391 bytes with no tracts | same | measured (`-- 0 1 78`, 1,000-sample row) |
| 26.6 kB a sample / 26.6 MB at a thousand | same | union 26,595,880 B ÷ 1,000 = 26,595.88 |
| 1.14 GB at a thousand | same | peak 1,141,251,880 B |
| "about 3.35 kB a stratum, at both stratum counts" | same | (1,141,251 − 3,391)/338 = 3,366 B; (473,484 − 3,391)/141 = 3,334 B |
| 78 genotypes = 12 allele lengths at 5 copies | `ng_prepass_handover_footprint.rs:50` | `allele_support` clips at −`reference_repeats`, so −5..=+6 is 12; 12·13/2 = 78 |
| the plan's argument list names six | report §"Where the code went its own way" | plan D1, line 148–150: generic, repeat-tract, joint fit, read groups, batches, ploidy |
| four distinct calibration scales 0.5, 0.25, 0.75, 0.125 | `run_parameters.rs:2050` | rate/mean = .001/.002, .002/.008, .003/.004, .004/.032 |
| "6.06 differences per 10,000 … 0.152 of those" | `run_parameters.rs:667` | segregating 0.0040; 2·0.2·1.0/(1.2·2.2) = 0.15152; product 6.0606e−4 |
| "1.67 in 1,000 … mean of 0.167" | `run_parameters.rs:681` | 0.001 + 0.004·(0.2/1.2) = 1.6667e−3 |

### F1 — Minor — "43 parts in 44" is 42 parts in 43

`doc/devel/reports/implementations/ng_prepass_calling_handover_2026-08-27.md`, §"One sample and a
thousand": *"A driver that projects each sample's substitution rates as that sample finishes, and
releases the rest, never holds the other 43 parts in 44."*

At a thousand samples the driver would hold `union + assembly` = 26,595,880 + 24,000 =
26,619,880 B out of a peak of 1,141,251,880 B. That is **1 part in 42.87**, so the part released
is 41.9 in 42.9 — "42 parts in 43" when rounded. "43 in 44" implies peak/held = 44, which is 2.6%
more than the measurement. As a bare fraction the error is small (43/44 = 97.73% against the
measured 97.67%, and 42/43 = 97.67%), but the "parts in" form states the ratio directly and the
ratio is 42.9, not 44.

### F2 — Minor — the genotype table is 2.50 kB a stratum, not 3.35 kB

Same section: *"What the difference is made of is the allele-length genotype table each stratum's
fit carries: about **3.35 kB a stratum**, at both stratum counts measured here."* Repeated in
`PROJECT_STATUS.md:44`.

The example takes the genotype count as its third argument, so the split is measurable. Running
`-- 338 1 0` alongside the default (1,000-sample row, `inputs` column, per sample):

| shape | inputs a sample |
| --- | ---: |
| 338 strata, 78 genotypes | 1,114,632 |
| 338 strata, 0 genotypes | 270,984 |
| 0 strata | 3,232 |

- genotype table: (1,114,632 − 270,984)/338 = **2,496 B a stratum** — exactly 78 × 32, one
  `GenotypeFrequency` each;
- everything else a stratum carries: (270,984 − 3,232)/338 = **792 B**;
- the difference the sentence is about (on the peak, as the report measures it): 3,366 B.

So the genotype table is **74% of the 3.35 kB**, and 870 B a stratum is the rest of the
`StratumFit` record plus the `BTreeMap` node and `StratumKey`. The 3.35 kB figure itself is right
for the *difference*; the attribution of all of it to the genotype table is 34% too generous. The
paragraph's conclusion ("calling never asks for it") survives either way — calling asks for one
field of the 792 B and none of the 2,496 B.

### F3 — Minor — "99.7%" and "43 parts in 44" are two different fractions on two different bases

Same section: *"1.14 GB at a thousand samples, and **99.7% of it is nothing the seam reads**."*
Repeated in `PROJECT_STATUS.md:41`.

Two numbers are defensible and they are not the same one:

- **the share the repeat tracts add**: 1 − 3,391/1,141,251 = **99.70%**. This is what the next
  sentence measures ("a sample weighs 3,391 bytes instead of 1,141,251"), and it is what 99.7%
  is arithmetically;
- **the share the seam never reads**, i.e. what a driver could release: `inputs`/`peak` =
  1,114,632,000/1,141,251,880 = **97.67%**. This is what the phrase *"nothing the seam reads"*
  says, and it is what F1's "parts in" sentence uses one paragraph later.

The seam does read something out of the 99.7%: one `substitution` field per stratum record. So
99.7% cannot be the share it never reads. The two sentences should be on one base.

### F4 — Minor — "201 other tests in that file" is 90 (or 202 for the module tree)

Report §"A sample's tract substitution rates": *"the permutation left all 201 other tests in that
file green."*

`grep -c '^\s*#\[test\]' src/ng/parameter_estimation/ssr/mod.rs` = **91** at HEAD, 90 at
`629e84ff`, and `cargo test --lib -- --list` confirms `ng::parameter_estimation::ssr::tests::`
holds exactly 91 cases. So *that file* has **90** other tests.

If "that file" was meant as the whole `ssr` module tree, the list gives 203 cases across
`ssr::tests` (91), `locus_offsets` (40), `stratum_table` (37), `slippage` (29) and
`offset_bucket` (6) — **202** others. Neither reading gives 201.

### F5 — Minor — "the constructor had ten callers" — there were 29

Report §"What it is, in one paragraph" and `PROJECT_STATUS.md:27`: *"the constructor had ten
callers and all ten were in its own test module."*

At `629e84ff`, `git grep -n 'RunParameters::assemble(' 629e84ff` returns **29** call sites, all in
`src/ng/calling/run_parameters.rs` at lines 766–1876. The file's `#[cfg(test)]` starts at line
474, so the second half of the claim is right — all 29 are in its own test module — but the count
is off by a factor of about three.

### F6 — Note — "78 … is a floor" holds only from five repeats up

`examples/ng_prepass_handover_footprint.rs:50–52`: *"78 allele-length genotypes a stratum, which
is what a five-copy dinucleotide tract has at two genome copies … Longer tracts have more, so
this is a floor rather than a typical value."*

`allele_support` (`ssr/mod.rs:246`) clips the low end at `−min(reference_repeats, 6)`, so the
support is `min(n,6) + 7` lengths and the diploid genotype count is `A(A+1)/2`: 55 at three
copies, 66 at four, 78 at five, and 91 at six or more. 78 is therefore a floor only over tracts
of five copies or more — and the example's own generator starts at four
(`RepeatCount((index / 6 + 4) as u32)`, line 165), whose true count is 66. The count for a given
repeat number is right; "a floor" is right for the majority of a catalogue and wrong for its
shortest tracts. (The word "dinucleotide" carries nothing here: the support depends on the repeat
count, not the period.)

---

## 2. The footprint example

**The three measured quantities are what the headers claim, and the union window is exactly
right.** `before_union` is taken with nothing between it and `after_inputs`, and
`before_assembled` immediately after the loop, so the `union` reading brackets exactly the four
maps the loop builds — `error_rate`, `minted`, `substitution`, `inbreeding` — and the header names
exactly those four. The temporary `BTreeMap` that `substitution_rate_by_stratum()` returns each
iteration is dropped before `before_assembled`, so it is not double-counted.

**The parts sum to the peak by construction, not by agreement.** `peak = after_assembled −
before_inputs`, and since `before_union == after_inputs` the three differences telescope to
exactly that. The header's "three parts that add up to the peak" is an identity; it is not
evidence that nothing was missed.

**My numbers match the report's exactly** at all three shapes it quotes.

### F7 — Minor — two of the run-wide things `from_prepass` builds appear in no column

`examples/ng_prepass_handover_footprint.rs:296–306`. The example calls `assemble` with
`&BTreeMap::new()` for contamination and `StratumFits::over(&[], BTreeMap::new())` for the
slippage gather. `from_prepass` builds neither of those empty:

- it fills `contamination_by_read_group` from `joint.contamination`, one entry per read group
  (`run_parameters.rs:242–250`). That map is missing from the `union` column;
- with it empty, `assemble` takes the `views.iter().all(Option::is_none)` branch
  (`run_parameters.rs:387`) and stores an **empty** `Vec`, where a real run stores one
  `ContaminationView` per read group. So the `assembly` column measures only the calibration
  vector — visibly so: it is exactly 24 B per read group at all four cohort sizes (24 / 240 /
  2,400 / 24,000), while the header says "the dense per-library **vectors**", plural;
- the slippage gather is stored on the result too (seventh argument, pass-through). At the
  default grain — every read group pooled into one slippage group — it does not grow with the
  cohort and omitting it is harmless. At the **specified** grain, one slippage group per read
  group, which is what the seam's own fixture `a_slippage_gather()` uses, each `StratumRow` holds
  three `Vec`s with one entry per group (`joint/stratum_fits.rs:430–434`), so the gather grows as
  strata × cohort — the same shape as the substitution map.

None of this changes the report's conclusion; the missing parts are small beside 1,141 kB a
sample. But the sentence *"What a run must hold for the whole of calling is the 26.6 kB a sample
of run-wide maps"* is a **lower bound**, not the figure, and the run driver is the reader who
would act on it.

---

## 3. The seam's own arithmetic — `from_prepass`

**No finding on placement or on the seed.** Each of the seven inputs is read where it belongs:

| input | what is taken from it |
| --- | --- |
| `generic_by_sample[i]` | `error_rate`, `minted_errors`, `inbreeding` |
| `repeat_tract_by_sample[i]` | `substitution_rate_by_stratum()` |
| `joint` | `contamination[sample]`, and the two moments |
| `read_groups` | the sample order, and each sample's own read-group set |
| `sequencing_batches`, `ssr_slippage_fits`, `ploidy` | passed straight through |

**The seed's two moments are in the right order.** `run_parameters.rs:263–266` passes
`joint.fitted_alternative_frequency()` then `joint.fitted_diversity()` into
`seed_from_moments(expected_frequency, diversity)`. The two are different newtypes
(`ExpectedAlternativeFrequency`, `ExpectedHeterozygosity`), so a swap at that call site does not
compile — which is worth saying, because the report's defect table lists "the seed's two moments
swapped" as a mutation two tests caught. Whatever that mutation was, it was not a swap of these
two arguments.

**Nothing is silently overwritten.** Every `insert` into the four run-wide maps is preceded by
`its_own_read_group` (`run_parameters.rs:586`), and `ReadGroups` partitions the ids across
samples, so no two iterations can write the same key. `checked_read_group_count_of`
(`run_parameters.rs:611`) then refuses a gap in the ids and refuses a rate without its minted
total or the reverse, in both directions — which is exactly what the new doc on
`GenericSampleParameters::minted_errors` (`generic/mod.rs:465`) claims.

### F8 — Minor — the run's ploidy is never checked against the substitution keys' ploidy

`src/ng/calling/run_parameters.rs:252–256` inserts every `StratumKey` the pre-pass produced,
whatever `key.ploidy` it carries, and passes the run's `ploidy` alongside as a separate argument.
The lookup, `FrozenParameters::ssr_substitution_rate_at` (`src/ng/calling/mod.rs:1030–1039`),
rebuilds the key with `ploidy: self.ploidy`.

So if the repeat-tract pre-pass ran at a ploidy other than the one handed to `from_prepass`,
**every** tract substitution rate becomes unreachable — `ssr_substitution_rate_at` answers `None`
at every locus, the caller falls back, and the run finishes. Nothing says so. Compare the
read-group axis, where the identical mistake is caught loudly by `its_own_read_group` with a
sentence naming the sample. `from_prepass`'s `# Panics` list (`run_parameters.rs:170–180`) names
the two list lengths, the missing contamination row and the doubly-claimed read group; it does not
name this.

In practice both ploidies come from one run config (`SsrAccumulators::new(ploidy)` on one side,
the driver's argument on the other), so this is a latent hazard rather than a live bug — and the
run driver, which is out of scope here, is what would introduce the second source. Worth one
assertion at the seam, where the two meet.

### F9 — Note — a missing top read group is caught, but by the wrong assertion

`from_prepass` holds `read_groups` and never compares its read-group count to the count `assemble`
derives from the union of the per-sample maps. A sample that arrives missing its highest-id
library shrinks the run's read-group axis silently at that step; it is then caught — by
`assemble`'s batching assert (`run_parameters.rs:394`), whose message blames the declared
batching. That is only true as long as the batching was minted over the run's real read-group
table.

---

## 4. The two new accessors

**`SsrSampleParameters::substitution_rate_by_stratum` (`ssr/mod.rs:772`)** returns exactly what
its doc says: one clone of `fit.substitution` per entry of `by_stratum`, key for key. Empty map in,
empty map out — exercised by the one-sample test, which builds
`SsrSampleParameters::of_substitution_rates(&[])` and assembles. The heading says
"`(library, tract shape, ploidy)`" and `StratumKey` is those three fields, so the shape claim is
right too.

**`JointFit::fitted_alternative_frequency` (`joint/fit.rs:405`)** returns
`ExpectedAlternativeFrequency::try_new(density.expected_alternative_frequency()).ok()`.
`try_new` goes through `checked_probability`, a `(0.0..=1.0).contains` test
(`ng/types.rs:766–773`), so it rejects NaN as well as out-of-range — which matters, because
`expected_alternative_frequency` is `p_fixed_alt + p_segregating()·a/(a+b)` and a degenerate
`a = b = 0` gives NaN. The doc's "`None` means the fit produced something that is not a frequency"
covers that case correctly, and its unreachability argument is sound: `p_segregating()` is
`(1 − p_invariant − p_fixed_alt).max(0.0)`, so `p_fixed_alt + p_segregating ≤ 1` and the Beta mean
is in `[0, 1]`.

### F10 — Note — "a wrap, not a computation, exactly as `fitted_diversity` beside it" — the two read different places

`joint/fit.rs:392–393`. `fitted_diversity` wraps the **stored field** `JointFit::expected_heterozygosity`;
`fitted_alternative_frequency` **recomputes** from `self.density.value`. In the fitter's own output
they do come off one density — `fit.rs:1571–1579` sets `density: Estimate { value: parameters.density, … }`
and `expected_heterozygosity: parameters.density.expected_heterozygosity()` in the same struct
literal — so `seed_from_moments`'s "**Both moments come off the same fitted density**"
(`run_parameters.rs:124`) is true of every `JointFit` this crate produces.

But every field of `JointFit` is `pub`, and a test in that same file constructs a fit where they
are independent: `fn a_fit_carrying(density: FrequencyDensity, expected_heterozygosity: f64)`
(`fit.rs:3173`). So the pairing is a convention held by one constructor, not an invariant of the
type, and the accessor doc's "exactly as `Self::fitted_diversity` beside it" describes a symmetry
the code does not have. Nothing is wrong today; the sentence would mislead the next person who
builds a `JointFit` by hand.

---

## 5. The claim that drove the design decision — `StratumFit::substitution`

**Verified. The claim holds.** There is no path in the tree where `StratumFit::substitution`
differs from the value `substitution_rates` produced for that key.

- `assemble_sample_parameters` (`ssr/mod.rs:2235`) is the only route that builds an
  `SsrSampleParameters` from evidence. It looks the rate up by the same `StratumKey` it is
  building the record for (line 2246), and hands it to `stratum_fit` (line 2270), which stores it
  in the `substitution` field **by move, unmodified** (line 2331). Nothing between reads it,
  scales it, or replaces it.
- all four in-tree call sites pass `&substitution_rates(&accumulators)` over the same accumulator
  they hand in — `ssr/mod.rs:6026`, `:6248`, `:6358`, and the new test at `:6157`.
- there is no production caller of `assemble_sample_parameters` at all yet; the run driver is out
  of scope on this branch.
- the only other constructions are `SsrSampleParameters::of_substitution_rates` (new, test-only,
  `ssr/mod.rs:789`), which stores the rates it is handed, and the footprint example
  (`ng_prepass_handover_footprint.rs:208`), which synthesises a constant.

So a new field would indeed have been a second copy. **No finding against the decision.**

One structural caveat, not a defect: `assemble_sample_parameters` takes `accumulators` and
`substitutions` as *independent* arguments and never checks that the second was computed from the
first. A caller that passed a map from another accumulator would get records carrying those
values. Every caller today does the right thing.

### F11 — Note — the new doc's "refuses to build a record" is a whole-sample abort, and one sentence it rests on is contradicted by two tests

The new doc on `substitution_rate_by_stratum` (`ssr/mod.rs:765–770`) says:

> *"[`substitution_rate_of`] answers `None`, [`substitution_rates`] leaves the key out, and
> [`assemble_sample_parameters`] refuses to build a record for a stratum with no rate. So a key
> that reaches this map has been measured, and there is no absence left for it to represent."*

The conclusion is right, and the "refuses" wording is the house term — the pre-existing test that
pins the behaviour is called `a_stratum_with_no_compared_bases_is_refused_rather_than_given_a_zero`
(`ssr/mod.rs:6708`) and argues the choice deliberately. Two things are nevertheless worth having
straight, because the run driver is the reader:

1. **The refusal is a `panic!`, not a skip** (`ssr/mod.rs:2246–2257`): *"no substitution rate was
   fitted for {key}, which holds {n} loci over {m} compared bases"*. The whole sample's parameter
   assembly dies; the stratum is not dropped and the rest kept. "Refuses to build a record" reads
   as the second.

2. **"Nothing observed can reach it" is false, and is asserted twice.** `substitution_rate_of`'s
   own doc (`ssr/mod.rs:1161–1162`, pre-existing) says *"a stratum with loci and no compared bases
   … Nothing observed can reach it: a read reaches a table only through a complete witness, which
   compares its bases."* Two tests in the same file construct exactly that state, from a tract
   every read shows as entirely deleted — `a_stratum_with_loci_and_no_compared_bases_is_absent_rather_than_zero`
   (`:3686`, asserting `loci() == 1` and `bases_compared() == 0`) and the `should_panic` test at
   `:6708`, **whose own doc repeats the sentence in the act of building the case it denies**:
   *"a tract every read shows as entirely deleted files a shape and compares no bases. Nothing
   observed reaches it."* `SsrAccumulators::strata()` (`:1113`) iterates the same `by_stratum` map
   `table_for` reads, so the key does reach `assemble_sample_parameters`.

So the state is reachable from data and the panic is the designed answer to it. All pre-existing;
this branch's contribution is a doc comment that leans on the unreachability claim. The number to
carry forward is that a single fully-deleted repeat tract aborts a sample's pre-pass assembly.

---

## Nothing found

- Section 3's placement and ordering questions: nothing wrong (F8 and F9 are about checks that are
  absent, not about values that are wrong).
- Section 5's central question — can the two ever differ: no, on every path in the tree.
