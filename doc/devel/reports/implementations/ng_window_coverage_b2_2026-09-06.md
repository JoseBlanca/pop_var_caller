# window coverage — B2: does a one-base record's head say what its evidence says?

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone B, step B2
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.1, §6 trap 7
**Review:** [ng_window_coverage_b2_2026-09-06.md](../reviews/ng_window_coverage_b2_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [B1](ng_window_coverage_b1_2026-09-06.md)

## The answer

**Yes, at every one of 8,784,182 records checked on two tomato stores, and the mechanism that
could have broken it does not occur there at all.** The rule B1 built takes a sample's depth at a
one-base record from the count the record's head already carries, rather than decoding the
evidence. That is only correct if the two numbers agree, and the argument for their agreeing —
that at a one-base locus no read's evidence can stop inside the locus, so every observation covers
the whole of it — is an argument about how the walk emits records, not a theorem. Spec §6's
seventh trap says so and asks for a measurement.

**Zero disagreements, and zero one-base records carrying a read whose evidence stops short.** The
second number says why the first is what it is — a read whose evidence stops inside the locus is
the only way the two counts can differ — but the two are **not independent findings**: at a
one-base locus such a read necessarily covers the only base, so "no disagreement" and "no such
read" are one fact seen twice. What makes it evidence rather than a tautology is that the shape is
constructible and the probe catches it: a test builds a one-base record whose head says 5 and
whose evidence says 7, and the measurement reports it.

Spec §3.1 stands as written; the rule does not become "build every body".

**And the cheap branch decides almost every position — 992 in every 1,000 on the six-accession
slice, 994 on the wider store.** That is the positions figure and not the records one, and they
are different questions: a wide **generic** record reports depth at one position, a wide **tract**
at every position of its span. In records, 1 in 871 spans more than one base on the slice and 1 in
1,221 on the wider store.

## What was measured, and on what

| store | ground | mean head depth | records | one-base records | disagreeing | with a read stopping short |
|---|---|---|---|---|---|---|
| six tomato accessions (`SRS3394606`, `SRS3394711`, `SRS3394712`, `SRS3394712_SRR7279484`, `SRS3394713`, `SRS3394714`) | the first two 100 kb intervals of `benchmarks/tomato1/regions.bed`, on `SL4.0ch01` | 14.4 | 1,164,208 | 1,162,872 | **0** | **0** |
| `SRS3394712` — the sample of `SRR7279481.p1` | all 80 intervals of that BED: 8 Mb over all 12 chromosomes | 10.3 | 7,627,555 | 7,621,310 | **0** | **0** |
| both | | | 8,791,763 | **8,784,182** | **0** | **0** |

*Mean head depth* is the mean of `reads_compared_with_reference` over the one-base records — the
reads whose whole sequence over the base was compared against the reference. It is below read
depth, and it is what the rule under test actually reports. **These are the benchmark's slice
CRAMs, which are deeper than the 63-accession tomato cohort's three reads a position**, so the
figures here describe a 10–14× corner of the depth range and not the 3× one.

What the two branches decide, in records and in positions:

| store | wide generic | wide tract | bases in wide records | of them, on tracts | positions reported | decided by the build branch | decided by the head |
|---|---|---|---|---|---|---|---|
| six accessions | 430 | 906 | 10,415 | 8,940 | 1,172,242 | 9,370 | **0.9920** |
| `SRS3394712`, 80 intervals | 2,081 | 4,164 | 50,647 | 43,897 | 7,667,288 | 45,978 | **0.9940** |

A wide generic record contributes one position and a wide tract its whole span, which is why the
bases are split by kind: pooling them makes the positions figure underivable from the output. No
repeat-cluster (`SsrBundle`) record appeared in either store, as
[`locus_generation.md`](../../ng/spec/locus_generation.md) §11 says. **No one-base record was a
tract in either store**, which is spec §3.1's "a single-base record is generic by construction",
measured rather than assumed.

## The deviation from the plan, and what was substituted

**The plan asks for the six-accession slice and "the whole-genome `SRR7279481` store". No ng store
this build can read, and none over a whole genome, exists on this machine.** Under
`/Users/jose/devel` there are ng-format `.psp` files above 1 MB — 330 of them — but the only large
one outside this branch's own scratch is a 51 MB store in another worktree written to an older
record layout, which this build refuses (`UnsupportedRecordEncoding … expected
"non-reference-reads", found "locus-kind"`). Everything large under
`/Users/jose/devel/pop_var_caller` — 1,568 files — is production's format, which ng's reader
cannot open at all. And the tomato CRAMs in `benchmarks/tomato1/` are cut to the benchmark's 80
intervals, so no whole-genome walk can be produced from them here.

**Substituted: the same accession over all 80 intervals** — 8 Mb across all 12 chromosomes,
against the slice's 200 kb on one — generated for this step with `generate-psps`. It is 40 times
the ground of the slice and 6.5 times its records, and it reaches every chromosome, so it answers
the "at scale, on varied ground" half of what the plan wanted. What it does not answer is whether
some corner of the genome the benchmark does not sample behaves differently. **The claim is
measured over 8 Mb of tomato at 10–14 reads a position, not over a genome and not at 3×.**

## Changes made

- **[`examples/ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)** —
  new. `walk` decodes every record of a store once and shows it to whatever `RecordMeasurement`s
  the caller passed; `SingleBaseEquality` is the one measurement this step needs, and at every
  record covering one base it compares the head's `reads_compared_with_reference` against the
  depth `num_obs_along_locus()` gives at that record's only position — the number the rule would
  have used had it built the body. It prints per-store and total counts, up to twenty
  counter-examples in full, and then asserts both that something was checked and that nothing
  disagreed. Plan step C3's whole-store window recomputation is a second implementor of that
  trait and one more entry in the slice; nothing in the walk or in this measurement moves.
- **[`src/ng/window_coverage/depth.rs`](../../../../src/ng/window_coverage/depth.rs)** — the
  rule's own documentation now carries the measurement, with its subjects and their depth, in
  place of an unquantified "great majority".

## What the review changed, and it was the difference between a reading and a check

**The probe as first written could not tell "nothing disagreed" from "nothing was checked".** Its
only assertion was that the disagreement count was zero, which holds vacuously over an empty set:
a reviewer wrote a psp with no records, ran the probe on it, and got a block of zeros and **exit
0** — the same verdict as the 7.6-million-record run. It now asserts, per store, that a record
covering one base was seen, naming the file if none was.

**And the discrimination check this report first offered did not exercise the mechanism the claim
is about.** Pointing the comparison at the wrong head field fires on records that carry no partial
witness at all, so it shows the operand is live and nothing more. Two mutations prove the gap:
disabling the partial-witness detector, and comparing the head against a recomputation of its own
quantity, **leave every printed number identical on both tomato stores**. Both are killed by the
fixture now shipped as a test — a one-base record with five whole reads and two whose witness is a
reach, where the head says 5 and the evidence says 7. The test asserts both quantities off the
record itself, so it pins that the shape is a genuine counter-example and not merely one the probe
labels as such.

Smaller, and each applied: a `first().unwrap_or(0)` that could only ever turn a missing answer
into an agreement is now a slice pattern that panics; every disagreement was retained to print
twenty (7,605,320 of them, measured, on a run where the rule failed wholesale) and is now capped;
the partial-witness test is an exhaustive match rather than `!= Complete`; the totals line prints
whatever the store count, beside a `stores-walked` line, so a run that died part-way through is
not mistakable for a complete report; and the counter-example line carries the two values the
comparison actually used rather than re-reading one of them.

**Two numbers in this report's first draft were wrong**, both in the direction that flatters the
rule: "about 1 record in 900" is 1 in 871, and "999 positions in every 1,000" was the *record*
share relabelled — the positions figure is 992 and 994 per 1,000. The probe now prints the
positions row itself, so the claim comes from the output rather than from arithmetic on it.

## Is the check able to fail?

**Three ways, all shown.**

1. **On real data:** with the comparison pointed at the wrong field of the head
   (`non_reference_reads`), the probe reports 197,051 of `SRS3394606`'s 197,141 one-base records
   disagreeing and aborts on its own assertion. That shows the operand is live.
2. **On the shape the claim is about:** the shipped test
   `a_partial_witness_at_a_one_base_locus_shows_up_as_a_disagreement` builds the counter-example
   and the measurement reports it — one disagreement, one partial witness. Two mutations that
   both tomato stores cannot distinguish from correct code fail this test.
3. **On an empty walk:** a store holding no records now fails the run naming the file, rather
   than printing zeros and exiting 0.

## Validation results

In the container, on this worktree:

- `cargo test --release --example ng_window_coverage_probe` — **5 passed, 0 failed**.
- The probe over the six-accession slice and over the 80-interval store — the tables above; exit
  0 both times, so neither assertion fired.
- `cargo build --release --example ng_window_coverage_probe` — clean;
  `cargo clippy --release --example ng_window_coverage_probe` reports nothing for this example.
- `rustfmt --check --edition 2024` clean on both changed files.
- `cargo test --lib --all-features` — unchanged at 6,343 by this step, which touches one example
  and one doc comment.

**The standing calling oracle was not re-run** and cannot have moved: this step adds an example
and edits a doc comment, and nothing calls the module yet.

## Tradeoffs and follow-ups

- **One link in the chain is verified but not pinned.** The probe builds bodies through the psp
  reader; a psp encoding that normalised a full-span partial witness to a whole one would hide the
  counter-example shape from any future run. A reviewer round-tripped such a record through
  `PspWriter`/`PspReader` and the probe caught it, so it holds today. Shipping that as a test costs
  about thirty lines of header construction in an example, and the claim belongs to
  `src/ng/psp/`; it is recorded here rather than pinned.
- **The probe decodes bodies against a live set the calling run will not have.** `RecordIter`
  applies each head's chain-id changes before decoding; the run's own `build` decodes against an
  empty live set, which is right only while the encoder writes no chain ids. The two agree today
  and stop agreeing at the psp path's Milestone E, at which point this measurement has to be
  retaken rather than cited. The probe's module doc says so.
- **The probe builds every body, which is the expensive walk** — 7.6 million of them on the wider
  store. That is deliberate: the whole question is what the evidence says. It is why this is an
  example rather than a test, and the stores it walks are generated rather than checked in.
- **`num_obs_along_locus()` allocates a vector per record** and only its first entry is read at a
  one-base one. Plan step C3 needs every position's depth anyway, so whether to ask
  `src/ng/locus_generation/` for a borrowing form is a question for that step.
