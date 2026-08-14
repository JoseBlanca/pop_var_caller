# ng census — implementation plan 1: the vocabulary and the encoding

**Status:** draft, 2026-08-14. The build order for a set of changes to code that already runs on real
reads: the rename of the joint route's records to **census evidence**, the deletion of the
coverage-by-window summary, and four changes to how a position's evidence is encoded. Design is
settled in [`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md),
[`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md),
[`parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md) and the two architecture
companions. **This plan turns that design into build order; it is not a place for new design.**

**Nothing here is new machinery.** Every step edits code that exists and is under test, which is why
the plan is short and why its oracles are unusually strong: milestone A must not move a single fitted
number, and milestone B must move them only where the specification predicts.

---

## Scope

**In:**

- the module and type rename to the census vocabulary, and the deletion of `coverage.rs`;
- the depth ladder extended from twenty bins to thirty, as a **refinement** of the existing one;
- the depth recorded as the position's **true** depth rather than the subsampled one;
- the allele count narrowed to one byte, with the per-position depth cap refusing a value a byte
  cannot hold;
- the repeat-tract record's two per-locus counts moved to the stratum, and its walked flag to a bit.

**Out (later plans):**

- **the census file and reading it by section** — [plan 2](census_file.md). Nothing here writes a
  file; the direct run holds everything, which is what the comparison against the per-sample route
  needs and all it needs.
- **the duplicated-copy class in its coupled form** — no plan yet. The class ships off
  ([`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §2.2) and the form that
  replaces it is a coupling rather than a factor, which is design work and not build order.
- **the repeat-tract estimator's borrowing on the tomato cohort**, and the one-sided borrow at the
  ends of the repeat-count axis — measurements, not code.

---

## Principles (how the order was chosen)

- **The answer-preserving work first, alone.** Milestone A renames and deletes and changes no
  arithmetic, so its oracle is byte-identical output. Bundling it with a step that legitimately moves
  a number would destroy that oracle for both.
- **Isolate a step whose failure is silent.** Three steps here produce a quietly-wrong number rather
  than a crash — the true depth, the one-byte count, and the stratum-grain counts. Each lands as its
  own commit with the fitted numbers checked before and after, so a `git bisect` can find it.
- **Types first, then implementation**, within each milestone (project rule).
- **Reuse over rewrite.** The ladder is extended by adding rungs to the existing generator, not by
  writing a second ladder; the per-sample route's twenty bins are the first twenty of the thirty.
- **Verify against ground truth.** The oracle throughout is the tomato cohort and the human trio on
  real reads, whose current numbers are recorded in
  [`joint_fit_against_truth_2026-08-13.md`](../reports/joint_fit_against_truth_2026-08-13.md) and
  [`trio_heterozygosity_excess_2026-08-14.md`](../reports/trio_heterozygosity_excess_2026-08-14.md) —
  not self-consistency.

---

## Preconditions (already in place)

- `ng-joint-fit` carries the merged work of five branches and `main`; `cargo check --all-targets` is
  green and the library suite passes 3,581 tests.
- The census records, the locus selection, the ordinary-position estimator, the contamination
  estimator and the repeat-tract estimator all exist and run on the tomato cohort and the human trio.
- `examples/ng_joint_records_walk.rs` walks alignments, fills the records and fits the cohort — it is
  the parity harness every milestone below is measured with.
- Baseline numbers to compare against exist for both cohorts, and the run scripts are in
  `tmp/run_records_{tomato_cohort,trio,hg002}.sh`.

---

## The steps

### Milestone A — the vocabulary, and one deletion

Nothing in this milestone changes an answer.

✅ **A1 — record the baseline.** Run the tomato cohort and the trio, and keep the full output. Every
step below compares against these two files.
*Depends:* —. *Source:* preconditions.

✅ **A2 — `records.rs` → `census.rs`, and the seven type names.** `KeptLoci` → `CensusLoci`,
`KeptLociDigest` → `CensusLociDigest`, `SampleRecords` → `SampleCensusEvidence`, `GenericRecords` →
`GenericEvidence`, `SsrRecords` → `SsrEvidence`, `RecordWriter` → `CensusWriter`, `RecordError` →
`CensusError`, and `RecordIdentity`/`SelectionIdentity` → `RecordingTerms`/`SelectionTerms`.
*Depends:* A1. *Source:* [arch records](../arch/parameter_prepass_joint_records.md) §4's table.

✅ **A3 — delete the coverage-by-window summary.** `coverage.rs` (about 560 lines), the `coverage`
field, `CensusWriter::finish`'s parameter, and the harness's `COVERAGE_DISCRIMINATOR` path.
*Depends:* A2. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §4.

✅ **A4 — assert the milestone.** Re-run both cohorts; the output must be **byte-identical** to A1.
*Depends:* A2, A3. *Source:* this plan's first principle.

> **Checkpoint A:** a rename and a deletion, with byte-identical output on two real cohorts. Pause for
> review.

### Milestone B — the encoding

Every step moves numbers. Each says where, and each is its own commit.

✅ **B1 — extend the ladder to thirty bins, as a refinement.** Ten more rungs at the existing ratio
carry the top from 124 to about 1,500; the first twenty bins keep their exact edges, so a census code
maps to a per-sample-route bin by collapsing everything above 124. **Own commit, do not bundle.**
Assert the refinement property in a test: every edge of the twenty-bin ladder is an edge of the
thirty-bin one.
*Depends:* A4. *Source:* [generic spec](../spec/parameter_prepass_generic.md) §4;
[records spec](../spec/parameter_prepass_joint_records.md) §2.2.

☐ **B2 — record the position's true depth, not the subsampled one.** The per-position cap still thins
the allele counts proportionally; the depth code stops being clipped by it. A consumer needing the
counts' own denominator computes `min(depth, cap)`, the cap being in the recording terms.
**Own commit, do not bundle** — this is a silent-failure step: it changes what a code means, and a
consumer that keeps reading the old meaning produces a wrong rate rather than an error. The oracle is
the trio, where 999 in 1,000 positions previously sat at exactly the cap.
*Depends:* B1. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §2.2, §4.1.

☐ **B3 — the allele count becomes one byte, and `DepthCap` refuses a value above `u8::MAX`.**
**Own commit, do not bundle** — a saturating count is a wrong allele fraction, not a panic. Assert
the refusal at construction and assert that a position at the cap round-trips.
*Depends:* B2. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §2.1, §2.2;
[arch records](../arch/parameter_prepass_joint_records.md) §1.2, §3.

☐ **B4 — the repeat-tract record's grain.** `covering_not_crossing` and `bases_compared` become one
value per (read group × stratum); `walked` becomes a bit per locus. The estimator's accumulator moves
into the writer. **Own commit, do not bundle** — the counts are a denominator, and a wrong denominator
is a plausible rate.
*Depends:* B3. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §3;
[arch records](../arch/parameter_prepass_joint_records.md) §1.4.

☐ **B5 — assert the milestone, and report what moved.** Re-run both cohorts and report every fitted
number against A1, **per parameter, never pooled**. The specification predicts change confined to
positions above 124 reads a position: the trio should move and tomato, at 2.4 to 30.6 reads, should
barely. Anything else that moves is a defect in one of B1–B4 and the per-step commits are what find
it.
*Depends:* B1–B4. *Source:* this plan's second principle.

> **Checkpoint B:** the encoding changed, with a per-parameter account of what moved on two real
> cohorts and why. Pause for review.

### Milestone C — one defect the milestone-A review found

**ADDED 2026-08-14 (owner), after Checkpoint A.** It is not encoding work and it was not in the
original scope; it is here because the review that milestone A required found it, it is a wrong
answer rather than an untidiness, and leaving it in a list of open items is how it gets lost.

☐ **C1 — a repeat tract's mismatching bases are attributed to the wrong reads.** When several reads
cross a tract, `CensusWriter::add_ssr` numbers each mismatch's read from zero **within one
observation** rather than within the locus, and the walk folds reads carrying identical bases into
one observation. So two reads that differ from the reference in different places are both written
down as read 0. **The distinction that destroys is the only reason the field exists**: the same
substitution at the same place on two reads is an allele, and two substitutions on one read is a bad
read ([records spec](../spec/parameter_prepass_joint_records.md) §3, §7.3;
[arch records](../arch/parameter_prepass_joint_records.md) §1.4). Carry a per-locus read counter
across the observation loop. **Own commit, after milestone B** — it changes what the census records,
which milestone B's own oracle forbids while B is running.
*Depends:* B5. *Source:* the A2 review,
[census_rename_a2_2026-08-14.md](../../reports/reviews/census_rename_a2_2026-08-14.md) finding B1,
which carries a test that fails against the code as it stands.

✅ **C2 — a mismatch outside the tract is not recorded, and the documentation says so.**
`TractDifference::offset` was documented as negative in the sequence before the tract and past the
end in the sequence after it, while the writer only ever compared a read against the tract itself.
**Settled 2026-08-14 (owner): the promise is withdrawn, and the reason is what the sequence either
side of a tract is for.** It sits in the locus only because the aligner needs it to anchor a read;
it is ordinary non-repetitive sequence and is not part of the locus, so what happens there belongs
to the generic path at the positions it already keeps. Landed **before** milestone B rather than
after, because it is a doc comment and a test and changes no recorded value. The test that stood
for the property asserted a hand-written `-2` was below zero and called no production code; it is
replaced by one that drives the writer.
*Depends:* —. *Source:* the same review, finding B2.

> **Checkpoint C:** two defects the census's own tests could not see, fixed. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | byte-identical output on the tomato cohort and the human trio, against the A1 baseline |
| B1 | a test asserting the twenty-bin ladder's edges are a subset of the thirty-bin one's; the per-sample route's cell count unchanged at 583 |
| B2 | the trio's fitted numbers move and tomato's barely do — the split the specification predicts |
| B3 | construction refuses a cap above 255; a position at the cap round-trips; allele fractions unchanged |
| B4 | the repeat-tract estimator's fitted numbers unchanged on both cohorts, the counts being a denominator it already summed |
| B | a per-parameter table of what moved, on real reads, against A1 |
| C1 | a test giving one tract two reads that differ in different places, and asserting they come back as two reads — it fails against the code as it stands |
| C2 | a test that drives the writer rather than comparing two values it wrote itself, whichever way the owner settles it |

---

## Out of scope (next plans)

- **The census file, its directory, and reading one section at a time** — [plan 2](census_file.md),
  which is also where `CohortCensusEvidence` and the scoped lending live.
- **The duplicated-copy class in its coupled form** — needs a specification first
  ([`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §2.2 says why the
  factorised form was withdrawn).
- **The fourth shape** — positions where a quarter of the reads disagree in every sample, which the
  benchmark trio showed and no class models
  ([`trio_heterozygosity_excess_2026-08-14.md`](../reports/trio_heterozygosity_excess_2026-08-14.md)).
