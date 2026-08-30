# Review — ng parameters file, C2 first half: `validate`

**Date:** 2026-08-30
**Scope:** the working-tree diff of C2's first commit — `validate.rs` (new), the `Meaningless`
variant, and the inverted test.
**Verdict:** Request changes — **0 Blockers, 4 Majors, 14 Minors**, from two agents in isolated
worktrees at `fa293d2a` with the step applied as a patch: one on correctness and design fidelity,
one reading the refusal messages as the geneticist and fact-checking every claim in the new prose.

---

## 1. Confirmed rather than found

- **Every refusal the plan names is implemented and every one fires.** The correctness agent built
  the table and matched each to its test.
- **The four predicates re-derived independently all hold**: `[0, 1)` for an inbreeding coefficient
  against `InbreedingF::try_new`, which rejects exactly one as a separate error because the prior
  multiplies its heterozygote branch by `1 − F`; an alternative concentration of exactly zero
  accepted, per spec §3.6; the NaN reasoning; and `i64::MAX` as the writer's real saturation value,
  against `a_toml_integer`.
- **The eight refusals the plan did not name are each defensible**, including the dense
  read-group-id check, which C2's own sentence in the plan asks for.
- **No panic path**: no indexing, no `unwrap`, no arithmetic that can overflow, outside tests.
- **Every measured figure in the tolerance's justification reproduces**, all six.

## 2. The Majors

**Ma1 — a refusal named a real key holding a different number.** The share-smoothing paths were
`shares_origin.shorter_share` and `shares_origin.fall_off`; the file's keys are
`shorter_share_smoothing` and `fall_off_smoothing`. Worse than an unfindable path: `shorter_share`
is a sibling key of that same row holding a perfectly good number, so a reader searching it stops
at the wrong thing. The level side was already correct, so this was an oversight rather than a
convention. **The test written to prevent this could not**: it compared only the last segment, and
`curve_weight` exists. It now checks every segment of every path, and the two nested families are
in its list.

**Ma2 — the dense axis was declared and nothing measured against it.** Four tables are keyed by the
read-group id — calibration, contamination, sequencing batches, slippage groups — and none was
checked for covering `0..n`, nor for duplicates. This writer cannot produce such a file; a hand edit
can, and the symptom is silent, because the projection reads these by *keyed lookup*: the missing
id's slot becomes a defaulted scale of one and the file's claim that the library was fitted
vanishes. Fixed with one helper over all four, plus the same treatment for the two per-sample tables
against `fitted_from.samples`.

**Ma3 — "every float that is not finite" overstated the walk by most of the repeat-tract section.**
Eleven floats were never reached: both `expected_slipped_reads`, and every numeric field of both
curve kinds. All are reachable — the writer emits `nan` and the reader takes it back. Now checked,
with a test, and the doc says what the curves are (provenance, not terms in a score) and why they
are refused anyway.

**Ma4 — two prose claims were false.** That recovering a line would require a second parse: the
`toml` crate re-exports `serde_spanned::Spanned`, which carries a span out of the parse already run.
The decision not to use it stands — the wrapper would sit on the *shape*, where it reaches
`Serialize`, `PartialEq` and the round-trip equality goal 1 rests on — but the stated reason was
wrong and is replaced by the real one. And "the projection is the one caller that must, and it calls
this first" describes a caller that does not exist; written in the present indicative it told a
reader that runs are protected by these checks today. They are not, and the header now says so.

## 3. The Minors, grouped

**Applied — correctness.** `ContaminationView::was_measured` is a conjunction, so *not measured* is
`markers == 0 || reads == 0`; the check refused only the conjunction, and the row it missed is the
worse one — zero markers with 90,233 reads says *measured, 3.1 in 100* in the file and reads back as
never measured. The contamination share is now half-open at one, matching
`FrozenContamination::new`'s own assert, which was otherwise reached as a panic several frames later
naming a read group rather than a file.

**Applied — the messages.** Every row path now spells its keys as the file spells them —
`[period = 2, reference_repeats = 6]` rather than `[period 2, reference_repeats 6]` — which makes
each path a literal substring of the file a reader can paste into a search. The two bare index forms
became `[read_group = 0]` and `[sample = "TS-1"]`, since there is no literal `by_read_group[0]`
anywhere in the file and the code keys on the id rather than the array position. The read-group gap
message quoted Rust's `0..3` range syntax and pointed at a *healthy* row — "with 2 where 1 should
be" — rather than naming the missing id; it now names it, and lists the ids it did find, capped
because that axis is cohort-sized. Three messages that stated a fact now say what to do: a future
format version says the file came from a newer build, the inbreeding bound is given in words rather
than as `[0, 1)`, and a substitution rate is called a probability rather than "a share".

**Applied — the tolerance's prose.** "A million times" was 1.5 million; "there is no value between
the two that would be a closer call" was simply false, with thirteen orders of magnitude between
them.

**Recorded, not applied.** The saturation marker is refused on the two count kinds and not on the
nine other integers the writer can saturate — `reference_repeats`, `cells`, `strata`, and the two
`fitted_*_repeats` on each curve kind. Those are provenance rather than evidence and the header now
says which are checked. Also: a warranted value may still carry `observations = 0`, where the
shape's own doc says a count is absent rather than zero; consequence is a wrong number in a report
rather than a wrong score.

## 4. Mutation testing

Nine mutations, all nine fail a test. One survived a first pass — see the implementation report §8;
the fixture for the spectrum-parity rule was also too short, so the length rule caught it and the
parity rule never ran.
