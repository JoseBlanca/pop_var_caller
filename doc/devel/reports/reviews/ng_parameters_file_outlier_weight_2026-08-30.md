# Review — the repeat-tract outlier weight, wired from the file to the scoring row

*2026-08-30. Two agents in isolated worktrees at `ab67df12`, each handed the step's diff as a
patch: one a correctness pass, one a design-fidelity pass against
`doc/devel/ng/spec/parameters_file.md` and `read_likelihoods.md` §4.5. **1 Blocker, 3 Majors and
8 Minors applied; 1 Major recorded for step D3.***

Both passes agreed the wiring is complete rather than half done — the number reaches the scoring
row, the run's report and the file from one field, so the three can no longer disagree — and both
found the same range defect at the file's edge.

## What was wrong, and what the fixes were

**B — a broken intra-doc link, which this crate denies rather than warns about.** The step moved
`DEFAULT_OUTLIER_WEIGHT` out of `repeat_tract_parameters.rs`'s imports and into its test module,
and the file's `//!` header still linked the bare name. `Cargo.toml` sets
`broken_intra_doc_links = "deny"`, and **neither `cargo clippy` nor `cargo test` runs rustdoc**,
so the step's own green validation could not see it. Fixed with the full path, and **measured
rather than reasoned about**: `cargo doc --no-deps` reports **25 unresolved links on this tree,
which is the pre-existing baseline** — it was 26 with the broken link, and the second one this
step had added (`inference::TractScoringFits`, whose type is one module deeper) is fixed too.

**M1 — deleting the forwarding call from one arm of `view()` failed no test.** `RunParameters::view`
has two arms, an uncontaminated one and a contaminated one, and the fixture that checked the
supplied weight survived the trip had an empty contamination map. So removing
`.with_repeat_tract_outlier_weight(...)` from the *other* arm left the whole library suite green,
and every run whose fit found a contamination fraction would have scored under 0.01 while its
report and its file said otherwise. That is the failure the step exists to prevent, one arm over.
A contaminated run now makes the same assertion.

**M2 — `validate` accepted two values the scoring row panics on, and the new doc said it did not.**
The weight went through `a_share`, which accepts the closed `[0, 1]`; the row asserts strictly
inside 0 and 1, and so does `RepeatTractOutlierWeight::supplied`. A file saying
`repeat_tract_outlier_weight = { value = 0.0, warrant = "supplied" }` therefore passed validation
and panicked several frames later naming a locus rather than the key — **and zero is the single
most likely edit to this key**, because turning the junk term off is what a number called an
outlier weight invites. The same shape as the contamination share of exactly one that step C2a
already moved earlier. The key now has its own open-interval check, and both endpoints have a
test.

**M3 (design) — the file could spell two warrants the memory shape cannot hold, and a third state
that is simply false.** `Warrant` in the file has four states because most of its numbers need
four; nothing fits this one, so `fitted_here` and `borrowed` are claims no run could make, and
`RepeatTractOutlierWeight` has nowhere to put them. **The state worth catching is the third**: a
person who takes spec §1.2 goal 3 at its word — copy the file your run wrote, change one line —
changes the *number* and leaves `warrant = "defaulted"` above it, which then says the run
inherited 0.01 beside a value that is not 0.01. `validate` now holds this one key to two of the
four, and the refusal for the edited-number case names the fix the file's own header already
teaches: change the warrant beside it to `supplied`.

**Eight Minors, all of them prose that the change made untrue.** Among them: the field doc on
`SsrLocusParameters::outlier_weight` still told a caller to pass the constant; `StatedConstants`'s
doc still said "the run's own value is `defaulted`", which is the thing the step undid; the two
module headers arguing that folding this number into the per-cell warrant would mark every tract
`Defaulted` now had to name `Supplied` too (the ruling survives — `Supplied` is one rung up a
four-rung ladder — but the sentence named the wrong warrant for half the cases the step made
reachable); a dropped line-continuation left fourteen spaces inside an assertion message; and one
count in a doc comment said thirty where `grep -rn "FrozenParameters::new(\|FrozenParameters::uncontaminated("`
returns **27**, of which this step forwards through 1.

## Recorded, not fixed — one Major that belongs to step D3

**D3's wholesale demotion is a *promotion* for this one number.** Spec §2.1 demotes every number
of a file whose census binding does not match to `Supplied`. `Provenance::strength` ranks
`Defaulted` 0 and `Supplied` 1, so a run that inherited 0.01 and read a mismatched file would come
out claiming a *better* warrant than it went in with — and `supplied(0.01)` is a legal
construction, so nothing would complain. The fix is one function
(`Provenance::weaker_of(file_warrant, Supplied)`, a no-op for every other number), and it belongs
where the demotion is written. Carried into `PROJECT_STATUS.md` as D3's.

## Two things both passes confirmed, and one they could not

**No production path reads the constant any more.** Every surviving mention of
`DEFAULT_OUTLIER_WEIGHT` outside `RepeatTractOutlierWeight::defaulted()` is inside a `#[cfg(test)]`
module or a doc comment, and `src/ssr/` — the separate, frozen caller — has none at all.
`RunParameters::view()` is the only production builder of a `FrozenParameters`, and both arms now
forward the value under test.

**The new guard's range is the scorer's exactly**, neither looser nor tighter: `supplied`'s
`> 0.0 && < 1.0` and the row's are the same comparison, and `is_finite()` is redundant beside them
rather than tightening.

**What no file-read check can close**: the row also asserts `weight + contamination fraction < 1`
per read group, and the fractions are the fit's rather than the file's. A supplied 0.9 passes
validation and fails at the first locus of a contaminated run. `supplied`'s doc now says so
instead of claiming a validated file cannot reach its panic.

## Tests

5,474 → **5,475** library tests, 0 failed — one new test
(`an_outlier_weight_whose_warrant_no_run_could_mean_is_refused`), with the rest of the step's
coverage added as assertions inside the three existing tests that already owned those seams.
`cargo clippy --all-targets --all-features -- -D warnings` clean; `cargo doc --no-deps` back to
its 25-link baseline.
