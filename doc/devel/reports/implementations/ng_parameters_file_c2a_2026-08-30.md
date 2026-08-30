# ng parameters file — C2, first half: refusing a file that parses and means nothing

**Date:** 2026-08-30
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone C, step C2 — **the
first of two commits**, see §2
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §3, §5, §6, §9
**Code:** [src/ng/calling/parameters_file/validate.rs](../../../../src/ng/calling/parameters_file/validate.rs) (new), and the `Meaningless` variant in [mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)

---

## 1. Plan

C2 is "the file shape → `RunParameters`, and the reader's `validate`", and the second half is the
owner's ruling of 2026-08-28: no step owned refusing a file that parses and says nothing, so C2
does. The plan lists the refusals; the owner added one more, an evidence count sitting at the
writer's saturation marker.

## 2. Why this is two commits, recorded rather than assumed

The plan's C2 is one step. It is landing as two because the projection needs **constructors on two
types outside this module** — `RunParameters`, whose fields are private and whose only constructor
takes the fit's raw inputs rather than assembled values, and `StratumFits`, whose constructor takes
the fit's own outcome types. `validate` needs neither and has its own failing test waiting to be
inverted, so it stands alone cleanly. The projection follows in its own commit.

**Nothing in a run calls `validate` yet**, and the module header says so. The projection is the
caller that must run it first; until that lands, reading a parameters file through this module's
public entry point does not validate it.

## 3. What it refuses

Every refusal the plan names, plus eight the plan did not — checked by review and each found
defensible: a format version this build cannot read, a ploidy of zero, an empty read-group table,
a gap in the read-group ids, a calibration multiplier at or below zero, any float that is not
finite, a negative `fall_off`, and non-positive concentrations.

**`fall_off` is checked for being a number and not bounded above**, because neither this file's
shape nor the fit that produces it documents an upper bound. A value above one would say a
two-repeat slip is likelier than a one-repeat slip — implausible, but nothing has established it is
impossible, and refusing a fit nobody has bounded would reject real data to enforce a guess.

## 4. Three things measured

- **The fixture's own length spectra sum to exactly one**, as do `[0.15, 0.7, 0.15]`,
  `[0.05, 0.1, 0.7, 0.1, 0.05]` and three thirds. **An earlier version of the tolerance's
  justification said they missed by one unit in the last place; that was recalled, not measured,
  and wrong.** What the tolerance is actually for is the fit's own normalisation: over geometric
  spectra normalised to one, 3, 5, 9 and 41 offsets sum exactly, 21 offsets miss by 1.1e-16 and 101
  by 6.7e-16. The tolerance of 1e-9 sits 1.5 million times above the worst of those and ten million
  below the smallest edit a person makes to a share.
- **`ContaminationView::was_measured` is `markers > 0 && reads > 0`**, so *not measured* is a
  disjunction. The first draft refused only both counts zero; a row with zero markers and 90,233
  reads says *measured, 3.1 in 100* in the file and reads back in memory as never measured. Now
  refused.
- **`FrozenContamination::new` asserts a half-open `[0, 1)`** — "a whole library of another
  individual's DNA is not a sample of this one". A share of exactly one was accepted here and
  became a panic several frames later, naming a read group rather than a file.

## 5. Changes made

- **New** `validate.rs`: `ParametersFile::validate`, nine per-section checks, seven helpers, and
  24 tests.
- **`mod.rs`**: the `Meaningless { field, problem }` variant, carrying the key's path in the file's
  own spelling rather than a line — and `line()` returns `None` for it, which is the first thing to
  make that `Option` reachable. The test A3 left waiting was inverted: both shapes still *parse*,
  which is the half it now pins, and both are then refused.

## 6. Tests

**79 → 103 in the module.** 24 new.

## 7. Validation

All in the dev container.

- `cargo test --lib ng::calling::parameters_file` — **103 passed, 0 failed, 2 ignored**.
- `cargo test --all-features --lib --tests` — **5,026 lib tests plus every integration binary, 0
  failed**.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean. `cargo fmt --check` — clean.

## 8. Mutation testing

Nine mutations, **all nine fail a test**: inbreeding closed at one; the spectrum parity check
removed; the tolerance loosened to 0.5; the saturation marker off by one; the all-unmeasured table
allowed; both-counts-zero allowed; a calibration of zero allowed; the parity check removed again
against its clippy-rewritten form.

**One survived a first pass and it is the interesting one.** Removing the parity rule — a length
spectrum must have an *odd* number of shares — left all 100 tests green, because the fixture for
that half was `[0.5, 0.5]`, which is also below three, so the length rule caught it and the parity
rule never ran. The fixture is now four equal shares.

## 9. What the review found

Two agents. **0 Blockers, 4 Majors, 14 Minors**, written up in the review report. The two that
mattered most:

- **A refusal named a real key holding a different number.** The share-smoothing paths were built
  as `shares_origin.shorter_share`, where the file's key is `shorter_share_smoothing` — and
  `shorter_share` is a *sibling key of the same row*, a plain float, in range. A reader would have
  searched, found a healthy number, and stopped. The test that exists to prevent exactly this
  compared only the path's last segment, so `curve_weight` passed it; it now checks every segment.
- **The dense read-group axis was declared and then nothing was measured against it.** Four tables
  are keyed by that id and none was checked for covering it. A hand edit deleting one calibration
  row passed, and the projection reads these by *keyed lookup* — so the row does not shift, it
  becomes a defaulted scale of one, and the file's claim that the library was fitted disappears
  with no message. That is spec §5's third row arriving through the back door.

Also applied: eleven floats in the repeat-tract provenance were unchecked though the header claimed
"every float"; the read-group gap message quoted Rust range syntax and pointed at a healthy row
rather than the missing id; every row path now spells its keys as the file does (`[period = 2,
reference_repeats = 6]`), which makes each one a literal substring a reader can paste into a search;
and three messages that stated a fact now say what to do.

**Two prose claims were wrong and are corrected**: that recovering a line number would need a second
parse (the `toml` crate re-exports `serde_spanned::Spanned`, which carries a span out of the parse
already run — the decision not to use it stands, but the reason given was false), and that the
projection "is the one caller that must, and it calls this first" (there is no such caller yet).
