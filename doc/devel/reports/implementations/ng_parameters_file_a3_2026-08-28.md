# ng parameters file — A3: absence is a missing key, never a sentinel

**Date:** 2026-08-28
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone A, step A3 — the last step of Milestone A
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §5
**Code:** [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)

---

## 1. Plan

Spec §5 lists five states a reader that collapses them will get wrong, and says of itself that it
is "the one most likely to cost a day". A3 expresses all five in the types, as **missing keys
rather than values standing in for one**.

## 2. What was already true, and what A3 changed

Three of the five were settled by A1 and A2 and needed nothing:

| state | how the file says it | since |
|---|---|---|
| a read group's error rate could not be fitted | a multiplier of 1.0 whose `warrant` is `defaulted` | A1, A2 |
| a stratum was furnished from its period's curves | no row in `length_spectrum_by_stratum` | A1 |
| a slippage group put no read in a stratum | no row for that pair | A1 |

The two contamination states are A3's work:

- **`ParametersFile::contamination` is `Option<Contamination>`.** An uncontaminated run writes no
  section at all, where a table of zeros would say every library was measured and found clean.
- **A row's fraction, its two evidence counts and its source moved into a
  `ContaminationMeasurement`, held as `Option`.** A read group that identified nothing has **no
  measurement**, where in memory it has a zero fraction beside two zero counts and
  `ContaminationView::was_measured` has to be asked which it is.

**That second change also drops a wart the in-memory type documents and cannot fix.** A read group
that identified nothing still has to carry a `ContaminationSource` there, and neither variant is
true of it — `run_parameters.rs` says so in as many words. Here it has none, because it has no
measurement to carry one on.

## 3. Assumptions

- **Two shapes remain writable that no run should mean**, and neither can be made unspellable: a
  contamination table whose every row is unmeasured (the uncontaminated run, longhand), and a
  measurement whose two counts are both zero (the in-memory `UNMEASURED_READ_GROUP`). No type can
  say a `Vec` is non-empty, or that not every row of one is absent. **Both are step C2's to
  refuse**, and `the_shape_accepts_two_things_step_c2_must_refuse` records that they are accepted
  today so the gap is visible rather than implied.
- **`NonZeroU64` on the two evidence counts was considered and not taken.** It would make the
  second shape unwritable, but a fit that returns an estimate with no evidence behind it would
  then have no file to be written to at all — and whether that can happen is a question about the
  contamination estimator, which step C4's round trip on a real fit is what answers.
- **The empty-list refusal is C2's, not C1's.** C1 is parsing, and both unwanted shapes parse.

## 4. Tests

Thirteen became fifteen, one ignored. The step's own test is
`each_of_the_five_states_is_a_missing_key_and_not_a_value`, one part per row of §5.

**Its first draft could not fail on the hazard its own message named**, and a review proved it by
writing the collapse: a writer that fills the (stratum × slippage group) axis with zero rows —
exactly what §5's fifth row forbids — added three rows to the file and left the test green. Part 5
had been asserting a property of `Vec::pop` over data the test itself had just edited, and never
looked at the emitted document.

**Both parts 4 and 5 now assert against the emitted document, on the gap the fixture already
has** — three strata crossed with two slippage groups is six pairs, of which three carry a row,
and one of the three strata carries a length spectrum. Re-running the same zero-row writer against
the fixed test fails it on the assertion that names the hazard. Emptying a table and checking it
came out empty is a weaker claim, and it is the one the first draft made.

## 5. Validation

| command | result |
|---|---|
| `cargo fmt -- --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | **14 passed, 0 failed, 1 ignored** |
| `cargo test --lib` | **4,934 passed, 0 failed, 12 ignored** |
| `cargo doc --no-deps` | zero unresolved links in this module |

## 6. Two things Milestone A cannot settle, raised at Checkpoint A

Both come from the design-fidelity review, which read spec §3 section by section against the types
and found them complete, and then looked for a sixth state the spec distinguishes and the types
collapse.

- **Spec §2.1's wholesale demotion has nowhere to write itself.** A file whose census binding does
  not match is demoted — "every number in it is marked `Supplied`" — and only five numbers carry a
  warrant that can say so. The slippage numbers, the contamination fraction and the prior's
  concentrations carry a different kind of warrant, and none of those vocabularies has a word for
  *a person handed me this*. The in-memory `LevelProvenance` has the same gap.
- **One axis of the file grows with the cohort and §9 does not price it.** The repeat-tract
  substitution rate is keyed by (read group × stratum × ploidy); §3.7 says "per stratum" and §9
  prices three axes, none of them this one.
