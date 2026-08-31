# ng parameters file — F1: what the three reviews changed

**Date:** 2026-08-31
**Step:** [parameters_file.md](../../ng/impl_plan/parameters_file.md) F1
**Reviews:** [correctness](ng_parameters_file_f1_correctness_2026-08-31.md),
[design fidelity](ng_parameters_file_f1_design_fidelity_2026-08-31.md),
[reader](ng_parameters_file_f1_reader_2026-08-31.md)

**After the fixes:** `cargo fmt --check` clean, `cargo clippy --all-targets --all-features -D
warnings` silent, `cargo test --lib` **5,599 passed, 0 failed, 13 ignored** (228 of them
`ng::calling::parameters_file`, against 194 at `ede29317`), `cargo doc --no-deps` at exactly the
48-line pre-existing baseline.

---

## 1. The one all three found, and how far it could be fixed

**All three reviewers found the same defect from three seats** — the correctness agent by running
the demoted write path, the design agent from spec §2.1, the reader by reading a produced file. The
new opening line said *"N of the 7 groups of numbers in this file were fitted from **your data**"*,
and a file demoted under §2.1 — every warrant `supplied` — still counted five groups as fitted.
`WhatTheRunFitted`'s own doc claimed the summary *"cannot come to disagree with the warrants"*.

**Three things done, and the third is the one that matters.**

- **One predicate was simply wrong.** `SubstitutionRateRow.rate` **is** a `WarrantedValue` and
  `demoted_to_no_better_than_supplied` demotes it; the predicate read `is_empty()`. It reads the
  warrant now, which moves the gap from five groups to four.
- **The doc no longer claims what it cannot deliver.** `GroupOfNumbers::states_whose_reads` is the
  new place the limit is written down, and `WhatTheRunFitted` points at it.
- **The file discloses the gap rather than papering over it.** The headline says *fitted from
  reads*, not *from your data*, and a derived paragraph names the three groups that can say **whose**
  reads and the four that cannot, ending: *"So in a file whose numbers were fitted over a different
  cohort those 4 still read as fitted."* Derived from `states_whose_reads`, so a group that gains a
  warrant leaves the sentence by itself.

**⚑ The residual is the owner's and is a stop-and-ask before F2.** Contamination measurements,
slippage rows, the two length-spectrum tables and the prior's rung carry no *handed-over* state, so
§2.1's demotion cannot reach them. This is the same gap `demoted_to_no_better_than_supplied` already
records for `SeedRung::FittedCurve`. F1 now states it truthfully; closing it is a change to what the
file records.

Two tests hold both halves: `a_demoted_file_stops_counting_the_groups_that_carry_a_warrant` pins
exactly which three move and which four stay, and
`a_demoted_file_says_a_fit_produced_its_numbers_and_not_that_this_run_did` pins the file's own
sentence.

## 2. The falsehood F1 introduced

**The empty-census note promised a demotion that does not happen** (reader B3). It said *"a run that
reads this file will find a disagreement at the first line"*; `census_disagreement` zips two empty
term lists and finds nothing, and a run with no census passes `None` and compares nothing — so the
sentence was false for exactly the reader the line above it names. The note now gives both answers
and says which run gets which, and two assertions pin them.

## 3. The two census seams

- **`fitted_under_another_census: Option<String>` conflated *they agree* with *nothing was
  compared*.** Replaced by a three-state `CensusAgreement`, because §13 test 5 wants the demotion
  visible in a run's report and F2 is about to read this.
- **`Some(CensusIdentity::of_a_run_with_no_census())` and `None` gave opposite outcomes** — a driver
  holding one identity would have got everything demoted. An identity naming no terms is treated as
  no identity, with a test.

## 4. The artefact's own prose

- **A `defaulted` multiplier is not 1.0**, by the owner's E1 ruling, and the row comment said the
  reads were *"taken at face value"* above a legitimate 4.0 — six Phred apart from what the run
  scored them at. The explanation is in the **section** now (a defaults run at the top of the cohort
  range would otherwise pay seven comment lines a read group), the row note is one line, and the
  header stops listing the multiplier among the numbers with a built-in value.
- **The four warrant words are defined** where the headline introduces them. `borrowed` carries the
  count and was explained nowhere.
- **The all-fitted arm names the seven groups**, which was the arm making the strongest claim and
  giving the reader nothing to check it against.
- **"A zero means it was measured"** is qualified — a defaults run's own
  `inbreeding_coefficient = 0.0, warrant = "defaulted"` broke it.
- **The absent-contamination note is separated by a blank line** so it stops reading as a remark
  about the calibration table above it. TOML cannot head an absent section; this is as far as it goes.
- Lists in the derived prose take an *and* before the last item.

## 5. Code robustness

- **`beside_the_vcf` never goes through a `String`.** `to_string_lossy` mapped every invalid byte to
  U+FFFD, so `/data/\xff\xfe.vcf` and `/data/\xfe\xff.vcf` shared one parameters file named after
  neither. It is `file_stem`/`extension` on `OsStr` now. Side effect: `.vcf` is a hidden file called
  `.vcf` and becomes `.vcf.parameters.toml` rather than being eaten.
- **`write_beside_the_vcf` is atomic** — a temporary file in the same directory, renamed over the
  destination. `fs::write` truncates first, so a failure part-way left a truncated parameters file
  beside a *complete* VCF, after every locus had been called.
- **`GroupOfNumbers::EVERY`'s doc** stops claiming a compiler guarantee it does not have: adding a
  variant breaks four `match`es but not the array, and every test iterates the array.
- **`RepeatTractLengthSpectra::key()`** names both its tables.
- The module header's over-broad claim about two VCFs, the degenerate-path behaviour, and the
  ~79 MB single allocation at 3,000 samples are all stated.

## 6. The surviving mutation

Sixteen mutations, fifteen killed, **one survived** — the no-count arm of `calibration_rows`,
replaced by the old `warranted_value(…, Reads(0))`, passed all 220 tests. The shape this project's
mutation testing keeps finding: every fixture with an absent count is a `defaulted` row, where both
branches write the same row. `a_fitted_multiplier_with_no_count_does_not_gain_a_count_of_zero` is
the fixture that tells them apart.

## 7. Considered and declined, with reasons

- **Re-adding the fitted-with-no-count guard for all three sources** (correctness Minor). It would
  panic on a legal file — a `fitted_here` multiplier with no `observations`, which `validate`
  accepts and which the same review's own analysis shows the new path round-trips and the old one
  did not. The guard is about *the fit's* rate set and belongs where it is.
- **The partial-coverage fallback note** (reader Major). E3 decided deliberately that a note firing
  on a partly-covered table would be true of almost every run and tell a reader nothing about
  theirs, and recorded why. Raised at Checkpoint F rather than overturned.
- **Conditional section prose** (reader Major) — a defaults run's `[repeat_tracts]` defines fifty
  lines of vocabulary for keys that are not there. Structural, and beyond F1.

## 8. One reviewer number that did not survive checking

The design agent refuted the header's *"eleven lines of prose"* correctly — it is **39**, and the
source had already been corrected to 39 before the review returned. Its replacement figure, *106
identical lines, first difference at 107*, was measured against
`ng_parameters_file_e2_defaults_run_as_written_2026-08-31.toml`, an artefact written at **E2**
(`6e434561`) whose preamble E3 later changed; that comparison first differs at line 23. Measured on
the base commit's own golden (`git show ede29317:…/every_shape_as_written.toml`), the first
`warrant` value is on **line 105**, which is what the source says.
