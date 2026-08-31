# ng parameters file — E3: the slippage slot, and what a run does without it

**Date:** 2026-08-31
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone E, step E3
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §8, §12 question 1
**Code:** `RepeatTractFitsUsed` in [run_report.rs](../../../../src/ng/calling/run_report.rs), built by
`RunParameters::report`; two conditional notes and one absent-section note in
[to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs)

---

## 1. The trace, run rather than recalled

`StratumFits::at` looks the read group up **before** the stratum and answers with one of
`NoSlippage`'s four variants. `inference::repeat_tract_parameters`'s gather loop then does this, per
`(read group, candidate)`:

- **no slippage** → `StutterModel::hipstr_shipped()` and `Provenance::Defaulted`, counted in
  `slippage_defaulted`, and counted *again* in `slippage_defaulted_by_an_unknown_read_group` for the
  two absences that mean the run is not what it claims;
- **no substitution rate** → `DEFAULT_SSR_SUBSTITUTION_RATE` and `Provenance::Defaulted`, counted in
  `substitution_defaulted`;
- the cell's warrant is `weaker_of` the two, so `Defaulted`.

**So the traced behaviour scores the tract rather than refusing it.** E3's brief says to stop for a
ruling in that case. **The owner ruled on 2026-08-31**: reasonable numbers stand now and will be
fitted from GIAB HG002 once the caller works; §8's measurement stays owed and does not block the
plan. That ruling is recorded in `PROJECT_STATUS.md` — it authorises proceeding, and a reader of the
code that cites it could not otherwise check it.

## 2. Where the gap is now visible

**Two places, because they answer different questions.** The per-locus counts answer *how much of
this tract fell back*; neither can answer *did this run fit any slippage at all*, which is what
decides whether a reader should distrust every repeat-tract call in the file.

- **`RepeatTractFitsUsed` on the run's report** — strata with slippage, fitted substitution rates,
  and the read groups the slippage declaration does not name. That last list is empty on a defaults
  run, and deliberately: being told nothing about slippage and being unable to look it up are
  different failures, and only the second means the parameters and the reads came from different
  runs.
- **Two notes in the produced file**, fired only where a table is **empty**. A partially fitted run
  gets neither, which is not an omission: a note saying *some tracts fell back* would be true of
  almost every run and would tell a reader nothing about theirs.

## 3. The note is derived, and checking it found something

**Every number in the note is read off `StutterModel::hipstr_shipped()`**, so the sentence a user
reads cannot come to disagree with the model their tracts were scored under. It is written in the
section's own three words — `share_of_reads_that_slip` = 0.10, `shorter_share` = 0.50,
`fall_off` = 0.05 — because the section teaches that parameterisation forty lines earlier and a note
quoting four direction shares instead lines up neither with the table it sits on nor with a later
fitted run's rows.

**⚑ The reviewer's own re-expression was wrong, and checking it is what turned up the fact worth
writing down.** They summed all four direction shares to 0.12. `share_of_reads_that_slip` is
`Slippage::level`, which `stutter_rates_for` splits into the **whole-repeat** pair alone; the
part-repeat mass is `PART_REPEAT_SHARE_OF_WHOLE` (0.05) times it, **added on top**. So the figure is
0.10, and — the consequence — at a slip share of 0.10 a fitted row would carry 0.005 of part-repeat
mass where the shipped model carries **0.02, four times as much**. The shipped model is not a point
the fit could produce, so a reader comparing the defaults against their own fitted rows would find a
term that does not add up. The note says so.

## 4. What the reader could not act on, and what was fixed

The geneticist generated a defaults run's file — nothing in the tree had shown one before E2's
review — and read it. Their Blocker and Majors, all applied:

- **one pair of numbers stands in for every stratum**, so a 20-base mononucleotide run and a 5-copy
  tetranucleotide are scored alike. They read "every repeat tract" as *scope*, not as that. It is
  the fact that decides whether they drop mononucleotide calls, and the file did not supply it;
- **a direction share covers slips of any size**, so *a whole repeat short* means one repeat or more;
- **"part repeat"** appeared once in 264 lines and was defined nowhere — two of the four numbers were
  therefore uninterpretable;
- **"a PCR library slips more than this"** was overstated: it conflates a few-cycle library prep with
  25–35-cycle targeted amplification, which differ by close to an order of magnitude. Qualified by
  cycle count now;
- **HipSTR replaces these constants by fitting**, so borrowing them is weaker than "HipSTR's shipped
  constants" alone implies;
- **the preamble and the note appeared to contradict each other** — "nothing here defaults that
  number in this file" against "every tract's cells took the caller's stated 0.001". Both are true,
  because the default is taken at the tract and never written as a row, and the preamble says that
  now rather than leaving the reader to work it out;
- **`fall_off` printed as `0.050000000000000044`.** The file's `Debug` float format is right for a
  *value* that must round-trip and wrong for prose, so the three numbers render to two decimals.

**And one the design review found in the same place**: the note explaining what an absent
`[contamination]` section means lived *inside* the section, so the only reader who needed it — the
one holding a file with no such section — never saw it, and the word *contamination* appeared
nowhere in a defaults run's whole file. Absence says its own name now.

## 5. Tests

Six added. Two of them exist because a mutation survived:

- `a_defaults_runs_report_says_every_tract_falls_back`, and
  `a_defaults_runs_tract_finds_no_stratum_and_falls_to_the_shipped_model` — the trace, asserted
  through `StratumFits::at` rather than through the caller, because what the gather does with that
  answer is `repeat_tract_parameters`'s own to test and it does.
- `a_defaults_runs_file_says_what_its_empty_repeat_tract_tables_mean` and
  `a_fitted_runs_file_carries_no_such_note` — the notes fire on the state and are not a paragraph
  the writer always emits.
- **`a_run_that_fitted_a_stratum_does_not_report_every_tract_falling_back`.** Measured: hard-coding
  `strata_with_slippage: 0` — which makes the predicate answer *true* for every run in the project —
  passed all 5,563 library tests. Three tests asserted it true and none asserted it false, so the
  field said nothing.
- **`each_empty_table_note_answers_for_its_own_table`.** Measured: pointing the substitution-rate
  note's guard at the *slippage* table passed all 192 tests of the module, because every fixture had
  both tables empty or both full.
- `a_run_holding_substitution_rates_and_no_slippage_reports_both` — added in the author's own
  mutation pass, because reporting a constant zero for that count passed everything else.
- `every_share_the_note_quotes_is_a_whole_percent`, in `to_toml`, because the note rounds to whole
  percents and a shipped share that stopped being one would be shown to a reader as a number the
  model does not hold. A comment claimed this test existed before it did.

**Mutations: eight by the author, thirty-two by the correctness review, two more by the
design-fidelity review.** Of the correctness review's, 22 were killed and 10 survived; one is
provably equivalent (exchanging two format arguments that are both 0.05, so the string is
byte-identical), and its three Majors were all test gaps rather than wrong code. **Two of the three
were already closed by the time it reported** — the same two the design-fidelity review found, which
is the pair worth remembering because both were a field and a guard whose *only* fixtures happened
to agree: the field was asserted true three times and false never, and the guard's two fixtures had
both tables empty or both full.

**The third was open and is now fixed: the read-group walk's lower bound.** `(1..len)` survived the
whole library suite, because the only fixture dropped read group **2** from the declaration. A run
whose slippage fit never named read group **0** reported nothing. The fixture drops both ends now,
and that mutation fails it.

**⚑ And one Minor is recorded rather than covered.**
[`GroupNotInTheFit`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs) — a read group
declared into a slippage group past the end of the fit's own rows — is the *other* absence meaning
*the run is not what it claims*, and `TractScoringFits` counts the two together for that reason. The
run-level field covers only the first. It is a property of each `(read group, stratum)` row rather
than of the declaration, so reporting it means walking the strata; and **neither production path can
produce it**, because the file's reader densifies every stratum row to `max(group) + 1` and
`gather_strata` sizes its groups from the map it is handed. The field's doc says which of the two it
answers and why, rather than claiming both.

## 6. What E3 did not do, and where it went

`PROJECT_STATUS.md` listed three things a geneticist could not see in a defaults run's file. E3 took
the repeat-tract one and the contamination one. The other two are about the file as a whole rather
than about slippage and belong with F1, which owns the writer and is where a run knows which of §7's
three sources it came from: **`[fitted_from]` still heads a file where nothing was fitted from
anything**, and **no line says no fit ran** — the reader asked for one above `format_version`, and
*"0 of 9 parameter groups were fitted from your data"* would serve a partially fitted run too.

## 7. Validation

In the container, on the committed tree:

- `cargo test --lib` — **5,565 passed, 0 failed, 13 ignored** (5,556 before this step).
- `cargo test --lib ng::calling::parameters_file` — **194 passed** (185 before).
- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo doc --no-deps` — 25 unresolved-link errors and 23 redundant-target warnings, both unchanged.
- `cargo test --all-targets` still exits 101 on the pre-existing panic in `benches/psp_writer_perf.rs`.
