# ng parameters file — E2: a run with no fit assembles from the defaults

**Date:** 2026-08-31
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone E, step E2
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §3.4, §3.5, §5, §7, §8
**Code:** `DEFAULT_INBREEDING_COEFFICIENT`, `DeclaredInbreeding` and `RunParameters::of_defaults` in
[parameters_file/defaults.rs](../../../../src/ng/calling/parameters_file/defaults.rs); two
assertions in [from_run_parameters.rs](../../../../src/ng/calling/parameters_file/from_run_parameters.rs)

---

## 1. Which door, and why not the other two

`RunParameters` has three constructors now and E2 adds the third.

- **`assemble`** takes the *fit's* raw outputs and derives the run's read-group axis from them. A
  run with no fit hands it two empty maps and it refuses a run with no read groups — the wrong
  complaint about the right situation.
- **`of_gathered_values`** takes nine already-assembled arguments, which is the whole of a run's
  parameters. It is what E2 builds on, and it is not the door a caller should use: nine arguments
  is nine chances to leave a default out.
- **`of_defaults`** takes what a run with no fit actually has — its read groups, its ploidy, and
  what it was told about inbreeding — and fills in the rest from the module's list.

## 2. The inbreeding coefficient, and the ruling that unblocked the step

`RunParameters` requires one coefficient a sample (§3.5: "At least one is required"), spec §8's
three cases have no slot for this parameter, and
[`generic::fallback`](../../../../src/ng/parameter_estimation/generic/fallback.rs)'s header forbids
the **fit** a default for it. **Owner's ruling, 2026-08-31: the run takes zero, and a user may state
one value for the whole run or a different value for any sample.**

`DeclaredInbreeding` carries the three states. A sample nothing was said about takes
`DEFAULT_INBREEDING_COEFFICIENT` marked `Defaulted`; a sample a statement reaches is `Supplied`. The
join is **by name**, which is spec §3.5's own rule for this quantity — a per-sample table joined by
row order is silently wrong against a re-ordered sample list.

**The argument for why the fit may not and the run may was wrong in its first draft, and the
geneticist reviewer is what found it.** It said the two are *different acts* — the fit *infers*, the
run *declares* — which does not survive the case that produces the file: under `nothing_said()`
nobody declared anything, and the zero is exactly the default the fit was forbidden. What is true is
**how far a wrong constant travels**: a fitted diversity divides by `1 − F` and carries the mistake
into every number the fit emits, where a defaulted coefficient at calling time reaches the calls and
stops. Zero is not harmless there either — it is Hardy–Weinberg, so a selfing cohort's heterozygotes
are over-called, and a landrace at `F = 0.9` scored at zero has that branch **ten times** what it
should be. Both `fallback.rs`'s header and the module here say which of the fit and the run they are
about now.

## 3. The slippage group is declared and empty, which are two different things

`of_defaults` declares every read group into slippage group 0 — the run's own declaration, and the
default the joint walk makes ([`ng_joint_records_walk.rs`](../../../../examples/ng_joint_records_walk.rs)) —
and gives it no strata.

**Declaring nothing would have been wrong in a way nothing would say.** `StratumFits::at` looks the
read group up *before* the stratum, so an undeclared read group answers `NoSlippage::UnknownReadGroup`
— which the type's own documentation calls *"the run is not what it claims"* and which
`TractScoringFits` counts apart from the ordinary absences for exactly that reason. Every cell of
every tract of a defaults run would have reported it. Declared, the answer is `NoSuchStratum`, which
is ordinary. Both reviewers re-derived this and it holds.

## 4. ⚑ Absorbed: the writer could not write a run that fitted nothing

Two assertions in `from_run_parameters.rs` refused a state a legal run produces, and E2 relaxed
both — narrowly, and each is still a *tightening* for every warrant but `Defaulted`:

- **at the door**, the rates map must cover every read group **or be empty** — the two-state
  contract the contamination axis already carries. A *short* map is the failure worth refusing,
  because it writes some other read group's count beside a multiplier;
- **at the row**, a missing rate is legal exactly where that read group's calibration is
  `Defaulted`, which is the one warrant that writes no count at all.

**Proved rather than argued, by the correctness reviewer**: on the module's mixed-warrant fixture
(read group 0 fitted, 1 defaulted, 2 borrowed), a rate set of the right cardinality that drops read
group 1's entry and carries a stranger at key 9 — refused before, admitted now — writes a
**byte-identical** file. The lookup is keyed by `ReadGroupId`, so a row can only read the rate filed
under its own id.

**⚑ And it leaves a gap F1 will hit.** §7 is "one writer, three sources", and two of the three now
work. A run scoring from a **supplied parameters file** has `Supplied` calibrations and no rate map,
and `of_run` panics on it — because the file carries the *counts* behind each multiplier
(`RunParametersFromFile::reads_behind_each_calibration`) and never the rates, which `of_run`'s
signature asks for. That mismatch predates E2 and is F1's to resolve; recorded in
`PROJECT_STATUS.md`.

## 5. A defaults run's file, and one defect in it that was E2's own

`a_defaults_run_writes_a_file_that_reads_back_as_the_same_run` goes the whole way — parameters →
file → TOML → file → parameters — so nothing here is proved about a file nobody could write.

**A coefficient an operator declared wrote `observations = { covered_positions = 0 }` beside it**,
which says the number was measured over no genome, while the file's own editing rule three lines
from its top tells that reader to *delete* `observations` on a value they supplied. There was an
existing rule here with a written rationale — *a count of zero is a count*, so that "fitted from
nothing" stays distinguishable from "a stated constant" — so the change is scoped rather than
sweeping: a **`supplied`** value with a zero count writes no key; a fitted one still writes its zero,
because a fit that produced a number from no reads is alarming and the count is what says so.

## 6. Tests

Eight added — seven in `defaults.rs`, one in `from_run_parameters.rs`. The ones that carry the step:

- `a_run_with_no_fit_takes_every_default` — all nine fields, field by field.
- `a_defaults_runs_tracts_find_no_stratum_rather_than_an_unknown_read_group` — §3 above.
- `a_stated_coefficient_is_supplied_and_an_unstated_one_is_defaulted` and
  `a_per_sample_statement_lands_on_that_sample_and_overrides_the_run_wide_one` — the three states
  and the name-keyed join. The second would pass every count-based check under a positional join,
  which is D2's Blocker one level down.
- `a_defaults_run_writes_a_file_that_reads_back_as_the_same_run` — §5.
- `a_rate_set_covering_some_of_the_runs_read_groups_is_refused_at_the_door` — added after review,
  because relaxing the door guard to `<=` survived all 184 tests otherwise.

**Mutations: mine nine, the reviewer's twenty-five.** Of the reviewer's, 22 killed, 3 survived: one
real-but-safe weakening (the `<=` above, now pinned), and two provably equivalent — one of which is
the fallback evidence count whose comment already carried the proof.

**⚑ One of my own mutations survived and the reason is worth keeping.** The first draft asserted
every coefficient equalled `DEFAULT_INBREEDING_COEFFICIENT`, which moves both sides together: a
build shipping the constant at 0.5 passed all 184 tests. Three assertions compare against a literal
zero now, and that build fails all three.

## 7. Review

Three agents in isolated worktrees: correctness with mutation testing, design fidelity, and the
produced file read as a geneticist. **1 Blocker, 12 Majors and about a dozen Minors**, and — as on
E1 — every finding but the Blocker was prose that said something untrue.

**The Blocker was the file telling its reader the caller is broken.** Every inbreeding row a defaults
run writes carried `origins::INBREEDING_COEFFICIENT`'s pre-ruling text: *"inbreeding has no default:
a run should not be able to write this line."* The geneticist's first action would have been to file
a bug rather than change `0.0` to `0.9`. It now says nobody stated one, what `1 − F` does, what a
landrace at 0.9 scored at zero costs, and what to type. **The file's header contradicted it a second
time**, with E1's own sentence *"do not mark an inbreeding coefficient or a substitution rate
`defaulted`"* standing above a file that marks every one of them so.

**Three of the Majors were wrong claims of mine**, all caught by running the thing:

- the mutation comment said a build at 0.5 "passed every test in this module until this line was
  written"; deleting only that line still fails three tests, because the next line kills the mutant
  on its own. What I had measured was the first draft.
- `of_defaults`'s `# Panics` credited `of_gathered_values` "one frame later"; that call is never
  reached, because `SequencingBatches::all_together` is an argument and panics first.
- `names_not_in` promised names "in the order they were given" and returns them sorted — it iterates
  a `BTreeMap`, and the one test had a single name.

**What the geneticist could not answer from the file is left to E3**, whose brief is that gap: an
empty slippage table reads as *no read group put a read in any stratum* where in fact every tract
was scored under HipSTR's shipped constants at 5 reads in 100 each way; `[fitted_from]` heads a file
where nothing was fitted from anything; and no line anywhere says no fit ran.

## 8. Validation

In the container, on the committed tree:

- `cargo test --lib` — **5,556 passed, 0 failed, 13 ignored** (5,548 before this step).
- `cargo test --lib ng::calling::parameters_file` — **185 passed** (177 before).
- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo doc --no-deps` — 25 unresolved-link errors and 23 redundant-target warnings, both unchanged.
- `cargo test --all-targets` still exits 101 on the pre-existing panic in `benches/psp_writer_perf.rs`.
