# ng parameters file — F1 correctness review

**Date:** 2026-08-31
**Step:** [parameters_file.md](../../ng/impl_plan/parameters_file.md) F1
**How it was run:** one agent, in a worktree detached at `ede29317` with the step applied as a
patch, building through the **worktree's own** `scripts/dev.sh` (every log carries
`Compiling pop_var_caller v0.1.0 (/Users/jose/devel/pop_var_caller-f1-correct)`). Baseline before
any change: 220 passed, 0 failed.

**Verdict: no blocker. One Major, seven Minors, and one surviving mutation out of sixteen.**

---

## Major — a file whose every number is `supplied` opened by claiming a fit

`what_was_fitted.rs`, `was_fitted`. Four of the seven groups answered from a table's **presence**
rather than from the warrant, so after spec §2.1's demotion the file a run writes opened with
*"5 of the 7 groups of numbers in this file were fitted from your data"* over rows whose every
warrant is `supplied`. **Run**, on the exact write path `to_run_parameters_for` takes when a census
disagrees; the same probe printed the substitution-rate rows it had just written as `[Supplied]`.

The same thing happened with no demotion at all, on spec §7's *first* source: a run scoring from a
hand-written file whose warrants are `supplied` wrote back the same sentence.

**The design-fidelity and reader agents found this independently**, from the spec's side and from
the artefact's. See the fixes report for what was done: the substitution-rate half is fixed, the
other four are the format's limit and are now disclosed in the file rather than papered over.

## Minor — the fitted-with-no-count guard stopped running on two of the three sources

`from_run_parameters.rs`. Before F1, `calibration_rows` asserted on **every** `of_run` call that a
read group with no count must be `Defaulted`. It now lives inside
`ReadsBehindEachCalibration::of_the_fits_rates`, so the other two constructors bypass it. **Run**:
`of_run` with `nothing_was_fitted(3)` against a run whose calibrations are
`FittedHere`/`Defaulted`/`Supplied` wrote three rows and `validate()` accepted them, where before
F1 it panicked.

**Declined, with a reason** (see the fixes report): the guard is about *the fit's* rate set, and
re-adding it for the file source would panic on a legal file — the same file this reviewer's own
"checked, and sound" section shows the new path round-trips and the old one did not.

## Minor — two different VCFs could be given one parameters file

`written_beside_the_vcf.rs` used `to_string_lossy`. **Run**: `/data/\xff\xfe.vcf` and
`/data/\xfe\xff.vcf` both produced `/data/<U+FFFD><U+FFFD>.parameters.toml`, so the second run
silently overwrote the first's parameters and neither file's name matched its VCF. **Fixed** — the
derivation is `file_stem`/`extension` on `OsStr` now, and `two_names_that_are_not_text_do_not_collide`
pins it.

## Minor — degenerate paths

**Run**, all of them: `""` → `.parameters.toml`; `"."`/`".."` → `./.parameters.toml`,
`../.parameters.toml`; `"/"` → `/.parameters.toml`; `".vcf"` → `.parameters.toml`; `"/data/"` (a
directory) → `/data.parameters.toml`, *beside* the directory rather than in it. **Partly fixed and
partly documented**: `.vcf` is now `.vcf.parameters.toml` (a hidden file called `.vcf` is not an
empty name with an extension, and `Path::extension` agrees), and the module header states that the
argument is a VCF *file's* path and what the rest give.

## Minor — `calls.vcf.gz` and `calls.bcf` share one parameters file

Reasoned only. The module header claimed *"Two VCFs of the same cohort in one directory keep two
parameters files"*, which holds for two stems and not for two formats of one. **The sentence is
fixed**, and says which case it does not cover and why that is usually right.

## Minor — `GroupOfNumbers::EVERY` claimed a compiler guarantee it does not have

Reasoned only. Adding a variant breaks the four exhaustive `match`es, not the `[Self; 7]` literal —
and every test in the module iterates `EVERY`, so a variant missing from it drops out of the
denominator with the suite green. **The doc is fixed** to say exactly that, and to say *add the
variant here first*.

## Minor — `write_beside_the_vcf` truncated in place

Reasoned only. `fs::write` truncates before writing, so a write that fails part-way leaves a
truncated parameters file beside a *complete* VCF, after every locus has been called. **Fixed** —
a temporary file in the same directory, renamed over the destination.

## Minor — the real defaults run's opening line was never asserted

`defaults.rs` builds a file from an actual `RunParameters::of_defaults` and checked nine phrases of
its prose, none of them the new headline, which was pinned only against a hand-built approximation
in `to_toml.rs`. The reviewer checked the real thing does land on 0 of 7. **Fixed** — the headline
is now one of the phrases that test asks for.

---

## Checked, and sound

**`warranted_value`'s two rules hold on every path.** The reviewer tabulated old against new for
every (warrant, count) pair. One row differs: `fitted_here`/`borrowed` with an **absent** count,
which the old path wrote as `observations = { reads = 0 }` and the new one writes as no count.
**Run**: such a file (which `validate` accepts) now round-trips file → run → file unchanged, where
the retired helper turned it into `reads = 0`. The new path round-trips strictly more files and
loses none.

**`census: None`.** The order is unchanged — validate, the three refusing bindings, then the census
— and `None` short-circuits before the demotion, so no path assuming a census is reached and
`demoted_to_no_better_than_supplied` is strictly *less* reachable than before.

**`what_the_run_fitted` on degenerate files.** Run: empty calibration table, empty inbreeding table,
a contamination section present with no rows, and one present with every `measurement` absent — all
match the documented rule.

**The length assertion in `of_run`** is not reachable from a legitimate caller: `of_run` has already
matched the read-group table against the run, the file source's length is the file's read-group
count which the refusals have already matched, and the defaults source is built from
`read_groups.len()`.

---

## Mutation testing — 16 applied, 15 killed, 1 survived

Every arm of `was_fitted` (7), `Warrant::was_measured_from_the_runs_reads`'s `Borrowed` arm, the
`census.and_then` in both directions, the two suffix lists in `beside_the_vcf`, and the
`nothing_was_fitted()` branch of all three prose builders. Mutations of the seven arms were run in
five separate builds so that none could mask another, with attribution confirmed from the assertion
message, which names the group.

**One survived, and it is the shape this project's mutation testing keeps finding: a branch whose
only fixtures happened to agree.** Replacing the no-count arm of `calibration_rows` with the old
`warranted_value(…, Reads(0))` passed all 220 tests — every fixture whose count is absent is a
`defaulted` row, where `warranted_value` drops the count anyway, so the two branches wrote the same
row. **The missing test is written**:
`a_fitted_multiplier_with_no_count_does_not_gain_a_count_of_zero`, a `fitted_here` multiplier with
no `observations`, round-tripped and compared whole.

Every mutation was restored, verified by `git diff` against the patch, and the last full run over
the restored tree passed 221 tests.
