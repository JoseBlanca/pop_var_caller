# ng parameters file — F1: one writer, three sources

**Date:** 2026-08-31
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone F, step F1
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §7, and §2.1 / §6 for the census
**Code:** `ReadsBehindEachCalibration` and the new `of_run` signature in
[from_run_parameters.rs](../../../../src/ng/calling/parameters_file/from_run_parameters.rs);
`CensusIdentity::of_a_run_with_no_census` and `to_run_parameters_for`'s optional census in
[bindings.rs](../../../../src/ng/calling/parameters_file/bindings.rs); two new files,
[what_was_fitted.rs](../../../../src/ng/calling/parameters_file/what_was_fitted.rs) and
[written_beside_the_vcf.rs](../../../../src/ng/calling/parameters_file/written_beside_the_vcf.rs);
three derived notes in [to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs)

---

## 1. The gap the step inherited, run before it was built on

Spec §7 is *one writer, three sources — a file the user supplied, the defaults in the binary, or
the fit*. **Two of the three worked.** A run scoring from a supplied file has `Supplied`
calibrations and no rate map, and `ParametersFile::of_run` took the fit's own
`BTreeMap<ReadGroupId, Estimate<ErrorRate>>`. Measured on `a_file_using_every_shape`, read into a
run and written back out:

```
thread '…::probe_supplied_run_writes_a_file' panicked at from_run_parameters.rs:282:13:
read group 0's calibration is Supplied and no rate was offered for it; only a `Defaulted`
calibration can have none …
```

So the source the format exists for — direct mode's user-facing input
([`run_streaming.md`](../../ng/spec/run_streaming.md) §2) — could not write the file it used.

**The minimal fix would have been wrong.** Relaxing the assertion to admit `Supplied` compiles and
passes, and drops every `observations` count on the way out: a file read and written back loses the
evidence behind its multipliers, which is the one thing `RunParametersFromFile` exists to carry back
(*"a run that read a file and then wrote one has to write back what it read"*).

## 2. What replaced it

`ReadsBehindEachCalibration` — a dense `Vec<Option<EvidenceCount>>` with **three constructors, one
a source**, and nothing else public:

| source | constructor | where the count lives |
|---|---|---|
| the fit | `of_the_fits_rates(rates, calibrations)` | on the `Estimate<ErrorRate>` assembly read and did not store |
| a supplied file | `as_a_file_recorded_them(counts)` | in the file's own `observations` |
| the defaults | `nothing_was_fitted(read_group_count)` | nowhere — nothing counted anything |

**The fit's two checks moved with the rates rather than being deleted.** *A rate set covering some
of the run's read groups* and *a rate whose warrant disagrees with its calibration's* are both about
the rates, so they are now `of_the_fits_rates`'s, and their four tests call it directly. What `of_run`
keeps is the one check the three sources share: the counts cover exactly the run's read groups, since
they are joined to the calibration axis by position.

**One behaviour changed, and it is a fidelity gain rather than a wash.** The retired test helper
`the_rates_the_projection_out_reads` wrapped a file's counts in a placeholder rate of 1e-3, turning an
*absent* count into `observations = 0`; under a `fitted_here` warrant that wrote `reads = 0` back into
the file, which says *this fit produced a number from no reads at all*. `as_a_file_recorded_them`
passes the `None` through.

## 3. The census, and the run that has none

`to_run_parameters_for` took `&CensusIdentity`. **Direct mode has none** — no pre-pass, no psp, no
census — and it is the mode this file format is for. The argument is `Option<&CensusIdentity>` now.

**`None` keeps the file's warrants, and spec §2.1 settles it rather than taste.** §2.1 considered
demoting on every read and rejected it *because it breaks the two-mode oracle*: the same cohort called
in direct mode from a file and in psp mode from the fit in memory must report the same warrants for
identical genotypes. Direct mode is exactly the mode with no census, so demoting whenever there is
nothing to compare against **is** demoting on every read under another name.

What that gives up is stated rather than hidden: a file fitted under another census of this cohort
reads into a direct-mode run with its `fitted_here` warrants intact, where the same file in psp mode
is demoted. That is a difference in what a run *reports* and never in what it *computes* — §2 is
explicit that consumers combine warrants and do not branch on them — and the alternative trades it for
a difference in what two modes report about the same call.

**The three refusals are untouched.** A missing census reaches only the fourth binding, the one that
demotes; the reference, the sample list and the read-group table still refuse, and the test asserts
they do with `None` in hand.

## 4. `[fitted_from].census` on a run that has none, and on a run that was demoted

Two cases, and they are the ones Milestone D handed back.

**A run with no census writes an empty list of terms**, through
`CensusIdentity::of_a_run_with_no_census()` rather than by constructing `Vec::new()` at a call site.
Empty rather than an absent section, because `census_disagreement` already treats *a term one identity
has and the other does not* as a disagreement: a **psp-mode** run reading such a file finds one at the
first term and demotes, which is the right answer, and the demotion costs the file nothing, since
`weaker_of` is a no-op on numbers that are already `supplied` or `defaulted`.

**Two things about that had to be corrected after review.** *Two empty term lists agree* — so a run
that itself has no census finds no disagreement, and the file's own prose had claimed otherwise. And
the reader now treats `Some(an identity naming no terms)` exactly as it treats `None`, because a
driver holding one `CensusIdentity` would otherwise hand it the writer's spelling of *this run had no
census* and get every number demoted — the outcome the `None` arm exists to avoid.

**A demoted run writes back the terms it read.** Three answers were possible and two say something
false — this run's own census claims the numbers were fitted under terms they were not; no census at
all claims they were fitted under none. The file that comes out is **stable**: read back by the same
run it is no longer demoted, because its census is now the file's own, and every warrant is already no
stronger than `supplied`, so writing it a second time gives the same file. That is spec §7's *a run is
reproducible from its own output* on the one path where a run's numbers change on the way in.
`a_demoted_run_writes_back_the_census_its_numbers_came_from` asserts all of it.

## 5. What a reader can now see, and the numbers behind it

**⚑ Everything in this section is what the step landed *after* three reviews changed it** — see
[the fixes report](../reviews/ng_parameters_file_f1_fixes_applied_2026-08-31.md). The first draft
said *"fitted from your data"*, which a file demoted under spec §2.1 cannot support, and its
empty-census note promised a demotion that does not happen.

**Measured, on this module's own two fixtures.** Before this step a fitted run's file and a defaults
run's file opened with **the same 39 lines of prose**, and the first thing in either that said which
run it was is the **`warrant` on line 105**, in the first row of `base_quality_calibration`.
Everything above that line was identical but for the cohort's own names.

Three notes now, all **derived from the numbers rather than recorded beside them** — the rule E3 set
for the missing-slippage note, and sharper here, because this file invites its reader to edit a value
and its warrant (§1.2 goal 3) and a recorded count would then be a sentence at the top contradicting
the rows below it.

- **Above `format_version`**: how many of the file's seven groups of numbers were fitted from reads,
  naming them. A defaults run's file opens *"Nothing in this file was fitted from reads — 0 of its 7
  groups of numbers"*; a fitted run's says *"All 7 groups … were fitted from reads"* and names them
  too, because that is the arm making the strongest claim and it gave the reader nothing to count;
  a partly fitted run gets the fraction and the names, which is the commoner case. **Then, on every
  file, which groups can say *whose* reads and which cannot** — see below.
- **`[fitted_from]`**: its first clause is conditional. It headed a defaults run's file with *"What
  these numbers were fitted from"*, which a geneticist read as a claim that they had been. The key name
  is the format and does not move; what it means is said in the note.
- **`[fitted_from.census]`**: an empty term list is explained where it stands.

**Whose reads is a second question, and only three of the seven groups can answer it.** A `warrant`
is the word that answers it, and the base-quality calibration, the inbreeding coefficients and the
repeat-tract substitution rates are the only groups carrying one on every number. A slippage row
carries a smoothing origin, a length spectrum a rung, a contamination row which reads it was fitted
from — each says *how* a number was arrived at and none has a state meaning *somebody handed this
over*. **So spec §2.1's demotion cannot reach four of the seven**, and the file says so rather than
letting a reader take a demoted file's slippage for this run's own fit. That is the same gap
`demoted_to_no_better_than_supplied` already records for `SeedRung::FittedCurve`, and closing it is
the owner's.

**A *group* is one thing a fit either did or did not do**, and there are seven: the base-quality
calibration, contamination, the inbreeding coefficients, the ordinary-site prior's seed, repeat-tract
slippage, repeat-tract length spectra, repeat-tract substitution rates. Not a count of *numbers* —
the substitution-rate table alone grows with `(read group × stratum × ploidy)`, so that would be a
report about the cohort's size rather than about the run. `[sequencing_batches]` and
`[stated_constants]` are excluded from the denominator because no fit can produce them, and including
them would put the full count out of every run's reach.

**The headline counts groups and says so**, since a group counts as fitted where the run measured any
part of it: a cohort in which one plant's coefficient could be fitted and another's could not did fit
the group, and the sentence after the count is what sends a reader to the individual `warrant`.

## 6. Where the file goes

`beside_the_vcf(vcf)` derives the path — the VCF's own name with its compression suffix and then its
format suffix taken off, plus `.parameters.toml`, so `calls.vcf.gz` gives `calls.parameters.toml` —
and `write_beside_the_vcf` writes it there and says where it went. Spec §7 says *beside its VCF* and
nothing more; the stem is used rather than a fixed `parameters.toml` so that two cohorts called into
one directory keep two files rather than one silently overwriting the other.

**What it does not give is a way to tell two files apart by their contents.** The file records no run
date, no caller version and no command line. **That is the command surface's and not this file's** —
what a run stamps into its own output is a property of the invocation, and spec §11 already puts the
neighbouring question (`dump-parameters`) there. Raised at Checkpoint F rather than built here.

**When it is called is still the run driver's**, and so is the one decision left at that call site:
whether a driver whose parameters cannot be projected should log and keep its VCF rather than panic
after the last locus.

## 7. What was left, and why

**`of_run` still requires a `ReferenceDigest`.** `ReferenceDigest::of` refuses a reference read from a
`.fai` alone, which holds no bases to digest, so one run still cannot satisfy the signature. **Left
required and recorded**, because the two missing bindings are not the same kind of thing: a run with no
census cannot check a binding that *demotes*, and §2.1 says what to do then; a run with no reference
digest cannot check one that **refuses**, which is §6's strongest guarantee — *a parameters file fitted
against a different assembly gives a plausible VCF with wrong repeat strata*. Making the digest
optional would weaken that, so it is the run driver's to raise with an owner rather than this writer's
to decide.

**`RunParameterReport` is untouched.** Making it a view over the file is F2.

## 8. Deviations from the plan, absorbed

- **A shared test helper moved.** `unwrapped_comments` existed twice — in `defaults.rs`'s tests and,
  written again, in the new prose tests. It is in the module's shared `tests` fixture now, and
  `defaults.rs`'s copy is gone.
- **`to_run_parameters`'s test module is `pub(super)`**, and two of its helpers
  (`the_files_read_groups`, `the_counts_the_projection_out_reads`) are visible across the module, so
  that the *supplied file* source is built one way rather than three.
- **One sentence of written prose changed** where this step's own artefact showed it false: the
  reference digest's note said *"the reference these numbers were fitted against"*, which is untrue of
  a defaults run's file. It says *"the reference this run ran against"*.
