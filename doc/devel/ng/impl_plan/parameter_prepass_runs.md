# ng — the fit stage: from stored pileups to a parameters file

**Status:** plan, 2026-09-04. **No code yet.** It follows
[`run_driver_psp_mode.md`](run_driver_psp_mode.md), which finished the walk stage, and picks up
the three things that plan handed on: reading a census back, building one *from* a stored
pileup, and the byte-for-byte census-equality check
([`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §7.12).

**This plan turns settled design into build order. It is not a place for new design.** Where it
meets a question the specs do not answer, it says so and stops at a checkpoint rather than
deciding.

---

## 1. What this closes

**No program in this tree has ever produced a fitted parameters file.** A calling run has two
sources for the numbers it scores with, and neither is a fit:

- `--defaults`, the constants compiled into the binary; or
- `--parameters <file>`, a file somebody hands it, which the run checks against its own
  reference and read groups and then writes back out
  ([`calling_run.rs`](../../../../src/pop_var_caller_exp/calling_run.rs) `run_parameters`).

The assembly point exists and is complete —
[`RunParameters::assemble`](../../../../src/ng/calling/run_parameters.rs) takes the nine groups
of numbers a run scores with — but **it is called 36 times in this tree and 35 of them are
inside test modules.** The thirty-sixth is `examples/ng_prepass_handover_footprint.rs` l.321,
which assembles from invented inputs to measure how much memory the result occupies. Nothing
joins a real fit to it.

So a cohort called today is called under the assumption of no base-quality calibration, no
contamination, no inbreeding, and a population whose allele-frequency density is a stated
constant. **What this plan delivers is a run whose numbers came from its own data**, and a file
that says, quantity by quantity, which of them did.

### 1.1 The route

    generate-psps          alignments  ->  <sample>.psp   (+ <sample>.census, see §3.1)
    generate-census        psps        ->  <sample>.census
    estimate-parameters    censuses    ->  cohort.parameters.toml
    call-from-psps         psps + that file  ->  the VCF

Each arrow is a file on disk, so any stage can be re-run without repeating the one before it.
That matters most at the last two: building a census is the expensive half and fitting is the
half that gets re-run while its knobs are chosen.

**Two words for one object, and this document uses both.** A **pileup** is what a walk over the
reads produces — every covered position with the reads seen at it — and a **`.psp`** is that
pileup written to a file. The prose below says *pileup* wherever it means the object and *psp*
only where the name is literal: a command, a flag, a path.

---

## 2. Scope

**In:**

- a second producer for the census — the same `CensusWriter`, driven from a stored pileup's
  decoded records instead of from alignments;
- §7.12's byte-for-byte agreement between the two producers;
- `generate-census`, the command, over a cohort of stored pileups;
- the wall time and peak memory of the two routes to a census, measured against each other;
- reading a cohort of census files back and fitting them — the generic half
  ([`fit_jointly`](../../../../src/ng/parameter_estimation/joint/fit.rs)) and the repeat-tract
  half (`strata_of_kept_loci` → `gather_strata` → `fit_strata`);
- assembling a `RunParameters` from that fit and writing the parameters file;
- `estimate-parameters`, the command;
- the base-quality calibration, which needs a per-read sum the joint fit does not produce
  (§3.4).

**Out:**

- **the inbreeding coefficient**, which stays what the run declares. It is fitted from a
  sample's windowed genome histogram, which is the other pre-pass route, not this one
  ([`fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs) l.308-315). A file this
  plan writes carries the declared value and says it was declared.
- **which knobs become flags.** The census selection's seed and two counts, the read filters
  and the five locus-generator knobs stay compiled-in constants — the owner's ruling of
  2026-09-04, whose stated reason was that nothing could read a census yet. **That reason
  expires at Checkpoint C**, and the question is put again there rather than answered here.
- **how many samples the fit holds at once, and running it inside a memory ceiling** —
  measurements against the object this plan builds
  ([`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §11 q8, q10).
- **psp-mode performance**, which belongs to the psp-mode performance plan
  (`cohort_merge_psp_path.md`, on branch `ng-psp-calling` and not yet merged here).

---

## 3. Decisions taken before the steps

### 3.1 Both routes to a census stay

**The owner's ruling, 2026-09-04.** `generate-psps` already writes `<sample>.census` beside
`<sample>.psp` from one pass over the alignment files, and that shipped. It is not removed when
the second producer lands. The two are **measured against each other** — wall time and peak
memory of one walk that writes both files, against a walk that writes the pileup and a second
pass that reads it back — and the measurement is recorded rather than acted on (B3).

Keeping both is also what makes §7.12 runnable at all: the check is that one sample's census,
built during the walk and built again from the pileup that walk wrote, is identical byte for
byte. **It is the only thing that says the pileup holds everything a census needs**, and it
fails on precisely the fields that do not survive the round trip.

### 3.2 A census stays one file per sample

The format names the pileup it was built from, by a digest of that pileup's header and its
record count. One file for a whole cohort would have nothing to name. Making the census
cohort-wide is a format redesign and is not in this plan.

### 3.3 Two commands rather than one

`generate-census` stops at a file; `estimate-parameters` starts from one. A single command
would re-read every pileup each time a fit is repeated. The cost of the split is one more
name on the command line.

### 3.4 The base-quality calibration is the one group of numbers this route cannot fit yet

`RunParameters::assemble` takes a per-read-group **error rate** and a per-read-group **minted
read-error total** — Σ over reads of `ln P(this read is wrong)`, with the read count it ran
over. It uses them together: where either is missing the read group takes
`ReadGroupCalibration::defaulted`, scale one, and the parameters file says `Defaulted` against
those numbers.

The joint fit produces the rates and not the totals. The totals come from
[`minted_error_by_read_group`](../../../../src/ng/parameter_estimation/generic/calibration.rs),
which sums one locus's complete observations, at generic loci only, reading each observation's
`q_sum` and `num_obs`.

**A stored pileup carries exactly those observations.** `src/ng/psp/record.rs` encodes
`SampleLocusObservations` to bytes and back with nothing dropped, `q_sum` is held in integer
steps rather than as a float, and a locus's kind is in the record body. So the totals can be
accumulated on a pass over the pileups that is already being made — which is what Milestone E
does. **Until it lands, every parameters file this route writes carries a defaulted
base-quality calibration and says so.** That is a smaller claim than the file could make, and
it is the honest one.

---

## 4. Preconditions — what is already built

Confirm each still holds before starting.

- **The census file.** Its byte layout, the directory of sections, reading one section without
  decoding the rest, and the freshness check against the pileup's identity —
  [`census_file.md`](census_file.md) milestones A and B, complete.
- **The fit reads a cohort of them.** `CohortCensusEvidence::new` refuses a cohort whose
  samples recorded under different terms before a section is decoded; `fit_jointly` takes it.
- **The repeat-tract half.** `strata_of_kept_loci`, `gather_strata`, `fit_strata` in
  [`ssr_fit.rs`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs).
- **All of it runs end to end today, from alignments, in one program.**
  `examples/ng_joint_records_walk.rs` builds each sample's census with `CensusWriter`
  (l.1077), fits the generic half (l.422), fits the repeat-tract half (l.771-869), and
  **already refits the same cohort from census files on disk** (l.1156-1159). That program is
  this plan's prototype and its oracle: what the commands must reproduce.
- **The walk-time census producer**, fed at the gatherer's yield point rather than by the
  pileup writer, so it records what the walk saw rather than what was stored.
- **A cohort of pileups is already opened by a command** — `call-from-psps` takes `--psp` once
  per sample or a directory of them, and the ground comes from the files' own headers. The two
  new commands take their cohort the same way and share the refusals.

---

## 5. Principles (how the order was chosen)

- **The prototype first, the command second.** Every step below has a working equivalent in
  `examples/ng_joint_records_walk.rs`. A step that cannot reproduce the example's answer is a
  step that has changed the arithmetic, and that shows up before any command exists to hide it.
- **Isolate a step whose failure is silent.** A census built from a pileup that quietly lacks
  one field produces a plausible fit, not a crash. §7.12 lands immediately after the producer,
  on a fixture carrying a repeat tract, whose per-read length is the field most likely to go
  missing.
- **Types first, then implementation** (project rule).
- **A number a run writes about itself is measured** — the two-route comparison (B3) is run and
  read, not predicted.

---

## 6. The steps

### Milestone A — a census built from a stored pileup

✅ **A1 — the producer.** Drive `CensusWriter` from a pileup's decoded record stream instead of
from alignments: the same writer, the same selection, fed downstream of every filter and cap
the walk applied. One sample.
*Depends:* —. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.1;
[`census_file.md`](census_file.md) C1.

✅ **A2 — the two producers agree byte for byte.** One sample's census built during a walk and
again from the pileup that walk wrote, identical. Run it on a fixture carrying every corner
state, **including a repeat tract**: a tract's per-read length is the field most likely not to
survive the round trip, and a census that lost it fits a stutter model on nothing.
*Depends:* A1. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §7.12;
[`census_file.md`](census_file.md) C2.

> **Checkpoint A: the pileup holds everything a census needs, demonstrated rather than assumed.**
> Pause for review.
>
> **Reached 2026-09-04**
> ([report](../../reports/implementations/ng_fit_stage_a_2026-09-04.md)). One sample's census
> built during its walk and built again from the psp that walk wrote is the same file byte for
> byte, on both samples of the varying cohort — the fixture with a ten-copy `GT` tract and three
> read groups. Four deliberate defects were run against the comparison: skipping tract loci fails
> all three tests, losing one read at every locus fails two, crediting every read to read group 0
> fails one (the sample with a single read group cannot see it), and **changing a read's minted
> error fails none** — because a census records depth codes and allele counts and no per-read
> quality at all. That last one confirms §3.4 from the code rather than from reasoning: the
> minted-error totals cannot come out of a census as the format stands, which is what step E2
> has to settle.

### Milestone B — `generate-census`, and the two routes measured

✅ **B1 — the command.** A cohort of pileups in, one census beside each. `--psp` once per sample
or a directory, as `call-from-psps` takes it; `--output-dir` as `generate-psps` takes it. A
census that cannot be written fails that sample rather than leaving a short one behind, and the
file appears whole or not at all.
*Depends:* A2. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.1;
[`census_file.md`](census_file.md) C3.

✅ **B2 — the report, and what it says about each sample.** Per sample: the pileup read, the
census written, both sizes, and the loci that went into it. A pileup holding no census-selected
locus is named as contributing nothing rather than omitted.
*Depends:* B1. *Source:* the run-level reporting rule
[`run_driver_psp_mode.md`](run_driver_psp_mode.md) Milestone F.

✅ **B3 — the two routes, measured.** Wall time and peak resident memory of **one walk writing
both files** against **a walk writing the pileup, then `generate-census` reading it back**, on
the six tomato accessions over the two 100 kb intervals and on one cohort large enough for the
second pass's read to dominate. Recorded in the milestone's report; it decides nothing on its
own.
*Depends:* B2.

> **Checkpoint B: two producers, one command, and a number for what each route costs.**
> Pause for review.
>
> **Reached 2026-09-05**
> ([report](../../reports/implementations/ng_fit_stage_b_2026-09-05.md)). `generate-census`
> builds each stored psp's census without opening an alignment file, and on the six tomato
> accessions over the two 100 kb intervals all six are byte-identical to the ones the walk wrote.
>
> **Building the census during the walk is the cheaper route: 1.28 s against 1.40 s over the
> work, and 192 MB peak resident against 188 MB**, across three repetitions that moved by 0.02 s
> and 1 MB. **That is on ground where the selection keeps 198,182 of 200,000 bases** — the budget
> is two million positions and this BED is 200 kb — so the census carries a share of the walk here
> that it would not carry on a whole genome, where the same budget keeps about 1 base in 400.
>
> The comparison's first run reported all six censuses different, and the cause was the harness:
> it recorded the route word into each psp's provenance, so the two psps' headers differed by one
> character and each census correctly named a different file. Sixteen bytes, in the digest. The
> script compares the files as well as timing them for exactly this reason.

### Milestone C — the fit, from census files

**Two steps were added on 2026-09-05, before the four this milestone was written with**, because
a cohort of censuses could not be assembled at all. Every census numbers its read groups from
zero, since a walk sees one sample — so on a two-sample cohort both censuses claim read group 0
and `CohortCensusEvidence::new` refuses them, correctly, as libraries that would be fitted as
one. That is psp mode's normal state, not a corner: the advertised way to walk a cohort is one
invocation a sample.

**The owner's ruling, 2026-09-05: the census records who its read groups are** — the `@RG ID` and
the library, per group — so a cohort merges on names alone and the fit's input stays the
censuses. *Rejected: having the fit open the psps too and renumber from their headers, the way
`call-from-psps` does.* It needs no format change, but it makes the fit's real input "the
censuses **and** their psps", which takes away the thing that made a separate census file worth
having — that the expensive half is done once and the fit re-runs off it alone. *Rejected without
being offered: renumbering by the order the censuses are named*, since the identity would then
depend on the argument order, and the parameters file still could not be written — it names a
read group by its `@RG ID`, its library and its sample.

✅ **C1 — the census records who its read groups are.** `DeclaredReadGroup` — the `@RG ID` and the
library — carried per group in walk-local order, written into the census header and read back,
with the format version bumped. Both producers supply it: the walk from its read-group table, the
psp-driven one from the psp's own header. **The two must still agree byte for byte**, which is
what says they are naming the groups the same way.
*Depends:* B1. *Source:* [`parameters_file.md`](../spec/parameters_file.md) §4's read-group table,
which is what these names have to reach.

✅ **C2 — a cohort of censuses merges on those names.** `CohortCensusEvidence::new` stops refusing
an index collision and instead **renumbers**: run-wide identifiers are assigned in (sample order,
the sample's own group order), which is the rule `ReadGroups::of_merged_tables` already uses for
alignment files, and each sample's section keys are relabelled to match. What it refuses instead
is **two samples declaring the same `@RG ID`**, which is the run-wide uniqueness rule and a real
error rather than a normal state.
*Depends:* C1.

✅ **C3 — read a cohort of censuses.** Open them, build `CohortCensusEvidence`, and surface its
refusals as command errors that name the sample and the field that differs — samples that
recorded under different terms, and a census whose named pileup identity does not match the
pileup beside it. **A refusal, never a panic**: this is the door a stale census arrives at.
*Depends:* C2. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.1;
[`fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs) `CohortCensusEvidence::new`.

✅ **C4 — both halves of the fit, giving the prototype's answers.** `fit_jointly` for the
generic half; `strata_of_kept_loci` → `gather_strata` → `fit_strata` for the repeat-tract half.
**The fitted numbers must equal the ones `examples/ng_joint_records_walk.rs` gets on the same
cohort** — this step changes where the evidence is read from, not the arithmetic.
*Depends:* C3.

⚠ **C5 — assemble a `RunParameters`, and write the file. BLOCKED, 2026-09-05.**

> **`assemble` refuses a read group that has a fitted error rate and no minted read-error
> total**, and says why: the two come from one pass over one set of reads, so one without the
> other means they saw different data (`checked_read_group_count_of`). **§3.4 assumed the pair
> would fall back to a defaulted calibration. It does not — it panics.**
>
> Three ways out, and the choice is the owner's. **Bring Milestone E forward**: accumulate the
> minted totals while `generate-census` reads the psps and store them, at which point the
> calibration is fitted and the pair is whole. **Or hand `assemble` a total that saw no reads**
> beside each fitted rate, which satisfies the check by defeating it. **Or supply no rates
> either**, at which point the run has no read-group axis and `assemble` refuses it outright, so
> this route could write no parameters file at all.
>
> Everything else of this step is written and compiles — the read-group table built from what the
> censuses declare, the nine arguments assembled, the file constructed. Its three tests are
> `#[ignore]`d with this reason.

☐ **C5 — assemble a `RunParameters`, and write the file.** The fit's outputs into
`RunParameters::assemble`, then the parameters file. **Each number carries the warrant it has
earned**: fitted where the fit produced it, defaulted where it did not, and the base-quality
calibration defaulted throughout (§3.4). The file names the reference it was fitted against and
the read groups it names, because that is what a calling run checks it by.
*Depends:* C4. *Source:* [`parameters_file.md`](../spec/parameters_file.md) §6, §7.

☐ **C6 — `estimate-parameters`, the command.** Censuses in, one parameters file out. It refuses
to write over a file it was handed, as a calling run does. Written twice from one cohort it
produces the same bytes.
*Depends:* C5.

> **Checkpoint C: a parameters file produced from data, for the first time in this tree.**
> **The knobs question is re-opened here** — the census selection's seed and two counts, the
> read filters and the locus-generator knobs were held as constants because nothing could read
> a census; something can now. Pause for review.

### Milestone D — the four stages, end to end

☐ **D1 — the whole route on the tomato slice.** `generate-psps` → `generate-census` →
`estimate-parameters` → `call-from-psps --parameters`, and the calling run's report names the
file's numbers as fitted rather than defaulted. A script beside the existing oracles.
*Depends:* C6.

☐ **D2 — what the fitted numbers change.** The same cohort called with `--defaults` and with
the fitted file: how many records each writes, and how many genotypes differ. **A measurement,
not a pass/fail** — it is the first look at what the fit is worth, and a run where nothing
moves is as informative as one where much does.
*Depends:* D1.

> **Checkpoint D: the four commands compose. Pause for review.**

### Milestone E — the base-quality calibration

☐ **E1 — accumulate the minted read errors while the census is built.** `generate-census`
already decodes every record of every pileup; `minted_error_by_read_group` takes one locus's
observations. Per read group, Σ `ln ε` and the read count, over complete observations at generic
loci, before the depth cap.
*Depends:* C6. *Source:*
[`calibration.rs`](../../../../src/ng/parameter_estimation/generic/calibration.rs).

☐ **E2 — where the totals live, and it is a design question.** Two numbers per read group have
to reach `estimate-parameters`. Putting them in the census changes the format and obliges the
walk-time producer to record them too, or §7.12 stops holding; putting them in a small file
beside it does not. **Bring both to the owner at the checkpoint rather than choosing here.**
*Depends:* E1.

☐ **E3 — the calibration is fitted, and the file says so.** The rates from the joint fit and
the totals from E2 into `assemble` together, so `ReadGroupCalibration::from_fitted_rate` is
reached and the warrant moves off `Defaulted`. A zero rate still takes the defaulted arm, and
still says so.
*Depends:* E2.

> **Checkpoint E: every group of numbers in the file is fitted or honestly labelled.**

---

## 7. Verification summary

| milestone | proven by |
|---|---|
| A | one sample's census built two ways, byte for byte, on a fixture carrying a repeat tract (§7.12) |
| B1-B2 | command-level fixtures; an unwritable census fails its sample; a pileup holding no selected locus is named |
| B3 | wall time and peak resident memory of both routes, on the tomato slice and on one larger cohort |
| C3 | every refusal provoked and named — mismatched recording terms, a census whose pileup identity has moved |
| C4 | the fitted numbers equal `examples/ng_joint_records_walk.rs`'s on the same cohort |
| C5-C6 | the file read back gives the same `RunParameters`; written twice, byte-identical; warrants match what was fitted |
| D | the four commands compose on the tomato slice, and the calling run reports fitted numbers |
| E | the calibration's warrant moves off `Defaulted`, and the fitted scale reproduces the rate |

---

## 8. Traps

- **A census is not a pileup and its staleness is silent.** The freshness check exists and names
  the pileup by a digest of its header **as written to the file**, not as the writer held it —
  the walk stage found that the pileup writer amends its header before writing it, so a census
  built from the in-hand copy names a file that does not exist and every check answers *rebuild*
  for ever. Any new producer takes the digest from the file.
- **A read group's identifier is unique across the whole run**, not merely within a sample.
  Refused when alignment files are opened and again when stored files are; a fixture that
  reuses one is refused before anything is measured.
- **The parameters file identifies samples by name and read groups by the sample and the
  identifier together**, never by position. A command that pairs a census to a sample by the
  order its flags were typed will pass every fixture and be wrong on a real cohort.
- **`--threads` builds rayon's global pool, which a process may build once.** A test sweeping
  thread counts runs every later count at the first one's width while reporting a sweep. Thread
  invariance is a script, as it is for the walk stage.

---

## 9. Out of scope (next plans)

- **The inbreeding coefficient on this route** — §2.
- **The fit's memory ceiling and how many samples it holds at once** —
  [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §11 q8, q10.
- **Which walk knobs become flags** — re-opened at Checkpoint C, decided by the owner.
- **psp-mode performance** — `cohort_merge_psp_path.md`, on branch `ng-psp-calling`.

### 9.1 Open question: does skipping the loci nobody varied at speed the fit up too?

**Raised by the owner, 2026-09-04.** The psp-mode performance work skips loci where no sample
carries a difference from the reference, and the question is whether the fit can skip them the
same way.

**What is settled is what the census needs**, and it does not point that way. A census records a
**depth code for every kept position and every read group, including the positions where nothing
varies** — the zero is the denominator, and a spectrum built without it has a shape and no scale
(`parameter_prepass_census_sites.md` §2). Building that depth means walking the record's
observations one read at a time: `CensusWriter::add_generic` adds each read's `num_obs` into
every position its witness covers, per read group. The record head carries one count for the
whole record and no read-group split, so **a body skipped on its head alone leaves the census
with no depth at that position at all** — which is the state that reads as *never walked*, and
is exactly the confusion `mark_walked` exists to prevent.

So the two stages want different things from the same skip: calling can drop an invariant locus
because it emits no record; the fit needs it because it is the denominator. **What might still
transfer is the cheaper decode rather than the skip** — reading a body far enough to total its
depths without materialising its observations. That is a measurement to make once the fit runs
end to end, not a design decision to take now.
