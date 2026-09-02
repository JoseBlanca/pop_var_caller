# ng — the run's parameters as a file

*Design spec, 2026-08-28. **No code yet — this settles the design.** It fills the entry
[`run_streaming.md`](run_streaming.md) §10 defers when it says "the parameters file's format —
what the user supplies in direct mode and the fit writes in psp mode. Its own spec, **on direct
mode's critical path**".*

*Reads on: [`run_streaming.md`](run_streaming.md) (the two modes and where in each the parameters
arrive); [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) and its siblings (what
each number means and how it is estimated); [`read_likelihoods.md`](read_likelihoods.md) §3.6 and
§4.4 (what contamination and a parameter's warrant do to a score). Read by: the emission step,
which prints what the run used beside what it called.*

---

## 1. What this is

Before ng can score a single read it needs a set of numbers that no locus can supply: how far to
trust each library's own base qualities, how much of each library's DNA came from somebody else,
how inbred each sample is, how variable the population is, and how often a read slips a repeat in
each kind of tract. Today those numbers exist only in memory — `RunParameters`
([`run_parameters.rs:97`](../../../../src/ng/calling/run_parameters.rs)) assembles them from the
pre-pass and calling reads them through a borrowed view. **Nothing writes them down, and nothing
can read them in.**

This document settles the file that does both: a **TOML** file, one per run, holding every number
calling runs on, together with what each number's warrant is and which inputs it was fitted from.

### 1.1 Why a run needs one at all

Three things follow from having it, and the first is a mode:

- **A run can be handed its parameters instead of fitting them.** That is the whole of direct
  mode: the alignment files go straight to a VCF, with no pre-pass and no psp, because the numbers
  the fit would have produced were supplied. Without this file that mode cannot run
  ([`run_streaming.md`](run_streaming.md) §2), which is why it sits on that mode's critical path.
  Running a caller this way is normal practice — freebayes and GATK are usually run against a
  published parameter set whose chemistry does not change week to week.
- **Every run can say what it ran on.** A genotype called at a contamination fraction of 3 in 100
  and one called at zero are different claims about the same reads, and nothing in the genotype
  says which it was ([`read_likelihoods.md`](read_likelihoods.md) §3.6). The file is where the run
  records that, in a form a person can read and a later run can re-use.
- **A cohort can be re-called without re-fitting.** The fit is the expensive stage; a parameters
  file lets a re-run skip it.

### 1.2 Goals

1. **Round-trip exactly.** What the fit produced, written and read back, gives the same genotypes.
   Not "close enough": the two-mode oracle ([`run_streaming.md`](run_streaming.md) §12) compares
   VCFs, and a parameter that survives a write to five decimal places will show up there as a
   changed call.
2. **Every number carries its warrant.** Fitted here, borrowed from a neighbouring grain, supplied,
   or defaulted — the four states `Provenance`
   ([`parameter_estimation/mod.rs:60`](../../../../src/ng/parameter_estimation/mod.rs)) already
   names. A value without its warrant is a number nobody can judge.
3. **A person can read it and change one line.** A user who wants to raise one library's error rate
   should not need a tool, a schema, or this document.
4. **It cannot be silently paired with the wrong inputs.** A parameters file fitted against a
   different reference, or against a different set of samples, must fail by name rather than
   produce a plausible VCF.
5. **Degrade across the committed range.** One sample to several thousand (`CLAUDE.md`). **Two of
   the file's axes grow with the cohort** — one row per sample, and one row per (read group ×
   stratum × ploidy) for the repeat-tract substitution rate — and §9 prices both. The second is
   the larger by two orders of magnitude at 3,000 samples, and §9 says what is and is not settled
   about it.

### 1.3 Non-goals, and what this document does not do

- **It does not define how any number is estimated.** That is step 4's eight documents, starting at
  [`parameter_prepass.md`](parameter_prepass.md). This document only carries their results.
- **It does not replace the census.** The census holds *evidence*; this file holds *the fitted
  numbers that evidence produced*. Where the census lives and how it is rebuilt is
  [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.
- **It does not fix key names or the TOML tree.** It fixes what is stored, in what unit, with what
  warrant, and what a reader must refuse. The exact spelling belongs to the companion architecture
  doc.
- **It does not decide what the VCF header prints.** What an output stage does with the run's
  parameters is the emission step's decision, not this one's.
- **It does not cover per-locus provenance.** Which rung of the tract ladder a repeat tract's
  prior came from, and how many of its scoring cells fell back to a stated constant, are properties
  of that tract and ride on `RepeatTractProvenance`, not on the run.

### 1.4 Vocabulary

- **read group** — one `@RG` line: reads from one library on one lane. The grain the base-error
  and contamination numbers are fitted at, because chemistry is a property of the preparation
  rather than of the plant.
- **slippage group** — the set of read groups whose reads are drawn under one set of repeat-slip
  numbers. Declared by the run, not inferred; `StratumFits::slippage_group_of`
  ([`stratum_fits.rs:370`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)).
- **stratum** — a (motif period, reference repeat count) class of repeat tract. Tomato SL4.00 holds
  462,701 STR loci in 141 of them at the STR path's calling floors
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6).
- **warrant** — what a number entitles a score to claim: fitted, borrowed, supplied, or defaulted.
- **the fit** — psp mode's middle stage, which turns the cohort's census files into these numbers.

---

## 2. The spine: a number, its warrant, and what a write does to it

**Every number in this file is a value plus a warrant plus a count of what was behind it.** That
shape already exists: `Estimate<T>` is `{ value, provenance, observations }`
([`parameter_estimation/mod.rs:122`](../../../../src/ng/parameter_estimation/mod.rs)), and the
warrant has four states, ranked weakest-last:

| warrant | what it means |
|---|---|
| `FittedHere` | estimated from this grain's own data |
| `Borrowed` | too little data at this grain, so the mean of the sample's other read groups was taken |
| `Supplied` | the run was handed this value rather than fitting it |
| `Defaulted` | nothing could be fitted and nothing was supplied, so a stated constant was used |

The ladder is not this document's invention and `Supplied` is not added for it — both are already
in the enum, and `Supplied` already sits *below* `Borrowed`, on the stated grounds that a number
the run was handed says nothing about this data where a borrowed one is at least a measurement of a
neighbouring grain.

**Consumers combine warrants; they do not branch on them**
([`read_likelihoods.md`](read_likelihoods.md) §4.4). A score resting on one fitted and one borrowed
parameter is a borrowed score. So the warrant changes what a run *reports*, never what it
*computes* — which is what makes the next decision safe.

### 2.1 Decision: the file preserves the warrant, and mismatched files are demoted to `Supplied`

A number the fit called `FittedHere` was fitted from *some* cohort's data. Read back into a run
over *that same* cohort it is still a fitted number; read into a run over a different one it is a
number somebody handed over. The file can tell those apart, because §6 binds it to the inputs it
was fitted from.

**So: a parameters file records the warrant the fitting run assigned. A reader keeps that warrant
where the file's binding matches the run's inputs, and demotes every number to `Supplied` where it
does not.** A number a person typed in has no fitted warrant to begin with and is `Supplied`
already.

*The alternative — demote on every read, on the grounds that a file is always something handed
over — was rejected because it breaks the two-mode oracle.* The same cohort called in direct mode
from a file and in psp mode from the fit in memory would then report different warrants for
identical genotypes, and the run-parameter block is exactly the thing an output stage prints. The
oracle would have to be told to ignore a difference that is real.

**Trap for the coder:** demotion is per-file, not per-number. A file whose binding does not match
is demoted wholesale — there is no state in which some of its numbers are fitted and others are
not because of the binding. Per-number differences come only from what the fit itself recorded.

---

## 3. What the file holds

Section by section, and every one of them is a field of `RunParameters`
([`run_parameters.rs:97`](../../../../src/ng/calling/run_parameters.rs)) except the identity block,
which is new here.

### 3.1 Identity and binding

What the numbers were fitted from: the reference's content digest, the ordered list of sample
names, the read-group table (id, `@RG ID`, library name, sample), and the census recording terms
the fit ran under. §6 says what each refuses.

### 3.2 The run's ploidy

One integer. It is a property of the run rather than of the fit, and it is written so that a
supplied file cannot be paired with a run at a different ploidy without saying so.

### 3.3 Per read group — the base-quality calibration

**One multiplier on each read's own reported error probability**, per read group, with its warrant
and the number of observations behind it.

**A read group with no usable rate gets a scale of one, marked `Defaulted` — never a fitted zero.**
A zero scale charges every read of that library the error floor, which is maximal confidence about
every base, drawn from a fit that found no errors at all. `ReadGroupCalibration::from_fitted_rate`
refuses it and the assembly takes `ReadGroupCalibration::defaulted` instead
([`run_parameters.rs`](../../../../src/ng/calling/run_parameters.rs), module documentation). The
file must be able to express that state, which is why the warrant travels with the value rather
than being inferred from it: a scale of exactly 1.0 is a legitimate fitted answer as well as the
default.

### 3.4 Per read group — contamination, and the batching it was drawn against

**Contamination is absent or measured, and never a fitted zero.** Three states, and the file has to
keep them apart:

- **the whole run is uncontaminated** — no read group identified any fraction. The read likelihood
  computes its plain formula, which is the simple case for that model rather than the weak one
  ([`read_likelihoods.md`](read_likelihoods.md) §3.6). Expressed as **absence**: the contamination
  table is not written at all.
- **measured and found clean** — a fraction of zero with non-zero evidence counts.
- **measured and non-zero** — a fraction, with its evidence counts.

Only the counts tell the second from the third; `ContaminationView::was_measured` is the predicate,
and where *some* read group identified a fraction, **every** read group needs an entry. `Option<T>`
is absence and never a sentinel — a missing table means the first state, and a zero never means the
first.

Beside it goes **the sequencing batching**: who was sequenced alongside whom, which is the
population a contaminating read is drawn from. It is written even where no contamination was
fitted, because it is a fact about the run rather than about the fit.

**The grain is the read group, and a row names both its read group and its library.** A
preparation sequenced over several lanes gives several read groups sharing one library name, and
two lanes of one preparation can carry genuinely different fractions, because index hopping happens
on a flowcell and not in a tube. A per-sample row would have to pick one fraction or average two,
and both throw away the distinction the grain exists for
([`run_report.rs:172`](../../../../src/ng/calling/run_report.rs)).

### 3.5 Per sample — the inbreeding coefficient

One number a sample, in the run's sample order. The fit keys them by sample *name*; the file writes
the name beside the value, because the order is the run's and a file that carried only an order
would be silently wrong against a re-ordered sample list. At least one is required.

**This is the file's only cohort-sized axis** — one row a sample, so 3,000 rows at the top of the
committed range (§9).

### 3.6 The ordinary-site prior's seed

Three values: the reference concentration, the total alternative concentration, and which regime
they were derived under (`SpectrumSeed`,
[`genotype_prior/mod.rs:497`](../../../../src/ng/calling/genotype_prior/mod.rs)).

**Written as the seed, not as the moments it came from.** The seed is built once per run by
`RunParameters::seed_from_moments` from two integrals of the fitted frequency density, and what
varies per locus is only how it is spread across that locus's alleles
([`calling_priors.md`](calling_priors.md) §2.3). Writing the moments instead would mean the reader
re-deriving the seed, and any change to that derivation would silently re-interpret every existing
file. The moments may be written *beside* it for a human, marked as informational.

**An alternative total of exactly zero is a real answer** — a fully invariant cohort — and is not
floored on the way in or out; the flooring belongs to the per-locus expansion.

### 3.7 Repeat tracts — slippage, length spectra, and the substitution rate

The largest section, and the one whose axes matter:

- **which slippage group each read group's reads are drawn under** — the run's own declaration.
- **per (stratum × slippage group): the slippage numbers**, with the warrant on each. The axis is
  the *slippage group*, not the read group (`StratumRow`,
  [`stratum_fits.rs:355`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)), which is
  what stops this section growing with the number of libraries.
- **per stratum: the fitted length spectrum and its concentration** — the tract ladder's top rung.
  **Only a stratum fitted on its own tracts has one.** A stratum furnished from its period's
  slippage curves carries no length spectrum at all, by construction, and that absence is what the
  middle rung exists to answer. So absence here is data, not a hole.
- **per motif period: one pooled length spectrum** — the middle rung, present only where the run
  asked for it.
- **the run's stated concentration** — the bottom rung: the median fitted concentration where any
  stratum was fitted, and `STATED_FLAT_CONCENTRATION`
  ([`stratum_fits.rs:342`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs), 1.0)
  where none was.
- **per stratum: the substitution rate** inside the tract.

### 3.8 The constants no fit produces

**Written explicitly, so that a run's inherited numbers are visible and editable rather than buried
in the binary.** The one that exists today is the share of repeat-tract reads that came from
nowhere the model can explain: `DEFAULT_OUTLIER_WEIGHT`
([`likelihood/ssr.rs:82`](../../../../src/ng/calling/likelihood/ssr.rs)), 0.01, inherited from the
existing caller and never measured here.

**Marking it soft is the point of writing it down.** A number that appears in a file the user can
edit is a number the project has admitted is a guess; one that only appears in the source reads as
a decision.

### 3.9 What the run counted as a repeat

**`[repeat_routing]` — the thresholds the run used to decide which stretches of the reference are
repeat tracts and which are ordinary sequence** ([`run_ssr_observations.md`](run_ssr_observations.md)
§2). A repeat catalog is built below every floor a caller routes on, deliberately, so a run picks
its own line inside that gap; two runs over the same reference and the same catalog can therefore
analyse different ground, and nothing else in this file would say so.

**A property of the run, not of the fit — so it is a section of its own and not part of
`[fitted_from]`.** `[fitted_from]` answers *where did these numbers come from*, and every mismatch
in it either refuses the file or demotes it (§6). This answers *what did this run treat as a
repeat*, which is the same kind of fact as the ploidy of §3.2: written so that a supplied file
cannot be paired with a run that routed differently without the difference being on the record.

**Eight values, spelled as the flags that set them.** Five have flags on
`call-from-alignments` — `min_copies` (six integers, one per period 1 to 6), `min_period`,
`max_period`, `max_str_len`, `min_purity`. Three do not, and are written anyway because a record
that omitted them would not say what the run actually asked the catalog for: `min_flank_bp`, which
is pinned at the catalog's own floor because the rows below it were never written; `min_score`,
which gates a scanner's output and so gates nothing in a run that reads a file; and
`bundle_threshold`, the distance within which two tracts are one bundle rather than two loci.

```toml
[repeat_routing]
min_copies = [8, 6, 6, 6, 5, 4]
min_period = 1
max_period = 6
max_str_len = 100
min_purity = 0.8
min_flank_bp = 15
min_score = 0
bundle_threshold = 15
```

**Absent means the file does not say**, which is §5's rule and not a new one: a file written by a
build older than this section, or one a person wrote by hand, records no routing, and a reader must
not read that as *the defaults*. It is a sixth row of §5's table.

#### A difference is reported and never refused

**Decided 2026-09-02 (owner)** — the ruling
[`run_ssr_observations.md`](run_ssr_observations.md) §2.3 records: *"The user could supply any
priors to the caller; it is their decision to use the same routing criteria or not."*

So a run handed a parameters file whose `[repeat_routing]` differs from its own **says so and calls
on**. Nothing refuses, and nothing is demoted to `Supplied` — unlike the census mismatch of §6,
which demotes because the numbers were fitted elsewhere. These numbers were not: they are as
warranted as they ever were, and only the ground they are applied to has moved.

**What the user is taking on, in one sentence**, because the report cannot be read without it: a
tract this run admits but the fit's own selection did not is scored from strata fitted over other
loci, or from the stated defaults — and the per-cell warrants already label which.

---

## 4. Decision: TOML, not JSON

**TOML.** Three reasons, and the first is that the repo already made this choice:

- **Production's psp header is TOML** — a length-prefixed plain-text body, parsed with the `toml`
  crate ([`src/psp/header.rs:10-12,42`](../../../../src/psp/header.rs)) — and ng's psp inherits
  that deliberately: its header "stays plain text so that `head` and a TOML parser can read it, as
  production's does" ([`psp_record_encoding.md`](psp_record_encoding.md) §1.3). A run's two
  user-readable artefacts would otherwise be in two different formats, for no gain.
- **Comments.** Goal 3 is a person changing one line, and what a person needs beside a number is
  where it came from and what moving it costs. In TOML that is a comment. In JSON it has to become
  another field, which then has to be parsed, validated, and kept in step with the value it
  describes — and a stale annotation is worse than none.
- **The dependency already exists.** No second parser, no second failure mode.

**What it costs.** The per-sample axis (§3.5) is an array of tables, which is TOML's most verbose
shape. Write each sample's row as a single inline table on one line, and the numeric rows of §3.7
as arrays of arrays rather than arrays of tables; both stay readable and neither needs a custom
encoder.

*JSON was the alternative considered.* It round-trips floats with the same care and has better
tooling outside this repo. It loses on both of the first two points, and the second is decisive:
this file's whole design puts a warrant beside every value, and a format that cannot annotate makes
every annotation into schema.

**Floats must round-trip exactly** (goal 1). TOML's float is an IEEE 754 double, so the values fit;
whether the crate's serializer emits enough digits to recover each one **has not been checked here**
and is a thing to establish before trusting it — §13's first test is what holds it, not this
sentence. If it does not, the fix is a serializer that formats floats for round-trip and not a
different file format.

---

## 5. Absent, zero, and default are three different claims

The failure this section exists to prevent is a reader that collapses them, and it is the one most
likely to cost a day. In one place:

| what is true | how the file says it | what a reader must not do |
|---|---|---|
| no read group identified any contamination | the contamination table is absent | write zeros for every read group |
| a read group was measured and is clean | fraction 0, evidence counts non-zero | read it as unmeasured |
| a read group's error rate could not be fitted | scale 1.0, warrant `Defaulted` | read the 1.0 as a fitted answer |
| a stratum was furnished from its period's curves | no length spectrum entry | fall back silently to the flat rung without saying so |
| a read group has no reads in a stratum | no row for that (stratum × slippage group) | write a zero slip rate |
| the file does not say what its run counted as a repeat | `[repeat_routing]` is absent (§3.9) | read it as the defaults |

**The rule underneath all five rows:** `Option<T>` is absence, never a sentinel, and a warrant is
carried rather than inferred from the value.

---

## 6. What the file is bound to, and what it refuses

The file names the inputs its numbers were fitted from, and a run reading it compares. This is the
same shape the census already uses one level down — a census file names the psp it was built from
by digesting that psp's header fields, and mismatches fail naming the field that differs
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1) — pointed at a
different object.

Four bindings, and each has a failure attached:

- **the reference's content digest.** A parameters file fitted against a different assembly gives a
  plausible VCF with wrong repeat strata. **Refuse.**
- **the sample list, in order, by name.** The inbreeding coefficients are per sample and the file
  carries names for exactly this reason (§3.5). A file listing samples the run does not have, or
  missing ones it does: **refuse**, naming the samples.
- **the read-group table.** The calibration and contamination axes are dense over `0..n` with
  nothing missing — a gap drops the highest read group entirely and surfaces as a panic at
  whichever locus first carries one of that library's reads
  ([`run_parameters.rs`](../../../../src/ng/calling/run_parameters.rs), module documentation).
  Checked here, at read, so the message is about the run rather than about a locus.
- **the census recording terms the fit ran under.** These identify *which* census produced the
  numbers. They are recorded in the census file rather than in the psp header, because several
  different censuses can be built from one psp
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1) — a psp has no
  single census to name. **A mismatch here is not a refusal**: the numbers are still numbers, and
  the run demotes the whole file to `Supplied` (§2.1) and carries on.

**The first three refuse and the fourth demotes**, and the line between them is whether the file
can still be *interpreted*. A file against another reference cannot: its strata mean something
else. A file fitted from a different census of this same cohort can — it is simply less warranted.

---

## 7. When it is written, and where

**Every run writes the parameters file it used, beside its VCF, whatever the numbers came from.**
One writer, three sources — a file the user supplied, the defaults in the binary, or the fit — and
after assembly the run cannot tell them apart, which is the point.

Three things this buys, and the third is what makes §8 workable:

- **A run is reproducible from its own output.** The file beside a VCF re-runs that VCF.
- **A defaults run is auditable.** The numbers a run guessed are on disk in the same form as
  numbers it fitted.
- **Editing starts from something.** A user who wants to change one number copies the file their
  run just wrote and changes a line, rather than composing one from this document.

**Writing is unconditional.** An earlier proposal made it a flag that could be switched off; there
is no case for that flag — the file is small beside a VCF (§9) and the run that most needs its
parameters recorded is the one whose operator did not think to ask.

---

## 8. Decision: the defaults live in the binary, not in a shipped file

**Owner's decision, 2026-08-28.** A user choosing to run without a fit should not have to find a
file on disk. So the default for every parameter is a named `pub const` in the source with its
origin recorded beside it, in this repo's existing convention, and "run with defaults" is a flag
rather than a path.

*The alternative — ship a defaults file and make it the only source — was rejected on ergonomics:
it makes the simplest run depend on locating an installed artefact.* What that alternative was
protecting is kept by §7 instead: because every run writes out what it used, a defaults run still
produces the file, so the defaults are inspectable and editable without ever being something the
user has to find first.

**Not every parameter has an honest default, and the file must not hide which.** Three cases:

- **has one**: the base-quality calibration scale (1.0 — no calibration), the repeat-tract outlier
  weight (0.01, inherited, §3.8), the flat concentration (1.0). All marked `Defaulted`.
- **absence is the default**: contamination. A run told nothing about contamination is
  uncontaminated, which is a real model state and not a guess (§3.4).
- **has one, and it has to be measured before it exists** — decided 2026-08-28: the
  per-(stratum × slippage group) slippage numbers. **The defaults are fitted from the GIAB HG002
  alignments** and compiled in like the others.

**That last default is the softest number in the file, and it should be read that way.** It is one
individual, one library preparation, one genome, standing in for every chemistry a run might use —
where the others are either "no calibration" (a scale of one) or a real model state (no
contamination). A tomato PCR library taking a human PCR-free slip rate is a guess in a way that a
scale of one is not. Two consequences:

- the default is marked `Defaulted` like the rest, and its origin — which alignments, at what
  depth, on which date — is written into the file as a comment beside it, so a user can see what
  they are inheriting without opening this document;
- **the measurement does not exist yet.** Nothing in this repo has fitted slippage from HG002 for
  this purpose. Until it does, a run without slippage numbers has no defaults to fall back on, and
  the fallback behaviour is whatever `StratumFits` already does with an absent
  `(stratum, slippage group)` row — see §12, question 1.

---

## 9. Cross-cutting concerns

**Size. Corrected 2026-08-28, on rows measured from the built shape** — an earlier version of this
paragraph counted three axes and missed the largest one.

**Four axes, and two of them grow with the cohort.** One row per read group (§3.3, §3.4), one row
per sample (§3.5), one row per (stratum × slippage group) (§3.7) — and **one row per (read group ×
stratum × ploidy)**, which is the grain the repeat-tract substitution rate is fitted and stored at
(`StratumKey`, [`ssr/mod.rs`](../../../../src/ng/parameter_estimation/ssr/mod.rs)). §3.7's phrase
"per stratum: the substitution rate" understates it: the rate is per read group as well, because
how often a base misreads inside a tract is a property of the chemistry.

Tomato SL4.00 has 141 strata at the STR path's calling floors, and a run usually has one slippage
group. Measured on the one-row-a-line inline form: **146 bytes an inbreeding row, 146 a
substitution-rate row.** At the top of the committed range — 3,000 samples, one library each:

| axis | rows | size |
|---|---|---|
| per sample (§3.5) | 3,000 | **0.44 MB** |
| per (read group × stratum × ploidy) (§3.7) | up to 3,000 × 141 | **up to 62 MB**, and 6 MB where one stratum in ten carries a fitted rate |

So the per-sample axis is what the old paragraph said it was — a few hundred kilobytes, negligible
beside the VCF, still openable in an editor. **The substitution rate is not**, and how large it
gets is set by how many (read group × stratum) pairs a cohort actually fits a rate for, which
nothing here has measured. A row exists only where one was fitted.

**What follows for a design, and it is not settled here.** At a few dozen samples the file is a few
megabytes and nothing needs doing. At a thousand it is the largest artefact a run writes that is
not a VCF. Two ways out if it matters — pooling the rate across read groups, which throws away a
chemistry distinction the fit makes, or writing this one table in a form that is not one line a
row — both cost something §1.2's goals ask for, and neither should be chosen before somebody counts
the fitted pairs on a real cohort.

**Memory.** Read once at run start into `RunParameters`, which is "assembled once, before any locus
is called, and never written afterwards"
([`run_parameters.rs:97`](../../../../src/ng/calling/run_parameters.rs)). The file's parsed form is
dropped after assembly; nothing holds the TOML tree during calling.

**Concurrency.** None. The file is read before any worker starts and written after the last locus
is emitted. It is the one part of a run with no parallelism to get wrong.

**Errors.** Every refusal in §6 names the field and the two values that differ, in the shape the
census's own refusal already uses. A malformed file fails at read with a line number, which is what
using an existing parser buys.

---

## 10. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the numbers themselves | `RunParameters` ([`run_parameters.rs:97`](../../../../src/ng/calling/run_parameters.rs)) | the file is its serialised form; the four assembly rules become §3 and §5 unchanged |
| value + warrant + count | `Estimate<T>`, `Provenance` ([`parameter_estimation/mod.rs:60,122`](../../../../src/ng/parameter_estimation/mod.rs)) | written as-is; `Supplied` already exists for this file's sake |
| what a run used, for an output | `RunParameterReport` ([`run_report.rs:54`](../../../../src/ng/calling/run_report.rs)) | becomes a *view over* the parameters file rather than a parallel structure; today it covers contamination, batching and the outlier weight only, and its only callers are tests |
| TOML parse/emit | the `toml` crate, already a dependency of the psp header ([`src/psp/header.rs:42`](../../../../src/psp/header.rs)) | no new dependency |
| binding by digest, and refusing on mismatch | the census file's identity block ([`census_file.rs`](../../../../src/ng/parameter_estimation/joint/census_file.rs)) | the same shape one level up (§6) |

**Parity oracle:** the two-mode comparison already required by
[`run_streaming.md`](run_streaming.md) §12 — the same cohort, called in psp mode from the fit in
memory and in direct mode from the file that fit wrote, must give the same VCF.

---

## 11. Deferred, with a recommended home

- **What the VCF header prints from this file.** The emission step's document. This one settles
  what is available, not what is shown.
- **A command that writes the defaults without running a caller.** Useful and cheap, but it is a
  command-surface question and belongs with the rest of `pop_var_caller_exp`'s subcommands, whose
  names kebab-case from their enum variants
  ([`cli.rs:21`](../../../../src/pop_var_caller_exp/cli.rs)).
- **Merging two parameters files** — for instance taking one cohort's chemistry with another's
  population numbers. No use for it has arisen; if one does, it belongs in this document rather
  than in a tool.
- **A version field's migration policy.** The file carries a version from the start; what a reader
  does with an older one is deferred until there is an older one.

---

## 12. Open questions

1. **Where the default slippage numbers come from — RESOLVED 2026-08-28 (owner): fit them from the
   GIAB HG002 alignments and compile them in.** This replaces the earlier framing, which asked
   whether a run without slippage should refuse to call repeat tracts or fall back to the tract
   ladder's bottom rung; with a measured default there is a third answer, and it is the one that
   lets a defaults run call tracts at all.

   **What is still owed, and it is work rather than a decision:** the measurement. Which strata it
   covers, at what depth, and what happens at a stratum HG002 does not populate — a repeat class
   that is rare in one human genome but common in a plant — are the fit's questions and belong in a
   report beside the other parameter measurements
   ([`../research/`](../research/)). Until that report exists, §8's third bullet has no numbers
   behind it, and a run with no slippage falls through to whatever `StratumFits` does with an
   absent `(stratum, slippage group)` row today
   ([`stratum_fits.rs:355`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)) —
   behaviour that already exists for partially-fitted runs and that nobody has traced. **Trace it
   before writing the reader**, because if it silently produces a score rather than refusing, the
   gap between "no default yet" and "default compiled in" is invisible in the output.
2. **Whether the ordinary-site prior's moments are written beside the seed — OPEN.** §3.6 writes
   the seed because it is what calling reads. The moments are what a human can interpret.
   **Leaning: write both, with the moments marked informational and ignored on read.** The risk is
   a user editing the moments and expecting the seed to follow. Confirm before code.
3. **Whether the read-group table belongs in this file at all — OPEN.** It is not a fitted number;
   it is run identity, and the alignment files already carry it. It is here because §6 needs
   something to check the dense read-group axis against. **Leaning: keep it, as identity rather
   than as parameters**, in its own section so nobody mistakes it for something to edit.

---

## 13. How we know it works

1. **Round-trip.** A `RunParameters` assembled from a real fit, written and read back, is equal to
   the original — every float, every warrant, every count. This is goal 1 and it is the test the
   whole design rests on.
2. **The two-mode oracle.** The same cohort called in psp mode from the fit in memory, and in
   direct mode from the file that fit wrote, gives the same VCF
   ([`run_streaming.md`](run_streaming.md) §12).
3. **Each of §5's five rows is a test.** A file with an absent contamination table gives an
   uncontaminated run; one with a zero fraction and non-zero counts gives a measured-and-clean read
   group; a `Defaulted` scale of 1.0 does not read back as fitted. These are the states a reader
   collapses, so each needs a fixture where collapsing them changes an answer — not merely a
   fixture where they differ.
4. **Each of §6's three refusals fires, and names the field.** A file against another reference,
   one with a sample the run does not have, one with a gap in the read-group ids.
5. **The fourth binding demotes rather than refuses**, and the demotion is visible in the run's
   report: same genotypes, every warrant `Supplied`.
6. **A hand-written minimal file runs.** Someone writes the smallest file that calls a cohort — no
   contamination, defaults where allowed — from this document alone, and it works. That is goal 3,
   and it is the only test here that cannot be automated.
