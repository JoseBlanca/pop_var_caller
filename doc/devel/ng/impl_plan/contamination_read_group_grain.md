# ng — contamination is a property of a library: moving the estimate to the read group

**Status:** plan, 2026-08-19. Nothing built. **Written to be run in its own conversation**, so it
states its own context rather than assuming the reader was in the one that raised it.

**The grain is ruled and is not what this plan investigates** (owner, 2026-08-19). What gets
contaminated is a **library**, not an individual: a second sample's DNA enters at library
preparation or at sequencing, so two libraries made from one plant can differ, one carrying stray
reads and the other none. **The contaminating source is a sample; the contaminated thing is a read
group.** ng's estimator currently produces one fraction per sample, which cannot express that.

**What this plan investigates is how to fit it at that grain without the estimate falling apart**,
and then changes the estimator. The reality is fixed; only the method is open.

---

## 1. What this decides

1. **What splits and what does not** when the fraction moves to the read group. §6 argues that only
   one parameter splits and that the expensive half of the model stays where it is — check that
   before designing anything, because it decides how bad the data problem is.
2. **How much data one fraction needs**, and what the estimator does for a read group that has less.
3. **How it is measured at all**, given that no benchmark cohort here has a multi-library sample
   (§7). This is the hard part and it is why the plan leads with it.

**The output** is a change to
[`joint/contamination.rs`](../../../../src/ng/parameter_estimation/joint/contamination.rs), its
specification [`../spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md)
§3.1 and §3.4, and a report carrying the numbers. **The owner decides on the numbers**; no threshold
is fixed here.

---

## 2. Why, and what is already agreed

**The primary specification has said read group all along.**
[`../spec/parameter_prepass.md`](../spec/parameter_prepass.md) §1's table of everything the
parameter fit produces lists **contamination at read-group grain**, under the principle its §1.1
states: *"Noise is chemistry… their unit is the read group."* The implementation and
[`../spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §3.1 say one per
sample. **This is a divergence to repair, not a new direction.**

**It matters on real data at about one sample in ten.** From a 2,085-file, 68-project tomato archive
survey, **157 of 1,707 samples carry more than one library** — 133 carry two, 20 carry three, and
four carry 7, 16, 16 and 42 ([`../spec/read_groups.md`](../spec/read_groups.md) §1, quoted at
[`../spec/parameter_prepass.md`](../spec/parameter_prepass.md) §5). At those samples a per-sample
fraction is an average of libraries that need not resemble each other.

**Nobody else does this, which is worth knowing before assuming it is easy.** ng's estimator, the
production caller, and `verifyBamID2` all put the fraction on the sample. Production splits only the
*contaminant's allele composition* to a sequencing batch (`q_b` per batch, `c_s` per sample,
[`contamination_estimation.rs`](../../../../src/var_calling/contamination_estimation.rs)). So there is
no reference implementation to port and the literature will not answer §1's second question.

---

## 3. Scope

**In.** The grain of the contamination fraction; what the fit has to change to produce it; what a
thin read group gets; how the change is measured; and what the emitted parameters must carry so a
consumer can tell a fitted fraction from a substituted one.

**Out.**

- **Whether to model contamination at all**, and whether it is on by default. Settled elsewhere
  ([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §3.6).
- **The read likelihood's use of the number.** It consumes a fraction per read group already — the
  evidence keys on read group and §2.3 there makes that a requirement — so this change *removes* an
  approximation from that document rather than asking anything of it.
- **The contaminant's allele frequencies.** A separate half of the model (§6), and this plan moves
  only the fraction.
- **Contamination at repeat tracts** ([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md)
  §4.5.1), which consumes whatever grain this produces.

---

## 4. Principles

1. **The grain is not the experiment.** No measurement can make a per-sample fraction correct for a
   sample whose two libraries differ. Measurements here decide *how* to fit it and what to do when
   the data is thin — never whether to.
2. **Do not reach for the panel-partitioning result as an objection without checking it applies**
   (§6). It is the obvious thing to cite and §6 argues it is about a different split.
3. **A change that cannot be measured on the data we hold must say so in its first paragraph**, and
   §7 says what to do instead. Reporting a benchmark result that could not have moved is worse than
   reporting nothing.
4. **A substituted fraction must never be indistinguishable from a fitted one** (§9).

---

## 5. What is already known, so nobody re-measures it

From [`joint/contamination.rs`](../../../../src/ng/parameter_estimation/joint/contamination.rs)'s module
documentation and the two reports behind it
([`../reports/joint_contamination_2026-08-12.md`](../reports/joint_contamination_2026-08-12.md), and a
2026-08-13 measurement it cites):

- **Structure hides contamination rather than inventing it.** On a panel of four subpopulations at
  `F_st` 0.20, a pooled allele frequency returns **0.005 for a sample truly contaminated at 3%**, and
  **exactly zero at 1%**. Both pass any threshold as clean.
- **Partitioning the panel to fix that is worse than ignoring it.** Estimating frequencies within
  groups of twelve adds about **0.015 to every sample's estimate** and puts **41 to 47 of 50 clean
  samples over a 1% threshold**. The fix used instead is a frequency per sample predicted from its
  ancestry coordinates, borrowing the *slope* along each axis and never a neighbour's allele counts.
- **A sample sitting alone at the end of an axis is refused**, because its fitted frequency is mostly
  its own echo — on a panel of 40, 5, 3 and 2 accessions with nobody contaminated, the group of two
  returned a spurious **0.031**. How much of its own frequency a sample supplies depends **only on
  its coordinates**, tracking the damage at 0.027, 0.307, 0.429 and 0.857 across those four groups.
- **Leaving a sample out of the frequency it is judged against changes nothing at 30 reads a
  position** on panels of 40 and of 12, **and at three reads it lifts the worst clean sample above
  the contaminated one**. Built, and off.
- **The estimator is built, wired into the joint fit, and off**, with an `estimate-contamination`
  subcommand exposing it
  ([`estimate_contamination.rs`](../../../../src/pop_var_caller_exp/estimate_contamination.rs)).

---

## 6. The structural question, and why it decides how hard this is

**The model has two halves and the plan's first job is to establish that only one of them splits.**

- **The frequency half** predicts, for each sample, its own expected allele frequency at each
  position, as a linear function of that sample's coordinates in a principal-component space, fitted
  across the whole panel. **Ancestry is a property of the individual**, and two libraries of one plant
  are the same individual — so this half does not split, and the panel it is fitted from does not
  shrink.
- **The fraction half** is one number saying what share of reads came from elsewhere. **That is the
  half that becomes per read group**, each fitted from that read group's own reads.

**If that decomposition holds, §5's 0.015 hazard does not apply to this change**, because that
measurement is about partitioning the *panel* to estimate *frequencies* from twelve samples — and
nothing here partitions the panel or the frequency model. What shrinks is the read count behind one
fraction, which is a different and smaller problem. **Establish this first: it is the difference
between a contained change and a dangerous one, and it is the argument most likely to be waved at
this plan without being checked.**

**Two things do follow and both need answering.**

- **The refusal criterion is per sample and should stay there.** It depends only on ancestry
  coordinates, so a sample whose frequency is its own echo is unfittable at every grain. Confirm it
  transfers unchanged rather than assuming it.
- **Read length is chemistry too.** Two libraries of one sample may have been sequenced at different
  read lengths, so they do not cover the same positions
  ([`../spec/parameter_prepass.md`](../spec/parameter_prepass.md) §5 records the same trap costing
  the stutter fit a spurious difference). Check whether the position set each fraction is fitted over
  has to be held common, and what it costs if it is.

---

## 7. The data problem, stated before any step depends on it

**No cohort in `benchmarks/` can tell the two grains apart.** Every sample in tomato1, tomato2 and
the GIAB trio carries **one library**
([`../spec/parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md), the same archive
survey: 1,550 of 1,707 samples carry one, and every benchmark sample does). At one library per
sample, per-read-group and per-sample are the same number, and a benchmark run would show **no
difference by construction**.

**So the measurement comes from two places and neither is the usual bench.**

- **Simulation, which is where the method is decided.** The harness behind §5's numbers already
  builds a panel with planted contamination
  ([`ng_joint_contamination_harness.rs`](../../../../examples/ng_joint_contamination_harness.rs)).
  Extending it to give a sample several libraries with **different** planted fractions is the
  instrument this plan needs, and it gives exact truth.
- **The tomato archive, which is where the reality is.** 157 samples with more than one library, four
  of them with 7 to 42. **There is no truth there**, so it cannot say whether an estimate is right —
  but it can say whether libraries of one sample come out **different**, which is the claim the whole
  change rests on and which nothing has yet checked.

---

## 8. The steps

### A — Confirm the structure before designing

**A1. Establish what splits.** Read the estimator and write down, parameter by parameter, which are
per sample, which per position, which per panel. **Produce the decomposition of §6 as a table with
the code behind each row**, or refute it. **Gate:** if the frequency half turns out to depend on the
grain of the fraction, this plan's cost estimate is wrong and the owner should hear that before
anything is built.

**A2. Confirm the refusal criterion transfers**, and say what a refusal means when one library of a
sample is refused and another is not.

### B — Build the instrument

**B1. Give the harness multi-library samples**, each library with its own planted fraction, its own
read count, and optionally its own read length. **Gate:** it must reproduce §5's existing numbers
when every sample has one library, or the extension has changed something it should not have.

**B2. Ask the archive whether libraries of one sample differ.** Run the existing per-sample estimator
separately over each library of the 157 multi-library samples and report the spread within samples
against the spread between them. **This is descriptive and it is the single most informative cheap
thing in this plan**: if libraries of one sample are indistinguishable in this archive, the change is
still correct and its urgency is different.

### C — Fit at the new grain

**C1. Split the fraction, share the frequency model**, per §6, behind a switch so both grains run in
one build.

**C2. Find the floor.** Sweep reads per library and positions per library on the simulated panel;
report the fraction's error against both. **What is wanted is the shape of the curve, not a
threshold** — where the estimate stops tracking truth, and how gracefully.

**C3. Decide what a thin read group gets.** Candidates, and the plan should measure rather than
choose: fit it anyway and report the evidence behind it; substitute the sample's other libraries'
pooled fraction; substitute the sample-wide fraction; or refuse it. **Note that refusing is not free**
— a refused library gets `c = 0`, which is a claim that it is clean, so the refusal must be visible
in the output rather than silently equivalent to a measurement of zero.

**C4. Check the read-length confound** (§6) on a simulated pair of libraries differing only in read
length, with equal planted fractions. Any difference the estimator reports there is an artefact.

### D — Land it

**D1. Change the estimator and the subcommand's output** to the new grain.

**D2. Amend [`../spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §3.1
and §3.4**, which say per sample, and its architecture sibling.

**D3. Tell [`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §3.6** that the approximation
it currently records — a per-sample estimate applied to every read group of that sample — is gone.

**D4. Provenance** (§9).

**D5. Re-run the existing benchmarks to show nothing moved.** Every benchmark sample has one library,
so **the correct result is no change at all**; anything else is a defect in the split.

---

## 9. Provenance — what must not be lost

A fraction emitted per read group can now be four different things, and a consumer cannot act on them
alike:

- **fitted from that read group's own reads**, with how many reads and positions stood behind it;
- **substituted** from the sample's other libraries or from a sample-wide fit, saying which;
- **refused** — not identified, which is not the same as zero and must not be encoded as zero without
  a flag beside it;
- **inherited** because the sample has one library, where the two grains coincide and no claim is
  being made about libraries at all.

---

## 10. What could make this exercise worthless

- **Measuring it on the benchmark cohorts.** They are single-library; the result is "no change" by
  arithmetic and says nothing (§7). It is also the right *regression* check — D5 — and reporting it
  as evidence of anything else would be wrong.
- **Citing the 0.015 panel-partitioning result as an objection without checking §6.** It is about
  splitting the panel to fit frequencies, not about splitting one sample's reads.
- **Encoding a refusal as a zero.** A library nobody could measure and a library measured as clean
  are different claims, and the caller would treat them alike.
- **Letting the frequency model follow the fraction to the read group.** That *would* partition the
  panel, and §5 has the measurement of what that costs.
- **Stopping at the estimator.** The number exists to change genotypes; the last check is what the
  caller does with it at a sample whose libraries differ.
