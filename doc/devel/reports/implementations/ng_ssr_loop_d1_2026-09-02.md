# D1 — the instrument for the tract QUAL experiment

**Date:** 2026-09-02. **Plan:** [`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md)
Milestone D step D1. **Design:** [`../../ng/spec/calling_loop_ssr.md`](../../ng/spec/calling_loop_ssr.md) §3.3.
**No change to the caller** — everything here is measurement machinery.

---

## What was built, and the two things it already settles about the runs

Three pieces:

- **[`benchmarks/lib/tract_qual_experiment.py`](../../../../benchmarks/lib/tract_qual_experiment.py)** —
  the scorer. Given a truth VCF, a caller's VCF, a confident-region BED and a BED of repeat
  tracts, it writes two tables: records binned by QUAL against the share that really are at a
  variant tract, and precision and recall as a QUAL threshold sweeps. Both split by motif
  period, homopolymer against period 2 and above.
- **[`examples/ng_tract_simulator.rs`](../../../../examples/ng_tract_simulator.rs)** — repeat
  tracts whose genotypes we chose, sequenced under a slippage we set: a reference, a BAM a
  sample, a truth VCF, the tract ground as a BED, and the stutter model as parameters-file rows.
- **[`benchmarks/lib/run_tract_qual_experiment.sh`](../../../../benchmarks/lib/run_tract_qual_experiment.sh)**
  and **[`benchmarks/ssr_hg002/src/run_ng_coverages.sh`](../../../../benchmarks/ssr_hg002/src/run_ng_coverages.sh)** —
  the drivers: build the tract ground, run the arms, score them.

**Two findings shape what D2 can report, and both were measured while building this.**

**The benchmark the plan names cannot answer the calibration question.** On the GIAB trio's
per-sample benchmark — 100 confident intervals a sample — ng writes **149 repeat-tract records
pooled over the three samples at 30×** (62, 41 and 46) and 150 at 50×. HG002's share of that
ground is 335 repeat tracts over 4,167 bases. A calibration asks how often a record at QUAL 30
is wrong, which is a claim about one in a thousand; 149 records cannot see it. So the experiment
also runs on **GIAB's HG002 tandem-repeat benchmark** (`benchmarks/ssr_hg002`): 50,000 Tier
intervals over 6,089,411 bases with 36,497 assembly-based truth records, at the same 30× and
50×. ng types 20,204 repeat tracts on it and writes **6,351 tract records at 30× and 6,408 at
50×** — 43 times the per-sample ground. Both grounds are kept: the small one for continuity with
the numbers already recorded, the large one for the calibration.

**The fitted-against-defaulted split has one side on any benchmark, and two only on the
simulator.** §3.3 asks for the cells split by where the model's numbers came from. No command
fits a parameters file (§3.4 — deferred), and every GIAB run's own report says so in the same
words: *"numbers behind the calls: 0 of 7 groups the file says were fitted … taken from
constants: … repeat-tract slippage, repeat-tract length spectra, repeat-tract substitution
rates"*. So on GIAB every tract call is `Defaulted` and the split is one cell. On the simulator
the truth is ours, so the simulator writes the stutter model it drew the reads under as the rows
a parameters file states slippage in; a run handed that file reports **"1 of 7 groups the file
says were fitted"**, and the same reads scored under the shipped default are the other arm.

## The two scoring rules, and the two defects the first version had

**Calibration is at the tract; the sweep is at the allele.** QUAL claims the samples here are
not all homozygous reference *at this tract*, so a record counts as truly variant when any truth
record falls inside the tract it sits at. The sweep is `score_ng_recall.sh`'s rule unchanged —
contig, position, REF and ALT equal after both sides are left-aligned and split — so its numbers
are read against the ones already recorded for these callsets.

Both defects were found by reading the records the first version called false at high QUAL, and
both are fixed:

- **The region masks were applied before left-alignment**, and the two orders do not treat the
  two sides alike at a confident interval's first base. At `chr1:69,233,430` the truth set
  carries `TATAATAATA → T`; the Tier interval starts at 69,233,431, so the truth record was
  dropped and ng's identical call — still at 69,233,431 before normalisation, and so still
  inside — was kept and scored a false positive at QUAL 922. Left-aligning first makes the two
  agree: at that site both are now dropped, which is right, because the event's own position is
  outside the confident region and neither side should speak for it.
- **A record was scored against its own span rather than its tract's.** At the homopolymer
  `chr1:52,776,219–52,776,230` ng writes an insertion at 52,776,219 and the truth set describes
  the same length change at 52,776,229, ten bases away at the other end of a compound repeat.
  Scored on the record's own span, padded by one, that read as a false positive.

**Together they moved 105 records into the wrong column.** On the tandem-repeat ground at 30×,
the first version called 105 records above QUAL 200 false where no truth variant was found; the
scorer as it stands calls **7**.

**One column was added for something the simulator made obvious.** A caller may list an
alternative allele and then give no sample a copy of it; after `bcftools norm -m -any` that
allele is its own record and the scorer's rule makes it a false positive. On the simulator at
30× over three samples, **302 of the 303 false homopolymer alleles at QUAL 0 and above have
`AC=0`** — nobody was called with them. That is a different failure from a sample genotyped
wrong, so the sweep carries `fp_with_no_called_copy` beside `fp`.

## The simulator, and the oracle that says its fixture is the one it claims

The forward model is the one §3.5 and the read-likelihood spec already describe: a diploid
genotype a sample a tract, then each read draws an allele, slips with probability
`slip_share` (shorter with `shorter_share`, by one repeat with `one_step_share`), and mis-reads
each base with `substitution_rate`. The defaults are the shipped stutter model's own numbers, so
a run at the defaults is the case where the caller's model is exactly right.

**The fixture checks itself, and the check found a real defect.** The driver types the
simulated reference with `ng_typed_region_dump` and stops the run unless the `ssr_locus` regions
are *exactly* the tracts the simulator laid down — one more is a repeat the flanks grew by
accident, one fewer is a tract the routing swallowed. On the first 4,000-tract fixture **5
tracts came back typed as ordinary sequence**: the flank ended in a near-copy of the motif
(`…CCAAGGT` before a tract of `AAGGC`), which the catalog reads as one more, impure copy at 80%
purity, so the tract it found was not the tract the truth file claimed. The flanks now carry a
written boundary — `period` bases each differing from the motif base opposite — and the check
passes on 4,000 tracts.

Beside that, three end-to-end agreements on a 60-tract fixture at 30× over two samples: the
typed regions reproduce the injected tracts exactly; the run reports 60 tracts built and 60
called; and it writes 43 records, against the 43 tracts the truth file says at least one sample
varies at.

**And the instrument is sensitive to the axis it was built for.** On a 3,000-tract fixture at
10× with one sample and a true slippage of 0.30 — three times what the shipped model assumes —
the defaulted arm writes 1,978 tract records and the arm handed the true slippage writes 1,644,
and the low-quality bins separate: of the period-2-and-above records below QUAL 1, the defaulted
arm has 79 of 428 at a truly variant tract and the fitted arm 151 of 305.

## What this step does not do

**No number here is the experiment's answer.** The figures above are the instrument's own
validation — ground sizes, a defect count before and after, a sensitivity check. The runs and
the report are D2, and the decision they feed belongs to `calling_quality_ssr.md`.

**Arm C is wired but not yet run.** The driver scores the production caller on both GIAB
grounds — the `high-recall` preset on the per-sample benchmark, `ssr-call` on the tandem-repeat
one — and neither has been scored yet.

## Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- The simulator's own oracle, run on every fixture the driver builds: the typed `ssr_locus`
  regions equal the injected tracts, or the run stops.
- The scorer runs on all three grounds: the per-sample benchmark (66 records on HG002's tract
  ground at 30×), the tandem-repeat benchmark (6,797 records on tract ground, 15,175 off it) and
  the simulator (3,528 on, 10 off, at 4,000 tracts and three samples).
