# Is ng's quality score worth believing at a repeat tract?

**Date:** 2026-09-02. **Asked by** [`calling_loop_ssr.md`](../ng/spec/calling_loop_ssr.md) §3.3,
which forbids designing a tract-specific quality correction before the failure has been
observed. **Run by** [`calling_loop_ssr.md`](../ng/impl_plan/calling_loop_ssr.md) Milestone D.
**This report takes no decision** — that belongs to `calling_quality_ssr.md`, written from it.

---

## 1. The answer, in one paragraph

**Ship the site quality as it stands at repeat tracts.** The failure §3.3 warned about — a
quality score that *grows confidently wrong* at a tract because the stutter model under-prices
slip products — **does not appear.** On GIAB's tandem-repeat benchmark at 30×, a tract record
written above QUAL 200 sits where the truth set carries no variant **5 times in 2,882**; on the
ordinary sequence of the same benchmark, scored by the same instrument against the same truth,
ng's corrected SNP and indel quality is wrong **19 times in 10,272**. Those are one in 580 and
one in 540 — the same number. In the middle of the range the tract quality is the *better* of
the two: between QUAL 30 and 50 a tract record is wrong 9 times in 429 against ordinary
sequence's 71 in 483. §3.3's decision rule asks whether the tract quality reaches the standard
the corrected SNP quality reaches on the same benchmark. It does.

**Two things it is not, and both matter to whoever writes the design.** It is not Phred: QUAL
200 claims to be wrong once in 10²⁰ and is wrong about once in 500 — but so is the SNP path's,
so that is a property of the fold and not of tracts. And **it is not gateable**: on ordinary
sequence, raising the threshold from 0 to 200 buys 14 precision points for 17 recall points; at
a homopolymer tract it buys **2 precision points for 28 recall points**, and then **raising it
further takes the 2 points back** — precision peaks at 0.850 at QUAL 50 and falls to 0.831 by
QUAL 200, within a point of the 0.825 it started at, having shed more than half the recall. At
period 2 and above the gate never buys more than 0.2 points at all.

---

## 2. What was run

| ground | what it is | arms | depths |
|---|---|---|---|
| **tandem-repeat benchmark** | GIAB's HG002 tandem-repeat set: 50,000 Tier intervals over 6.09 Mb, 36,497 assembly-based truth records. ng types 20,204 repeat tracts in it and writes 6,351 tract records at 30× | ng; the existing repeat-tract caller (`ssr-call`) | 30×, 50× |
| **per-sample benchmark** | The GIAB trio's 100 confident intervals a sample — where every standing ng number was measured | ng; the existing caller, `high-recall` preset | 30×, 50× |
| **simulator** | Tracts whose genotypes we chose, sequenced under a slippage we set. 4,000 tracts, one sample | ng at `--defaults`; ng handed the slippage the reads were drawn under | 10×, 30×, 100× × slippage 0.02, 0.10, 0.25 |

Instrument: [`benchmarks/lib/tract_qual_experiment.py`](../../../benchmarks/lib/tract_qual_experiment.py),
built as step D1 ([report](implementations/ng_ssr_loop_d1_2026-09-02.md)); driver
[`benchmarks/lib/run_tract_qual_experiment.sh`](../../../benchmarks/lib/run_tract_qual_experiment.sh).
**Calibration is at the tract** — a record counts as right when the truth set carries a variant
in the tract it sits at. **The threshold sweep is at the allele**, on
`score_ng_recall.sh`'s rule: contig, position, REF and ALT equal after both sides are
left-aligned and split.

**Why a second GIAB benchmark.** The per-sample benchmark holds about 4,200 bases of repeat
tract a sample and ng writes **149 tract records on it pooled over the three samples at 30×**.
Every one of them is at a tract the truth set does carry a variant in, at every QUAL — so its
calibration table is a column of zeroes and says nothing. A claim about one error in a thousand
needs more records than that. The tandem-repeat benchmark has 43 times as many.

---

## 3. Calibration: the shape, and it is not the shape the risk predicted

ng's tract records on the tandem-repeat benchmark, binned by QUAL. **Observed error** is the
share of the bin's records sitting at a tract with no truth variant; **claimed** is what QUAL
says that share should be.

| QUAL | 30×, records | observed error | claimed | 50×, records | observed error |
|---|---:|---:|---:|---:|---:|
| 0–1 | 163 | 0.258 | 0.95 | 58 | 0.276 |
| 1–3 | 69 | 0.087 | 0.67 | 30 | 0.167 |
| 3–10 | 127 | 0.031 | 0.26 | 50 | 0.200 |
| 10–20 | 201 | 0.020 | 0.039 | 64 | 0.078 |
| 20–30 | 196 | 0.036 | 0.0040 | 61 | 0.049 |
| 30–50 | 429 | 0.021 | 0.00021 | 168 | 0.036 |
| 50–100 | 923 | 0.017 | ≤ 10⁻⁵ | 688 | 0.022 |
| 100–200 | 1 253 | 0.0032 | ≤ 10⁻¹⁰ | 1 161 | 0.012 |
| **200+** | **2 882** | **0.0017** | ≤ 10⁻²⁰ | **4 017** | **0.0022** |

Three readings, in the order they matter.

**The curve is far too flat, in both directions.** Below QUAL 10 the score is *too timid* — at
QUAL 0–1 it claims to be wrong 95 times in 100 and is wrong 26. Above QUAL 20 it is
*overconfident* — at QUAL 20–30 it claims one error in 250 and delivers one in 28. Above QUAL
100 the observed error stops falling at about one in 500 whatever QUAL says.

**But the SNP path does exactly the same thing, on the same benchmark's ordinary sequence.**
Scored by the same instrument against the same truth set, ng's corrected SNP and indel quality
on the Tier intervals' non-repeat sequence:

| QUAL | 30×, records | observed error | tract ground, observed error |
|---|---:|---:|---:|
| 0–1 | 1 574 | 0.440 | 0.258 |
| 10–20 | 348 | 0.279 | 0.020 |
| 20–30 | 282 | 0.231 | 0.036 |
| 30–50 | 483 | 0.147 | 0.021 |
| 50–100 | 865 | 0.057 | 0.017 |
| 100–200 | 1 442 | 0.025 | 0.0032 |
| **200+** | **10 272** | **0.0018** | **0.0017** |

**At the top the two are the same number and everywhere else the tract path is better.** That
is the opposite of what §3.3 feared, and it is the answer to the decision rule.

**Depth does not make it worse in the way the risk describes.** From 30× to 50× the top bin
moves from 0.0017 to 0.0022 on tract ground and from 0.0018 to 0.0030 on ordinary sequence —
both slightly up, the generic path more so. Nothing here is a tract-specific inflation with
depth.

**What the nine remaining high-quality errors are.** Of the 4,135 tract records above QUAL 100
at 30×, **9 sit at a tract the truth set carries nothing in** — 3 at homopolymers, 3 at period
2, 3 at period 4. They are not stutter noise: they are whole length changes several
repeats long — eight copies of `AC` inserted at `chr5:14,207,463`, four copies of `TATT` at
`chr18:71,537,338`. A mispriced slip product would look like a one-repeat difference; these do
not.

---

## 4. Gateability: the tract quality is not a gate, and the SNP quality is

Precision and recall as the threshold sweeps, tandem-repeat benchmark, 30×, ng:

| ground | QUAL ≥ 0 | QUAL ≥ 30 | QUAL ≥ 100 | QUAL ≥ 200 |
|---|---|---|---|---|
| tract, homopolymer | 0.825 / 0.812 | 0.849 / 0.718 | 0.848 / 0.536 | 0.831 / 0.360 |
| tract, period 2+ | 0.835 / 0.690 | 0.841 / 0.645 | 0.836 / 0.549 | 0.834 / 0.427 |
| ordinary sequence | 0.759 / 0.582 | 0.844 / 0.551 | 0.877 / 0.491 | 0.900 / 0.415 |

*(precision / recall)*

**On ordinary sequence the gate works: from 0 to 200 it buys 14.2 precision points and costs
16.7 recall points, a shade over one for one.** At a homopolymer tract, going from 0 to 100
buys **2.3 precision points and costs 27.6 recall points** — twelve recall points a precision
point — and going on to 200 *loses* precision as well, 0.848 back to 0.831, for another 18
recall points. At period 2 and above the gate buys 0.2 points and costs 14.1 by QUAL 100, and
by QUAL 200 it is back where it started at 0.834 with a quarter of the recall gone.

**Precision at a tract is flat from QUAL 10 upwards**: 0.849, 0.850, 0.850, 0.848 at thresholds
10, 30, 50 and 100. There is no threshold that separates ng's true tract calls from its false
ones, because QUAL is not what distinguishes them.

**What does distinguish them is already in the record.** Of the 807 false homopolymer alleles at
30× with no threshold, **107 are alleles no sample was called with** — listed in the ALT column
with `AC=0`. The other 700 are alleles a sample was given and the truth set does not have.

---

## 5. The existing repeat-tract caller, on the same ground

`ssr-call` at 30×, tandem-repeat benchmark, against ng:

| | ng, QUAL ≥ 10 | existing caller, QUAL ≥ 10 |
|---|---|---|
| period 2+, precision | 0.841 | 0.854 |
| period 2+, recall | **0.671** (2 952 of 4 398) | **0.234** (1 028 of 4 398) |
| homopolymers | 0.849 / 0.772 | *not in its catalog* |

**At the same precision ng recovers 2.9 times as many true repeat-tract alleles.** The
comparison is of what each caller emits on this ground and not of their genotypers: the existing
caller's own catalog holds motif periods 2 to 6 and no homopolymers at all
(`benchmarks/ssr_hg002/README.txt`), so its homopolymer row is 33 true alleles against 4,684 and
is a fact about its catalog rather than about its calling. Its low recall at period 2+ is the
emission gap already recorded for it
([`ssr_emission_drop_attribution_2026-07-08.md`](ssr_emission_drop_attribution_2026-07-08.md)).

**Its QUAL is not a continuum.** 1,091 of its 1,939 records on tract ground sit in the QUAL 0–1
bin, where its observed error is 0.058, and it writes nothing at all above 100. So there is no
threshold sweep to compare — its emission decision is made elsewhere, which is precisely what
§3.3 describes arm C as.

---

## 6. The simulator: slippage is not what the residual error is

The simulator is the only place the true slippage is known, so it is the only place the risk can
be tested directly. 4,000 tracts, one sample, ng at `--defaults` — which assumes 10 reads in 100
slip.

**Above QUAL 10 the observed error is zero at every setting tried**, including at 100× and
including when the reads really slip 25 times in 100, two and a half times what the model
assumes:

| true slippage | depth | records above QUAL 200 | observed error |
|---|---|---:|---:|
| 0.02 | 100× | 2 002 | 0.0000 |
| 0.10 | 100× | 2 004 | 0.0000 |
| 0.25 | 30× | 1 437 | 0.0000 |
| 0.25 | 100× | 2 004 | 0.0000 |

**So mispriced slippage does not produce the one-in-500 floor seen on real reads.** Whatever
those nine records are, the simulator says they are not the stutter model being wrong — which
is a direct instruction to whoever designs `calling_quality_ssr.md`: a slippage correction
would be aimed at a failure this measurement does not find.

**What mispriced slippage does cost is alternative alleles nobody carries.** At 30×, going from
a true slippage of 0.10 to 0.25 takes period-2+ precision from 0.883 to 0.665, and **546 of the
846 false alleles have `AC=0`** — listed but not genotyped into any sample. Recall is untouched
(0.999 against 0.996).

**And knowing the true slippage recovers about a fifth of that.** Handed the model the reads
were drawn under, precision at slippage 0.25 and 30× goes from 0.665 to 0.720 at period 2+ and
from 0.654 to 0.710 at homopolymers, for half a recall point (0.996 to 0.992 at period 2+).
At 100× it recovers nothing (0.723 to 0.725). **So the fitted-against-defaulted axis is worth about five precision points at a
tract, and only where the slippage is badly wrong and the depth is moderate.** That is the whole
of what this experiment can say about it: no command fits a parameters file
(§3.4), so on both GIAB grounds every tract call is `Defaulted` — each run's own report says so,
*"0 of 7 groups the file says were fitted"*.

---

## 7. The per-sample benchmark, for continuity

The ground C4's numbers were measured on. Pooled over the three samples, tract ground only:

| depth | arm | homopolymer | period 2+ |
|---|---|---|---|
| 30× | ng, QUAL ≥ 0 | 0.977 / 0.962 | 0.849 / 0.804 |
| 30× | existing caller | 0.992 / 0.947 | 0.978 / 0.804 |
| 50× | ng, QUAL ≥ 0 | 0.977 / 0.962 | 0.855 / 0.839 |
| 50× | existing caller | 1.000 / 0.947 | 0.980 / 0.875 |

*(precision / recall)*

**132 truth alleles at homopolymers and 56 at period 2+, so read the differences as small
counts**: ng's homopolymer precision of 0.977 against 0.992 is 3 false calls against 1. The
existing caller here is the SNP and indel caller at its `high-recall` preset, not `ssr-call`.

---

## 8. What this report does not settle

- **The decision.** `calling_quality_ssr.md` owns it; §1 states the input the decision rule
  asked for and nothing more.
- **Whether QUAL should be Phred at all.** Both paths claim one error in 10²⁰ at QUAL 200 and
  deliver one in 500. That is a property of the site-quality fold, measured here on two grounds
  at once, and it is a bigger question than tracts.
- **The false alleles nobody is called with.** 107 of 807 false homopolymer alleles at 30× and
  546 of 846 on the simulator at high slippage. They are candidate selection's, not quality's,
  and nothing here says what to do about them.
- **Cohorts.** Every ground here is one sample or three. QUAL is a cohort claim, and what it
  does at fifty samples or a thousand is unmeasured.
