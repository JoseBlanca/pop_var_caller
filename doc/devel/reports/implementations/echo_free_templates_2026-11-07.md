# A sample does not read its own reads back (2026-11-07)

Branch `echo-free`, from main f1097328. Spec changed: `neighbour_evidence_phasing.md` §2.2, §2.3, §5.

## 1. What was wrong

`fixed_f_2026-11-06.md` §5 left a loss untraced. In a 40-sample outbred cohort with every sample at
2×, refinement got fewer heterozygotes right than calling alone: 5,007 of 13,709 against 7,521.
Counting the rounds averaged showed the loss growing after round 2:

| rounds averaged | heterozygotes right |
|---|---|
| 2 only | 7,111 |
| 2–3 | 4,636 |
| 2–5 (the pass's) | 5,007 |

The heterozygotes lost were called homozygous **evenly on both sides, at every allele frequency**.
Where the alternative allele is carried by 10–20% of the cohort's copies, 429 were called
homozygous alternative, against 45 with refinement off. A pull towards the cohort's frequencies
would push rare alleles towards the reference. An even split is what a sample's two copies settling
on whichever allele its one or two reads showed would look like.

**The cause is an echo.** Sample A copies sample B in round `r`, so A's haplotype there carries B's
alleles. In round `r + 1`, B copies A — the template choice ranks haplotypes by how little they
disagree with B's genotype, and A's resembles it most — so B's haplotype now carries A's alleles,
read from A's reads. In round `r + 2`, A reads B's haplotype back as a relative's and counts its own
reads a second time. Each round tightens it. Where deep samples hold the panel their reads dominate
and the echo does little; with none it decides the genotype.

## 2. What was tried

All on the 40-sample all-2× cohort and the 40-sample outbred cohort at 30×, 10×, 5× and 2× (ten
each), with switches in a diagnostic build (tag `het-loss-diag`). Heterozygotes right; for the
mixed cohort at 2×, and its phase.

| change | all 2× | mixed, 2× | mixed, phased right, wrong links |
|---|---|---|---|
| main | 5,007 | 4,378 | 14,584, 9 |
| stored haplotypes drawn hard, one allele pair from each joint | 5,124 | | |
| stored haplotypes the most likely pair | 4,326 | | |
| stored haplotypes' joint squared | 4,204 | | |
| the panel frozen after round 1 | 7,215 | 4,258 | 15,048, 82 |
| the panel frozen after round 2 | 3,890 | 4,326 | 14,964, 21 |
| a sample never copies one that copied it the round before | 6,836 | 4,276 | 12,531, 4 |
| **a copy that followed the reader is read without it** (adopted, §3) | 7,042 | 4,383 | 14,750, 8 |

Hardening the haplotypes does not help: softness is not the cause. Freezing the panel after round 1
removes the chain and most of the loss, but costs the mixed cohort its phase. Barring whole samples
removes the echo but also a deep sample's best relatives, which are often mutual. Reading only the
echoed copy without the reader keeps both.

## 3. The change

| what | where | spec |
|---|---|---|
| Each round stores, per sample, site and copy, the template the path followed, and per sample the samples its templates came from | `DrawnPath::followed` (`model.rs`), `HaplotypeStore::followed`, `sources` (`refine.rs`) | phasing §2.3, §5 |
| A copy that followed the reader, where the reader's reads leave its genotype unsure, is read with the followed template replaced by the cohort's frequencies — in the template choice and in the walk | `HaplotypeStore::template_into`, `ECHO_SETTLED_LIKELIHOOD` = 0.01 | phasing §2.2, §2.3 |

The replaced copy is `p(a) · Σ_b h(b) · L_t({a, b})`, normalised, with `p` the cohort's frequencies,
`h` the copy-owner's other copy and `L_t` its likelihood (flat at its hidden records); where both
copies followed the reader, `p(a) · L_t({a, a})`. Recomputing the copy exactly needs a second stored
haplotype per copy; the diagnostic build did so and gave the same genotypes within noise (all 2×:
6,991 against 7,105; mixed: 19,568 against 19,553 heterozygotes right over all depths).

**Two choices, measured.**
- *Only where the reader's reads leave its genotype unsure* (second most likely genotype above 1/100
  of the most likely). Removing the echo everywhere gave the 22-sample F2 6,483 heterozygotes right
  and 6,242 phased; with the limit 6,545 and 6,358. Limits of 1/1,000 and 3/100 were within noise of
  1/100 everywhere.
- *A sample's hidden records are judged by its reads* (the records hidden from its own walk for the
  fill check). Judged unsettled, the F2 fell to 6,413 and 6,164: a 30× parent's hidden records read
  its offspring's copies without its reads, and its haplotypes there went soft for every offspring
  copying it. The copy-owner's hidden records do read a flat likelihood, so its held-out reads cannot
  reach its own fill check through the reader.

## 4. Measurements

The branch against main, refinement on; genotypes right/wrong summed over the cohort's depths;
phased right, links, wrong links.

| cohort | main: heterozygotes | branch | main: phase | branch |
|---|---|---|---|---|
| outbred, 40, all 2× | 5,007/7,045 | **7,042/5,007** | 801, 348, 5 | 1,030, 421, 1 |
| outbred, 40, 30×/10×/5×/2× | 19,562/532 | 19,575/529 | 14,584, 9,842, 9 | 14,750, 10,022, 8 |
| inbreeding 0, 40, 30×/20×/10×/5× | 20,485/165 | 20,485/165 | 16,227, 11,763, 7 | 16,192, 11,730, 2 |
| inbreeding 0, 12 of those | 5,996/73 | 6,000/70 | 3,920, 2,437, 8 | 3,865, 2,378, 7 |
| inbreeding 0.5, 40 | 9,613/42 | 9,614/39 | 8,487, 7,457, 7 | 8,500, 7,469, 5 |
| selfing (0.9), 40 | 2,002/42 | 2,009/38 | 1,819, 1,638, 0 | 1,813, 1,627, 0 |
| F2, 22 (2 parents at 30×, offspring at 1×) | 6,507/97 | 6,516/105 | 6,264, 6,048, 2 | 6,320, 6,068, 2 |

In the all-2× cohort the homozygotes lose a little: alternative 7,953 → 7,887 right, reference
16,637 → 16,449; all three classes together 29,597 → 31,378 right, against 29,381 with refinement
off. The mixed cohort's homozygous references called wrong rise from 139 to 159.

**With the echo removed the rounds no longer drift** — heterozygotes right in the all-2× cohort by
rounds averaged (diagnostic build, the same rule):

| rounds | main | echo removed |
|---|---|---|
| 2 | 7,111 | 7,453 |
| 3 | 4,636 | 7,412 |
| 5 | 5,007 | 6,991 |
| 8 | 5,295 | 7,109 |

The later rounds still help where deep samples are: the mixed cohort's 5× samples called wrong 189
heterozygotes over rounds 2 alone and 141 over rounds 2–5.

**Small cohorts.** Subsets of 4 and 12 samples, one to three at each depth, each with its own
parameters file:

| subset | main: phased right, links, wrong | branch |
|---|---|---|
| outbred, 12 | 3,183, 1,776, 7 | 3,171, 1,807, 15 |
| outbred, 4 | 490, 218, 12 | 474, 201, 15 |
| inbreeding 0, 4 | 711, 351, 13 | 687, 339, 19 |
| inbreeding 0.5, 4 | 488, 281, 11 | 487, 281, 12 |
| inbreeding 0.8, 4 | 310, 238, 12 | 299, 228, 9 |

Genotypes there move by under 1 in 200. Wrong links rise in three of five, by 3 to 8; the 4-sample
cohorts were already wrong 3.7 to 5.5 times in 100 links (`phase_background_gate_2026-11-06.md`).

**Time** on the 40-sample inbreeding-0 cohort, 6 threads: 9.4 and 9.8 s on main, 9.8 and 9.9 s on
the branch. **Memory:** 4 bytes a sample-record more (§5 of the spec), about 0.5 GB at 3,000 samples.

**refine_sim on rick, 400 samples:** pending — rick was running the benchmark suite.

## 5. Tests

- `refine::tests::a_copy_that_followed_the_reader_is_read_without_it` — 3 samples at 3×: a copy
  that followed the reader is read as the formula of §3 where the reader's reads are unsure and
  whole where they settle it; the other copy, and any other reader, read it whole; both copies
  following the reader give `p(a) · L({a, a})`.
- No test reproduces the drift itself: the crate's simulated cohort (`test_cohort`) does not drift
  at 2× with 8, 16, 32 or 64 founders — heterozygotes right rose from 2 rounds to 8 in each, with
  the echo left in.
- The two bit-fingerprint tests of a phased window re-recorded; the test window's phased columns
  rise from 2,782 to 2,834.
- `templates.rs`: the landing weight and the stretch scale now use `crate::float::powi` and
  `crate::float::ln` — clippy's disallowed methods, since std's round differently by platform; no
  fingerprint moved here. A separate commit.
- The whole crate's tests pass in release, except the two `should_panic` tests outside the pass that
  panic only through debug assertions.

## 6. Reproducing

Cohorts as in `fixed_f_2026-11-06.md` §7 and `phase_agreement_default_2026-11-06.md` §4; scored
with `benchmarks/suite/lib/refine_score.py` for phase and a by-depth genotype count against the
truth. The diagnostic build's switches (tag `het-loss-diag`): `PVC_DIAG_HAP`, `PVC_DIAG_FREEZE`,
`PVC_DIAG_NOECHO`, `PVC_DIAG_ECHO` (`1`, `walk`, `choice`), `PVC_DIAG_ECHO_GL`,
`PVC_DIAG_ECHO_APPROX`, `PVC_DIAG_ROUNDS`.
