# Why thin data loses every stutter error: diagnosis (2026-10-15)

Branch `str-penalised-fit`, head `53a3e4f9`. Diagnosis only: no committed code was changed. The
measurements below were made on an instrumented private copy of the tree
(`tmp/diag_errors/src_copy/`), with the remedies prototyped in Python on the matrices it wrote out.

## The answer

**The joint inversion fails because of the allele lengths, not the stutter curves.** The penalised
fit computes every error of a motif length from one matrix, H: the curvature of the fit's objective
(how sharply the fit's score falls when a number moves away from its fitted value). If H cannot be
inverted, commit `0b12941d` withholds every error of that motif length. In all 56 failing runs I
re-ran, H has a direction along which the curvature is negative, so H cannot be inverted in any
basis. In 54 of the 56, that direction is the trade between the two outermost allele lengths on one
side of the allele window (8 and 9 repeats shorter than the reference, or 8 and 9 longer) at one
thin class. The stutter numbers alone always have positive curvature.

**The cause is a term left out of H.** The fit adds a prior to every class's allele-length shares:
the floor of spec §2.4, one locus's worth of allele draws spread over the 19 lengths plus one draw a
class. The fit maximises the reads' likelihood plus this prior, but H contains only the reads' part
(`penalised_fit.rs:1428-1429`: "The half allele draw a class adds is left out of `H` and `P`"). At a
thin class the prior holds a rare length's share well above where the reads alone would put it. At
that point the reads' likelihood can curve upward along the trade between two rare lengths. Only
the prior's own curvature makes the total curve downward there. With the prior's curvature added, H
inverts in all 384 runs, and the errors it then gives are about as honest as the errors of the runs
that already succeed.

**Recommendation:** add the allele-length prior's curvature to H and to the middle term, as §2.10.4
already does for the stutter prior. It is about 20 lines of code. It returns every error lost in the
56 runs. It leaves the stutter errors of runs that already succeed unchanged in the median, and it
moves the reported share errors of classes with 8 or more loci by a median of 1%. Keep `0b12941d`'s
rule as the fallback for anything that still fails. This changes which errors are reported, so it is
the owner's decision.

## Terms used below

- **H**: the curvature matrix that the errors are read from, `A + P` plus a ridge of 10⁻⁶ on the
  diagonal (`penalised_sandwich`, `penalised_fit.rs:2683-2697`). One row and column for each fitted
  number at each repeat count. Two-base runs have 21 numbers a repeat count (slip rate, shorter
  share, fall-off, and 18 allele-length shares measured against the reference length). Over the 22
  repeat counts of the two-base designs that makes 462 rows.
- **A**: the reads' curvature, one block for each class. It is computed numerically, by moving each
  number 0.001 on its logit scale either way and measuring how the expected counts' slope changes
  (`penalised_fit.rs:1518-1547`).
- **P**: the curves' bend penalties plus the stutter prior's curvature `a·p(1−p)`
  (`penalised_fit.rs:1612-1639`).
- **S**: the curvature of the allele-length prior that the fit maximises but H leaves out. The fit
  adds `allele_pseudo_draws + allele_prior_tracts × (class draws ÷ loci with reads) ÷ 19` draws to
  every length (`penalised_fit.rs:1318-1338`). With `Q` the total of those prior draws over the 19
  lengths and `s` the fitted shares, S is `Q·(diag(s) − s sᵀ)` over the class's 18 share logits.
- **Inverting H.** `dense_inverse` (`penalised_fit.rs:2701-2716`) uses a Cholesky factorisation. It
  stops at the first row whose pivot is zero or negative (`:2710`), and the motif length then gets
  `UndeterminedErrors::EveryNumberButBetweenLocusShare` (`:1690-1695`).
- **Honest error**: one that holds the truth within one error about 68 times in 100, and within two
  errors about 95 times in 100. Coverage below is written "within one / within two", for example
  0.65/0.88. Errors are compared raw, before any withholding rule, as in D1's tables.
- **Thin class**: a class with fewer than 8 loci with reads, the threshold below which spec §2.6
  withholds every error.

## The evidence

### 1. The failures reproduce exactly

I re-ran all 384 of D1's coverage jobs with at most 500 samples (`tmp/d1/cov.txt`, the same
arguments; 8 threads). I did not re-run the 96 jobs at 2,000 samples. In the re-run, the same 56
runs as in D1 got no stutter error, run for run (`tmp/diag_errors/failing_rerun.txt` against D1's
list `tmp/diag_errors/failing.txt`). The 33 failing runs at 2,000 samples were not examined.

My Python reconstruction of the production errors from the matrices written out agrees with the
bench's own errors to a relative difference of 8×10⁻¹⁴, over 924 errors in the three
`p2_n60_d10_pcr_free` runs. So the remedies below are measured on the production matrices, not on
an approximation of them.

### 2. Where H fails, and along what

**Two-base thin-tail design, 60 samples × 3 reads, seed 1** (`tmp/d1/mixed_p2.toml`; repeat counts
6–14 with 33–100 loci, 15–27 with 1–7):

- H has 462 rows. The Cholesky factorisation fails at row 234 (repeat count 17, the share 9 repeats
  shorter than the reference), with a pivot of −9.07.
- The smallest eigenvalue of H is −0.385. The next ones are 0.69, 1.15 and 1.24, and the largest is
  1.3×10⁵. So exactly one direction is negative, and it is not near zero.
- That direction's weight: repeat count 16's shares 8 and 9 shorter (−0.66, +0.62), with smaller
  weights on the same two shares at repeat counts 15 and 17. It is the trade "more of length −9,
  less of length −8" at a 2-locus class, spread a little to its neighbours by the shares' bend
  penalty.
- Four classes have a negative eigenvalue in their own A block, all along the same −8/−9 trade:
  - repeat count 16 (2 loci): −2.24;
  - repeat count 21 (7 loci): −0.05;
  - repeat count 24 (3 loci): −0.24;
  - repeat count 26 (5 loci): −1.10.
  With S added, each is positive: 1.37, 1.70, 0.94 and 0.61.
- The stutter numbers alone (slip rate, shorter share and fall-off at all 22 repeat counts, 66 rows)
  have smallest eigenvalue 19.8. Nothing there is uninformed.
- With S added, H's smallest eigenvalue is 1.76 and the factorisation succeeds.

Seeds 2 and 3 of the same cell look the same. Seed 2 fails at repeat count 17 (3 loci), smallest
eigenvalue −0.60. Seed 3 fails at repeat count 25 (4 loci), −0.77. Both negative directions lie on
the −8/−9 trade. The stutter block's smallest eigenvalues are 18.6 and 17.5. With S added, H's
smallest eigenvalues are 1.61 and 1.41.

**Two-base, 60 samples × 10 reads, PCR-free, seed 2** (`p2_n60_d10_pcr_free`, the run of the
milestone-C fix review): the factorisation fails at row 361 (repeat count 23, share −8), pivot
−0.117. The smallest eigenvalue is −0.206, carried by repeat count 23's shares −8 (+0.88) and −9
(−0.43). That class has 16 loci. Its expected allele draws are 72 at −9 and 33 at −8, where the
truth is about 0.003 of each, so a few loci carry these lengths by chance. Its A block's smallest
eigenvalue is −1.16, and +1.35 with S. The shorter share's curve is stiff here (strength 8.5×10⁹, so
H's largest eigenvalue is 1.35×10¹¹), but the stutter block's smallest eigenvalue is 18.6, so the
stiffness is not the cause. This confirms the milestone-C review's suspicion (M1) that this run's H
is genuinely not positive definite. Seeds 1 and 3 of the cell succeed. Each has one or two classes
whose A block is negative along the same kind of trade (−0.22 at repeat count 13, and −0.62 and
−0.04 at 15 and 25). There the bend penalty from the neighbouring repeat counts is just enough to
cover it: seed 3's smallest eigenvalue of H is 0.255.

**All 56 failing runs** (`tmp/diag_errors/tally.txt`): the negative direction's two heaviest numbers
are:

| what carries the negative direction | runs |
|---|---|
| shares −8 and −9 of one class | 37 |
| share −8 of one class alone | 11 |
| shares +8 and +9 | 4 |
| share +8 alone | 2 |
| share +1 (1 sample × 3 reads, a 1-locus class) | 1 |
| share −4 (1 sample × 3 reads, a 2-locus class) | 1 |

- The smallest eigenvalue of the stutter block, over the 56 failing runs, is 0.09, at 1 sample × 3
  reads. It is above 1 in every run with 10 samples or more.
- **The smallest eigenvalue of H with S added, over all 384 runs, failing or not, is 0.005** (1
  sample × 3 reads, a 1-locus class). It is positive in every run.

So the brief's hypothesis, a stutter-curve direction with no information, is not what happens. The
failing direction carries negative curvature of 0.04 to 36 units. It is not a near-zero one. It
lies in the allele-length shares. The stutter errors are lost only because `0b12941d`'s rule blanks
the whole motif length when any direction fails.

### 3. Negative curvature, not a point on a bound or a numerical artefact

- **Not numerical.** The negative eigenvalues are up to −36 at 500 samples, against a ridge of
  10⁻⁶ and a Cholesky rounding of about 10⁻⁵ even at stiffness 10¹⁰ (milestone-C review, M1). A
  change of basis or a rescaling cannot help: a symmetric matrix keeps its number of negative
  eigenvalues in every basis (Sylvester's law of inertia).
- **Not a bound.** The shares are a softmax of unbounded logits and have no clamp. The two clamps in
  `numbers_of` (fall-off at 0.99, F's concentration) touch only the fall-off and F, and neither
  carries the failing direction.
- **Not a missing stutter prior.** The stutter prior is in H (`penalised_fit.rs:1627-1639`, spec
  §2.10.4).
- **A missing allele-length prior.** A is the curvature of the reads' likelihood alone. The fit's
  point maximises likelihood + allele-length prior + bend penalties, so the matrix that must be
  positive there is A + S + P. Away from the likelihood's own maximum, a likelihood over unseen
  genotypes can curve upward. At a thin class the prior's 3 to 4 draws a length hold every rare
  length's share well above what the reads alone would give. The measured result fits this
  reading: A + P is negative in 56 runs and A + S + P is positive in all 384.
- **Why the window's outer lengths in particular** is not established. One guess, unchecked: alleles
  beyond the ±9 window must be explained by the edge lengths, and reads separate the −8/−9 split
  poorly. Either way, the remedy does not depend on the answer.

## The options, measured

Each remedy was applied to the production matrices of all 384 runs (`tmp/diag_errors/remedy.py`).
The errors were mapped to each number's own scale exactly as `standard_errors` does, and compared
with the truth (`tmp/diag_errors/summary.py`, output `tmp/diag_errors/summary.txt`). Coverage is
within one / within two errors, with n the number of class-numbers. The shares count only lengths
that the owner's rare-length rule reports (share × loci ≥ 3).

### (c) Add the allele-length prior's curvature S to H and to the middle term (recommended)

**Errors come back in all 56 failing runs.** H + S inverted in every one of the 384 runs.

Coverage in the 56 failing runs, with the remedy:

| number | 1–2 loci | 3–4 loci | 5–7 loci | 8+ loci |
|---|---|---|---|---|
| slip rate | 0.53/0.83 (330) | 0.59/0.86 (313) | 0.65/0.88 (239) | 0.56/0.91 (57) |
| shorter share | 0.55/0.88 | 0.50/0.78 | 0.54/0.84 | 0.58/0.82 |
| fall-off | 0.66/0.90 | 0.57/0.88 | 0.61/0.88 | 0.82/1.00 |
| shares (reported) | – | 0.17/0.58 (12) | 0.52/0.77 (100) | 0.69/0.95 (172) |

For comparison, the 328 runs that already succeed, under today's code:

| number | 1–2 loci | 3–4 loci | 5–7 loci | 8+ loci |
|---|---|---|---|---|
| slip rate | 0.49/0.80 (1,507) | 0.56/0.86 (1,522) | 0.56/0.87 (2,441) | 0.65/0.93 (327) |
| shorter share | 0.52/0.83 | 0.58/0.91 | 0.55/0.89 | 0.60/0.91 |
| fall-off | 0.59/0.86 | 0.61/0.89 | 0.61/0.89 | 0.63/0.90 |
| shares (reported) | – | 0.12/0.68 (34) | 0.54/0.82 (1,132) | 0.63/0.93 (1,045) |

The recovered errors are about as honest as those of the runs that never failed, cell for cell. The
largest shortfall is the shorter share at 3–4 loci: 0.50/0.78 against 0.58/0.91. Both tables fall
short of honest (0.68/0.95) at under 8 loci, as D1 already found. That shortfall is the thin-class
question, which the 8-loci rule handles, not this defect.

**Effect on the runs that already succeed** (error with S divided by today's error):

- Slip rate, shorter share and fall-off at classes with 8+ loci: median 1.000. The 5th to 95th
  percentiles lie within 0.97–1.001, so the change is under 3%.
- The same three numbers at thin classes: median 1.000. The 5th to 95th percentiles lie within
  0.84–1.02. A few move more (from 0.06 to 1.75), all in classes whose errors the 8-loci rule
  withholds.
- Shares at classes with 8+ loci: median 0.99, 5th–95th percentile 0.96–1.12. Coverage of the
  reported ones goes from 0.63/0.93 to 0.63/0.92 (1,045).
- Shares at thin classes: median 0.89, 5th–95th percentile 0.57–3.4. Coverage at 5–7 loci falls
  from 0.54/0.82 to 0.48/0.75 (1,132). Today these errors are withheld by the 8-loci rule. This
  matters only if the owner lowers that rule for shares.
- The stutter coverage of these runs moves by at most 0.03 in any cell.

**Variant: S in H only, not in the middle term.** The coverage is the same to within 0.01
everywhere. Putting S in both follows the stutter prior's convention (§2.10.4), which keeps a number
with no reads from getting an error near zero.

**Code.** About 20 lines in `standard_errors`. For each class, compute `Q` as `count` does
(`penalised_fit.rs:1318-1338`; best moved into one helper that both call, so the two cannot drift),
and add `Q·(diag(s) − s sᵀ)` into `penalty` over the class's share block, next to the stutter prior
(`:1627-1639`). Two doc comments need changing:

- `:1428-1429` should say the prior is included;
- `UndeterminedErrors` (`:540-570`) should drop "a share pushed to nearly 0 or 1" as a typical
  cause, or keep it as the residual case.

The fitted values do not change. Only the errors change, which moves the oracle's pinned checksums.

### (a) Drop only the directions that are not positive

Eigen-decompose H, invert only on its positive directions, and withhold the error of any number
with more than 0.1% of its weight on a dropped direction. I tried two thresholds:

- **Drop eigenvalues ≤ 10⁻⁸** (only the non-positive ones):
  - Runs that already succeed are unchanged, since they have no such direction.
  - In the failing runs, 927 of 939 slip-rate errors and 928 of 939 shorter-share errors return,
    against all 939 under (c). The withheld ones are mostly at 1–2-locus classes, the class carrying
    the bad trade and its neighbours.
  - For the three stutter numbers, coverage is within 0.07 of variant (c) in every cell, for
    example slip rate at 5–7 loci 0.66/0.89. For the reported shares it is higher within two errors
    (0.89 against 0.77 at 5–7 loci, 95 against 100 shares), because the shares of the class carrying
    the bad trade are withheld and keep the errors of a matrix without S.
- **Drop eigenvalues ≤ 10⁻⁹ × the largest** (a relative threshold, the usual choice): this
  withholds about 6 errors in 10 even in runs that already succeed. A stiff curve makes the largest
  eigenvalue 10¹¹, so the threshold reaches 100 and drops well-informed directions. **It is
  unusable here.** This is the tolerance risk the milestone-C review predicted.

**Code.** Production needs a portable eigen-decomposition (Jacobi rotations, about 100 lines, about
a second a motif length at 462 rows), an absolute tolerance constant, and a rule that says which
numbers a dropped direction touches. **Weakness:** it reports errors conditional on a direction
that the fit's own curvature calls wrong, and it treats the symptom: the matrix it inverts is still
missing the prior the fit used.

### (b) Invert each curve in its own basis

**Cannot help here.** That remedy fixes rounding along a stiff curve's straight line. These H's
are not positive definite in exact arithmetic, and no change of basis changes the sign of an
eigenvalue. The stutter curves, the only stiff ones, are well informed in every failing run
(smallest stutter-block eigenvalue 0.09 at one sample, above 1 elsewhere). It may still be worth
having for the review's cause 1 (a stiff curve through one repeat count with slipped reads), but no
run in this set shows that case.

### (d) Narrow `0b12941d`'s rule

The rule could withhold only the failing block, here the shares of one class, instead of the whole
motif length. With (c) in place nothing fails in these 384 runs, so this would be a safety net with
no measured use yet. It is not recommended now.

## Recommendation

1. **Adopt (c).** It corrects H to the curvature of the objective the fit actually maximises, so it
   is a bug fix rather than a rescue. It returned every lost error in the 56 runs at 500 samples or
   fewer, with coverage matching the runs that never failed. It changed today's stutter errors by
   under 3% at classes with 8+ loci, and today's reported share errors by a median of 1%. It is
   about 20 lines.
2. **Keep `0b12941d`'s rule** as the fallback. Log, as today, any motif length that still fails.
3. **Before merging, re-run the 33 failing runs at 2,000 samples × 3 reads.** The same mechanism is
   likely there but was not checked, by instruction. Re-run the oracle-like check as well, since
   the errors move.
4. **Correct the stale statements about the cause** in the handoff (§4), `tmp/d1/RESUME.md` (open
   observation 1) and the `UndeterminedErrors` comment: thin classes lose their stutter errors
   because of an allele-length trade at the window's edge, not because of an uninformed stutter
   direction. `LOCI_FOR_AN_ERROR` is unrelated, as D1 found.

## How to reproduce

The instrumented copy adds one block before `penalised_sandwich` in `standard_errors`. When
`PVC_DUMP_H` is set, it writes `<PVC_DUMP_H>.p<period>.json` with each class's A and B blocks, P as
triplets, and the fitted shares, counts and curves.

```
cd /Users/jose/devel/pop_var_caller-str-penalised-fit
mkdir -p tmp/diag_errors/src_copy && git archive 53a3e4f9 | tar -x -C tmp/diag_errors/src_copy
# apply the dump block (tmp/diag_errors/src_copy/src/parameter_estimation/joint/ssr_fit/penalised_fit.rs,
# search PVC_DUMP_H), then:
(cd tmp/diag_errors/src_copy && CARGO_TARGET_DIR=$PWD/../target \
   cargo build --release --features bench-fixtures --example ng_str_fit_bench)
# one run:
RAYON_NUM_THREADS=8 PVC_DUMP_H=tmp/diag_errors/cov/mixed_p2_n60_d3_s1 \
  tmp/diag_errors/target/release/examples/ng_str_fit_bench --period 2 --samples 60 --depth 3 \
  --seed 1 --truth tmp/d1/mixed_p2.toml --methods D_tract_prior \
  --truth-out tmp/diag_errors/cov/mixed_p2_n60_d3_s1.truth.tsv \
  --fixture-out tmp/diag_errors/cov/mixed_p2_n60_d3_s1.json
uv run --no-project --with numpy python tmp/diag_errors/remedy.py diag \
  tmp/diag_errors/cov/mixed_p2_n60_d3_s1.p2.json
# all coverage jobs at <= 500 samples, then the tables:
uv run --no-project python tmp/diag_errors/drive.py tmp/d1/cov.txt tmp/diag_errors/cov
uv run --no-project --with numpy python tmp/diag_errors/remedy.py rows tmp/diag_errors/cov_rows.tsv tmp/diag_errors/cov
uv run --no-project python tmp/diag_errors/failing.py tmp/diag_errors/cov > tmp/diag_errors/failing_rerun.txt
uv run --no-project --with numpy python tmp/diag_errors/summary.py tmp/diag_errors/cov_rows.tsv tmp/diag_errors/failing_rerun.txt
uv run --no-project --with numpy python tmp/diag_errors/tally.py tmp/diag_errors/cov
```

The `p2_n60_d10_pcr_free` runs are in `tmp/diag_errors/grid/`, made by `drive.py` from
`tmp/diag_errors/extra.txt`. These runs omit `--both-fallbacks`, which affects only the bench's
fallback methods.
