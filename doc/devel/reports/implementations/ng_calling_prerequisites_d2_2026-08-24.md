# ng — calling prerequisites, D2: the error-rate scale's denominator

**2026-08-24**, branch `ng-calling-prerequisites`. Step D2 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.2 — **which this step corrected
twice, on the owner's ruling, before any of it could be built.**

---

## 1. The step as written could not be built, and why

The caller will charge each read the error probability the walk minted for it, rescaled by one
number per library so that the average comes out at the rate the pre-pass measured. The pre-pass
fits that rate. The step was to add the denominator: the average minted error over the same reads.

**§3.2 asked for "a running sum of the per-read error probability, and the count of reads it was
summed over", and nothing carries the first of those.** The fold derives one number per read — the
worse of its window's base quality and its mapping quality — sums the *logarithms* into an
observation's `q_sum`, and throws the individual reads away. `Σ ε` is not recoverable from `Σ ln ε`,
so supplying it meant a second accumulation at fold time and a new field on every observation,
touching about 89 construction sites.

**The owner's ruling: take the geometric mean, and the whole difficulty disappears.** The deciding
fact is what the scale is applied *to*. The model charges an observation `exp(q_sum / num_obs)` — a
geometric mean — and so does production, clamped
([`posterior_engine.rs:1536`](../../../../src/var_calling/posterior_engine.rs)). **Production has no
recalibration at all**, so "what does production do here" has no answer to copy; what it has is the
quantity, and it is this one. A scale built from an arithmetic mean and applied to a geometric one
would not make the calibrated property hold in the model's own terms, so paying for the arithmetic
sum would have bought an inexactness rather than removed one.

So the two numbers §3.2 asks for are `Σ q_sum` and `Σ num_obs` per read group — **already carried**,
and the step becomes adding up numbers the walk has produced rather than minting anything.

**The second correction: only one of the two routes can carry it.** The step asked both. The census
route's per-position record is a depth code and a sparse list of allele counts with no quality in it
at all, so it cannot supply either number under either definition of the average. Its accumulator now
waits on the open comparison between the two routes: if it wins, its records gain a quality field; if
the histogram route wins, nothing is owed.

## 2. Changes made

**[`src/ng/parameter_estimation/generic/calibration.rs`](../../../../src/ng/parameter_estimation/generic/calibration.rs)**, new.

- `MintedReadErrors` — one read group's summed log error and read count, with the two means it
  implies. `None` at zero reads, because a scale needs a denominator and there is none.
- `minted_error_by_read_group` — one locus's totals per read group, over the same observations under
  the same gate as the error-rate histogram.
- `fold_into` — a locus's totals into a running table.

**[`generic/accumulators.rs`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)** —
`GenericAccumulators` accumulates them in `add_locus` beside the read-group histogram, merges them
across region shards, and exposes them with `minted_errors()`. **Beside the histograms rather than
inside them**, because nothing is fitted from this: the fit reads the tables and the scale divides
this into the fit's answer afterwards.

**[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.2** — both corrections, with
their reasons.

Six tests, plus two assertions on existing accumulator tests.

## 3. The one thing in the implementation that is not obvious

**The running sum is a fixed-point integer, in units of 2⁻²⁰, and not an `f64`.**
`GenericAccumulators::merge` promises order-independent merging across region shards — its own test
merges three shards in all six orders and asserts the results agree — and that holds for its tables
because they sum integer counts. Floating-point addition is not associative, so a running `f64`
would have broken that contract in a way no existing test could have seen, and the symptom would
have been a run whose genotypes depended on the order its shards happened to finish in.

**The cost is bounded on the number that is used.** Each conversion rounds by at most half a unit
and each covers at least one read, so however many loci and shards a run has, the error on the
*mean* stays under 2⁻²¹ ≈ 4.8 × 10⁻⁷ — against a mean log error of order 5 to 20. The accumulator
test asserts exactly that bound, and the miss it measures is 1.2 × 10⁻⁹.

## 4. What checking my own work changed

**The server-side outage that hit this step's reviews is why this section exists**: three of four
review agents died on API errors, so the check that the denominator runs over the same reads as the
fit — the requirement §3.2 spends most of its length on — I made myself. What it found:

**The two site sets are identical, by construction rather than by luck.** Both paths run behind the
same `LocusKind::Generic` gate in `add_locus`, both iterate `locus.complete_observations()`, and both
count `num_obs`. So a repeat-tract locus, a partial read and an observation no read is behind are
excluded from both, and no gate exists on one path and not the other.

**The grain matches the fit's answer and not the histogram's key.** The histogram is keyed
`(read group, ploidy)`; these totals are keyed by read group alone. That looked like a mismatch and
is not: the fit gathers a group's cells across every ploidy it covered into **one** scan and returns
one rate per group, because chemistry does not know about chromosomes. A denominator keyed by ploidy
would have been the mismatch.

**The depth cap is quality-blind, which is what the argument needed.** `CountedSite::capped` draws
hypergeometrically from `(depth, alt_reads, cap)` — three counts, no quality anywhere — so the reads
it notionally keeps are a uniformly random subset, and the mean log error over all reads is an
unbiased estimate of the mean over the kept ones. **Unbiased and not equal**, which is the honest
statement: the cap keeps a count, never an identity, so the kept subset's own mean is unobservable.

**One real defect, and widening the sum removed it.** The running total was an `i64`. One read group
of one human sample reaches about 2.2 × 10¹⁶ scaled units, which an `i64` holds four hundred times
over — but `add` is documented to fold across *samples* as well as shards, because a read group is a
library and a library can hold more than one plant. Past about four hundred human-scale samples in
one read group an `i64` saturates, and `saturating_add` pins rather than panicking, so the symptom
would have been a mean that was merely wrong. It is an `i128` now: 7 × 10²¹ such samples, beyond any
run, and the map holds one entry per read group so the width costs bytes.

### What the reviews changed — nothing, because they could not run

**Seven review agents died on `529 Overloaded` across this step**, in two waves: all four of the
first launch, then three of four relaunches. The two that got furthest reported one line each before
failing. So this step has **no independent review**, and §4 above is what I checked myself in their
place — the site-set comparison, the grain, the cap's quality-blindness, and the overflow bound that
turned out to be a real defect.

**What that leaves unchecked, named rather than glossed:** whether any of the seven new tests is
redundant with another; whether a genuinely `f64` running sum is distinguishable at fixture scale, so
whether the integer choice rests on a test or only on the argument in §3; and an independent pass
over the prose, which in five of the last six steps found a claim that was wrong. **The measurement I
promised the owner — how far the geometric and arithmetic means actually sit apart on real reads — is
also not done**, because it was one of the agents that died. It is the first thing to run when
capacity returns, and until it does, §3.2's choice rests on the argument that the geometric mean is
what the model consumes, not on a number.

**One near-miss worth recording, because the rule that saved it exists for exactly this.** A
`perl -0pi -e` substitution truncated `calibration.rs` to zero bytes — silently; `cargo fmt`
succeeded on the empty file and only a test run showing *0 tests* caught it. The saved review diff
was what recovered the 388 lines. The session's standing rule is never to `git checkout` a file to
undo an edit; **it should say never to run `perl -i` or `sed -i` on a source file either**, and the
next handoff says so.

## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::parameter_estimation::generic::calibration` | **7 passed, 0 failed** |
| `cargo test --lib ng::parameter_estimation::generic::accumulators` | **16 passed, 0 failed** |
| `cargo test --lib` | **4,181 passed, 0 failed, 14 ignored**, 735.38 s |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 redundant-explicit-link-target warnings — the recorded baseline, unchanged |

**Every mutation quoted below was applied to this tree and run.** Dropping the summed log error, and
counting observations instead of reads, each redden the tests written for them — the second in three
places. **Two mutations turned out to be behaviour-neutral and are reported as such rather than as
gaps:** routing the fold through an `f64` round trip changes nothing at fixture scale, because
`i128 → f64 → i128` is exact below 2⁵³ scaled units, and the same is true of the merge. That is why
§3's claim is argued rather than measured, and why the review that was to test it harder is the one
still owed.

## 6. Follow-ups

- **The two averages, measured on real reads.** Owed to the owner and not done: the review agent that
  was to run it died. Until it is, the geometric choice rests on the argument that it is what the
  model consumes.
- **An independent review of this step**, which nothing has had. The three categories that would
  have run are in `/Users/jose/devel/d2_review_brief.md`.
- **The census route's accumulator** waits on the comparison between the two error-rate routes. Its
  records carry no quality, so if that route wins they gain a field; if the histogram route wins,
  nothing is owed.
- **Nothing consumes the accumulator yet.** The consumer is the read likelihood's error-rate scale,
  in [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md).
