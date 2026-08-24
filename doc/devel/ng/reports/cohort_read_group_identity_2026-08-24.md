# ng — sixty-three libraries merged into one, and the refusal that stops it

*A defect in how the cohort route identified its libraries, the guard that now refuses it, and what
the guard cost the fixtures it was checked against. 2026-08-24.*

---

## 1. What was wrong

**A run that assembled several samples into one cohort could merge every one of their libraries into
a single library, silently, and then fit one sequencing-error rate where it should have fitted one
per library.**

A library — one `@RG` record, one DNA preparation of one plant — is identified by the position its
record takes in the run's flat list of every input file's `@RG` records
([`build_read_groups`](../../../../src/ng/read/input/read_groups.rs)). That makes the identifiers
unique across *the paths that call was given*, which is what its own documentation promises. Called
with one file at a time, the promise has nothing to act on: every single-library sample comes back
holding identifier `0`.

Two callers did exactly that — the shipped `estimate-contamination` subcommand
([`estimate_contamination.rs`](../../../../src/pop_var_caller_exp/estimate_contamination.rs)) and
the tool that drives the cohort fit on real alignments
([`ng_joint_records_walk.rs`](../../../../examples/ng_joint_records_walk.rs)) — and the cohort took
them without a word: it collected the samples' identifiers into a set, so sixty-three samples all
claiming `0` became one entry.

**What it cost, and why nobody saw it.** The fit keys its sequencing-error rates on the library
alone, so sixty-three libraries claiming one identifier were fitted as one and reported as one rate.
The contamination half survived intact, because it is keyed by sample *and* library; and the one
shipped consumer of this path emits only contamination and discards the rates. So nothing wrong
reached a file. What it cost is that **the cohort route's per-library error rate had never been
produced by any driver that existed** — and that is the number the read likelihood would consume if
this route is the one kept.

The merged answer does not look wrong, which is what made it dangerous. Blend the ordinary and
mismapped rates by the share of positions that are mismapped and you get one number a reader can
compare across the two routes; on the 63-accession tomato panel that number was 3.9813 × 10⁻³ from
the merged cohort route against 3.9721 × 10⁻³ from the per-sample route's read-weighted mean over
the same reads — **two parts in a thousand apart**, inside the precision the tool prints (both
measured 2026-08-24, before this branch). A reader comparing them would have concluded the two
routes agree. What they could not see is that one of the two is sixty-three libraries averaged into
one, whose own rates are measured below to span **16.45-fold**.

---

## 2. The rule the guard rests on

**A read group belongs to exactly one sample.** The run's read-group table files every identifier
under the single sample its `@RG SM` names, so two samples holding one identifier is never a cohort
that happens to overlap — it is always a cohort whose samples had their identifiers minted
separately. There is no legitimate case to carve out.

One near-miss that looks like a counter-example and is not: two files may declare the *same library
name* and that is ordinary — one preparation sequenced across two lanes — but they are two `@RG`
records and get two identifiers. A repeated **name** is fine; a repeated **identifier across
samples** is not.

**Why the check cannot be satisfied by a cohort assembled wrongly.** Every call to the identifier
minter numbers from zero. So a run that mints once over every alignment gives every library a
different identifier and cannot collide; a run that mints a file at a time gives *every* sample's
first library the identifier `0` and always collides, at `0`, as soon as there are two samples.
There is no arrangement of separate mintings that produces disjoint identifier sets. That is why an
index is all the evidence has to carry for the question to be answerable at the cohort's door — the
census needed no new field. The one case the check cannot see is a cohort of **one** sample, where
there is nothing to merge and nothing to catch.

---

## 3. What was built

**The guard.** `CohortCensusEvidence::new` now makes two checks before a single section is decoded:
that the twelve recording terms agree, and that no two samples claim the same read group. Its
refusal names both samples and the identifier they share, in the voice the terms refusal already
used. Asserted in
[`a_cohort_refuses_two_samples_whose_read_groups_were_identified_separately`](../../../../src/ng/parameter_estimation/joint/census.rs),
whose two samples are named `a` and `b`:

> samples a and b both claim read group 0; a read group belongs to one sample, so these censuses had
> their read groups minted separately and their libraries would be fitted as one

The two refusals are one type, `CohortRefusal`, and the fit turns either into its own error —
`JointFitError::IdentityMismatch` as before, and `JointFitError::SharedReadGroup` for the new one.

**The two callers.** Both now identify their read groups once, over every alignment at once, and
walk one sample per entry of the by-sample view — which is the pattern `SampleReads::open_only_sample`
already documented for cohort tools. With that in place the guard has nothing to fire on.

---

## 4. What the guard caught in this repository's own fixtures

**Every multi-sample fixture in the cohort modules was built the shape the guard forbids** — all
samples on read group `0`. That is not a detail of the fixtures; it is the unphysical assumption
written down. One of them said so out loud, in
[`ng_joint_contamination_control.rs`](../../../../examples/ng_joint_contamination_control.rs):
*"Library `k` of every plant is read group `k` … Real read groups are unique to a plant, but sharing
them here keeps the error rate fitted from the whole panel exactly as it is today."*

All of them now give plant *s*'s library *k* an identifier of its own. **This changes what the fit is
asked to do**, and the change is the real one: because a library belongs to one plant, a cohort of
*N* samples has at least *N* libraries, and **every per-library rate is always fitted from one
sample's reads**. The precision the drawn-cohort tests used to assert against — one rate from ten
samples' reads — was a regime the caller can never be in. Those assertions now check the mean over
the libraries, which carries the same precision and is an honest quantity.

**A fixture where every library misreads at the same rate cannot tell a fit that keeps them apart
from one that pooled them**, so one was added that does not. Six libraries drawn over a 6-fold range,
8,000 positions at 12 reads:

| library | drawn | fitted |
|---|---|---|
| 0 | 0.00200 | 0.00207 |
| 1 | 0.00300 | 0.00304 |
| 2 | 0.00400 | 0.00386 |
| 3 | 0.00600 | 0.00653 |
| 4 | 0.00900 | 0.00889 |
| 5 | 0.01200 | 0.01180 |

Each within 25% of its own rate; the fitted spread **5.70-fold against the drawn 6.00**. A pooled
fit would return one number for all six and a spread of 1.

---

## 5. The end-to-end check

`ng_joint_records_walk` over the 63-accession tomato panel (`benchmarks/tomato1`, all 63 CRAMs,
1,999,404 kept ordinary positions and 4,164 STR loci over 8 Mb of analysed regions).

**Before this change the run printed one library line** — measured on this same cohort and tool on
2026-08-24, before this branch:

```
read group ReadGroupId(0): a read misreads at 0.00334 at an ordinary position and 0.0237 at a mismapped one
```

**After it prints sixty-three**, one per accession, the fit converging in 30 passes over 724 s:

```
read group ReadGroupId(0):  a read misreads at 0.00353 at an ordinary position and 0.0285 at a mismapped one
read group ReadGroupId(1):  a read misreads at 0.00185 at an ordinary position and 0.0214 at a mismapped one
read group ReadGroupId(2):  a read misreads at 0.00304 at an ordinary position and 0.0233 at a mismapped one
…
read group ReadGroupId(62): a read misreads at 0.00430 at an ordinary position and 0.0292 at a mismapped one
```

| | lowest library | highest library | spread | mean over the 63 |
|---|---|---|---|---|
| ordinary position | 0.00051 | 0.00839 | **16.45-fold** | 0.00299 |
| mismapped position | 0.0030 | 0.0398 | **13.27-fold** | 0.0209 |

**That spread is what the single line was hiding.** The merged rate of 0.00334 sits in the middle of
a range whose ends differ by a factor of sixteen — high enough to look like a plausible cohort-wide
number, and wrong by a factor of six and a half for the cleanest library and by two and a half for
the dirtiest. Nothing in the merged output said so.

The cohort itself was accepted without a refusal, which is the other half of the check: all 63
samples agreed on the twelve recording terms and no two claimed one library, so the guard passed
silently on a correctly assembled run. Whole run 2,308 s: 63 walks, the ordinary-position fit, the
contamination fit and the repeat-tract fit.

**The repeat-tract half is unaffected and stays pooled on purpose.** It reports *"63 read groups in
1 slippage group, pooled"* — slippage is fitted over one group by design, so splitting the
identifiers changes nothing there.

---

## 6. Gates

- `cargo fmt --check`: clean.
- `cargo clippy --lib --tests --all-features -- -D warnings`: clean.
- `cargo test --lib`: 4,190 passed, 0 failed, 14 ignored — three more than the 4,187 on `main`, which
  are the three tests added here.
- `cargo doc --no-deps`: 23 unresolved links and 12 redundant explicit link targets, the same counts
  as `main`.
- `cargo clippy --example` on the six examples touched: four clean. `dhat_ng_joint_fit` (12 errors)
  and `ng_joint_duplicated_in_fit` (2) fail, and fail identically on `main` with these changes
  reverted — unused imports behind a feature gate, and a `needless_range_loop`, none of them on a
  line this work touched.
