# Code review — ng cohort merge, C2: the serial driver and the oracle

*2026-08-17. Two category passes (reliability; extras + smells), each in its own isolated git
worktree detached at the commit under review. Audit trail in
`tmp/review_2026-08-17_ng-cohort-merge-c2/`.*

## 1. Scope

- **What:** the working-tree diff of step C2, as commit `ef3fed1e` (parent: the C1 commit
  `502086e1`). The new `src/ng/run/cohort_merge/serial.rs` (+467) and one line in `mod.rs`.
- **Out of scope:** `build.rs`, `close.rs` (earlier commits); frozen production; the cache and
  the organiser.

## 2. Verdict

**Approve with changes** — no Blockers, 5 Majors, and the two reviewers converged on four of
them. Mutation pass: **8 run, 1 survived, 0 changed no behaviour**.

## 3. Execution status

- `cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 111 passed at review time.

## 4. The Majors

**M1 — the end-to-end fixture's quality assertion could not fail, and both reviewers proved it
the same way.** Every one of the substitution sample's thirty minted records carried the
identical `q_sum` of −20.723265836946418, because every base was quality 30 and `q_sum` is
Σ ln P(base-call error) whether or not the base matched. So "the weakest of the six means"
equalled the strongest, the mean, and any one of them: flipping `max` to `min` in the
composition — the exact divergence from production the code calls deliberate — left the test
green. Fixed by giving one read a bad base (quality 5) inside the locus; the mutation now
fails it, `left: -20.72 … right: -6.91`.

**M2 — my own doc had both fixture positions wrong by five bases.** It said a substitution at
107 and a deletion over 105–109 where the reads produce 112 and 110–114 — consistent with each
other, and contradicted by the test's own assertion eight lines below, which is why nothing
caught them. The prose is now the data's.

**M3 — no partition-invariance test, in the file whose header calls itself the oracle.** Spec
§15 names that property this component's regression anchor and the plan's C2 line cites it as
its source. Both reviewers checked that it holds — one at region widths 1/3/7/20, the other at
1/20/200/1000 — so it was cheap to pin, and it is now pinned at 1, 6, 60 and 120 regions over
sixty loci.

**M4 — overlapping or repeated analysed regions duplicated observations silently.** With
`[1–60, 40–100]` the locus at 45 came back twice; with the same region twice, everything did;
descending regions came out in the order given. The doc priced a violation as wrong *order*
where it is a wrong *answer* — and nothing downstream can tell one locus carried twice from a
cohort that varied at two places. Now a release-level check with two tests, and the doc says
why a bare slice is the wrong shape to trust.

**M5 — the driver is quadratic in the number of analysed regions, and its doc claimed
otherwise.** Every call closes the loci from the beginning of the observations it is given and
discards those before its own ground, so the same 20,000 observations cost **5.4 ms in one
region and 184 ms in a thousand** (release, one sample), growing linearly in region count.
`build_region`'s doc said the serial driver hands over whole analysed regions "where the prefix
is empty by construction" — true only when there is exactly one. The doc now states the
measurement and says to hand it the run's own regions, with the cache named as what makes
short ones affordable.

Two coverage findings of the same shape as C1's: `min_alt_obs` was never passed at a value of
its own, and the failed-spans test put both its failed loci in the first of two regions, so
stopping after one region or walking them backwards left it green. Both closed, and the
stop-after-one mutation now fails four tests.

## 5. What held up

- **Sample 0 is the sample the prose calls sample 0**, the substitution sample really is the
  six-record case inside the locus, and the deletion sample really has one — checked by
  dumping every minted record.
- Swapping the two samples, moving the substitution outside the deletion's span, and splitting
  the deletion across two placements are all caught by the existing assertions.
- Dropping the failed-span gathering, ignoring the region bounds, and building one region over
  everything are all caught.

## 6. Cross-category, and it is the fixture's most valuable finding

**The generic mint writes a record at every covered position** — thirty records for a
thirty-base read, twenty-nine of them an all-`A` observation against an all-`A` reference. So
`alleles_of_sample`'s explanation of the gap between a sample's records ("where this sample
minted nothing, because none of its reads departed from the reference there") names the wrong
cause: on that path a gap is ground **no read covered**. Stronger, the case is unreachable
there — a read named on both sides of a gap spans it, so a record exists. The doc is corrected
in this commit, with the caveat the reviewer attached: one fixture under the default generator
configuration, and whether a depth cap can suppress a record at a covered position was not
audited.

## 7. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh cargo test --lib
```
</content>
