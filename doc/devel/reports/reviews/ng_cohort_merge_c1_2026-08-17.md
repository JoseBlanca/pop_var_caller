# Code review — ng cohort merge, C1: one builder's job over one region

*2026-08-17. Three category checklists, each in its own isolated git worktree detached at the
commit under review. Per-category audit trail in
`tmp/review_2026-08-17_ng-cohort-merge-c1/`.*

## 1. Scope

- **What:** the working-tree diff of step C1, as commit `39ebb05e` (parent: the
  sparse-support-row commit, branch `ng-cohort-merge`). One file,
  `src/ng/run/cohort_merge/build.rs`, +322/−1: `RegionOutcome`, `build_region`, seven tests.
- **Out of scope:** `close.rs`, `mod.rs`; frozen production; the cache and the organiser.
- **Categories:** reliability (the mutation pass), extras (intent, ownership beyond this
  function, hot path), smells + naming.

## 2. Verdict

**Approve with changes** — 3 Blockers, 2 Majors, and a naming pass. Every Blocker was a test
that could not fail; the code was right in all three.

## 3. Execution status

- `cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 104 passed at review time.
- Mutation pass: **15 run, 6 survived, 3 of those changed no behaviour** — so three real
  survivors, each proved by a probe that passes on the code and fails on the mutant.

## 4. The three Blockers

**B1 — no locus opened on a region's first base.** Every fixture's loci opened strictly
inside their region, so widening `start < region.start` to `<=` passed all 104 tests — and
that mutation makes *every* builder skip such a locus: the earlier one breaks past it, the
owner skips it, and it leaves the run with nothing panicking and no span in the failed list.
Building regions are dealt out end to end at twenty bases, so it is about **one locus in
twenty** that opens on one. Two categories found this independently.

**B2 — the contig terms in the break were untested.** The other-contig fixture put both loci
at position 12, so positions never compared in a way that could trip it, while its own doc
claimed coverage "however its positions compare". Reducing the break to
`locus.region.start > region.end` truncated a whole region's output — proved with a locus on
the earlier contig at position 80 against a region ending at 50.

**B3 — `min_alt_obs` was never passed at a value of its own**, so hardcoding the default
inside `build_region` survived. Its sibling `max_cohort_locus_span` was covered twice at a
non-default value; the threshold that decides whether a locus is built at all was not.

All three are closed, and each mutation was re-run against the new test and now fails.

## 5. The two Majors, both about claims rather than code

**M1 — the failed spans' stated purpose does not hold under this function's own input
contract.** The architecture says a failed locus's span is needed so the organiser can
displace a locus a neighbouring builder built inside that ground. But every builder is handed
everything overlapping its region, so a later locus overlapping an earlier one would have
chained into it and been skipped as the earlier region's: **two loci owned by different
regions cannot overlap at all** under this contract. What remains load-bearing is spec §3.3's
count. The doc now says that, and says the displacement rule is the organiser's to keep or
retire — **which is a finding E2 should have before it builds its overlap resolution**, since
the fixture the plan prescribes for E2 is the very one that demonstrates the exclusion.

**M2 — the cost of handing a builder more than its region, measured.** The walk starts at the
beginning of the slices it is given, so every locus opening before the region is closed and
then discarded: **3.3 µs per prefix base at 63 samples, 40 µs at 250** — 63 µs with no prefix
against 16.7 ms with five thousand bases of it. Handing every builder the whole stretch is
therefore quadratic in its length: about **23 hours for a megabase at 63 samples**, against
seconds when each builder gets a window. Recorded at `build_region`, and it is what the
observation cache (milestone D) exists to fix; the serial driver avoids it by handing over
whole analysed regions, where the prefix is empty.

## 6. Naming, applied

`RegionOutcome::observations` → `cohort_observations` and `failed` → `failed_locus_spans`
(a `Vec<GenomeRegion>` spelled exactly as a boolean would be, beside a `Verdict::Failed`), and
`build_region`'s `region` parameter → `builder_region`, because *region* was carrying three
meanings in one file and two of them met in one condition. The type keeps the architecture's
name; the field renames are a recorded departure from arch §4's spelling.

Three loose words corrected: an 8-base record called a *chain*, a sample said to cover "the
tail" of a locus where it covers the eighth base of thirteen, and `RegionOutcome`'s doc
crediting the organiser with summing counts the struct does not hold.

## 7. What's good

- The partition test earns its place: it killed 6 of the 15 mutants on its own, including
  both wrong ownership rules — a locus claimed by two regions and a locus claimed by none —
  because its deletion opens on the **last** base of one ten-base piece and ends in the next.
- Three of the six surviving mutations were shown to change no behaviour rather than being
  reported as findings: `break` → `continue` and dropping either disjunct return byte-identical
  outcomes, because the walk's yield order makes the break condition monotone.

## 8. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
```
</content>
