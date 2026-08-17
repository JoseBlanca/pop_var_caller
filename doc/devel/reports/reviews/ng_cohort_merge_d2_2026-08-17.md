# Code Review: ng cohort merge — D2, the merge read through the cache
**Date:** 2026-08-17
**Reviewer:** rust-code-review skill (orchestrator), three category sub-agents in isolated worktrees
**Scope:** the uncommitted D2 diff — a second driver in `serial.rs`, plus one accessor on the D1 cache
**Status:** Approve-with-changes

---

### 1. Scope

- **Reviewed:** the working-tree diff of step D2, as the stash commit `1af2c517` over `284b680b`.
- **In scope:** everything D2 adds to
  [serial.rs](../../../../src/ng/run/cohort_merge/serial.rs) —
  `merge_cohort_through_cache`, `building_regions_of`, the extracted region guard and the new
  tests — and `held_observations` on
  [organise.rs](../../../../src/ng/run/cohort_merge/organise.rs), plus this step's report.
- **Out of scope:** `merge_cohort_serially` and its tests (the oracle, committed at C2 and
  deliberately untouched), `build.rs`, `close.rs`, the rest of `organise.rs` (committed at D1).
- **Categories dispatched:** reliability alone (the step's whole claim is a correctness one),
  errors + idiomatic together, and naming + module_structure + smells together. Three agents
  rather than seven: the diff is one function, one helper and one accessor.
- **Audit trail:** `tmp/review_2026-08-17_ng-cohort-merge-d2/`.

### 2. Verdict

**Approve-with-changes.** Byte-identity holds well past the fixtures — a reviewer's randomised
differential over **600 layouts found no disagreement** — but the claim itself had a hole no
fixture covered, and the two properties the step was built to make testable were both pinned
more weakly than the report said.

### 3. Execution status

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `154 passed; 0 failed` at review time.
- `cargo test --lib` — `3777 passed; 0 failed; 11 ignored` (593 s).
- `cargo clippy --all-targets` — not run; red on this branch for 49 pre-existing reasons
  (standing item).

Findings labelled "Needs verification": none.

**Mutation totals across the fan-out:** 28 mutations run; 6 survived; 1 of those changed no
behaviour on legal input.

### 4. Open questions and assumptions

1. **May a source be polled after it yields `Err`?** The driver's documented retry depends on
   it; `Iterator` grants nothing. Affects M5 — resolved by stating the requirement and what a
   source that cannot honour it must do instead.
2. **Do malformed analysed regions reach this signature?** The oracle's own doc says yes ("a
   user-supplied BED is exactly the shape that arrives overlapping"), which is what makes B-class
   finding M1 worth fixing rather than assuming away.

### 5. Top 3 priorities

1. **M1** — an inverted analysed region breaks the byte-identity claim outright: the oracle
   builds nothing, the cached driver builds everything, and with a second region the same locus
   comes back twice.
2. **M2** — the test that pins "the ground is divided" tests the divider, not the driver;
   mutating the driver's call site left it green.
3. **M3** — `cover(the whole analysed region)` — the cache's entire purpose undone — survived
   all 154 tests.

### 6. Findings

#### Majors

**M1: serial.rs — an analysed region with its ends inverted walks past the guard, and the two
drivers then disagree.** **Categories:** reliability, errors (convergent). **Confidence:** High,
demonstrated on unmutated code.
`building_regions_of` orders the two ends; `build_region` compares them raw. On a single
analysed region `50-1` over records at 12 and 45 the oracle emits **nothing** and the cached
driver emits **both loci**. Followed by a second region, the guard's own comparison
`(earlier.contig, earlier.end) < (later.contig, later.start)` reads `50-1` as ending at 1 and
accepts the pair — and the locus at 45 comes back **twice**, which is the corruption that guard's
panic message says it exists to prevent.
*Fix:* refuse a region whose ends are inverted, in the guard both drivers share.

**M2: serial.rs — the division test tests the divider, and nothing said the driver used it.**
**Categories:** reliability. **Confidence:** High.
The brief's first "right answer, wrong machine" property is *that the ground is divided at all*.
Mutating the driver's call site to `std::iter::once(*analysed_region)` — the shape the defect
takes — left the tile test green; the only failure was the eviction test, and only because that
fixture happens to be 600 bases at 20-base regions. At width 5 or 13 the same fixture ends
holding nothing, so a plausible fixture change would remove the only pin on *both* properties at
once.
*Fix:* assert that the width changes what the cache holds — 60 held at one region for the
stretch, 2 at twenty-base regions.

**M3: serial.rs — covering the whole analysed region instead of each building region survives
every test.** **Categories:** errors. **Confidence:** High, behaviour proven.
The existing eviction test reads the window only *after* the merge, where the final eviction
hides the difference. Under the mutation the first cover draws 30 observations where the shipped
code draws 2.
*Fix:* read the window at the one moment the driver is caught mid-stretch — a source failure.

**M4: serial.rs — the overlap guard's decisive comparison has no fixture on it.**
**Categories:** reliability. **Confidence:** High.
Both refusal tests overlap by 21 bases. Relaxing `<` to `<=` keeps the module green and makes
*both* drivers emit a locus on a shared base twice — so neither the byte-identity fixtures nor
the randomised differential can see it.
*Fix:* two regions sharing exactly one base.

**M5: serial.rs — a failed merge advances the cache, and the obvious retry comes back short with
`Ok`.** **Categories:** errors. **Confidence:** High, measured.
On a source failure the driver returns `Err` and drops what it built, but the cache keeps the
readers' position and the eviction that has already happened. A caller that retries the same
ground with the same cache is not refused: measured on a 60-locus fixture failing after the
thirtieth, **32 of 60 loci came back and the return value was `Ok`**.
*Fix:* say so in the doc now; making it unrepresentable belongs with the organiser, which owns
the cache.

**M6: organise.rs — `held_observations` is exercised by one single-sample test.**
**Categories:** reliability. **Confidence:** High.
Replacing the sum with the first sample's length passes everything. A memory report that counts
one sample of a thousand under-reports the cache by three orders of magnitude, and the assertion
the accessor exists to support would stay green if eviction reached only the first sample.

**M7: serial.rs — a source out of coordinate order aborts the process, though the signature
promises `Result`.** **Categories:** errors. **Confidence:** High that it is reachable; Medium on
the remedy. This driver is the first thing to take observations from an arbitrary source, and
`[profile.release]` sets `panic = "abort"`.
*Fix now:* name it on the driver, so the promise the signature makes is accurate. The durable
fix is the `RunError` migration `organise.rs` already records that it owes.

#### Minors

**Mi1 — the eviction test's stated mechanism is not what produces its number.** (reliability,
smells, convergent) "The last region's own record and the one draw past it" — measured, both
held records lie *inside* the last building region and the source is spent, so nothing was drawn
past the ground; at five-base regions the same fixture ends holding none.

**Mi2 — comparing whole `Debug` renderings gives a 36 kB failure with no pointer at the
difference.** (smells) Entry by entry gives 888 bytes naming the first entry that differs. The
sibling file already does it that way.

**Mi3 — `both_drivers_agree` reads as a claim, asserts nothing, and returns three values.**
(naming, idiomatic, convergent) Five call sites each write the comparison themselves; one that
forgot would compile, run both drivers and pass.

**Mi4 — `building_regions_of` is the organiser's geometry sitting in the oracle's file.**
(module_structure) Its second caller is the milestone-E organiser, which hands regions out; left
here, that organiser either imports from the reference implementation's file or re-derives the
clamp.

**Mi5 — the file header describes one driver and calls the file "the oracle"; it now holds
two**, and describes the windowed view in the future tense in a file that contains it.
(module_structure)

**Mi6 — "building region" carries the step's whole argument and is defined only in a sibling
file's header.** (naming) The clamp argument and the eviction argument are both statements about
building regions.

**Mi7 — the doc gives what the cache saves and never what it costs.** (naming) The saving is
quantified (3.3 µs per prefix base at 63 samples); the price — one cover per building region,
`sweeps × (samples + held)` — is measured ten lines away in another file and not mentioned.

**Mi8 — `held_observations` returns a count and is named for the observations.** (naming) The
crate's convention is `_len`, and this accessor's name will spread to the run's memory report.

**Mi9 — the coordinate-ceiling test collects an unbounded iterator**, so the regression it
guards against is a SIGKILL of the test binary rather than a diff. (idiomatic, measured)

**Mi10 — `SourceFailed` is now defined twice, identically**, beside four existing copies of
`region`/`region_on` across the module. (module_structure)

**Mi11 — the doc's "eviction happens before each cover" has no regression test**; swapping the
two lines leaves the module green. It is a cost, not a correctness property: measured, the
pre-cover eviction takes the window from three records to one at every region.

**Mi12 — the report's "14 new tests" is 12**, its byte-identity table names one test that never
runs the oracle, and one mutation row is true only of a narrower mutation than it describes.

#### Nits

The width has three names in one step; `refuse_overlapping_analysed_regions` names one of the
two things it enforces; `source_of` against the sibling's `source_over`; one test builds its
fixture twice; the `width - 1` is the one arithmetic operator on its line without a comment
saying why it is safe.

### 7. Out of scope observations

- **`with_observations`'s left trim is inert in this caller** — the driver evicts at exactly the
  same base immediately before, proven a no-op over every fixture and 600 random layouts. Two
  mechanisms doing one job; milestone E's organiser should know which it is relying on.
- **Gathering one region's outcome into the accumulator is written out twice**, once per driver,
  in two functions whose whole relationship is that they must produce the same bytes. A method on
  `RegionOutcome` would make that structural; it lives in `build.rs`.
- **The four copies of `region`/`region_on`** in `build.rs` and `close.rs` remain.

### 8. Missing tests to add now

All added; see the fixes-applied report. In brief: the inverted analysed region; two regions
sharing one base; the width-changes-the-window pair; the window at the moment of failure; the
held count across two samples and over none; no analysed regions; no samples; the division's
reading of an inverted region; and the randomised differential between the two drivers.

### 8a. The diff's own quantitative claims

| claim | verdict |
|---|---|
| "14 new tests" | **WRONG** — 12 by `git diff \| grep -c '#\[test\]'` |
| "every one of them can fail" | **not what the evidence shows** — nine driver mutations, each killing at least one test, is a different claim |
| "five fixtures assert byte-identity" | **count correct, one row wrong** — the deletion-boundary test never runs the oracle |
| "the window holds 2 of 60" | **number correct, mechanism WRONG** — both are the last region's own records; at five-base regions the count is 0 |
| "without the eviction call it holds all sixty" | CHECKED-CORRECT |
| "the width ignored kills 2 tests" | **true of a mutation inside the divider, false at the driver's call site** |
| "one fixture defect — `build_region` refuses two members disagreeing on the reference" | CHECKED-CORRECT, with the panic and its pinning test cited |
| "`154 passed; 0 failed`" | CHECKED-CORRECT |

**Five wrong claims, every one the author's own about the author's own fixtures.** Same pattern
as D1.

### 9. What's good

- **The randomised differential came back clean at 600 layouts**, with 600 of 600 containing a
  record that straddles an analysed region's right edge — the regime the fixtures cannot reach.
- **`building_regions_of` as `impl Iterator` from `std::iter::successors` allocates nothing**,
  and the two alternatives were tried: `from_fn` needs mutable state to say the same thing, and
  a stepped range can express neither the clamp nor `u64::MAX`.
- **The successor is computed from the *clamped* end**, which is what makes the division stop at
  the analysed region's last base rather than one width past it.
- **Four arithmetic decisions were each mutated and each killed** by an existing test.
- **`with_observations`'s covered-ground check is unreachable from this driver** — every call is
  preceded by a cover on the same region — which is the right relationship between a guard and
  its caller.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh bash tmp/mutate_d2.sh      # the mutation battery, 13 mutations
```
