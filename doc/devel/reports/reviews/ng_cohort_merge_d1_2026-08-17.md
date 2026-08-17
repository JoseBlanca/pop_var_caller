# Code Review: ng cohort merge — D1, the observation cache
**Date:** 2026-08-17
**Reviewer:** rust-code-review skill (orchestrator), seven category sub-agents in isolated worktrees
**Scope:** the uncommitted D1 diff — `src/ng/run/cohort_merge/organise.rs` (new) and one line in `mod.rs`
**Status:** Approve-with-changes

---

### 1. Scope

- **Reviewed:** the working-tree diff of step D1, as the stash commit `15ade0f1` over `3732cbff`.
- **In scope:** [organise.rs](../../../../src/ng/run/cohort_merge/organise.rs) (new, code and
  tests), [mod.rs](../../../../src/ng/run/cohort_merge/mod.rs) (the added declaration and the
  header paragraph), and this step's implementation report.
- **Out of scope:** `build.rs`, `close.rs`, `serial.rs` (milestones A–C, committed), everything
  outside `src/ng/run/cohort_merge/`.
- **Categories dispatched:** reliability, errors, naming, idiomatic, refactor_safety, smells,
  and module_structure + defaults together (one file, both checklists). `unsafe_concurrency`
  was skipped — the file has no `unsafe`, no threading primitives and no shared state yet —
  and `tooling` was skipped, nothing in `Cargo.toml` changed.
- **Audit trail:** the seven per-category files are in `tmp/review_2026-08-17_ng-cohort-merge-d1/`.

### 2. Verdict

**Approve-with-changes.** The window `cover` builds was proved *sufficient* — a reviewer drove
a real `build_region` through the cache over 600 random cohort layouts and compared every
region against the same builder handed the whole stretch, with **no disagreement anywhere**.
The defects are all in what the tests could see: the file's one real algorithm, the
sample-sweeping fixpoint, was pinned by nothing, and four claims in the implementation report
about its own fixtures were wrong.

### 3. Execution status

Run by the orchestrator in the container, and reproduced independently by each sub-agent in its
own worktree:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `test result: ok. 131 passed; 0 failed; 0 ignored;
  0 measured; 3634 filtered out`.
- `cargo test --lib` — `test result: ok. 3754 passed; 0 failed; 11 ignored` (660 s).
- `cargo clippy --all-targets --all-features` — **not run against this diff**; it is red on this
  branch for 49 pre-existing reasons in `examples/`, `benches/` and other modules' test code,
  and was red identically before this work (standing item).

Findings labelled "Needs verification": none. Every finding below rests on a mutation that was
run, with the two outputs quoted in the category file.

**Mutation totals across the fan-out:** 30 mutations run by four categories; 7 survived; 2 of
the survivors changed no behaviour on legal input.

### 4. Open questions and assumptions

1. **Does the `ObservationSource` that lands later permit a source to be polled after it
   yields `Err`?** `cover`'s documented retry depends on it and `Iterator` grants nothing.
   Affects B2. Resolved in this step by *stating the requirement* rather than assuming it.
2. **Is the sweep worst case reachable on real data?** A chain of overlapping observations
   running through the cohort in decreasing sample order costs `samples` sweeps — 28 ms for one
   11-base region at 3,000 samples. Affects M6. Documented, not optimised.

### 5. Top 3 priorities

1. **B1** — nothing pinned the fixpoint: a cover capped at two sweeps passed all 131 tests and
   disagreed with a whole-stretch builder on 410 of 600 layouts.
2. **B2** — nothing pinned what happens after a source fails: latching the sample as spent on
   `Err` passed all 131 tests and would make one sample in k silently uncovered for the rest of
   the run.
3. **M1** — the test *named* for the fixpoint could not fail for it, and its own docstring said
   the opposite.

### 6. Findings

#### Blockers

**B1: organise.rs — the fixpoint is pinned by nothing; a two-sweep cap passes every test.**
**Categories:** reliability. **Confidence:** High.
Replacing `loop { … }` with `for _ in 0..2 { … }` passes all 131 tests. It is not a no-op: a
three-sweep chain (sample 2 widens the reach to 70, sample 1 then to 120, sample 0 only then
reachable and widening to 300) gives a window one observation short, and run through
`build_region` **410 of 600 seeded layouts disagree with the whole-stretch oracle** — the first
turning a locus the oracle refuses into one the cache-fed builder emits. Every existing fixture
reached its fixpoint in two sweeps.
*Fix:* `a_chain_that_needs_a_third_sweep_is_drawn_whole`, plus the differential test of M3.

**B2: organise.rs — nothing pins the reader's state after a failure; a one-line change makes a
sample vanish silently.** **Categories:** errors. **Confidence:** High.
`pull` deliberately did not latch `spent` on `Err`, and `cover`'s doc rested on that ("a run
that retries must not have lost it"). No test made a second cover after a failure. Latching on
`Err` too passes all 131 tests; a probe fails under it and passes without. A sample wrongly
retired keeps being handed out as an empty slice, so every later cohort locus records it as
uncovered — a legitimate-looking fact, no panic, wrong genotypes for one of k samples across
the rest of the genome.
*Fix:* the test, **and** state in the doc the requirement the retry places on a source, which
`Iterator` does not grant.

#### Majors

**M1: organise.rs — the test named for the fixpoint cannot fail for it, and its docstring says
the opposite.** **Categories:** reliability, refactor_safety (convergent). **Confidence:** High.
`the_frontier_follows_a_widening_in_a_later_sample` claimed "a cover that swept the samples once
would leave it out". It would not: the far sample's one record is drawn during the first sweep —
the draw is what discovers it is beyond the reach — and then held, so it is in the window either
way. Under the single-sweep mutation that test stays green; the only failure is the *boundary*
test. It is also the test the implementation report credited with killing two mutations it does
not kill.
*Fix:* give the far sample a **second** record, so a single sweep is observable, and correct the
docstring.

**M2: organise.rs — `with_observations` takes a `GenomeRegion` and reads two of its three
fields.** **Categories:** refactor_safety, smells, module_structure/defaults, idiomatic
(convergent, four categories). **Confidence:** High.
Replacing `span.end` with `Position(0)` at every call site left the whole suite green: the field
reached no behaviour. Two consequences — the milestone-E organiser's span errors would be
accepted in silence, and `cover` and `with_observations` disagreed about an inverted region
(`cover` takes `end.max(start)`, the trim took `start` raw), which a reviewer's probe confirmed
hands back an **empty** window over ground `cover` had just drawn.
*Fix as applied:* keep the region and **make the end load-bearing** — it is checked against the
ground `cover` reached (M3), and both ends are read with the same `min`/`max` defence.

**M3: organise.rs — the cache records nothing about how far it is covered, so a failed cover and
a successful one are indistinguishable to a builder.** **Categories:** errors, reliability
(convergent). **Confidence:** High on the fact, Medium on impact (the organiser that would
mishandle it is milestone E).
`cover` returns `Result` and leaves no mark. A caller that catches the error and carries on — or
one that simply asks for a wider region than it covered — gets a short window and closes a locus
over ground the reader never reached. This is the exact failure the file exists to exclude; every
other path is defended structurally.
*Fix:* a `covered_to` field, set only by a successful cover, and a release check in
`with_observations`.

**M4: organise.rs — `evict_before` was never tested with more than one sample.** **Categories:**
reliability. **Confidence:** High.
Both eviction tests used a one-sample cache, so `for sample in self.samples.iter_mut().take(1)`
passes all 131 tests. Output is identical — at 63 samples that is 62/63 of the cache never
released, and only a heap profile would notice.

**M5: organise.rs — the coordinate-order check's contig half was pinned in one direction only.**
**Categories:** errors. **Confidence:** High.
Weakening the check to compare only *within* a contig passes: no fixture had a source going back
a contig while positions rise, which is the multi-file / psp merge case.

**M6: organise.rs — the quoted cost of the sweep loop omits the term that grows with the
cohort.** **Categories:** smells. **Confidence:** High.
`draw_to`'s doc priced the loop as "one comparison per held observation per sweep, over a window
short by construction". Two costs were missing. Each sweep visits **every** sample; and the sweep
count is unbounded — measured, a chain running backwards through the cohort costs 11.8 µs for one
11-base cover at 63 samples against **28,430 µs at 3,000**. Separately, "short by construction"
names a construction that does not exist yet: with eviction lagging 50 regions (≈200 held per
sample) the same walk at 1,000 samples costs 1,028 µs a cover against 616 µs when eviction keeps
pace.

**M7: organise.rs — the ordering assertion belongs on the "must become a `RunError`" list and
did not say so.** **Categories:** errors. **Confidence:** High.
Coordinate order is a *producer's* guarantee, which is the test `build.rs` uses to sort its three
assertions onto that list; each of those carries the migration sentence and this one did not.
`[profile.release]` sets `panic = "abort"`, so on a corrupt psp this aborts the process rather
than returning an error — and this file is the **first** such check the psp path will reach.

#### Minors

**Mi1 — `frontier` is already this module's word for something else.** (naming) The design
documents use *frontier* for how far output has been **released** (arch §4 "evicted behind the
released frontier"; spec §6.4). The growing quantity here is what `close.rs` calls the **reach**.
Both meanings would be live in this one file once the organiser lands.

**Mi2 — one concept, two names, in two files of one module.** (naming, module_structure,
refactor_safety, convergent) `start_of` here is character-for-character `close.rs`'s private
`key_of`, and `reach_of` is the same shape one field over.

**Mi3 — three verbs for taking the next observation from a source** (draw / pull / poll), with
no distinction between them; the assertions read `pulled.get()` in tests named for *drawing*.
(naming)

**Mi4 — `held` is half a name**, and `Organiser` in arch §4 already uses `held` for *outcomes*.
`first_reaching` promises an observation and returns an index; `dropped` names a count of things
not yet dropped. (naming)

**Mi5 — the header cites two measurements by internal label, one without its cohort size.**
(naming) *The C1 review* and *the C2 review* have no paths, and the 184 ms figure was measured at
**one sample** while the 3.3 µs beside it is at 63 — tied together by a "which is why".

**Mi6 — `spent: bool` reimplements `std::iter::Fuse`.** (idiomatic) The field's own doc is the
one-sentence description of `Fuse`.

**Mi7 — the "stragglers stay" rule is untestable on legal input and did not say so.**
(reliability, refactor_safety) Within one sample records are disjoint and ascending, so reach is
monotone and no fixture can separate a prefix drain from a filter. One reviewer proposed a
fixture; the other proved that fixture would be illegal input.

**Mi8 — `mod.rs`'s "what has landed" paragraph is two milestones out of date.**
(module_structure, errors) Neither `serial` — the module's own oracle — nor `organise` appeared.

**Mi9 — no eviction test crosses a contig**, and no test covers a source that is empty from the
first poll. (reliability)

**Mi10 — `expect` in non-test code without the repo's `// PANIC-FREE:` comment**, where a
restructure removes it outright. (errors, idiomatic)

**Mi11 — the reason the cache lives in the organiser's file is not recorded**, so a later
reviewer seeing an 800-line file named for a job it half does is likely to split it — which would
convert a guarantee privacy can enforce into a comment. (module_structure)

#### Nits

`spent` → a predicate name; `grew` / `folded` / `from` are each one word short of what they
hold; *fold* invents a second verb for `close.rs`'s "extends the reach"; `k` is used undefined;
`the_window_still_works_after_an_eviction` names no claim; the test helper `at` should parallel
`region_on`; "this is the module's dominant memory" asserts a size the spec supplies (320 bases
at 16 builders); the rustdoc link to a design document is a filesystem path and renders broken.

### 7. Out of scope observations

- **Four copies of the `region` / `region_on` test helpers**, one per file of this module, plus
  four near-identical observation fixtures. Past the three-copy line; the fix is a `#[cfg(test)]`
  fixture module under `cohort_merge/`, and it touches three files this step did not.
- **The concurrency shape milestone E needs an answer for:** `with_observations(&self)` against
  `cover(&mut self)` means the borrow checker forbids covering while any builder holds a window.
  That is the right refusal for one thread; several builders reading while the organiser advances
  the readers needs a design (a per-round split, an `RwLock`, or owned windows).
- **`GenomeRegion` has public fields and no `start <= end` invariant**, which is why three
  separate places in this module carry a `min`/`max` defence.

### 8. Missing tests to add now

All were added; see the fixes-applied report. In brief:
`a_chain_that_needs_a_third_sweep_is_drawn_whole`;
`a_builder_fed_from_the_cache_closes_the_loci_a_whole_stretch_would` (the differential);
`a_cover_can_be_made_again_after_a_failure`;
`a_source_that_goes_back_a_contig_is_refused`;
`eviction_drops_from_every_sample_not_only_the_first`;
`eviction_at_a_later_contig_drops_the_previous_contigs_window`;
`a_survivor_of_an_eviction_still_widens_the_next_reach`;
`two_evicted_at_once_leave_the_window_sound`;
`covering_a_region_with_an_empty_source_hands_out_nothing`;
`a_window_over_ground_no_cover_reached_is_refused`;
`a_window_reaching_past_the_covered_ground_is_refused`.

### 8a. The diff's own quantitative claims

Re-derived by the reliability agent against the implementation report:

| claim | verdict |
|---|---|
| "Tests — 16" | CHECKED-CORRECT |
| "twelve mutations, all twelve fail at least one test" | headline holds for the four re-run; **the table under it had eleven rows** |
| "single sweep → killed by the widening test" | **WRONG** — that test passes under the mutation; the boundary test is the killer |
| "never grow the frontier → the widening test" | **WRONG** — same |
| "the widening fixture's sample order is chosen so a single sweep fails" | **WRONG** — the overshoot keeps the far record either way |
| "the boundary case needed three samples to be visible" | **WRONG** — two suffice, if the widened sample carries a second record |
| "the overshoot is what a boundary-stopping cache loses" (the mechanism) | CHECKED-CORRECT |
| "`drawing_stops…` pins the count exactly" | CHECKED-CORRECT |
| "`cargo test --lib` green" | CHECKED-CORRECT |
| §2's timing figures (3.3 µs/base, 5.4 ms vs 184 ms) | quoted from the C1 and C2 reviews, not re-measured; the naming agent found the second is a **one-sample** figure quoted beside a 63-sample one |

**Five wrong claims, every one the author's own about the author's own fixtures**, against
figures quoted from other documents that were all correct — the pattern the plan-driven skill
names.

### 9. What's good

- **Every struct literal spells every field** — `SampleWindow`, and the test fixtures' six-field
  `SampleLocusObservations` and ten-field `SequenceObservation` — so a field added to any of the
  three is a compile error here (refactor_safety).
- **The closure-passing accessor enforces its own contract structurally**: the elided lifetimes
  in `f: impl FnOnce(&[&[…]]) -> R` are higher-ranked, so a builder *cannot* return a slice it
  was handed. "A builder may not hold observations" is the compiler's rule, not a convention
  (idiomatic, with the error quoted).
- **The generic-over-`E` shape is the right one**, and its two alternatives were tried: a second
  struct parameter needs `PhantomData`, and a local trait would be a second `ObservationSource`
  to delete when the run's own lands (idiomatic).
- **`draw_to` keeping no resume mark is load-bearing, not laziness** — the reviewer built the
  mark and it goes out of bounds on a window that lost two entries at once (refactor_safety).
- **The prefix-drain reasoning survived being attacked from both directions**: one agent asked
  for a fixture, the other proved the fixture would be illegal input, and the second is right
  (reliability).

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh bash tmp/mutate_d1.sh      # the mutation battery, 16 mutations
```
