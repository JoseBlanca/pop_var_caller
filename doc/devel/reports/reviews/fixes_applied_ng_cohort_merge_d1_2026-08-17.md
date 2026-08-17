# Fixes applied — ng cohort merge D1 (the observation cache)

*2026-08-17. Against [the D1 review](ng_cohort_merge_d1_2026-08-17.md); step D1 of
[the plan](../../ng/impl_plan/cohort_merge.md).*

**Everything actionable was applied.** Two Blockers, seven Majors, eleven Minors and the nits;
nothing deferred to the owner, nothing disputed. The module's tests go from 131 to 142, the
cache's own from 16 to 27.

## The two Blockers, and what now catches them

- **B1, the unpinned fixpoint.** Added `a_chain_that_needs_a_third_sweep_is_drawn_whole` —
  three samples widening each other in turn, where the third sweep is what puts the last
  observation in the window — and the review's differential test (below). A cover capped at two
  sweeps now fails both.
- **B2, the state after a failure.** Added `a_cover_can_be_made_again_after_a_failure`: a source
  that fails on its second observation and yields a third afterwards, with a second sample beside
  it so the retry has to rebuild the reach from the held window. **And the doc now states the
  requirement the retry places on a source** — that it may be polled after yielding `Err`, which
  `Iterator` does not grant — and says what a source that cannot honour it must do instead
  (yield `None`, which reads as a spent sample).

## The Majors

- **M1** — the widening test now gives its far sample a **second** observation, so a single sweep
  is observable, and its docstring says what actually pins what: one observation beyond the reach
  is held whatever the cover does, because the draw that discovers it keeps it. Under the
  single-sweep mutation this test now fails; before, it passed.
- **M2 + M3 together, and the resolution is one change.** Four categories found `span.end`
  unread; two found that nothing recorded how far the cache was covered. Both are fixed by making
  the end *do that job*: `ObservationCache` gains `covered_to`, set only by a **successful**
  cover, and `with_observations` checks the span's last base against it at release level. The
  parameter is kept rather than narrowed to a `GenomePosition`, because keeping it is what makes
  the check possible — and the inverted-region disagreement the reviewers found goes with it,
  since both ends are now read with the same `min`/`max` defence `cover` uses. Two tests:
  `a_window_over_ground_no_cover_reached_is_refused` and
  `a_window_reaching_past_the_covered_ground_is_refused` (the second was added after the first
  mutation run showed a check against the span's *start* still passing).
- **M4** — `eviction_drops_from_every_sample_not_only_the_first`, two samples each with something
  to lose and something to keep.
- **M5** — `a_source_that_goes_back_a_contig_is_refused`, positions rising while the contig falls.
- **M6** — `cover`'s doc now prices the loop as `sweeps × (samples + held)`, names the worst case
  (a chain running backwards through the cohort, one sweep per sample), and carries the review's
  five measured figures. It also says that the `held` term stays short **only while the organiser
  evicts at the pace it releases ground**, and that this module cannot enforce that.
- **M7** — the ordering assertion carries the migration sentence its three siblings in `build.rs`
  carry, and says it is the first such check the psp path will reach.

## The Minors and nits

- **Mi1** — `frontier` renamed **`chain_reach`** throughout, in the code and in the test names.
  The design documents use *frontier* for how far output has been released, and the organiser
  brings that meaning into this same file at milestone E.
- **Mi2** — `start_of` and `reach_of` are gone. `SampleLocusObservations` gains
  `start_position()` and `reach_position()` beside its existing `reach()`, and `close.rs`'s
  private `key_of` is deleted in favour of the same accessor — one definition where there were
  two, in the two files that walk by it.
- **Mi3** — one verb: everything **draws** (`draw_next`, `last_drawn`, the test counter `drawn`),
  matching the spec's own word. *Pull* and *poll* are gone.
- **Mi4** — `held` → `held_observations`, `first_reaching` → `first_reaching_index`,
  `dropped` → `first_survivor`, the local `windows` → `observations_per_sample` (which is what
  every consumer already calls it).
- **Mi5** — the header cites the two reviews by path, and says the 184 ms figure was measured on
  **one sample** rather than tying it to the 63-sample one with a "which is why". The
  implementation report carries the same correction.
- **Mi6** — the `spent` flag is gone; the source is held as `Fuse<S>`. A failure stays
  `Some(Err(_))` and so leaves the source live, which is what B2's retry needs.
- **Mi7** — the "stragglers stay" paragraph now says no legal input can separate a prefix drain
  from a filter, and why (a sample's records are disjoint and ascending, so reach is monotone).
  No fixture was added: the one proposed would have been illegal input. The mutation is recorded
  as a survivor by design.
- **Mi8** — `mod.rs`'s "what has landed" paragraph now names all five files, `serial` included.
- **Mi9** — `eviction_at_a_later_contig_drops_the_previous_contigs_window` and
  `covering_a_region_with_an_empty_source_hands_out_nothing`.
- **Mi10** — the `expect` is gone: the ordering check binds the previous position with
  `if let Some(previous)` and formats it inline.
- **Mi11** — the header records why the cache lives in the organiser's file: when the organiser
  lands, `cover` and `evict_before` can become **private to the file**, which a cache in a file
  of its own could not be.
- **Nits** — `at` → `position_on`; `the_window_still_works_after_an_eviction` →
  `a_cover_after_an_eviction_draws_forward_and_keeps_the_survivor`; `grew` → `reach_grew`;
  `folded` → `considered`; `k` replaced by "one slice reference per sample"; the memory claim
  carries its number (320 bases at 16 builders on 20-base regions); the `cover` loop is split
  into a named `sweep`, so the fixpoint is one line; the rustdoc filesystem link is now a plain
  path, as in the sibling files.

## Two things the reviewers asked for that were **not** done, with the reason

- **`held_observations()` as a public memory probe.** No consumer exists, and `-D warnings`
  would reject it; the memory question it answers is milestone E's, where the organiser will
  want it. Recorded as an open item rather than added dead.
- **A shared `#[cfg(test)]` fixture module** for the four copies of `region` / `region_on` across
  this module's files. It is a real duplication and the fix touches three files this step did not
  write. Recorded as an open item.

Two design-document items are the owner's and are recorded, not made: the arch's file tree
describes `organise.rs` as the organiser without mentioning the cache, and arch §4's
`with_observations` sketch now differs from the code in what `span` does.

## Mutation re-run

The battery was rewritten for the new code and re-run in the container: **sixteen mutations
written, one does not compile, fifteen run, fourteen killed.** The single survivor is the
prefix-drain-versus-filter one, which Mi7 establishes cannot be separated by any legal input.
Both Blocker mutations — a two-sweep cap and a spent-on-`Err` latch — now die, as do the
five that the review's new tests target.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` → `142 passed; 0 failed`; the whole library suite
re-run after the fixes → `3765 passed; 0 failed; 11 ignored` (559 s), against `3754 passed`
before them.
