# ng — read filtering in stages, C1b: `container.rs` gets a test module

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `c718a1c` (C1)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) — an
**owner-added step** between C1 and C2, not one of the plan's original four.
**Origin:** the deferred findings of
[`fixes_applied_2026-08-03_v5.md`](../reviews/fixes_applied_2026-08-03_v5.md) §5 (B2's review).

---

## 1. Plan

`src/ng/read/input/aligned_reads_reader/container.rs` had **no test module at all**, and B2's
review found four input classes on `DecodedContainer::fill_raw_read` that nothing reached:

1. a record with **no name** — `PackedReadEntry::name` is `None`, and `None` is not an empty name;
2. **empty** sequence / quality / CIGAR spans;
3. the **clear-and-refill** claim the function's own doc makes — no test served a long read and
   then a short one through one buffer;
4. `out.data_mut().clear()` — whose doc names its own silent failure mode, and **whose deletion
   left the whole suite green**.

B2 deferred them with a reason and a trigger: *"worth taking **before C2**, which moves the
filtering loop and is exactly the change that could hand this function a buffer with a history."*

**That reason is wrong, and C1b's review is what established it.** The buffer already has a
history: `ReadFilter::next` refills **one** `NoodlesRawAlignedRead` for a whole pass, so every read
after the first arrives carrying the previous one, and has since B2. Measured by instrumenting
`fill_raw_read` and running the CRAM cursor walk — read 0 arrives empty, every read after it
arrives with 30 bases in place; asserting freshness turns three *existing* tests red.

So nothing was waiting on C2. What was missing is not the condition but the check: **no test
anywhere compares a served read's content against an independent expectation**, which is why
deleting `sequence.clear()` grows sequences past 300,000 bases while the rest of the suite — the
CRAM-versus-BAM oracle included — still passes. A regression in these clears would corrupt
production reads *today*, silently, and only this module would catch it. §6 records the correction
against the original deferral.

**Almost test-only.** The nine-test diff was; the review then added four production doc comments —
`fill_raw_read`'s buffer-history statement and panic contract, and `push`'s note that its error
path is 4 GiB-gated. No production *behaviour* changes.

## 2. Assumptions

**The tests build a `DecodedContainer` directly rather than through `decode_container_at`**, which
is the only thing that builds one in production and which needs a real CRAM. That is deliberate,
not a shortcut: the decode is already covered by the CRAM walks in `open_bam.rs`, and the part
these findings are about is the *packing round trip* — `push` flattening a record into the two
shared buffers and `fill_raw_read` rebuilding it. Both are reachable from a test module in the
same file, so the tests exercise the same path `decode_container_at` takes, minus the file. That
keeps them fast, deterministic, and able to construct shapes a real CRAM will not conveniently
produce (an unnamed record, an empty one).

**No new production API was added to make this testable.** The private fields are reachable
because the test module is a child of the module that owns them.

## 3. Changes made

One new `#[cfg(test)] mod tests` in `container.rs`, with a module doc recording why it drives the
round trip directly, which four things had never been checked, and — after the review — that the
buffer already has a history.

Helpers: `container_of(records)` (packs through the real `push`), `full_record(name, bases)` —
every scalar given a **distinct** value, so a packing that transposed two of them cannot survive —
and `raw_read_from(container, i)`.

## 4. Tests added

| test | class | what it pins |
|---|---|---|
| `a_record_and_its_read_group_survive_the_round_trip_through_the_packed_form` | — | whole-`RecordBuf` equality, so a scalar added to `PackedReadEntry` and forgotten in either direction fails without anyone extending an assertion; plus the read group, the other half `fill_raw_read` sets |
| `an_unnamed_record_does_not_inherit_the_name_left_in_the_buffer` | 1 | `None` and an empty name are different states, and one buffer serves both records — the only arrangement in which this can fail |
| `an_empty_name_stays_empty_rather_than_becoming_absent` | 1 | the other direction of the same distinction |
| `an_empty_record_inherits_no_bases_qualities_or_cigar_from_the_one_before_it` | 2 | a zero-length byte range copies as a no-op, so the *clears* are all that stand between this record and the last one's contents |
| `a_shorter_record_keeps_no_tail_of_the_longer_one_before_it` | 3 | the clear-and-refill promise, stated as a test |
| `a_buffer_carrying_auxiliary_tags_comes_back_with_none` | 4 | the clear whose deletion left the suite green |
| `each_record_is_stamped_with_its_own_read_group` | — | asserting one group would pass on a `fill_raw_read` that ignored `i` |
| `a_container_counts_the_records_it_packed_and_charges_none_elsewhere` | — | `len` and `other_sample_records` are separate quantities; packing charges nothing to another sample |
| `records_read_their_own_slices_of_the_shared_buffers` | — | entry *i* reads its own slices of a shared buffer — invisible to a single-record fixture, and it survives a mutation serving every read as the *next* one, which the two-record tests pass |

**Four more were added by the review**, each closing a mutation that had survived the whole suite:

| test | what it pins |
|---|---|
| `every_scalar_is_read_from_the_entry_asked_for` | six of the seven scalars could be read from the wrong entry with the suite green — the fixture gave every record identical scalars, so the module pinned span per-entry-ness and not scalar per-entry-ness |
| `a_span_past_the_index_width_is_refused_rather_than_truncated` | `Span::new`'s documented refusal; truncation would hand back another record's bytes |
| `shrinking_gives_back_the_slack_the_buffers_grew_by` | `shrink_to_fit`, which nothing could observe and which looks like dead code to a reader who has not met `decode_container_at`'s comment |
| `serving_a_second_read_reuses_the_first_reads_allocations` | that reuse *happens*, not merely that it is safe — the property the packed form exists for, and behaviourally invisible |

**Suite: 2,847 → 2,860 (+13)** — 9 as built, 4 from the review. Fully accounted.

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,860 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

The dumps cannot move — the only production changes are doc comments — but they were run rather
than reasoned about, because the plan checks every step against them.

### Mutations run

Each marker `grep -c`-confirmed present before the run, the file byte-restored after. The review
ran **47** more in an isolated worktree and reproduced all four of the first group exactly, kill
set for kill set.

| mutation | killed |
|---|---|
| `None => *out.name_mut() = None` → `None => {}` | the unnamed-record test, **alone** |
| `out.data_mut().clear()` deleted | the aux-tag test, **alone** |
| `sequence.clear()` removed (append rather than replace) | the empty-record and shorter-record tests |
| `&self.index[i]` → `&self.index[0]` | 5 tests |
| **`flags` read from `self.index[0]`** | `every_scalar_is_read_from_the_entry_asked_for`, **alone** — nothing at all before the review |
| **`Span::new` truncates instead of refusing** | `a_span_past_the_index_width_…`, **alone** — nothing before |
| **`shrink_to_fit` emptied** | `shrinking_gives_back_the_slack_…`, **alone** — nothing before |
| **name replaced wholesale rather than cleared and refilled** | `serving_a_second_read_reuses_…`, **alone** — nothing before, and behaviourally invisible |

The aux-tag row is the one the step was asked for; the bottom four are ones the review found.
Before C1b, every one of these eight left the suite green.

## 6. Tradeoffs and follow-ups

- **The deferral's stated reason was wrong and is corrected at its source.**
  [`fixes_applied_2026-08-03_v5.md`](../reviews/fixes_applied_2026-08-03_v5.md) §5 now carries the
  correction: nothing was latent, and nothing was waiting on C2. That entry also named the third
  deferred item differently from §4 of the same report ("buffer-shrink" against "clear-and-refill");
  C1b covers **both**, so the deferral is discharged under either reading.
- **`push`'s error path stays unreachable** — it needs a 4 GiB payload, so building the input costs
  4 GiB of memory. The refusal it propagates is pinned where the truncation would happen instead,
  on `Span::new`, and `push`'s doc now says so rather than leaving a coverage audit to re-derive it.
- **A property test over the round trip is deferred.** `push` + `fill_raw_read` is a round-trip law
  over a structured domain, and a proptest would subsume the scalar gap automatically and reach
  orderings the fixed fixtures do not. It needs an `arb_record_buf()` generator and is its own
  piece of work; the five hand-written tests already kill every mutation it was proposed to catch.
- **The decode path stays covered where it was** — `decode_container_at` and `owner_of_cram_record`
  are exercised by the CRAM walks in `open_bam.rs`, including the multi-container and
  multi-read-group fixtures B2 added.
