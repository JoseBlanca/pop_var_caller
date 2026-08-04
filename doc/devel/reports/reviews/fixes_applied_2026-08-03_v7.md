# Fixes applied — ng read filtering in stages, C1b

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `c718a1c`
**Review:** [`ng_read_filtering_stages_c1b_2026-08-03.md`](ng_read_filtering_stages_c1b_2026-08-03.md)
**Impl report:** [`ng_read_filtering_stages_c1b_2026-08-03.md`](../implementations/ng_read_filtering_stages_c1b_2026-08-03.md)

Applied: 13 · Applied with adaptation: 0 · Already fixed: 0 · **Deferred: 1** · Disputed: 0

---

## 1. Findings table

| id | severity | subject | decision | status | validated |
|---|---|---|---|---|---|
| M1 | Major | "every caller passes a fresh buffer" is false | Apply | Applied | prose; claim re-verified against `filtering.rs` |
| M2 | Major | six of seven scalars readable from the wrong entry | Apply | Applied | Pass — new test kills all seven, alone |
| M3 | Major | `Span::new`'s overflow refusal uncovered | Apply | Applied | Pass — kills the truncation mutation, alone |
| M4 | Minor | `shrink_to_fit` unobservable and unpinned | Apply | Applied | Pass — kills the no-op mutation, alone |
| M5 | Minor | allocation reuse untested | Apply | Applied | Pass — kills the wholesale-replace mutation, alone |
| M6 | Minor | `filled` is a bare participle | Apply | Applied | renamed `raw_read_from`, 9 call sites |
| M7 | Minor | "span" collides with "reference interval" | Apply | Applied | test renamed, doc reworded |
| M8 | Minor | "short read" is a charged domain term | Apply | Applied | sibling pair renamed to read as a pair |
| M9 | Minor | "four input classes" mislabels two items | Apply | Applied | rewritten as a bulleted list |
| M10 | Minor | `push`'s error path unreachable | Apply | Applied | recorded on `push`'s doc |
| M11 | Minor | the deferral's own text disagreed with itself | Apply | Applied | v5 §5 annotated; both readings now satisfied |
| nits | Nit | understated names, `expect` labels, redundant lead, panic contract | Apply | Applied | Pass |
| — | — | a proptest over the round trip | Defer | Deferred | see §3 |

## 2. M1 — the correction that changes what the step is worth

Three places claimed the four properties were latent because `fill_raw_read` is only ever handed a
fresh buffer, and that C2 would be the change that first hands it a used one. The claim came from
B2's deferred finding and was repeated without checking.

It is false. `ReadFilter::next` refills **one** `NoodlesRawAlignedRead` for a whole pass, so every
read after the first has arrived with a history since B2. The reviewer measured it by instrumenting
`fill_raw_read` and running the CRAM cursor walk; asserting freshness turns three *existing* tests
red.

**This makes C1b more important than its own doc claimed, not less.** The reason the gap survived
is not that the condition is absent but that nothing checks it: no test anywhere compares a served
read's bases against an independent expectation, so deleting `sequence.clear()` grows sequences past
300,000 bases with the pre-existing suite — CRAM-versus-BAM oracle included — still green. A
regression in these clears would corrupt production reads today.

Corrected in the module doc, in the two test docs that repeated it, on `fill_raw_read`'s own doc,
and in the origin report (`fixes_applied_2026-08-03_v5.md` §5), whose stated reason for deferring
was the same false claim.

## 3. Deferred, with a reason

**A property test over the round trip.** `push` + `fill_raw_read` is a round-trip law over a
structured input domain, and the reviewer specified a proptest that would subsume the scalar gap
automatically and reach orderings the fixed fixtures do not — long→short, short→long,
unnamed→named, and repeats of one index through one shared buffer.

It was **specified but not compiled** by the reviewer, and it needs an `arb_record_buf()` generator
covering name-`None`/empty/arbitrary-bytes and zero-length sequence, quality and CIGAR. That is its
own piece of work, and the five hand-written tests added here already kill every mutation it was
proposed to catch. Recorded rather than rushed into a step whose subject is a deferred finding from
two steps ago; it belongs with whoever next touches this module's coverage.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,860 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**Suite 2,847 → 2,860 (+13)**: the original 9, plus 4 added by this review pass. Fully accounted.

### Mutations run after the fixes

Each marker `grep -c`-confirmed present before the run, file byte-restored after.

| mutation | killed by |
|---|---|
| `flags` read from `self.index[0]` | `every_scalar_is_read_from_the_entry_asked_for`, **alone** (was: nothing, whole suite green) |
| `Span::new` truncates instead of refusing | `a_span_past_the_index_width_is_refused_rather_than_truncated`, **alone** (was: nothing) |
| `shrink_to_fit` emptied | `shrinking_gives_back_the_slack_the_buffers_grew_by`, **alone** (was: nothing) |
| name replaced wholesale instead of cleared and refilled | `serving_a_second_read_reuses_the_first_reads_allocations`, **alone** (was: nothing — behaviourally identical) |
| `None => *out.name_mut() = None` → `None => {}` | the unnamed-record test, alone |
| `out.data_mut().clear()` deleted | the aux-tag test, alone |
| `sequence.clear()` removed | the empty-record and shorter-record tests |
| `&self.index[i]` → `&self.index[0]` | 5 tests |

Four mutations that previously survived the entire 2,856-test suite are now each killed by exactly
one test.

## 5. Disputed findings

None.

## 6. Failed-validation findings

None. The one finding whose suggested code the reviewer had not compiled — the allocation-reuse
test — compiled and passed as written once the capacities were read through the `_mut()` accessors,
which the reviewer flagged as the open question.
