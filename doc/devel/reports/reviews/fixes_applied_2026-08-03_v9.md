# Fixes applied — ng read filtering in stages, C3

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `0b5f5e7`
**Review:** [`ng_read_filtering_stages_c3_2026-08-03.md`](ng_read_filtering_stages_c3_2026-08-03.md)
**Impl report:** [`ng_read_filtering_stages_c3_2026-08-03.md`](../implementations/ng_read_filtering_stages_c3_2026-08-03.md)

Applied: 9 · Applied with adaptation: 0 · Already fixed: 0 · Deferred: 0 · Disputed: 0

Every finding applied. Nothing deferred, because nothing needed a decision.

---

## 1. Findings table

| id | severity | subject | status | validated |
|---|---|---|---|---|
| B1 | Blocker | the held record's read-group replay is untested | Applied | Pass — new test kills both mutations, alone |
| B2 | Blocker | the early stop's `on_this_contig` guard is untested | Applied | Pass — new test kills the mutation, alone |
| M1 | Major | the report's `decode_fails` claim is false, twice | Applied | report corrected; the record moved onto the variant |
| M2 | Major | prose introduced or made stale by this patch | Applied | five sites |
| Mi1 | Minor | `AlignmentCursor`'s doc severed (from C2) | Applied | doc reattached, `cargo doc` still 12 |
| — | — | arch §3.3's `header` / `other_sample_reads` | Applied | recorded for the checkpoint; code unchanged, which the review confirmed is right |

## 2. The two Blockers

Both are lines C3 **moved without changing**, in the `read_next` whose bodies the step's own safety
argument called byte-identical. That argument is right, and it is exactly why the review mutated
rather than trusted it: identical bodies stay as covered as they were, and these two were not
covered at all.

**B1 — the held record's read group.** The early stop fires on a record already handed over, so it
is put back with the group it arrived with; on CRAM that group is decided at container decode and
travels attached, with no `RG` tag to recover it from. Replaying `None` kills a valid
multi-read-group CRAM run at every region boundary; replaying `Some(ReadGroupId(0))` attributes a
read to the **wrong library**, silently. Both survived all 2,856 tests.

`the_held_record_carries_its_own_read_group_into_the_next_region` uses a two-group resolution and
puts the held record in the *second* group, so a replay under the first group's id — or under none
— is visible.

**B2 — the early stop's contig guard.** A sorted file groups contigs, so a record of a later contig
always begins "past" this region's end by position alone. Stopping on it loses every remaining read
of the region. Dropping the guard passed all 2,856; the existing fixture's foreign-contig record
sits *inside* the region, so the position test never fires on it.

`a_record_of_another_contig_does_not_trigger_the_early_stop` puts one past the end, between two
in-region reads.

## 3. M1 — a correction that moved a record as well as fixing a sentence

The claim was that the doubles took "the only construction of `ReadFilterError::Decode`". Verified
against `git show 0b5f5e7`: `FakeRecord` never constructed a `ReadFilterError`, and
`decode_fails: true` appeared nowhere — both constructors set `false`, the two tests that set it
having died at **C2**. Dead code, deleted.

`Decode` *is* unreachable and unpinned, from C2. The note C1 wrote about it lived in the test block
C3 deleted, so the record has been moved onto **the variant itself**, where it survives further
deletions — together with a correction to its doc, which named "the unmapped flag clear yet no
position" as a cause the layer below discards before the conversion is reached.

## 4. Mi1 — a C2 defect found by a C3 reviewer

`AlignmentCursor`'s doc comment was severed by C2's insertion of the error and tally types, which
landed *inside* it. The cursor's rustdoc was the fragment *"is no stream object to give back."*
Repaired by moving the inserted block above the doc. `cargo doc` stays at the 12-link baseline.

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,858 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**Suite 2,857 → 2,858 (+1):** one deleted, two added. Fully accounted.

### Mutations re-run after the fixes

| mutation | before | after |
|---|---|---|
| held record replayed as `Some(ReadGroupId(0))` | **survived 2,856** | killed by the held-record test, alone |
| held record replayed as `None` | **survived 2,856** | killed by the same test |
| the `on_this_contig &&` guard dropped | **survived 2,856** | killed by the foreign-contig test, alone |

## 6. Disputed findings

None. The review contradicted the impl report on `decode_fails` and was right.
