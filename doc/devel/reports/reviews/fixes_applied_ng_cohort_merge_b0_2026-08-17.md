# Fixes applied — ng cohort merge B0

*2026-08-17, against [the review](ng_cohort_merge_b0_2026-08-17.md) of plan step B0
([impl report](../implementations/ng_cohort_merge_b0_2026-08-17.md)). Every finding is answered
below.*

## The two that mattered, and both are now measured

**The Blocker: "every read is named" is asserted as a count, not as presence.** The differential's
new check read `!chain_ids.is_empty()`, which a walk naming one read in a hundred satisfies. It now
reads

```rust
2 * observation.chain_ids.len() >= observation.num_obs as usize
```

— not `== num_obs`, because a read *pair* whose two mates both reach one observation collapses onto
a single id, and no more than two can. Re-run against the mutation that survived the review:
truncating each observation's ids to the first one now fails, in each mint path separately.

**The Major: the walk comparison accepted a mid-region renumbering.** The helper kept one renaming
map per locus; it now keeps one per call, and the tiling test calls it **twice — once per region**,
which is what pins the renumbering to the join rather than merely tolerating it somewhere. The
helper is now tested itself: a renaming that holds throughout is accepted, one that changes at the
second position is refused, and two reads collapsing onto one identity is refused from the other
side.

Mutants re-run after the fixes, one at a time from a pristine copy
(`./scripts/dev.sh cargo test --lib ng::locus_generation`):

| mutation | before the fixes | after |
|---|---|---|
| general fold keeps only each observation's first id | survived (whole suite green) | **FAILED**, 359 passed / 3 failed |
| fast lane keeps only each observation's first id | survived | **FAILED**, 361 / 1 |
| fast lane withholds a reference-matching read's id (the revert) | caught | **FAILED**, 361 / 1 |
| the helper's map rebuilt per locus (the weakening) | not detectable | **FAILED**, 363 / 2 |

## Every finding

| id | what it said | what was done |
|---|---|---|
| **B1** | naming pinned as presence, not as count | **Applied**, as above. |
| **M1** | the helper accepts a mid-region renumbering | **Applied**: one map per call, the tiling test compares region by region, and three tests now pin the helper's own behaviour. |
| **M2** | a comment states the deleted rule | **Applied**: it now says which rule each side follows, and records that this is the third rule the comparison has lived under — the earlier two are why it is not simply a subset check. |
| **Mi1** | the public `chain_ids` doc predates the ruling | **Applied**: it carries both questions the id answers, that an id names a read *within one walk*, and that the STR path has none and needs none. |
| **Mi2** | the two new reference-matching filters drop the complete-observation guard | **Applied** in both places, with the reason `matches_reference` gives. |
| **Mi3** | `assert_eq!(len, len)` left behind | **Applied**: both removed, their messages moved onto the helper's call. |
| **Mi4** | stale "the walker drops REF chain ids" comments in the copied suite | **Applied** where the comment is about ng; production's own copied text is left alone, since it is still true of production. |
| **Nits** | five | **Applied**: the clone became a `mem::take`, the failure dumps the offending observation only, the doc link points at the field, and the unreachable disjunct is gone with the new form of the assertion. The two maps keep their names, which say which direction each carries. |

## Claims corrected in the report

- **The cost is not only eight bytes a read.** An observation that used to hold an empty `Vec` now
  holds a populated one, so it is an allocation per observation per position as well.
- **Its size had already been measured, by production.** The ids production drops are ~96.6% of all
  chain ids on real cohorts (`pileup/walker/open_record.rs:155`), so keeping them all is about
  thirty times the ids either caller used to carry. ng's old rule was per read rather than per
  bucket, so the two do not withhold exactly the same set; the report says so.
- **The segment rule is stated for observations, not loci**, and there is a **second dependency of
  the same kind**: §6.2 makes segment independence conditional on every sample sharing one
  segmentation. Both are now in the report.

## Validation

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::locus_generation` — 365 passed, 0 failed, 1 ignored.
- `cargo test --lib` — **3,682 passed, 0 failed, 11 ignored** (3,679 at review time; the three added
  are the helper's own).
