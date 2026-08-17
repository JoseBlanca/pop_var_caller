# Fixes applied — ng cohort merge, step A2

*2026-08-17, branch `ng-cohort-merge`. Input:
[the A2 review](ng_cohort_merge_a2_2026-08-17.md) — 4 Major, 6 Minor, 8 Nits over four
category checklists. Every finding is accounted for below.*

## Findings table

| ID | Title | Severity | Decision | Status |
|---|---|---|---|---|
| M1 | no guard against bases of a different length — indels | Major | Apply | **Applied** |
| M2 | two more open-coded copies in library code; doc claimed exclusivity | Major | Apply | **Applied with adaptation** |
| M3 | the complete-only rule's recorded cost had no test | Major | Apply | **Applied** |
| M4 | `reach` doc claimed agreement "everywhere" — false at the ceiling | Major | Apply | **Applied** |
| Mi1 | `saturating_add` untested at its boundary | Minor | Apply | **Applied** |
| Mi2 | `reach` pinned at three points where a property was available | Minor | Apply | **Applied** |
| Mi3 | inverted-region test claimed agreement in prose | Minor | Apply | **Applied** |
| Mi4 | `GenomeRegion::len` defect recorded only in a test comment | Minor | Apply (doc only) | **Applied with adaptation** |
| Mi5 | two crate names for one quantity (`alt_reads`) | Minor | Apply (doc) | **Applied with adaptation** |
| Mi6 | `reach` doc said "locus" where the arch says "observation" | Minor | Apply | **Applied** |
| Nits | eight, listed in the review | Nit | Apply / Dispute | **6 Applied, 2 Disputed** |
| — | arch §2.3 vs arch cohort_merge §2 wording collision | — | Ask | **Raised at Checkpoint A** |
| — | `GenomeRegion::len`'s arithmetic; the census's untested multi-base path | — | Defer | **Out of scope, recorded** |

## The four Majors

**M1 — three assertions, and they are the ones that matter.** The predicate had no
fixture comparing bases of a *different length* from the reference, so three containment
mutants passed the whole suite while silently reclassifying indels as reference.
`matches_reference_compares_the_bases_it_is_given` now covers a trailing deletion (`AC`
against `ACGT`) and insertions on both sides (`ACGTT`, `TACGT`). **Verified after
applying:** `reference_bases.starts_with(&self.bases)` and
`self.bases.starts_with(reference_bases)` each now fail exactly that test — `28 passed;
1 failed` in both cases, against `29 passed; 0 failed` on the real code.

**M2 — the claim and the code were made true together, rather than one being chosen over
the other.** All four surviving open-coded spellings now call the predicate:
`depth_and_alt_reads.rs:82` and `:234` in library code, and two helpers in
`locus_generation/pileup/tests.rs`. The doc no longer claims to be "the one definition of
non-reference in the codebase" — it claims to be the one place the *comparison* is
written, and says in as many words that which observations to ask about (the subset, the
depth cap, the read-group grain) is still the generic pre-pass's, citing that module and
`arch/parameter_prepass_generic.md` §2.3. That makes both architecture documents true as
written, which choosing a side would not have. **The adaptation:** the review offered
"move the sites" or "weaken the claim"; the fix does both, because the collision was
between two different senses of one word, not between two designs.

**M3 — `a_variant_seen_only_by_partial_reads_is_not_counted`.** A partial that disagreed
with the reference over the two bases it saw, carrying 7 reads, contributing 0. The
existing partial test could not see this: a partial that *agreed* answers 2 either way.

**M4 — the sentence was wrong and the code was right.** At `u64::MAX` production's
expression gives `18446744073709551614` and this one `18446744073709551615`; the last
base of a one-base region at `u64::MAX` is `u64::MAX`. The claim is narrowed to "every
region below the top of the coordinate space", the divergence is stated with its
direction, and the ceiling test now asserts `production_reach(u64::MAX, 1) ==
u64::MAX - 1`.

## The Minors

`non_reference_reads_saturates_rather_than_wrapping` pins the boundary — **verified:
`wrapping_add` now fails exactly that test.** A proptest checks `reach` against
production's arithmetic over the whole well-formed domain, bounded below the ceiling
because that is the one documented exception. The inverted-region test calls
`production_reach` instead of asserting a literal. `types.rs` carries the
`GenomeRegion::len` warning where a reader meets it, naming the release behaviour (a
length of 0, not a panic) and the test that pins the consequence — **doc only; the
arithmetic is untouched**, being a shared type outside this step's blast radius.
`non_reference_reads`'s doc names the pre-pass's `alt_reads` as the same quantity, so a
grep works from either side, and explains why `alt` survives in `min_alt_obs`. `reach`
says "observation" where spec §1.3 reserves *locus* for the cohort locus.

## Disputed

- **The `..locus(region, observed)` functional update in the new fixture.** The reviewer
  filed it as a nit against itself: `locus` spells all six fields, so a field added to
  `SampleLocusObservations` is a compile error there, and the signal is not lost.
- **A destructuring rebuild of the same fixture** — same reason.

## Validation

Re-run in the container after the fixes: `cargo fmt --check` clean;
`cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::locus_generation::tests` 29 passed, 0 failed;
`cargo test --lib ng::parameter_estimation::generic` 262 passed, 0 failed, 5 ignored —
the pre-pass suite, which is what vouches for M2's two moved call sites. Full-suite
figures are in the commit message.

**Three mutants re-run after the fix, each now killed:** the two containment mutants
(M1) and `wrapping_add` (Mi1).
