# Fixes applied — ng read filtering in stages, C4

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `5e3b22f`
**Review:** [`ng_read_filtering_stages_c4_2026-08-03.md`](ng_read_filtering_stages_c4_2026-08-03.md)
**Impl report:** [`ng_read_filtering_stages_c4_2026-08-03.md`](../implementations/ng_read_filtering_stages_c4_2026-08-03.md)

Applied: 7 · Applied with adaptation: 1 · Already fixed: 0 · **Deferred: 1** · Disputed: 0

---

## 1. Findings table

| id | severity | subject | status | validated |
|---|---|---|---|---|
| B1 | Blocker | the window untested on the arm every real walk takes | Applied | Pass — test D kills the mutant, alone |
| B2 | Blocker | "nothing else moves" unpinned for the region in flight | Applied | Pass — tests A and B kill three mutants |
| M1 | Major | `last_region_start = None` survives | Applied | Pass — test C kills it |
| M2 | Major | the `saturating_sub` rationale points the wrong way | **Applied with adaptation** | §2 |
| Mi1 | Minor | the `None`-arm comment is falsified by C4 | Applied | rewritten |
| Mi2 | Minor | the replay property deserves its own test | Applied | test E |
| — | — | six stale `reset_counts` references in the design docs | **Deferred** | §3 |

## 2. M2 — the rationale corrected, the code kept

The review is right that the justification pointed the wrong way. Underflow is **unreachable**:
every contribution to the narrowing's other-sample count is a `+=`, and the baseline is sampled from
that same number, so the baseline can never exceed it. `saturating_sub` guards nothing reachable.

The code is kept — it is free, and a plain `-` would put a panic in an accounting helper for no
benefit — but the comment now says what is true, and names the hazard that **is** real and runs the
other way: `CramAlignedReadsReader` adds a container's foreign-record count each time it decodes
one, so a container decoded twice contributes twice. A window therefore *starts* at zero honestly
and may **over**-report afterwards on CRAM — exactly as the unwindowed number already does, which
that field's own doc has always called container-granular.

**Not fixed, deliberately:** it is the reader's accounting, and it pre-dates the window by two
milestones. Recorded at the code.

## 3. Deferred, and it is the checkpoint's

**Six live `reset_counts` references across the spec, the plan and the arch**, plus the arch sketch
still calling the tally field `counts`. The review is right that the rename was executed by halves.

This skill does not edit the spec, the architecture or the plan — a design-doc change is the owner's.
So it goes to **Checkpoint C** as a concrete either/or: amend the three documents to
`reset_read_group_counts`, or rename the method back to `reset_counts` and accept that it names one
of two tallies by the other's word. The impl report's §4 already records the deviation and its
reason.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,865 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**Suite 2,858 → 2,865 (+7):** two as built, five from the review.

### The six survivors, re-run after the fixes

| mutation | before | after |
|---|---|---|
| the `first_mut` arm ignores the window | **survived 2,860** | killed by `read_group_counts_scopes_the_other_sample_rider_when_the_new_window_met_its_own_reads`, alone |
| reset sets `region = None` | **survived 2,860** | killed (2 tests) |
| reset sets `examined = 0` | **survived 2,860** | killed (2 tests) — the mutant served a read twice |
| reset sets `last_emitted = None` | **survived 2,860** | killed by the order-guard test, alone |
| reset sets `last_region_start = None` | **survived 2,860** | killed (2 tests) |
| `saturating_sub` → `-` | **survived 2,860** | still survives — unreachable by construction, §2 |

## 5. Disputed findings

None. The review contradicted the impl report's rationale for the saturating subtraction and was
right.
