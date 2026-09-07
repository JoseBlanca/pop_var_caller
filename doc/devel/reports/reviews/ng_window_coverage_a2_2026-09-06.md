# Review — window coverage A2: the floor under a window

**Date:** 2026-09-06
**Reviewed:** commit `711c06a5` on `ng-window-coverage` (plan step A2 of
[window_coverage.md](../../ng/impl_plan/window_coverage.md))
**Impl report:** [ng_window_coverage_a2_2026-09-06.md](../implementations/ng_window_coverage_a2_2026-09-06.md)
**Per-category files (audit trail):** `tmp/review_2026-09-06_window-coverage-a2/`

## 1. Scope

The A2 diff only — about sixty lines of source and five tests inside `src/ng/window_coverage/`.
The module itself landed at A1 and was reviewed then; two agents in isolated worktrees split the
categories, one on `reliability`, one on `naming`/`idiomatic`/`defaults`/`errors`/
`refactor_safety`/`smells`.

Out of scope: `src/sample_summary/` (frozen; read as the oracle); the four checks that are red on
`main` before this branch.

## 2. Verdict

**Approve with changes, and one of them was a real defect rather than a missing test: the floor
counted records, not positions.**

## 3. Top three

1. **The floor counted observations, not distinct covered positions.** `observe`'s order rule is
   non-decreasing, so one base can be observed several times; probed, one base observed five
   times cleared a floor of five and was folded into the histogram — which is precisely the
   "as confident as five hundred positions" window the floor exists to silence. Every doc, and
   spec §3.3, say *covered positions*. **Fixed:** the accumulator now carries a second counter,
   the number of different coordinates the running sums sit on, and the floor is judged against
   that. The divisor of the two means is unchanged, so the differential against production still
   matches bit for bit.
2. **No fixture exercised the floor's selectivity.** All five of A2's original tests used streams
   whose verdict was uniform — all absent or all speaking — so a mutant that made absence
   *sticky* (once one window is silenced, every later one is) passed the whole suite. **Fixed:**
   a stream holding one isolated position and a run of four, at a floor of three, so both
   verdicts occur in one stream from one accumulator.
3. **`assert_valid` accepted a floor no window could reach.** A floor of 100 over a 10-base
   window silences every window of every sample and reports nothing about it — in a module whose
   own configuration doc argues that a bad configuration must fail loudly, and at exactly the
   number plan step D1 sweeps by measurement. **Fixed:** the assertion rejects a floor above what
   the window is wide enough to hold, naming both numbers.

## 4. What's good

- **The decision A2 made that the spec does not state — an absent window is not folded — is
  right, and the tests pin both halves of it.** The reliability agent checked the arithmetic
  behind the stated reason (both the GC and the depth cast of `NaN` come back 0, so a folded
  absent window lands in one specific cell) and confirmed that a "count it but do not bin it"
  mutation fails the `windows_folded` assertion while a "bin it" mutation fails the cell-sum
  assertion.
- **The floor set to 1 in every pre-existing test hides nothing, and is stronger than "off".**
  Since a window always holds at least its own centre, a floor of 1 sits exactly on the boundary,
  so those tests and the differential are themselves what catch a `<` widened to `<=`.
- **Refactor safety is actively good here.** The configuration has no `Default`, so the new field
  is a compile error at every literal until each one decides about it, and `finish` destructures
  it exhaustively with a named ignore rather than `..`.

## 5. Findings

### Major

**M1 — the floor counted records rather than distinct coordinates.** *Category: reliability.*
Described above. **Applied**, with `a_repeated_position_does_not_clear_the_floor_on_its_own`.

**M2 — no fixture straddled the floor.** *Category: reliability.* Described above. **Applied**,
with `a_stream_holding_both_a_silent_window_and_a_speaking_one_is_scored_per_window`.

**M3 — a floor no window could reach was accepted.** *Category: reliability.* Described above.
**Applied**, with `new_panics_on_a_floor_no_window_could_reach`.

**M4 — `is_absent` decided a two-field invariant from one field.** *Categories: idiomatic,
reliability (convergent).* Both fields are `pub`, so `WindowCoverage { gc_fraction: NAN,
mean_depth: 12.0 }` is writable crate-wide and `is_absent` called it present — a real depth over
a missing GC fraction scored as a measurement. The first builder outside this module arrives at
plan step C2. **Applied** as fail-closed on either field, with
`is_absent_is_true_when_only_one_field_is_missing`. The reviewer also proposed a `debug_assert`
that the two fields agree; **not taken** — `is_absent` is a question, and a question that
panics is a worse contract than one that answers safely. The stronger lock (a private field and
a constructor) is a follow-up for C2, when a builder outside this module exists; today a
two-`f32` positional constructor would buy a transposition hazard for a lock nothing needs.

### Minor

**Mi1 — the configuration's own doc listed two constants and A2 added a third.** *Category:
defaults.* The type-level doc — the site a caller reads when deciding what to put in each field —
said "the settled values live in `WINDOW_BP` and `GC_BINS`", so a reader meets
`MIN_WINDOW_POSITIONS` beside them with nothing at that site saying it is of a different kind.
**Applied:** the doc now says which two are production's and inherited, and that the third is
provisional until D1 defends it with data.

**Mi2 — seven doc claims the floor made stale.** *Categories: naming, reliability (convergent).*
The private method's doc was updated at A2 and the public ones were not: the accumulator's type
doc still promised a fold per finalised window, its `windows_folded` field doc still said "one
per covered position", the module doc still said ng's version was "about to diverge", and the
differential's doc still listed the floor as a future departure. **All applied.**

**Mi3 — two exits from `finalise_centre`, each advancing the centre cursor.** *Category: smells.*
`finalise_all` loops on that cursor, so an exit that forgot it would hang the run rather than
fail a test — and A3 adds a third path through this method. **Applied:** one branch produces the
value, one push and one increment at the bottom.

**Mi4 — the floor's effect leaves no trace a consumer can read.** *Category: defaults.* `finish`
consumes the accumulator, and the histogram carries neither the floor nor how often it fired,
while `windows_folded` now depends on it — so a consumer cannot tell a thinly-covered sample from
one whose windows were mostly silenced. The reviewer proposed
`CoverageByGcHistogram::windows_under_the_floor`. **Deferred, and raised with the owner at
Checkpoint A**: `CoverageByGcHistogram` is the interface between this plan and the hidden-paralog
filter's, which another branch is building from spec §3.5's stated shape, so adding a field to it
is not this step's to decide. The test-side benefit was taken without the type change — the
differential now asserts that no window it compares is absent, which makes "the floor is switched
off in this comparison" a checked claim rather than a comment.

### Nits

Applied: a section heading that ended up labelling the floor tests instead of the panic tests it
was written for; an assertion message carrying eighteen spaces from a hand-wrapped string; a
positive control inside `an_n_position_does_not_count_toward_the_floor`, so its assertion is about
the `N` not counting rather than about the fixture being thin; and
`is_absent_accepts_any_nan_payload_but_equality_does_not`, which records that a `NaN` arriving
through a codec is absent but does not compare equal to this module's own absent value — so
`is_absent` is the safe question and `== absent()` is not.

## 6. Verification

- **Mutation testing.** The reliability agent ran 13 mutations against A2 as committed: 11
  caught, 2 survived, one of those (`is_absent` reading `gc_fraction` instead of `mean_depth`)
  indistinguishable on any value the module can produce. The real survivor was sticky absence,
  which M2's new test now catches.
- **Four mutations re-run by the orchestrator on the fixed tree**, each failing exactly one test
  (36 passed, 1 failed, every time): the floor judged against records rather than distinct
  coordinates; absence made sticky; the distinct count left un-reset at a contig change; and the
  unreachable-floor assertion removed.
- `cargo test --lib --all-features` — **6,312 passed, 0 failed, 15 ignored** (6,306 before the
  fixes; 6,301 at A1; 6,275 without the module). 37 in this module.
- `cargo clippy --lib --all-features` — no diagnostic in `src/ng/window_coverage/`.
  `rustfmt --check` clean.
- The calling oracle was not re-run: nothing calls the module, so no VCF byte can move.

## 7. For the owner, at Checkpoint A

Two of these are spec sentences and one is a shared type, so none is the implementer's to change:

- **Spec §3.4 now contradicts the code.** It says "Every finalised window is folded into a
  per-sample count matrix"; a window under the floor is finalised and is not folded, for the
  arithmetic reason in §5's Mi4 above. §3.3 would carry the sentence.
- **`windows_under_the_floor` on the histogram** — recommended, but it changes the type the
  filter plan reads.
- **The 12-byte deque entry in spec §5**, carried over from A1's review: the entry §3.6 actually
  specifies is 24 bytes.
