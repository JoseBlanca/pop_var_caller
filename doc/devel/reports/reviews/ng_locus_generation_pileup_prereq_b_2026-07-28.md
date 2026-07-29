# Code Review: ng generic locus generator — prerequisites, Milestone B

**Date:** 2026-07-28
**Reviewer:** rust-code-review skill (orchestrator), 4 category sub-agents (6 checklists)
**Scope:** the Milestone B diff — the shared locus type: `ReadCoverage` reshaped, `ObservedSequence`
extended
**Status:** Request-changes → **all applied**; see §11

---

### 1. Scope

- **What was reviewed:** a diff — `d80bbe2..e76fe2a`, four commits (`6750d73` B1, `5054515` B2,
  `45dab6b` the fixture strengthening, `e76fe2a` B3 docs).
- **In-scope files:** [locus_generation/mod.rs](../../../../src/ng/locus_generation/mod.rs),
  [locus_generation/ssr.rs](../../../../src/ng/locus_generation/ssr.rs), and the six
  `examples/ng_ssr_*.rs` the diff touches; the three specs for the intent check.
- **Out of scope:** production (`src/pileup/`, `src/psp/`, `src/var_calling/`, `src/vcf/`);
  `src/ng/read/input/` (Milestone A, reviewed separately).
- **Categories dispatched:** `reliability` (this step is flagged silent-failure-prone, so test
  discrimination is the whole question), `refactor_safety` (a type reshape across a shared type,
  done partly by search/replace), `idiomatic` + `naming` (new constructors and predicates), `smells`
  + `extras` (hot path; and "diff matches stated intent" against the plan's two named oracles).

One review over the milestone rather than one per step, as at Milestone A and for the same reason.

### 2. Verdict

**Request-changes, and the reviewers earned it.** One Blocker and eight Majors, of which **six were
"a test that cannot fail"** — including one I introduced by weakening an existing test in a
mechanical rewrite, and one where my own doc comment claimed coverage the test did not have. All are
applied in `934ea0f`.

The reliability agent did not judge by reading: it **mutated the source and ran the suite** for every
claim, and quoted the output. That is what turned three of these from opinion into fact.

### 3. Execution status

Container (`./scripts/dev.sh`), verbatim:

- `cargo fmt --all --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics.
- `cargo test --all-features` — `2495 passed; 0 failed` at review time; **2503 after the fixes**.
- `cargo test --examples` — green. Host-native: 2503.
- **Not run:** `cargo test --all-targets --all-features` (pre-existing bench panic), `cargo doc
  --no-deps` (11 pre-existing bad intra-doc links) — PROJECT_STATUS *Standing project-wide items*.

**Mutations run to verify the fixes** (each applied, suite run, reverted):

| mutation | before the fix | after |
|---|---|---|
| `segment.start()` → `locus.margin_start.get()` | 61 passed, 0 failed | `placed_left_is_counted_against_the_tract_anchor_not_the_margin` FAILS |
| both `Complete` arms of the flushness predicates → `false` | all passed | `flushness_is_derived_from_where_the_run_sits` FAILS |
| delete the read-group `then_with` | passed 6 runs in 10 | `rows_sort_by_group_within_an_allele_deterministically` FAILS 3/3 |
| swap `from_left`/`from_right` at the mint | (B1: dump byte-identical) | dump FAILS since `45dab6b` |

### 4. Open questions

1. **Should `Observed` be a private-field newtype?** (M6) Its fields are `pub`, so the constructors'
   clamping is convention, not an invariant — and `num_obs_along_locus` clamping "defensively" with a
   comment saying the producer enforces it is the tell. **Not applied:** a type-shape change to the
   type Checkpoint B exists to freeze.
2. **Should `locus_len` be a `LocusLen` newtype?** (M7, convergent ×3) `from_left(10, 4)` and
   `from_left(4, 10)` both compile and the clamping hides the transposition; and the length fed to
   the mint (`segment.tract_len()`) comes from a different expression than the one every dump feeds
   the predicate (`region.len()`). **Not applied**, same reason.
3. **Should the expanded-allele merge be prevented at all?** (M1) It is arguably correct — identical
   constraints are one cell — but it is not the pre-reshape answer. Preventing it means keeping a
   side tag, which is the encoding the spec rejected. **Recorded and tested, not changed.**
4. **Should the two dashboard dumps gain a read-group column?** (M8) They now emit one row per cell
   with no way to tell the cells apart, and the cohort stutter dump is the analysis that motivated
   the field. **Not applied:** it changes an artifact the marimo dashboards parse.

### 5. Top 3 priorities

All applied. 1: the `placed_left` Blocker. 2: the expanded-allele merge — a behaviour change
documented as a labelling one. 3: the three predicates/tests that could not fail.

### 6. Findings

#### Blocker

**B1: [ssr.rs](../../../../src/ng/locus_generation/ssr.rs) — `placed_left`'s anchor was pinned by
nothing.** **Category:** reliability. **Confidence:** High (mutation-verified). *Applied.*

`tally` takes `locus_start` as a bare `u64`; its unit test supplies that argument itself, so it
pinned `<` versus `<=` and nothing about which coordinate the generator passes. A 0-based/1-based
slip, the margin start, or a widened region's start would all land identically unnoticed — on the
term production subtracts from QUAL. Fixed by two generator-level tests that bracket the anchor from
both sides on a fixture where the two candidate anchors give 4 and 0.

#### Major

**M1 — the expanded-allele collapse is a *merge*, not just a label.**
**Categories:** reliability, refactor_safety, idiomatic, smells, extras (convergent, 5).
*Applied (test + doc).* Because the STR reach is in read bases, a saturating run makes
`from_left` and `from_right` return the same value, so two opposite-sided partials with equal bases
share a bucket key and become one row. Unreachable on any fixture whose reads are exact reference
slices. The plan's stated equivalence `PartialRight(n) ⇔ Observed { len - n, n }` stops holding at
`n = len`.

**M2 — `partial_reach_beyond_locus_is_clamped` was weakened by the mechanical rewrite.**
**Category:** reliability. *Applied.* Its two halves construct the same value post-reshape, so it ran
one case twice while its doc claimed "on **both** ends". Split into a clamp test plus two
value-level constructor tests.

**M3 — `is_flush_left` / `is_flush_right` had no test.** **Categories:** reliability, idiomatic.
*Applied.* Both `Complete` arms are unreachable from every call site, so inverting them changed
nothing.

**M4 — the read-group tie-break's determinism test was probabilistic.** **Category:** reliability.
*Applied.* 6 green runs in 10 with the tie-break deleted.

**M5 — `num_obs_along_locus` had no interior-run case.** **Category:** reliability. *Applied.* The
case the reshape exists to represent, and the only one where a wrong window neither panics nor
clamps.

**M6 — exhaustive matches traded for guard chains ending in `_`.** **Category:** refactor_safety.
*Applied.* Five label sites now destructure `Observed`, so the compiler can force them again.

**M7 — the clamping invariant is unenforced / the two `u16`s are transposable.**
**Categories:** idiomatic, naming, smells (convergent, 3). *Deferred to §4.1–4.2.*

**M8 — row-splitting silently changed the two dashboard dumps.** **Category:** extras. *Documented,
column deferred to §4.4.*

#### Minor and nits (applied)

Six stale comments (the `(bases, read_coverage)` dedup banner; "six reads / two partials / one left,
one right" prose `45dab6b` invalidated); the spec's flush-right predicate written `==` where the code
uses `>=`; the two STR specs `locus_generation_ssr.md` / `read_preparation_ssr.md`, which B3 skipped
though spec §10 named them as fold-in homes.

Not applied, recorded: `placed_left` reads like a boolean beside `num_obs`/`num_fwd` (kept — it is
production's name and the spec's); `tally`'s `locus_start: u64` where `Position` exists; the dump's
`render()` no longer injective on rows; `partial:interior` is dead on the STR path; `coverage_label`
now triplicated verbatim across three examples.

### 7. Out of scope observations

The "byte-identical dump" oracle both milestones rest on is a **manual** check with no committed
golden artefact — worth a snapshot file if a third milestone leans on it.

### 8. Missing tests added now

`placed_left_is_counted_against_the_tract_anchor_not_the_margin`,
`a_read_starting_on_the_tract_anchor_is_not_placed_left`,
`an_expanded_allele_merges_the_two_sides_into_one_row`,
`rows_sort_by_group_within_an_allele_deterministically`,
`flushness_is_derived_from_where_the_run_sits`,
`depth_over_an_interior_run_raises_only_the_witnessed_stretch`,
`from_right_places_the_run_against_the_right_border`,
`from_left_and_from_right_agree_once_the_reach_covers_the_whole_locus`.

### 9. What's good

- **B1's fixture repair is real**, and the reviewers confirmed it independently: pre-`pl2` the side
  swap was byte-identical, post-`pl2` it fails.
- **`coverage_order` is total and injective**, and the sort key equals the bucket key — so the
  determinism claim holds as written.
- **Every `locus_len` supplied at a mint or predicate site is the right quantity**, checked
  independently by two categories.
- **The clamp-then-derive order** in `from_right` is right, and its stated reason is right.

### 10. Commands to re-verify

`./scripts/dev.sh cargo fmt --all --check`; `… clippy --all-targets --all-features -- -D warnings`;
`… cargo test --all-features` (expect 2503); `… cargo test --examples`. Mutations: §3's table.

### 11. Author response

- **B1, M1–M6** — fixed in `934ea0f`, each mutation-verified where a mutation applies.
- **M7, M8** — deferred to Checkpoint B as owner decisions (§4); both are changes to a type or an
  artifact this checkpoint exists to settle.
- Minors/nits — applied except the five recorded above.

Per-category audit trail: `tmp/review_2026-07-28_ng-shared-locus-type/` (gitignored).
