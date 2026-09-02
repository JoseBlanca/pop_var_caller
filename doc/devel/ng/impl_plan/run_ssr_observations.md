# ng STR observations — implementation plan

**Status:** draft, 2026-09-02. Turns the settled design of
[`../spec/run_ssr_observations.md`](../spec/run_ssr_observations.md) — the routing policy,
the per-sample tract slot, and the merge carrying the locus kind — into build order. **Not a
place for new design**; every step cites the spec section it realizes.

**Runs in parallel with** [`calling_loop_ssr.md`](calling_loop_ssr.md), on its own branch and
worktree (§Branch). The two plans meet at two seams, both owned here: the
`CohortObservation.kind` field (Milestone A — merged to `main` first, because everything in
the loop plan's selection work reads it) and the run driver's temporary set-aside guard
(Milestone C — which the loop plan later replaces with its dispatch).

---

## Scope

**In:** `pop_var_caller_exp/call_from_alignments.rs` (routing flags, criteria into
`segments_over`, the temporary driver guard's counts surfaced);
`ng/run/cohort_merge/{close,build}.rs` (the kind field); `ng/run/walker.rs`
(`generic_path_generators` fills the tract slot); `ng/locus_generation/ssr.rs` (the
`read_filter_counts` override); `ng/run/report.rs` (wording that survives a filled slot);
`ng/calling/parameters_file/` (the run's routing criteria recorded in the file it writes);
one additive section in [`../spec/parameters_file.md`](../spec/parameters_file.md) settling
that record's spelling before it is coded.

**Out (with owners):**

- **Everything that calls a tract** — selection, the driver's real dispatch, the record —
  [`calling_loop_ssr.md`](calling_loop_ssr.md).
- **Bundles and satellites** — unchanged, still counted refusals (spec §8).
- **The routing-frontier measurement** (where the line *should* run) — spec §8's deferral.
- **The fit-mode command** — deferred future work (owner, 2026-09-02).

## Principles (how the order was chosen)

- **The cross-plan seam first, and alone.** Milestone A is one field on a merged, tested
  module, touches no other change of this plan, and unblocks the parallel plan — so it lands
  and merges before anything else, exactly as the selection plan originally ordered it.
- **Behaviour change before new machinery.** The routing fix (Milestone B) changes which
  ground existing code sees and is measurable on its own; the new machinery (Milestone C)
  comes after, so a benchmark shift can be attributed to one milestone.
- **Isolate the silent steps.** Two steps fail quietly rather than loudly — the
  `read_filter_counts` override (drop rates under-report with every number plausible) and
  the routing swap (a wrong criteria value routes wrongly genome-wide with nothing crashing).
  Each is its own commit with its oracle green before and after.
- **Container builds**: all `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place)

- The spec, with its decisions closed — the owner's catalog-is-candidates ruling, the CLI
  spelling (spec §9, resolved), the kind-as-a-field decision.
- [`SsrGenerator`](../../../../src/ng/locus_generation/ssr.rs) built and tested;
  [`StrRepeatCriteria`](../../../../src/ng/repeat_catalog/criteria.rs) with its admissibility
  check; the catalog segments reader parameterised by criteria.
- The routing oracle:
  [`examples/ng_typed_region_dump.rs`](../../../../examples/ng_typed_region_dump.rs)
  reproduces a run's routing byte-for-byte (verified 2026-09-02: 539,460 / 32,577 bases on
  HG002 against the run's own report).
- The measured baseline to improve on:
  [`ng_str_path_losses_2026-09-02.md`](../../reports/ng_str_path_losses_2026-09-02.md).

## Branch, worktree, and the parallel plan

- **Branch** `ng-ssr-observations`, from `main`. **Worktree**
  `../pop_var_caller-ssr-observations` (`git worktree add ../pop_var_caller-ssr-observations
  ng-ssr-observations`) — the loop plan works in its own worktree at the same time, so
  neither builds in the other's tree. Use absolute paths in commands; a relative path from a
  reset shell can silently hit the main checkout. Remove the worktree when the plan
  completes.
- **File ownership while both plans are in flight** — do not edit the other plan's files;
  changes cross only through `main`:

  | this plan owns | the loop plan owns |
  |---|---|
  | `run/cohort_merge/`, `run/walker.rs`, `run/segments.rs`, `run/report.rs` | `calling/` (selection, inference, quality) |
  | `pop_var_caller_exp/call_from_alignments.rs` | `run/callers.rs` (after this plan's guard merges), `run/records.rs` |
  | `locus_generation/ssr.rs` | `examples/` probes, benchmark harness scripts |
  | `calling/parameters_file/` (criteria record) | |

- **Merge to `main` at every checkpoint**, Milestone A's immediately — the loop plan rebases
  on it. The one shared touch point is `call_one_generic_locus`'s call site: this plan adds
  the set-aside guard (C2); the loop plan replaces it with the real dispatch after this
  plan's Checkpoint C has merged.

---

## The steps

### Milestone A — the merge carries the locus kind (the seam; merge first)

**A1. `ClosedLocus` and `CohortObservation` carry `LocusKind`.**  ✅
As [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) A1 wrote it, now delivered here
(its ownership note points back): the closer has the kind in scope where it builds a
`ClosedLocus` ([`close.rs:713-721`](../../../../src/ng/run/cohort_merge/close.rs)) and drops
it; carry it through, and have `CohortObservation::over`
([`build.rs:996`](../../../../src/ng/run/cohort_merge/build.rs)) clone it onto the assembled
locus. The clone is per *built* locus and boxes two flanks only on the `Ssr` variant — say so
in the doc comment; `Generic` clones nothing. No verdict moves; no test changes outcome
beyond constructors naming a kind.
*Depends:* none. *Source:* spec §4; [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §2.

> **Checkpoint A:** the merge's tests unchanged in outcome, the field populated on both
> kinds. **Reviewed and merged to `main` at once** — the loop plan's selection work reads
> this field from its first step. Pause for review.

### Milestone B — the routing policy is the caller's

**B1. The criteria from flags, and the refusals.**  ✅
`call-from-alignments` gains `--min-copies`, `--min-period`, `--max-period`,
`--max-str-len`, `--min-purity`, spelled as `type-regions` spells them, defaults the calling
floors; `segments_over` builds its `StrRepeatCriteria` from them instead of `default()`
([`call_from_alignments.rs:845`](../../../../src/pop_var_caller_exp/call_from_alignments.rs)).
Asking below the catalog's floors surfaces the existing typed refusal with the flag's name in
the message. **Own commit — a wrong value here routes wrongly genome-wide with nothing
crashing**; the oracle is B3's parity check.
*Depends:* none (parallel to A). *Source:* spec §2.1, §2.3, §9.

**B2. The parameters file records what the run routed with.**  ✅
Documentation first: one additive section in
[`../spec/parameters_file.md`](../spec/parameters_file.md) naming the routing-criteria
record's TOML spelling, beside `[fitted_from]`'s conventions. Then the run's written file
carries it, and a supplied file's differing criteria are *visible, never blocking* — the
owner's ruling, spec §2.3. No comparison refuses anything.
*Depends:* B1. *Source:* spec §2.3.

**B3. Routing parity, and the no-regression pin.**  ✅
Two tests against external oracles: the run's ground partition equals
`ng_typed_region_dump`'s at the same criteria (the spec's standing oracle, kept as a fixture
test on the committed synthetic reference); and on ground that is generic under both the old
and the new default, the VCF is byte-identical to a pre-change run — the change is *which*
ground is generic, never what happens on it.
*Depends:* B1. *Source:* spec §10.

**B4. The recovery, measured.**  ✅
Re-run the GIAB per-sample benchmark at 30× and 50× with the new default routing and
re-score (the loss report's own scripts). Expected from the upper bound: overall recall from
0.935 toward ≈ 0.97 (SNPs) and 0.818 toward ≈ 0.94 (indels); record what is actually
measured beside the prediction in a short report, and update the benchmark dashboards' ng
rows.
*Depends:* B1. *Source:* spec §2.2, §10;
[the loss report](../../reports/ng_str_path_losses_2026-09-02.md) §5.

> **Checkpoint B:** flags in, criteria recorded, parity and no-regression green, the
> recovery measured against its prediction. Merge to `main`. Pause for review.

### Milestone C — the tract slot is filled, and the accounting is paid

**C1. `SsrGenerator` reports its read-filter counts.**  ✅
Override [`read_filter_counts`](../../../../src/ng/locus_generation/mod.rs) (default: empty)
with the generator's own reader tallies, per read group. **Own commit, do not bundle** — the
failure is silent under-reporting with every number plausible; the oracle is a test that
walks tract ground and asserts non-zero per-group drops where reads were filtered.
*Depends:* none. *Source:* spec §3.2.

**C2. The slot filled, and the driver's temporary set-aside.**  ✅
`generic_path_generators` ([`walker.rs:1588`](../../../../src/ng/run/walker.rs)) builds an
`SsrGenerator` into the `Ssr` slot — unit-robust aligner, `SsrGeneratorConfig::default()`,
same `WalkReference` accessors as the pileup generator. And the driver, on meeting a
non-generic cohort observation, **sets it aside and counts it** — never hands it to
`select_generic` — until the loop plan's dispatch replaces this guard. The bundle slot stays
unfilled.
*Depends:* A1, C1. *Source:* spec §3.1, §5.

**C3. The run report survives both states.**  ☐
The ground partition still sums to 100% with tract regions handled; the *not called* line
now names what is actually unbuilt (bundles, satellites) rather than every repeat; the
set-aside tract loci of C2 are a counted line of their own. Wording per spec: a smaller,
honest *not called* line, not a zero and not the old sentence.
*Depends:* C2. *Source:* spec §3.2, §5.

**C4. Invariance, end to end.**  ☐
The E2 oracle re-run with the slot filled: byte-identical VCF at pools of 1–16 across two
building-region widths on the concurrency fixture, and the report's counts identical too.
*Depends:* C2. *Source:* spec §3.3, §10.

> **Checkpoint C:** slot filled, guard in place, report honest, invariance green. Merge to
> `main` — **this unblocks the loop plan's driver milestone.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the merge's existing tests unchanged in outcome; the field populated on both kinds |
| B | **external, twice**: `ng_typed_region_dump` parity on the ground partition; byte-identical VCF on ground generic under both defaults; the GIAB recovery measured against the report's prediction |
| C | non-zero per-group drop tallies on tract ground (the silent-failure oracle); partition sums to 100% in both states; E2 byte-identity at 1–16 threads |

## Out of scope (next plans)

- **The dispatch that calls a tract** — [`calling_loop_ssr.md`](calling_loop_ssr.md)
  Milestone C replaces C2's guard.
- **Bundles** — their own spec first (spec §8).
- **The routing-frontier measurement** — cheap once both paths run; spec §8's home.
- **The fit-mode command** — deferred future work (owner, 2026-09-02).
