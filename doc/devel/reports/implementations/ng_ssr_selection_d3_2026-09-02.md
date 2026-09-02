# D3 — what the genotype prior takes, and the one it no longer does

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone D step D3, and **Checkpoint D**. **Design:**
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §2.3 (and its §4 and §5
retractions); [`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §3.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs),
`SsrLocusSelection` and `repeat_count_of_bases`.

---

## What landed

**`SsrLocusSelection`** — the shared selection bundle, plus `repeat_counts`: one entry per
surviving candidate, in candidate-id order, the reference at index 0.

**And one producer of that integer.** `repeat_count_of_bases` is now the only floor division in
this path: the ladder keys its rungs through it, `select_ssr` fills `repeat_counts` through it, and
the refused-tract early return uses it too rather than writing the division out a second time. The
genotype prior's `fill_ssr_seed` takes exactly this slice and must not recompute it — a second
producer only has to disagree once for a candidate's prior mass to land on the wrong length, and
nothing downstream could see it.

## The mode is not returned, and D2 is why

Arch §2.3 lists `modal_repeat_count` on this type. **Arch §4 and §5 then retract its reason**: on
2026-08-27 the prior's seed was re-indexed by offset from the **reference** tract length, which
every locus already knows, so the cohort's commonest length stopped being an input to it — *"nothing
consumes it now"*.

**D2 removed the last consumer inside this module.** The periodicity grid was the one remaining
reader of the ladder's mode, and it is now anchored on the reference tract's length as well. So
returning the mode would be a field with no reader anywhere.

It is therefore not on `SsrLocusSelection`. **The ladder still computes it**, so restoring it is one
line, and the architecture's own open question — whether selection should carry the mode at all —
is left where it is rather than answered by shipping a field nobody reads. **This is a departure
from the plan's D3, which says to return both**, and it is recorded here and in the type's own
documentation.

## Tests — 5 new

| test | what it pins |
|---|---|
| `the_repeat_counts_run_parallel_to_the_surviving_candidates` | parallel to the **survivors**, not to the merge's table, with the reference at 0 |
| `a_part_copy_carries_the_floored_count_the_ladder_keyed_it_by` | 7 bases → 3 copies, 9 → 4; and the tract is periodic because the two are one whole copy apart |
| `a_refused_tract_returns_one_repeat_count_for_the_reference` | a `NotPeriodic` locus is still indexed for the prior |
| `a_homopolymers_repeat_counts_are_its_lengths` | period 1, where a wrong period is invisible in the ladder |
| `every_candidates_count_is_the_rung_the_ladder_put_it_on` | the coupling itself: every survivor walked back to its merge index and compared against the ladder's own rung key |

## What the mutations found

Three deliberate defects, all caught:

| mutation | outcome |
|---|---|
| the counts taken over the merge's table rather than the survivors | caught — `the_repeat_counts_run_parallel_to_the_surviving_candidates` |
| the count rounds up instead of down | caught — 2 tests, one of them the ladder's own from B1 |
| a refused tract returns no counts at all | caught — `a_refused_tract_returns_one_repeat_count_for_the_reference` |

The second is the one that matters for the "one producer" claim: because the ladder and the prior's
slice now come from the same function, a single mutation to it breaks tests in both, which is
exactly the coupling the field exists to guarantee.

## Checkpoint D

**`select_ssr` is complete and proven on hand-built loci**, including the `NotPeriodic` and
reference-only outcomes. The four passes in order: the periodicity verdict, the fold, the ladder,
nomination with its `±1` rescue and cohort union, and the admission of the sequences on a promoted
rung with the cap and the leftover. Nothing outside the module calls it; the driver's dispatch to
it is the STR loop plan's own Milestone C, and the differential against production and HG002 is
Milestone E.

## Validation

All in the container (`./scripts/dev.sh`):

- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **5,980 passed, 0 failed, 14 ignored**;
  `ng::calling::allele_candidates` at **158**, from 153 at D2 and **93 before this work began**.
- `cargo doc --no-deps` — 26 unresolved-link errors, unchanged from the pre-change tree, none in
  these files.
