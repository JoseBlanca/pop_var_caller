# B2 — the tract path's settings, and one sample's reads counted onto the ladder

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone B step B2, and **Checkpoint B**. **Design:**
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §2.2, §3.1;
[`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §5, §7, §12 Q1–Q2.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs),
with one buffer added to the shared scratch in
[`mod.rs`](../../../../src/ng/calling/allele_candidates/mod.rs).

---

## What landed

**`SsrSelectionConfig`, and it has no `Default`.** Its three fields are the shared support rule
and cap, the periodicity gate's share, and the run's ploidy. `at_ploidy(ploidy)` is the only
constructor of the defaults, so there is no way to reach a configuration without naming how many
copies of the genome a sample carries — which is the design's point (arch §2.2: "from
`FrozenParameters::ploidy`, never a constant here"), since a polyploid region is in scope and
ploidy is the one thing that changes how many rungs a sample promotes.

**Two constants, both settled elsewhere and both soft.**
`DEFAULT_MAX_CANDIDATE_ALLELES_SSR = 32` against the ordinary path's six — the owner's decision of
2026-08-24, unmeasured and bounded by cost: 528 diploid genotypes a sample a locus against six's
21. `DEFAULT_MAX_OFF_GRID_SHARE = 0.10` — production's `max_out_of_frame_frac`
([`candidate_set.rs:85`](../../../../src/ssr/cohort/candidate_set.rs)), inherited and never
measured by them or by us, and every number spec §4.1 and §5 report was taken with that gate
switched off entirely.

**`fill_sample_reads_per_rung` — one sample's spanning reads counted onto the ladder's rungs.**
It comes back one entry per rung, in the ladder's own shortest-first order, so nomination can walk
rungs by index. The read-group rows are pooled through the *same* `one_run_per_allele` the ordinary
path's fold uses — a read is a read whichever lane produced it, and asking the rule of each row
separately would be a stricter rule applied to exactly the samples carrying more than one library.

## Two departures from the architecture sketch, and one from the spec

**The off-grid share is a validated newtype, where arch §2.2 writes a bare `f64`.** The failure it
guards is silent in both directions: above one, no sample is ever non-periodic and the
`NotPeriodic` verdict can never be reached; below zero, every sample is, and every repeat tract in
the run comes back as the reference alone. That is the same class of defect this module's own A1
review found in `MinAltReadShare` — "a negative share does not crash, it deletes half the rule" —
so `MaxOffGridShare` follows that type's shape and delegates its range check to it, rather than
writing a second spelling of *is this a fraction of one*.

**The ladder gained a reverse index**, `rung_of_table_index`, filled during the B1 build. A sample's
support rows name merge-table indices in the merge's order; the ladder's own buffer is in rung
order. Without the reverse direction a per-sample fold would have to re-derive each sequence's
repeat count from its bases and search the rungs for it — a second producer of the one integer the
ladder exists to have one producer of.

**And the support share is the ordinary path's 10 in 100, where spec §5 writes 5.** This is the
one departure that changes a number, so it is set out in full below.

### Why 10 and not the spec's 5

**The spec's own reason for its number is what points at the other one.** §5 sets the tract share
to 5 in 100 *"so that one number governs both paths"*, citing the ordinary path's value as it stood
when that sentence was written. **That value moved to 10 in 100 the same day**, by the owner's
decision (commit `2a7170f2`, 2026-08-24), taken deliberately against what recall alone would say:
recall says 5 is free and 10 costs two true alleles on the trio at 300×, and it ships at 10 because
the candidate count is the other side of the trade and nothing had measured it. Setting 5 on the
tract path now would create exactly the second number §5's sentence exists to avoid.

**What it costs, from spec §5's own sweep at 300× on HG002.** On the class this path was designed
for — a heterozygote whose two copies are the same length spelled differently, 296 of HG002's 695
heterozygous tracts — the two rules give:

| share | both spellings offered | candidate sequences per tract |
|---|---|---|
| 5 in 100 | 86.1% | 1.26 |
| 10 in 100 | 85.8% | 1.22 |

**Three tracts in a thousand, and fewer candidates.** At 30× and below the two are the same rule,
because the floor of two reads decides — the share is inert below 21 compared reads a sample.

**It shifts what Milestone E checks against.** That milestone reproduces spec §4.1's and §5's HG002
numbers through the shipped module and says "a difference is a defect and must be traced". The
targets are now **85.8% and 1.22**, not 86.1% and 1.26 — recorded here so the difference reads as
this decision rather than as the defect being hunted.

## Tests — 11 new

| test | what it pins |
|---|---|
| `the_tract_cap_is_thirty_two_where_the_ordinary_paths_is_six` | both caps in one assertion, and that the config carries the tract one |
| `the_support_rule_is_the_ordinary_paths_share_and_floor` | the share above, as a line someone has to edit rather than a drifting default |
| `the_ploidy_is_the_callers` | a triploid run reaches the config as three |
| `the_off_grid_share_refuses_anything_that_is_not_a_fraction_of_one` | negative, above one, `NaN`, `+∞`; and the default at 0.10 |
| `a_samples_reads_land_on_the_rungs_its_sequences_sit_on` | a rung the sample showed nothing at is a zero **in place**, not a missing entry |
| `two_spellings_of_one_length_add_on_their_shared_rung` | 4 and 5 reads at five repeats are 9 at five repeats |
| `one_samples_two_read_groups_land_on_one_rung` | 3 from one lane and 4 from the other are 7 |
| `a_sample_with_only_partials_counts_at_no_rung` | 40 partial reads name no length |
| `refilling_for_a_second_sample_leaves_none_of_the_firsts_counts` | the buffer is one per worker |
| `the_histogram_totals_the_samples_compared_reads` | the numerator's rungs and the denominator's total have one producer |
| `a_row_naming_an_allele_the_tract_does_not_hold_is_refused` | the merge bug spec §8 makes an assertion |

## What the mutations found

Three deliberate defects, applied to the tree, run, and copied back:

| mutation | outcome |
|---|---|
| the rung's count overwrites instead of adding | caught — `two_spellings_of_one_length_add_on_their_shared_rung`, and **only** that one |
| the histogram is not emptied between samples | caught — `refilling_for_a_second_sample_leaves_none_of_the_firsts_counts` |
| every allele recorded as sitting on the shortest rung | caught — 4 tests |

**The first mutation is why `two_spellings_of_one_length_add_on_their_shared_rung` exists.** Every
other histogram fixture puts each of a sample's alleles on a rung of its own, or reaches one rung
through two read groups — which `one_run_per_allele` has already pooled into a single row group
before the accumulation runs. So an overwrite in place of the add gives the identical answer on all
of them. The one input that separates the two is a sample carrying **two distinct sequences of the
same length**, which is the interrupted repeat this whole path was designed to offer both spellings
of. The test was written because the mutation was predicted to survive, and it did until the test
existed.

## Validation

All in the container (`./scripts/dev.sh`), on the tree as committed:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::calling::allele_candidates` — **118 passed**, from 107 at B1 and 93 before
  this milestone.
- `cargo test --all-targets --all-features` — library **5,933 passed, 0 failed, 14 ignored**, every
  other target green; the run exits 101 on the pre-existing index-out-of-bounds at
  `benches/psp_writer_perf.rs:386`, in production's psp writer bench.
- `cargo doc --no-deps` — 26 unresolved-link errors, the same 26 as on the pre-change tree, none in
  these files.
