# Fit precision, step D1 — the fixed order a repeat-tract stratum takes its samples in

**Date:** 2026-10-02. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step D1.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §4.4 (the order). **Branch:** `fit-precision`.

## 1. What was built

Step D2 will fit each repeat-tract stratum of a large cohort on a subset of its samples: the first ones in a
fixed order, grown until the stratum's slippage level is measured precisely enough. This step builds the order.

`parameter_estimation::joint::sample_order` (a new module):

- **`sample_order(names, seed)`** returns the indices of `names` ranked by a hash of each name and `seed`,
  lowest first. The hash is the census's, xxh3 with a seed (`loci::hash_position`), but taken over the name's
  UTF-8 bytes directly: the census goes through Rust's `Hash` for `str`, which appends a byte the standard
  library has not settled on, and a toolchain upgrade could move the order (review Mi1).
- **A sample's place depends on its name alone**, so the order of the names is the same whatever order they
  arrive in. Two different names with equal hashes are ranked by name. A run never holds two samples of one name
  (`RunError::SampleAppearsTwice`).
- **`SAMPLE_ORDER_SEED`**, a fixed constant. Changing it changes which samples every large cohort's subsets hold.
- Nothing calls it yet; D2 does, from `run::census_fit`, where the cohort's sample names are in hand.

No fitted number moves.

## 2. Tests

| test | what it shows |
|---|---|
| `the_order_does_not_depend_on_the_order_the_names_arrive_in` | every rotation and the reversal of thirty names give the same sequence of names |
| `every_sample_is_ranked_once` | a permutation, and not the arrival order; no names, no order |
| `the_seed_sets_the_order` | another seed, another order |
| `the_order_of_known_names_is_pinned` | eight names' order, recorded in the Linux dev container (arm64), which must hold on every platform; a change of hash or of how a name reaches it fails here |

- Review ([fit_precision_d1_2026-10-02.md](../reviews/fit_precision_d1_2026-10-02.md)): seven mutations, five
  killed; the two survivors change nothing a run can reach (a 64-bit hash collision; duplicate names).
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean. Module: 4 passed.
  Full suite: 5,015 passed, 3 failed (the pre-existing `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 8 ignored.

## 3. For step D2

The fit sees no sample names today. D2 computes the order once in `run::census_fit` from
`cohort.sample_names()` — the indices `SampleTractReads::sample` and the homozygote excess use — and hands the
order to the fit, with each sample's rank for testing whether a row's sample is in the subset.
