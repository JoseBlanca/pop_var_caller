# B3 — the allele count is one byte, and the cap refuses what a byte cannot hold

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md) step B3.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§2.1 and §2.2; [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.2 and §3.
**Date:** 2026-08-14.

---

## 1. What changed

**`AlleleObservation::reads` is a `u8` where it was a `u32`.** It is the count of a kept
position's reads that showed one non-reference allele, and it is the field the sparse list is
made of.

**`DepthCap` refuses a value above 255, and that is what makes the byte safe.** A count cannot
exceed the depth of the position it sits at, and the walk thins a position's counts by the
same ratio it thins its depth — so the cap bounds the count exactly. Nothing else stopped the
two drifting apart: the cap is a run-time value, and one set above 255 would have let the
counts saturate silently while the depth field beside them said otherwise.

The refusal is a **`const fn` constructor with an assertion**, `DepthCap::new`, not a fallible
one. The tuple field is now private, so the only ways to make a cap are that constructor and
`DepthCap::MAX`. Two consequences, both deliberate:

- **a cap written as a constant is refused while the run is being built.** The walk harness
  holds `const DEPTH_CAP: DepthCap = DepthCap::new(124);`, and a constant above 255 there is a
  compile error rather than a panic thirty seconds into a genome. That is the same move the
  ladder already makes with its five-bit `const _: () = assert!(…)`;
- **there is no longer an "uncapped" cap.** Two sites answered `DepthCap(u32::MAX)` for a
  cohort with no samples in it; they now answer `DepthCap::MAX`, the widest a byte holds. Both
  clamp nothing, because both are reached only when there is nothing to clamp.

**`thin_to_cap` returns the byte the record stores** rather than a `u32` a caller must narrow.
It panics rather than truncating if a surviving count will not fit — which would mean an
allele carried more reads than its position held. A truncated count is a wrong allele
fraction with no symptom; the panic names the count, the depth and the cap.

## 2. What moved on real reads: the entries are identical and the list is a third smaller

Both oracles, before and after, at four container CPUs — the 8-accession tomato subset (2.4 to
30.6 reads a position) and the GIAB trio (a few hundred).

**No fitted number moved, on either cohort. Not one line of either fit's output differs.** The
entire diff of both runs is the size of the sparse list and the per-sample total that contains
it:

| | entries | before | after | bytes an entry |
|---|---:|---:|---:|---:|
| HG002 | 29,798 | 0.358 MB | **0.238 MB** | 12.0 → **8.0** |
| HG003 | 32,762 | 0.393 MB | **0.262 MB** | 12.0 → **8.0** |
| HG004 | 47,552 | 0.571 MB | **0.380 MB** | 12.0 → **8.0** |
| tomato, deepest of the eight | 7,550 | 0.091 MB | **0.060 MB** | 12.0 → **8.0** |
| tomato, shallowest of the eight | 422 | 0.005 MB | **0.003 MB** | 12.0 → **8.0** |

**The entry count is unchanged everywhere**, which is the assertion that matters: the same
positions carry the same alleles at the same counts, and only the width they are held in
changed.

**Eight bytes and not the six the specification prices.** Spec §2.2 says an entry goes "from
about nine bytes to six" — that is the packed wire format. In memory the struct is a `u32`
index, a one-byte allele and a one-byte count, and Rust aligns it to the `u32`: six bytes of
content in eight. So the measured saving is **a third of the sparse list**, not the 44% the
byte format will give when the census is written to a file (plan 2). Both figures are about
the same change; they are the same entry counted twice, in memory and on the wire.

**What it is worth where.** The trio's HG004 holds 47,552 entries over 59,737 kept positions —
0.19 MB saved on a 0.62 MB sample. At the two million positions and 100× the specification
prices, the same rate is about **1 MB a read group**, and there the sparse list is the larger
half of the generic record rather than a rounding error.

## 3. Tests

Three added to `joint::census`, and one existing test widened:

| test | what it pins |
|---|---|
| `a_cap_a_byte_cannot_hold_is_refused_at_construction` | `DepthCap::new(256)` panics, naming the reason |
| `the_largest_cap_a_byte_holds_is_accepted` | the boundary from the other side — 255 is a cap, and `DepthCap::MAX` is it |
| `a_position_at_the_cap_round_trips_at_the_widest_count_the_byte_holds` | 255 reads on one allele at a cap of 255: the code stands for 255, the count comes back 255, and the fold across read groups does not truncate it |
| `a_thinned_share_rounds_to_nearest_and_never_loses_the_last_read` | unchanged in what it asserts; `thin_to_cap` now takes a `DepthCap` and returns a `u8` |

**The five synthetic-draw fixtures now say what they assume.** `fit.rs`, `contamination.rs` and
three examples draw allele counts from a Poisson depth without thinning them to any cap, so
each one converts with an `expect` rather than a cast: a draw above 255 stops the fixture
instead of wrapping. None of them draws that deep today — the deepest default is 8 reads a
position — but the sweep example takes its depth from the command line.

## 4. One place the invariant is not structural, and what it does there

`thin_to_cap` leaves a count untouched when the position's depth is zero, which happens only
for an observation whose read group the writer was not told about: the depth loop skips such a
group, so its depth stays zero while its reads still arrive. That is a pre-existing
inconsistency — the group ends up with a sparse entry and a never-walked depth array — and B3
does not change it. What B3 changes is the failure: above 255 reads it now panics naming the
count rather than writing the count modulo 256. Neither oracle reaches it; the walk declares
its read groups from the alignment header the observations come from.

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo test --lib ng::parameter_estimation::joint::census` | `29 passed; 0 failed` (26 before) |
| `cargo test --lib` | `3,578 passed; 0 failed; 11 ignored` (3,575 before) |
| the 88-second tomato oracle | §2 — no fitted number differs |
| the 74-second trio oracle | §2 — no fitted number differs |

**The two red gates are the two that were red before this branch's first commit**, neither in
code this plan touches: `cargo clippy --all-targets --all-features -- -D warnings`, and
`cargo test --all-targets`, which panics in `benches/psp_writer_perf.rs:386` with an index out
of bounds.
