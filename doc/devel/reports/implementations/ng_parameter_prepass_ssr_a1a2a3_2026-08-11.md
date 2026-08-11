# ng step 4, the STR path — A1+A2+A3: the module, the stratum, and the two widths

*Implementation report, 2026-08-11. Steps A1, A2 and A3 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), run under the
plan-driven-implementation skill. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) and
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md).*

## Plan

Three plan steps in one loop, because A1 is scaffolding and would otherwise be a commit of four
module doc comments. Named here rather than merged silently, and the same bundling the sibling
generic path used for its own A1+A2+A3.

- **A1** — `src/ng/parameter_estimation/ssr/` with its four files, wired into step 4's surface.
- **A2** — `SsrPeriod` in the shared vocabulary; `RepeatCount` and `Stratum` beside the path that
  fits per stratum.
- **A3** — the offset scalars, the bucket rule, and the two widths that are not the same width.

No logic beyond the bucket rule and the allele support; the mathematics is Milestone D.

## Assumptions and deviations

Seven, none of which changes the design. Each is a choice the plan or the architecture left open, or
a place where following the architecture's sketch to the letter would have made an illegal state
representable.

1. **`RepeatCount` and `Stratum` live in `parameter_estimation/ssr/mod.rs`, not in `types.rs`.**
   Architecture §2.1 lists all three types under a heading that says "extend `types.rs`", but its
   own prose immediately says *"`SsrPeriod` is STR domain vocabulary with consumers in steps 6 and
   7 … The rest are step-4's own"*, and [`module_layout.md`](../../ng/arch/module_layout.md) names
   only `SsrPeriod` as expected in the shared file. Took the narrower reading: the one type another
   step will name goes in the shared vocabulary, the two that only step 4 speaks stay with step 4.
   Moving them later is a re-export, not a redesign.
2. **`OffsetBucket`'s value is private, where the architecture sketches `pub u8`.** A bucket indexes
   an entry's counts, so `OffsetBucket(200)` would be an out-of-bounds index that compiles.
   [`bucket_of`](../../../../src/ng/parameter_estimation/ssr/mod.rs) is the only constructor and it
   is total. `RepeatCount` and `WholeRepeatOffset` keep their public fields, because every value of
   each is a legal thing to observe.
3. **The offset vocabulary lives in `mod.rs`, not in `locus_offsets.rs`.** The architecture's file
   table gives `locus_offsets.rs` one job — turning a locus into an entry — which is Milestone C.
   The offsets themselves are spoken by all four files, so they sit where the module's other
   vocabulary does.
4. **`allele_support` takes a repeat count and reads `ALLELE_OFFSET_LIMIT` from the constant**,
   where the measurement harness's version takes the limit as an argument. Nothing in the shipped
   path varies it; the harness keeps its own copy and its own sweep.
5. **`Motif::ssr_period()` returns an `SsrPeriod`, not a `Result`.** `Motif::new` already rejects a
   length outside `1..=MAX_MOTIF_LEN`, so the conversion cannot fail; the `expect` states that
   invariant rather than guarding a fallible path.
6. **Four `Display` impls the architecture does not list** — on `SsrPeriod`, `RepeatCount`,
   `Stratum` and `WholeRepeatOffset`. A6's two error messages and Milestone E's summary both name a
   stratum, and the alternative is each of them formatting the pair itself. `WholeRepeatOffset`
   renders signed, because `+2` and `-2` are different observations and a bare `2` says neither.
7. **The three sibling files carry a module doc and no empty `#[cfg(test)] mod tests`**, where A1
   asks for the block. An empty test module is a placeholder that asserts nothing; the tests land
   with the code they test, in B, C and D.

## Changes made

**[`src/ng/types.rs`](../../../../src/ng/types.rs)** — shared vocabulary:

- `SsrPeriod`, a checked `u8` rejecting zero and anything past `MAX_MOTIF_LEN`. Zero is the one that
  has to be unrepresentable: a tract's length becomes a repeat count by dividing by its period.
- `DomainError::SsrPeriod(u8)`, whose message quotes the range from `MAX_MOTIF_LEN` rather than a
  literal.
- `Motif::ssr_period()`, beside the existing `period()` rather than in place of it, so the three
  modules calling that accessor are untouched.

**[`src/ng/parameter_estimation/ssr/`](../../../../src/ng/parameter_estimation/ssr/)** — the new
folder, four files, wired into `parameter_estimation/mod.rs`:

- `mod.rs` — `RepeatCount`, `Stratum` (ordered by period then repeat count, which is the order the
  monotonicity rule walks), `WholeRepeatOffset`, `OffsetBucket`, `bucket_of`, `allele_support`, and
  the three constants: `OFFSET_HALF_RANGE = 4`, `OFFSET_BUCKETS = 9`, `ALLELE_OFFSET_LIMIT = 6`.
- `locus_offsets.rs`, `stratum_table.rs`, `slippage.rs` — module docs only, each saying what lands
  there and in which milestone.

**Two widths, and the doc comments carry which one is load-bearing.** The recorded offset range can
be narrow — with an end bucket scored by summing over what it absorbs, ±1 still returns the slippage
level to within 0.05% against alleles reaching ±3. The allele support is what decides the answer,
and its cost curve is a cliff rather than a slope: 2.5% of loci outside costs 0.1% of the level,
19.3% outside costs +499% with the direction asymmetry destroyed.

## Tests added

Fourteen, all in the two files' own `#[cfg(test)]` blocks.

In `types.rs` (4): every period in scope round-trips; zero, seven and 255 are rejected by name;
`Motif::period` and `Motif::ssr_period` agree at all six lengths — the live path mints a period only
through a motif, so a disagreement there moves every stratum a locus is filed under; a period
renders as the bare number.

In `ssr/mod.rs` (10):

- **strata sort by period and then by repeat count** — the ordering the monotonicity walk depends
  on. Sorting by repeat count first would compare a dinucleotide's fit against a hexamer's.
- **a stratum names itself** as "period 2, 6 repeats", and cannot be built at a period of zero.
- **the bucket count is read against the range it derives from**, not against a literal, so moving
  one constant without the other fails here.
- **every offset inside the range gets its own bucket, in order**, with the reference length in the
  middle.
- **offsets past the range saturate**, at −5, −12, −40 and `i8::MIN` and at the positive extremes
  (5, 12, 40 and `i8::MAX`): a read 40 repeats short is a real observation and has to land
  somewhere.
- **exactly the two ends are saturating** — the property the scoring rule of Milestone D must know,
  where plugging in the edge instead of summing over what it absorbs costs +33% of the level.
- **the allele support clips below the reference and not above**: 10 lengths at 3 repeats, 13 at 6
  and 13 at 20, which is the arithmetic behind 55, 91 and 91 genotypes.
- **two boundary cases the copy floors make unreachable but the arithmetic must survive**: a repeat
  count of `u32::MAX`, which does not fit the offset type, and a tract of no repeats at all.

## Validation

All commands run in the dev container (`./scripts/dev.sh`), on this branch.

*(As first written, before the review. The numbers below are corrected — four of them were wrong,
and every one was a claim about this work's own gate run rather than a figure quoted from a design
document. What the review changed in the code is in the fixes-applied report.)*

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --all-features` | **3,386 passed, 0 failed, 10 ignored**, across 12 test binaries — 3,306 of them the library's, which is the single line an earlier draft of this table quoted as the whole run |
| `cargo test --lib parameter_estimation::ssr` | 10 passed |
| `cargo test --lib ng::types` | 33 passed |

**Two gates are red on this branch and neither is this work's.** `cargo clippy --all-targets` fails
with **8** errors across four files under `examples/` — `ng_str_stutter_rate.rs`,
`ng_generic_loci_dump.rs`, `shared/stutter_model.rs`, `shared/stutter_table.rs` — though a given run
may print only 7, because `-D warnings` aborts the build at the first failing example and whether
`ng_generic_loci_dump` finishes checking varies between runs. The lints are `needless_range_loop`
(×2), `needless_borrows_for_generic_args` (×2), `manual_is_multiple_of`, `manual_checked_div`,
`collapsible_if` and `nonminimal_bool` — not all of one kind, as an earlier draft of this paragraph
said. **All 8 reproduce at `b96fdcd0` with this patch not applied**, which is what makes them
inherited rather than introduced; none of those files is touched here. `cargo doc --no-deps --lib`
reports 13 unresolved intra-doc links, **none of them in the files this step wrote** — the one
inside `types.rs` is at line 513, in `Motif`'s pre-existing doc comment.

## Tradeoffs and follow-ups

- **`allele_support` allocates a `Vec` per call.** It is called once per stratum per fit, against a
  search that re-walks a stratum's whole table per candidate, so it is not on any hot path. If
  Milestone D's inner loop ever wants it allocation-free, the shape is a `RangeInclusive<i8>`
  mapped rather than collected.
- **Nothing yet reads any of this**, which is the milestone's own rule: the mathematics is proven
  against the exact-bias harness before a locus is read.
- **A5 will need `Estimate<T>` and `Provenance`** from the step's surface, both of which exist.
