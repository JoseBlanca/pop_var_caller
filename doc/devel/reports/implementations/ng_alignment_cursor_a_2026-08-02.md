# ng — the alignment cursor, Milestone A: implementation report

*Plan: [impl_plan/alignment_cursor.md](../../ng/impl_plan/alignment_cursor.md).
Design authority: [spec](../../ng/spec/alignment_cursor.md) and
[arch](../../ng/arch/alignment_cursor.md). Evidence:
[perf_ng-generic-pileup_2026-07-31.md](../reviews/perf_ng-generic-pileup_2026-07-31.md).
Branch `ng-generic-perf`. One section per plan step; the milestone's steps land as
separate commits and this file grows with them.*

**Milestone A changes no behaviour.** It commits the instrument every later milestone is
verified against, and adds the types the cursor will be built from. The probe's real-data
numbers must be **identical** at the end of it.

---

## A1 — the probe, its tests, and the fixture builders it shares with the bench

### Plan

`examples/ng_generic_walk_probe.rs` existed in the working tree, uncommitted. A1 lands it
**with a test module**, and with the fixture builders factored so
`benches/ng_generic_pileup_perf.rs` shares them rather than carrying its own copy of the
synthetic contig and the BAM writer.

Three moves:

1. Lift the bench's fixture builders into `examples/shared/synthetic_alignment.rs`,
   parameterised by a `SyntheticGeometry` instead of the bench's file-level `SPAN` and
   `READ_LEN` constants. Included by `#[path = …] mod`, which is this repo's existing way of
   sharing a body between example targets (`examples/shared/witness_side.rs`).
2. Point the bench at it, deleting its copy.
3. Give the probe a test module over that same fixture, and register the example with
   `test = true` so the tests actually run.

### Assumptions

- **Where shared code lives.** The plan does not say. Rust gives no import path from one
  example to another or from a bench to an example, so the choices were a `pub` module in
  the library or a `#[path]`-included file. The library is the wrong home — a synthetic BAM
  writer is not part of what `pop_var_caller` offers — and `examples/shared/` is the
  precedent already in the tree.
- **What the probe's tests should pin.** Spec §11 names the real anchor (HG002 chr21) and
  says the unit suite is not the bar. So the tests here do not attempt to be that bar; they
  pin the two properties the cursor work will be judged against, at a size that runs in the
  suite: the walk covers its span, and **the answer does not depend on how the span is cut
  into regions**.

### Deviations from the plan or the file as it stood

**One, and it changes a decision recorded in the probe's own comments.** The file carried a
`compile_error!` making both `alloc-mimalloc` and `dhat-heap` a build failure, argued as
"the exclusion has to be stated, not implied". Committing that turns CI red: this project
lints and tests with `--all-features`
([ci.yml:35,47](../../../../.github/workflows/ci.yml#L35)), which enables both.

The guard's *intent* is kept and its instrument moved: mimalloc is now gated on
`not(feature = "dhat-heap")` (matching `examples/dhat_ng_open_files.rs:58`), and
`allocator_is_ambiguous()` reports the contradiction and `main` refuses the run. The thing
worth stopping was never a build — it was a **measurement** that reports mimalloc's name and
dhat's behaviour — and that is a run. The check uses `cfg!` rather than `#[cfg]` so the
condition is type-checked in every build.

*Raised at Checkpoint A rather than absorbed silently, because it edits a decision the perf
review wrote down.*

Two smaller ones, absorbed:

- `parse_env_u64` keeps its shape but delegates the rule to a pure `parse_count(name,
  value)`. The zero-rejection rule was untestable otherwise: `std::env::set_var` is `unsafe`
  in edition 2024 and this crate forbids `unsafe`.
- `SyntheticSample` does **not** carry the geometry it was built from. Both callers pass one
  in and still hold it; a second copy is a second thing that can be read after the first has
  moved on.

### Changes made

| file | change |
|---|---|
| [examples/shared/synthetic_alignment.rs](../../../../examples/shared/synthetic_alignment.rs) | **new.** `SyntheticGeometry` (span, read length, coverage, and the read count they imply) and `SyntheticSample::build` (seeded contig → FASTA, tiling reads with one mismatch each → coordinate-sorted BAM). The builders below `build` are private: both consumers use the same two-item surface. |
| [benches/ng_generic_pileup_perf.rs](../../../../benches/ng_generic_pileup_perf.rs) | 170 lines deleted; includes the shared module. `check_the_walk_covers_the_span` now takes the geometry instead of reading file-level constants, so its two assertions name the numbers they check. Its behaviour is unchanged. |
| [examples/ng_generic_walk_probe.rs](../../../../examples/ng_generic_walk_probe.rs) | **committed**, plus: the allocator guard above, the `parse_count` split, and a test module (nine tests as written, nineteen after review). |
| [Cargo.toml](../../../../Cargo.toml) | `[[example]] ng_generic_walk_probe, test = true`. |

**A cargo behaviour worth recording, and a measurement of it I got wrong.** An example's
`#[test]`s do not run under a plain `cargo test` unless the target says `test = true`, which
is what this entry buys — and it matters for the same reason the Milestone D witness review
gave: a guard nobody runs is not a guard. My first version of this paragraph also claimed
`--all-targets` *builds* examples without running their tests, and said it had been checked.
It does run them. The measurement behind the claim had been aborted by a **pre-existing panic
in an unrelated criterion bench** (`benches/psp_writer_perf.rs:386`) that ends an
`--all-targets` run before it reaches any example; three review sub-agents measured the truth
independently. The manifest comment now states what holds.

### Tests added

Nine as first written; **nineteen after the review**, which found the nine inadequate to
their own stated job — see [the fix report](ng_alignment_cursor_a_fixes_2026-08-02.md). The
original nine:

| test | what it would catch |
|---|---|
| `the_whole_contig_walk_yields_one_locus_per_reference_position` | a walk that generates less than it claims — which reads as a *fast* run, not a wrong one. Pins depth too: one read per position would still cover the span. |
| `cutting_the_span_into_regions_does_not_change_what_the_walk_produces` | **the property the cursor must not break.** Same bases as one region, as five, as a dozen → same loci and same observations. The unit-scale form of spec §11.2. |
| `the_typed_region_walk_reaches_the_generator` | the default path (typed regions, not the whole-contig bound) never reaching the generic slot. |
| `the_locus_ceiling_stops_the_walk_exactly` | `PVC_PROBE_MAX_LOCI` stopping near a count instead of at it. |
| `splitting_a_generic_region_tiles_it_exactly` | a gap or an overlap between pieces at five chunk sizes, including chunk 1 and a chunk wider than the region. |
| `only_generic_regions_with_a_chunk_are_split` | the grain knob cutting up a region typed for another generator. |
| `a_failed_region_passes_through_the_splitter` | a failed region swallowed on its way past — which would un-call the rest of a contig silently. |
| `a_count_knob_rejects_zero_and_anything_that_is_not_a_count` | a knob's zero taken as a value; each knob's zero is a different silent no-op. |
| `the_report_prints_the_keys_the_recorded_baselines_name` | a renamed printed key. Spec §11's anchor is a **line of text**, so a rename breaks every comparison against every recorded run and nothing else would notice. |

### Validation

Real output, host-native (`cargo` is allow-listed; the container was not needed):

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --tests --examples --all-features` (CI's command) — `2770 passed; 0
  failed; 5 ignored`, plus the probe's `9 passed; 0 failed`.
- `cargo test --lib ng::` (debug) — `1471 passed; 0 failed; 2 ignored`, unchanged from the
  base commit. *Release is red on a clean tree for reasons predating this branch: four tests
  assert on `debug_assert!` messages that release compiles out.*
- `cargo bench --bench ng_generic_pileup_perf -- --test` — six cases, all `Success`, so the
  bench's own `check_the_walk_covers_the_span` assertions survived the fixture move.
- **The real-data anchor, re-run after the change**, HG002 30× chromosome 21:
  `loci=236081 observations=251786 reads_admitted=54709` — the recorded baseline, digit for
  digit.

### Tradeoffs and follow-ups

- The probe's test geometry is 5,000 bp at 10×. It cannot reach the region size the first
  retention attempt first diverged at (100 kb) — that check stays the real-data one, by
  design.
- `refuse_an_ambiguous_allocator` is a startup guard on a path CI never runs (`main`), so
  nothing tests it. Testing it would mean a subprocess run under a feature combination that
  exists only to be rejected.
