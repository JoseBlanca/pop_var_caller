# ng — read filtering in stages, D1: the first filter runs with no reference

**Date:** 2026-08-04 · **Branch:** `ng-generic-perf` · **Base:** `f5630f8` (Checkpoint C's fix)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **D1** —
the first of Milestone D, *the two tests output identity cannot see*.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §5, §8 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.2.
**Review:** [`ng_read_filtering_stages_d1_2026-08-04.md`](../reviews/ng_read_filtering_stages_d1_2026-08-04.md)
(4 Major / 8 Minor / 8 nits) · **Fixes:**
[`fixes_applied_2026-08-04_v1.md`](../reviews/fixes_applied_2026-08-04_v1.md).

---

## 1. Plan

> The first filter runs **with no reference at all** — construct it and drive it without a
> `RawRefSeq` in scope. Untested, spec §5's capability quietly stops being true.

Checkpoint D's mutation: *give the first filter a reference requirement and D1 must not compile.*

## 2. What shipped

`src/ng/read/reference_free_first_filter.rs` — a `#[cfg(test)]` sibling file under `read/`,
declared from `read/mod.rs` the way `left_align_parity` is. It holds a signature coercion, a config
built field by field, and one test that drives a whole reference-free pass: reader → narrowing →
first filter → tally.

**No production code changed.** The whole step is `#[cfg(test)]`.

## 3. It was built twice, and the first build was worthless

The first build followed the plan's sentence literally: a `#[cfg(test)]` module inside
`filtering.rs` whose mechanism was its import list — nothing in scope that could serve a reference
base. It passed, and the named mutation broke it.

**Three reviewers measured its unique detection power at zero, and one of them broke spec §5's
property with it green.** The review report has the detail; the two findings that decided the
rebuild:

- The checkpoint mutation was already caught by `cursor.rs:706` — **production code**, so
  `cargo build` fails on it with the module deleted. The module's alarm was a fourth copy of an
  alarm already ringing.
- `ReadFilterConfig::default()` launders a reference added *to the config*. A reviewer put a working
  `&dyn RawRefSeq` on the config behind `Default`, had the first filter read from it, and got
  2,867 passed, 0 failed with the module green. Filtering on flag and MAPQ now genuinely needed a
  reference, and the one test written to notice did not.

**Deviation from the plan, recorded.** The plan's D1 says *"construct it and drive it without a
`RawRefSeq` in scope"*. What shipped keeps that intent and drops the mechanism: scope is not what
pins the property, because a scope is a convention and `use super::*` erases it silently (a
reviewer ran that and everything stayed green). The property is pinned by a **signature coercion**
and by **how the config is constructed**, and the reference-free scope is where those two live so
that neither can be repaired by reaching for something already to hand. Intent delivered more
strongly; mechanism changed. Absorbed rather than escalated because it changes no design and stays
inside the step — but flagged at Checkpoint D so the owner can overrule.

## 4. The one thing here nothing else catches, measured

The coercion is the right *statement* of the property but is not a unique alarm: the same mutation
fails `cursor.rs` and `filtering.rs`'s `pre` helper with `E0061`. Nor is the sibling-module
placement unique — the visibility mutation breaks `cursor.rs` too.

The config is. Reproducing the reviewer's escape and then repairing it the way its author would —
supplying the new field in `Default` **and** in `filtering.rs`'s `post_config`, whose module already
holds an `InMemoryRefSeq`:

```text
error[E0063]: missing field `reference` in initializer of `filtering::ReadFilterConfig`
   --> src/ng/read/reference_free_first_filter.rs:100:5
error: could not compile `pop_var_caller` (lib test) due to 1 previous error
```

One error, in the new file, nowhere else.

## 5. Mutations run

| mutation | result |
|---|---|
| `_reference: &impl RawRefSeq` on `verdict_on_raw_read` | the coercion fails with **`E0308`** — the only type-level error; the three call sites give `E0061` |
| a field on `ReadFilterConfig`, `Default` **not** supplying it | fails at `Default`, at `post_config`, and here — three sites |
| a field on `ReadFilterConfig`, `Default` **and** `post_config` supplying it | **fails here alone** — §4 |
| `verdict_on_raw_read` narrowed back to private | fails here (`E0603`) and at `cursor.rs:75` — co-detected, not unique |
| the secondary record's flag set to `SUPPLEMENTARY` | the pass test fails; the exhaustive `ReadFilterCounts` equality is what catches it |

## 6. Two tests I wrote and then deleted

Both were in the rebuilt file before the measurement, and both were padding:

- **`the_flag_bits_the_pass_sets_are_the_bits_the_filter_reads`** — claimed to guard against
  noodles' `Flags::*` drifting from production's `FLAG_*`. Setting the secondary record's flag to
  `SUPPLEMENTARY` fails the pass test and leaves this one **green**: the pass's exhaustive count
  equality already catches every misrouting.
- **`a_reference_free_pass_still_knows_which_read_group_a_drop_came_from`** — its justification was
  that a read-group resolution needing a reference would leave the pass green. False: the pass
  drives the same narrowing and would not compile either. Its one assertion moved into the pass's
  loop.

## 7. Tests

| test | what it pins |
|---|---|
| `reference_free_first_filter::a_pass_over_raw_reads_filters_and_charges_every_drop_with_no_reference` | a whole pass — reader, narrowing, first filter, tally — over six records, keeping one and charging each of five flagged drops, from a scope that names no reference and constructs the config field by field |
| `const _FIRST_FILTER_TAKES_NO_REFERENCE` | the first filter's whole signature, in one line no import can repair |

**Suite: 2,866 → 2,867 (+1).** Fully accounted: one test added, and the first build's one test
replaced rather than kept.

## 8. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,867 passed**, 0 failed, 5 ignored |
| `cargo test --examples` | 52 passed, 0 failed |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` against the `f5630f8` baseline | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709`, every counter line identical |

The change is entirely `#[cfg(test)]` and cannot reach a release binary; the dumps were run because
the plan checks every step against them.

## 9. Open, for Checkpoint D

- **Spec §5 overclaims and should be amended.** *"The reference stops being a precondition for
  filtering at all"* is true of the filter, the reader and the narrowing — and false of
  `AlignmentCursor<R: RawRefSeq>` and `AlignmentFile::cursor`, whose bound is unconditional. So §5's
  three named callers still cannot reach a *file's* reads without producing a reference. Recorded in
  the new module's doc; no design doc edited.
- **The `const` witness is `#[cfg(test)]`.** Beside `verdict_on_raw_read` it would make
  `cargo build` carry the signature check too. That is a line in a shipping file, and Milestone D is
  scoped as tests only.
- **Plan D1's "construct it"** predates C2 making the first filter a free function; there is nothing
  left to construct.
