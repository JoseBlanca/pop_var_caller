# ng — read filtering in stages, D2: the conversion is asked for nothing when the first filter drops

**Date:** 2026-08-04 · **Branch:** `ng-generic-perf` · **Base:** `5320fe4` (D1)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **D2** —
the last step of the plan.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §2, §3, §6, §8 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.4.
**Review:** [`ng_read_filtering_stages_d2_2026-08-04.md`](../reviews/ng_read_filtering_stages_d2_2026-08-04.md)
(4 Major / 12 Minor / 12 nits) · **Fixes:**
[`fixes_applied_2026-08-04_v2.md`](../reviews/fixes_applied_2026-08-04_v2.md).

---

## 1. Plan, and why it could not be followed

> **D2.** The conversion is asked for nothing when a read fails the first filter. **Use a read
> that would fail to convert:** unmapped, with no alignment start. Filter #5 drops it before the
> conversion; if the conversion were hoisted it would raise a fatal decode error instead.

**That fixture cannot exist**, and the reason is structural.
`RegionRawAlignedReads::read_next` yields a record only after proving three things:

1. `reference_sequence_id() == Some(this contig)`
   ([`region_raw_aligned_reads.rs:189`](../../../../src/ng/read/input/region_raw_aligned_reads.rs#L189));
2. `overlaps()`, which returns `false` unless **both** `alignment_start()` and `alignment_end()`
   are `Some` ([`:259`](../../../../src/ng/read/input/region_raw_aligned_reads.rs#L259));
3. a read group, stamped by the CRAM arm or resolved to `RecordOwner::Mine`
   ([`:230`](../../../../src/ng/read/input/region_raw_aligned_reads.rs#L230)).

`NoodlesRawAlignedRead::decode` refuses exactly three things: no read group, no reference sequence
id, no alignment start ([`aligned_read.rs:258`](../../../../src/ng/read/aligned_read.rs#L258),
[`:134`](../../../../src/ng/read/aligned_read.rs#L134),
[`:140`](../../../../src/ng/read/aligned_read.rs#L140)). **The three guarantees and the three
refusals are the same three.** A read with no alignment start has no footprint, so it is discarded
a layer below the filters, uncounted — not charged to `unmapped` either.

So the plan's D2 would pass whether or not the conversion had been hoisted, which is the one thing
it exists to catch. **Measured rather than argued: with the conversion hoisted above
`verdict_on_raw_read`, `cargo test --lib` gives 2,866 passed, 0 failed.**

**The owner chose the replacement** (2026-08-04) from a written set of options, and also lifted the
bar on editing the design docs.

## 2. What shipped

**The step-1 loop's inline `self.buffer.decode()` became `self.convert_buffered_read()`**, a
private method that increments a counter and makes that same call. The counter is
**`CursorCounts::reads_converted`** — a real field on the struct that already exists to say what
the cursor did, beside `reads_decoded`.

The tests are two, plus two fixture helpers (`flagged_read_at`, `read_at_mapq`).

**The counter is inside the callee, and that is the entire mechanism.** Write the increment as its
own statement beside the call in `next_filtered_read` and a hoist that moves the call leaves it
behind — whereupon it counts the survivors and reports **exactly what the test expects**, so the
test passes while every read is being converted. I measured that before writing it down. Inside
`convert_buffered_read` there is no such shape: the count goes wherever the call goes.

**It ships rather than hiding behind `#[cfg(test)]`, and that changed during the step.** The first
build used a `#[cfg(test)]` thread-local with a reset protocol. Two reviewers independently built
the field version, measured identical detection power, and pointed out that the thread-local's
justification (*"a field would have to be threaded through every constructor"*) is refutable by one
grep — `AlignmentCursor` has exactly one struct literal, and it already holds a `CursorCounts`.
The field needs no reset, is per-cursor rather than thread-global, and is folded across a sample's
files for free by the exhaustive-destructure `AddAssign`. It also avoids introducing the crate's
first `#[cfg(test)]` static and first `#[cfg(test)]` statement inside a production function body.

**Not a test double** (spec §6, which the deleted `RecordSource` doubles cost four steps to be rid
of). The reader, the narrowing and the conversion are all the real ones; a real conversion is
counted as it happens. Nothing is stubbed and no layer is bypassed.

**A compile-time pin was possible and was not taken.** A reviewer built a zero-sized witness minted
by the first filter's `Keep` arm and required by the conversion, which turns the hoist into a build
error. Rejected for two reasons, both recorded on the method: spec §1 says the design "does not …
add a type", and the witness pins only half the property — no type can forbid someone writing a
second copy of the length rule and checking it early, which is what the second test catches.
Carried to Checkpoint D.

## 3. The property is stated as an accounting identity

**conversions = the reads that reached the second filter = kept + second-filter drops.**

Two tests, one for each direction, because the ordering can be broken both ways and neither shows
in any output:

| test | breaks if |
|---|---|
| `a_read_the_first_filter_drops_is_never_converted` | the conversion **rises** above the flag/MAPQ checks |
| `a_read_the_second_filter_drops_has_already_been_converted` | one of the second filter's checks **sinks** below the conversion |

The second is not symmetry for its own sake. Filter #7 compares a sequence length, which a raw
record can answer without being converted, so moving it down looks like a free saving — and it
changes no output at all. Spec §2 keeps #7, #9 and #8 together on the decoded read because writing
any of them against noodles' types is a second copy of a rule, which is what `filtering.rs` guards
against hardest.

The first test's script is **mixed** — two clean reads among six dropped ones, **one per first-filter
reason** — so the assertion says *which* reads were converted rather than only that the counter is
small. Covering all six rather than a sample is not padding: review converted just the `Unmapped`
and `LowMapq` drops and the whole suite stayed green, and #2 is the highest-volume drop on real
data.

## 4. Mutations run

Every one re-run against the final version:

| mutation | result |
|---|---|
| **the conversion hoisted above the first filter** (the plan's named mutation) | test 1 FAILED — **2,868 passed, 1 failed**, `left: 8, right: 2`: killed by that test alone |
| **filter #7 hoisted above the conversion** (a length check on the raw record, then `continue`) | test 2 FAILED — **2,868 passed, 1 failed**, killed by that test alone, and **no output moved** |
| **the increment removed** from `convert_buffered_read` | **both** tests FAILED — **2,867 passed, 2 failed** |
| **two drop reasons convert before rejecting** (`Unmapped`, `LowMapq`) | test 1 FAILED — `left: 4, right: 2`. **This one survived the three-reason fixture the step was first built with** |
| **the increment left at the call site while the call is hoisted** | **survives** — both tests pass. This is why the count lives in the callee, and the reason is not the one first written down: the abandoned increment reports the *expected* number, not zero |

Each test kills a distinct mutation on its own; removing the instrument kills both.

## 5. Speed

The counter **ships**, so this is a real per-conversion increment in the release build and a
measurement is owed. `ng_generic_walk_probe` on HG002 chr21, **six runs a side, one machine, one
session**, nothing else running:

| | runs (s) | median |
|---|---|---|
| before | 1.839, 1.846, 1.826, 1.827, 1.819, 1.826 | 1.8265 |
| after | 1.853, 1.848, 1.834, 1.840, 1.830, 1.838 | 1.8390 |

**Noise, not a cost.** The medians differ by **0.7 %** and the ranges overlap — the fastest *after*
run (1.830) beats the slowest *before* run (1.846) — which is the opposite of the separation
pattern B1 used to call its 1.4 % consistent. One `u64` increment per conversion, ~55k of them on
chromosome 21.

*(An earlier six-run pair measured the `#[cfg(test)]` thread-local version, which compiled out in
release; it is superseded and not reported. A single `seconds=2.481` reading taken with a
`cargo clippy --all-targets` build running alongside it is discarded, recorded here only so the
number is not mistaken for a measurement.)*

## 6. Tests

| test | what it pins |
|---|---|
| `a_read_the_first_filter_drops_is_never_converted` | of eight records — two clean and **one per first-filter reason** (duplicate, low MAPQ, supplementary, secondary, unmapped-but-placed, QC-fail) — exactly **two** conversions happen, the two clean reads are served, and the whole `ReadFilterCounts` matches field for field |
| `a_read_the_second_filter_drops_has_already_been_converted` | a 20-base read clears the first filter, is converted **once**, and is then charged to `too_short` |

**Suite: 2,867 → 2,869 (+2).** Fully accounted.

## 7. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,869 passed**, 0 failed, 5 ignored |
| `cargo test --examples` | 52 passed, 0 failed |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` against the `f5630f8` baseline | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

The dumps *had* to be byte-identical here in a way they did not at D1: this step changes production
code, and now adds a shipped counter. That they are is the evidence the change is behaviour-free.

## 8. Design documents amended

The owner lifted the no-edit bar (2026-08-04), so the two claims this milestone found to be wrong
are corrected at the source rather than left as notes in code:

- **spec §5's overclaim**, found at D1 — *"the reference stops being a precondition for filtering
  at all"* is true of the filter, the reader and the narrowing, and false of
  `AlignmentCursor<R: RawRefSeq>` and `AlignmentFile::cursor`.
- **spec §8's second test and arch §3.4's loop** — the plan's D2 fixture was impossible, and what
  replaced it is a counted conversion step.

They travel in their own commit, so this step's diff stays about the code.
