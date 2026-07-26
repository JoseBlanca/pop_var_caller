# ng read preparation — C1–C3: parity, the soft-mask measurement, and the `F1` check

**Date:** 2026-07-26 · **Plan:** [read_preparation.md](../../ng/impl_plan/read_preparation.md) steps
C1, C2, C3 · **Spec:** [read_preparation.md](../../ng/spec/read_preparation.md) §5, §6, §9

Combined implementation + review + fixes report. Milestone C turns three claims the spec review left
as arguments into measurements.

## 1. Plan

One `#[cfg(test)]` parity module, `src/ng/read/left_align_parity.rs`, following the pattern of
[delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs) — ng's test, production as the
yardstick, no test-only dependency in shipping ng code.

## 2. Changes made

**[src/ng/read/left_align_parity.rs](../../../../src/ng/read/left_align_parity.rs)** (new): a FASTA
writer that preserves case (the soft-masked arm depends on it), a `Repository` builder, both sides'
`prepare_read`, a field-by-field `PreparedRead` comparison, and one fixture of eight reads.

The reference is `chr1 = CCTAAAAGGTTTTTACGCATGCGCG` — an `A` run at 4–7, a `T` run at 10–14, unique
sequence from 16. The reads cover: no indel; a deletion in the `A` run spelled at its right end; an
insertion after the `T` run; the same deletion already leftmost; a deletion in unique sequence; a
soft-clipped read with a shiftable deletion; a reverse-strand read; a paired second-of-pair.

## 3. Tests added (5)

| test | what it establishes |
|---|---|
| `ng_matches_production_on_an_uppercase_reference` | **the port anchor** — every field of every read byte-identical to production's `--no-baq` output |
| `the_parity_fixture_actually_exercises_left_alignment` | the fixture is not vacuous: exactly `del_in_a_run`, `ins_in_t_run`, `reverse`, `soft_clipped` are rewritten, **by name** |
| `the_shiftable_reads_land_on_their_leftmost_spelling` | the shifted spellings themselves, so an indel moved to the *wrong* place fails here, not only against production |
| `production_left_alignment_does_nothing_on_a_soft_masked_reference` | **the measurement** — see §5 |
| `left_alignment_cannot_change_the_mismatch_fraction_verdict` | discharges spec §5's step-1/step-2 ordering argument, across six thresholds |

## 4. The vacuity check earned its keep immediately

The first draft of the fixture **passed byte-parity while proving almost nothing**: only one read was
actually rewritten. Three reads were malformed, and each failed silently:

1. `del_in_a_run` was named for the `A` homopolymer but its `M5 D1 M4` deleted the **`G` at 8**, whose
   left neighbour is `A` — nothing to shift.
2. `unique_context`'s `M3` covered reference `GCA` while the read supplied `CAT` — every base
   mismatched.
3. `soft_clipped`'s CIGAR consumed nine read bases while its `seq` had ten. That is a **malformed
   read**, and `left_align_cigar` bails safe on exactly that, handing the CIGAR back untouched — so
   both sides "agreed" while normalizing nothing.

All three passed the parity assertion. The lesson is in the fix: the check now asserts the rewritten
reads **by name**, because a floor ("at least N moved") would still have admitted two of the three.

A second error was caught the same way: the expectation said three reads shift when the answer is
four — `reverse` carries the same shiftable spelling as `del_in_a_run`, being the same deletion on the
other strand.

## 5. The soft-mask measurement — mechanism confirmed, exposure small

`production_left_alignment_does_nothing_on_a_soft_masked_reference` asserts three things; the middle
one is the finding, and it is stated positively rather than as a disagreement:

1. ng's answer on masked sequence is **identical to its own answer on unmasked sequence**;
2. production's answer is **the mapper's original CIGAR, byte for byte** — it normalized nothing;
3. so the two disagree, on exactly the reads with somewhere to shift.

**Then the exposure was measured on the real references** (full scans, not samples):

| reference | lowercase bases |
|---|---|
| GRCh38 no-alt analysis set | **0** |
| tomato SL4.0 | **227,170** — ch01 507, ch03 13,008, ch05 2,416, ch06 30,813, ch07 19,733, ch08 41,118, ch09 33,895, ch10 46,229, ch11 27,961, ch12 11,490 |

So the defect is **real but latent**: ~0.03% of the tomato genome, none of GRCh38. It **cannot**
explain production's indel deficit, which was measured on GRCh38. It would bite a run against UCSC
`hg38.fa` or a RepeatMasker-masked plant genome.

**A related correction, recorded in `PROJECT_STATUS.md`:** the 2026-07-24 normalizer screen does not
bear on this. It fetched the **canonical** view ([ng_normalizer_screen.rs:228](../../../../examples/ng_normalizer_screen.rs))
and compared ng's three normalizers **to each other** — never against production's raw-byte
behaviour. "Which normalizer" and "does production's normalization work at all" are different
questions, and only the first was tested. The screen's conclusion had been compressed into
"normalization placement is not the lever behind the indel deficit", which does not follow.

## 6. Validation

In the container: `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
and `cargo test --lib --all-features` — results recorded in the commit message. `cargo test
--all-targets` additionally runs bench mains and hits a **pre-existing** panic in
`benches/psp_writer_perf.rs:386` (`index out of bounds: the len is 3300000 but the index is 3300000`),
proven unrelated by re-running that bench with this work stashed.

## 7. Review and fixes applied

Reviewed inline. No Blocker, no Major. The substantive findings were the fixture defects in §4, all
fixed in the same step; the rest were import tidying (fully-qualified paths pulled up into `use`).

## 8. Follow-ups

- **The magnitude on real data is still unmeasured.** This is a synthetic 25-base fixture: it proves
  the mechanism, not the size. The honest next step, if anyone wants it, is the normalizer screen run
  raw-vs-canonical over tomato's 227 kb of masked sequence.
- The two `OPEN:` items from arch §6 are untouched: `Ok(None)` carries no decline reason, and
  call-vs-port of `prepare_passthrough` stands at "call it".
