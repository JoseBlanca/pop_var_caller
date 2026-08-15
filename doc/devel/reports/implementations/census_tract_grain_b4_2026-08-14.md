# B4 — a repeat tract's two counts move to the stratum, and its walked flag to a bit

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md) step B4.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§3; [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.4.
**Date:** 2026-08-14.

---

## 1. What changed

**Two of a tract record's counts are now kept once per (read group × stratum) instead of once
per tract.** A stratum is a *(motif period, reference repeat count)* pair — the group the
slippage numbers are fitted within.

- **reads that reached a tract and crossed none of it** — the censored ones, which report no
  length. On the trio 37 reads in every 100 that reach a tract never cross it;
- **the bases the substitution rate is fitted against** — the denominator of a rate that is
  itself fitted per stratum.

**And the flag saying the walk reached a tract is a bit rather than a byte** (`WalkedBits`, the
STR half's `PackedDepthCodes`).

**The estimator's accumulator moved into the writer.** `gather_strata` used to add the
per-locus censored count into a per-stratum running total the moment it read it, once per
(locus × sample × read group). It now reads one number per (sample × read group × stratum).
The total is unchanged, and §2 measures that it is.

## 2. What moved on real reads

Both oracles, before and after, at four container CPUs.

**No fitted number moved, on either cohort, and the censored-read totals are identical**:
**102,308** reads reached a tract without crossing it on the trio and **19,008** on the eight
tomato accessions, before and after. That is the assertion this step needs — the count is a
denominator, and a wrong one is a plausible rate rather than a crash.

The whole diff of both runs is the size of the repeat-tract record:

| | tracts | strata | tracts a stratum | before | after |
|---|---:|---:|---:|---:|---:|
| tomato, 8 accessions | 349 | 26 | 13.4 | 25.0 bytes a tract | **21.0 to 21.3** |
| the GIAB trio | 216 | 32 | 6.75 | 25.0 bytes a tract | **25.2** |

**On the trio the record grew, and the reason is arithmetic rather than a defect.** Per tract
the change saves 7 bytes (a `u16`, a `u32` and a `bool`) and spends an eighth of one; per
stratum it spends a pair of map entries, 48 bytes. So the record shrinks wherever a stratum
holds more than about **seven** tracts. Tomato's eight accessions sit at 13.4 tracts a stratum
and save a sixth; the trio sits at 6.75 and pays 0.2 bytes a tract. It is a property of a
216-tract census, not of the encoding.

*(The spread across tomato's eight accessions — 20.9 to 21.3 — is the number of strata each
sample actually put a censored read in: 20 to 23 of the 26 the selection holds. A stratum a
sample never charged costs it nothing. The trio's 25.2 is all 32 charged.)*

**At the scale the census actually runs at, the saving is the one the specification prices.**
Tomato SL4.00 holds **462,701 kept tracts in 141 strata**, 3,281 a stratum: 18 bytes of offset
buckets, an eighth of a byte of walked flag and 0.015 bytes of per-stratum pair, so **25.0 →
18.1 bytes a tract**, and **11.57 MB → 8.39 MB a read group**. Spec §3 predicts "about 3.2 MB
of 11.57 MB a read group — a quarter of the repeat-tract record"; this is 3.18 MB and 27%.
**That figure is arithmetic from the two measured inputs** (the per-tract bytes above and the
locus and stratum counts from `ng_joint_loci_probe`), not a measurement — the full-cohort run
in B5 is where it is read off a real census.

*What the byte figures count.* The harness prices a vector by its element size and a map by
its entries' content, so the "bytes a tract" above are content bytes and exclude a `BTreeMap`
node's own overhead — the same convention every other row of that table already used.

## 3. What this step gives up, and it is in the specification

**A tract can no longer say, on its own, that reads reached it and stopped inside it.**
`SsrLocusState` had four states and now has three: never walked, no length reported, crossed.
The fourth — *reached, not crossed* — was the per-locus count that moved, and per locus the fit
never acted on it: it summed it on sight.

Spec §7.4 anticipated exactly this. It asks for the state to be planted deliberately, "so the
test either shows a lower bound surviving or shows that it does not". It does not survive per
locus; it survives per stratum, which is the grain the loss actually runs along — a tract
longer than the reads is never crossed, in any sample at any depth, so the censoring follows
repeat count and repeat count is half of what a stratum is.

**`Stratum` moved from `ssr_fit.rs` into `census.rs`** and is re-exported where it was. The
census is now keyed by it, and the data module must not depend on the module that does the
mathematics — the arch document already places it in the census (`SectionKey::Ssr(ReadGroupId,
Stratum)`, §1.1a).

**One widening, recorded because it changes a stored number in a case nothing reaches.** Reads
that reached a locus and produced no observation at all used to be added into a `u16` that
saturated at 65,535; the per-stratum count is a `u64` and does not. No locus in either oracle
comes near it.

## 4. Tests

| test | what it pins |
|---|---|
| `one_locus_marked_walked_leaves_every_other_locus_alone` | the bit array: twenty loci in three bytes, set at every offset within a byte, no neighbour touched |
| `the_three_states_at_an_str_locus_are_distinguishable` | what a locus can still say — and the doc comment says which state left and where it went |
| `a_read_that_reached_a_tract_without_crossing_it_is_counted_against_its_stratum` | the censored reads land in their own stratum and no other, driven through the writer |
| `two_strata_are_charged_separately_for_their_censored_reads_and_compared_bases` | **the one thing a per-locus count could not get wrong**: two tracts of twelve bases differing only in motif period, each stratum charged its own censored reads and its own compared bases |
| `a_slipped_read_contributes_a_length_and_no_base_comparison` | unchanged in what it asserts; the denominator is now read by stratum |

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo test --lib ng::parameter_estimation::joint::census` | `32 passed; 0 failed` (29 before) |
| `cargo test --lib` | `3,581 passed; 0 failed; 11 ignored` (3,578 before) |
| the 88-second tomato oracle | §2 — no fitted number differs |
| the 74-second trio oracle | §2 — no fitted number differs |

**The two red gates are the two that were red before this branch's first commit**, neither in
code this plan touches: `cargo clippy --all-targets --all-features -- -D warnings`, and
`cargo test --all-targets`, which panics in `benches/psp_writer_perf.rs:386`.
