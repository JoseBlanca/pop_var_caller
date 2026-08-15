# B4 — the same cohort fitted from memory and from files, on real reads

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step B4 — implementation plan 2,
the last step of milestone B.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§7.15 (both halves).
**Date:** 2026-08-14.

---

## 1. The answer

**Every fitted number is identical, to the last bit, on both cohorts** — the eight tomato
accessions and the GIAB trio, each fitted once from a census held in memory and once from one
read off files written during the same run.

| | files written | read back to fit | agreement |
|---|---:|---:|---|
| tomato, 8 accessions | 0.545 MB, **0.068 MB a sample** | 0.465 MB | every number identical |
| the GIAB trio | 0.819 MB, **0.273 MB a sample** | 0.773 MB | every number identical |

**The fit read 85% of what tomato's files hold and 94% of the trio's, and the shortfall is the
repeat tracts** — the ordinary-position half never asked for a tract section, so those bytes were
never touched. That is the by-section design doing the one thing it exists for, measured rather
than asserted.

*(Measured with `CENSUS_FILES=<dir>` on `tmp/oracle_tomato_files.sh` and `tmp/oracle_trio_files.sh`,
which are the two oracles with that variable set. At a 60,000-position census, a tenth of the real
target; spec §6.1 prices a two-million-position census at about 6 MB a read group, and these are
0.068 and 0.273 MB at a thirtieth of the positions.)*

## 2. The counting reader, which is the half worth having

Spec §7.15 asks for two things, and the second is the one that separates a real
section-by-section reader from a decoder that reads the file and hands back a slice: **the bytes
actually read**. Both give the same values; only the count tells them apart.

`census_file::bytes_read()` counts, per thread, the bytes that leave a census file — incremented
where the read happens, so what is asserted is the read itself and not a promise about it.

| assertion | measured |
|---|---|
| a call for one stratum reads that stratum's extent and no other byte | exactly the extent, and it is smaller than the file |
| a band of two sections reads exactly those two | the sum of their extents |
| a second call for the same band reads them again | twice the sum — **nothing was retained between the calls** |
| a resident census reads nothing at all | 0 |

## 3. Tests

Three new, on top of B2's parity of values:

| test | what it pins |
|---|---|
| `a_cohort_fitted_from_files_gives_the_parameters_it_gives_from_memory` | three samples over 400 drawn positions, fitted both ways: the log-likelihood, the density's two shapes, the mismapped share, and every sample's heterozygosity, homozygote excess and error rate, all by exact equality |
| `asking_for_one_section_reads_that_section_and_no_other_byte` | §2's four rows |
| `a_resident_census_reads_nothing` | the other end of the same property |

And in the walk over real alignments, `CENSUS_FILES=<dir>` writes every sample's census, opens
them again as files, refits, and prints **the largest disagreement anywhere** rather than the
second answer — on both cohorts it is nothing.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo clippy --lib --all-features -- -D warnings` | clean; the example this step touched is clean too |
| `cargo test --lib ng::parameter_estimation::joint::census_file` | `15 passed; 0 failed` (12 before) |
| `cargo test --lib` | `3,607 passed; 0 failed; 11 ignored` (3,604 before) |
| the 88-second tomato oracle, no files | byte-identical to B2's |
| the 74-second trio oracle, no files | byte-identical to B2's |

**Checkpoint B is met**: a census on disk, read a section at a time, with the fit's answers
unchanged on both cohorts.
