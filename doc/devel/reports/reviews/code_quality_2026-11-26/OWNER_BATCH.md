# Owner-decision batch (2026-11-28)

The owner approved these decisions from `PLAN.md` §2B. Each package agent reads
`PACKAGE_BRIEF.md` first; this note **overrides it where they differ**.

## What the owner decided

| # | Decision |
|---|---|
| B2 | Fix the homozygote-excess M-step. Then rerun the ten-sample cohort named in `JUMP_SLACK`'s comment and report whether the slack can go. Do not remove it in this batch. |
| B3 | Fix the two STR-fit defects: reset the change measures after a jump; restore the fallback when the round limit ends on a jump. Leave out the `settled` flag (no format change). For the clamp, only **measure** how often the level passes 0.95, and report it. |
| B4 | Seed the census writer's STR map with every declared read group. |
| B5 | Reads that stop the walk today — deletions longer than `max_record_span`, and reads with no base qualities — are dropped before the walk and counted under a new read-filter reason (tally field, and psp provenance key). |
| B6 | Refuse `--max-locus-span`, `--max-str-len` and `--max-str-bundle-len` above 65,535; check first whether more code assumes `u16`. Document that `--threads 0` on `repeat-catalog` means one thread. |
| B7 | Include `--max-str-len` in the catalog reader's padding. |
| B8 | Test the STR read's "outside the tract" condition before the read cap. Rename `reads_discarded_by_cap` or add a second counter, so neither changes meaning silently. |
| B9 | Store each CRAM index entry's running maximum reach. Add a CRAM fixture with one long read. |
| B13 | Remove the `lowDepth` FILTER line from the VCF header. |
| B16 | Delete the per-locus slippage re-fit. |
| B17 | Delete contamination's five research options; keep their measurements in the module doc. |
| B18 | Delete `CensusMomentSums`. First check that no `Debug` dump of `JointFit` reaches a user file. |
| B19 | Remove the reading of parameters files written before 2026-10-16: `share_curve.rs`, most of `slippage_curve.rs`, the two provenance variants. Do it on both sides, writer and reader. |
| B20 | Delete `PspWriter::append` and its tests. Record "not implemented" in the psp spec §6.4. |
| B22 | Add an `fsync` to the catalog writer before its rename. |
| B23 | `--max-memory` sizes are **binary**, the way RAM is sold: `4G` = 4 × 2³⁰ bytes. Make the parsing and every log line that prints a memory size agree with that, including `peak.rs`. Use the same spelling for every log line that prints a size. |

Not in this batch: B10, B11, B12, B14, B15.

## Rules for this batch

- **Deletions and B22 change no output.** `scripts/code_quality_oracle.sh --identity` must still pass on those commits.
- **The fixes change output by design** (B2, B3, B4, B5, B7, B8, B9, B13). For each fix commit:
  1. Add a test that fails before the fix and passes after, as usual.
  2. Run `--identity`. It will report `DIFFERS`.
  3. Check that only the outputs your fix should move have moved, and say why each one moved. For a VCF, give the number of records whose genotype, QUAL or GQ changed. For a parameters file, give the fitted numbers before and after.
  4. Put that evidence in the commit message and in your reply.
  5. Do not re-record the baseline: the manager does that once, after the whole batch is merged.
  6. An output that changes for a reason your fix does not explain is a stop: report it.
- **Keep each fix in its own commit**, separate from any refactor or deletion. Do deletions before fixes in the same file.
- **Tests run without `--release`** (see `PACKAGE_BRIEF.md`, How to build and test).
- **Text the gate cannot see.** Error text, log text and `--help`: show before and after in the commit message.
</content>
</invoke>
