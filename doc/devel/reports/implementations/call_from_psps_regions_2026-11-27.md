# `call-from-psps --regions` — implementation report, 2026-11-27

Plan: [`call_from_psps_regions.md`](../../ng/impl_plan/call_from_psps_regions.md), steps 1 to 7.
The owner's decisions:
- a locus belongs to the run whose regions hold its first base;
- the psps are opened and the BED checked before the reference is read;
- the hidden-duplication filter is fitted per job, and this is documented.

## What a person can now do

Write one psp per sample over the whole genome, run `estimate-parameters` once over all of them,
then call each chromosome as its own job with
`call-from-psps --regions chrN.bed`, and join the VCFs with `bcftools concat`.

If the BED asks for any base the psps do not cover, the run stops before the reference is read.
Example:

```
error: --regions asked.bed asks for regions the psps do not cover, so nothing was called:
  chr1:101-150   on chr1 the psps cover chr1:1-100
The psps cover chr1:1-100 (100 bases). Positions are 1-based, both ends included.
Remove these regions from the BED, or write the psps again over regions that include them.
```

`pvc` exits with status 1 and writes no VCF. A contig the reference lacks, or an empty BED, is
refused with the same messages `--regions` gives on every other command.

## Changes, by step

1. **The set difference.** `GenomeRegions::parts_outside` returns the requested regions minus
   the regions the psps cover, in one pass over the two sorted lists. `GenomeRegions::on_contig`
   was added for the message. Both are in `src/region_typing/mod.rs`.
2. **The flag and the refusal.**
   - `--regions` on `CallFromPspsArgs`.
   - `run_ground::regions_inside_the_psps` and `GroundError::RegionsOutsideThePsps`. The refusal
     lists at most 10 pieces, then counts the rest.
   - `OpenPspCohort::contig_table`.
   - The cohort is now opened before the reference is read (`src/cli/call_from_psps.rs`).
3. **The caller accepts a part of what the files cover.**
   - `SegmentationInputs::first_difference_in_how_repeats_were_sorted` checks the setting, the
     catalogue and the criteria, but not the regions.
   - `PspVariantCaller::open` uses it, then refuses a run whose regions reach outside a file's
     with the new `RunError::RegionsOutsideThePsp`. That check is defensive: the command refuses
     first.
   - The segmentation is cut over the requested regions.
4. **The sources skip what the run does not call.**
   - New `ObservationSource::nothing_wanted_before`, a no-op by default. The cache calls it at
     every eviction, from `SampleWindow::evict_before`, so all three merge drivers get it.
   - New `PspSummarySource::over_part`. It enters the file at the first region's start minus the
     observation reach ceiling (`records_from`), then drops, unread, every record that ends before
     the floor. A record is dropped the way the depth ceiling drops one: its bytes are truncated
     and `HeldLiveSets::skip` is applied.
   - **The source is used only when the run's regions differ from the files'.** A run without
     `--regions`, or with a `--regions` naming everything the psps cover, builds the same sources
     as before. Its output is unchanged by construction, not by argument.
5. **Proof** (`src/cli/mode_equivalence.rs`), on a new two-contig varying cohort
   (`test_fixtures::a_two_contig_varying_cohort_on_disk`). Each contig has a SNP in each sample
   and an STR length variant.
6. **Report and log.** For a `--regions` run, the report's ground line names the BED and what the
   psps cover. The debug log names both too. Code: `RunTail::regions_asked_for` and
   `RunReport::with_regions_asked_for`.
7. **Docs.**
   - Developer docs: the module comment, `PspVariantCaller::open`, the help text, spec
     `run_streaming.md` §5.3 and §11 question 5, and `command_line.md` §4.2 and §5.2.
   - User docs: `doc/user/calling.md` (the psp paragraph and "Splitting a run by region"),
     `options.md`, `limits.md`, and the cluster recipe in `recipes.md`.

## Step 5 results

All on the fixture: two samples, about 20 reads a position, contigs of 600 bases. This shows
the rule holds on records that carry calls. It says nothing about rates on real data.

| Oracle | Result |
|---|---|
| (a) `--regions chrW` on whole-genome psps against psps walked over chrW only, filter on (0.01), refinement on its default | same VCF, every line; 3 records, all on chrW |
| (b) chrV job then chrW job, joined, against one whole run, filter off | same records, 6 of 6 |
| (c) cuts inside chrV after bases 60, 121 (on sample one's SNP) and 210 (inside the STR) | the two jobs' records, joined, are the whole run's, every record compared whole, including chrV:121 on the cut |
| (d) `--regions` naming everything the psps cover, against no flag, filter on | same VCF |
| A whole-genome run's parameters file used for a one-contig job | binds and calls |

**The design allows differences near a cut, and none can occur on this fixture.** It allows
three:
- a phase set split in two;
- an STR record widened over its flank that is not merged with a generic record across the cut;
- window depth measured without the records before the cut, which only the hidden-duplication
  filter reads.

None can occur here. No record is phased (`PS` is `.`), no generic record lies within a flank of
the tract, and the filter is off. So (c) compares every record whole. The record on the cut,
chrV:121, is written by the first job, whose regions hold its first base. It is identical in the
split and whole runs at all three cuts:

```
chrV 121 . T A 232.0 PASS AF=0.250000;AC=1;AN=4;DP=40;ABPEN=0.0;SPPEN=42.4;BQPEN=0.0;MQREF=60.00;MQALT=60.00;MQDIFF=0.00 GT:GQ:DP:AD:PL:GP:PS 0/1:99:20:10,10:307,0,307:0,1,0:. 0/0:62:20:20,0:0,60,735:1,0,0:.
```

The allowed differences are still untested by any fixture.

**Mutation check.** With the floor disabled (`nothing_wanted_before: None` in `over_part`), (a)
fails. The first contig's records reach the hidden-duplication filter's depth model, and
`PARALOG_LR` moves, for example from −5.53 to −5.65 at chrW:121.

## Blocks read

On `psp_source`'s multi-block fixture (two contigs, four blocks each), a source over the second
contig began 5 of the file's 8 blocks. Those are the second contig's four, plus the first contig's
last block, which the index cannot rule out. Test:
`a_source_over_the_second_contig_reads_only_the_blocks_that_can_hold_it`.

## Tests added

- `region_typing`: 6 example tests and a property test against the set difference taken base by
  base.
- `run_ground`: 5 tests (inside and merged, outside with the full message, a gap, more than 10
  pieces, BED-reader refusals).
- `call_from_psps`: the flag parses; refused before the reference is read, with the exact
  message and no VCF; a missing contig; equal to everything the psps cover; part of it is
  called.
- `psp_source`: the seek, and the floor.
- `psp_caller`: a run over part of the files opens, and one reaching outside is refused.
- `report`: the ground line.
- `mode_equivalence`: the five oracles above.
- `tests/call_from_psps_regions_refusal.rs`: the `pvc` binary, which checks exit status 1, the
  printed message and that no VCF is written.

## Validation

- `cargo test --lib`: 5,323 passed, 0 failed.
- The integration test passes.
- `cargo clippy --all-targets`: no warning in a changed file; the remaining warnings are older.
- rustfmt was run on the changed files only.

## Left open

- **Gaps between regions are read through, not skipped.** A BED of many intervals spread over
  the genome reads each file from its first interval onwards. It never reads more than a run
  without `--regions` would. A second seek is a follow-up, if a measured BED needs one.
- **Not measured on real data**: wall time and memory of per-chromosome jobs, and how many calls
  the per-job hidden-duplication fit changes.
- `PROJECT_STATUS.md` was not updated.
