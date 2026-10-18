# `--call-satellites` and `--generic-path-only` — implementation plan

**Status:** draft for the owner's review, 2026-10-18. Step 6 of
[`command_line.md`](command_line.md), which says this step gets its own short plan before any
code. The design is spec [`command_line.md`](../spec/command_line.md) §4.4.1 (owner's rulings,
2026-10-16); this plan answers the three questions that section leaves to implementation, says
where each change goes, and stops where a choice is the owner's.

**Branch:** `call-satellites`, worktree `../pop_var_caller-call-satellites`, from `main` at
d963b776.

**Why a file of its own** rather than a section of `command_line.md`: the other steps of that plan
are interface work that must not change a call, and this one changes calls, adds a psp field and
has its own accuracy check. It is as long as the rest of that plan together. `command_line.md`'s
Step 6 points here.

---

## What the two options do, in one paragraph

A run sorts the reference into four kinds of ground (spec §4.4.1's table): ordinary sequence,
which the SNP/indel model calls (the *generic path*); STR tracts and STR bundles, which the STR
model calls; and *satellites* — a tract longer than `--max-str-len` or a bundle longer than
`--max-str-bundle-len`, both 100 bp — which a default run does not call. `--call-satellites`
hands satellites to the generic path. `--generic-path-only` hands all repeat ground to it and
needs no catalogue. A run given both is refused.

## The three questions the spec left open, answered from the code

### 1. Can the ground be sorted with no catalogue? Not today; it needs a new branch, not an empty file

Every route that sorts the ground refuses a missing catalogue before it starts:
`refuse_a_catalog_that_is_not_there` (`src/cli/run_ground.rs:294`, again at `:375`) and the
reader itself (`src/repeat_catalog/reader.rs:63-67`). The catalogue is also part of what a psp
records: `SegmentationInputs.catalog` is a required `RepeatCatalogHeader`
(`src/segmentation_inputs.rs:31`), written to every psp as `[segmentation.catalog]`
(`src/psp/segmentation_section.rs:121-140`) and compared across a cohort
(`segmentation_inputs.rs:67-79`).

A catalogue with zero rows *would* give one generic stretch per contig, because the walk loops
over the header's contig table, not over the rows (`reader.rs:238-251`;
`parquet_file.rs:593` returns no rows for a contig with none, and `fill_generic_gaps`,
`src/region_typing/mod.rs:1405-1432`, fills the whole contig). But an empty catalogue is still a
file the user must build or be handed, its header must match the reference
(`reader.rs:593-640`), and the psp would record a catalogue that played no part. The spec says
no catalogue is needed, so:

**The plan: a new branch in `run_ground`, beside `segments_cut_from` (`run_ground.rs:399-429`),
that builds the segmentation as one `Generic` stretch per analysed span and opens no file.**
`Segmentation::build` (`src/run/segments.rs:71-97`) already takes any stream of typed regions; what
changes is what it records (next section). Two other users of the catalogue need the same branch:

- the census plan `generate-psps` builds, which opens the catalogue to choose the STR loci it
  samples (`src/cli/generate_psps.rs:652-669` → `select_kept_loci`,
  `src/parameter_estimation/joint/loci.rs:827-856`). Under `--generic-path-only` the generic
  sample is drawn from the analysed regions directly and the STR sample is empty;
- `estimate-parameters`, which compares a given `--catalog` with the one the psps were walked
  with (`src/cli/estimate_parameters.rs:843-863`). Under `--generic-path-only` psps it refuses
  `--catalog`, as `generate-psps` does.

### 2. How satellite ground joins its neighbours: always ordinary sequence, so it merges into one stretch

**A satellite never touches an STR tract or bundle.** When the ground is sorted, anything within
`bundle_threshold` (15 bp) of a satellite is absorbed into it and the satellite grows to cover it,
repeated until nothing more is in reach (`src/region_typing/mod.rs:908-941`). A bundle that becomes
a satellite already holds everything within that distance of its members (`mod.rs:806-810`). So
each satellite has ordinary sequence on both sides — at least 15 bp of it — or a contig end. The
junction the spec asked about is therefore between ordinary sequence and ordinary sequence, and
the STR edges further out are the existing tract-to-generic junctions, unchanged.

**The handed-over satellite must merge with the ordinary sequence either side into one generic
stretch, not sit beside it as a second one.** A generic stretch is where the SNP/indel model builds
its loci, and an indel cannot span two of them: a deletion running past a stretch's last base is cut
there (`src/locus_generation/pileup/open_record.rs:2994-2998`), and an insertion on its last base is
left for the next stretch (`open_record.rs:2595-2634`), which only claims it when that next
stretch follows repeat ground (`src/locus_generation/mod.rs:1099-1107`). The code states the rule:
"Split a run at *p* and an indel spanning *p* is callable by neither half"
(`src/region_typing/mod.rs:1398-1402`). One merged stretch has no edge there to lose anything at.

**Where the merge goes: when the run assembles its stretches, not inside the catalogue reader.**
The reader's stretches are also what the census samples the SNP/indel fit's positions from
(`loci.rs:836-846` keeps only `Generic` stretches). Dropping satellites inside the reader would
add satellite positions to that sample, and the spec says `--call-satellites` changes nothing in
the fit. So the reader stays as it is, and a new step between the reader and `Segmentation::build`
in `segments_cut_from` (`run_ground.rs:411-428`) does three things:

1. **clips each satellite to the analysed spans.** A satellite is emitted whole even past a BED
   edge, as a tract is (`src/repeat_catalog/segments.rs:101-119`; `clips_at_a_bed_edge`,
   `mod.rs:688-693`); as generic ground it is clipped, as all generic ground is;
2. **merges it with the generic stretches either side** into one;
3. **keeps the list of satellite spans it handed over**, for the census and the run report below.

What follows from it, by design:

- **STR tracts a satellite absorbed stay inside it** and are called by the generic path too. The
  absorption happens when the ground is sorted (`mod.rs:941`), and undoing it would make a
  different sorting, not this option.
- **The census.** `generate-psps` marks every generic stretch as walked
  (`src/run/gatherer.rs:315-318`). The positions it records are the ones `select_kept_loci` chose,
  which do not include satellite ground, so marking it walked should change nothing. The test
  below pins that: one fit from default psps and one from `--call-satellites` psps of the same
  data must give byte-identical parameters files. If they do not, the census is handed the
  satellite list and leaves those spans out.
- **The run report.** Today a satellite is counted as "not called — tandem arrays longer than this
  run types as callable" (`src/run/report.rs:284-288`, set at `locus_generation/mod.rs:1109-1112`).
  Under `--call-satellites` that line reads zero and a new one counts the satellite bases handed
  to the generic path. The width refusal line that the spec asks to keep is already there
  (`report.rs:355-364`, "loci the merge declined to assemble for being too wide").

### 3. Does the parameters file accept a missing stutter section? No — but it accepts an empty one

**There is no stutter section that can be left out.** Every STR number sits in `[repeat_tracts]`,
which is a required field of `ParametersFile` (`src/calling/parameters_file/mod.rs:584`, not an
`Option`, no serde default), the writer always emits it (`to_toml.rs:426`), and a file without it
is refused at parse time as malformed (`from_toml.rs:88-97`). Every table inside it is required
too (`to_toml.rs:2226-2227`).

**An empty section is accepted everywhere, and is what a fit over zero repeats writes today.**

- The stutter fit with no repeat classes logs "there was nothing to fit" and returns nothing,
  without an error (`src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:846-859, 900-903`); tests
  already write and read such a file (`src/run/census_fit.rs:1357+`, read back at `:1524`).
- Validation passes it (`validate.rs:428-435`: "A run with no repeat tracts at all declares
  none"), and the file says in words that nothing was fitted (`to_toml.rs:556-563, 782-832`).
- At calling time each tract then uses HipSTR's stutter model and the 0.001 substitution rate,
  counted as defaulted (`src/calling/inference/repeat_tract_parameters.rs:378-412, 549-583`).
- Under `--generic-path-only` there are no tracts, so no calling run reads any of it.

Leaving the section out, as the spec's ruling words it, means format 6 of the parameters file:
`repeat_tracts` becomes an `Option`, version-gated (as `validate.rs:822-836` gates format 5's key),
and about ten readers learn to handle `None` (`to_run_parameters.rs:369, 510`, `bindings.rs:711,
733`, `what_was_fitted.rs:136-163`, `validate.rs:436-529`, `from_run_parameters.rs:419, 688`). It
also moves all four parameters-file and VCF digests in `src/cli/cross_platform_digests.rs` (`:267`,
`:395`, `:502`, `:608`), since every file's version line changes. **Decision A** below.

---

## The psp header

**A new required key in `[segmentation]`, `repeat-ground`, with one of three values:
`str-paths` (default), `call-satellites`, `generic-path-only`.** On the typed side it replaces
the two repeat fields of `SegmentationInputs` with one enum, so a value that cannot happen — a
catalogue recorded for a run that used none — cannot be built:

```rust
pub enum RepeatGround {
    /// Tracts and bundles to the STR paths; satellites not called.
    StrPaths { catalog: RepeatCatalogHeader, repeat_tract_criteria: StrRepeatCriteria },
    /// As `StrPaths`, and satellites to the generic path.
    CallSatellites { catalog: RepeatCatalogHeader, repeat_tract_criteria: StrRepeatCriteria },
    /// Every repeat to the generic path; no catalogue was read.
    GenericPathOnly,
}
```

Under `GenericPathOnly` the header has no `[segmentation.catalog]` and no
`[segmentation.repeat-tract-criteria]` table; under the other two it has both, as today.

**Mixing settings in a cohort is refused by the check that already refuses mixed criteria.**
`first_difference` (`segmentation_inputs.rs:67-79`) gains the setting as the first field it
compares, so the existing `RunError::CohortWalkedUnderDifferentSettings`
(`src/run/mod.rs:611-622`, raised from `psp_caller.rs:1039-1062`) names it: "samples A and B do
not agree on the setting for repeat ground". `call-from-psps`, `estimate-parameters` and
`regenerate-census` all open psps through that check.

**The format version: 2.1 → 3.0, and a 2.x psp still reads, as `str-paths`.**

- *Why a major bump.* A 2.x reader ignores keys it does not know (`header.rs:2206` pins that). A
  2.1 build handed a `call-satellites` psp would read it as a default one, sort the ground with
  satellites, and meet generic records inside what it thinks is a satellite. A major bump makes
  it refuse the file by version instead (`header.rs:725-733`), as 2.0 was bumped for the same
  reason (`header.rs:483-490`).
- *Why 2.x still reads.* Every 2.x psp was written by a build that had only the default, so
  reading it as `str-paths` is a fact about that version, not a guess about a missing key.
  `OLDEST_READABLE_MAJOR` (`header.rs:511`) stays 2. In a 3.x psp the key is required, and a
  missing one is refused by name, as `max-str-bundle-len-bp` is (`segmentation_section.rs:394-400`).
  Nobody has to write their psps again. **Decision B** below.

The oracle does not compare psp bytes (the header holds a timestamp), so the default runs stay
byte-identical in everything it does compare.

## What does not change

- **The parameters file**, under Decision A's recommendation: no format change, so a default run's
  file is byte-identical. Its `[repeat_routing]` section records the six repeat settings and only
  reports a disagreement, never refuses (`bindings.rs:261`); a file fitted under
  `--generic-path-only` and used by a default run makes every tract fall back to HipSTR's
  stutter, which the run report already counts. That is the same as a fit that found no tracts
  today, so it is not refused.
- **The VCF header.** `##commandline` already records the option; no new line.
- **`--max-locus-span`** stays in force under both options, as the spec says.
- **`type-regions`**, the developer tool that prints the sorted ground, keeps printing satellites;
  it is what the accuracy check uses to find satellite ground.

## A fix this needs first: reference memory along one long stretch

The SNP/indel walk releases the reference behind it only when a new stretch begins
(`release_reference_behind`, `src/locus_generation/pileup/generator.rs:1223, 1286-1301`). Its
three readers — the walk's, the read preparer's and the read filter's — each keep everything from
the stretch's start. The code's own note sizes it: "On human chromosome 1 that is around 250 MB,
against a walk that otherwise peaks near 25 MB" (`generator.rs:1262-1263`). Today that never
happens because repeats cut human sequence into stretches of a few kilobases. Under
`--generic-path-only` a whole contig is one stretch, so a run would hold about 750 MB of reference
**per sample** — 750 GB at a thousand samples. `--call-satellites` creates stretches of megabases
around human centromeric satellites (HSat2/3, period 5).

**The fix: also release as the walk advances inside a stretch**, keeping the same margin
(`max_record_span`) behind the walk's position. Releasing is a hint; a base asked for after it was
released is read again (`generator.rs:1282-1284`), so no call can change. `resident_reference_bases`
(`generator.rs:1310`) gives the test its measure.

## Build order

Each step is built, its tests run, and committed; the owner looks before the next starts.

1. **Release the reference inside a long stretch.** No call changes; the identity oracle is run
   once here to show it.
2. **The setting as a type, in the psp header.** The `RepeatGround` enum in
   `SegmentationInputs`; the psp header's key, format 3.0, 2.x read as `str-paths`; the cohort
   check naming the setting. Every run still writes `str-paths`, so no call changes.
3. **`--call-satellites`.** The flag on `generate-psps` and `call-from-alignments` (*Advanced*),
   refusing it together with `--generic-path-only`; the clip-and-merge step in
   `segments_cut_from`; `call-from-psps` and `estimate-parameters` following the psps; the report
   line.
4. **`--generic-path-only`.** The flag; the refusals of `--catalog` and the six repeat settings
   alongside it; the segmentation with no catalogue; the census with no STR sample;
   `estimate-parameters` writing the empty stutter section (or, under Decision A's other
   option, none).
5. **Help text and spec.** Both options' help carries the spec's warning: on satellites they add
   mostly SNPs and short indels, and the satellite's own length variation mostly stays uncalled.
   The spec's §4.4.1 gets the answers settled here.
6. **Checks**, below.

## Tests added

| Step | Test | What it pins |
|---|---|---|
| 1 | a walk over one 1 Mb generic stretch holds no more reference than `max_record_span` plus a read's window per reader | the memory bound, measured with `resident_reference_bases` |
| 2 | each of the three settings round-trips through the psp header | the encoding |
| 2 | a 2.1 header reads as `str-paths` | old psps keep working |
| 2 | a 3.0 header without `repeat-ground` is refused naming the key | no silent default |
| 2 | a cohort of one `str-paths` and one `call-satellites` psp is refused naming the setting | the spec's mixing rule |
| 3 | both options together are refused before anything is read, on both subcommands | the spec's rule |
| 3 | a satellite between two generic stretches comes back as one generic stretch, and no two generic stretches are adjacent anywhere | the merge |
| 3 | a satellite crossing a BED edge is clipped to the span | the clip |
| 3 | a fixture with a SNP inside a satellite: not called by default, called under `--call-satellites`; a deletion across the satellite's old edge is called as one allele | end to end |
| 3 | default and `--call-satellites` psps of the same data fit byte-identical parameters files | the spec's "changes nothing in the fit" |
| 3 | direct mode and psp mode give identical VCFs under `--call-satellites` | the existing mode-equivalence test, run under the new setting |
| 4 | `--generic-path-only` runs with no catalogue file on disk | the spec's "no catalogue needed" |
| 4 | `--generic-path-only` with `--catalog` or any of the six repeat settings is refused | the spec's refusal |
| 4 | `estimate-parameters` on `--generic-path-only` psps writes a file the calling run reads, with an empty stutter section | the parameters-file answer |
| 4 | direct mode and psp mode give identical VCFs under `--generic-path-only` | mode equivalence |

## The checks before merging

**The identity oracle** (`scripts/promote_ng_oracle.sh`, four tomato2 accessions at low depth
over 20 regions, fitted run included): every checksum it compares must match the baseline. It
compares VCFs, window coverage and parameters files, not psp bytes, so the psp header change does
not reach it. Under Decision A's other option the parameters file's version line moves and the
baseline is re-recorded with a note.

**The accuracy check: `--call-satellites` against GATK on GIAB HG002, over satellite ground only.**

- *Satellite ground* is what `type-regions` types `Satellite` under default settings, intersected
  with the truth set's confident regions. The size of that intersection is reported first: the
  GIAB v4.2.1 confident regions leave out much repetitive sequence, so it may hold few truth
  variants, and the count says how much the rest of the numbers can carry.
- *For each caller* — pvc default (expected: no calls there), pvc `--call-satellites`, GATK — the
  number of calls in satellite ground, and how many match the truth, split into SNPs and indels;
  and the truth variants there that each caller misses. For pvc, also how many loci were refused
  for width there.
- *Data and depth*, each number labelled with them (**Decision C** below):
  - HG002 at about 300 reads a position over 1,000 regions (5.08 Mb,
    `benchmarks/human_genome_bottle/`), against GIAB v4.2.1. A GATK 4.6 VCF of these reads
    already exists, so this costs one pvc run.
  - HG002 at 30× and at 10× over the GIAB Tandem-Repeat benchmark's 50,000 regions
    (`benchmarks/ssr_hg002/`), against that benchmark's assembly-based truth, which is built to
    cover repeats. GATK has to be run, but only over the satellite intervals there.
  - Both are one sample. Cohort size is not exercised by this check; it is an accuracy check
    on ground the default does not call, not a check of the cohort machinery.
- *Scoring.* In repeats one variant has many spellings, and an exact match on position and alleles
  (`benchmarks/lib/compare_to_truth.sh`, bcftools) counts a differently spelled correct call as
  both a false call and a miss. RTG `vcfeval` compares haplotypes instead and is what the
  2026-10-03 GIAB review recommended (`doc/devel/ng/research/giab_unexplained_fn_2026-10-03.md`).
  It is not installed on the host or in the container. **Decision D** below.

## Decisions for the owner

**A. The stutter section under `--generic-path-only`: write it empty (recommended), or leave it out.**
Empty needs no change to the file format or its readers, keeps every digest, and is the file a fit
over zero repeats already writes, which already says in words that nothing was fitted. Leaving it
out is what the 2026-10-16 ruling says; it costs format 6, about ten readers, and moves all four
pinned digests, and gains a file one section shorter.

**B. psp format: 3.0 that still reads 2.x as default (recommended), 3.0 that refuses 2.x, or 2.2.**
2.2 lets an older build read a `call-satellites` psp as a default one, which is wrong silently.
Refusing 2.x makes every user write their psps again for an option they did not use.

**C. Accuracy data: both sets above (recommended), or only the 300× set.** The 300× set is one pvc
run and the GATK calls exist; but it is one depth, and its truth may hold few variants in
satellites. The Tandem-Repeat set adds a truth built for repeats and two lower depths, for a
GATK run over a few satellite intervals per depth.

**D. Scoring: fetch RTG Tools 3.13 into the project's `tmp/` and run `vcfeval` with the
container's Java (recommended), or score with bcftools exact matching.** The fetch is a download
from GitHub (the RTG release archive, about 50 MB). Exact matching undercounts both callers in
repeats by an amount nobody knows, which is what the comparison is about.
