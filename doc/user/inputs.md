# Preparing the inputs

The required input files are:

- the reference genome the reads were aligned to,
- its repeat catalogue, and
- one or more alignment files.

## The reference

A FASTA file. If it has no `.fai` index beside it, pvc builds one.

## The repeat catalogue

Build it once per reference:

```
pvc repeat-catalog --reference genome.fa
```

This writes `genome.fa.repeats.parquet` beside the reference, and every later command looks for
it there. If the reference's directory is read-only, you can write it elsewhere with `--output`
and pass `--catalog` to every later command.

`repeat-catalog` uses one thread by default, because each thread holds a whole chromosome in
memory. `--threads` sets more; the catalogue is the same at any thread count. Measured:

| genome | threads | time | peak memory |
|---|---|---|---|
| human, GRCh38 | 1 | about 40 s | about 1.4 GB |
| human, GRCh38 | 6 | about 12 s | about 4.7 GB |
| tomato, SL4.0 (782 Mb, 13 chromosomes; 6,441,172 repeats) | 1 | 17 s | 484 MiB |

## The alignments

BAM or CRAM files, with these properties:

- **Coordinate-sorted**, as declared in the header (`@HD SO:coordinate`).
- **Indexed**: `.csi` or `.bai` for BAM, `.crai` for CRAM. A `.bai` cannot index chromosomes
  longer than about 537 Mb; use `.csi` for those genomes. `--build-index-if-missing` builds a
  missing index, but it can crash on a CRAM whose slices span chromosomes, so build CRAM indexes
  with `samtools index` instead.
- **Aligned to the same reference**: the `@SQ` lines must match the FASTA in chromosome names,
  lengths and order, and in MD5 checksum where both carry one.
- **Read groups with sample names**: every file needs at least one `@RG` line with an `SM` tag.
  `samtools addreplacerg` adds one. The `SM` value is the sample's name in the VCF.
- **Read-group IDs unique across the run**: two files whose read groups share an `ID` are refused
  together. This catches cohorts where one pipeline gave every file the same ID, such as `ID:1`.
  Rename them with `samtools addreplacerg`.

Two files with the same `SM` are taken as the same sample, and their reads are pooled. One file
may hold several samples. If a read group has no `LB` (library) tag, pvc names its library from
the sample, the read-group `ID` and the file name.

In psp mode the sample name becomes a file name, so a name that is empty, is `.` or `..`, or
contains `/` is refused. Dots elsewhere in a name, as in `S.1`, are fine.

## Regions

`--regions` takes a BED file: the first three columns are the chromosome, as named in the FASTA,
and the start and end, counted from 0 with the end excluded (the BED convention). Further columns,
blank lines, `#` comments and `track` or `browser` lines are ignored. Intervals that overlap or
touch are merged, and an end past the chromosome's end is cut back to it. A chromosome not in the
reference is refused.

## Read filters

pvc ignores these reads:

- duplicates, QC failures, secondary, supplementary and unmapped reads;
- mapping quality below 20;
- reads shorter than 30 bases;
- reads whose alignment is malformed: an insertion next to a deletion, or a deletion at either
  end of the read;
- reads where more than 1 base in 10 disagrees with the reference. Only disagreeing bases of
  quality 10 or more count as disagreements.

## Open-file limit

A one-step (direct mode) run keeps every alignment open for the whole run: 4 open files per BAM,
11 per CRAM, 4 more per sample, plus 32, whatever the number of threads. The default limit is 256
open files on macOS and 1,024 on most Linux systems. With one BAM per sample, that is reached at 28
samples on macOS and 124 on Linux; with one CRAM per sample, at 14 and 66.

pvc checks this before starting, and if the limit is too low it says which value to set:

```
ulimit -n 4096
```

A psp-mode run keeps one file open per sample, so it reaches the limit later; it checks the same
way before starting. See
[Calling § Large cohorts](calling.md#large-cohorts).
