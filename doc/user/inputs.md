# Preparing the inputs

pvc needs three things: the reference genome the reads were aligned to, its repeat catalogue,
and one or more alignment files. It checks all of them before it reads any data, and a run that
cannot use an input refuses to start and says why.

## The reference

A FASTA file. If it has no `.fai` index beside it, pvc builds one. pvc reads the whole FASTA once
at the start of each run and checks it against its index.

## The repeat catalogue

Build it once per reference:

```
pvc repeat-catalog --reference genome.fa
```

This writes `genome.fa.repeats.parquet` beside the reference, and every later command looks for
it there. If the reference's directory is read-only, write it elsewhere with `--output` and pass
`--catalog` to every later command. A catalogue built on another reference is refused.

`repeat-catalog` uses one thread by default, because each thread holds a whole chromosome in
memory: six threads on the human genome need more than 16 GB. One thread builds the human
catalogue in under two minutes. The catalogue is the same at any thread count.

`repeat-catalog` also has options that change how repeats are found (under *Advanced* in its
`--help`). Leave them alone. A catalogue built with stricter settings than a calling run asks
for, such as a higher `--min-copies`, is refused. But a catalogue built with other scoring
weights (`--match-reward`, `--mismatch-penalty`) is accepted without complaint, and the calls
quietly change.

## The alignments

BAM or CRAM files, with these properties:

- **Coordinate-sorted**, as declared in the header (`@HD SO:coordinate`).
- **Indexed**: `.csi` or `.bai` for BAM, `.crai` for CRAM. A `.bai` cannot index chromosomes
  longer than about 537 Mb; use `.csi` for those genomes. `--build-index-if-missing` builds a
  missing index, but on a CRAM whose slices span chromosomes it can crash, so build CRAM indexes
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

## Reads pvc ignores

These filters are fixed:

- duplicates, QC failures, secondary, supplementary and unmapped reads;
- mapping quality below 20;
- reads shorter than 30 bases;
- reads whose alignment is malformed: an insertion next to a deletion, or a deletion at either
  end of the read;
- reads where more than 1 base in 10 disagrees with the reference. Only disagreeing bases of
  quality 10 or more count as disagreements.

## Open-file limit

A one-step (direct mode) run keeps every alignment open for the whole run: 2 open files per BAM,
3 per CRAM, 2 more per sample, plus 32. The default limit is 256 open files on macOS and 1,024 on
most Linux systems. With one BAM per sample, that is reached at 56 samples on macOS and 248 on
Linux; with one CRAM per sample, at 44 and 198.

pvc checks this before starting, and if the limit is too low it says which value to set:

```
ulimit -n 4096
```

A psp-mode run keeps one file open per sample, so it reaches the limit later. See
[Calling § Large cohorts](calling.md#large-cohorts).
