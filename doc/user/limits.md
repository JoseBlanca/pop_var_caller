# Limits and refusals

## What pvc does not call

- **Structural variants**, and deletions longer than 100 bases (`--max-locus-span`).
- **Satellites**: repeats, or STR bundles, longer than 100 bases. The log counts how many were
  skipped.
- **Repeats with units longer than 6 bases.** They are called as ordinary sequence.
- **Sites where every sample matches the reference.** There is no gVCF.
- **Loci deeper than `--max-depth`** (1,000 reads a position in any read group), because such
  depth usually comes from collapsed duplications.

## Limits

| Limit | Value |
|---|---|
| ploidy | 1 to 16 |
| alleles at an SNP/indel site | 6, reference included (`--max-generic-alleles`) |
| candidate alleles at an STR | 32 |
| chromosome length in the VCF | 2,147,483,647 bases |
| chromosome length with a `.bai` index | about 537 Mb; use `.csi` above that |

## Refusals

pvc checks its inputs and outputs before reading any data, so a run that cannot finish usually
stops at once and says why. The common ones:

| Message about | Cause | Fix |
|---|---|---|
| an existing output or log | the file is already there | `--force`, or another `--output` |
| the repeat catalogue | missing, or built on another reference | `pvc repeat-catalog --reference genome.fa` (the message gives the command) |
| sort order | the alignment is not coordinate-sorted | `samtools sort` |
| `@SQ` lines | the alignment was made against another reference, or its chromosomes are in another order | align against this reference, or pass the right one |
| read groups or `SM` | no `@RG`, no `SM`, or the same read-group `ID` in two files | `samtools addreplacerg` |
| the index | no index beside the alignment | `samtools index` |
| open files | the limit is too low for this many alignments | `ulimit -n` with the value the message gives, or psp mode |
| psp files that disagree | written over different `--regions`, catalogues or repeat settings, or the same sample twice | write them again with the same settings |
| a stale psp | its summary for `estimate-parameters` is missing, or was written under other settings | the message gives the `generate-psps … --force` command for each sample |
| a parameters file | estimated on another reference, other samples or other read groups | estimate it again on this cohort |

One thing pvc does not catch: two `generate-psps` commands writing the same sample at the same
time overwrite each other. Run one command per sample.

## When a run fails part-way

- A VCF being written is named `<output>.tmp`; a failed run leaves it behind for you to delete.
- A psp being written is named `<sample>.psp.<number>.partial`; a failed run deletes it and leaves
  any earlier psp for that sample untouched.
- The scratch directory (`pvc_tmp_…` inside `--tmp-dir`) is removed at the end of a run, whether it
  succeeded or failed. A run that was killed leaves it behind.
- With `--output -`, a failed run has already sent part of a VCF downstream.
- The log records what happened. A crash inside pvc is recorded there too, and the run exits with
  status 101. Please report it, with the log, as an
  [issue](https://github.com/JoseBlanca/pop_var_caller/issues).
