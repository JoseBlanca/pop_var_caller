# pvc user guide

For what pvc is and how to install it, see the [project README](../../README.md). For the
statistical models behind the calls, see [doc/methods/](../methods/) *(in preparation)*.

## Contents

- [Preparing the inputs](inputs.md): the FASTA reference and how to build its repeat catalogue,
  what pvc requires of the BAM or CRAM files, which reads it ignores, and the open-file limit.
- [Calling](calling.md): the two ways to call, direct mode and psp mode, and how to choose
  between them; adding samples to a cohort; large cohorts; outputs, logs and scratch space.
- [Estimated parameters](parameters.md): the parameters pvc estimates from your reads and uses as
  priors during calling, their built-in values, and how to read, edit and reuse the parameters
  file.
- [Options](options.md): ploidy, inbreeding, which sites are written, depth and locus-size limits,
  and what counts as a repeat, with when to change each.
- [The VCF](vcf.md): the header, QUAL and GQ, every INFO, FORMAT and FILTER field, and how STRs
  are written.
- [Limits and refusals](limits.md): what pvc does not call, its hard limits, why a run can refuse
  to start, and what a failed run leaves behind.
- [Estimating contamination](../help/estimate_contamination.md).

`pvc <command> --help` lists every option of a command with its default.

## Definitions

- **psp file**: a per-sample pileup file that holds what one sample's reads showed at every
  position. It is written once from the sample's alignments, and lets you estimate parameters and
  call without reading the alignment files again.
- **Direct mode** and **psp mode**: the two ways to run pvc. Direct mode calls straight from the
  alignments; psp mode goes through psp files. [Calling](calling.md) compares them.
- **Parameters file**: a text file (TOML) with the numbers pvc estimates from your data and uses
  as priors during calling, such as each sample's inbreeding coefficient. See
  [Estimated parameters](parameters.md).
- **Read group**: one `@RG` line of a BAM or CRAM header, usually one library sequenced on one
  lane. pvc estimates base-quality calibration and contamination per read group.
- **STR** (short tandem repeat): a stretch of DNA where a unit of 1 to 6 bases repeats, such as
  `ACACACACAC`. pvc calls an STR as one locus whose alleles are whole repeat sequences.
- **STR bundle**: repeats lying within 15 bases of each other, close enough that one read spans
  several. pvc calls them together, as one locus.
- **Satellite**: a repeat or bundle longer than 100 bases. Reads of ordinary length cannot span
  one, so pvc does not call it.
- **Repeat catalogue**: the list of every tandem repeat in the reference, made once per reference
  by `pvc repeat-catalog`.
- **Stutter**: when DNA is copied inside a repeat, the copy often gains or loses a whole unit, so
  a sample carrying 12 copies might show some reads with 11 or 13.
