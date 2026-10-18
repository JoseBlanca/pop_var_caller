# pvc user guide

This guide is for people running pvc: how to prepare the inputs, how to call a cohort, which
options to change and when, and how to read what pvc writes. For what pvc is and how to install
it, see the [project README](../../README.md). For the statistical models behind the calls, see
[doc/methods/](../methods/) *(in preparation)*.

| I want to… | Read |
|---|---|
| check my reference and alignments will be accepted | [Preparing the inputs](inputs.md) |
| call one sample or a cohort | [Calling](calling.md) |
| add samples to a cohort I have already called | [Calling § Adding samples](calling.md#adding-samples-to-a-cohort) |
| spread a large cohort over several machines | [Calling § Large cohorts](calling.md#large-cohorts) |
| understand or edit the numbers pvc estimates from the data | [Estimated parameters](parameters.md) |
| change ploidy, filters, depth limits or what counts as a repeat | [Options worth changing](options.md) |
| understand a VCF field, a FILTER value or an STR record | [Reading the VCF](vcf.md) |
| know what pvc does not call, or why a run refused to start | [Limits and refusals](limits.md) |
| estimate contamination on its own, to compare with another tool | [Estimating contamination](../help/estimate_contamination.md) |

`pvc <command> --help` lists every option of a command with its default. This guide does not
repeat that list. It says which options matter and when to change them.

## Words this guide uses

- **STR** (short tandem repeat): a stretch of DNA where a unit of 1 to 6 bases repeats, such as
  `ACACACACAC`. pvc calls an STR as one locus whose alleles are whole repeat sequences.
- **Stutter**: when DNA is copied inside a repeat, the copy often gains or loses a whole unit, so
  a sample carrying 12 copies shows some reads with 11 or 13. pvc expects this and does not call
  such reads as new alleles.
- **STR bundle**: repeats lying within 15 bases of each other, close enough that one read spans
  several. pvc calls them together, as one locus.
- **Satellite**: a repeat or bundle longer than 100 bases. Reads of ordinary length cannot span
  one, so pvc does not call it.
- **Repeat catalogue**: the list of every tandem repeat in the reference, made once per
  reference by `pvc repeat-catalog`.
- **psp file**: a per-sample file that holds what one sample's reads showed at every position.
  It is written once from the sample's alignments, and lets you estimate parameters and call
  without reading the alignments again.
- **Direct mode** and **psp mode**: the two ways to run pvc. Direct mode calls straight from the
  alignments; psp mode goes through psp files. [Calling](calling.md) compares them.
- **Parameters file**: a text file (TOML) with the numbers pvc estimates from your data before
  calling, such as each sample's inbreeding coefficient. See [Estimated parameters](parameters.md).
- **Read group**: one `@RG` line of a BAM or CRAM header, usually one library sequenced on one
  lane. pvc estimates base-quality calibration and contamination per read group.
