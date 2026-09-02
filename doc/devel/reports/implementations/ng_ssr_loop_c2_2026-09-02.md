# C2 — the record says the locus is a repeat tract

**Date:** 2026-09-02. **Plan:** [`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md)
Milestone C step C2. **Design:** [`spec/calling_loop_ssr.md`](../../ng/spec/calling_loop_ssr.md)
§3.2; [`spec/vcf_output.md`](../../ng/spec/vcf_output.md) §6–§8.
**Modules:** [`src/ng/run/records.rs`](../../../../src/ng/run/records.rs),
[`src/ng/run/callers.rs`](../../../../src/ng/run/callers.rs) (the verdict's route),
[`src/pop_var_caller_exp/cli.rs`](../../../../src/pop_var_caller_exp/cli.rs) (a sentence that had
become false).

---

## What landed

**The motif reaches the record.** `evidence_for_output` reads the kind off the called locus's own
candidate table — a tract's table is `LocusKind::Ssr` and carries the motif — so the field that
was `None` at every locus is now `Some(TractAnnotation::new(motif))` at a tract. Everything the
file says about the repeat is written from that one value: the `STR` flag, `RU`, `PERIOD`, and
each called allele's `REPCN`.

**Selection's verdict reaches the `FILTER` column.** It used to end inside the driver; it now
travels to the record with the remapping and the leftover, and `filter_for` maps it:
`NotPeriodic` → `notPeriodic`, `Truncated` → `tooManyAlleles`, anything else → `PASS`.

## Three decisions, each with its reason

**A loop that did not settle outranks everything selection could say.** Not a preference:
`assemble_record` asserts that the filter is `EMNoConv` exactly when the loop failed to converge,
so the two cannot be carried together. A tract that was truncated *and* did not converge says
`EMNoConv`, and the truncation is visible in the record anyway — the alternatives it kept are the
alternatives it kept.

**`lowDepth` is declared and never written.** The plan lists it among the verdicts selection
mints. It is not one: production refuses a tract whose *cohort-summed* depth is under ten, and ng
does not port that gate — depth is asked once, upstream, per sample, by the merge's keep rule, and
spec §6 is explicit that there is no depth verdict on this path. The vocabulary stays in the
header because a file written by an older caller can carry it.

**A truncated SNP/indel locus is still a `PASS`, unchanged.** `tooManyAlleles` is spec §8's
*tract* filter; the ordinary path's cap cuts the lowest-ranked alternatives and calls the locus
over the rest. What changed is only that a tract's truncation now reaches the column.

## `REPCN` needed no wiring, and that is worth stating rather than passing over

The plan asks for "the per-allele repeat counts" to reach `REPCN`. They already do, by a different
route: the encoder derives a called allele's repeat copies from the record's own bases and the
motif, so it cannot drift from the sequence the same record writes.

**That leaves two producers of one integer**, and D3's report is where the hazard was named: the
genotype prior indexes its length spectrum by selection's count, and a disagreement would put a
candidate's prior mass at one length while the file reported another, with nothing failing. They
agree because both are the same floor division, and
`the_records_repeat_count_is_the_one_the_prior_is_given` is what holds them to it — over a whole
unit, a part unit, and a sequence shorter than one unit, which is where a floor division and a
rounded one part company.

## The bcftools round-trip, measured

A real run: HG002 at 30× over the first 200 intervals of the Tier BED, through ng's own catalog.
**78 records written, 17 of them repeat tracts.** Read and rewritten by `bcftools view` 1.16:

- **the tract fields are byte-identical over all 17 records** — `bcftools query -i 'STR=1' -f
  '%CHROM %POS %RU %PERIOD [%REPCN] %FILTER\n'` gives the same output before and after, so
  `bcftools` parses `RU`, `PERIOD` and the per-sample `REPCN` as declared;
- **the only column that moves is `INFO`, and only its Float fields** — `MQALT=70.00` comes back
  as `MQALT=70`, and the same for `AF`, `ABPEN`, `SPPEN`, `MQREF`, `MQDIFF`. Every one of those is
  on SNP/indel records too. It is `bcftools` normalising a float's trailing zeros, not anything
  the tract fields introduced.

One of the 17 is a `1/2` tract — two alternative lengths, `REPCN=15,16` — which is the case that
would break first if `REPCN` were built from the unsorted candidate order rather than from the
genotype's.

## The command's own description had become false

`call-from-alignments --help` opened with *"Repeat tracts are analysed and NOT called: every locus
goes down the SNP/indel path"*. It is now the opposite, so the sentence is replaced rather than
left to mislead: tracts are called through their own model, and what is still set aside is a
repeat *cluster* too tangled to have clean flanks.

## Tests — 4 new

| test | what it pins |
|---|---|
| `a_tract_record_carries_its_motif_and_an_ordinary_one_carries_none` | the motif reaches a tract's record, and an ordinary locus over repeat-looking bases is not a tract record |
| `a_tracts_filter_comes_from_selection_unless_the_loop_did_not_settle` | both tract verdicts reach `FILTER`; non-convergence outranks them; a truncated SNP/indel locus is still a `PASS` |
| `the_records_repeat_count_is_the_one_the_prior_is_given` | the two producers of the repeat count agree, including on a part unit and a sub-unit sequence |
| (the encoder's existing `REPCN` tests) | unchanged and still green — the encoding half was already built |

## Mutation testing — three run, three killed

| mutation | outcome |
|---|---|
| the motif never reaches the record | killed |
| the filter ignores selection's verdict | killed |
| a tract verdict outranks non-convergence | killed |

## Validation

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` 6,011 passed / 0 failed / 14 ignored, all in the container. The `bcftools`
round-trip above ran in the container too, at `bcftools 1.16` / `htslib 1.16`.
