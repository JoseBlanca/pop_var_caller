# Plan — no `*` allele: an STR record widened over a flank absorbs the generic records it overlaps

*2026-10-29. Branch `str-edge-no-star`. Implements the amendment of
[`../spec/str_edge_deletion_alleles.md`](../spec/str_edge_deletion_alleles.md) §6 (2026-10-29).
The owner's direction: pvc never writes `*`.*

## What changes, in one paragraph

Today `src/run/edge_deletion_gathering.rs` rewrites each generic record that overlaps a widened STR
record: copies the STR allele deletes become `*`, and a generic record left with only the
reference and `*` is dropped. After this work the same pass **merges** each such generic record
into the STR record: one record over the union of their spans, whose alleles are haplotypes, the
flank allele of each sample paired with its STR allele by the STR allele's own sequence, by
homozygosity, or else by the pairings the cohort shows (spec §6.2). No other module calls
differently; the change is confined to the pass, the record it builds, and the documents.

## Size of the change, measured on main's build (5ffb2444)

| run | generic records written with `*` | samples whose pairing comes from the cohort |
|---|---|---|
| stuttering benchmark, 1 sample, 8 runs | 0 to 2 a run | 0 |
| stuttering benchmark, 50 samples, 8 runs | 9 to 36 a run | 30 to 40 a run |
| tomato, 62 accessions, about 3× | 97 | 15 |

Records the pass drops are not counted by the run today; step 1 adds the count.

## Steps — one commit each, touched tests run at each

1. **Count before changing anything.** Write the pass's counts to the run log (`pvc.log`): copies
   taken in agreement, copies overridden, generic records dropped. Run the oracle and the
   stuttering benchmark's 50-sample PCR run at 15× on this build to have the "before" numbers,
   dropped records included. *Stop: report the counts.*
2. **The merge** (spec §6.1–§6.7), replacing `reconcile` in `edge_deletion_gathering.rs`:
   pairing (steps 1–3 of §6.2), the allele table, `GQ`, `AD`/`DP`/mapping-quality pools, `REPCN`
   through `TractWidening::flank_bases_per_allele` (each merged allele's flank bases are its length
   less its STR spelling's), `QUAL`/`FILTER`/`AF`. New counts: generic records merged, samples
   paired from the cohort, samples no-called by the merge. The three existing tests are rewritten
   to the merged expectation; new tests:
   - a sample heterozygous at both records, paired from three other samples' homozygous pairings;
   - the same at one sample: the first pairing, `GQ` 3;
   - a generic no-call keeps the STR genotype; an STR no-call is `./.` with its `AD`;
   - `AD` of an STR allele split over two merged alleles goes to `DP − ΣAD`;
   - `REPCN` equal for two alleles differing only in the flank;
   - a generic record on each flank of one STR record;
   - a generic record extending left of the widened span, with an empty allele padded by the
     generic record's own left base (the only case where the merged span starts before the STR
     record and needs a padding base; a generic record with an empty allele always carries it).
   *Stop: report.*
3. **`*` becomes unwritable.** Delete `SPANNING_DELETION`; the encoder's `*` branch
   (`src/vcf/encode.rs:189`) becomes a refusal, so a future path that produced one fails a test
   instead of writing it.
4. **Documents**: `vcf_output.md` §3 (records never claim the same base, now without exception)
   and wherever it names the gathering; `doc/user/vcf.md` "STR records";
   `doc/user/overlapping_variants.md` "What pvc writes" and the README's comparison-table cell for
   pvc. The last two live on branch `docs-overlapping`: rebase onto main once it is merged, or
   edit them there.
5. **Checks before merging** (the brief's list):
   - no `*` in any VCF the oracle and the benchmarks write (`bcftools view -H -i 'ALT="*"'`);
   - the identity oracle (`scripts/promote_ng_oracle.sh`): every changed record lies at a widened STR
     record, each listed with why it changed;
   - the stuttering benchmark grid, 16 runs, true stutter, against main's binary, overall and for
     the pairs carrying a deletion across a segment's edge and the pairs with a SNP beside the
     segment, each number labelled with its run (samples, PCR or not, depth);
   - tomato, 62 accessions: records changed, and the 111 flank genotypes §6.2 expects lost.
   Run one job at a time on the laptop's 6 fast cores.

## Decisions taken in the spec that the owner should confirm

These are written into §6 as recommendations; each is reversible inside step 2.

- **Pairing from the cohort where the reads cannot say** (§6.2, step 3), with `GQ` at most 3 when
  the cohort has nothing to offer. Alternatives: no-call every such sample (loses its `REPCN`
  too), or always pick the first pairing at `GQ` ≤ 3 (no cohort counting).
- **A sample no-called at the STR record is no-called in the merged record**, losing a flank
  variant's call where the generic record had one (111 on tomato, 2 on the stuttering benchmark).
- **The cap of 32 alleles does not apply to the merged record** (largest measured: 12).
- **`AD`, `DP` and mapping qualities from the STR record**; the generic record's reads are not
  written.
