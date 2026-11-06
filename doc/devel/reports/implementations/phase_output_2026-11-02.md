# Writing the phase the neighbour-evidence pass estimates: the build

*Implementation report, 2026-11-02, branch `neighbour-pass`, steps `1ab612e2` to `e31e568f`. Spec:
[`neighbour_evidence_phase_output.md`](../../ng/spec/neighbour_evidence_phase_output.md); the
prototype check that set its rule:
[`phase_blocks_check_2026-11-02.md`](../../ng/research/neighbour_evidence/phase_blocks_check_2026-11-02.md).
Scored with `doc/devel/ng/research/neighbour_evidence/phase_vcf_score.py`; outputs in the
worktree's `tmp/phase_check/rust/`.*

## Summary

**The pass now writes `GT` with `|` and a `PS` field where it is confident, and on every inbred or
mixed-depth simulated cohort at most about 1 link in 200 inside a phase set is wrong.** Genotypes and
sites are unchanged; the cost is 10% at 200 samples.

- **Inbred and mixed-depth cohorts (16 or 64 founders, 50 or 200 samples):** 44% to 77% of the true
  heterozygous SNPs phased, 0.05% to 0.51% of the links inside phase sets wrong.
- **Outbred cohorts:** 50 samples of 16 founders, 24% phased, 0.39% wrong. 50 samples of 64
  founders, 13% phased, 4.7% wrong. 10 samples of 64 founders, 21 SNPs phased, 7 of 13 links wrong.
  The last two are the confidently wrong case the prototype found; the spec keeps it open.
- **Phase sets are short:** a median (N50) of 2 to 12 records, 14 to 77 kb.

## What was built

| step | commit | what |
|---|---|---|
| 1 | `1ab612e2` | `model::walk_and_draw_paths`: the round's own path as before, and N more from the same forward probabilities on a stream of their own (`StreamPurpose::ExtraPaths`) |
| 2 | `7686b276` | `phase.rs`: the consensus of the averaged rounds (two-state alignment, `f32`), each path's two bits at every allele pair, the link test over all 36 paths, blocks; `RefinedRecord::phases` |
| 3 | `804cf7be` | `line::rewrite_line` writes `a|b`, `PS` last in `FORMAT`, the repeat counts in `GT` order |
| 4 | `4964bcc3` | header (`GT`'s description, `PS`, `extra_paths`, `lean_min`, `path_agreement`); `--phase-min-path-agreement` on both calling commands and `refine-genotypes`; `refine-genotypes` refuses phased input |
| 5 | `e31e568f` | joins across windows: the stage keeps each sample's last written block; `phase::written_in_core` joins on 3 agreeing shared records |

## Measured

`pvc refine-genotypes` on each cohort's calling VCF (`calls_off.vcf.gz`), 6 threads. SNP records
only: an indel or a repeat can be spelled differently in the call and in the truth, and a match by
sequence then counts a right genotype as wrong. On the first cohort that mismatch made 398 of 408
"wrong" phased non-SNP columns.

| simulated cohort | true het SNPs phased | phased with a wrong genotype | links inside phase sets wrong | N50 records (kb) |
|---|---|---|---|---|
| 50 (10 at 30×, 40 at 3×), inbreeding 0.5, 16 founders | 5,534 of 7,532 (73%) | 10 | 4 of 4,375 (0.09%) | 8 (51) |
| same, inbreeding 0.8 | 2,572 of 3,323 (77%) | 10 | 2 of 2,174 (0.09%) | 11 (67) |
| same, inbreeding 0.8, 64 founders | 1,521 of 3,476 (44%) | 7 | 6 of 1,184 (0.51%) | 7 (52) |
| 200 (40 at 30×), inbreeding 0.8, 64 founders | 11,082 of 15,411 (72%) | 50 | 5 of 9,311 (0.05%) | 12 (77) |
| 50 at 3×, inbreeding 0.5, 16 founders | 4,763 of 7,510 (63%) | 12 | 5 of 3,505 (0.14%) | 6 (37) |
| 10 at 3×, inbreeding 0.95, 64 founders | 0 of 207 | 0 | — | — |
| 50 at 3×, outbred, 16 founders | 3,958 of 16,416 (24%) | 8 | 8 of 2,028 (0.39%) | 2 (14) |
| 50 at 3×, outbred, 64 founders | 2,237 of 17,244 (13%) | 43 | 55 of 1,170 (4.70%) | 2 (17) |
| 10 at 3×, outbred, 64 founders | 21 of 3,161 (1%) | 0 | 7 of 13 | 3 (19) |

**Against the prototype** with the same rule (the check's 36-path rule): Rust phases more on seven
of the nine cohorts — 72% against 46% at 200 samples, 13% against 7% on the 64-founder outbred
cohort — fewer on the 10-sample outbred one (21 SNPs against 48), and none on the 10-sample inbred
one, as the prototype. Its links are wrong less often on four cohorts (0.05% against 0.28% at 200
samples, 4.7% against 7.8% on the 64-founder outbred one) and more often on four: three under the
1% limit (0.09% against 0.05% on the first cohort, 0.51% against 0.35%, 0.14% against 0.10%) and
the 10-sample outbred one far over it (7 of 13 against 15 of 42). Its phase sets are
about half as long (N50 8 against 13 records on the first cohort). The prototype fixed the copying length
at 300 kb; Rust fits it per sample from {30, 300} kb, and a 30-kb fit breaks links over 30 kb.
Which of the two causes each difference was not separated.

**Deep and shallow samples:** of the wrong links on the mixed-depth cohorts, 1 is between deep
samples (in 3,808 links) and 16 between shallow ones (in 13,236).

**Unchanged:** every record's first eight columns — sites, `ALT`, `INFO` with `AC`/`AN` — are
identical to the build before this work (`d31e4477`) on all nine cohorts. One thread and six write
the same bytes on two cohorts (one inbred, one outbred). `bcftools view` and `query` read the
output, including multi-allelic phased genotypes such as `0|3`.

**Joins:** the stage test (20 outbred samples of 4 founders at 4×, 1 Mb in 200-kb windows) has 29
links crossing a window boundary, none wrong, and 9 wrong in 5,494 links in all. The simulated
cohorts' contigs are 1 Mb and the default window is 2 Mb, so joins did not occur in the table above.

**Time** (single runs, 6 threads): 200 samples 15.9 s → 17.5 s; 50-sample cohorts 3.2–3.7 s →
3.4–4.0 s.

## Not done

- **Memory**: the paths' bits take 16 bytes at every allele pair, not the spec's 5 — 1.9 GB at 3,000
  samples in a floor-sized window. Packing them is straightforward; left until a cohort that size
  runs.
- **Joins on real windows**: measured only in the stage test; the simulated contigs are shorter
  than a window.
- **The outbred failure**: unchanged from the prototype check; the spec leaves it open (§6, §8 read
  evidence through chain ids).
- **`whatshap compare`**: not installed here; the scorer above stands in for it.
