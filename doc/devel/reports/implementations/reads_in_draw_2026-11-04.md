# Reads in the pass's haplotype draw — implementation report

*2026-11-04, branch `reads-in-draw`. Builds
[`neighbour_evidence_reads_in_draw.md`](../../ng/spec/neighbour_evidence_reads_in_draw.md), whose
§8 the owner ruled the same day: the draw given the path first, prototyped behind a switch and
measured before the tests, `HET_MASS_MIN` at ½ with ¼ and ¾ measured beside it.*

## What was built

Along each path a round draws, each sample's read fragments now set the relative orientation of the
heterozygotes they link (spec §2.2). The step sits after the path is drawn, in `draw_haplotypes`
(`src/neighbour_evidence/model.rs`):

1. every site's joint of the two copies' alleles is computed as before (`joint_along`);
2. `orient_by_reads` keeps the sites whose heterozygous mass is at least `HET_MASS_MIN` = ½ and
   that are not hidden, and walks them in order; consecutive ones the sample's fragments link form
   a group;
3. `orient_group` draws each group's orientations from the two-state chain of spec §2.2 and moves
   each grouped site's heterozygous mass onto the drawn order;
4. each copy's haplotype is the joint's margin, as before (`copies_of`).

`refine.rs` hands the window's fragments to `solve_round`, which gives the draw a closure returning
`score_read_link`'s `L` for two sites. A new random stream, `StreamPurpose::Orientation`, one a
sample a round, is read by the round's own path and then its extra paths. The switch is
`Phasing::with_reads_in_draw`, set from `NeighbourEvidenceDecision::reads_in_draw`: on wherever the
pass runs, off only in the one test that checks the link test's own `GQ` rule (below). No option or
header value changes.

**Without fragments nothing changes.** The draw step is skipped when the window has no fragments,
and a sample with fragments that link nothing keeps its joints bit for bit and reads nothing from
the new stream (tests `reads_that_link_nothing_change_no_bit`, `an_unlinked_site_keeps_its_joint`).
The margins are summed in the same order as before, so the haplotypes are bit-identical.

## Measurements

All on the founder simulator's cohorts of the read phase report
(`read_phase_2026-11-03.md`), 3× unless marked, called with refinement on by the prototype, once
with the draw switched off — main's behaviour, the prototype's environment switch, since removed —
and once with it on; the timings compare main's build (`0019f331`) with the finished branch. Links are scored with `read_phase_score.py`
against the simulator's truth; genotypes with `founder_grid.py`'s scoring. Scripts and outputs:
`tmp/reads_in_draw/` (`run.sh`, `genotypes.py`, `timing.py`). One seed each.

### The cohorts this change is for: outbred, SNPs 1 in 100 bases apart

50 samples at 3×, 2 Mb. A *link* is two consecutive phased records of one sample in one phase set;
it is wrong when their written orientations disagree with the truth's.

| founders | wrong links, main → this branch | true heterozygotes phased | wrong genotypes where the sample has reads | calls without reads: made, wrong |
|---|---|---|---|---|
| 64 | 1,059 of 30,237 (3.5%) → 916 of 45,443 (2.0%) | 41,020 → 61,688 of 151,380 | 35,681 → 33,656 of 468,796 | 3,816, 164 → 3,967, 168 |
| 16 | 40 of 44,563 (0.09%) → 30 of 57,168 (0.05%) | 59,147 → 74,989 of 145,477 | 17,157 → 14,476 of 419,457 | 8,967, 40 → 9,340, 45 |

**Phase and genotypes both improve.** At 64 founders the branch writes half as many links again,
and fewer of them are wrong in number, not only in share. Genotypes improve by 6% and 16%: other
samples now copy haplotypes that carry their neighbours' read phase, so the templates are better
beyond the records reads link — the hypothesis of spec §1.1. **The 1-in-100 ceiling is not met at
64 founders**: 2.0% remains. Phased heterozygotes with a wrong genotype rise with the phasing,
from 504 to 950 at 64 founders (1.2% to 1.5% of those phased).

### `HET_MASS_MIN`

| threshold | 64 founders: wrong links | heterozygotes phased | wrong genotypes with reads | 16 founders: wrong links |
|---|---|---|---|---|
| ¼ | 917 of 43,824 (2.1%) | 59,650 | 34,026 | 31 of 55,806 |
| **½** | **916 of 45,443 (2.0%)** | **61,688** | **33,656** | **30 of 57,168** |
| ¾ | 967 of 44,136 (2.2%) | 59,915 | 34,007 | 40 of 56,049 |

½ is best on every column; it stays.

### No harm elsewhere: the nine cohorts of the phase output report

SNPs 1 in 1,000 bases apart, the simulator's default.

| cohort | wrong links, main → branch | wrong genotypes with reads, main → branch |
|---|---|---|
| 50 at 3×, F 0.5, 16 founders | 3 of 3,669 → 4 of 3,681 | 615 → 607 of 47,443 |
| 50 (10 at 30×), F 0.5, 16 founders | 4 of 4,407 → 3 of 4,366 | 409 → 427 of 48,084 |
| 50 (10 at 30×), F 0.8, 16 founders | 1 of 2,139 → 1 of 2,149 | 197 → 192 of 47,579 |
| 50 at 3×, outbred, 16 founders | 9 of 2,186 → 5 of 2,311 | 2,756 → 2,697 of 46,915 |
| 200 (40 at 30×), F 0.8, 64 founders | 4 of 9,367 → 6 of 9,388 | 1,296 → 1,321 of 221,346 |
| 50 (10 at 30×), F 0.8, 64 founders | 6 of 1,187 → 6 of 1,228 | 572 → 572 of 52,798 |
| 10 at 3×, F 0.95, 64 founders | none phased either way | 165 → 165 of 7,912 |
| 10 at 3×, outbred, 64 founders | 6 of 23 → 6 of 35 | 2,174 → 2,176 of 8,170 |
| 50 at 3×, outbred, 64 founders | 64 of 1,320 → 60 of 1,475 | 4,528 → 4,486 of 52,206 |
| **all nine** | **97 of 24,298 → 91 of 24,633** | **12,712 → 12,643** |

Two inbred mixed-depth cohorts gain wrong genotypes, 18 of 409 and 25 of 1,296 (4% and 2%), inside
the spread between seeds of about 5% (phasing spec §1); every other cohort is level or better. Calls
where the sample has no reads move by at most 6 a cohort.

### Cost

Wall time of the whole `call-from-psps` run, 6 threads, three runs each alternating the builds,
on an M-series laptop:

| cohort | main | this branch |
|---|---|---|
| 50 samples, dense, 64 founders | 10.7, 11.0, 11.1 s | 11.3, 11.4, 11.3 s |
| 200 samples (40 at 30×), F 0.8 | 29.3, 29.4, 29.3 s | 29.5, 29.4, 29.4 s |

About 3% on the dense cohort, where groups are many, and nothing measurable at 200 samples with
SNPs ten times sparser.

## Tests

In `src/neighbour_evidence/model/tests.rs`, the spec §6 list:

- `reads_alone_orient_two_sites_alike_in_their_templates` — equal own weights and one fragment's
  `L` give the same orientation with probability 0.98, the first site a fair coin;
- `three_fragments_overturn_confident_templates_and_two_do_not` — against templates favouring the
  other order by about 39,600 to 1, two fragments turn it about 6 times in 100, three about 75;
- `hidden_and_homozygous_sites_stay_out_of_a_group` — both keep their joints, and the sites either
  side are still linked across them;
- `the_written_joint_keeps_homozygous_mass`, `a_one_two_heterozygote_is_read_as_alleles_one_and_two`;
- `an_unlinked_site_keeps_its_joint`, `reads_that_link_nothing_change_no_bit`,
  `reads_change_the_haplotypes_and_not_the_walk`.

**Thread invariance** is the stage test `the_pass_with_reads_writes_the_same_bytes_at_any_thread_count`
(one thread against eight, refinement on, fragments at every heterozygote), which now also checks
that its output differs from the same run with reads kept out of the draw. Spec §6 point 5 named
`mode_equivalence.rs`; its cohort has two samples, too few for refinement, so the stage test is the
one that exercises the draw.

**One test changed its premise.** `a_column_reads_newly_phase_has_gq_30_or_more` asserted that every
column phased with fragments and not without them reached `GQ` 30, the link test's rule for links
resting on reads. With reads in the draw the population's own phase uses them, and may phase a
heterozygote below 30 (the first the test met was at `GQ` 14), which the population is allowed to
do. The test now keeps the fragments out of the draw, so it checks the link test alone.

## What is left

- **The 64-founder dense cohort is at 2.0%**, above the ceiling. The next step is spec §2.4's exact
  draw, reads in the path; an owner decision, not part of this branch.
- **Real data** (spec §6 point 6): the tomato read check and HG002's segmental duplications, after
  merging.
- The same measurement at a second seed was not run: no decision here rests on a difference small
  enough to need it.
