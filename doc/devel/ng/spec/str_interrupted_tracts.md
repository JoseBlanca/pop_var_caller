# Interrupted STR segments: left out of the stutter fit

*Status: settled, 2026-11-07. The owner's ruling of 2026-11-06 was "leave interrupted tracts out of
the stutter fit; at calling, score each interrupted tract with the rate of its longest pure run". The
calling half was then built as a prototype and measured on genotypes (§3), changed nothing that
could be told from noise, and on 2026-11-07 the owner ruled: "we keep the fit change, but we leave
calling as in main". The measurements behind the fit change are §1; what this document changes is
§2 and §4; how we know it works is §7.*

**What this does, in words.** The repeat catalog accepts an STR segment with a few bases that break
the repeat, such as `TTTTTTTAAATTTTTTTTT`, and gives it the repeat count of its whole length: a
homopolymer of 19. Such an *interrupted* segment slips much less than a pure run of 19 bases,
because a polymerase slips along the uninterrupted stretch it is copying, and here the longest one is
9 bases. Until this change the repeat-tract fit put pure and interrupted segments of one length into
the same class and fitted one slip rate between them. From this change **the fit learns slippage
from pure segments only** (§2). Calling is unchanged: every segment, pure or interrupted, is scored
with its own class's fit (§3).

Nothing changes in the catalog, the census or the psp files, and pure segments are fitted and
called exactly as before.

## 1. Why: two kinds of segment in one class, measured

**From 17 repeats up, a homopolymer class is partly interrupted segments.** On the tomato reference
(SL4.0), at the tracts tomato2's census kept:

| repeats | segments | interrupted |
|---:|---:|---:|
| 12 to 14 | 907 | 0 |
| 17 | 40 | 4 |
| 18 | 23 | 6 |
| 19 | 30 | 17 |
| 20 | 22 | 15 |
| 21 | 16 | 8 |

The scanner holds one segment together down to a purity of about 0.78 (`ssr_repeat_scanner.md`
§3.3); a 3-base interruption leaves a homopolymer of 15 bases at 0.80, so from about that length
on a class can hold both kinds.

**The two kinds slip at different rates.** On kimura's cohort of 2,176 tomato accessions, counting at
each segment the reads at neither length of the genotype each sample's reads favour (the stutter-tail
rule's count, `parameter_prepass_ssr.md` §4.6.1), the median segment at 19 repeats has 12.0 such
reads in 100 when pure and 1.8 in 100 when interrupted. At 20 repeats it is 11.9 against 3.0.

**What the mixture does to the fit, and then to the stutter-tail rule.** The first fit of the
19-repeat class lands between the two kinds, at 8.2 slipped reads in 100. A pure segment at 12 in 100
is then nearly twice its class's rate, and the stutter-tail rule leaves out a segment at twice its
class's rate once its reads are many enough. Each round lowers the class's rate and brings more pure
segments over the line:

| round | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| slipped reads in 100, 19 repeats | 8.2 | 7.5 | 6.8 | 5.9 | 4.8 | 3.6 | 3.0 | 2.8 | 2.7 |

By the last round the rule had left out 255 of the 299 pure segments and 5 of the 324 interrupted
ones; at 20 repeats, 126 of 179 pure and 5 of 359 interrupted. Calling then scored every pure
homopolymer of 19 at 2.7 slipped reads in 100, about a quarter of what such a segment slips, and a
read one repeat short reads as a second allele rather than as stutter. *Measured with every read
group at one stutter, on branch `tail-rule-probe`; the fit with a stutter per read group left out
261 at 19 repeats against this run's 260.*

**Tomato2 shows the same fall without the rule's help**: its 19-repeat class fits at 6.3 slipped
reads in 100 against 25.9 at 16 repeats, and the rule left out none of its segments. The rounds made
the fall worse at kimura's depth; the mixture is what starts it.

**Elsewhere the rule removes what it was built for** (same run): one to five segments in 100 at
dinucleotides of 6 and 7 repeats, moving their rates by less than one part in seven; and at
hexanucleotides of 4 repeats 78 of 373 segments, all but one of them pure and of ordinary depth, at
which 8 samples in 100 show three or more lengths against 4 in 1,000 at the segments kept — the
signature of copies elsewhere in the genome. This document does not change the rule.

## 2. The fit: interrupted segments are left out before anything else

**A segment is interrupted when the catalog does not record it as perfect**:
`SsrSegment::is_perfect()` is false, meaning its bases are not the motif tiled end to end
(`repeat_catalog.md` §3.2). The catalog trims every segment to start and end on two whole copies,
so a pure segment always reads as perfect.

**Where it happens.** When the fit reads the census (`strata_of_kept_loci`, `gather_strata`), a
segment that is not perfect is skipped and counted, before the guard
(`parameter_prepass_joint_records.md` §3.3) and before the two rules of `parameter_prepass_ssr.md`
§4.6.1. Its reads teach the fit nothing about slippage.

**What it still counts toward:**

- **the depth rule's crossing share** — of the reads that reached a class's segments, the share that
  crossed one. Interrupted segments' reads count, as the guard's segments' reads already do, so the
  share describes the class's length and not which of its segments the fit kept;
- **the substitution rate.** Its two counts are kept per census section, one read group's segments
  of one class, and not per segment, so an interrupted segment's bases cannot be taken out of them.
  That leaves them as they are today: `parameter_prepass_ssr.md` §4.6 already records that an
  interruption inside a tract is charged to this rate.

**A class every segment of which is interrupted has no fit of its own**, like any class with no
reads; the run log's existing warning for an emptied class names it, with the count as the reason.
On tomato2 that is the homopolymer classes of 23 and 25 repeats, two segments each. Calling needs
them only for the candidates of a pure segment of that length, and those fall back as a missing class
does today (`str_fit_to_calling.md` §3).

**Nothing upstream changes.** The census keeps the same segments, a psp holds the same sections and
counts, and no digest moves: the kept-loci digest covers the generic path's positions only, and the
selection digest covers the catalog's settings, not how the fit groups segments.

## 3. Calling is unchanged, and why

**Calling scores every segment with the fit of its own class `(period, R)`**, as before this
document: an interrupted homopolymer of 19 gets the 19-repeat class's stutter and allele-length
prior, now fitted from pure segments only.

**The alternative, built and measured.** A prototype looked each candidate's stutter up at the class
of the segment's *longest pure run* `P`, moved by the candidate's length (`P + c − R`), and,
separately, the allele-length prior at `(period, P)`. The longest pure run is the longest stretch in
which every base equals the one a motif length further on, in whole copies; for
`TTTTTTTAAATTTTTTTTT` it is 9. Measured on genotypes, against main and against the fit change alone:

- **Tomato2**, 59 inbred accessions at about 5 reads a position, calling restricted to the 309
  single-base segments of 15 to 26 repeats in its regions. Heterozygous calls, which in inbred lines
  are mostly errors:

  | build | pure segments, 17–26 repeats (5,251 genotypes) | interrupted, 17–25 (3,304 genotypes) |
  |---|---:|---:|
  | main | 71 | 20 |
  | fit change (§2) | 63 | 20 |
  | + stutter at the pure-run class | 64 | 20 |
  | + allele-length prior at the pure-run class | 64 | 22 |

  The 15- and 16-repeat segments, which none of the changes should touch, stay at 75 heterozygous
  calls in every build. The fit change removes 8 heterozygous calls, all at pure segments; the
  stutter lookup changes nothing at interrupted segments, and the prior adds 2.
- **HG002**, GIAB tandem-repeat benchmark at 5×, 15× and 30×, 20,183 segments of which 1,954
  interrupted: every build is within 3 right genotypes and 3 false calls of main at each depth.

**So the calling change is not made.** It adds a second notion of a segment's class to calling for no
effect that the two data sets can tell from noise. The caveat is size: tomato2 holds only 56
interrupted segments of 17 repeats or more, and a cohort with many more could show an effect these
cannot. Reopen it with such a cohort's measurement, not on the argument above. The prototype is
`tmp/calling_proto_final.patch` on branch `pure-run-stutter`'s worktree, not committed.

**What this leaves open.** `read_likelihoods.md` §4.6's purity adjustment on the slippage level stays
deferred (§10 there): an interrupted segment's candidates are still scored at their own class.

## 4. The parameters file

**`[repeat_tracts_left_out].by_stratum` gains a count, `interrupted`**: how many of the class's
segments §2 left out. It goes into format 9, which is not yet on any pushed main and so has never
been written by a released build; a format-8 or older file has no such count and is read as having
left none out. Nothing at calling reads it, like the rest of the section.

**The file's note says** that interrupted segments are left out of the fit, because they slip along
their longest pure run rather than as a segment of their length.

## 5. What each end of the range gets

The change depends on the reference and not on the reads, so which segments it leaves out is the
same at every depth and cohort size.

- **One sample, few reads.** The fit has fewer segments to learn from at long homopolymers, where
  segments were already few (30 at 19 repeats in tomato's selection); the curves across classes
  carry those classes, as they carry any thin class (`str_penalised_fit.md` §2.11).
- **Thousands of samples.** The stutter-tail rule now meets classes of one kind of segment, so its
  rounds have no second population to strip away. Measured on kimura's 2,176 accessions (§7).
- **Genomes with few interrupted segments** (human at the lengths the catalog keeps): the fit loses
  few segments, and its rates barely move; HG002's 19-repeat rate went from 6.3 to 6.7 slipped reads
  in 100 at 30×.

## 6. The allele-length prior stays at the segment's own class

**Settled with the owner, 2026-11-06, and confirmed by measurement on 2026-11-07**: an interrupted
segment keeps the prior of its own class `(period, R)`. Moving it to the pure-run class changed 2 of
3,304 genotypes at tomato2's interrupted segments (§3).

## 7. How we know it works

1. **Kimura's fit with §2 in place** (2,176 accessions, every read group at one stutter): the
   homopolymer slip rates rise steadily with length, 17.1, 18.4, 20.3, 21.9, 23.4, 23.9 and 25.5
   slipped reads in 100 from 16 to 22 repeats, with no fall at 19 and 20, and the stutter-tail rule
   keeps all 299 pure segments at 19 repeats, which it had cut to 44.
2. **Tomato2's fit**: the 19-repeat class's rate moves from 6.3 slipped reads in 100 to 22.7.
3. **Genotypes** (§3): tomato2's heterozygous calls at pure homopolymers of 17 to 26 repeats fall from
   71 to 63; HG002 within ±3 of main.
4. **Unit tests**: the fit skipping a segment that is not perfect and counting it, while its reads
   still count toward the crossing share; a format-8 file with an interrupted count refused.
