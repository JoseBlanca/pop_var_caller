# Interrupted STR segments: left out of the stutter fit, scored by their longest pure run

*Status: proposal, 2026-11-06. Owner's ruling of the same day: "leave interrupted tracts out of the
stutter fit; at calling, score each interrupted tract with the rate of its longest pure run". The
measurements behind it are §1; what this document changes is §2 to §5; how we know it works is §7.*

**What this does, in words.** The repeat catalog accepts an STR segment with a few bases that break
the repeat, such as `TTTTTTTAAATTTTTTTTT`, and gives it the repeat count of its whole length: a
homopolymer of 19. Such an *interrupted* segment slips much less than a pure run of 19 bases,
because a polymerase slips along the uninterrupted stretch it is copying, and here the longest one is
9 bases. Today the repeat-tract fit puts pure and interrupted segments of one length into the same
class and fits one slip rate between them, and calling scores both with it. From this change:

- **the fit learns slippage from pure segments only** (§2);
- **calling scores an interrupted segment's reads at the slip rate of a pure segment as long as its
  longest uninterrupted run** — for `TTTTTTTAAATTTTTTTTT`, a homopolymer of 9 (§3).

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

## 3. Calling: an interrupted segment is scored as its longest pure run

**The longest pure run** of a segment is its longest stretch of bases in which every base equals the
one a motif length further on, counted in whole copies: a stretch of `L` such bases holds `⌊L / k⌋`
copies of a motif of `k` bases. It is read off the reference bases of the segment, which calling
already holds as allele 0 (`candidate_alleles_ssr.md`). Written this way it does not depend on where
the motif's phase starts, so an interruption that shifts the phase (`CACACAT` then `ACACAC`) leaves
two runs of the same motif rather than one run and a mismatch. For a perfect segment the stretch is
the whole segment and the run is the segment's own repeat count.

**The stutter of each candidate is looked up at the pure-run class, moved by the candidate's
length.** Today a candidate of `c` repeats at a segment is scored with the slippage of class
`(period, c)` (`read_likelihoods.md` §4.4). At a segment of reference count `R` whose longest pure
run is `P`, the candidate is scored with the slippage of class `(period, P + c − R)`:

- the reference allele itself, `c = R`, slips as a pure run of `P`;
- a candidate one repeat longer, `c = R + 1`, slips as a pure run of `P + 1`: the extra copy is taken
  to sit in the longest run, which is where slippage most likely put it.

**Why the candidate's length and not its own longest run.** A candidate's bases are the reference's
with whole copies added or removed, so its longest run is `P + c − R` exactly when the change sits in
the longest run, and calling does not know where a length change sits. The formula takes the most
likely place, and it needs no new sequence work per candidate.

**When `P + c − R` is below one**, or names a class with no fit, the candidate falls back as a
missing class does today: the shipped stutter, marked `defaulted` (`str_fit_to_calling.md` §3).
At tomato's interrupted homopolymers, `P` is 7 to 10 copies shorter than `R`, so only a candidate
that has lost most of the segment reaches below the fitted classes.

**What else the candidate's class sets is unchanged.** The substitution rate is still looked up at
`(period, c)`, because it was fitted over interrupted and pure segments alike (§2). The allele-length
prior is still the segment's own class `(period, R)` (§6). A
segment's read-group multipliers are the read group's at the motif length, so they apply as before.

**Bundles are not touched.** A bundle locus already looks each of its pure reference runs up at its
own run's class (`str_bundle_loci.md` §7.1).

**Where it goes in the code.** One place: the slippage lookup in `TractScoringFits::gather_for_locus`
(`calling/inference/repeat_tract_parameters.rs`), which today asks `StratumFits::at(read_group,
period, candidate.repeat_count)`. It asks with `P + c − R` instead, `P` computed once a locus from
the reference allele's bases. Candidates' repeat counts, the prior's centring and every other use of
`c` stay as they are.

**This fills the purity adjustment `read_likelihoods.md` §4.6 left open** ("the level a candidate
gets is the stratum's, unadjusted for that candidate's own purity"), deferred there to "the
interrupted-repeat work" (§10). It is narrower than that work: it says nothing about alleles that
differ by an interruption rather than by length, which `doc/devel/specs/ssr_interrupted_repeat_recall.md`
addresses for the older caller.

## 4. The parameters file

**`[repeat_tracts_left_out].by_stratum` gains a count, `interrupted`**: how many of the class's
segments §2 left out. It goes into format 9, which is not yet on any pushed main and so has never
been written by a released build; a format-8 or older file has no such count and is read as having
left none out. Nothing at calling reads it, like the rest of the section.

**The file's note says** that interrupted segments are left out of the fit and that calling scores
them at their longest pure run.

## 5. What each end of the range gets

The change depends on the reference and not on the reads, so it acts the same at every depth and
cohort size: an interrupted segment is scored as its longest run whether one sample covers it at 3
reads or 2,000 at 15.

- **One sample, few reads.** The fit has fewer segments to learn from at long homopolymers, where
  segments were already few (30 at 19 repeats in tomato's selection); the curves across classes
  carry those classes, as they carry any thin class (`str_penalised_fit.md` §2.11).
- **Thousands of samples.** The stutter-tail rule now meets classes of one kind of segment, so its
  rounds have no second population to strip away (§7 measures it).
- **Genomes with few interrupted segments** (human at the lengths the catalog keeps): the fit loses
  few segments, and calling changes only at the segments that are interrupted.

## 6. The allele-length prior stays at the segment's own class

**Settled with the owner, 2026-11-06: an interrupted segment keeps the prior of its own class
`(period, R)`.** The prior says how a class's chromosomes spread over lengths around the reference
(`str_fit_to_calling.md` §3). After §2 the 19-repeat class's prior is fitted from pure segments only,
and pure segments may vary in length more than interrupted ones, so an interrupted segment could get a
prior somewhat wider than its own. That is an expectation, not a measurement, and moving the prior to
the pure-run class `(period, P)` would give it the spread of much shorter segments, which need not
suit it better. The stutter difference is measured and large (§1); the prior weighs most at a few
reads a sample and little at depth. So the prior is left alone, and §7's first check measures what
§2 does to it: if the 19- and 20-repeat classes' spread widens materially once interrupted segments
are out, this is reopened with that number.

## 7. How we know it works

1. **The kimura probe again**, with §2 in place: the homopolymer slip rates at 17 to 22 repeats lie
   between their neighbours', with no fall at 19 and 20, and the stutter-tail rule leaves out at
   those lengths about as few segments as at 14 to 16 repeats (0 to 1 in 1,000 there). **And the
   allele-length prior's spread at 19 and 20 repeats**, against today's fit of the mixed classes
   (§6).
2. **Tomato2's fit**: the 19-repeat class's rate moves from 6.3 slipped reads in 100 toward the 16-
   to 18-repeat classes' 18 to 26.
3. **Calls, HG002 chr20 at 30× and 5×**: no loss of true calls at interrupted segments, and the
   counts of false calls and wrong genotypes at STR segments reported against the current main.
4. **Unit tests**: the longest pure run of perfect, interrupted and phase-shifted segments; the
   class a candidate is looked up at, including below one; a perfect segment's lookup unchanged; the
   fit skipping a segment that is not perfect and counting it.
