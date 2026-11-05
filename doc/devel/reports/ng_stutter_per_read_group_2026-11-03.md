# Stutter per read group: measured on the simulation, HG002 and tomato2

*2026-11-03, branch `stutter-per-read-group` (commits `0f3e6b05` to `3d0279fa`). Spec:
[`str_library_stutter.md`](../ng/spec/str_library_stutter.md); the thin-class rule it needed:
[`str_penalised_fit.md`](../ng/spec/str_penalised_fit.md) §2.11.*

## 1. What the branch does

Before this branch every read group of a cohort was scored at one slip rate per class. A **class** is
every repeat tract with one motif length and one reference repeat count. Now each read group gets a
**stutter multiplier** per motif length: its reads slip at the class rate times its multiplier, capped
at 0.95. A multiplier of one is a read group of typical stutter; the multipliers of a motif length
have a geometric mean of one, weighted by each read group's reads. The fit places the class rates and
the multipliers together, the parameters file (format 8) writes them in
`repeat_tracts.stutter_multiplier_by_read_group`, and calling scores each read group at its own rate.

A multiplier with fewer than about 20 slipped reads behind it is held near one by a pull toward one,
and is written without a standard error. At 3 reads a position that is every multiplier of motifs of
three bases or more: a tomato2 accession has about 1,160 slipped reads at homopolymers, 160 at
dinucleotides, 9 at trinucleotides and under 2 beyond.

The branch also brings one change that acts on every run, with one read group or many: **a class of
fewer than 3 loci with reads no longer shapes its motif length's stutter curves**, and past the last
class that does it is held at that class's numbers (§4 shows what this moves on HG002).

## 2. The simulation: 50 inbred samples at 3 reads, each with its own true stutter

The stuttering simulation draws each sample's library with its own multiplier, log-normal with a
standard deviation of 0.4 on the log scale (0.57 to 1.77 between the 10th and 90th sample, 0.43 to 3.19 at the extremes). One
multiplier per sample, the same at every motif length.

**The multipliers recovered**, against the truth rescaled to a geometric mean of one:

| motif length | median miss | within two standard errors of the truth |
|---|---|---|
| 1 | 4.4% | 43 of 50 |
| 2 | 8.6% | 39 of 50 |
| 3 | 13.5% | 35 of the 46 that have an error |
| 4 to 6 | 17% to 25% | no error written: they lean on the pull |

For 50 draws from a correct error about 47 would fall within two; homopolymers are close, and the
shortfall grows with the motif length because the pull draws every multiplier toward one (a fitted
log multiplier is 0.91, 0.82 and 0.60 of the true one at motif lengths 1, 2 and 3).

**The calls.** Three runs on the same psps, scored per sample and tract against the simulated truth
(`scripts/ng_str_stutter_score.py`). A *cell* is one sample at one tract.

| parameters | true variants, genotype right | true variants, genotype wrong | reference tracts called variant |
|---|---|---|---|
| fit with one stutter for every read group (the branch's start, `92b00dd8`) | 41,855 | 991 | 387 |
| fit with each read group's multiplier (this branch) | 41,880 | 929 | 324 |
| the same fit with its stutter replaced by the true one | 41,941 | 920 | 318 |

Out of 46,231 cells at true variants and 103,769 at reference tracts. The two kinds of error together
are 1,378 with one stutter, 1,253 with the multipliers and 1,238 with the true stutter: **the
multipliers remove 125 of the 140 errors the true stutter removes.** The genotypes gained are fewer,
25 of the 86 the true stutter gains; most of the change is fewer wrong calls rather than more right
ones. The gain sits at homopolymers and dinucleotides, where the multipliers are placed (wrong plus
false calls 350 → 296 and 347 → 303, against 279 and 311 with the true stutter); at three to six
bases, where the multipliers stay near one, the three runs make 681, 654 and 648 such errors.

*Not the benchmark's own true-stutter run.* `ng_true_stutter` holds the true stutter with every other
number at the built-in defaults, including the allele-length prior and the inbreeding coefficient of
an inbred cohort; it gets 27,560 genotypes right, so it measures the defaults rather than the stutter.
The third row above is built from the branch's fitted file with only the stutter tables replaced
(`tmp/true_on_fit.py` on the branch).

The one-sample simulation at 30× has one read group and is not rerun: a run of one read group has
every multiplier one.

## 3. tomato2: 59 accessions at about 3 reads, one sequencing project

**The multipliers are narrower than the per-accession figure that motivated the work.** Between the
10th and the 90th accession they run 0.88 to 1.13 at homopolymers and 0.84 to 1.18 at dinucleotides
(extremes 0.70 to 1.30 and 0.59 to 1.42), with median standard errors on the log scale of 0.028 and
0.069. The fivefold spread of §1.1 of the spec (0.37% to 1.98% of a sample's reads off its most common
length) was measured on tomato1, 51 accessions from several experiments; tomato2's accessions all come
from one project (PRJNA353161).

**The motif lengths do not behave alike.** An accession's homopolymer and dinucleotide multipliers
correlate at r = 0.56 on the log scale, and 13 of the 59 differ by more than two of their combined
standard errors, where about 3 would by chance if each accession had one multiplier. This supports
the owner's ruling of one multiplier per motif length.

**Read groups stuttering differently is not why tomato2's tracts vary more than a binomial allows.**
At tracts where every accession's reads favour the same length, off-length reads vary from tract to
tract more than a binomial at the class's rate predicts. Allowing each read group its own rate leaves
that unchanged:

| motif length | tracts | observed over binomial variance, one rate a class | each read group's own rate |
|---|---|---|---|
| 1 | 6,148 | 2.01 | 2.01 |
| 2 | 422 | 2.64 | 2.65 |
| 3 to 6 | 157 | 3.80 | 3.79 |

With multipliers this close to one, a tract read by one mix of accessions expects nearly what another
does, so the mixing cannot produce a twofold excess. What is left is variation between tracts of one
class, which the planned per-locus slip rate is meant to model (`str_penalised_fit.md` §6).

**Time.** The repeat-tract fit took 245 seconds, with two other runs sharing the laptop, against 38
seconds with one stutter for every read group; on an idle laptop at the build step it took 3 min 08 s.

## 4. HG002 chr20, 30× and 5×: one read group

Every HG002 run has one read group, so its multipliers are written `defaulted` at one and calling is
not touched by them. **The thin-class rule does move its fitted stutter**, at the long ends of several
motif lengths, where classes hold one or two loci. Built without the rule, this branch's fit equals
the starting commit's line for line except for the new table, so the rule is the whole of the move.

| depth | what moved more than 10% | largest moves |
|---|---|---|
| 30× | 46 of 153 classes | dinucleotides of 35 to 45 repeats 0.08–0.12 → 0.09–0.26; pentanucleotides of 12 to 20 repeats 0.015–0.020 → 0.018–0.039; tetranucleotides of 23 to 25 repeats 0.15–0.48 → 0.10 |
| 5× | 47 of 146 classes | trinucleotides of 12 to 29 repeats 0.019–0.63 → 0.017–0.044; dinucleotides of 33 to 50 repeats 0.11 → 0.12–0.18; hexanucleotides of 6 to 13 repeats 0.009–0.17 → 0.008–0.027 |

**The calls barely move.** Scored with vcfeval against GIAB on chr20's confident regions:

| depth | true positives | false positives | false negatives |
|---|---|---|---|
| 30×, start / branch | 81,339 / 81,338 | 373 / 374 | 841 / 842 |
| 5×, start / branch | 57,504 / 57,504 | 4,346 / 4,346 | 24,913 / 24,913 |

The rule removes the worst end-of-curve rates, the 5× trinucleotide at 29 repeats at 0.63 and the 30×
tetranucleotide at 25 at 0.48, and also raises some dinucleotide ends at 30× to 0.26, with nothing on
chr20 to say which is right. Few tracts on chr20 are that long, which is why the calls do not move.

## 5. The fit's time at many read groups (2026-11-04)

**On kimura's cohort (2,655 read groups) the multipliers example ran 23 hours without finishing and
was stopped**; without multipliers the same repeat-tract fit took 1 h 05 min there. The cause was
measured on drawn homopolymer cohorts (9 classes of 100 loci at 3 reads, one read group a sample, 6
threads; the ignored test `measure_the_fit_time_against_the_read_groups`):

| read groups | one stutter for every read group | each read group's multiplier |
|---|---|---|
| 60 | 1.5 s | 7.1 s |
| 960 | 4.3 s | 115–126 s |

The time follows the distinct rows a round scores: a row keyed by read group cannot merge with
another read group's, so one class grows from 3,455 rows to 60,502. Building a read model per
multiplier is 0.3 s of the 87 s spent weighing genotypes. Rounding the multipliers to bins merges
too little (2% bins: 25,479 rows; 10%: 12,774).

**What is built: each motif length starts where one stutter for every read group settles**, and the
multipliers then enter and alternate with the class rates until both settle (spec §3.1, amended).
At 960 read groups the fit takes 82 s against 115–126 s; on tomato2 the multipliers example takes
79 s against 108 s. The multipliers agree with the old path to 0.02% on tomato2; slip rates to 0.1%,
except at the long ends of the dinucleotide and trinucleotide curves (classes of 1 to 52 loci), up
to 1.3%, because the curves' stiffness is now chosen during the pooled rounds. On the simulation the
calls move by three errors in 46,231 variant cells (41,877 right, 934 wrong, 322 false calls against
41,880, 929 and 324).

*Tried and dropped*: extrapolating a steady approach to its end (fewer rounds on drawn cohorts, but
on tomato2 the tetranucleotide fit never settled and the homopolymer fit took 763 rounds where it
takes 199); alternating the slip-rate curve and the multipliers inside a round (142 keyed rounds
either way); one read model per class in the standard-error pass (13 s either way).

**Then four more, owner 2026-11-04** (spec `str_penalised_fit.md` §2.6, §2.8;
`parameter_prepass_ssr.md` §4.6.1). On tomato2's repeat-tract fit, each added to the ones before:

| added | fit time | largest move against the fit before them |
|---|---|---|
| the pooled start (committed `6b07b4fd`) | 80 s | — |
| tolerance 10⁻⁵ instead of 10⁻⁷ | 39 s | multipliers 0.0004%, slip rates 0.06% |
| the loop jumps along its own path, guarded by the likelihood | 32 s | dinucleotide multipliers 0.11%, slip rates 0.25% |
| one read model per class in the standard errors' per-tract passes | 23 s | none (every fixture's bytes unchanged) |
| a tail-rule refit starts from the fit before | 20 s | a homopolymer multiplier 0.01%; a slip rate 0.43%, class of 3 loci |

Against the pooled start alone, 37 of tomato2's 106 classes move their slip rate by more than 0.1%
and none by more than 0.43%; a multiplier's error is 3% (homopolymers) to 7% (dinucleotides). On the
drawn cohort of 960 read groups and 1,000 loci a class the fit takes 87 s against 303 s. On the
simulation the repeat-tract fit takes 12 s against 34 s and the calls are those of the pooled start
(41,877 right, 934 wrong, 322 false calls). *Tried and dropped*: the jump without its guard (on
tomato2 the tetranucleotide fit never settled); recomputing alternate classes each round (more
rounds in proportion).

**What is left**: the standard errors' steps in the three stutter numbers, each scoring the
class's keyed rows, are about a third of the drawn cohort's 87 s.

### 5.1 Kimura, 2026-11-05: memory freed on one thread, and a fit that circled

**Kimura's run of the four changes kept one thread busy, not 96.** After 11 hours its main thread
had used 9.6 hours of processor time and each of the 96 pool threads about 1.5. kimura's processor
(Xeon 6527P, 48 cores on two sockets) allows no sampling profiler to a user, so the fit's parts were
timed on the calling thread, on the drawn cohort of 960 read groups and 200 loci a class:

| part of the fit, per-read-group, seconds on the calling thread | 1 thread | 32 threads |
|---|---|---|
| scoring the rows | 56.0 | 7.8 |
| each row's counts | 96.8 | 11.8 |
| merging the row table's blocks (one thread) | 12.3 | 10.2 |
| freeing the scored rows and their read model | 7.9 | 20.9 |
| freeing each row's counts | 4.2 | 3.2 |

**The largest part was freeing memory**: about 200,000 small vectors a class and round, made on the
pool's threads and freed on the calling one, which Linux's allocator does slowly and one at a time,
longer at 32 threads than at one. Four changes, none of which moves a number:

- each read model's rows and every scored row in one array (21 s → 0.3 s at 32 threads);
- a keyed table built once a run of rounds, given each round's multipliers (tables 13.2 s → 2.6 s;
  peak memory 1.14 → 1.26 GB on the drawn cohort);
- the read groups' multiplier solves and the read models side by side.

At 32 threads the drawn cohort's fit went from 78 s to 37 s, and tomato2's on the laptop from 20 s
to 14.9 s, every output byte the same. **Then, owner 2026-11-05, each round's counts summed in fixed
chunks side by side** (the per-row counts' freeing, 4 s → 0.3 s; 49.7 s → 45.5 s at 16 threads); this
sums in another order, and no number on tomato2 moves by more than 10⁻⁸ of itself.

**The jumps of §5 circled on kimura's homopolymers.** From round 226 to round 331, when it was
stopped, the fit's largest change was 6.6·10⁻⁴, 3.1·10⁻⁴ and 2.3·10⁻⁴ in every cycle of three
rounds. The guard compares the reads' log-likelihood, while a round climbs the penalised objective,
so a jump back passes it. The jumps now stop after 30 rounds without a new lowest change
(`str_penalised_fit.md` §2.8; tomato2 16.0 s against 14.7 s). With the first version of that rule
the homopolymer fit settled in 157 rounds, 18 minutes.

**The stutter-tail rule's refits were each a whole fit.** A refit that starts from the fit before
let the between-locus share move again and read it off pairs three times, the first fit's 30, 6 and
2 rounds on tomato2. It now keeps the share read off pairs (owner 2026-11-05;
`parameter_prepass_ssr.md` §4.6.1): tomato2's homopolymer refit takes 18 rounds against 80, moving a
multiplier by at most 0.002% and a slip rate by at most 0.03%. On kimura the refits went from about
15 minutes each to 16, 9.5 and then about 7.

**A separate finding: the stutter-tail rule removes a large share of kimura's long homopolymers.**
It used all 10 of its rounds and left out 925 tracts, 261 of the 624 at 19 repeats and 123 of the
539 at 20; on tomato2 it left out 3. A tract read by 2,176 samples carries enough reads for the
rule's test to flag far smaller departures than at 59. Whether that is right is not measured here.

*Kimura's time for the whole fit*: pending (the run started 15:40, with the machine at a load of 75
to 85 on 48 cores, about 54 of them this run's).

## 6. Not measured here

- **Kimura's cohort** (2,176 accessions from many projects): the first run stopped after 23 hours
  (§5); the run with §5.1's changes had settled homopolymers and was fitting dinucleotides when this
  was written. Its multipliers are the measurement the owner's extrapolation idea needs (spec §12).
- **The PCR-free arm of the simulation**, which the spec lists: not rerun.
- **The call-time cost of the multipliers**: calling looks a read group's rate up as before, with one
  more map lookup; not timed.
