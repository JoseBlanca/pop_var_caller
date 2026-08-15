# What was asked, what I tried, what broke, and what I now believe

*A narrative account of the work of 2026-08-12 and 2026-08-13 on ng's second route to the calling
parameters. **Written to be read start to finish by someone who has read none of the specifications**,
so §1 sets the scene before anything is measured. The four technical reports it summarises are
[`duplicated_locus_probe_2026-08-12.md`](duplicated_locus_probe_2026-08-12.md),
[`joint_str_estimator_2026-08-12.md`](joint_str_estimator_2026-08-12.md),
[`joint_contamination_2026-08-12.md`](joint_contamination_2026-08-12.md) and the earlier
[`joint_fit_estimator_2026-08-12.md`](joint_fit_estimator_2026-08-12.md); each carries the numbers in
full and the commands that produced them.*

---

## 1. The setting: why there are two routes at all

Before a caller can genotype anything it needs to know the things it will assume. How often does a
sequencer misread a base? How heterozygous is each plant? How inbred? How diverse is the population
the plants came from? At a microsatellite, how often does the polymerase slip? And how much of each
sample's DNA came from a different plant altogether?

None of those can be assumed. They differ between species, between sequencing runs, and between
individual accessions, and a caller that guesses them wrongly makes systematic errors rather than
random ones. So ng measures them first, in a pass that reads every alignment once. That pass is what
these documents call **step 4**.

**There are two ways to do it, and this project is building both so it can compare them.**

- **The route that is already built** summarises each sample as it walks — how many positions had one
  read disagreeing with the reference, how many had two, and so on — and then fits that sample's
  parameters from its own summary. It is cheap and it needs no other sample.
- **The route being designed** keeps the raw evidence, but only at a bounded set of positions — about
  two million, **the same positions in every sample** — and then fits every parameter once, across the
  whole cohort at once.

**The second route exists for one reason.** A summary has forgotten which position was which, so it
cannot weight a plant's genotype by *that position's own allele frequency in the cohort*. The joint
route can. Everything else follows from that: it can tell a position where every sample looks
heterozygous (which is an artefact) from one where a quarter do (which is a real variant); it can ask
which alleles the population actually carries, which is what contamination is identified by; and it
can fit one microsatellite's behaviour using every sample that has reads there.

**The state when I started.** The second route had six specification documents, a program that chooses
the positions, and a program that checks the arithmetic of the estimator against truths it makes up
itself. Several of its decisions had been written by analogy with the first route and had **no
measurement behind them at all**. My task was to take those decisions and find out whether they are
right, in this order: a question about duplicated stretches of genome, then two modules of code, then
the microsatellite half of the estimator, then contamination.

---

## 2. Duplicated stretches: does the population exist, and can coverage find it?

### What was unknown

The estimator has two kinds of position. At an **ordinary** one, reads disagree with the reference
only when the sequencer misreads a base — about 2 bases in 1,000. At a **noisy** one they disagree
more often, about 1 in 20, which is what mismapped reads look like.

Neither describes a position where **half** the reads disagree. Yet those exist, and there is a
concrete reason: **if a plant carries two copies of a stretch that the reference genome has only
once, both copies' reads pile up at the same position.** Wherever the two copies differ from each
other, half the reads say one base and half say the other — which is exactly what a heterozygote
looks like. The model's only home for such a position is *heterozygous*, and heterozygosity is one of
the quantities the whole pass exists to estimate. So the artefact inflates the answer.

The proposal was a third kind of position, and — this is the important part — one told apart not by
its own reads but by **the read depth around it**, because two copies collect twice the reads. The
specification could not say whether that population is worth modelling, because **nobody had looked**.
It named the measurement and stopped: *on one tomato sample, what fraction of positions sit in a
window of doubled coverage and also read about half non-reference?*

### What I did

Wrote a program that walks one alignment and, for every 500-base-pair window of the genome, records
the average read depth. Two details mattered.

- **Coverage genuinely varies with GC content**, and by a lot: on tomato the median window depth runs
  from 16.2 reads a base at 20% GC to 29.0 at 36%, a factor of 1.79. That is larger than the doubling
  being looked for, so each window's depth is divided by what its own GC content predicts, using a
  curve fitted from the sample itself.
- **Depth is compared to the sample's own median**, so "twice the coverage" means twice *this*
  sample's normal, not an absolute number.

Then, for every position, whether it reads between 35% and 65% non-reference.

### What came back

**The population is there, and coverage separates it cleanly.** On a 25× tomato accession:

| | share of positions | reading near half non-reference |
|---|---:|---:|
| windows at normal coverage | 86% | **0.033%** |
| windows at about twice coverage | 0.86% | **1.26%** |

That is a 38-fold difference at matched read depth, and the two conditions together land on 1 position
in 8,600 — 24 times what they would if coverage and read composition were unrelated. All eight
accessions I walked, from 2.5× to 28.7×, separate the same way.

**The 1.26% is itself a biological number.** It is the share of positions at which the two duplicate
copies differ from each other — sequence divergence between recent duplicates. Nothing in the caller
sets it, and it is why the whole population is small: a duplication is invisible everywhere its copies
still agree.

### Three things I did not expect

**First, the specifications were wrong about the size, by a factor of thirty.** They put this
population at 1,700 to 8,400 positions per two million and concluded that duplication contributes more
false heterozygotes than sequencing error and mismapping put together. Measured, it is **150 to 590**.
The old figure had read a fitted mixture weight — *how many positions sit in such a window* — as
though it were a count of positions actually showing the signal. Both quantities are real; they differ
by exactly that 1.26% divergence. The conclusion those documents build on it survives, because the two
large terms were always the other two, but the paragraph asserting duplication is the largest is now
corrected.

**Second, and this changes what gets built: a 500 bp window does not work at the depth this project
actually runs at.** A window tells one copy from two only once it has collected enough reads. At 25×
a 500 bp window has plenty. At 3.6× it has almost none, its average depth is scatter, and the
separation disappears entirely:

| accession's mean depth | 2.5× | 3.6× | 5.2× | 9.9× | 13.3× | 25.2× | 28.7× |
|---|---:|---:|---:|---:|---:|---:|---:|
| how much better than chance | 1.6× | **1.3×** | 1.5× | 2.5× | 7.7× | 24.0× | 24.9× |

Widening the window restores it: the 2.5× accession goes from 1.6× at 500 bp to 15× at 10 kb. What
matters is depth times width — about 12,000 aligned bases — and the tomato archive runs near three
reads a site, so the fixed 500 bp grid in the specification would have found nothing in exactly the
samples this was adopted for. **The fix is to store fine windows and let the fit add neighbouring ones
together**, which is exact; the reverse is impossible.

**Third, whose duplication is it — the reference's or the plant's?** This decides whether the coverage
summary has to be kept per sample, which costs 1.6 MB each. I put eight accessions on one window grid.
Of the 84 windows that at least one of them reads at doubled coverage, **40 are read that way by
exactly one accession** and 11 by seven or eight. So both kinds are present and the per-plant kind is
the larger: some windows are the reference's own collapsed repeats, which everyone sees, and more are
copy number that varies between accessions. Keeping it per sample is therefore necessary — a single
shared list would force one accession's duplication onto samples that do not carry it.

---

## 3. Two modules of code, and the one idea in them worth telling you about

The next task was ordinary construction: finish the module that decides which positions are kept, and
write the one that defines what each sample records at them. Both are done and tested.

One design idea in them is worth a geneticist's attention, because it is about a failure that is
otherwise **silent**. Every sample is walked independently, on a different machine, possibly months
apart. If two samples disagree about which positions they kept — a different random seed, a different
region file, a different repeat catalogue — then pooling them is meaningless, and **nothing in the
data would look wrong**. So each sample carries thirteen values describing what it was asked for and
in what units it wrote things down, and the fit refuses to proceed if any two disagree, naming the one
that does. Five of those thirteen say only *in what units*: two samples can have kept exactly the same
positions and still have written down depths on different scales, and every other check passes.

Beside that there is a checksum of the positions actually written, blocked per megabase, so a
disagreement says *chromosome 3, megabase 41* rather than *somewhere in two million*.

---

## 4. Microsatellites: the strongest argument the route has produced

### What was unknown

At a microsatellite the observation is not which base a read showed but **how long a repeat tract it
reported**. The polymerase sometimes adds or drops a whole repeat unit while copying, so a read
reports a tract one unit longer or shorter than the DNA it came from. Three numbers describe that:
how often it happens, whether it more often shortens or lengthens — at tomato dinucleotides shorter
outnumbers longer **4.9 to 1** — and how quickly two-unit slips fall off against one-unit slips.

The built route fits those three numbers for a whole **stratum**: all tracts sharing a motif length
and a reference repeat count. The joint route proposed instead to let **each tract's own distribution
of lengths across the cohort** weight the fit, controlled by one further number saying how
monomorphic tracts tend to be — small meaning most tracts are fixed at one length, large meaning every
tract carries the whole stratum's range.

Two specification sections said this, they had been written in one sitting by analogy with the
ordinary-position half, and **not a line of code existed for either.**

### What came back, and the headline is one number

I built a program that draws a stratum with known slippage numbers, draws each tract's length
frequencies, draws each plant's genotype, draws reads, and then fits — comparing the per-stratum model
against the per-tract one.

**At tomato's three reads a site, the per-stratum model does not lose accuracy so much as lose the
parameters.** Against a true slippage rate of 0.0800:

| | slippage rate | shorter-vs-longer | fall-off |
|---|---:|---:|---:|
| truth | 0.0800 | 0.830 | 0.250 |
| per-stratum — **what is built today** | **0.0233** (−71%) | pinned at 1.000 | collapsed to 0 |
| per-tract joint fit | 0.0803 (+0.3%) | 0.828 | 0.253 |

The direction split pinning at 1.000 means the fit concluded that reads *never* lengthen a tract,
which is false; the fall-off collapsing to zero means it concluded that two-unit slips never happen.
Those are the signatures of a parameter that is not identified, not of one estimated badly. At six
reads a site the gap narrows (0.0569 against 0.0800), so **the shallower the cohort, the more this
route buys** — the right direction for a 63-accession archive sequenced at three reads.

**The proposed control number also works.** The "how monomorphic" number comes back at 0.487 against a
truth of 0.500, and the per-stratum model really is its extreme case: its error runs from −37.1% where
87% of tracts carry one length to −0.9% where none does. So *per tract or per stratum* genuinely
reduces to one fitted number, which is what the specification hoped and had no way to know.

### What failed, and I withdrew it

The specification also said **how** to compute the sum over a tract's possible length frequencies:
enumerate the cases where a tract is fixed at one length, or segregating exactly two.

**That cannot represent a tract carrying three lengths, and the consequence is not subtle.** The fit's
only remaining way to explain three lengths among a tract's reads is to say the reads slipped, so the
fitted slippage rate absorbs it:

| share of tracts carrying three or more lengths | error in the fitted slippage rate |
|---:|---:|
| 0% (two length classes — the enumeration is complete) | **+0.9%** |
| 18% | +23.7% |
| 40% | +63.0% |
| 99.9% | **+722%** |

The first row is the control that makes this a statement about the enumeration rather than about my
coding of it: where there is no third length to miss, it is accurate. And the last row is where the
other specification section needed it to work, because "every tract carries every length" is exactly
the per-stratum limit — so the enumeration destroys the very comparison the design was adopted for.

Withdrawn, and replaced with a fixed 256-point numerical integration whose cost does not grow with the
number of lengths at all. At the ±4 units the records store it is *smaller* than the enumeration it
replaces, and it stays accurate when the truth is drawn from the enumeration's own family — the
reverse is not true.

---

## 5. Contamination: the one where I had the direction backwards

### What was unknown

Contamination here means the fraction of a sample's reads that came from a different plant — a second
seedling in the tube, a neighbouring library on the same run. It is invisible wherever the two plants
carry the same allele, and where they differ it shows as a **small** share of reads carrying the other
allele. That smallness is what tells it from a heterozygote, whose two alleles are balanced.

Nothing about it had been measured. And the test that matters is not that a contaminated sample
returns its own fraction — it is that **a clean panel returns zero**, because a false positive means
telling somebody to repeat a sequencing run.

### The first result, and it reversed something written

The specification worried that on a panel of landraces from several regions, using **one pooled
allele frequency** for the whole panel would make population structure look like contamination: an
accession from a diverged region carries alleles the pooled frequency calls rare, and rare alleles
turning up in a sample is the contamination signature.

I drew exactly that panel — 50 accessions in 4 subpopulations, divergence set by `F_st`, three reads a
site — and it is the other way round.

| divergence | true 1% | true 3% | true 10% |
|---|---:|---:|---:|
| none (`F_st` 0) | 0.0103 | 0.0325 | 0.1065 |
| `F_st` 0.10 | 0.0037 | 0.0209 | 0.0829 |
| **`F_st` 0.20** | **0.0000** | **0.0050** | 0.0584 |

On an unstructured panel the estimate is exact. On a diverged one **a genuinely 3%-contaminated
accession comes back at 0.5%, and a 1% one at exactly zero** — both pass any sensible threshold as
clean. **Structure does not invent contamination; it hides it.**

The mechanism, in words: a pooled frequency on a structured panel says every position is genetically
uncertain — the panel-wide frequency is intermediate where each subpopulation is nearly fixed. A model
that already expects uncertainty explains odd reads by genotype uncertainty rather than by
contamination, so the contamination parameter goes to zero. It has lost its power, not its
calibration.

**I misread my own first table because of this.** The clean-panel run showed the pooled frequency
returning zero for every accession and flagging nobody, and read alone that looks like the best of the
three methods. It was not; it was an estimator with nothing left to say. Only the run with a genuinely
contaminated accession in it distinguished the two, which is why that control now sits in the
specification as the test that catches a broken frequency.

**The direction also matters for how the pipeline should behave.** A false-positive story says watch
for panels that flag too much. The measured failure says the opposite: **watch for a panel that flags
nothing**, because that is what a broken frequency looks like.

### The second result: the parameter is fine, the frequency is the problem

Handed the **correct** per-subpopulation frequencies — which no real fit can have, so this is a
ceiling rather than a method — a clean panel returns **0 of 50** accessions above 1% at every
divergence, and a 3% contaminant comes back at 3.0–3.5%. So nothing about population structure makes
contamination unmeasurable. The entire difficulty is in obtaining the frequency.

### The third result, which is a warning about how

I also tried the obvious thing: split the panel into its four groups and estimate each group's
frequencies from its own twelve members. **That is worse than ignoring structure altogether** — it
adds about +0.015 to every accession's estimate and puts 41 to 47 of 50 clean accessions over a 1%
threshold. The reason is simply that twelve samples is too few to estimate an allele frequency from,
and a *noisy* frequency inflates contamination for the mirror-image reason that a *biased smooth* one
deflates it: contamination is the parameter that absorbs reads which do not fit the genotype the prior
expected, so any error in the frequency manufactures some.

**Two error modes, opposite signs.** That pair is the most useful thing this measurement produced, and
it is why splitting a panel into groups is the wrong instinct.

### The fourth result: how many markers, and the budget sits on the threshold

| usable markers | worst spurious estimate in a clean panel | a truly 3% accession |
|---:|---:|---:|
| 3,400 | **1.85%** | 2.90% |
| 13,700 | **0.86%** | 3.00% |
| 55,000 | 0.53% | 3.20% |
| 219,000 | 0.21% | 3.08% |

**A contaminated accession's estimate is right from 3,400 markers up.** What more markers buy is the
*noise floor on the clean ones* — and that floor, not the estimate, is what a threshold has to clear.
The two-million-position budget yields about 10,000 usable markers, which puts the floor at about 1%,
which is the threshold itself. Rather than raise the budget five-fold, the specification now asks for
the panel's own spread of estimates to be reported beside each accession's, so a sample is judged
against its cohort instead of against a constant. That costs nothing — every sample is fitted anyway.

---

## 6. Your two corrections, and what they changed

**On clustering.** You warned that k-means behaves badly when groups differ greatly in size. Agreed,
and the general form of it was already in my measurement: any scheme that assigns samples to groups
and then estimates within a group divides the data, and the estimate degrades with the smallest group.
The answer is not a better clustering algorithm but **no clustering**. Each accession gets a few
coordinates from a principal-component decomposition of the cohort's genotypes — the scatter plot you
already make — and then at each locus the allele frequency is fitted as a **straight line in those
coordinates, using all fifty accessions**. An accession alone at one end of an axis gets its own
frequency, and what it borrows from the panel is the *slope* — how fast frequency changes along the
axis — not any neighbour's allele counts. There is no threshold anywhere on the plot, so an admixed
accession simply gets an intermediate frequency.

**On what a contaminant actually is.** You pointed out that the risk is the samples **sequenced
together**, not the biological population. That is a different and better question, and it applies to
one specific half of the model: the *contaminant's* genotype prior. The *sample's own* prior is still
about ancestry. Splitting those two apart is the change your remark caused, and it removes the harder
problem from the contaminant's side entirely, because a sequencing batch is a list rather than
something to infer.

**So I checked whether we can know the list, and mostly we cannot.** The tomato accessions' read-group
headers carry no platform unit, and SRA rewrote the read names to `SRR7279481.37559618:TTAGGC:37559618`
— the multiplexing barcode survived, the flowcell and lane did not. The human benchmark kept its read
names intact (`HISEQ1:23:H9UD5ADXX:2:…`, flowcell `H9UD5ADXX`, lane 2) although its header says
`PU:unknown`. So one cohort has it hidden in the reads and the other has lost it, and neither has it
where the file format puts it. Your decision — **user-supplied, defaulting to "everything was
sequenced together", refusing a partial list rather than completing it** — is now in the specification
with that evidence beside it.

It also makes the answer better where the list *is* known: a batch is a few dozen samples, all of
which the fit already holds genotype estimates for, so it can ask **which** batch-mate the stray reads
resemble. A real event favours one donor consistently; a spurious estimate favours none.

---

## 7. What I got wrong along the way

Worth recording, because each one nearly went into a document as a finding.

1. **I read a null result as a success.** The clean-panel table showed the pooled frequency flagging
   nobody, and that is what a good method looks like. It was an estimator with no power. Caught only
   because I ran the positive control, which I had nearly skipped as a formality.
2. **My first version of the "fix" was a strawman.** I implemented per-individual frequencies by
   splitting the panel into groups, which is what `verifyBamID2` explicitly does *not* do — it borrows
   across the whole panel. So my first table said the fix was worse than the problem. Adding a third
   arm using the frequencies the genotypes were actually drawn at is what separated *the frequency is
   wrong* from *the frequency is right but noisily estimated*, and only then did the picture make
   sense.
3. **My first microsatellite program could not be run.** It recomputed every read's likelihood for
   every candidate value of every parameter, which put one fit at several hours. Restructuring it so
   the read likelihoods are computed once per slippage value and reused across the hundreds of
   candidate frequency vectors brought it to about a minute. No result changed; it simply became
   possible to have one.
4. **I trusted a specification's arithmetic over a measurement, twice** — the duplication population's
   size and the direction of the structure effect. Both were stated confidently, both had a citation,
   and both were wrong in a way only a measurement would show. The second one contradicted a number
   quoted correctly two paragraphs above it in the same document.

---

## 8. What I now believe, and what is still unknown

**Believed, with numbers behind it.**

- Duplicated stretches are a real population, coverage finds them, and they are about thirty times
  rarer than the design assumed — real enough to model, not the dominant term.
- The coverage window has to be sized by depth, and the tomato archive needs about ten times the width
  the specification fixed.
- The microsatellite half of this route recovers slippage where the built route cannot, and the gap is
  largest at exactly the shallow depth this project runs at. **This is the strongest case for the
  route so far.**
- Contamination is measurable and population structure hides rather than fakes it; the whole
  difficulty is getting a per-accession allele frequency without partitioning the panel.

**Still unknown, and I would not guess.**

- **Everything above except the duplication measurement is against made-up truths.** That grades the
  arithmetic and the model's own consistency, not whether the model describes tomato. The comparison
  the whole route exists for — both routes on real reads, on the human trio and on the 63 accessions —
  needs the estimator itself, which does not exist yet. That is the next substantial build.
- **Whether the per-accession frequency is worth adopting** is being measured as this is written: the
  early numbers say it recovers a 3% contamination that pooling misses entirely, at the cost of a
  higher noise floor on clean samples, and the run in flight says whether that floor falls with more
  markers and whether accessions in a group of two still get their contamination back.
- **How much structure the 63 tomato accessions actually have** is unmeasured, and it decides which row
  of §5's table they sit in.
