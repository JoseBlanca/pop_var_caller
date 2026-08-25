# ng — the calling loop: turning evidence and parameters into genotypes

*Design spec draft, 2026-08-19. **No code yet — this settles the design.** Third and last of the
three documents on variant calling. The other two define functions: the **genotype prior**
([`calling_priors.md`](calling_priors.md)) says how likely each genotype is before any read is
looked at, and the **read likelihood** ([`read_likelihoods.md`](read_likelihoods.md)) says how
probable a sample's reads are given each genotype. **This document says when those two are called,
with what, and how their answers become a call.***

*Reads on: [`cohort_merge.md`](cohort_merge.md) — the evidence this consumes;
[`run_streaming.md`](run_streaming.md) §3.5 — the run this sits inside. Production's equivalent is
[`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs) (SNP/indel) and
[`ssr/cohort/em.rs`](../../../../src/ssr/cohort/em.rs) (STR). Everything said about those files is a
record of what they do, not a proposal to change them — `src/ssr/` and `src/var_calling/` are frozen
production.*

---

## 1. What this is

**At one locus the caller has to answer two questions at once, and neither can be answered first.**
What genotype does each sample carry? And how common is each allele in this cohort? A sample's
genotype is easier to judge once you know which alleles are common — that is the whole point of the
prior — but the only way to learn which alleles are common is from the samples' genotypes.

**The way out is to guess, then improve.** Start from a guess at the allele frequencies; work out
each sample's genotype probabilities under that guess; add those up to get a better estimate of the
frequencies; repeat until the frequencies stop moving. That loop is this document. It is
expectation-maximization, and the name is worth stating once and then not leaning on: what matters
below is which quantities move on each pass and which do not.

**When the loop stops, one more pass over the samples produces the output**: for each sample, the
most probable genotype, and how confident that is.

**The allele frequencies are not the only thing the loop can move.** At a repeat tract, how often a
read gains or loses a whole repeat — and which way, and by how much — is measured before calling
starts. Three of those numbers are pooled over every tract that shares a motif length and a repeat
count ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4); the fourth, **how often a read
slips at all, is no longer a per-cell number at all** — it is a curve in repeat count, fitted once
per motif period and read off at each cell
([`str_slippage_level_curve.md`](str_slippage_level_curve.md), settled 2026-08-20). **§5.1 and §12's
Q2 were written against the per-cell level and have not absorbed that**, which changes what a
per-locus re-fit is pulled back *toward* and therefore what "how far a locus sits from its stratum"
measures. Q2 carries it. A single tract can behave unlike the
rest of its class: an interruption, a nearby indel, somatic instability. So the loop can also be
built to re-fit those numbers from the locus's own reads and run itself again on them.
**Production's STR caller does this, HipSTR does it differently, and ng as specified here does
neither** — §5.1 lays the three side by side, and §12's Q2 is the measurement that
decides which one ng keeps.

### 1.1 Goals

1. **Call the whole cohort at a locus in one pass of the pipeline**, from the evidence the merge
   hands over and the parameters the pre-pass froze.
2. **The same evidence gives the same calls**, at any thread count and any division of the genome
   into work. The merge already carries this requirement and this step is where it is easiest to
   lose ([`run_streaming.md`](run_streaming.md) §12).
3. **Compute the expensive part once for each set of parameters it depends on** — never once per
   pass. Which of the two costs at a locus is the larger is unmeasured (§12's Q5); what is settled
   is that the read likelihood must not be rebuilt needlessly. With every parameter frozen it is computed once
   per locus and every pass reuses it; if the slippage numbers are re-fitted at the locus, it is
   rebuilt once per re-fit, which is at most four builds at production's round cap against the 50
   a per-pass rebuild would cost (§5, §8).
4. **Degrade across the committed range** — one sample to several thousand, a few reads a position
   to several hundred (`CLAUDE.md`, *What this caller has to work on*). §7 answers both ends.
5. **A locus that will not converge is emitted and flagged, not dropped and not fatal.**

### 1.2 Non-goals, and what this does not do

- **It does not define the prior or the read likelihood.** Both are siblings. This document calls
  them.
- **It does not decide whether a variant is emitted.** That is site filtering and emission, step 11
  of [`ng_proposal.md`](ng_proposal.md), and production's STR work found emission and genotyping
  behave independently. What this produces is genotypes and their confidence; what is done with
  them is a separate document.
- **It does not define the site quality, but its last pass computes one** *(amended 2026-08-25)*.
  [`calling_quality.md`](calling_quality.md) owns both quality numbers. Two of the three things it
  needs exist only inside this loop — the per-sample posterior row, which is one reused buffer, and
  the `samples × genotypes` likelihood table, which is per-worker scratch — so that document's §3
  places their arithmetic in this loop's final pass rather than downstream, and §9 below records
  what comes out. Nothing about the loop's own arithmetic changes; it calls a function and fills
  two more fields.
- **It does not select the candidate alleles from what the merge unified** — that is step 6's,
  and it has its own spec ([`candidate_alleles.md`](candidate_alleles.md)). What this document does
  own is whether the loop may **add** to the set once calling starts (§4.1, §12's Q3).
- **It fits nothing the parameter pre-pass fits.** Every error rate, contamination fraction and
  inbreeding coefficient arrives frozen and leaves unchanged. What the loop may move is §5's, and it
  is two things: this locus's allele frequencies, always; and this locus's slippage numbers, if
  §12's Q2 says they earn their cost.
- **It does not change how wide a locus is, and this is excluded rather than deferred.** A locus's
  extent is settled before calling starts, by the merge chaining together the positions whose
  observations overlap ([`cohort_merge.md`](cohort_merge.md) §4.1). **freebayes does the opposite** —
  it grows a haplotype window until the window stops changing, so what is being called is decided by
  the calling ([`ng_proposal.md`](ng_proposal.md) step 6). ng will not do that. The loop is handed a
  locus and calls it; §4.1's discovery adds *alleles* within that extent and never widens it. It is
  worth saying plainly because "discover alleles from the calling" and "discover the locus from the
  calling" sound alike and only the first is on the table.
- **It does not phase, and it does not link loci.** Each locus is called on its own evidence.

### 1.3 Vocabulary

Three terms, and only the third is this document's.

- **locus** and **cohort observation** — one stretch of genome, with every sample's evidence over it,
  grouped and collated by the merge ([`cohort_merge.md`](cohort_merge.md) §4).
- **candidate genotype** — one assignment of alleles to a sample's copies. At a locus with `A`
  candidate alleles and ploidy `P` there are `C(A + P − 1, P)` of them: 21 at six alleles and a
  diploid, 126 at a tetraploid.
- **expected allele copies** — for one sample, how many copies of each allele it carries *on
  average* over its own genotype probabilities. It is a fractional number, not a call, and using it
  rather than a called genotype is what lets the loop work at low coverage where no genotype is
  certain. The prior consumes exactly this ([`calling_priors.md`](calling_priors.md) §6).

---

## 2. The loop

**In words first.** Every sample's reads are scored against every candidate genotype, once. Then:
guess the allele frequencies; combine them with each sample's read scores to get that sample's
genotype probabilities; add the probabilities up across samples to get better frequencies; repeat.
Stop when the frequencies stop moving. Then score every sample one last time and take its best
genotype.

**In steps**, at one locus, for a cohort of `n` samples:

**Three loops, one inside the next, and ng ships with the outer two switched off** — so by default
the whole of the outer two is one pass through their bodies and only the innermost repeats.

```text
candidate alleles := from candidate selection                            (§4)

repeat                                          the discovery round — §4.1. OFF by default
  slippage numbers := the parameter fit's per-stratum values

  repeat                                        the slippage round — §5.1. OFF by default
    for each sample s, each candidate genotype g:
        Lg_table[s][g] = the genotype likelihood       (read_likelihoods.md §2.1)
                        — built from the per-allele emissions Lr(observation | allele),
                          which is where the cost is (§4, §8).
                          On a discovery round after the first, only the NEW alleles' Lr
                          values are computed; on the SNP/indel path the surviving
                          Lg rows carry over too, on the STR path they do not (§4)

    initialise
        one E-step with NO prior — reads only — then sum: the cohort's expected copies (§3)

    repeat                                      the frequency loop — the innermost
        E-step   for each sample s:
                     for each genotype g:
                         score = Lg_table[s][g] + log prior(g | s)   (calling_priors.md §3)
                     posterior over g := softmax of those scores
                     that sample's expected copies := Σ_g posterior(g) × copies of each allele in g

        M-step   cohort's expected copies := Σ over samples, in a fixed order

        stop when the cohort's expected copies have stopped moving (§6)

    if the slippage re-fit is off: leave the slippage round
    slippage numbers := re-fit from this locus's reads, weighted by the genotype posteriors,
                        pulled back toward the per-stratum values                     (§5.1)
    stop when every re-fitted number has stopped moving, or at the round cap (§6)

  if discovery is off: leave the discovery round
  candidate alleles += tract lengths the converged posteriors say are being explained
                       as slippage but recur too often in one sample to be    (§4.1)
  stop when a round adds nothing, or at the allele cap

finally
    if discovery ran: drop alleles no sample's best genotype used, and re-run
                      the frequency loop on what is left                       (§4.1)
    for each sample: report the genotype with the highest posterior, and its confidence
```

**Every sample gets its own prior, and no two are the same.** Each sample is judged against what the
*other* samples showed: its own expected copies are subtracted from the cohort total before its
prior is built, so its own reads cannot arrive twice — once through the likelihood, and once through
the frequency they helped set ([`calling_priors.md`](calling_priors.md) §6). The subtraction leaves
every sample with a different concentration, so there is no single prior a locus could work out once
and hand round. **The E-step therefore builds `n` priors per pass, where a shared prior would build
one.**

**What that is worth, in the work one pass actually does.** Take a locus with six candidate alleles
and a diploid cohort, so 21 candidate genotypes (§1.3). For one sample, one pass does:

- **the prior over those 21 genotypes** — 42 `lgamma` calls and 6 `logsumexp` calls. The `lgamma`
  count is production's primitive as written: one hoisted per allele, then one more for each allele
  a genotype actually carries, which is one for each of the 6 homozygous genotypes and two for each
  of the 15 heterozygous ones ([`genetics.rs:149`](../../../../src/genetics.rs),
  [`:164`](../../../../src/genetics.rs)). The `logsumexp` calls are the inbreeding mixture's, one per
  homozygous genotype ([`calling_priors.md`](calling_priors.md) §3.2);
- **everything else** — 21 additions to combine each genotype's read likelihood with its prior, 21
  exponentials to normalise them, and 36 multiply-adds to turn the result into expected allele
  copies.

**So the prior is not a small term that the leave-one-out rule multiplies — it is the part of the
E-step that carries the expensive function, and the rule multiplies it by the cohort size.** Take a
thousand samples and five passes, of which the first builds no prior (§3): one locus costs
`42 × 1000 × 4` = **168,000 `lgamma` calls**, where a shared prior would cost 168. **How much
wall-clock that is has not been measured**, here or in production; what is measured is the count.

**Whether it dominates the locus is a different question, and cohort size is not what decides it.**
Both costs are linear in the number of samples — the likelihood table is `samples × genotypes` values
built once from `samples × observations × candidates` emission evaluations, the loop's work is
`passes × samples × genotypes` — so `n` cancels out of the
comparison. What decides it is the allele count, because genotypes grow as `C(A + P − 1, P)` while
candidates grow as `A`, and the pass count, which multiplies only the loop's side. §12's Q5 open
question is where the two cross; neither side has been timed for ng.

**What cohort size does decide is whether the `n` priors are worth building at all.** §11 records the
fast path: past some cohort size the leave-one-out subtraction is too small to change a genotype, and
one prior would serve every sample. The size at which that becomes true is the number nobody has.

### 2.1 freebayes has no per-sample prior, so it needs no leave-one-out subtraction

**freebayes never asks how likely one sample's genotype is.** It asks how likely the *whole cohort's*
genotypes are, together: one genotype for every sample, scored as a single object it calls a
**genotype combination** — *"a combination of genotypes for the population of samples in the
analysis"* ([`Genotype.h:165`](../../../../freebayes/src/Genotype.h)). Its prior is a probability of
that whole assignment. It sums four terms, all on by default
([`Genotype.cpp:1595`](../../../../freebayes/src/Genotype.cpp);
[`Parameters.cpp:438-444`](../../../../freebayes/src/Parameters.cpp) for the defaults), and none of
them belongs to any one sample. Two carry the population genetics:

- **how likely that arrangement of alleles is, given the allele counts the assignment implies** —
  one over the number of ways those counts could be dealt out across the cohort's chromosomes, times
  the number of orderings its heterozygotes allow
  ([`Genotype.cpp:1403`](../../../../freebayes/src/Genotype.cpp));
- **how likely those allele counts are at all** — Ewens' sampling formula, evaluated on the pattern
  the assignment implies: how many alleles appear on exactly one chromosome, how many on exactly
  two, and so on ([`Genotype.cpp:612`](../../../../freebayes/src/Genotype.cpp) builds that pattern,
  [`Ewens.cpp:32`](../../../../freebayes/src/Ewens.cpp) scores it). Its only parameter is the
  mutation rate `θ`, fixed on the command line at 0.01 per base and multiplied by the length of the
  haplotype being called ([`freebayes.cpp:309`](../../../../freebayes/src/freebayes.cpp)). **Nothing
  in the calling estimates it.**

The other two are quality terms rather than population genetics, and this document does not propose
porting them: a binomial penalty per allele for lopsided strand and read-placement counts, and a
multinomial term charging an assignment for allele support that does not match the frequencies it
implies ([`Genotype.cpp:1530-1568`](../../../../freebayes/src/Genotype.cpp)). ng puts that class of
signal in the site filter ([`read_likelihoods.md`](read_likelihoods.md) §3.7), not in the prior.
**They are named here because they are on by default**, so a comparison against freebayes that
leaves them enabled is not comparing the two population-genetic terms alone.

**A sample's own genotype probability is then a sum, not a term.** Every assignment that gives that
sample that genotype contributes its joint posterior; the contributions are added and normalised
([`Marginals.cpp:41`](../../../../freebayes/src/Marginals.cpp)). That is marginalising over what the
other samples might be, rather than conditioning on an estimate of what they are.

**Why this belongs in this document.** *No allele frequency is ever estimated*, so no sample is ever
judged against a number its own reads helped produce, and **there is nothing for a leave-one-out
subtraction to repair.** The double counting §2 subtracts away is not a hazard of cohort calling; it
is the price of the shortcut this loop takes — replacing the other samples' genotypes with a summary
of them, so that each sample can be scored on its own. A coder who reads only the paragraphs above
could reasonably conclude that the subtraction is what makes a cohort prior correct. It is not. It
is what makes *this* cohort prior correct.

**And freebayes pays for its version.** There are `genotypes^n` assignments at a locus, so it cannot
enumerate them: it starts at the assignment the read likelihoods alone prefer and climbs, stopping
when the best assignment repeats or a per-locus iteration cap is reached
([`Genotype.cpp:1102`](../../../../freebayes/src/Genotype.cpp)). **What it finally sums over is a
narrow neighbourhood of the winner, and how narrow depends on a setting a reader will not guess.**
Two paths exist:

- **the default.** `-W/--posterior-integration-limits` is `1,3`
  ([`Parameters.cpp:477`](../../../../freebayes/src/Parameters.cpp)), so the marginalisation runs
  over `bandedGenotypeCombinations` — for each sample only its next-best-ranked genotypes, in
  combinations that move at most a bounded number of samples at once — plus the all-homozygous
  assignments;
- **the exhaustive-local path, which the code itself labels a temporary hack.** With bandwidth and
  band depth both zero it rebuilds every assignment differing from the winner in exactly one
  sample's genotype — `n × (genotypes − 1)` of them
  ([`Genotype.cpp:1197-1201`](../../../../freebayes/src/Genotype.cpp) is the branch,
  [`:903`](../../../../freebayes/src/Genotype.cpp) the one-step rule).

Either way the sum is over a neighbourhood, not the space
([`Marginals.cpp:41`](../../../../freebayes/src/Marginals.cpp) does the summing). Under the default,
most of a sample's alternative genotypes appear in **no** assignment at all, so their posterior is
whatever the normaliser leaves them.

**The trade, in one line: ours replaces the other samples with a summary and must subtract a
sample's own contribution back out; freebayes keeps the other samples and must give up on ever
seeing more than a sliver of what they could be.**

**Both are built, and §12's Q1 is the comparison.** ng's default is the loop of §2, which is linear
in cohort size at both ends of the committed range (§7); whether a whole-cohort search can be made
linear too is Q1's to answer rather than this paragraph's. But the choice has
never been measured against the alternative on this project's data, and "we inherited the
expectation-maximization shape" is not a reason to keep it. §12's Q1 says what to build, in what
order, and what to report.

---

## 3. Starting the loop: the first pass runs on the reads alone

**The prior of [`calling_priors.md`](calling_priors.md) is used on every pass but the first, and
nothing about it is discarded.** The first pass is a special case for a mechanical reason: **the
prior needs a number that does not exist yet.** That document builds each sample's concentration as
`α_seed + (the cohort's expected allele copies − the sample's own)`
([`calling_priors.md`](calling_priors.md) §6), and expected allele copies are what a *previous* pass
produces. On the first pass through a locus there is no previous pass. So the first pass has to be
told what to use instead, and there are only two candidates: the seed concentration on its own — the
prior with its cohort term set to zero — or no prior at all, every genotype equally likely.

**Decision: no prior on the first pass.** Its one job is to turn the reads into a first estimate of
the cohort's expected allele copies, and from the second pass on the prior runs in full, seed
included, and the loop converges under it. **What the first pass settles is where the iteration
starts, not what the model is.**

**Why not the seed: the seed on its own can talk the loop out of a variant that is really there.**
The seed says a locus is almost certainly invariant — the hom-ref genotype carries about `1 − 3θ/2`
of the prior mass and a heterozygote about `θ` ([`calling_priors.md`](calling_priors.md) §4), a pull
of roughly 30 Phred at `θ = 0.001` and 20 Phred at `θ = 0.01`. Apply that on the first pass and, at a
locus where the reads are thin, every sample that carries the variant can be scored hom-ref. Their
expected copies of the alternative allele then come out near zero, so the cohort term for that
allele is near zero on the second pass, so the prior is still just the seed — and so on. **The loop
converges, and it converges to no-variant, having never let the reads speak.** GATK names this in its
own allele-frequency calculation, which starts its frequencies flat and only then switches to the
Dirichlet posterior mean:

> *"first iteration uses flat prior in order to avoid local minimum where the prior + no pseudocounts
> gives such a low effective allele frequency that it overwhelms the genotype likelihood of a real
> variant … we want a chance to get non-zero pseudocounts before using a prior that's biased against
> a variant"*
> ([`AlleleFrequencyCalculator.java:188`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/afcalc/AlleleFrequencyCalculator.java))

Production does the same and marks it as a step of its own
([`posterior_engine.rs:2586`](../../../../src/var_calling/posterior_engine.rs) — *"a flat genotype
prior (likelihood only) so the cohort's reads set an honest initial frequency before the
leave-one-out prior engages"*).

**Inherited from both, and never measured here — soft.** The reasoning above is sound arithmetic
about the prior's size, not a count of calls that moved. It should bite hardest where the read
likelihood is weakest against a 20-to-30 Phred prior, which is the tomato panel's corner at 3 reads
a position, and hardly at all at 300. **Q7 is the measurement.**

**It runs at the start of every outer round, not only the first, and that is a decision rather than
an accident of indentation.** §2's pseudocode puts the prior-free pass inside the slippage round and
inside the discovery round, so a second round throws away the expected copies the previous one
converged to and re-derives them from the reads. The reason above — *the number does not exist yet* —
is true only of the very first round; the reason for the later ones is different: **a round starts
with new parameters, and the expected copies it inherits were converged under the old ones.** Carrying
them in would seed the new round at the old round's answer, which is the same self-reinforcing start
the flat pass exists to prevent, with the seed replaced by the previous round. Production restarts
its frequency loop from the seed on every re-fit round for the same reason
([`em.rs:592-602`](../../../../src/ssr/cohort/em.rs)).

**At one sample the first pass has nothing to do**, because the cohort term is zero at every pass:
the prior is the seed on pass 1 whether the first pass is flat or not, so a flat pass 1 followed by a
seeded pass 2 reaches exactly what a single seeded pass would have reached. It costs one pass and
changes no genotype. **Whether to spend a branch on skipping it is Q6, not a rule taken here** —
§7 explains why no branch is needed for correctness.

---

## 4. The candidate set: fixed by default, grown as a measured option

**What ships: the alleles a locus is called over are settled before the first pass and do not change
during it. That is a default, not a settled decision** — the loop is built so it can grow them, and
§12's Q3 is the measurement that will keep the default or overturn it.

The version to build is HipSTR's: after the loop has converged, look for tract
lengths the model is currently explaining as slippage but which recur too often in one sample to be
slippage, add them, and run the loop again — **between whole runs of it, never between two of its
passes** (§2's outermost `repeat`; the table below is why). It is built because it is the one
discovery mechanism that plausibly finds alleles our selection cannot (§4.1).

**Why it is off by default, and it is a cost-against-benefit judgement rather than a mechanical
objection.**

- **The cost is certain and is paid at every locus.** A discovery round is a whole extra run of the
  loop to convergence, because the alleles it looks for are the ones the *converged* posteriors say
  are being explained as slippage (§4.1). A locus where nothing is found still pays a second run to
  establish that. Expected one or two rounds, so roughly a doubling of the loop's own arithmetic —
  not of the emission table, which is appended to rather than rebuilt.
- **The benefit is unmeasured and expected at a minority of loci.** Nothing on this project's data
  says how often an allele is actually hidden under stutter. That figure is the first thing §12's Q3
  reports, and until it exists the option cannot be defaulted on honestly.
- **And it can make calls worse where this caller is weakest.** An allele admitted in error is not
  free: it takes prior mass from every other genotype at that locus, and it gives a homozygous
  sample a second allele to be heterozygous for. HipSTR's bar is a conjunction — at least 2 reads
  **and** at least 15% of one sample's spanning reads — and at low depth **only the count binds**: at
  3 reads a position, 2 reads is already 67%, so the fraction is satisfied automatically and the whole
  bar is "2 of this sample's 3 reads". Two reads is what a single stutter product looks like. So the
  corner of the committed range where the evidence is thinnest is where the weaker half of the bar is
  the only half doing anything (§4.1).

**A separate question is where a discovery round may sit, and there the objection is mechanical.**
Three things react differently to a candidate being added, and a coder needs the map before assuming
either that growth is cheap or that it throws work away:

| what | survives a candidate being added? |
|---|---|
| the emission, `Lr(observation \| allele)` — the expensive part, computed once per (sample, observation, candidate) | **yes.** It depends on the observation, the allele and the frozen parameters, and on nothing else. Growing the set **appends** one column; every value already computed stays right |
| the genotype likelihood `Lg` the loop consumes, per (sample, genotype) | **on the SNP/indel path, yes** — a genotype carrying none of the new allele has an unchanged likelihood, so the existing rows survive and new rows are added. **On the STR path, no**: the term for reads no allele explains is spread over the tract lengths the stutter model can reach *from the candidate set* ([`read_likelihoods.md`](read_likelihoods.md) §4.5), so a wider set changes that spread and every genotype's likelihood with it |
| the prior, per (sample, genotype) | **no, on either path.** The alternative concentration is shared out across however many alternative alleles the locus carries ([`calling_priors.md`](calling_priors.md) §4), so adding one changes `α_alt` for all of them, and with it every genotype that carries any alternative |

**Read that table as the placement rule, not as the decision.** The expensive table survives and the
cheap ones do not, so nothing in it argues against discovery — it argues that a candidate may not be
added *between two passes of the frequency loop*, because the third row plus the convergence test
would leave the loop comparing a different quantity from one pass to the next. Add alleles between
whole runs of the loop, as §2's outermost `repeat` does, and every row of the table is satisfied.

**This section recorded that no document said how the candidate alleles are chosen. One does
now** ([`candidate_alleles.md`](candidate_alleles.md), with its architecture and a shipped
`select_generic`; the repeat-tract half is
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) and is written but not built). The gap was
found here — [`cohort_merge.md`](cohort_merge.md) §13 passed *"choosing candidate alleles from the
table"* to "the calling steps' spec", all three calling documents passed it on to
[`ng_proposal.md`](ng_proposal.md) step 6, which is a row in a table comparing five callers, and
nobody had written the design. **What follows is kept because it is what the design had to
answer**, not because it is still open.

**One consequence was blocking a formula, and that is how the gap was noticed.** The SNP/indel
emission needs `q_sum_other`, the pooled error mass of reads matching no candidate
([`read_likelihoods.md`](read_likelihoods.md) §3.3), and **nothing upstream produces it**: the merge
unifies *every* sequence a sample showed into the allele table ([`cohort_merge.md`](cohort_merge.md)
§4.2), so there is no leftover until something narrows the table. **Selection is what creates that
pool, so whoever specifies selection owes the pool** — the alleles it dropped, with their summed
per-read error mass, not merely a count of them. **It does**: `LocusSelection::unmatched` carries
the reads and the mass per covering sample, summed from the merge's own rows rather than
re-derived from a count and a rate.

**What is true of both paths today: neither invents a sequence.** A coder coming from HipSTR will
expect otherwise. Production's STR candidate assembly admits a repeat length only where some
sample's reads reached it (`occupied`,
[`candidate_set.rs:221`](../../../../src/ssr/cohort/candidate_set.rs)), including the ±1 neighbours
it nominates for a sample that under-resolves. What the two paths really differ in is the *shape* of
the selection, and that is what sets how large `A` gets:

- **the generic path selects flat.** The merge has already unified every sample's observations into
  one allele table ([`cohort_merge.md`](cohort_merge.md) §4.2), so what is left is a cap and a
  support bar over that table.
- **the STR path selects with structure.** Sequences are pooled into **rungs** keyed by repeat
  count, each sample's peaks nominate which rungs to promote, and same-length spellings clear a
  recurrence bar to be promoted alongside the rung's representative
  ([`candidate_set.rs:194`](../../../../src/ssr/cohort/candidate_set.rs)). **That structure is not
  incidental** — [`calling_priors.md`](calling_priors.md) §5.1 lays the prior's mass over the rung
  ladder, and its Q3 asks how a rung's weight divides between two alleles of the same length. The
  ladder and the STR prior are one design.

### 4.1 How HipSTR builds its candidate set — the mechanism Q3 would switch on

**Production picks its alleles once and lives with them; HipSTR picks a set, genotypes, looks at what
the genotyping implied, and picks again.** This section sets HipSTR's version out in full, because it
is what ng builds behind the switch and a coder needs the mechanism, not a summary of it.

| | production's STR path | HipSTR |
|---|---|---|
| when the set is fixed | **before genotyping, once** ([`candidate_set.rs:194`](../../../../src/ssr/cohort/candidate_set.rs)) | **grown during genotyping, to a fixpoint, then pruned** ([`seq_stutter_genotyper.cpp:568`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)) |
| what a candidate is | the tract sequence | the tract **plus five reference bases either side**, so an indel just outside the tract travels with the allele instead of being charged as noise (`LEFT_PAD`/`RIGHT_PAD`, [`HaplotypeGenerator.h:59`](../../../../HipSTR/src/SeqAlignment/HaplotypeGenerator.h)) |
| how a sequence gets in | its repeat count is a clear peak in some sample — or a ±1 neighbour of one, where some sample reached that length — and it clears a bar within its length: minimum reads, minimum distinct samples, minimum share of that length's reads | any of three routes ([`HaplotypeGenerator.cpp:157`](../../../../HipSTR/src/SeqAlignment/HaplotypeGenerator.cpp)): **strong in one sample** — at least 2 reads and at least 20% of that sample's reads; or **common overall** — above 5% of samples or 5% of reads; or **named in a reference panel VCF**, which is the only way in with no read support at all |
| **alleles hiding under stutter** | **never sought.** A repeat count that only ever appears as a slip product of a called allele cannot become a candidate | **sought explicitly, and this is the substantive difference.** After aligning and scoring, each read's alignment is retraced; where the trace says *this read slipped*, the tract sequence it implies is counted per sample. Any such sequence not already a candidate, with at least 2 reads and at least 15% of that sample's spanning reads, is admitted ([`seq_stutter_genotyper.cpp:863`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)) |
| what happens after an allele is added | — | realign, rescore, **and look again** — the loop repeats until a round finds nothing new |
| alleles that earned nothing | kept | **pruned twice**, with the posteriors recomputed after each: those no sample's best genotype used, then those no read spans ([`:648`](../../../../HipSTR/src/seq_stutter_genotyper.cpp), [`:657`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)) |
| runaway | too many alleles → the locus is a no-call | too many haplotypes → **the locus is abandoned** |
| with a reference panel supplied | not a mode production has | **all of the above is switched off** and the panel's alleles are the set ([`:640`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)) |

**The one sentence that separates them: HipSTR asks whether something it is calling stutter is
actually an allele, and production never asks.** That is a real question at a repeat tract, because a
short allele in a sample that also carries a long one is exactly what a contraction slip looks like.

**And HipSTR is the existence proof for §4's table above.** When it adds an allele it does not rebuild
its read-to-allele scores: it maps every surviving haplotype to its new index, **copies its scores
across unchanged**, marks only the genuinely new haplotypes for realignment, and scores those
([`:361`](../../../../HipSTR/src/seq_stutter_genotyper.cpp),
[`:380`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)). Growing the set costs the new columns and
the recomputed posteriors, and nothing else. **So growing the candidate set does not throw the
expensive work away**, and any argument for the fixed default that rests on rebuilt tables is wrong.
§4's default rests on the extra convergence each round costs and on the risk of minting an allele.

#### Where ng's version has to differ from HipSTR's

**HipSTR does run an expectation-maximization loop with allele frequencies in it, but it is not the
loop that produces the genotypes.** Per locus it does three things in order:

1. **A length-based EM.** The reads are reduced to tract lengths in base pairs, and one loop
   re-estimates the allele frequencies **and** all six stutter parameters together on every pass,
   stopping when the log-likelihood stops rising or every parameter moves by less than `1e-4`
   ([`em_stutter_genotyper.cpp:170`](../../../../HipSTR/src/em_stutter_genotyper.cpp)). **Only the
   stutter model survives it. The allele frequencies it fitted are thrown away**
   ([`genotyper_bam_processor.cpp:148`](../../../../HipSTR/src/genotyper_bam_processor.cpp)).
2. **Sequence-level genotyping, with that stutter model frozen** and a genotype prior of two
   constants — one for a homozygote, one for a heterozygote
   ([`genotyper.cpp:34`](../../../../HipSTR/src/genotyper.cpp)). This is where the discovery loop of
   §4.1 sits, so what a discovery round wraps is **one scoring pass**.
3. **Optionally, all of it again.** With `--recalc-stutter-model` the maximum-likelihood alignments
   are retraced, the length EM is re-run from them, and the whole of step 2 repeats
   ([`seq_stutter_genotyper.cpp:1466`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)).

**So HipSTR is not avoiding the cost of an EM — it is putting the EM somewhere else.** It fits the
stutter parameters by expectation-maximization over lengths, and then calls genotypes without one,
against a prior that knows nothing about this locus. **ng merges the two:** one loop fits this
locus's allele frequencies, feeds them to the prior, produces the calls, and optionally re-fits the
slippage numbers around all of it (§5.1). That is why a discovery round costs a full convergence here
and a single scoring pass there.

**Decision: look for hidden alleles against the converged posteriors, not against a first scoring
pass.** The cheaper alternative — look once, early, as HipSTR effectively does — would ask the
question at the moment the answer is least informed, and a hidden allele is exactly the case where
the cohort's evidence decides it: one sample's two extra reads at a shorter length are slippage until
three other samples show the same length. That is what the loop supplies and a single pass does not.
Rounds are expected to be one or two, since a round that finds nothing ends the loop, so the expected
cost is about twice the loop rather than many times it — but **that expectation is a guess and Q3
measures it**, because a locus that keeps finding one more allele is the shape that would hurt.

**A middle setting exists and Q3 should measure it too, because it may get most of this for a
fraction of the cost.** Converge once; **hold the frequencies fixed at what that converged to** and
run the discovery rounds against them, which makes each round a scoring pass rather than a
convergence — HipSTR's cost with a prior that has actually seen the locus; then converge once more at
the end on the final allele set. It gives up only the frequencies' response to an allele *while it is
being decided on*, and that is what the final convergence puts back. **Three settings for Q3, then:
off; discover against converged-and-frozen frequencies; discover against a full convergence each
round.**

**Two things ng inherits from HipSTR unchanged, and one it must not.** Inherit the recurrence bar's
*shape* — a candidate must clear both a read count and a share of one sample's spanning reads, so a
single stray read cannot mint an allele — and inherit the prune, because an allele that no sample's
best genotype used has cost every other sample prior mass for nothing. **Do not inherit the
thresholds as measured constants:** HipSTR's 2 reads and 15% were set on high-depth human data.
**The two halves bind at opposite ends of the depth range:** below about 13 reads the count is the
only constraint, since 2 reads already clears 15%; above it the fraction takes over and makes the bar
stricter. So a single pair of numbers cannot be right at 3 reads a position and at 300. **They are
inherited and soft, and Q3 sweeps them.**

---

## 5. What moves during the loop, and what does not

**Two quantities can move, and they move on two different clocks.** This locus's allele frequencies
move on every pass — that is the loop of §2. This locus's slippage numbers move, at most, between
whole runs of that loop, because changing them invalidates every read likelihood. Everything else is
frozen before the first pass and stays frozen.

| | moves? | where it comes from |
|---|---|---|
| cohort's expected allele copies | **every pass** — this is the loop | the M-step |
| each sample's prior | every pass, but only through the line above | [`calling_priors.md`](calling_priors.md) §3, §6 |
| this locus's slippage level, direction split and fall-off | **between runs of the loop, and only where per-locus re-fitting is switched on** — ng ships with it off (§5.1) | re-fitted from this locus's own reads, pulled back toward the frozen per-stratum value |
| each sample's read likelihoods | **no** — rebuilt only when the line above changes | [`read_likelihoods.md`](read_likelihoods.md), computed once per set of slippage numbers |
| per-read-group error rate, contamination, STR substitution rate | **no** — frozen by the parameter fit | [`read_likelihoods.md`](read_likelihoods.md) §6.1 |
| each sample's inbreeding coefficient | **no** — frozen by the parameter fit | [`calling_priors.md`](calling_priors.md) §7 |
| the candidate alleles | **no while the loop runs** — a discovery round may add to them between whole runs of it, and ng ships with that off (§4.1) | candidate selection, plus discovery where switched on |

### 5.0 A sample the candidate step could not call, and what the M-step does with it

**This whole section is the SNP/indel path's, and the repeat-tract path sets no sample aside**
(owner, 2026-08-24). §5.0.1 is why the two differ; everything before it describes the generic path.

**Candidate selection can hand this loop a sample it has already declared uncallable at the locus.**
When the allele cap cuts a sequence that a sample's own reads had earned, that sample carries
something the locus is no longer called over, and
[`candidate_alleles.md`](candidate_alleles.md) §4.1 fixes its genotype as **missing** rather than
letting the caller invent one. The flag travels on the evidence as
`GenericSampleEvidence::genotype_must_be_missing`
([`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.1); this loop does not re-derive
it and could not, because the pooled error mass that would be the only trace is identical under
every genotype and cancels.

**The E-step is unaffected: such a sample is simply not scored, and emission writes its genotype as
missing.**

**Decision: it leaves the loop entirely, before the first pass** (owner, 2026-08-24). An uncallable
sample is set aside when the loop is entered and takes no part in either step: it is not scored, it
contributes nothing to the M-step's sums, and emission writes its genotype as missing because it
has no call rather than because a call was withheld. **The locus is still called, on the samples
that can be.**

**Why it cannot merely be left in with an uncertain posterior.** The M-step sums every sample's
expected allele copies into the cohort's, and that sum is what the next pass's prior is built from
(§2). A sample whose true allele is absent from the table does not have an *uncertain* posterior —
it has one over the **wrong set**, and it puts its mass on whichever surviving genotype its reads
mismatch least, which is usually the homozygous reference. Including it pulls the locus's allele
frequencies toward the reference by exactly the samples carrying the rarest alleles: the error is
systematic and in one direction, not noise that averages out.

**And this decision is what makes candidate selection's truncation defensible at all**, which is
why it belongs to that document as much as to this one.
[`candidate_alleles.md`](candidate_alleles.md) §4.1 prefers cutting the allele list over refusing
the locus, on the ground that refusing costs every sample a locus most of them were callable at.
That preference holds **only** if the samples that lost an earned allele stop contributing to the
locus's numbers — otherwise the locus is corrupted cohort-wide and refusing it would be the honest
answer. Setting them aside at the loop's entrance is what closes that.

**Excluding costs the loop nothing it does not already tolerate.** A sample with no reads at the
locus already contributes zero for every genotype and is decided by the prior alone (§7), so a
smaller cohort is a shape the loop is built for. What it does change is the denominator of the
fitted frequencies, which is the right denominator: the samples the locus was actually called on.

**The ruling has a producer and no carrier, and the gap is in the shared vocabulary rather than
in this loop** (found 2026-08-24, checking candidate selection against the modules it feeds).
Selection already decides the fact — `UnmatchedSupport::genotype_must_be_missing()` is true for a
sample whose own reads earned an allele the cap then cut. What does not exist is anywhere to put
the answer: `SampleGenotypeCall` is `{ genotype, genotype_quality }` with no absent variant, and
`Genotype::new` panics on an empty multiset, deliberately — *"an empty multiset is not a haploid
call, it is a sample with no genome"*. **So a missing genotype cannot currently be expressed**, and
this section's decision cannot be implemented until it can. The type is `calling/mod.rs`'s; whoever
builds the loop adds the variant there, and §9's hand-off is written against its existence.

**And the two per-sample lists are in different orders.** `LocusSelection::unmatched` is parallel
to the merge's covering samples, `LocusInference::per_sample` is one entry per sample of the run in
run order. A sample that does not cover the locus and a covering sample that lost nothing are
different facts that look identical if the join is by position rather than through the sample index
the merge records.

*What would still be worth measuring, though nothing turns on it:* the tomato panel at the loci
where the cap binds, comparing the fitted frequencies with the samples in and out — it would say
how large the bias would have been. Not measurable until selection is wired into the builder, which
is [`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md)'s work.

#### 5.0.1 Why the repeat-tract path does not need this

**At a repeat tract the loop can put back what selection cut, so no sample is set aside** (owner's
decision, 2026-08-24). A discovery round between whole runs of the loop looks at what the converged
posteriors are explaining as slippage, nominates the tract lengths that recur too often in one
sample to be slippage, and **adds them to the candidate set** (§4, §4.1). A length one sample's
reads earned and the cap removed is therefore not gone for the rest of the locus's calling, which is
the fact the generic path does not have: there the candidate set is settled before the first pass
and stays settled to the end (§4, §5's table). So the repeat-tract path scores every covering sample
on every pass, and `genotype_must_be_missing` is a flag the generic evidence carries and the
repeat-tract evidence does not.

**What separates the two paths is what happens to the cut allele's reads, not the cap itself.** On
the SNP/indel path a read whose sequence is no longer a candidate contributes only the pooled error
mass `q_sum_other` ([`read_likelihoods.md`](read_likelihoods.md) §3.3) — the same number under every
genotype, so it cancels, and the sample's posterior comes out confident over a set that cannot
represent what it carries with nothing in the arithmetic saying so. That is what the paragraphs
above answer. At a tract, a length off the candidate set is still **reached by the stutter model
from a candidate** ([`read_likelihoods.md`](read_likelihoods.md) §4, and §4.5 for what happens
beyond the slip cutoff), so those reads have a likelihood that differs between genotypes rather than
one that cancels.

**One condition, and it is stated because §4 makes it necessary.** The discovery round this rests on
**ships off**, and §12's Q3 is the measurement that would default it on. With it off, a repeat
tract's candidate set is as fixed during the loop as the generic path's, and the whole weight falls
on that path's selection admitting enough rungs to begin with — which is why its cap is 32 against
the generic path's six ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §12). **That
document's §12 argues for 32 partly from this section's missing-genotype rule and so needs
re-reading against this ruling**; the ruling is not in doubt, the sentence resting on it is.

**Why the slippage numbers get a clock of their own rather than joining the frequencies on the
inner one.** Building the read-likelihood table takes `candidates × Σ_s (observations in sample s)`
evaluations of a function whose per-entry cost against the loop's own per-genotype work is
unmeasured (§8, §12's Q5) — on the STR path one evaluation is a stutter term and a substitution
alignment, against the loop's two `lgamma` calls. The slippage numbers are inside that function. Re-fitting them on every pass would rebuild the table up to 50 times per
locus; re-fitting them between runs of the loop rebuilds it at most four times at production's cap
(§5.1). Same adaptation, an order of magnitude apart in cost.

### 5.1 Re-fitting the slippage numbers at the locus — three settings, one to be measured

**What is at stake.** How often a read gains or loses a whole repeat is measured before calling
starts, pooled over every tract sharing a motif length and a repeat count
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4). The case for letting one tract depart
from its class is that a tract can behave unlike it — an interruption, a nearby indel, somatic
instability ([`read_likelihoods.md`](read_likelihoods.md) §6.1). The case against is that the reads
a tract's numbers would be re-fitted from are the very reads being genotyped, so the numbers can end
up describing this locus's noise rather than its chemistry. **Nobody has measured which effect is
larger**, which is why this document specifies the frozen behaviour, builds the machinery for the
other two, and hands the choice to a measurement (§12, Q2).

**The two reference implementations differ in five ways, not one**, and the differences are what the
measurement has to hold apart.

| | ng as specified | production's STR caller | HipSTR |
|---|---|---|---|
| where the numbers start | the per-stratum fit | the per-stratum fit | six fixed constants, identical at every locus and on every dataset — there is no genome-wide fit at all ([`em_stutter_genotyper.cpp`](../../../../HipSTR/src/em_stutter_genotyper.cpp), `init_stutter_model`) |
| how far the locus may pull them | not at all | part way — 50 pseudo-counts on the shape and 20 slipped reads on the level, so a tract with no slips collapses back to its stratum's value ([`stutter.rs:184`](../../../../src/ssr/cohort/stutter.rs)) | all the way — its pseudo-counts are one count and `log(1.1)`, there only to keep the geometric proper |
| how many numbers move | none | three: direction split, fall-off, level ([`em.rs:572`](../../../../src/ssr/cohort/em.rs)) | six: those three for whole-repeat changes, and three more for part-repeat ones |
| what the re-fit reads | — | each read attributed to the **called** genotype's nearest allele; a no-call contributes nothing ([`em.rs:1192`](../../../../src/ssr/cohort/em.rs)) | every read weighted by the **whole genotype posterior**, no genotype called |
| loop shape | one loop | **nested** — an outer round re-fits, rebuilds the likelihood table, and runs the frequency loop to convergence; at most 3 rounds, settable to 0 ([`em.rs:140`](../../../../src/ssr/cohort/em.rs), reached through [`driver.rs:319`](../../../../src/ssr/cohort/driver.rs)) | **flat, and upstream of genotyping** — frequencies and slippage numbers re-estimated in the same M-step of one loop over read lengths, every pass, until the likelihood stops rising or all six parameters move by less than `1e-4`. Only the slippage numbers leave it; the frequencies it fitted are discarded and the genotyper uses a two-constant prior (§4.1) |

**Decision: ng builds the nested shape.** It is what makes the three settings one code path rather
than three: frozen is that code at zero rounds, and HipSTR's setting is that code at zero pull-back.
The flat shape would make frozen a second code path *and* rebuild the likelihood table on every
pass, which is the cost the table exists to avoid.

**Decision: the re-fit reads the genotype posteriors, not the called genotypes.** Everything else in
this loop is built on expected allele copies rather than calls, because that is what lets it work
where no genotype is certain (§1.3), and the re-fit should read the same quantity. This is HipSTR's
choice and not production's, and the alternative would push a thinly-covered tract's numbers toward
whichever genotype won by a whisker and drop its no-call samples entirely. **It is not free, and
production's comment says why it did not build it:** attributing every read under the full genotype
posterior instead of to one called genotype multiplies the attribution work by the genotype count,
unless a per-allele attribution is cached — *"the soft per-read responsibility split is the deferred
refinement"* ([`em.rs:18`](../../../../src/ssr/cohort/em.rs)). The size is unmeasured and Q2
should report it alongside the rest.

**Two things the measurement must not bundle together.** *How many numbers move* and *how hard they
are pulled back* are independent, so a six-number unshrunk arm run against a three-number shrunk one
cannot say which difference did the work. **Vary the pull-back; hold the count at three.** The three
extra numbers HipSTR fits describe part-repeat changes, and ng holds those at production's
placeholders — a part-repeat share fixed at 5% with the two one-step shares tied
([`read_likelihoods.md`](read_likelihoods.md) §4.2) — which have no owner yet
([`read_likelihoods.md`](read_likelihoods.md) §10). A six-number arm would be re-fitting parameters
ng has not specified. **Six becomes askable once the part-repeat estimator has an owner**, and
[`read_likelihoods.md`](read_likelihoods.md) §6.1 records the same ordering from its own side.

**Nothing here applies to the SNP/indel path**, and it is worth saying why rather than leaving the
asymmetry to be noticed. That path has no slippage numbers, its per-read-group error rate is frozen
because re-fitting it from the loci the merge kept would measure the merge's own selection, and its
contamination fraction is frozen because how contaminated a sample is is a property of the sample
and not of a locus ([`read_likelihoods.md`](read_likelihoods.md) §6.1). On that path exactly one
quantity moves.

---

## 6. Stopping

**Converge on the quantity the loop actually feeds back — the cohort's expected allele copies —
divided by the number of chromosomes in the cohort.** Stop when the largest such change between two
passes falls below a threshold. **The division is not decoration and dropping it is the easy
mistake:** expected copies are a count, and the threshold is a fraction. `1e-3` in raw copies means a
different thing at 1 sample from what it means at 1,000, so a criterion written on counts tightens by
the cohort size across the range this caller commits to.

**Why that quantity and not the frequencies the loop reports.** It is what the M-step produces and
what the next E-step consumes, so its movement is what says the loop has not settled. Production
tests exactly this — the change in expected counts, scaled by `1/(ploidy × n_samples)`
([`posterior_engine.rs:2718-2722`](../../../../src/var_calling/posterior_engine.rs)) — and its own
comment gives both halves of the reason: the reported frequency estimate is *"a pseudocount-scaled
readout that does NOT feed back"*, so testing it let a larger pseudocount damp the delta and stop the
loop early; and dividing by the chromosome total *"keeps `q` on the same `[0, 1]` frequency scale …
so `convergence_threshold` and its validation range carry over unchanged in meaning."*

**The two constants are inherited and are marked soft.** Production uses a threshold of `1e-3` and a
cap of 50 passes ([`posterior_engine.rs:86`](../../../../src/var_calling/posterior_engine.rs),
[`:96`](../../../../src/var_calling/posterior_engine.rs)), and its own comment records that the EM
converges in 3 to 5 passes on the GATK reference data, so the cap is ten times the observed need
rather than a tuned value. **ng inherits both, as named constants, and neither has been measured
here.** What would set them is §12's Q4.

**Where the slippage numbers are re-fitted, the outer rounds need their own stopping rule**, and it
is not this one: they stop when every re-fitted number moves less than a threshold of its own, or
when the round cap is reached. Production stops on both, with a cap of three
([`em.rs:572`](../../../../src/ssr/cohort/em.rs)); HipSTR stops when the whole model's likelihood
stops rising or every parameter moves less than `1e-4`. **ng inherits production's rule**, because
its numbers are the ones being compared between rounds and a likelihood test would need the table
rebuilt to be read. Both thresholds are inherited and neither has been measured here.

**A locus that hits the cap is emitted with a flag, not dropped and not fatal.** Production learned
this the explicit way: its error variant for non-convergence is retired in favour of emitting with
`converged = false`, *"so a single hard site doesn't kill the whole cohort run"*
([`posterior_engine.rs:26`](../../../../src/var_calling/posterior_engine.rs)). ng carries the same
rule and the same reasoning. **The flag must reach the output**, because a genotype from a loop that
did not settle is a different claim from one that did, and nothing downstream can tell them apart
otherwise.

---

## 7. One sample, and a thousand

**At one sample the frequency loop is one pass and that is the correct answer, not a degraded one.**
The prior's cohort term subtracts a sample's own expected copies from the cohort total, and with one
sample those are the same number, so the term is exactly zero and the prior never moves
([`calling_priors.md`](calling_priors.md) §6). Nothing the M-step produces can change the next
E-step. **So the caller does one E-step against the frozen prior and reports.** §3's prior-free first pass
is pointless here rather than wrong: the pass after it reaches the same genotype, so it costs a pass
and changes nothing.

**The outer round is the exception, and it runs the other way round.** Re-fitting the slippage
numbers needs reads at the tract, not samples in the cohort, so one deep sample is where it has most
to work with and a wide shallow panel is where it has least
([`read_likelihoods.md`](read_likelihoods.md) §6.1). Switched on, a single 300× sample can therefore
run several outer rounds while the 63-accession tomato panel at 3 reads a position collapses to the
per-stratum values and runs one. **Cohort size is the wrong axis to reason about this on**, and it
is the axis the rest of this section is about.

**No branch on cohort size is needed to get this**, and none should be written: the subtraction is
zero by arithmetic, so a loop that checks convergence will stop after one pass on its own. What
*is* worth a branch is the wasted second pass, and that is an optimisation rather than a
correctness rule.

**At several thousand samples the loop is linear in cohort size and holds no shared state**, so a
locus can be called on one worker while its neighbours are called on others. What grows is the
per-pass cost: `n` priors built and `n × genotypes` rows scored (§2 counts both). **The prior is the term that
cannot be shared**, for the leave-one-out reason of §2, and §11 records the fast path that would
recover it.

**At three reads a position** the read likelihoods are nearly flat across genotypes, so the prior
carries the call and the loop settles quickly — there is little in the data for the frequencies to
move towards. **At three hundred** the likelihoods are overwhelming, the prior is irrelevant to the
call, and the loop settles quickly for the opposite reason. **The slow case is the middle**, where
reads and prior are comparable, which is worth knowing before anyone reads a convergence-failure
count as a data-quality signal.

---

## 8. Cost, memory and determinism

**What dominates.** Two costs, and they scale differently. The read likelihood is computed once per
`(sample, observation, candidate)` — or once per re-fit round where §5.1's per-locus re-fitting is
switched on, which is one build plus at most three more; the loop's own arithmetic is
`passes × samples × genotypes`.
Since genotypes grow as `C(A + P − 1, P)` while candidates grow as `A`, **the loop's arithmetic
overtakes the likelihood as the allele count rises**, and the crossover depends on how expensive
each side's per-entry work is. A candidate's read likelihood on the STR path costs a stutter term
and a substitution alignment; the loop's per-genotype work is not a multiply-add either, because the
prior puts an `lgamma` in it (§2). Neither has been measured for ng.

**What is held.** Per locus: the read-likelihood table, `samples × genotypes` floating-point values;
the current and previous expected copies, `2 × alleles`; and one sample's posterior row,
`genotypes`. **Nothing is allocated inside the loop** — the caller hands in scratch sized by the
locus's shape and the loop fills it. Production lifted exactly these buffers out of its own
iteration after a profile put the allocator's self-time at about 16% of cycles
([`posterior_engine.rs:1874`](../../../../src/var_calling/posterior_engine.rs)).

**Determinism has one rule and it is the M-step's.** The cohort's expected copies are a sum over
samples, floating-point addition is not associative, so **the sum runs in a fixed sample order**.
That is the contract [`calling_priors.md`](calling_priors.md) §8 names and it is the same
requirement the merge already carries for byte-identical output at any worker count
([`run_streaming.md`](run_streaming.md) §12). Everything else in the loop is a pure function of its
inputs: no RNG, no clock, no thread-dependent iteration.

**Errors.** A locus that does not converge is not an error (§6). The failure modes that *are* errors
are caller bugs — a genotype table that disagrees with the allele count, a non-finite likelihood
reaching the loop, a sample count that disagrees between the evidence and the parameters. These are
assertions, and the structural ones hold in release.

---

## 9. Where this runs, and what it hands on

**Inside the merge's builder, on the region the builder owns**, rather than after the whole region
stream has drained. [`cohort_merge.md`](cohort_merge.md) §14's sixth open question raises this and
leans the same way; this document takes it, for that document's reason — the buffer between the
workers and the ordered output then holds *called variants*, which are smaller than the observations
they came from, and the run's skeleton already collects results per region
([`run_streaming.md`](run_streaming.md) §3.5).

**It commutes either way**, so this is a memory decision rather than a correctness one, and it can be
revisited when emission fixes the shape of what a call is.

**What it hands on** is, per locus: the allele table, and per sample the most probable genotype and
its confidence, plus the cohort's expected allele copies — **plus, on the SNP/indel path only, for a
sample the candidate step declared uncallable (§5.0), the fact that its genotype is missing rather
than any of the above.** Emission writes that sample's `GT` as missing; it is a decision taken
before the reads were scored and not a low-confidence call, and the two must not be conflated in
the output.

**The posteriors themselves are not handed on** *(amended 2026-08-25; this sentence used to say they
were)*. `CallingScratch` holds one genotype-length posterior row that each sample in turn is scored
into, so when the last sample has been scored the earlier rows no longer exist; keeping them would
mean a second `samples × genotypes` buffer held for the whole locus. Nothing downstream asked for
them. What *did* ask — the site quality, which needs the likelihood table — is answered by computing
it here instead ([`calling_quality.md`](calling_quality.md) §3.2), not by widening the hand-off.

**So two more things leave with each locus, both owned by
[`calling_quality.md`](calling_quality.md)** *(added 2026-08-25)*: the site quality before its
artifact correction, computed from the likelihood table once the loop has stopped; and nine pooled
read counts — reference and primary-alternative reads, their forward-strand and placed-left shares,
the locus total, and how many alternative reads the called genotypes imply — which are all the
artifact correction needs and are the same nine numbers whatever the cohort size. That document's
§3.5 also pins one rule on whoever consumes them: **there is one quality field, and it is written
twice** — the baseline here, the corrected value at the first output stage — never two fields.

**Back to the uncallable sample, because it is an absence in all of this too: it is scored at no
point, so it has neither posteriors nor expected copies.** §5.0 sets it aside before the first pass,
so this is an absence rather than a value with a flag beside it, and the cohort's expected copies
are a sum over the samples the locus was called on — the same cohort the site quality's count axis
runs over ([`calling_quality.md`](calling_quality.md) §5.1). **A repeat tract never produces such a
sample** (§5.0.1): every covering sample there is handed on with a genotype and a confidence, and
the expected copies are a sum over all of them. **The cohort's expected allele
copies are not a by-product** — site filtering and emission read them, and recomputing them
downstream from the called genotypes would give a different number, because a called genotype has
thrown away the uncertainty the expected copies still carry.

---

## 10. Reuse map

| what | production code | how ng reuses it |
|---|---|---|
| the loop's shape — E-step, M-step, convergence delta behind one seam | [`posterior_engine.rs:2635`](../../../../src/var_calling/posterior_engine.rs) | **shape ported.** The seam is what lets the two paths share a loop while differing in their emission, which is the same split [`read_likelihoods.md`](read_likelihoods.md) §2 makes |
| the prior-free first pass | [`posterior_engine.rs:2586`](../../../../src/var_calling/posterior_engine.rs) | ported with its reasoning (§3) |
| convergence on expected counts, and the emit-with-flag rule | [`posterior_engine.rs:26`](../../../../src/var_calling/posterior_engine.rs), [`:2657`](../../../../src/var_calling/posterior_engine.rs) | ported; both constants inherited and marked soft (§6) |
| the frequency update | [`em.rs:816`](../../../../src/ssr/cohort/em.rs) | the STR path's version of the same M-step; ng writes one for both paths |
| the final pass that produces calls | [`em.rs:857`](../../../../src/ssr/cohort/em.rs) | ported as §2's last step |
| scratch lifted out of the iteration | [`posterior_engine.rs:1874`](../../../../src/var_calling/posterior_engine.rs) | ported, with its measured reason (§8) |
| the outer round that re-fits the slippage numbers and rebuilds the likelihood table | [`em.rs:572`](../../../../src/ssr/cohort/em.rs) | **shape ported, default not** — ng builds the same nesting and ships with the round count at zero, so the setting is a configuration rather than a code path (§5.1) |
| the re-fit itself — direction split, fall-off, pulled back toward the per-stratum value | [`stutter.rs:184`](../../../../src/ssr/cohort/stutter.rs) | **shape ported, input not** — ng weights each read by the genotype posteriors where production attributes it to the called genotype (§5.1) |

**Parity oracle, and it is one loop rather than two.** Production's **SNP/indel** loop is the oracle:
given the same likelihood table, the same prior parameters and the same candidate set, ng must
reproduce its genotypes, and any difference must be attributable to a decision one of the three
calling documents records.

**Production's STR loop cannot serve as one, and the reason is not a bug in either.** It converges on
a different quantity at a different scale — the normalised frequencies `π`, at a tolerance of `1e-6`
([`em.rs:137`](../../../../src/ssr/cohort/em.rs)) — where ng converges on expected copies over
chromosomes at `1e-3` (§6). Two loops that stop at different points on the same trajectory will
disagree on any genotype near a boundary, for no reason either document records. **What the STR side
gets instead is a differential:** run ng's loop against production's convergence rule and tolerance
on the same likelihood table, and require the genotypes to match; then change the rule back and
report what moved. That has a failing state, which parity-with-an-escape-clause does not.

---

## 11. Deferred, with a recommended home

- **A shared prior for large cohorts.** The leave-one-out correction to a sample's prior shrinks as
  the cohort grows, and past some size one prior would serve every sample, turning `n` evaluations a
  pass into one. [`calling_priors.md`](calling_priors.md) §10 homes this here as a perf item with a
  measured threshold, and the threshold is what is missing: nobody has measured at what cohort size
  the shared prior stops changing a genotype. **Home:** here, once there is a cohort large enough to
  measure it on.
- **Discovering alleles from the calling — no longer deferred.** It is built here, off by default,
  and measured (§4.1, §12's Q3). Note what it is *not*: it adds alleles to a locus whose extent is
  already settled. Growing the extent itself is excluded, not deferred — §1.2.
- **Re-fitting the part-repeat slippage numbers at the locus** — the three HipSTR moves and ng does
  not (§5.1). It cannot be measured until ng specifies them rather than inheriting production's
  placeholders. **Home:** the parameters fit, which
  [`read_likelihoods.md`](read_likelihoods.md) §10 names as the right owner and records as unclaimed;
  once it claims them, §12's Q2 gains a fourth setting.

---

## 12. Open questions

**Q1 — is the summarise-and-condition loop of §2 the right shape, and is the Dirichlet prior the
right prior?** Two questions, because the alternative changes both at once. **All the arms are to be
built and compared** (owner, 2026-08-21).
§2.1 sets out the two: ours estimates one cohort-wide quantity and scores each sample against it,
paying a leave-one-out subtraction to stop a sample's reads counting twice; freebayes scores one
genotype-per-sample assignment at a time under a joint prior — a combinatorial term for the
arrangement and Ewens' sampling formula for the allele counts — and needs no such subtraction
because it never estimates a frequency.

**Three arms, because scoring assignments and using Ewens' prior are two separate changes and
running them together answers neither question.**

| arm | how the cohort is handled | the prior over a cohort assignment |
|---|---|---|
| **A — the default** | summarise the others, condition each sample on the summary | Dirichlet-multinomial per sample, leave-one-out ([`calling_priors.md`](calling_priors.md) §6) |
| **B** | score whole assignments | the **same model, written jointly** — the Dirichlet-multinomial over the cohort's own allele counts, which has a closed form and needs no leave-one-out because nothing is conditioned on a summary |
| **C** | score whole assignments | freebayes' — the arrangement term times Ewens' sampling formula (§2.1) |

**A against B isolates the shape**; **B against C isolates the prior.** Run only A against C, as an
earlier version of this question said to, and a difference is unattributable: the two differ in both
at once, and by a lot — at one diploid sample with two alleles Ewens gives the two homozygotes
**equal** weight, `1/(2 + θ)` each, where the Dirichlet seed puts about `1 − 3θ/2` on hom-ref, which
is roughly 30 Phred apart on hom-alt at `θ = 0.001`.

**Run every arm at `F = 0`, and that is a property of the measurement rather than of the caller.**
Arms B and C have nowhere to put a per-sample inbreeding coefficient: both joint priors are written
over the *cohort's allele counts*, while arm A applies a two-branch identity-by-descent mixture to
each sample separately (§3.2 of [`calling_priors.md`](calling_priors.md)). **Nor does one compose
cheaply.** The identity-by-descent branch means a homozygous-by-descent sample contributes **one**
chromosome to the cohort's counts rather than two, so the exact joint prior is a mixture over which
homozygous samples took that branch — `2^k` terms at `k` homozygotes, which is not a closed form at
63 accessions. Left unstated, the comparison would run arm A with tomato's `F ≈ 0.8` against two arms
with no `F` at all on a panel that is homozygous nearly everywhere, and the B-minus-A difference
would be mostly the inbreeding term. **Setting `F = 0` in all three costs nothing and is what makes
the comparison mean what it says.**

**What that defers, stated so it is not lost: whether a whole-cohort scorer can carry inbreeding at
all.** If a joint arm wins at `F = 0`, that question has to be answered before it could replace arm
A, because a caller that cannot express `F` is not usable on a selfing crop — and the answer may be
that it cannot, which would decide against the joint arms on its own. **Home:** here, as the
follow-up Q1 hands on if it does not end in arm A.

**The fourth cell — summarise-and-condition with an Ewens prior — is not in the table because it does
not exist in closed form.** Ewens scores an unlabelled partition of the whole cohort's chromosomes;
there is no per-sample factor to hold one sample out of. **That is a finding, not an omission:** the
frequency-spectrum prior is only available to a joint scorer, so if arm C wins, the shape decision
follows from it rather than being separable.

**One property of arm C to read its result against.** The composition is not a distribution over
assignments: summed over every assignment at a two-allele, one-diploid-sample locus it gives
`(2 + θ)/(1 + θ)`, about **2**, because the arrangement term normalises within a labelled count
vector while Ewens is over unlabelled partitions, so each partition is counted once per labelling.
Different partition shapes admit different numbers of labellings, so the distortion is not a constant
that cancels when the posterior is normalised. **Report it rather than repairing it silently** — a
normalised variant is a fourth arm, not a bug fix.

**Build it in four pieces, in this order.**

1. **The two joint priors, as pure functions of an assignment's per-allele chromosome counts.** Both
   are small: freebayes' whole Ewens implementation is 51 lines
   ([`Ewens.cpp`](../../../../freebayes/src/Ewens.cpp)), and the cohort Dirichlet-multinomial is the
   primitive ng already ports ([`calling_priors.md`](calling_priors.md) §9) evaluated on cohort
   counts instead of one sample's.
2. **An exhaustive scorer for small cohorts** — enumerate every assignment, sum the joint posteriors,
   and read off each sample's genotype by marginalising. At eight samples and three genotypes that is
   6,561 assignments; at five samples and 21 genotypes, 4.1 million. **This is the piece that
   separates two things a search would confound:** whether a joint model disagrees with our loop at
   all, and whether the search's neighbourhood is too narrow. If arm B's exhaustive scorer agrees
   with arm A on small cohorts, every difference at large cohorts is the search, not the model.
3. **The search, for the cost comparison at realistic cohort sizes** — the local moves and the
   convergence rule of §2.1, shared by arms B and C.
4. **Arm C's prior into the same scorer.** Once arms A and B are settled this is a swap of one
   function, which is why it costs almost nothing to answer the second question too.

**What to report, per arm.** Peak memory per locus and wall-clock per locus, at one sample, at the
63-accession tomato panel, and on a cohort of a thousand; **plus the count of genotypes that differ
from arm A**, because cost alone cannot decide it where the arms disagree. **Report B-against-A and
C-against-B separately** — those are the shape difference and the prior difference, and a single
A-against-C number is the one figure that cannot be acted on. *Leaning: the arms will agree on the
easy majority and the comparison will turn on the hard minority, so report the disagreement counts
first and the timings against them.*

**One prediction worth stating so the measurement can refute it.** Freebayes copies the entire
assignment for every neighbour it scores ([`Genotype.cpp:916`](../../../../freebayes/src/Genotype.cpp)),
so its retained set is one `n`-long assignment per neighbour — **quadratic in cohort size on its
exhaustive-local path** (`n × (genotypes − 1)` neighbours), against this loop's
`samples × genotypes` likelihood table, which is linear. **That is a property of freebayes'
implementation and not of the method**: every neighbour differs from the winner in a bounded number
of samples, so a native implementation can store the moves rather than the assignment and stay
linear. **Build it natively** — the repo has been bitten before by porting an upstream shape instead
of writing one for its own data layout. If the measurement reports a quadratic curve, the first thing
to check is whether ng copied the copy.

**Q2 — should the slippage numbers be re-fitted at the locus, and if so how hard should the locus be
allowed to pull them?** **Settle first what the level is pulled back *toward*.** This question and
§5.1 were written against a slippage level fitted per `(motif period, repeat count)` cell;
[`str_slippage_level_curve.md`](str_slippage_level_curve.md) has since made it a curve in repeat
count fitted once per motif period. The re-fit's shrinkage target is therefore a point on a fitted
line, not a cell's own estimate, and the two differ most exactly where the old cell was thin — which
is where a per-locus re-fit was most likely to look useful. **The arms below are unchanged; what
changes is the baseline they are measured against and the meaning of "far from its stratum".**

Three settings to run, and they are one code path at three configurations
rather than three implementations (§5.1): **frozen** — the per-stratum numbers, no rounds, what this
document specifies; **pulled back part way** — production's 50 pseudo-counts on the shape and 20
slipped reads on the level, at most three rounds; **free** — the locus's own reads set the numbers
outright, HipSTR's setting, which in the nested shape is the middle one at zero pull-back. Hold the
number of re-fitted parameters at three and the re-fit's input at the genotype posteriors in all
three arms (§5.1), so the only thing varying is the pull-back.

*Leaning: ship frozen and let the measurement move it, and expect the answer to depend on reads at
the locus rather than on cohort size* — at 300 reads on one sample there is enough at a tract to
fit from, at 3 reads across 63 accessions there is not, which is what the pull-back exists to
refuse ([`read_likelihoods.md`](read_likelihoods.md) §6.1).

**Settled by:** the two benchmarks the STR path already has — the HG002 tandem-repeat bundle, where
the truth is an assembly and the depth is high, and tomato's recurrence-based standard
([`silver_standard.py`](../../../../benchmarks/ssr_tomato1/scripts/silver_standard.py)), where the
truth is weaker but the depth is 3 reads a position across 63 accessions. Between them they cover
both ends of the axis that matters. **Report three numbers, not one:** genotype accuracy given
detection under each setting; the share of loci whose re-fitted numbers land far from their
stratum's, which is what says whether adaptation can matter at all; and the wall-clock cost, since
each round rebuilds the whole read-likelihood table (§8). **This one is not cheap and it is not
blocked** — the data exists and the work is a measurement harness.

**Q3 — should the loop be allowed to discover alleles hiding under stutter, the way HipSTR does?**
**Build it and measure it on real data** (owner, 2026-08-21). §4.1 sets out the mechanism: after the
loop converges, retrace what the model is explaining as slippage, and admit a tract length that
recurs in one sample past a bar; realign only the new alleles' columns; run again; stop when a round
adds nothing; then drop the alleles no sample's best genotype used. ng ships with it off.

**Three settings, because the middle one may be the answer:** off; discovery against the converged
frequencies **held fixed**, so each round is a scoring pass; and discovery against a full convergence
every round (§4.1). The middle setting is HipSTR's cost with a prior that has seen the locus, and if
it matches the third on accuracy it is the one to keep.

**What to report, and the first number is not accuracy.** *How often it fires at all* — the share of
loci where a round admits anything, on each benchmark. If that share is a handful of loci in ten
thousand, everything after it is a rounding error and the option can be dropped rather than tuned.
Then, on the loci where it fires: genotypes that change, alleles that survive the prune, discovery
rounds per locus, and wall-clock against the loop with discovery off.

**Where it should pay, and where it should be dangerous — the two must be reported apart.** *Should
pay:* the HG002 tandem-repeat bundle, high depth, where a short allele masked by contraction slip is
exactly the case an assembly-based truth set can adjudicate and a caller usually misses. *Dangerous:*
the tomato panel at 3 reads a position, where the fraction half of HipSTR's bar is inert — 2 reads
out of 3 is 67% — so admission rests on a 2-read count alone, which is also what one stutter product
looks like. **Sweep the bar** — HipSTR's 2 reads and 15% were set on high-depth human data and are
inherited, not measured here (§4.1) — and report **the depth at which the fraction starts to bind**,
which is where the two halves change places (about 13 reads at these values).

*Leaning: expect it to fire rarely and to matter where it fires.* If the tomato panel shows it
minting alleles at low depth and HG002 shows it recovering real ones, the answer is not on/off but a
depth floor below which discovery does not run, and Q3 should return that floor rather than a verdict.

**Q4 — are the inherited threshold and iteration cap right for this caller?** `1e-3` and 50 passes
come from production, where the observed need is 3 to 5 passes on GATK reference data (§6). Neither
has been measured here, and this caller's range reaches places production's did not — a thousand
samples, three reads a position. *Leaning: inherit both and instrument rather than tune* — emit the
pass count per locus, and set the constants from the distribution once there is one. **Settled by:**
the pass-count distribution on the tomato panel and the GIAB trio, plus the count of loci hitting the
cap. **Cheap, and it needs no truth set.**

**Q5 — where does the loop's arithmetic overtake the read likelihood?** §8 says the crossover exists
and that neither term has been measured. It decides whether the allele cap matters for speed as well
as for memory. *Leaning: none — this is a measurement, not a judgement.* **Settled by:** timing one
locus at several allele counts on both paths, which needs no data beyond a fixture.

**Q6 — is the wasted second pass at one sample worth a branch?** §7 shows the loop terminates
correctly at one sample without any test of cohort size, at the cost of one extra E-step. *Leaning:
measure before branching* — a branch on cohort size is exactly the kind of thing that later grows a
second meaning, and the cost is one pass over one sample. **Settled by:** the single-sample profile
the STR benchmark already produces.

**Q7 — does starting the first pass without a prior change any call on this project's data?** §3
takes it from production and GATK and argues it from the size of the seed's pull, which is arithmetic
rather than a count of calls. *Leaning: keep it — the failure it guards against is a converged wrong
answer rather than a slow one, and that is the expensive kind.* **Settled by:** re-calling one tomato
interval and the GIAB single-sample bundle with the first pass seeded instead of flat, and counting
genotypes that move and variants that disappear. Cheap, needs no truth set to *detect* a difference,
and needs GIAB only if there is one to adjudicate. **Where to look if the answer is "none":** the
tomato panel at 3 reads a position, since the guard is worth least at high depth and §3 predicts it
is worth most there.

---

## 13. How we know it works

**Unit tests, each pinning a property rather than a value.**

1. **One sample reaches its fixed point on the loop's first pass, and the second pass only confirms
   it.** The prior-free initialisation and the first seeded pass give *different* expected copies —
   the seed is worth 20 to 30 Phred (§3), so they must — and from the first pass on nothing moves:
   pass 2's expected copies equal pass 1's **bit for bit**, and the loop stops. Assert that equality
   and the pass count of two; **do not assert that the copies equal their initial value**, which is
   false by construction, and do not assert a pass count of one. Whether a branch skips the second
   pass is Q6's; this test must pass either way, so it asserts on the *genotype* and on
   pass-1-equals-pass-2, not on the pass count when the branch is present.
2. **The M-step is order-independent in its inputs and fixed in its arithmetic.** Permuting the
   samples changes no genotype. And the sum runs in a fixed order: **the mutation check is on the
   summed expected copies, compared bitwise**, because reordering a floating-point sum moves them by
   a few units in the last place and an argmax over 21 genotypes will not flip on that — a mutation
   test whose observable is the genotype passes against a non-deterministic implementation and
   proves nothing.
3. **The prior-free first pass does what it is for.** On a cohort whose reads support a non-reference
   allele at a locus, the expected allele copies after the first pass reflect the reads and not the
   seed. **And the property §3 actually claims is the stronger one:** seeding the first pass instead,
   at a locus thin enough that the seed outweighs each sample's reads, must converge to no variant
   where the prior-free start converges to one. A test that only checks the first pass has not tested
   the trap.
4. **A locus that hits the cap is emitted with the flag set**, rather than dropped or raised as an
   error, and the flag reaches the output (§6). **Do not also assert that the convergence delta
   decreases every pass.** Expectation-maximization guarantees a monotone *likelihood*, not a
   monotone parameter delta, and §6 claims no such thing — HipSTR's own loop records the same caveat
   about its likelihood under pseudocounts
   ([`em_stutter_genotyper.cpp:196`](../../../../HipSTR/src/em_stutter_genotyper.cpp)). A test
   written that way is either flaky or pinned to one hand-picked fixture.
5. **The likelihood table is computed once per set of slippage numbers, never once per pass.**
   Instrument the emission call count for one locus and assert it equals
   `candidates × Σ_s (observations in sample s) × builds`, independent of the number of passes —
   **`Σ_s`, not a three-way product**, because samples do not have equal observation counts and a
   fixture built so they do is the one shape that hides the bug. With re-fitting off, `builds = 1`.
   With it on, **`builds` is the number of rounds that actually rebuilt**, which is one less than the
   number of rounds run whenever the last round detects convergence and leaves before rebuilding —
   production's loop breaks in exactly that order
   ([`em.rs:575-580`](../../../../src/ssr/cohort/em.rs)). Assert against an instrumented build
   counter, not against a round count. **This is the test that catches the whole design being paid
   for twice**, and it fails silently in production terms — a loop that recomputes gives identical
   answers, only slower. **Run it at both settings**: a per-pass rebuild hides behind a per-round one
   when only the re-fitting arm is tested.
6. **A re-fit that changes nothing changes nothing.** With per-locus re-fitting switched on at a
   locus whose reads carry no slips, the re-fitted numbers come back equal to the per-stratum ones
   and the genotypes are bit-for-bit those of the frozen setting. This pins the collapse production
   relies on ([`stutter.rs:184`](../../../../src/ssr/cohort/stutter.rs)) and catches a re-fit that
   drifts on an empty profile.
7. **Nothing allocates inside the loop.** The allocation count over a locus is independent of the
   pass count, and — with re-fitting on — of the round count.
8. **The parity oracle** (§10): the same inputs give production's genotypes, or the difference is
   traced to a recorded decision.
9. **On a cohort small enough to enumerate, the search reproduces the exhaustive scorer's
   genotypes** (§12, Q1) — or the difference is the search's narrowness and is reported as such
   rather than as a model difference. **Do not add "and at one sample the exhaustive scorer agrees
   with §2's loop".** It does not, and the reason is not marginalisation: Ewens' formula is
   symmetric in the alleles, so for one diploid sample it gives the two homozygotes equal weight —
   `1/(2 + θ)` each after normalising, against a heterozygote's `θ/(2 + θ)` — while this loop's seed
   puts about `1 − 3θ/2` on hom-ref and `θ/2`-scale mass on hom-alt
   ([`calling_priors.md`](calling_priors.md) §4). At `θ = 0.001` those differ by about 30 Phred on
   hom-alt. **The two approaches differ in their prior as well as in how they marginalise**, which is
   Q1's first problem to solve, not a test to write.
10. **The joint prior's normalisation is stated, and it is not one.** Summing the exponentiated joint
    prior over every assignment at a two-allele one-diploid-sample locus gives `(2 + θ)/(1 + θ)`,
    about **2** — because the arrangement term normalises within a labelled count vector while Ewens'
    formula is a distribution over *unlabelled* partitions, so each partition is counted once per
    labelling. freebayes never normalises this prior; it normalises the posterior over the
    assignments it kept ([`Genotype.cpp:1595`](../../../../freebayes/src/Genotype.cpp),
    [`Marginals.cpp:41`](../../../../freebayes/src/Marginals.cpp)). **Assert the closed-form total for
    a small enumerable locus, whatever it is** — a test asserting it sums to one fails against a
    correct implementation.
11. **Discovery finds a planted allele, terminates, and costs nothing when off** (§4.1). Three
    properties. On a made-up locus where one sample carries a short allele whose reads the frozen
    slippage numbers explain as contraction, discovery must admit that length and the sample must end
    up heterozygous, where the same locus with discovery off calls it homozygous. **Termination is the
    half that will bite:** a locus whose reads support no further allele must end after one round, and
    a locus that keeps admitting must stop at the allele cap rather than run — assert the round count,
    not just that the call is right. And with discovery off, the emission evaluation count and the
    genotypes must be **bit-for-bit** those of test 5, so the option cannot cost anything when unused.
12. **A discovered allele's columns are appended, not rebuilt.** After a discovery round admits `k`
    alleles, the emission evaluation count has risen by exactly the entries those `k` alleles need and
    by nothing else — the property HipSTR gets by mapping surviving haplotypes to their new indices
    (§4.1). Without this test the append is an intention rather than a fact.

**The end-to-end check, and the definition of done for the manager:** the two benchmarks the sibling
documents already use — GIAB single-sample genotype accuracy at true variants, and the HG002
tandem-repeat bundle scored on genotype accuracy given detection. **This document owns neither
number**; what it owes is that changing the loop moves them only for reasons it records. A loop
change that moves either number without a decision behind it is a defect in the loop.
