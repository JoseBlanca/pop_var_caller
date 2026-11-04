# ng — reads in the pass's haplotype draw

*Spec, 2026-11-04. No code yet. It realises the item deferred in
[`neighbour_evidence_read_phase.md`](neighbour_evidence_read_phase.md) §8, "the pass's sampler using
reads", whose condition is now met: on the outbred simulated cohort, links that rest on the
population are still wrong 3.6 times in 100 with reads in the link test, against a ceiling of 1 in
100 (`doc/devel/reports/implementations/read_phase_2026-11-03.md`). It changes one step of the
phasing model of [`neighbour_evidence_phasing.md`](neighbour_evidence_phasing.md) §2.1, and reads
the fragments that [`neighbour_evidence_read_phase.md`](neighbour_evidence_read_phase.md) carries to
the pass.*

---

## 1. What this is

### 1.1 The problem

The neighbour-evidence pass phases a sample by copying: along the window, each of the sample's two
copies follows a template, a haplotype of another sample. Each round, the pass draws one **path**
for the sample (which templates the copies follow, record by record) and then estimates the
sample's two haplotypes along that path. Other samples copy those haplotypes in the next round, and
the paths of the last four rounds make the written phase (phase output spec §2.1, §2.3).

Reads enter that model only at its very end. Since the read phase build, a link test after the
rounds lets reads break a population link they contradict, and join two records the population
could not orient (read phase spec §2.7). **Nothing earlier sees them**: the haplotypes a sample
passes on to the cohort are estimated from its templates and its genotype likelihoods alone. So
where a sample's reads show that two heterozygotes lie on one copy, the other samples copying it
still read "either allele, half and half" at both records.

The consequence, on the founder simulator with dense SNPs (one every 100 bases, 50 samples at 3×,
`read_phase_2026-11-03.md`):

| simulated cohort | links wrong, population only | with reads in the link test | links reads broke, of them wrong |
|---|---|---|---|
| outbred, 64 founders | 1,466 of 28,007 (5.2%) | 1,082 of 30,321 (3.6%) | 385 of 442 |
| outbred, 16 founders | 69 of 41,978 (0.16%) | 47 of 44,407 (0.11%) | 23 of 48 |

Reads break mostly wrong links (87% of those they broke at 64 founders), but they can only test the
links they reach, and 3.6% stays above the ceiling. The links they do not reach are wrong because
the population's haplotypes are poorly phased. **The hypothesis this spec tests** is that they are
poorly phased partly because each sample's reads never reached them.

### 1.2 What changes

**Where a sample's reads link two of its heterozygous records, the haplotypes drawn along its path
put the two alleles on the copies the reads favour.** The reads decide the *relative* orientation of
heterozygotes they link; they do not change any genotype, any path, or how the haplotypes are
estimated anywhere else. Every later step reads the changed haplotypes:

- **other samples' templates**, from the next round: a template now carries the sample's read phase;
- **the phase consensus and the path agreement**: the 36 paths the link test compares carry the
  read phase too, so a population link and a read link are no longer two separate tests at the
  records reads reach.

This is the method of SHAPEIT2's phase-informative reads (research note
[`05_read_backed_phasing.md`](../research/phasing_imputation/05_read_backed_phasing.md) §3.1): a
read likelihood multiplied into the phaser's draw of haplotypes, with a disagreeing read charged its
error, not vetoed. SHAPEIT2 multiplies it into the draw of the whole haplotype pair; this spec
multiplies it into the second half of the pass's draw only, the haplotypes given the path (§2.4).

### 1.3 Goals

- At the dense outbred cohorts above, fewer wrong links among those the population makes, at the
  same or a larger share of heterozygotes phased. The ceiling of 1 wrong in 100 is the target; §6
  says why it may not be reached by this change alone.
- No more wrong genotypes than now, on any cohort, beyond the seed spread (about 5%, phasing spec
  §1).
- **Byte-identical output wherever no sample has a fragment**: a run with only STR records, an
  input with no read pairs reaching two records, a cohort where refinement does not run.
- Deterministic, the same at any thread count, as the pass is now.

### 1.4 Non-goals

- **Reads in the path itself** — which templates a sample follows. Considered and deferred (§2.4,
  §7).
- **Writing anything new.** No VCF field, header value or option changes. The `##neighbourEvidence=`
  line already states the read settings (`read_reach`, `read_cap`).
- **Changing the link test** of read phase spec §2.7. Its read rows stay; §6 measures how often they
  still fire.

### 1.5 Vocabulary

The terms of the phasing and read phase specs, plus four:

- **Joint at a record**: for one sample, given its path, the probability that copy 0 carries allele
  `a` and copy 1 allele `b`, for every ordered pair `(a, b)` (phasing spec §2.1, `copies_along` in
  `src/neighbour_evidence/model.rs:677`). Each copy's haplotype is a margin of it.
- **Heterozygous mass**: the joint's total over `a ≠ b` — how likely the sample is heterozygous at
  the record, given its path.
- **Orientation** of a heterozygous record: *straight* when copy 0 carries the lower-numbered allele
  of each heterozygous pair, *swapped* when copy 1 does. Two records with the same orientation are
  in cis, in the read phase spec's sense; with different ones, in trans.
- **Read-linked group**: for one sample along one path, a run of consecutive orientable records
  (§2.2), each linked to the next by at least one fragment.

---

## 2. The method

### 2.1 In words

Today, at each record along the drawn path, the joint is computed on its own and each copy's
haplotype is its margin. Where the two templates are alike there, the joint gives the two orders
`(a, b)` and `(b, a)` the same weight, both copies get the same margin, and the record carries no
phase. Where the templates differ, one order outweighs the other, and that sets the orientation.

With reads, for each sample, along each drawn path:

1. **Compute every record's joint, as now.**
2. **Find the read-linked groups**: the records whose heterozygous mass is at least one half, in
   order; consecutive ones are linked when the sample has fragments at both. A record outside every
   group is left exactly as now.
3. **In each group, draw an orientation for every record**, from a two-state chain:
   - each record's own weight on *straight* is the joint's mass on the orders that put the lower
     allele on copy 0, and on *swapped* the mass on the reverse orders. This is the templates' say;
   - each link between consecutive records multiplies `e^{L/2}` into the two records having the same
     orientation and `e^{−L/2}` into their having different ones, `L` being the link log-odds of the
     read phase spec (§2.5 there), computed from this sample's fragments at the two records;
   - the orientations are drawn backwards from the group's last record, each given the one after it,
     as the path itself is drawn.
4. **Write each grouped record's haplotypes with the drawn orientation**: the joint's heterozygous
   mass moves onto the drawn order, and its homozygous mass stays where it is. Each copy's haplotype
   is the margin of that joint, as now.

The sample's genotype posterior is not touched: it comes from the walk, before any draw. What
changes is what other samples copy and what the phase step reads.

### 2.2 The chain, precisely

For one sample along one path, write `J_r(a, b)` for the joint at record `r`, normalised to sum to
one, and `H_r = Σ_{a≠b} J_r(a, b)` for its heterozygous mass.

**Orientable records.** A record is orientable when `H_r ≥ ½` (`HET_MASS_MIN`) and it is not one of
the sample's hidden records (phasing spec §4.1). Below one half the path says the sample is more
likely homozygous there, and the record's orientation would be drawn for alleles it probably does not
carry.

**The heterozygote a link reads.** At a biallelic record it is alleles 0 and 1. At a record with more
alleles it is the pair `a < b` with the largest `J_r(a, b) + J_r(b, a)`, ties to the lower pair in
the order the loops of `copies_along` visit them; that is the `Heterozygote` handed to
`score_read_link` (`src/read_phase/link.rs:192`).

**Links.** For consecutive orientable records `r` and `s` of the sample, `L_rs` is
`score_read_link` over the sample's fragments at `r` and `s`. The pair is linked when at least one
fragment counts there (`ReadLink::fragments > 0`). A fragment that reaches `r` and `t` but skips an
orientable record `s` between them adds nothing: the chain links only consecutive records, as the
link test does. §6 measures what this loses.

**Weights.** Within one group, with `o_r ∈ {straight, swapped}`:

- own weight: `w_r(straight) = Σ_{a<b} J_r(a, b)`, `w_r(swapped) = Σ_{a<b} J_r(b, a)`;
- link weight: `φ_rs(o_r, o_s) = e^{+L_rs/2}` if `o_r = o_s`, `e^{−L_rs/2}` otherwise;
- the chain's probability of a set of orientations is proportional to `Π_r w_r(o_r) · Π_links φ_rs`.

A forward pass over the group, normalised at each record, and a backward draw give an exact draw
from that distribution. One uniform a record of the group, from the stream of §4.

**The haplotypes written.** At a grouped record with drawn orientation `o`, the new joint `J′` is:

- `J′(a, a) = J(a, a)` for every `a`;
- for every `a < b`: if `o` is straight, `J′(a, b) = J(a, b) + J(b, a)` and `J′(b, a) = 0`; if
  swapped, the reverse.

Copy 0's haplotype is `Σ_b J′(a, b)`, copy 1's `Σ_a J′(a, b)`. For a biallelic record with
`J = (x, u, v, y)` over `(0,0), (0,1), (1,0), (1,1)`: straight gives copy 0 `(x + u + v, y)` and
copy 1 `(x, u + v + y)`.

**Why the reads change only the orientation.** A read's allele at each record is already in the
sample's genotype likelihoods, so already in `J`. What a fragment adds that no single record holds
is which two alleles came from one molecule. The chain uses exactly that and nothing else: `L` is
built from whether the fragment's two alleles agree with cis or with trans (read phase spec §2.5),
and the homozygous mass is never moved. Counting the reads again for *het or not* would count them
twice.

### 2.3 What the change does at three kinds of record

- **Templates alike, reads link it.** Today: no phase, both copies the same margin. Now: the group's
  first record takes a fair coin's orientation (its own weights are equal), and each later record's
  orientation follows the reads relative to it. The whole group's phase is relative, which is all a
  phase is (read phase spec §2.7): the coin decides nothing a reader can see.
- **Templates differ, reads agree.** The drawn orientation is the templates' almost always; nothing
  visible changes.
- **Templates differ, reads disagree.** The two pull against each other, each with its weight. One
  fragment at a base error of 1 in 100 is worth `L` ≈ 3.9 (read phase spec §2.5). Two confident
  templates carrying different alleles, each blurred by the mismatch rate `ε` = 0.01 to 0.995 and
  0.005, favour their order by `(0.995 / 0.005)²` ≈ 39,600, a log-odds of about 10.6. So one or two
  fragments do not overturn confident, differing templates; three do. That is SHAPEIT2's "charged,
  not vetoed", and it means **the change acts mostly where the templates are alike or unsure** —
  which is where the haplotypes carry no phase today.

### 2.4 Decisions, with what they beat

- **Reads in the draw of haplotypes, not in the path — recommended.** SHAPEIT2 samples the haplotype
  pair from the population model times the read likelihood, so reads also shape which templates are
  followed. Here the path is drawn without reads, then the haplotypes given the path and the reads.
  That is an approximation of the joint draw: a path the reads contradict is not made less likely.
  The exact version needs the walk's state to remember the orientation of the last heterozygote, two
  to four times the states, so two to four times the walk's cost, the pass's largest. The draw here
  costs a fraction of one walk (§4). Build the cheap one, measure, and take the exact one only if
  §6 falls short.
- **Orientation hard, alleles soft — recommended.** Inside a group the orientation is a drawn
  choice, as the path is; the alleles keep their probabilities. The alternative, keeping both orders
  at their averaged weight, is what loses the phase today: two copies with equal margins carry no
  link, whatever the reads say (phasing spec §2.2, "a drawn path, not the posterior's margins").
- **Only grouped records change.** A record no fragment links keeps today's haplotypes, so a cohort
  without fragments writes today's bytes (§6, point 1). Drawing an orientation everywhere would put
  a coin's phase on every record where the templates are alike, which other samples would then copy
  as if it were evidence.
- **`HET_MASS_MIN` = ½ — a starting value.** It decides which records join a chain. Too low, and
  homozygous records split groups apart; too high, and a heterozygote at 3× whose path is unsure is
  left out. §6 measures ¼ and ¾ beside it.
- **Every round, every path.** Round 1's haplotypes are what round 2's templates are; leaving reads
  out of round 1 would lose a round. The extra paths of the link test (phase output spec §2.3) are
  drawn the same way, so the paths' agreement includes the reads' uncertainty.

### 2.5 What does not change

The walk, its posterior and every genotype it gives; the parameter fit; template choice; the
averaging over rounds; the no-read rules and the fill check; the consensus and its alignment of
rounds (`SamplePhase::add_round`, `src/neighbour_evidence/phase.rs:237`); the link test and its read
rows; the read check; the phase writer of reads alone, where refinement does not run.

---

## 3. Where it sits in the code

| what | where now | change |
|---|---|---|
| the window's fragments | `refine_window_with_reads` (`src/neighbour_evidence/refine.rs:586`), handed only to `final_solution` | also handed to `solve_round`, one sample's slice a record |
| the draw | `draw_haplotypes` and `copies_along` (`model.rs:578`, `:677`) | `copies_along` returns the joint; a new step after the per-record joints applies §2.2, then takes the margins |
| a link's log-odds | `score_read_link` (`src/read_phase/link.rs:192`) | unchanged, called from the new step |
| the stream | `Uniforms`, `StreamPurpose` (`src/neighbour_evidence/random.rs:16`) | a new purpose, `Orientation`, keyed like `Path` plus the path's number within the round |

The draw lives in `model.rs` because the joints do; the walk does not know about reads, so the new
step takes a closure or a small trait that, given two record indices and their heterozygotes,
returns `L`. The window solve supplies it from the sample's fragments, returning nothing at hidden
records. That keeps `model.rs` free of the fragment types, as `walk` is free of the store's.

**Prototype first, in Rust.** The Python phaser reads no reads, and teaching it fragments would mean
writing fragments out, which nothing does. The change is small enough to build behind a constant,
measure on the simulated cells (§6), and only then add the tests and the review.

---

## 4. Cost, memory, determinism

**Time.** Per path, the new step visits each record once to find groups, and each linked pair once
to compute `L`: a merge of two lists of at most 32 read alleles (the read cap) each. The walk spends
`K² ≈ 289` products a record on its pair states. The new step is a few percent of a walk, and a
round draws at most nine paths a sample. Measured, not assumed, in §6.

**Memory.** None held. The fragments are already in the stage (`HeldSites`,
`src/neighbour_evidence/stage.rs:470`); a group's chain is two numbers a record, freed with the path.

**Determinism.** The orientations read uniforms from their own stream, so drawing them changes
neither the path's uniforms nor the extra paths' (the same reason those have separate streams,
phase output spec §2.3). A sample with no group draws nothing from it. `L` is summed in fragment
order (read phase spec §2.5). The draw is the same at any thread count, as the rest of the round is.

---

## 5. Traps

- **Hidden records keep the sample's reads out.** A hidden record is one where the sample is treated
  as having no reads, so the fill check can test the prediction there (phasing spec §4.1). Its
  fragments must not enter a link either, or its reads reach the record through the orientation.
  Hidden records are never orientable (§2.2).
- **"Lower" means the record's allele numbering.** Straight means copy 0 carries the lower-numbered
  allele at every heterozygous pair of the record, which is what `score_read_link` calls cis. At a
  `1/2` heterozygote the lower allele is 1, not the reference; a chain that assumed allele 0 there
  would read every link of that record with the wrong sign.
- **`L` is signed for the heterozygote it was computed on.** At a multi-allelic record the chain's
  orientation applies to every heterozygous pair, but `L` speaks only of the chosen one. Moving mass
  of the other pairs with the same orientation is a convention, not evidence; it is harmless while
  their mass is small, and a three-way heterozygote is impossible in a diploid.
- **STR records carry no fragments.** The STR path names no reads (read phase spec §2.1), so a group
  ends at every STR record that is orientable, and reads never link across one.
- **The first round's templates are calling's `GP`**, not haplotypes. Round 1's draw is still made
  with reads, because round 2 copies its haplotypes.
- **A collapsed duplication** — two copies of a sequence mapped to one place — gives a false
  heterozygote with fragments that confidently link it. The link test keeps such contradictions from
  being written (read phase spec §2.7); here they can shape templates. The hidden-duplication filter
  removes most before the pass; §6 checks HG002's segmental duplications.

---

## 6. How we know it works

1. **Byte identity without fragments.** With refinement on and every fragment removed, a simulated
   cohort's output is today's, byte for byte. The cross-platform digests do not move either, but
   that proves nothing here: refinement does not run on their cohort.
2. **Unit tests:**
   - the chain on hand values: two records with equal own weights and `L` = +3.9 between them draw
     the same orientation with probability `e^{3.9}/(1 + e^{3.9})` ≈ 0.98, checked over many seeds
     against the binomial bound;
   - two confident, differing templates at one record (log-odds about 10.6, §2.3) against two
     fragments and then against three: the templates win in the first case and lose in the second,
     with the probabilities the chain gives;
   - a hidden record is never in a group; a record with `H` just below ½ splits a group;
   - a `1/2` heterozygote at a three-allele record reads alleles 1 and 2;
   - the written joint keeps every homozygous mass and moves the whole heterozygous mass;
   - no fragments: the haplotypes equal today's, bit for bit.
3. **The simulated cells, against today's main** (`grid_table.py`, `read_phase_score.py`, the cells
   of `read_phase_2026-11-03.md`), one seed each, then a second seed where a difference decides:
   - the two dense outbred cohorts: wrong links, population only and with reads; heterozygotes
     phased; wrong genotypes where samples have reads; fills and wrong fills; how often the link
     test's read rows still fire;
   - the nine cohorts of the phase output report: the same counts, as the no-harm check;
   - the 200-sample cohort from 64 founders, 40 at 30×;
   - `HET_MASS_MIN` at ¼ and ¾ on the two dense cohorts;
   - the share of fragment-linked heterozygote pairs the chain loses because a fragment skips an
     orientable record (§2.2).
4. **Cost**: the pass's time with and without the step, at 50 and 200 samples.
5. **Thread invariance**: the end-to-end test of `mode_equivalence.rs` at one thread and four, on a
   fixture with a read-linked group.
6. **Real data, after merging**: the tomato read check (§2.8 of the read phase spec) with and without
   the change, once the cohort's psps are regenerated; HG002's segmental duplications, with the
   duplication filter on and off.

**What to expect.** The read phase build's measurement bounds the direct effect: reads reach the
records of a link in a minority of the population's links at SNP spacings of 1 in 100, fewer at
1 in 1,000. This change improves the haplotypes other samples copy, so its effect reaches beyond the
records reads link; how far is what point 3 measures. SHAPEIT2 reported 22% longer stretches between
switch errors from adding 5–20× reads (research note §3.1), on humans and with 1000 Genomes behind
it. If the dense 64-founder cohort stays above 1 in 100, the next step is the exact draw of §2.4.

---

## 7. Deferred, with a recommended home

- **Reads in the path** (§2.4): the walk's state extended by the last heterozygote's orientation.
  Home: if §6 point 3 leaves the dense outbred cohort above the ceiling.
- **Links across a skipped record** (§2.2): a chain cannot carry a factor between non-consecutive
  records; a group of three with a fragment from the first to the third would need the chain's
  state to hold two orientations. Home: if §6 shows the loss matters.
- **Reads deciding heterozygous or homozygous** — whether a fragment's alleles favour one genotype
  over another beyond what the likelihoods already say. Not planned: they would count the reads
  twice (§2.2).

---

## 8. Questions for the owner

1. **Build the draw given the path first, and the exact draw only if it falls short?** *Recommended:
   yes.* The cheap draw costs a few percent of a walk; the exact one doubles to quadruples the walk,
   the pass's largest cost, and would be built blind to how much the cheap one gives.
2. **Prototype in Rust behind a constant, measure on the simulated cells, then tests and review —
   or the full build with tests first?** *Recommended: prototype first*, as the light development
   loop does: the change is about 150 lines in two files, and whether it earns its place is a
   measurement.
3. **`HET_MASS_MIN` = ½ to start, with ¼ and ¾ measured?** *Recommended: yes.*
