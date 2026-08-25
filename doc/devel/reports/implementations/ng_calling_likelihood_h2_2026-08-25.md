# ng read likelihoods — H2: contamination at a tract, and a term that was reading a different question

*Implementation report, 2026-08-25. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step H2 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). **This completes
Milestone H and the plan — Checkpoint H.***

## 1. What it is

Two things, and the second is the one the plan says not to bundle with anything.

**The lengths the outlier weight is spread over** — `fill_reachable_lengths`, in
[`ssr_emission.rs`](../../../../src/ng/calling/likelihood/ssr_emission.rs). Production spreads
its outlier weight over `D`, the number of distinct sequences the **whole cohort** showed at the
locus ([`em.rs:393`](../../../../src/ssr/cohort/em.rs)), so a single sample showing two sequences
got a junk floor of 0.005 and a 63-accession panel showing twenty got 0.0005 — ten times lower.
A sample's own genotype likelihood moved when an unrelated sample joined the run, which is not a
property a per-sample likelihood may have. ng asks the candidate set and the two slip cutoffs.

**Contamination** — spec §4.5.1's third term:

```text
log Lg(g)  =  Σ_o  n_o · log[ (1 − λ − c) · Σ_a (k_a/P) · Lr(o | a)
                            +      λ      · U(o)
                            +      c      · seed(o) ]
```

A contaminating read at a tract is not junk: it shows a length that is a real allele in somebody,
often a common one. Today it must be explained as slippage where the model can reach it —
inflating the apparent slip rate — or falls to the outlier floor where it cannot, and the first is
how a contaminated sample gets called heterozygous for its contaminant's allele.

## 2. The open question the architecture recorded, and how it is resolved

`arch/read_likelihoods.md` §4.1 and `arch/calling_priors.md` §5 both record that the genotype
prior builds its seed with one entry **per candidate** while `c · seed(o)` needs a probability per
**length**, and list the cases the gap turns on.

**The distribution is keyed to the locus's reachable lengths**, the same support the outlier
weight uses. Two candidates spelling one length share one entry; a read at a length no candidate
reaches gets nothing and falls to the outlier floor, which spec §4.5.1 says explicitly is right
rather than a loss; and a read that ran out inside the tract gets the seed's mass **at or above**
what it witnessed. Converting the prior's per-candidate shape into that stays the calling loop's
job — it is the only place holding both tables.

**Both reviewers who were asked to challenge this agreed, and one gave the argument that settles
the censored case**: truncation is a property of read length against tract length, not of whose
DNA a read is, so zeroing a truncated read's contamination term would put un-modelled
contamination back exactly at long tracts — where truncation happens, and where the term is most
needed.

**Two things the resolution does not settle, now said so in the type's own contract.** Most
entries will be zero, because the reachable support is far wider than the candidate set — about 39
lengths against 5 candidates at a dinucleotide locus. And **how two candidates spelling one length
should share that length's mass is still the specification's open question 3**; keying to lengths
is what makes it answerable in one place rather than decided twice by accident.

## 3. What the review found, and the one finding that was a real defect

**Three agents, 17 mutations, two Blockers and seven Majors.**

### The junk term was answering a different question from the other two

**For a read that ran out inside the tract, the other two terms are tails and the junk term was a
point mass.** `censored_emission` asks *how likely is a tract at least this long*; the seed asks
the same of the contaminant; and the junk term asked *how likely is a junk read to show exactly
this length*. Measured on a locus of 31 reachable lengths, a read witnessing 8 bases got a junk
floor 25 times smaller than the matching tail form, and 91 times smaller than its own
contamination term — **so a truncated read was preferentially explained as somebody else's DNA by
about two orders of magnitude, and the gap widened as the tract lengthened.**

The reviewer filed it Medium and asked whether spec §4.5's *"uniform over the tract lengths the
model's support reaches"* means the probability of the observation or a flat constant. **Spec §2.1
settles it in its own notation: the term is written `λ · U(o)`, a function of the observation.**
So this was a defect, introduced at H1 and made visible by H2 converting one of the three terms.

**Fixing it broke the term's own purpose, and that is worth recording.** With a strict `U(o)`, a
read whose length is outside the candidates' reach collects nothing — and §4.5 exists precisely so
that *"one such read"* cannot drive every genotype to zero. The resolution is that the reachable
count is a **normaliser for how many lengths are in play, not a membership test**: the mass is
floored at one length, so a junk read always has somewhere to go. The seed is a real distribution
and does test membership, which is why the two are not the same expression.

### The two Blockers were both fixtures that could not fail

- **Nothing distinguished the three-term mixture from folding `c` into `λ`** — the wrong answer
  spec §4.5.1 names outright. Doing so left 19 of 19 tests green while moving a genotype by 4.8
  nats; so did dropping the seed entirely (49 nats on one entry); and so did `1 − λ` in place of
  `1 − λ − c`, which drives a row entry **positive**. The new test turns on the seed being
  *peaked*: with a flat seed the three-term form reproduces `λ + c` exactly, and with a peaked one
  the two part company.
- **Neither branch of the seed lookup had a value-pinning test.** The off-by-one at the witnessed
  length survives and halves the survival mass; the existing censored test could not see it,
  because its seed put zero mass there.

### The other Majors

- **The per-read-group fraction was never exercised across two groups** (17 nats under mutation).
  The trap the reviewer flagged is real: the ordinary two-group fixture gives each group its own
  slippage row, so the rows differ whether or not the fraction is read — a fixture sharing
  parameters was needed.
- **A seed's values were unchecked.** Unnormalised gives *positive* log-likelihoods; one negative
  entry gives `NaN` throughout — and summing to one does not catch that, since −0.5 and +1.5 sum
  to one. Both are now checked, and the sum matters arithmetically rather than as bookkeeping
  because a truncated read takes a **suffix**: a seed summing to two at 3 in 100 contamination
  moved a genotype by 1.3 on the Phred scale with nothing failing.
- **`λ + c ≥ 1` was accepted and the floor fired silently**, flattening the row to a part in 10⁹.
  Spec §4.5.1 asks for the floor and it stays, but as a guard on the arithmetic rather than a
  policy: the two shares come from different fits and neither knows about the other, so a sum past
  one is a fit that has gone wrong and is now refused aloud.
- **The support's ascending order was relied on and unchecked** — a reversed slice returns 0.0
  where the sorted one returns 0.335, with no panic.
- **The fraction table's width was checked per observation**, contradicting the row's own
  documented promise that everything is checked once up front.
- **The headline per-sample test tested nothing.** It called a pure function twice with identical
  arguments and killed 0 of 17 mutations. It is replaced by one that pins the junk floor's **size**
  — 0.01/31 from the locus, against the 0.005 production would have given this two-sequence sample.

### What the independent oracle confirmed

A reviewer transcribed spec §4.5.1 into an oracle sharing no code with the row and compared
**4,128 cells** — three motifs × four candidate sets × one and three read groups × ploidy 1/2/4 ×
four contamination fractions × seed on and off, every case carrying a read at an unreachable length
and a read that ran out. **Worst disagreement: zero units in the last place.** Five edge cases the
row's own tests skip came back the same.

## 4. What the repair is actually worth, measured

**The junk floor's swing falls from tenfold to about 2.2–2.6 fold**, not to nothing. What is
removed is the dependence on *what samples showed*; the candidate set is itself cohort-derived, so
a locus admitting one more candidate reaches more lengths — one extra candidate at a five-candidate
dinucleotide locus moves the floor by 5 parts in 100. That is the honest size: much smaller, and
no longer a function of the reads.

## 5. Validation

Dev container, rustc 1.98:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features --tests -- -D warnings` — clean.
- `cargo test --lib` — **4,532 passed, 0 failed, 14 ignored**; **219 in `ng::calling::likelihood`**,
  of which **24 are the row's**.

**Six mutations re-run against the repaired tests, all killed**: folding `c` into `λ`; both seed
boundary shifts; every read group charged the first one's fraction; and the junk term returned to a
point mass. Each restore was checksum-verified before the next ran.

## 6. What the plan leaves open

- **Nothing converts the prior's per-candidate seed into the per-length one the row takes.** That
  is the calling loop's, and no test exercises the conversion because it does not exist yet.
- **Whether a reachable length no candidate spells should get zero or the geometric evaluated
  there** is the live half of the architecture's open question — see §2.
- **`emission` scores 22 length changes the support calls unreachable** (G1's report), unchanged.
- **Spec §12's items 14–19** are measurement runs needing genotypes end to end, recorded in
  [`calling_bakeoffs.md`](../../ng/impl_plan/calling_bakeoffs.md).
