# ng — the neighbour-evidence pass: refining each sample's genotypes with the loci around them

*Spec, 2026-11-01; the owner ruled its four open decisions the same day (§9). It turns the evaluation in
[`../research/neighbour_evidence_hmm.md`](../research/neighbour_evidence_hmm.md) into a design, under
the owner's ruling recorded in its §11. It depends on `PL` and `GP` in every record
([`vcf_output.md`](vcf_output.md) §7.2, branch `pl-gp`) and is to be measured on the founder cohort of
[`str_stutter_benchmark.md`](str_stutter_benchmark.md) §12 (branch `cohort-founders`).*

---

## 1. What this is

pvc calls each locus on its own: a sample's genotype comes from its reads at that locus and a prior
from the cohort's allele frequencies. In a cohort of inbred lines the loci around it say more. The
lines descend from a few haplotypes, so a sample's chromosome is, stretch by stretch, a copy of
other samples'; and an inbred sample's two copies are one copy counted twice over long runs. The
pass uses both: for each sample, it asks along the chromosome *which other samples is this one
copying here, and are its two copies the same here?*, and re-calls each genotype from its own reads
together with the answer.

On tomato1 the prototype of this model removed 36–56% of genotype errors in a cohort at 1× and 3×,
and three in four true non-reference genotypes at sites where the sample had no read; in an outbred
simulated cohort it removed none (the report's §4–§5). This spec is that model, made a stage of pvc.

### 1.1 Goals

- Re-call every sample's genotype at every record pvc writes, SNPs, indels, multi-allelic records and
  STR records alike, from the sample's own likelihoods and the cohort's shared haplotypes.
- Give a genotype, with an honest quality, where a sample has no reads.
- Run on by default where it helps and off where it cannot, saying which and why.
- Work from what a record already carries — `PL`, `GP`, the alleles — so the same code can run on a
  finished VCF.

### 1.2 Non-goals

- **Phasing.** The pass infers which samples a sample copies, not which of its alleles sit on one
  chromosome; `GT` stays unphased.
- **New loci.** A locus pvc did not build stays unbuilt; the pass changes genotypes only. On one
  tomato accession alone at 3×, 42% of the true non-reference genotypes were at such loci (report
  §5.2); that loss belongs to candidate generation, not here.
- **Read linkage** — reads covering two sites. Left out by the owner's ruling (report §11).
- **Outbred cohorts.** Using shared haplotypes there needs phased templates and repeated sweeps over
  the cohort ([`07_decision`](../research/phasing_imputation/07_decision_separate_tool_and_route.md)
  §2); this pass does neither and is off there (§5).

### 1.3 Vocabulary

- **Template.** Another sample, as a source a copy can copy from. At each record a template is its
  *allele shares*: the expected fraction of its copies carrying each allele, from its `GP`. An inbred
  sample's shares are near 0 or 1, which is why a genotype can stand in for a haplotype.
- **Autozygous.** Of a stretch where a sample's two copies are one ancestral copy counted twice; the
  sample is homozygous along it.
- **Window.** A stretch of one contig whose records the pass holds and solves together (§2.2).

---

## 2. Where it sits

### 2.1 A stage over finished records

**The pass reads finished VCF lines and writes them rewritten.** That is the one form every run's
records share at the end: with the hidden-duplication filter on, records are parked on its spill as
finished lines and read back by its third pass; with it off they are encoded straight for the
writer. A stage in front of the writer that takes lines in genome order and gives lines back sits in
both, and the same code can run as its own command over a VCF a run already wrote (§9 Q4).

This is the "pass over the emitted records" [`run_streaming.md`](run_streaming.md) §4.3 requires for
any coupling between segments: nothing about calling changes, and the calling workers never wait
for it.

### 2.2 Windows

The pass cannot write a record until it has seen the records after it, because a sample's state at a
site depends on both sides. It holds records in **windows**: a *core* of `C` bases, written when the
window is solved, with *margins* of `M` bases either side that are read for context and written by
the neighbouring window. A record is written once the stream has passed the end of its window's
right margin. Records are written in the order they arrived.

**Defaults to be measured, starting at `C` = 2 Mb and `M` = 1 Mb.** The prototype's autozygous runs
and shared stretches were 10 to 300 kb on tomato, limited by tomato1's 100-kb regions; a margin of a
megabase holds several of either. On a run over short regions — a BED of 100-kb intervals — a
window is a region and its margins are empty; that is what the prototype measured.

A contig's end closes its last window. Nothing crosses contigs.

### 2.3 Parallelism

Within a window, every sample is solved on its own once the templates are fixed (§3.3), so the
samples of a window are shared among the run's threads. Windows are solved one after another, which
keeps the written order without reordering.

---

## 3. The model

### 3.1 In words

For one sample, walk along the window's records. At each, the hidden state says either *the two
copies are copying templates `k₁` and `k₂`*, or *the two copies are one, copying template `k`*.
Between neighbouring records a copy keeps its template with a probability that falls with the
distance between them, and the autozygous state persists the same way. At each record, a state is
scored by how well the sample's reads fit the alleles its templates carry. A forward and a backward
walk give, at every record, the probability of each state given the whole window; the genotype
posterior is the sample's own likelihood weighted by what the states say.

### 3.2 The inputs, read off the record

For a record with `A` written alleles and, at diploidy, `G = A(A+1)/2` genotypes:

- **The sample's likelihood** `L(g)` for each genotype, from `PL`: `L(g) ∝ 10^(−PL_g / 10)`.
- **Each template's allele shares** `q(a)`, from its `GP`: `q(a) = Σ_g GP_g · copies_a(g) / 2`.

A sample written `.` in both (set aside by the candidate step) is no template at that record and is
not re-called there. A record with no rows at all — one merged from an STR record and its flank's
generic records — is passed through unchanged and is not a site of any window's walk.

### 3.3 The templates

For each sample and window, the `K` other samples whose allele shares agree best with its own over
the window's records, measured as `−Σ_records Σ_a |q_self(a) − q_other(a)|` over the records where the
sample's `PL` is not flat (where its own reads said something). **The sample itself is never a
template.** `K` = 16 was what the prototype ran; how accuracy and cost move with `K` is a
measurement (§8).

**A known circularity, kept small and stated.** A template's `GP` came from a prior built from the
other samples' allele copies — the target's included. So a little of the target's own evidence
returns through its templates. With 56 others the prototype measured no sign of it; it grows as the
cohort shrinks, which is one reason the pass is off below a cohort size (§5).

### 3.4 The emission

A template's shares are blurred by a mismatch rate `ε` — a mutation since the shared ancestor, or the
template's own wrong call — `q̃(a) = (1 − ε)·q(a) + ε / A`, with `ε` = 0.01 as in the prototype. Then:

    not autozygous, templates (k₁, k₂):   e = Σ_a Σ_b q̃_k₁(a) · q̃_k₂(b) · L(a/b)
    autozygous, template k:               e = Σ_a q̃_k(a) · L(a/a)

The sums run over alleles, so a record with eight alleles, an STR among them, costs more arithmetic
and nothing else.

### 3.5 The transitions

Between records `d` bases apart:

- **A copy keeps its template** with probability `exp(−d / S)`, and otherwise draws one uniformly; in
  the non-autozygous state the two copies do this independently.
- **The autozygous state persists** with probability `exp(−d / R)`; the other state persists with
  `exp(−d / R′)`, `R′ = R · (1 − F) / F`, so a share `F` of the sample is autozygous in the long run.
  A change of state draws its templates afresh.

`S`, `R` and `F` are per sample. **`F` starts at the sample's fitted inbreeding coefficient** from the
parameters file and is re-estimated as the mean autozygous posterior; `S` and `R` are chosen by
likelihood from small grids (`S` ∈ {30, 300} kb, `R` ∈ {10 kb … 3 Mb}), as the prototype did. Refit per
window or once per contig is a cost-against-accuracy measurement (§8).

### 3.6 The posterior

At a record, with `P(state | window)` from the two walks and `e(state)` the record's emission:

    P(g | window) ∝ L(g) · Σ_states [ P(state | window) / e(state) ] · P(g | state)

where `P(a/b | k₁, k₂) = q̃_k₁(a) q̃_k₂(b) + q̃_k₁(b) q̃_k₂(a)` for `a ≠ b` (and the one term for `a = b`),
and `P(a/a | autozygous k) = q̃_k(a)`. Dividing by the emission removes the record's own reads from the
state probability, so they enter once, through `L(g)`.

---

## 4. What it writes

| field | after the pass |
|---|---|
| `GT` | the argmax of the new posterior |
| `GQ` | from the new posterior, as [`calling_quality.md`](calling_quality.md) §4 defines it |
| `GP` | the new posterior |
| `PL` | unchanged — it is the evidence |
| `DP`, `AD`, `REPCN`, `REPCNB` | unchanged (`REPCN`/`REPCNB` follow the new `GT`) |
| `INFO` `AC`, `AN` | recounted from the new genotypes |
| `QUAL`, `AF`, every other `INFO` field, `FILTER` | unchanged |

**A sample whose reads said nothing** — flat `PL`, written `./.` by calling (§7.1 there) — gets the
pass's posterior, and **is written with a genotype when the pass gives one**: the evidence is no
longer absent, it is the neighbouring loci's. Its `GQ` says how sure that is, and its flat `PL` says
none of it came from its own reads. *This reverses §7.1's rule for the pass's output; ruled by the owner, 2026-11-01 (§9 Q1).*

**The header** gains one line recording whether the pass ran and why, and its parameters:

    ##neighbourEvidence=<ran=yes,why="57 samples, median inbreeding 0.98",templates=16,mismatch=0.01,core=2000000,margin=1000000>

---

## 5. When it runs

**On by default when the cohort can supply templates and is inbred enough for genotypes to stand in
for haplotypes; off otherwise, saying so in the header and the run report.** The two conditions:

- **At least `N_min` samples with rows at the record's alleles.** Starting value 10: the simulation
  showed most of the gain at 10 inbred samples (errors 796 → 432 at 3×, report §4).
- **The cohort's median fitted inbreeding coefficient at least `F_min`.** Starting value 0.5. The
  outbred simulation gained nothing; between the two ends is unmeasured.

Both are starting values, to be set from §8's measurements. A command-line option overrides the
decision either way — `--neighbour-evidence auto|on|off`, `auto` the default — because a user who
knows their material is a selfing panel at F = 0.4 should be able to say so.

**At a ploidy other than two the pass is off**: the pair state is diploid.

---

## 6. One sample and three thousand, three reads and three hundred

- **One sample.** Off: there are no templates. The report measured what a lone sample would gain from
  the autozygous state alone — 5% at 3× on tomato — and the owner ruled that not a goal.
- **Three thousand samples.** Two costs grow with the cohort.
  - *Choosing templates* compares each sample with every other over a window: `N² × records` a window.
    At 3,000 samples and the 52,000 records of a 2-Mb tomato window that is about 5 × 10¹¹ operations
    a window, too many. **The first version compares on a fixed random sample of the window's
    informative records (starting at 500)**, about 5 × 10⁹ a window; at that scale a positional
    Burrows–Wheeler index is the tools' answer (07 §2.6) and is deferred (§10).
  - *The walks* cost `K² + K` per sample per record per pass: 272 at `K` = 16, so a genome of tomato's
    density (about 24 million records over 900 Mb, own arithmetic from tomato1's 26 a kilobase) costs
    about 6.5 × 10⁹ per sample per pass. Over 3,000 samples and the grid's passes that is the
    dominant cost of a run; §8 measures it and how far `K` can fall.
- **Three reads a position.** Where it helps: the prototype's gains are all at 1× and 3×.
- **Three hundred reads.** The sample's own likelihood dominates every emission and the posterior is
  the calling's; the pass costs time and changes nearly nothing. It still runs if the switch says so:
  depth does not turn it off, because a deep cohort can hold shallow samples.

**Memory.** A window holds its records' rows: at 57 samples and 52,000 records of 2 Mb plus
margins, about 100,000 records × 57 samples × about 6 numbers, 270 MB at 8 bytes; at 3,000 samples
14 GB. **The pass therefore holds rows as `f32` and windows shrink with the cohort** —
`C + 2M` chosen so `records × samples` stays under a budget (starting at 200 million rows, 4.8 GB);
the margin is kept, the core shortened. Both are measured (§8).

---

## 7. Traps

- **Genotype order.** `PL` and `GP` are in the shared table's order (`vcf_output.md` §7.2); the pass
  must read and write them through the same table, never by a hand-rolled index.
- **A record whose alleles the pass dropped nothing from but whose `GT` moved** must have `REPCN` and
  `REPCNB` recomputed from the new `GT`, and `AC`/`AN` recounted; a stale `REPCN` is a silent error.
- **A flat `PL` is not "no information" for the walk:** it contributes an emission of 1 to every
  state, which is correct, and must not be skipped, or the distances between the remaining records
  would stretch the transitions.
- **The window's right margin is not written by that window.** Writing a record twice, or not at all,
  at a window boundary is the one bookkeeping error the oracle in §8 checks for directly.

---

## 8. How we know it works

1. **Against the Python reference.** The prototype's `hmm.py` on the same input lines gives the same
   posteriors to 1e-6, on a tomato window and on simulated windows with multi-allelic and STR records.
2. **On the founder cohort** (`str_stutter_benchmark.md` §12), scored against `truth.vcf` and
   `ancestry.tsv`: founders 4, 16, 64; F 0, 0.5, 0.8, 0.95; samples 10, 50, 200; depth 1, 3, 10, 30.
   What it sets: `N_min`, `F_min`, `K`, `C` and `M`, whether `S` and `R` are refitted per window.
   The autozygous posterior is scored against the runs in `ancestry.tsv`.
3. **On tomato1 held out**, as the report measured: the nine deep accessions cut to 3× and 1×, their
   full-depth calls as truth. The prototype's numbers (report §5.2) are the floor the Rust version has
   to reach.
4. **Cost** at 57 and 500 samples on tomato1's regions, and extrapolated to 3,000 by §6's arithmetic
   checked against the 500.
5. **Invariants:** every record in is a record out, in order; the pass with every sample's `PL` flat
   and `K` = 0 returns its input's genotypes; a record without rows comes out byte-identical.

---

## 9. Decisions — **ruled by the owner, 2026-11-01**

All four as recommended.

1. **A sample whose reads said nothing: write the pass's genotype, or keep `./.`?** **Ruled: write
   it.** A call there is the aim the owner set ("where possible a call at a site with little or
   no evidence"), its `GQ` is honest about how sure it is, and its flat `PL` lets anyone filter it
   out. Keeping `./.` would throw away the three in four the prototype recovered.
2. **Mark the refined genotypes per sample?** **Ruled: no new per-sample field.** The header says
   the pass ran; `PL` against `GP` shows what it changed; a field per sample would add bytes to every
   column of every record for a fact one header line states.
3. **`--neighbour-evidence auto|on|off`, with `auto` deciding by §5's two conditions?** **Ruled:
   yes.**
4. **A stand-alone command over a finished VCF**, the same stage run on a file (`pvc refine-genotypes
   --vcf in.vcf.gz --output out.vcf.gz`)? **Ruled: yes, second, after the in-run stage** — it costs a
   thin command around the same code, and it lets a cohort called in parts, or called before the pass
   existed, be refined once its `PL` and `GP` are written.

---

## 10. Deferred, with a recommended home

- **Template choice at thousands of samples by a positional Burrows–Wheeler index** (07 §2.6) — here,
  once §8 measures the sampled-record choice at 500 samples and it does not scale.
- **Phased templates and sweeps over the cohort**, the outbred case — not this pass; 07 §2.
- **Ploidy above two** — a state per multiset of templates; here, when a polyploid cohort is a target.
- **Refining with the pass's own output** — templates from the refined genotypes, a second sweep —
  here, if §8 shows the first sweep leaves accuracy on the table.
