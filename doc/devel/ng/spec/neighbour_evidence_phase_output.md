# ng — writing the phase the neighbour-evidence pass estimates

*Spec, 2026-11-02. No code yet — this settles the design. It extends
[`neighbour_evidence_phasing.md`](neighbour_evidence_phasing.md) (the pass's phased templates),
which estimates each sample's two haplotypes and so far writes only genotypes, and
[`neighbour_evidence_pass.md`](neighbour_evidence_pass.md) (the pass). The owner set it as a goal
the same day, reversing both specs' non-goal "`GT` stays unphased", with one rule: **where the phase
is not confident it is not written — the output is blocks of confident phase, not a continuous
phasing.** The measurements behind it are `phase_accuracy.py` and `beagle_switch_errors.py` in
`doc/devel/ng/research/neighbour_evidence/`.*

---

## 1. What this is

After the pass's rounds every sample carries two haplotypes in every window: for each record, each
chromosome copy's probability of carrying each allele. They decide the written genotype but are
then thrown away. **This writes the phase they carry**: at a heterozygous genotype, which allele is
on which copy, as `GT` with `|` (`0|1` or `1|0`), and a `PS` field naming the block of records whose
phases belong together.

**Phase is written only where it is confident.** A heterozygous genotype whose phase is uncertain
is written unphased (`0/1`) with `PS` empty, and a link between two phased genotypes that is
uncertain ends one block and starts the next. Downstream tools read a block as "these alleles are
on one chromosome together" and nothing across blocks, so a short block costs little and a wrong
link inside a block costs a lot.

**What was measured** (Python prototype of the pass, founder simulator; the consensus of §2.1):
how often the phase flips wrongly between consecutive phased heterozygous records of a sample, and
the share of truly heterozygous records phased with confidence (`|h₁ − h₂| ≥ 0.5`), against Beagle
5.5 phasing calling's genotypes on the same cohorts:

| simulated cohort | the pass: wrong flips (true hets phased) | Beagle 5.5 |
|---|---|---|
| 50 samples, 10 at 30× + 40 at 3×, inbreeding 0.5, 16 founders | 1.5% (84%) | 6% (71%) |
| same, inbreeding 0.8 | 2.4% (83%) | 5% (68%) |
| same, inbreeding 0.8, 64 founders | 8.0% (54%) | 22% (66%) |
| 50 samples all at 3×, inbreeding 0.5 | 2.2% (82%) | 8% (65%) |
| 50 samples all at 3×, outbred | 6.7% (59%) | 21% (72%) |
| 200 samples, 40 at 30×, inbreeding 0.8, 64 founders | 3.5% (68%) | 5% (69%) |

These are flips counted over every phased record, with no block breaks: the blocks of §2.3 exist to
take the uncertain links out, so the flips *inside* blocks are to be lower — measured in §7. Every
simulated contig was one window, so joining windows (§2.4) is unmeasured. The prototype does not fit
each sample's copying lengths; the Rust pass does, and gave clearly better genotypes.

### 1.1 Goals

- Write `GT` phased and `PS` for every heterozygous genotype whose phase is confident, in blocks
  whose internal links are confident.
- Leave everything else as it is: an unphased heterozygote stays `a/b`, a homozygote stays `a/a`,
  `./.` stays `./.`.
- Join a sample's blocks across windows where the windows' overlap shows how, and start a new block
  where it does not.
- Deterministic and the same at any thread count, as the pass is.

### 1.2 Non-goals

- **A continuous phase** (owner): blocks may be short and many.
- **Phase of homozygotes.** A homozygote's phase carries no information; it stays `a/a` with `PS`
  empty.
- **Phase across `./.` or unphased records as such** — a block's records need not be consecutive;
  `PS` ties its heterozygotes together.
- **Read-backed phasing** — reads spanning two heterozygous sites. Out, as for the pass.
- **A per-record phase quality field.** `PS` blocks carry the decision; a quality is *Deferred*
  (§8).

### 1.3 Vocabulary

- **Haplotype** (as in the phasing spec): one copy's probability of carrying each allele, at each
  record.
- **Orientation** of a heterozygote `{a, b}`: whether copy 0 carries `a` and copy 1 `b`, or the
  reverse.
- **Phased record** of a sample: a heterozygous genotype written with `|`.
- **Link**: the relative orientation of two consecutive phased records of a sample — whether their
  alleles on copy 0 are on one chromosome together.
- **Block**: a run of a sample's phased records joined by confident links; written with one `PS`.

---

## 2. The method

### 2.1 A consensus over the rounds

Each round's haplotypes come from one drawn path, so a round's phase has random switches, and the
copy called 0 in one round may be copy 1 in the next — the walk is symmetric in its two copies. The
genotype averages rounds 2 to `R` freely (it does not depend on which copy is which); the phase
cannot be averaged until the rounds' copies are lined up.

**Each round's haplotypes are aligned to the consensus so far, then added to it.** For one sample
in one window, after round 2 the consensus is round 2's haplotypes. After each later round:

1. Over the records where both the consensus and the round *lean* — one copy towards some allele,
   the other copy away from it, `|h₀(a) − h₁(a)| ≥ 0.5` for some allele `a` — mark whether the round's
   orientation agrees with the consensus's.
2. Choose, over those records, a path of two states — the round kept or the round's copies swapped
   — that costs 1 at each record where the state disagrees and 2 at each change of state, the
   cheapest (Viterbi over two states; ties to "kept").
3. Each record takes the state of the last marked record at or before it (the first marked record's
   state before it). Swap the round's two copies where the state is "swapped", and add them to the
   consensus sums.

The consensus at a record is the sums over the rounds added, divided by their count. **Measured
with alignment to round 2 rather than to the running consensus** (`phase_accuracy.py`): it halved
the flips against the last round alone (3.5% → 1.5%, 4.5% → 2.2%, 6.4% → 3.5% on the cohorts of §1).
Aligning to the running consensus needs no second store of round 2 and is the natural streaming
form; the two are compared in §7 before the Rust build fixes one. *Soft, a starting value*: the costs
1 and 2.

### 2.2 Which heterozygotes are phased

A sample's record is phased when:

- the genotype written is heterozygous, `{a, b}` with `a ≠ b`, after the pass's allele subsetting
  (`line.rs`, `rewrite_line` renumbers the alleles it keeps; the orientation is carried through the
  renumbering); and
- **the consensus orients it confidently**: with `c₀` and `c₁` the consensus haplotypes,
  `w = c₀(a)·c₁(b)` and `w′ = c₀(b)·c₁(a)`, the larger of `w / (w + w′)` and `w′ / (w + w′)` is at
  least **`PHASE_MIN` = 0.9**. Copy 0 carries `a` when `w ≥ w′`.

`PHASE_MIN` = 0.9 is the prototype's `|h₀ − h₁| ≥ 0.5` made general: at `h₀ = 0.75`, `h₁ = 0.25` it
is `0.75² / (0.75² + 0.25²)` = 0.9. *Soft*: set from the prototype's rule, not tuned.

A record where the genotype is heterozygous but the consensus does not orient it is written
unphased: `a/b`, `PS` empty.

### 2.3 Blocks: which links are confident

**A link between two consecutive phased records of a sample is kept only when every averaged round
agrees on it** — each round, aligned (§2.1), orients the two records the same way relative to each
other as the consensus does. A link any round disagrees on ends the block; the next phased record
starts a new one. With `R` = 5 that is four rounds agreeing.

This is the cheapest test the rounds offer: a switch error inside one round's path shows as a
disagreement between rounds at exactly the link it crosses, and the consensus alone hides it.
**Whether four rounds catch enough** is measured (§7): the flips left inside blocks against the
flips at all links, and the blocks' lengths. *Soft*: if the flips inside blocks stay above about
1%, the next lever is more rounds for the phase (not for the genotype) or a stricter `PHASE_MIN`,
both measured before chosen.

A block also ends:

- where two consecutive phased records are farther apart than the sample's copying length `S`
  (pass §3.5): beyond `S` the model has no information linking them; and
- at a window's join where the join is not confident (§2.4), and at a contig's end.

### 2.4 Joining windows

Each window phases its samples from scratch, and its copy 0 has nothing to do with the last
window's copy 0. **Neighbouring windows overlap**: window `w + 1`'s left margin is the end of
window `w`'s core, already written. So for each sample the stage keeps, from window `w`, the
orientation it *wrote* at its phased records within the next window's left margin, and the block
they belong to. Window `w + 1` computes its own orientation at those records (its margin is phased
too); then:

- **at least 3 shared phased records, and the two windows' relative orientation the same at all of
  them** — window `w + 1`'s orientation is kept or swapped as a whole for that sample, to agree, and
  its first block continues `w`'s last;
- **otherwise** — window `w + 1`'s first block is new.

Swapping a sample's whole window is exact: the window's phase has no meaning except relative to
itself. *Soft*: 3 records and "all agree", unmeasured; §7 measures joins.

### 2.5 `PS`

A block's `PS` is **the `POS` of its first phased record** — the convention of WhatsHap and GATK's
read-backed phasing, so tools that read those read this. A block joined across windows keeps the
`PS` it started with. Two of a sample's blocks never share a `PS` (each starts at a different
record).

---

## 3. What the VCF says

- **`GT`**: a phased record `a|b` in copy order (copy 0's allele first); every other column as now.
- **`PS`**: a new `FORMAT` field, `Number=1,Type=Integer`, **appended after the existing fields**,
  so every field the pass already rewrites keeps its place; `.` where the column is not phased.
  Written on every site record the pass writes, so the `FORMAT` string is the same on every such
  record.
- **Header**: `##FORMAT=<ID=PS,...,Description="Phase set: the POS of the first phased record of
  this block; records of a sample with the same PS are phased together">`; `GT`'s description
  (`src/vcf/header.rs:610`, "Genotype, unphased") becomes "Genotype; phased (|) where PS is set,
  unphased (/) elsewhere" when the pass ran; `##neighbourEvidence=` gains `phase_min=0.9`.
- **`refine-genotypes`** writes the same, as it shares the stage.

---

## 4. Where it sits in the code

| what | existing code | change |
|---|---|---|
| the rounds and the store | `refine::refine_window_in_rounds`, `HaplotypeStore` | each sample keeps a consensus row (sums, `f32`, `2·(A − 1)` a record, as the store) and, at its core records, one orientation bit a round; the alignment of §2.1 runs inside the sample's round closure after its draw |
| a record's result | `RefinedRecord` | gains, per sample, the orientation (copy 0's allele of the pair) or none, and the block start (§2.5) |
| the window's join | `stage::NeighbourEvidenceStage` | keeps, per sample, the written orientation and `PS` at the last window's records inside the next window's margin; passes them into the solve; reads back the window's last |
| the line | `line::rewrite_line` | writes `|` and `PS` through the allele renumbering |
| the header | `vcf/header.rs` | `PS`, `GT`'s description, `phase_min` |

**Memory**, per sample and window: a consensus row the size of one store buffer (8 bytes a
biallelic record), and `R − 1` orientation bits a core record (a byte). At 3,000 samples and a
floor-sized window (phasing spec §5) that is about 0.9 GB more.

**Time**: the alignment is one pass over a sample's records a round, a few operations a record —
small beside a walk.

---

## 5. Traps

- **The copy labels are arbitrary in every round and every window** — never compare two rounds'
  or two windows' copy 0 without aligning them (§2.1, §2.4).
- **Allele subsetting renumbers alleles** after the genotype is chosen (`rewrite_line`): the
  orientation is "copy 0 carries allele `a` of the record as the pass read it" and is written in the
  renumbered alleles.
- **A record whose genotype is homozygous in the written `GT`** is never phased, even if the
  haplotypes lean apart there — the genotype is the averaged posterior's, the haplotypes are a
  consensus; they can disagree, and the `GT` wins.
- **A hidden record** (phasing spec §4.1) is written with a posterior like any other and is phased
  like any other.
- **The join reads what was written**, not what window `w` computed before writing — they are the
  same unless a later rule changes a written genotype, and reading the written one keeps it so.
- **Multi-allelic records**: `{a, b}` may be `1|2`; the lean test of §2.1 is "for some allele", and
  `w`, `w′` of §2.2 read the two alleles of the genotype.

---

## 6. When it runs

Whenever the pass runs, with no option (as phasing, ruled 2026-11-02).

---

## 7. How we know it works

1. **The prototype first** (before the Rust build): aligning to round 2 against aligning to the
   running consensus (§2.1); the flips inside blocks under §2.3's rule against all links, and the
   blocks' lengths (N50 and the share of phased records in blocks of at least 10); on the cohorts of
   §1. The rule is kept if the flips inside blocks are at most a third of the flips at all links.
2. **The Rust build against the truth** on the same cohorts and the 200-sample one: flips inside
   blocks, phased share, block lengths — against the prototype and Beagle.
3. **Joins**: a stage test with windows small enough that each contig holds several; the flips at
   joins against the flips inside windows; and the share of joins that continue a block.
4. **Unit tests**: the two-state alignment on a hand case (a round equal to the consensus with one
   stretch swapped is swapped back exactly there); `PHASE_MIN` at the boundary; a multi-allelic
   `1|2`; allele renumbering; a homozygote never phased; `PS` the first phased record's `POS`.
5. **Invariance**: the same bytes at 1, 4 and 30 threads.
6. **A real cohort**: the tomato test on kimura, where the deep accessions' read-backed phase (from
   WhatsHap, if run) or their consistency across overlapping windows can be compared.

---

## 8. Deferred, with a recommended home

- **A per-record phase quality** (`PQ`): the consensus's `w / (w + w′)` is one; whether users need it
  beside `PS` is for the owner once blocks exist. Home: a follow-up to this spec.
- **Phase from reads** spanning two heterozygotes, to join blocks the population cannot. Home: its
  own spec, if blocks prove short in deep samples.

---

## 9. Open questions — the owner's to rule, with leanings

- **Q1. `PHASE_MIN` = 0.9** (§2.2): the prototype's rule, generalised. Leaning: keep, and move it
  only on §7.1's measurement.
- **Q2. A link kept only when all averaged rounds agree** (§2.3). Leaning: keep; it is the one
  signal of a switch the rounds give for free.
- **Q3. A join needs 3 shared phased records, all agreeing** (§2.4). Leaning: keep; measured in
  §7.3.
- **Q4. `PS` = the `POS` of the block's first phased record** (§2.5). Leaning: keep — the WhatsHap
  and GATK convention.
