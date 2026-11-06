# ng — writing the phase the neighbour-evidence pass estimates

*Spec, 2026-11-02; revised the same day after an adversarial review
([`neighbour_evidence_phase_output_spec_review_2026-11-02.md`](../../reports/reviews/neighbour_evidence_phase_output_spec_review_2026-11-02.md)).
No code yet — this settles the design. It extends
[`neighbour_evidence_phasing.md`](neighbour_evidence_phasing.md) (the pass's phased templates),
which estimates each sample's two haplotypes and so far writes only genotypes, and
[`neighbour_evidence_pass.md`](neighbour_evidence_pass.md) (the pass). The owner set it as a goal
the same day, reversing both specs' non-goal "`GT` stays unphased", with one rule: **where the phase
is not confident it is not written — the output is blocks of confident phase, not a continuous
phasing.** The measurements behind it are `phase_accuracy.py` and `beagle_switch_errors.py` in
`doc/devel/ng/research/neighbour_evidence/`, and the review's probes (`blocks_probe.py` and its
outputs, in `tmp/phase_output_review/` of the `neighbour-pass` worktree).*

**What the review changed.** The first draft's block rule compared rounds *after* aligning them, and
the alignment hides exactly the switches the rule was meant to catch; the rule now compares each
round's own, unaligned phase (§2.3). The confidence threshold is now the prototype's measured rule,
not a looser one that only resembled it (§2.2). The join of windows is defined per sample over the
last written block (§2.4). Repeat counts follow the phased `GT` (§3), an already-phased input is
refused (§3), multi-allelic records keep their rounds' haplotypes (§2.1), and §7 is rebuilt so it
can decide §9. The owner then ruled §9: the stricter link rule, and no cohort gate — one ceiling on
flips inside blocks for every cohort instead (§6, §7.1).

---

## 1. What this is

After the pass's rounds every sample carries two haplotypes in every window: for each record, each
chromosome copy's probability of carrying each allele. They decide the written genotype but are
then thrown away. **This writes the phase they carry**: at a heterozygous genotype, which allele is
on which copy, as `GT` with `|` (`0|1` or `1|0`), and a `PS` field naming the block of records whose
phases belong together.

**Phase is written only where it is confident.** A heterozygous genotype whose phase is uncertain
is written unphased (`0/1`) with `PS` empty, and a link between two phased genotypes that is
uncertain ends one block and starts the next. Tools that read `PS` (WhatsHap; SHAPEIT4 with
`--use-PS`) read a block as "these alleles are on one chromosome together" and nothing across
blocks, so a short block costs little and a wrong link inside a block costs a lot. Tools that
ignore `PS` (`bcftools consensus -H`, for one) take every `|` at face value; for them a
single-record block would be a random orientation, which is why such records are written unphased
(§2.2).

**What was measured** (Python prototype of the pass, founder simulator; the consensus of §2.1, with
no block rule): how often the phase flips wrongly between consecutive phased heterozygous records of
a sample — a *flip* — and the share of truly heterozygous records phased with confidence
(`|h₀ − h₁| ≥ 0.5`), against Beagle 5.5 phasing calling's genotypes on the same cohorts:

| simulated cohort | the pass: flips (true hets phased) | Beagle 5.5 |
|---|---|---|
| 50 samples, 10 at 30× + 40 at 3×, inbreeding 0.5, 16 founders | 1.5% (84%) | 6% (71%) |
| same, inbreeding 0.8 | 2.4% (83%) | 5% (68%) |
| same, inbreeding 0.8, 64 founders | 8.0% (54%) | 22% (66%) |
| 50 samples all at 3×, inbreeding 0.5 | 2.2% (82%) | 8% (65%) |
| 50 samples all at 3×, outbred | 6.7% (59%) | 21% (72%) |
| 200 samples, 40 at 30×, inbreeding 0.8, 64 founders | 3.5% (68%) | 5% (69%) |

The comparison is not like for like: the pass's flips are counted over the records it chose to
phase, Beagle's over every heterozygote it called. Both are scored at truly heterozygous records.
Every simulated contig was one window, so joining windows (§2.4) is unmeasured. The prototype phased
biallelic SNPs only and did not fit each sample's copying lengths; the Rust pass does fit them, and
gave clearly better genotypes.

**The blocks of §2.3 lower the flips but do not remove them.** Measured by the review on the same
rounds, the flips left inside blocks are 0.54% on the first cohort and 5.4% on the 64-founder one
(§2.3). The rounds share a switch error that every round makes together a third to two thirds of the
time, and no comparison between rounds can see it. §6 and §7.1 deal with what that means on diverse
cohorts.

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
- **Keeping a phase the input already carries** (§3).

### 1.3 Vocabulary

- **Haplotype** (as in the phasing spec): one copy's probability of carrying each allele, at each
  record.
- **Orientation** of a heterozygote `{a, b}`: whether copy 0 carries `a` and copy 1 `b`, or the
  reverse.
- **Lean** of a copy pair at a record, for the pair `{a, b}`: `w − w′`, with `w = h₀(a)·h₁(b)` and
  `w′ = h₀(b)·h₁(a)`. Its sign is the orientation; its size is how sure. For a biallelic record it
  is exactly the difference of the two copies' alternative probabilities, `h₁ − h₀` for `a` = REF.
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

**Each round's haplotypes are aligned to the running consensus, then added to it.** For one sample
in one window, after round 2 the consensus is round 2's haplotypes. After each later round:

1. Over the records where both the consensus and the round *lean* by at least 0.5 — for a biallelic
   record, `|h₀ − h₁| ≥ 0.5` in alternative probability; for a multi-allelic record, on the allele
   `a` with the largest `|c₀(a) − c₁(a)|` in the consensus — mark whether the round's orientation
   agrees with the consensus's (the same sign of `h₀(a) − h₁(a)`).
2. Choose, over those records, a path of two states — the round kept or the round's copies swapped
   — that costs 1 at each record where the state disagrees and 2 at each change of state, the
   cheapest (Viterbi over two states; ties to "kept").
3. Each record takes the state of the last marked record at or before it (the first marked record's
   state before it). Swap the round's two copies where the state is "swapped", and add them to the
   consensus sums.

The consensus at a record is the sums over the rounds added, divided by their count. **Measured**
(`phase_accuracy.py`): the consensus halved the flips against the last round alone (3.5% → 1.5%,
4.5% → 2.2%, 6.4% → 3.5% on the cohorts of §1). Aligning to the running consensus or to round 2
makes no difference (1.49% against 1.44%, 2.40% against 2.33%, 7.79% against 7.77% on the three
mixed-depth cohorts; review m4); the running consensus is kept because it needs no second store.
*Soft, a starting value*: the costs 1 and 2. A swapped stretch is undone exactly only if it holds at
least five marked records inside the window, or three at either end (two changes of state cost 4).

**The alignment serves the consensus only.** The block rule (§2.3) does not use it.

**What each round also keeps, per sample, at every record of the window** (core and left margin —
the join of §2.4 reads the margin): for a biallelic record, two bits — the sign of the round's own
lean, *before* alignment, and whether it leans by at least 0.5. For a record with more than two
alleles, the round's own haplotypes, because the pair `{a, b}` whose orientation matters is known
only after the last round has set the genotype; such records are few.

### 2.2 Which heterozygotes are phased

A sample's record is phased when all of these hold:

- the pass called it (outcome *Called*). A record the pass left *Unchanged* — a sample with no
  informative record in the window, which keeps its `GP` shares as both haplotypes
  (`refine.rs:346–375`) — has `w = w′` and is never phased; the same holds for *Withheld*;
- the genotype written is heterozygous, `{a, b}` with `a ≠ b`, after the pass's allele subsetting
  (`line.rs`, `rewrite_line` renumbers the alleles it keeps; the orientation is carried through the
  renumbering);
- **the consensus leans on it by at least `LEAN_MIN` = 0.5**: `|w − w′| ≥ 0.5` with `c₀`, `c₁` the
  consensus haplotypes. Copy 0 carries `a` when `w > w′`. This is the prototype's measured rule
  exactly (for a biallelic record `|w − w′| = |c₀ − c₁|`), and it reads a multi-allelic genotype's
  two alleles;
- it is not a block of one record (§2.3): a phased record with no confident link to either
  neighbour is written unphased.

The first draft used `w / (w + w′) ≥ 0.9` and called it a generalisation of the prototype's rule. It
is strictly looser: `w / (w + w′)` is the probability of one orientation *given* a heterozygote, and
is high when both orientations are unlikely. On the first cohort it phased 4,425 truly homozygous
records against 100 (review I1), which matters wherever the written heterozygote is wrong — sharpest
at hidden records, re-estimated without the sample's reads.

A heterozygous genotype filled at a record without reads (pass, no-read rules) is phased like any
other. Its phase was measured no worse: links touching a fill flipped 4 times in 456 (0.9%) against
91 in 5,777 (1.6%) between records with reads, on the first cohort (review I7). Its genotype is held
to `GQ` 10, as every fill.

### 2.3 Blocks: which links are confident

**A link between two consecutive phased records of a sample is kept only when, in every averaged
round, the round's own relative orientation of the two records is the consensus's** — copy 0
carrying the same allele of each genotype at both records, or not. The round's relative orientation
does not depend on which copy it calls 0, so it needs no alignment: for biallelic records it is the
exclusive-or of the round's two stored sign bits (§2.1). A link any round disagrees on ends the
block; the next phased record starts a new one. With `R` = 5 that is four rounds agreeing.

**Why not the aligned rounds.** The first draft compared the rounds after §2.1's alignment. The
alignment puts its change of state exactly at a round's switch error, so after alignment the round
agrees with the consensus on both sides and the test sees nothing; with four rounds it can never
fire between records where the consensus leans by more than 0.75. It left 62% to 75% of the flips;
the unaligned test leaves 33% to 63% (review B1):

| cohort | flips, no block rule | flips inside blocks, aligned test | flips inside blocks, unaligned test | block N50 in records: none / aligned / unaligned |
|---|---|---|---|---|
| 10 at 30× + 40 at 3×, inbreeding 0.5, 16 founders | 93 of 6,228 (1.49%) | 65 of 5,837 (1.11%) | 31 of 5,753 (0.54%) | 128 / 51 / 39 |
| same, inbreeding 0.8 | 65 of 2,704 (2.40%) | 40 of 2,539 (1.58%) | 34 of 2,524 (1.35%) | 77 / 43 / 35 |
| same, inbreeding 0.8, 64 founders | 141 of 1,810 (7.79%) | 94 of 1,624 (5.79%) | 87 of 1,605 (5.42%) | 53 / 24 / 23 |
| 50 at 3×, inbreeding 0.5 | 131 of 6,093 (2.15%) | 84 of 5,632 (1.49%) | 58 of 5,563 (1.04%) | 127 / 41 / 34 |
| 50 at 3×, outbred | 635 of 9,796 (6.48%) | 408 of 7,802 (5.23%) | 322 of 7,601 (4.24%) | 107 / 14 / 11 |
| 200, 40 at 30×, inbreeding 0.8, 64 founders | 350 of 10,221 (3.42%) | 263 of 9,518 (2.76%) | 219 of 9,426 (2.32%) | 73 / 42 / 37 |

(Prototype rounds, running consensus, records phased at `|c₀ − c₁| ≥ 0.5`, no single-record or `S`
rule.)

**What no round comparison can see.** A switch that all four rounds make together leaves them in
agreement. The rounds are not independent — each copies the haplotypes the round before drew — and
such shared switches are a third to two thirds of the flips (31 of 93 on the first cohort, 87 of 141
on the 64-founder one). More rounds or a stricter link rule shorten the blocks without reaching
them.

**The rule adopted (owner, 2026-11-02, "for now, we'll see how it works") is the stricter one:
every averaged round also leans by at least 0.5 — its own, unaligned lean — at both records of the
link.** Measured by the review: first cohort 16 flips in 4,827 links (0.33%) at a block N50 of 16
records, against 31 (0.54%) at 39 with agreement alone; 64-founder 48 in 1,361 (3.5%). With 9 rounds
for the phase, 64-founder 23 in 1,149 (2.0%) at a block N50 of 8.

**Then made stricter by the prototype check (§7.1;
[`phase_blocks_check_2026-11-02.md`](../research/neighbour_evidence/phase_blocks_check_2026-11-02.md)):
each round draws 8 more paths from the same forward probabilities, and the test runs over all 36
paths (4 rounds × 9).** A link is kept only when every path gives the two records the consensus's
relative orientation and leans by at least 0.5 at both. The extra paths serve the link test only;
the round's own path still feeds the next round and the consensus, so genotypes are unchanged. The
rule passes the 1% ceiling on every inbred and mixed-depth cohort: 0.05% to 0.35% flips inside
blocks, with 28% to 61% of true heterozygotes phased. The round rule above passes on only two of
them. It fails on outbred cohorts: 1.2%, 7.8% and 36% (§6). Fewer extra draws are not enough: with
4, the 64-founder 50-sample cohort is at 1.43%.

As built, each extra path is a backward draw of its own, from the forward probabilities the
round already holds; drawing them from shared cumulative weights a record would be cheaper, but
the whole phase added 10% at 200 samples (§4), so it was not needed. The link test needs each path's orientation bit at every record whose paths all lean,
since only those can be linked: 36 bits, about 5 bytes a record a sample. (As built, 16 bytes at
every allele pair — §4.)

A stricter `LEAN_MIN` is **not** a lever on the block rule — above 0.75 no round can disagree — but
it is one on which records are phased.

A block also ends:

- where two consecutive phased records are farther apart than the sample's copying length `S`
  (pass §3.5). Beyond `S` a copy still keeps its template with probability `e⁻¹` ≈ 0.37, so the
  model links the records only weakly; measured, links longer than 300 kb flipped 23 times in 49
  over the six cohorts, about chance. `S` is fitted per sample from {30, 300} kb; the prototype ran
  at 300 kb only, so the break at 30 kb is unmeasured;
- at a window's join where the join is not confident (§2.4), and at a contig's end;
- between two phased records **at one `POS`**: there the link is not broken but the second record is
  written unphased, so two blocks never share a `PS` (§2.5). No two PASS records shared a `POS` in
  the first cohort's refined VCF (0 of 2,728) — the pass merges overlapping records — so this guards
  a case that may never occur.

### 2.4 Joining windows

Each window phases its samples from scratch, and its copy 0 has nothing to do with the last
window's copy 0. Neighbouring windows overlap: window `w + 1`'s left margin covers records an
earlier window has already written. That margin need not lie inside `w`'s core alone — the planner
shrinks a core down to the margin's length (`windows.rs`, `size_ladder`), so after a short core the
next margin can reach into `w − 1`'s core.

**The join state is kept per sample, not per window: the sample's last written block** — its `PS`,
and the orientation written at its phased records within the last `margin_bases` bases before the
next window's core. The stage keeps it across windows that write no phase for the sample (no site in
the core, or the sample not walked) and clears it at a contig's end.

Window `w + 1` phases its left margin like its core (its consensus and blocks run over the whole
window). For each sample, take the records of the last written block that `w + 1` also phases, and
of those the **`k` = 3 nearest the core boundary**. Then:

- **if there are `k` such records, they lie in one block of `w + 1` that runs on, with confident
  links, to `w + 1`'s first phased core record, and the two windows' relative orientation is the same
  at all `k`** — `w + 1`'s orientation is kept or swapped as a whole for that sample, to agree, and
  that block continues the last written block (keeps its `PS`);
- **otherwise** — `w + 1`'s first block in the core is new.

Swapping a sample's whole window is exact: the window's phase has no meaning except relative to
itself. Comparing only the `k` records nearest the boundary, inside one block of each side, is what
makes a join possible at all: the first draft compared the whole margin, about 62 links a side on
the first cohort at the default 1-Mb margin, which at 1.5% flips a link joins about one time in six
(0.985¹²⁴ ≈ 0.15, if flips are independent), and it crossed `w`'s own block breaks. *Soft*: `k` = 3,
unmeasured; §7.3 measures joins.

### 2.5 `PS`

A block's `PS` is **the `POS` of its first phased record** — the convention the VCF specification
recommends ("the position of the first variant in the set", VCF 4.4, genotype fields), and the one
WhatsHap writes. A block joined across windows keeps the `PS` it started with. Two of a sample's
blocks never share a `PS`: each starts at a different record, and §2.3 keeps a block from starting
at the `POS` of the record before it.

---

## 3. What the VCF says

- **`GT`**: a phased record `a|b` in copy order (copy 0's allele first); every other column as now.
- **`REPCN` and `REPCNB`** are declared "GT order" (`header.rs:616–617`) but are written today in
  ascending allele order (`line.rs`, `parse_genotype` and `rewrite_line`), which agreed while `GT`
  was ascending too. **They are written in the order of the written `GT`**, so a `1|0` STR record
  gives each haplotype its own copy number.
- **`PS`**: a new `FORMAT` field, `Number=1,Type=Integer`, **appended after the existing fields**,
  so every field the pass already rewrites keeps its place; `.` where the column is not phased.
  Written on every site record the pass writes, so the `FORMAT` string is the same on every such
  record.
- **Header**: `##FORMAT=<ID=PS,...,Description="Phase set: the POS of the first phased record of
  this block; records of a sample with the same PS are phased together">`; `GT`'s description
  (`src/vcf/header.rs:610`, "Genotype, unphased") becomes "Genotype; phased (|) where PS is set,
  unphased (/) elsewhere" when the phase is written. `##neighbourEvidence=` — built by
  `NeighbourEvidenceDecision::as_header_value` (`stage.rs:66–79`) — gains `lean_min=0.5`.
- *The two items below described `pvc refine-genotypes`, removed on 2026-11-03 (read phase spec
  §5); they are kept as the record of what was built.*
- **`refine-genotypes`** writes the same, as it shares the stage, but builds its own header from the
  input's lines (`refine_genotypes.rs:164–200`), so the `GT` and `PS` declarations are added there
  too.
- **An input already phased is refused** by `refine-genotypes`: a header that declares `PS`, or a
  `|` in any `GT` it reads, stops it with a message naming the field or the record (strip the phase
  first, e.g. `bcftools annotate -x FORMAT/PS`). Otherwise a `FORMAT` would carry `PS` twice, and
  the columns the pass leaves *Unchanged* would keep an input block whose `PS` follows the same
  convention and can collide with a pass block. pvc's calling never writes `|`, so a run is not
  affected.

---

## 4. Where it sits in the code

| what | existing code | change |
|---|---|---|
| the rounds and the store | `refine::refine_window_in_rounds`, `phase::SamplePhase` | each sample keeps a consensus row (sums, `f32`, `2·(A − 1)` a record, as the store) and, at every allele pair of every record, each path's two bits — its lean's sign and whether it reaches 0.5 — as two 64-bit masks; the extra paths come from `model::walk_and_draw_paths`; the alignment of §2.1 runs inside the sample's round closure after its draw (built: steps 1–2, `1ab612e2`, `7686b276`) |
| the blocks | new, after the last round | the phased records (§2.2), the links (§2.3), the window's blocks per sample |
| a record's result | `RefinedRecord` | gains, per sample, the orientation (copy 0's allele of the pair) or none, and the block's `PS` |
| the window's join | `stage::NeighbourEvidenceStage` | keeps, per sample, the last written block (§2.4) across windows; passes it into the solve; updates it from what the window writes |
| the line | `line::rewrite_line` | writes `|`, `PS`, and the repeat counts in `GT` order, through the allele renumbering |
| the header | `vcf/header.rs`, `stage.rs` (`as_header_value`), `cli/refine_genotypes.rs` | `PS`, `GT`'s description, `lean_min`; the refusal of a phased input |

**Memory**, per sample and window, as built: a consensus row the size of one store buffer (8 bytes
a biallelic record), and the paths' bits — two 64-bit masks, 16 bytes, at every allele pair of
every record (one pair at a biallelic record, three at a triallelic). At 3,000 samples and a
floor-sized window of 39,000 records (phasing spec §5) that is 39,000 × 3,000 × 8 bytes = 0.94 GB for
the consensus and 1.9 GB for the bits — beside the store's two buffers, 1.9 GB. The bits could be
packed into 9 bytes (36 paths × 2 bits) or kept only where every path leans; not done, since no
cohort that size has been run and the owner prefers the simpler code until one is.

**Time**, measured (`refine-genotypes`, 6 threads, this laptop): the extra paths and the phase add
10% at 200 samples (15.9 s → 17.5 s); on the 50-sample cohorts, single runs of 3.2 s to 4.0 s,
from 6% faster to 21% slower — differences of the size of one run's noise. The extra paths are drawn one by one, each a backward pass; drawing them from shared
cumulative weights is not needed at this cost.

---

## 5. Traps

- **The copy labels are arbitrary in every round and every window** — never compare two rounds'
  or two windows' copy 0 without aligning them (§2.1, §2.4). The link test of §2.3 compares relative
  orientations, which have no labels.
- **Do not test links on the aligned rounds** (§2.3): the alignment absorbs the switches.
- **Allele subsetting renumbers alleles** after the genotype is chosen (`rewrite_line`): the
  orientation is "copy 0 carries allele `a` of the record as the pass read it" and is written in the
  renumbered alleles; the repeat counts follow it.
- **A record whose genotype is homozygous in the written `GT`** is never phased, even if the
  haplotypes lean apart there — the genotype is the averaged posterior's, the haplotypes are a
  consensus; they can disagree, and the `GT` wins.
- **A record the pass did not call** (*Unchanged*, *Withheld*) is never phased (§2.2).
- **A hidden record** (phasing spec §4.1) is written with a posterior like any other and is phased
  like any other; `LEAN_MIN` on `|w − w′|` keeps one whose copies both lean to one allele unphased.
- **The join reads what was written**, not what window `w` computed before writing — they are the
  same unless a later rule changes a written genotype, and reading the written one keeps it so.
- **Multi-allelic records**: `{a, b}` may be `1|2`; the alignment's lean is on the consensus's most
  divergent allele (§2.1), and the phase on the genotype's own pair (§2.2) from the kept round
  haplotypes.

---

## 6. When it runs

Whenever the pass runs, with no option and **no gate on the cohort** (owner, 2026-11-02): pvc is a
general tool, and the phase is written only where it is confident, so the confidence rule — not a
judgement about the cohort — decides how much is written. A cohort that phases badly must come out
with few, short blocks, not with wrong ones.

**This is a requirement on the rule, and the rule does not meet it yet.** Round agreement measures
how much the draw varies, not whether the model is right, and a switch all rounds share passes it
(§2.3). Measured, links the rule calls confident are wrong 0.33% of the time on the first cohort and
3.5% on the 64-founder one; on 10 outbred samples at 3× (64 founders) 34 of 85 links inside blocks
flipped (40%) under the plain agreement test — every round agreed on a near-random phase. So
"confident" does not yet mean the same thing on every cohort. §7.1 sets one ceiling for all cohorts,
and the Rust build waits until the rule meets it on each.

**Measured with the 36-path rule (§2.3), the outbred cohorts still fail.** At 50 outbred samples
from 64 founders, 1,977 links have all 36 paths agreeing, and 244 of them (12%) are wrong. The
draws sample the model's uncertainty given its templates. When the templates do not hold the
sample's haplotypes, the model is confidently wrong, and no number of draws shows it. Fixing that
needs a check that does not come from the model's own draws (open, §9).

---

## 7. How we know it works

1. **The prototype first** (before the Rust build), on the six cohorts of §2.3's table plus 10
   samples outbred and 10 samples at inbreeding 0.8 (each at 3×, from 64 founders), and an
   all-outbred 50-sample cohort from 64 founders:
   - the link rule: unaligned test, + every round leaning, + 9 rounds — flips inside blocks, block
     N50 in records and in bases, share of phased records in blocks of at least 10;
   - `LEAN_MIN` at 0.5, 0.6 and 0.8: share of true heterozygotes phased and flips inside blocks;
   - joins: run `phaser.py` on two overlapping halves of a contig and apply §2.4 — the share of joins
     that continue a block, and the flips at joins against the flips inside one window;
   - deep and shallow samples reported apart, and fills apart.
   **The rule passes if the flips inside blocks are at most 1% on every cohort** — one ceiling,
   whatever the cohort, because "confident" must mean the same everywhere (§6). **The ceiling is
   meant to become the user's parameter** (owner, 2026-11-02), e.g. `--phase-max-link-error 0.01`:
   a user who wants longer blocks accepts more wrong links, and one who feeds a tool that trusts
   `PS` absolutely asks for fewer. That needs a link confidence that *is* an error rate — a number
   per link whose complement, over the links that pass it, matches the share that flip — so that
   the parameter is a threshold on it. The plain all-rounds-agree rule is a yes/no and cannot be
   set to a rate; the draw-based confidence below can (the share of drawn paths that disagree with
   the consensus at the link). §7.1 measures whether it is calibrated: links binned by it, the
   observed flip share against the bin's value, on every cohort. If it is, the option is added with
   1% as its default; if not, the rule stays fixed and the parameter waits. A cohort may phase
   few records; it may not phase them wrongly. Where a cohort misses the ceiling, the rule is made
   stricter by a signal that tightens on its own where the phase is poor, never by excluding the
   cohort. Candidates, each measured before chosen: several paths drawn per round from the same
   forward pass (independent given the round's templates, so a switch they share is the model's,
   not the draw's); a higher lean at the link's records; more rounds for the phase.
2. **The Rust build against the truth** on the same cohorts: the same figures against the prototype
   and Beagle, scored also with `whatshap compare` against the simulator's truth VCF (switch errors
   and block lengths per sample, read from `PS`) — which checks the format through an independent
   reader. Read the output with `bcftools` as well.
3. **Joins in Rust**: a stage test with windows small enough that each contig holds several,
   including a short core followed by a longer margin and a window with no site for a sample; the
   flips at joins against the flips inside windows.
4. **Unit tests**: the two-state alignment on a hand case (a stretch of five marked records swapped
   inside the window is swapped back exactly there); the unaligned link test catching a switch in one
   round that the alignment absorbs; `LEAN_MIN` at the boundary; a multi-allelic `1|2`; allele
   renumbering; an STR record written `1|0` with its repeat counts in that order; a homozygote never
   phased; an *Unchanged* column never phased; a single-record block written unphased; `PS` the first
   phased record's `POS`; two records at one `POS`; `refine-genotypes` refusing an input with `PS`
   and one with `|`.
5. **Invariance**: the same bytes at 1, 4 and 30 threads.
6. **A real cohort**: the tomato test on kimura, where the deep accessions' read-backed phase (from
   WhatsHap, if run) or their consistency across overlapping windows can be compared.

---

## 8. Deferred, with a recommended home

- **A per-record phase quality** (`PQ`): the consensus's lean is one; whether users need it beside
  `PS` is for the owner once blocks exist. Home: a follow-up to this spec.
- **Phase from reads, through calling's chain ids** (owner, 2026-11-02: now or in a following stage).
  A read or read pair that covers two heterozygotes shows which alleles sit together. Calling
  already gives every read a chain id, and the psp stores one per observation
  (`psp_chain_id_encoding.md`). But the pass reads only VCF lines, so the ids never reach it. This
  is the one source of evidence that does not come from the copying model, which is what the outbred
  cohorts need (§6). It would serve two ways: break a link the reads contradict, and keep a link the
  population cannot phase. How much it can reach depends on how far apart a sample's heterozygotes
  are, against fragments of about 350 bases:
  - in the founder simulator, 8–9% of consecutive true heterozygotes of a sample lie within 500
    bases, on every cohort, inbred or outbred (median gap 3.8 kb);
  - in the tomato cohort's calls (62 accessions, pvc 5ffb2444), the five most heterozygous
    accessions have a median of 117 bases between consecutive heterozygous SNPs, and 73% of the gaps
    are within 350 bases. Some of those heterozygotes may be clustered false calls, which was not
    checked.
  So on real diverse material, reads may reach most links, and the simulator understates that.
  Home: a following stage with its own spec. It needs a channel from calling to the pass: per sample,
  the read evidence linking consecutive heterozygous records. First measure the reach on real psps:
  the share of consecutive called heterozygotes with at least one read covering both, by depth.
- **A run-level check of the phase**, like the fill check for genotypes: none exists, because the
  pass cannot see a switch all rounds share. The deep samples' read-backed phase could supply one.
  Home: with phase from reads.

---

## 9. Questions ruled

Both questions were ruled by the owner on 2026-11-02:

- **Q1. How strict is a link?** The stricter rule — every round leans at both records — "for now,
  we'll see how it works" (§2.3).
- **Q2. Phase on a diverse cohort where blocks still flip?** No cohort gate: pvc is a general tool,
  and writing only confident phase must handle every cohort by itself (§6). The ceiling of §7.1
  enforces it.
- **The 1% ceiling** (owner): a parameter the user may want to set. Measured (§7.1), the share of
  disagreeing paths is not an error rate. The owner then ruled to expose it as a **quality proxy**
  instead: `--phase-min-path-agreement`, the share of the 36 paths that must agree, and lean, at a
  link. Lowering it gives longer blocks with more wrong links. The user page states the error rates
  measured for each kind of cohort, and says that 1 is the strictest setting and still leaves
  outbred cohorts at 1.2% to 36% (§6).
- **The default is 0.8** (owner, 2026-11-06; it was 1, all paths, until then): optimized for
  inbreeding 0 to 0.8 and samples at 5× to 30×. From 12 samples up it keeps wrong links under half
  a link in 100 and phases 3% to 27% more heterozygotes right than 1; at 4 samples no share
  keeps under 1 in 100, and 0.8 raises them from 1.5–3.2 in 100 to 3.3–4.8
  (`doc/devel/reports/implementations/phase_agreement_default_2026-11-06.md`).
- **Read evidence through chain ids** (owner, 2026-11-02): now or later. Placed in a following stage
  (§8), with its reach measured first.

The questions as they were put follow, for the record.

Settled by the review and adopted (no ruling needed): the link test on unaligned rounds (§2.3), the
threshold as the prototype's `|w − w′| ≥ 0.5` (§2.2), the join over the last written block with
`k` = 3 (§2.4), `PS` = `POS` of the first phased record (§2.5), single-record blocks unphased, the
repeat counts in `GT` order and the refusal of a phased input (§3), the running consensus (§2.1).

- **Q1. How strict is a link?** The unaligned test alone, or also "every round leans at both
  records". The stricter rule cuts the flips inside blocks by about a third (first cohort
  0.54% → 0.33%; 64-founder 5.4% → 3.5%) and the blocks to less than half their length (N50 39 → 16
  records on the first cohort). §7.1's criterion — flips inside blocks at most a third of the flips
  with no block rule — is met by the plain unaligned test on the first cohort only (31 of 93); the
  stricter rule meets it on both cohorts measured (16 of 93, 48 of 141). **Leaning: the stricter
  rule** — the owner's rule is that a short block costs little and a wrong link a lot. To confirm on
  all cohorts in §7.1 before the build.
- **Q2. Write phase on a diverse cohort where blocks still flip?** The 64-founder inbred cohort runs
  the pass by default and keeps 2% to 3.5% flips inside blocks at best (9 rounds, stricter rule, N50
  8 records). Options: write it (the 200-sample version of that cohort is at 2.3% with the plain
  unaligned test, so more samples help); write it only above a stated sample count; or require a
  stricter `LEAN_MIN`. **Leaning: write it, by §6's rule, and decide on a sample count only if §7.1
  on 64 founders at 50 and 200 samples, outbred included, shows the flips do not fall with samples.**
  (Not adopted: the owner ruled out any cohort gate, §6.)
