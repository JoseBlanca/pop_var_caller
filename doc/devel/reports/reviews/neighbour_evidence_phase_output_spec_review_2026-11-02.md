# Review: writing the phase the neighbour-evidence pass estimates (spec of 2026-11-02)

*Adversarial review of [`neighbour_evidence_phase_output.md`](../../ng/spec/neighbour_evidence_phase_output.md)
against the pass's code (`src/neighbour_evidence/`, `src/vcf/header.rs`,
`src/cli/refine_genotypes.rs`), the two specs it extends, the prototype
(`doc/devel/ng/research/neighbour_evidence/phaser.py`, `phase_accuracy.py`,
`beagle_switch_errors.py`) and the VCF 4.4 specification. The owner's rule — phase is written only
where it is confident, in blocks that may be short and many — is taken as given. Probe scripts and
their outputs are in `tmp/phase_output_review/` of the `neighbour-pass` worktree; the commands are at
the end.*

**Terms used below.** A *flip* is a wrong change of phase between two consecutive phased
heterozygous records of a sample, judged against the simulator's truth (a switch error). A *link* is
such a pair of consecutive phased records. *Flips inside blocks* counts only the links a block keeps.
All numbers are on the founder simulator's cohorts of the spec's §1, from the dumped prototype rounds
the spec's table was made from (`tmp/template_choice/phase_dump/`); like the spec's table, they score
truly heterozygous biallelic SNP records.

## Summary

**One blocker.** The block rule of §2.3 compares rounds *after* they have been aligned to the
consensus, and the alignment removes the very disagreements the rule looks for. Measured on the six
cohorts of §1, it takes out a quarter to a third of the flips (for example 93 flips become 65 on the
first cohort), never the two thirds §7.1 requires. A version of the test that reads each round
before alignment does better (93 become 31), and it needs no alignment at all. Even so, on the
64-founder and outbred cohorts no variant I tried brought the flips inside blocks below 2%.

**What holds.** The spec's §1 numbers reproduce from the dumps (1.5%, 2.4%, 8.0%, 2.2%, 6.7%, 3.5%
with `phase_accuracy.py`). Aligning to the running consensus instead of to round 2 changes nothing
measurable, so §7.1's first comparison can be settled now. The block break at the copying length `S`
is supported: the 49 links longer than 300 kb across the six cohorts flipped 23 times, about half,
which is no better than chance. The memory figure of §4 is right.

**Seven important findings** (I1–I7): the phasing rule of §2.2 is not the prototype's rule
"generalised" but a looser one (I1); joining windows as written would rarely join, and its inputs
are not defined (I2); phase stays unconfident on diverse cohorts whatever the block rule (I3); the
repeat-count fields are written in the wrong order once `GT` is phased (I4); an input that already
carries `PS` or `|` produces a broken or mixed file (I5); multi-allelic records have no defined
orientation in the rounds (I6); and §7 cannot decide the open questions of §9 as written (I7).

---

## Blocker

### B1. The link test compares aligned rounds, and alignment hides the switch errors the test is meant to see

**Section:** §2.3, with §2.1.

**The mechanism.** A switch error in one round is a place where that round's phase turns over
relative to the truth. The alignment of §2.1 lays the round's copies over the consensus with the
fewest changes of state. Where the round and the consensus differ by a switch, the cheapest
alignment puts a change of state exactly at that switch. After alignment, the round agrees with the
consensus on both sides, so its relative orientation across the link agrees with the consensus's,
and the test of §2.3 sees nothing. The test fires only where the alignment did *not* follow the
round: a stretch too short to pay for two changes of state (four marked records or fewer, at
costs 1 and 2), or a record where the round's haplotypes are near-equal and its orientation is a coin
toss.

There is also a bound that makes the test inert at confident records. The consensus is the mean of
the `R − 1` = 4 aligned rounds. Write `d = c₀ − c₁` for the consensus's difference in alternative
probability and `d_r` for one aligned round's. If one aligned round orients a record against the
consensus (`d_r ≤ 0`), then `d ≤ (0 + 1 + 1 + 1) / 4 = 0.75`. So at every link between two records
with `|c₀ − c₁| > 0.75`, every aligned round agrees with the consensus at each record, and the test
cannot fire.

**Measured** (`tmp/phase_output_review/blocks_probe.py`, running consensus, records phased where
`|c₀ − c₁| ≥ 0.5`, the prototype's rule). "Spec" is §2.3 as written. "Unaligned" is the alternative
proposed below: each round's own relative orientation of the two records, read before alignment,
must agree with the consensus's.

| cohort | flips, no block rule | flips inside blocks, §2.3 | flips inside blocks, unaligned | block N50 (records): none / §2.3 / unaligned |
|---|---|---|---|---|
| 10 at 30× + 40 at 3×, inbreeding 0.5, 16 founders | 93 of 6,228 (1.49%) | 65 of 5,837 (1.11%) | 31 of 5,753 (0.54%) | 128 / 51 / 39 |
| same, inbreeding 0.8 | 65 of 2,704 (2.40%) | 40 of 2,539 (1.58%) | 34 of 2,524 (1.35%) | 77 / 43 / 35 |
| same, inbreeding 0.8, 64 founders | 141 of 1,810 (7.79%) | 94 of 1,624 (5.79%) | 87 of 1,605 (5.42%) | 53 / 24 / 23 |
| 50 at 3×, inbreeding 0.5 | 131 of 6,093 (2.15%) | 84 of 5,632 (1.49%) | 58 of 5,563 (1.04%) | 127 / 41 / 34 |
| 50 at 3×, outbred | 635 of 9,796 (6.48%) | 408 of 7,802 (5.23%) | 322 of 7,601 (4.24%) | 107 / 14 / 11 |
| 200, 40 at 30×, inbreeding 0.8, 64 founders | 350 of 10,221 (3.42%) | 263 of 9,518 (2.76%) | 219 of 9,426 (2.32%) | 73 / 42 / 37 |

§7.1 keeps the rule "if the flips inside blocks are at most a third of the flips at all links". The
§2.3 rule leaves 62% to 75% of the flips (65 of 93, 40 of 65, 94 of 141, 84 of 131, 408 of 635,
263 of 350), so it fails on every cohort.

**The bound, checked.** Phasing only records with `|c₀ − c₁| ≥ 0.8` makes the §2.3 rule remove
nothing on the three cohorts I ran: the flips inside blocks equal the flips at all links (53 of
5,369 on the first cohort, 71 of 1,428 on the 64-founder one, 325 of 6,330 on the outbred one;
`blocks_strict.out`). So §2.3's fallback "a stricter `PHASE_MIN`" switches the block rule off. Its
other fallback, more rounds, helps a little: with 9 rounds on the 64-founder cohort the §2.3 rule
leaves 64 flips of 1,480 links (4.3%) against 94 of 1,624 with 5 (`blocks_rounds9.out`).

**Why the unaligned test is the right one.** Whether two records are on the same copy *within one
round* does not depend on which copy that round calls 0, so it can be compared with the consensus
without aligning anything. It catches a switch in any one round. It still cannot catch a switch that
all four rounds make together, and those are a third to two thirds of the flips (31 of 93 on the
first cohort, 87 of 141 on the 64-founder one): the rounds are not independent, since each copies the
haplotypes the round before drew.

**Proposed amendment.**

1. Replace §2.3's test by: *a link is kept only when, in every averaged round, the round's own
   relative orientation of the two records — copy 0 carrying the same allele of each genotype at
   both, or not — is the consensus's.* Store per round, per record, the sign of the round's
   orientation (one bit); a link's relative orientation is the exclusive-or of two bits. The
   alignment of §2.1 then serves the consensus only.
2. Change §2.3's fallbacks: a stricter `PHASE_MIN` is not a lever for the block rule (it is a lever
   on its own, see I1); "every round leans at both records" is one that works (first cohort: 16 of
   4,827, 0.33%, at an N50 of 16 records; 64-founder: 48 of 1,361, 3.5%).
3. State in §2.3 that round agreement cannot see a switch shared by all rounds, and that on diverse
   cohorts this is most of them; I3 takes up what the spec should do there.

---

## Important

### I1. `PHASE_MIN` = 0.9 on `w / (w + w′)` is not the prototype's rule; it phases far more records

**Section:** §2.2, §9 Q1.

**The mechanism.** `w / (w + w′)` is the probability of one orientation *given that the record is
heterozygous*. It is large whenever one orientation is much likelier than the other, even if both
are unlikely. With copy 0 at 0.01 alternative and copy 1 at 0.3, neither copy leans to the
alternative, `|c₀ − c₁|` = 0.29, and the prototype would not phase the record; but
`w / (w + w′)` = (0.3 × 0.99) / (0.3 × 0.99 + 0.01 × 0.7) = 0.977, and §2.2 phases it. The reverse
holds: `|c₀ − c₁| ≥ 0.5` implies `w / (w + w′) ≥ 0.9`, with equality only at 0.75 and 0.25 (the
example of §2.2). So §2.2's rule is strictly looser than the one measured.

The exact generalisation is the *unnormalised* difference. For a biallelic record with alternative
probabilities `c₀` and `c₁`, `w − w′ = (1 − c₀)·c₁ − c₀·(1 − c₁) = c₁ − c₀`. So `|w − w′| ≥ 0.5` *is*
the prototype's rule, and it reads the genotype's two alleles at a multi-allelic record as §2.2
wants.

**Measured** (round-2 alignment, no block rule, `blocks_f16_F0.5_mixed.out`,
`blocks_other_cells.out`):

| cohort | truly het phased, `|c₀ − c₁| ≥ 0.5` | truly het phased, `w/(w+w′) ≥ 0.9` | flips, each | truly homozygous records passing the rule |
|---|---|---|---|---|
| first cohort | 6,312 (84%) | 6,817 (91%) | 1.44% / 1.67% | 100 / 4,425 |
| 50 at 3×, outbred | 9,694 (59%) | 11,679 (71%) | 6.70% / 7.49% | 342 / 8,441 |

The last column matters wherever the written genotype is a wrong heterozygote: the looser rule then
phases it, where the prototype's would almost never. (Homozygous written genotypes are never phased,
so most of those 4,425 records are harmless; how many are written heterozygous was not measured,
because the dumps hold haplotypes, not posteriors.)

**Proposed amendment.** Define confidence as `|w − w′| ≥ 0.5`, call it what it is (the prototype's
rule), and drop the claim that 0.9 generalises it. Hidden records (phasing spec §4.1) make this
sharper: their haplotypes are re-estimated without the sample's reads, so the genotype can be
heterozygous from the reads while both copies lean to the reference — the case where the two rules
differ most.

### I2. Joining windows: comparing over the whole margin will rarely join, and the inputs are not defined

**Section:** §2.4, §2.5, §4 (stage row).

1. **"All agree" over the whole margin fails most of the time.** On the first cohort a sample has
   about 63 phased records a megabase (6,297 phased records over 50 samples and 1.99 Mb of contig).
   The default margin is 1 Mb (`windows.rs`, `WindowSizes::default`), so a join compares about 62
   links on each side. At the measured 1.49% flips a link, the chance that neither window flips
   anywhere in the margin is about 0.985¹²⁴ ≈ 0.15, assuming flips are independent and the margin
   phases as well as a core (it phases worse: a margin's far edge has context on one side only).
   So about one join in six would succeed on the easiest cohort.
2. **The comparison crosses window `w`'s own block breaks.** The records `w` wrote in the next
   margin belong to several of `w`'s blocks (about 4.6 a megabase under §2.3 on the first cohort:
   460 blocks over 50 samples and 1.99 Mb). Across a break, `w`'s relative orientation is exactly
   what `w` judged unconfident, so "all agree" mixes confident and unconfident information.
3. **The margin is not always `w`'s core.** The planner shrinks the core down to the margin's length
   per window (`windows.rs`, `size_ladder`), so the margin of `w + 1` lies inside `w`'s core only if
   `w`'s core was at least `w + 1`'s margin. After a dense stretch a short core can be followed by a
   longer margin, which reaches into `w − 1`'s core.
4. **Nothing says what carries over a window that writes no phase for a sample** — a window with no
   site in its core (`stage.rs`, `release`: `core_sites.is_empty()` skips the solve), or one where
   the sample has no informative record and is not walked (`refine.rs`, `SampleRounds::templates`
   is `None`).
5. **"Window `w + 1`'s first block" is undefined** once `w + 1` computes blocks over its margin:
   its own link test can break between the shared margin records and its first core record.

**Proposed amendment.** Make the join state per sample, not per window: *the last written block*
— its `PS`, and the orientation written at its phased records within the last `margin_bases` bases.
The stage keeps it across windows that write nothing for the sample, and clears it at a contig's
end. A join compares only the records of that block that `w + 1` also phases, taking the `k` nearest
the core boundary (`k` soft, 3 to start), and requires them to lie in one block of `w + 1` that runs
on, with confident links, to `w + 1`'s first phased core record. Measure in §7.3 the share of joins
that succeed against the share of `k`-record windows of a single window's phase that hold no flip.

### I3. On diverse cohorts no block rule tried makes the phase confident

**Section:** §1, §2.3, §6, §9 Q2.

The owner's rule is that unconfident phase is not written. On two kinds of cohort, the blocks stay
wrong at rates far above the 1% the spec names in §2.3:

- **The 64-founder inbred cohort**, which the pass runs on by default (median inbreeding 0.8): the
  best variant in B1, every round leaning at both records, leaves 48 flips of 1,361 links (3.5%);
  with 9 rounds 23 of 1,149 (2.0%) at a block N50 of 8 records.
- **A small outbred cohort**, 10 samples at 3× from 64 founders (`f64_F0_n10_d3`, run here with
  `phaser.py --rounds 5 --haplotypes 16 --one-per-sample`): only 104 of 3,161 truly heterozygous
  records are phased (3%), and inside blocks 34 of 85 links flip (40%) under both the §2.3 rule and
  the unaligned one (`blocks_small.out`). All rounds agree on a phase that is close to random. The
  pass does not run on this cohort by default (its median inbreeding is below 0.5;
  `stage.rs`, `decide`), but `--neighbour-evidence on` runs it.

This is the CLAUDE.md case of a number that comes back confident and wrong at one end of the range.
Round agreement measures how much the draw varies, not whether the model is right.

**Proposed amendment.** Add to §6: *phase is written only when the pass runs by its default
decision* (at least ten samples, median inbreeding at least 0.5); with `--neighbour-evidence on`
below that, genotypes are refined and `GT` stays unphased, and the header says so. Make the
64-founder cohort an explicit gate in §7: if the flips inside blocks stay above a stated limit there,
the owner decides between shorter blocks (every round leaning) and no phase. **I am not sure** the
default-decision gate is enough: the 64-founder cohort passes it and still gives 3.5% at best.

### I4. The repeat-count fields are written in sorted allele order, but their declaration says `GT` order

**Section:** §3, §5 (allele renumbering).

`REPCN` and `REPCNB` are declared "GT order" (`header.rs:616–617`). `rewrite_line` writes them from
`genotypes[sample]`, which is the two alleles ascending (`line.rs:115–124`, `parse_genotype`;
`line.rs:189–194`, `old_pairs[best]`; written at `line.rs:290–304`). Today `GT` is ascending too, so
they agree. Once a record is written `1|0`, the repeat counts must follow copy order, or an STR
record states the wrong copy number for each haplotype.

**Proposed amendment.** In §3 and the §4 row for `line::rewrite_line`: the repeat-count fields are
written in the order of the written `GT`; add an STR record written `1|0` to §7.4's unit tests.

### I5. An input that already carries `PS` or `|` gives a duplicate `FORMAT` key or mixed phase sets

**Section:** §3, *`refine-genotypes`*.

`refine-genotypes` refuses only a VCF whose header says the pass ran
(`refine_genotypes.rs:145–153`). A pvc VCF that went through WhatsHap has `|` and `PS` and is
accepted. Then:

- appending `PS` to a `FORMAT` that has it gives two `PS` keys, which no VCF reader accepts;
- columns the pass leaves unchanged (`Outcome::Unchanged`) keep their input `GT` and `PS`
  (`line.rs:187`, `237–239`), and the input's `PS` values follow the same convention (`POS` of a
  block's first record), so an input block and a pass block can carry the same `PS` and be read as
  one.

**Proposed amendment.** State a rule: on every site record the pass rewrites, input phase is
removed — `GT` rewritten with `/` and any input `PS` replaced — and `PS` is written in place where
the `FORMAT` already has it. Or refuse an input with `PS` in its header, naming it. Add the case to
§7.4.

### I6. A multi-allelic record has no orientation the rounds can store

**Section:** §2.1 step 1, §2.3, §4 (memory), §5 (multi-allelic trap).

- §2.1 marks a record where "both lean, for some allele `a`", then asks whether the round's
  orientation agrees with the consensus's. With three alleles the two can lean on different
  alleles, and "orientation" is not defined until the allele is fixed. Say which: the allele with
  the largest `|c₀(a) − c₁(a)|` in the consensus, the round agreeing if its difference there has the
  same sign.
- §2.3 and the memory figure store one bit a round a record. The bit that matters is the orientation
  of the written genotype's pair `{a, b}`, and the pair is known only after the last round, when the
  posteriors are averaged. One bit per round cannot hold it for `A > 2`.

**Proposed amendment.** For biallelic records keep one bit. For `A > 2`, keep each round's
haplotypes at those records (they are few; the store already holds `2·(A − 1)` numbers there) and
compute the pair's orientation at the end. Say that the prototype phased biallelic SNPs only, as the
phasing spec's review said of its own table (phasing review I9).

### I7. §7 cannot decide §9's questions, and leaves out the cases most likely to fail

**Section:** §7, §9.

| question | what §7 measures | what deciding it needs |
|---|---|---|
| Q1 `PHASE_MIN` | one value, 0.9 | a sweep, and the rule of I1: `|w − w′|` at 0.5, 0.6, 0.8 |
| Q2 every round agrees | the one-third criterion | already answered: it fails on all six cohorts (B1). Compare the unaligned test and "every round leans" |
| Q3 3 records, all agree | in Rust only, after the build | the rule of I2 first; joins can be prototyped by running `phaser.py` on two overlapping halves of a contig |
| Q4 `PS` = `POS` | nothing needed | — |

Missing from §7:

- **The range.** Small cohorts (10 samples; 2 and 3 samples with `--neighbour-evidence on`) and an
  outbred cohort, where I3 found 40% flips inside blocks.
- **Deep against shallow samples**, reported apart. The 30× and 3× samples phase differently, and a
  pooled figure hides which one carries the flips.
- **Fills.** A heterozygous genotype written without reads is phased like any other under the spec.
  Its phase is not worse (first cohort: links touching a fill flip 4 times in 456, 0.9%, against 91
  in 5,777, 1.6%, between records with reads; `fills.out`), but its genotype is held to `GQ` 10 only.
  Measure it apart and say in §2.2 that fills are phased.
- **Multi-allelic and STR records** (I4, I6).
- **A reader round trip**: score the output with `whatshap compare` against the simulator's truth
  VCF (it reports switch errors and block lengths per sample from `PS`), and read it with
  `bcftools`. That checks the format and gives an independent scorer.

---

## Minor

- **m1. `PS` is unique only if no two phased records of a sample share a `POS`** (§2.5: "each starts
  at a different record"). A block broken between two records at one position gives the next block
  the same `PS`. The pass merges overlapping records into one, and no PASS record shares a `POS` in
  the first cohort's refined VCF (0 of 2,728), so this may never occur; window cores never split a
  position (`windows.rs`, `plan_window`). Rule it out in the text: *a link between two records at
  one `POS` that is not confident leaves the second record unphased instead of starting a block.*
- **m2. §7.4's alignment test needs a length.** "A round equal to the consensus with one stretch
  swapped is swapped back exactly there" holds only for a stretch of at least five marked records
  inside the window, or three at an end: two changes of state cost 4, and ties go to "kept".
- **m3. §2.3's reason for the `S` break is wrong; the rule is right.** Beyond `S` a copy keeps its
  template with probability `e^(−1)` ≈ 0.37 (pass §3.5), so the model still links the records
  weakly. Measured, the links longer than 300 kb flip about half the time (23 of 49, all six
  cohorts), which supports the break. `S` is fitted from {30, 300} kb per sample; the prototype ran
  at 300 kb only, so a 30-kb sample's breaks are unmeasured.
- **m4. §2.1's comparison is settled.** Running consensus against round 2, flips with no block rule:
  1.49% against 1.44%, 2.40% against 2.33%, 7.79% against 7.77% on the three mixed-depth cohorts
  (`blocks_f16_F0.5_mixed.out`, `blocks_other_cells.out`). Adopt the running consensus and drop the
  comparison from §7.1.
- **m5. The comparison with Beagle in §1 is not like for like.** The pass's flips are over the 54% to
  84% of heterozygotes it chose to phase; Beagle's are over every heterozygote it called. Both are
  scored on truly heterozygous records, not on written ones. Say so under the table.
- **m6. Single-record blocks, and tools that ignore `PS`.** A block of one record carries no phase
  and exposes a random orientation to readers that ignore `PS` (`bcftools consensus -H` takes `|`
  at face value). VCF 4.4 defines `PS` and ignores it for unphased genotypes, so writing a lone
  heterozygote `a/b` loses nothing. Recommend it, and soften §1's "downstream tools read a block as
  …" to the tools that read `PS` (WhatsHap; SHAPEIT4 with `--use-PS`, which takes `PS` as a 32-bit
  integer and an error rate, documented at 0.0001 — a rate these blocks do not reach).
- **m7. §2.5's attribution.** The `POS`-of-the-first-record convention is the VCF specification's
  own recommendation ("The recommended convention is to use the position of the first variant in
  the set as the PS identifier", VCF 4.4, Genotype fields). Cite that; I did not verify GATK's.
- **m8. Code locations in §4 and §3.** The `##neighbourEvidence=` value is built in `stage.rs:66–79`
  (`as_header_value`), not in `vcf/header.rs`; `phase_min` goes there. `refine-genotypes` builds its
  own header from the input's lines (`refine_genotypes.rs:164–200`) and replaces only `GP`'s
  declaration, so `GT`'s new description and `PS`'s declaration must be added there too. Add both to
  the §4 table.
- **m9. Samples the pass did not walk.** A sample without an informative record in the window keeps
  its `GP` shares as both haplotypes (`refine.rs:346–375`), so `w = w′` and it is never phased; its
  outcome is `Unchanged`. Say so in §2.2, so a coder does not phase calling's heterozygotes there.
- **m10. Orientation bits at margin records.** §4 stores them "at its core records", but a joined
  block's links run through `w + 1`'s left margin (I2), so the bits are needed there too. The memory
  is unchanged in order: the window's 39,000 records at 3,000 samples, at a byte each, is 117 MB.

## Figures checked

- **§4's memory.** A consensus row of `2·(A − 1)` `f32` is 8 bytes a biallelic record, as the store
  (`refine.rs:351`). At the floor window of phasing spec §5, 39,000 records × 3,000 samples × 8 bytes
  = 0.94 GB, plus 13,000 core records × 3,000 × 1 byte = 39 MB of bits: the spec's "about 0.9 GB" is
  right. The time claim (one pass a round, a few operations a record) holds.
- **§1's table** reproduces with `phase_accuracy.py` on the dumps: 92 of 6,248 (1.5%), 64 of 2,707
  (2.4%), 145 of 1,807 (8.0%), 132 of 6,092 (2.2%), 644 of 9,594 (6.7%), 360 of 10,200 (3.5%); the
  last-round figures of §2.1 likewise (229 of 6,499, 3.5%, on the first cohort).

## Citations checked

Correct: `src/vcf/header.rs:610` (`GT`, "Genotype, unphased"); `refine::refine_window_in_rounds`
(`refine.rs:475`); `HaplotypeStore` (`refine.rs:337`), whose row is `2·(A − 1)` `f32` a record
(`refine.rs:351`); `RefinedRecord` (`refine.rs:158`); `stage::NeighbourEvidenceStage`
(`stage.rs:140`); `line::rewrite_line` (`line.rs:134`), which renumbers the alleles it keeps
(`line.rs:208–216`); pass §3.5 (the copying length `S`); phasing spec §4.1 (hidden records) and §5
(the floor-sized window). Wrong place: the `##neighbourEvidence=` line (m8). The two research scripts
named in the preamble exist and are the ones that made §1's table.

## Reproducing

From the worktree root, host Python:

```
D=tmp/template_choice/phase_dump
uv run --with numpy python tmp/phase_output_review/blocks_probe.py $D f16_F0.5_n50_d3+30x10
CONS=running PHASES=lean LINKS=none,spec,raw,raw_lean \
  uv run --with numpy python tmp/phase_output_review/blocks_probe.py $D <cell>
PHASE_LEAN=0.8 CONS=running PHASES=lean LINKS=none,spec,raw \
  uv run --with numpy python tmp/phase_output_review/blocks_probe.py $D <cell>
uv run --with numpy python tmp/phase_output_review/run_phaser4.py --rounds 9 --haplotypes 16 \
  --one-per-sample --dump tmp/phase_output_review/dump9 tmp/founder_grid/f64_F0.8_n50_d3+30x10
uv run --with numpy python tmp/phase_output_review/fills_probe.py $D tmp/founder_grid/<cell>
```

`run_phaser4.py` runs `phaser.py` unchanged with four worker processes; its first five rounds match
the dumps exactly (141 flips of 1,810 links on the 64-founder cohort both times).
