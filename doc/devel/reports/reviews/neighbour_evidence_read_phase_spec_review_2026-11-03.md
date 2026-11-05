# Review: `neighbour_evidence_read_phase.md` (draft of 2026-11-03)

*Adversarial review of the spec that carries read-backed phase evidence from the cohort merge to
the phase writers. Branch `neighbour-pass`, worktree `pop_var_caller-neighbour-pass`, spec
uncommitted. Every file:line below was opened and checked on 2026-11-03. The owner's six settled
points are taken as given; findings about them are only where the spec fails to honour them or where
they create a problem the spec does not address.*

## Verdict

**Not ready to build. It needs one design change and several corrections first.**

The plumbing is sound. The capture point, the four allele renumberings, the ordered thread, the
spill route and the command removal all check out against the code. A few citations are off by a
few lines and none is wrong in substance.

What fails is the statistics where reads meet the neighbour-evidence pass. A link made by reads
inside the pass has no floor on how sure the genotype is. A het (heterozygous genotype) the pass
writes with as little as a 50% chance of being right can be linked with full confidence. That alone
can take a cohort past the 1-in-100 ceiling the spec says it holds.

Four further problems need fixing in the text before anyone writes code:

- the error formula reverses its sign for a read allele whose error is above one half;
- the read cap picks reads by a chain id that depends on what ground was walked earlier;
- the hidden-duplication filter's spill grows by several times at the top of the range, and the
  spec does not size it;
- the testing plan cannot detect the errors `e_floor` exists for.

**Counts: 1 Blocker, 9 Major, 16 Minor.**

---

## Blocker

### B1. Read-made links inside the pass have no genotype-confidence floor (§2.7, §6 trap 7)

**Evidence.** §2.6 holds reads-only phase to `GQ` ≥ 20, because "a link to a homozygote miscalled
heterozygous is a genotype error the reads cannot see". §2.7 has no such floor. It relies on
`LEAN_MIN`, and §6 says "`LEAN_MIN` (pass) [is] what keep[s] it out". That does not hold.
`LEAN_MIN` = 0.5 bounds `w − w′`, with `w = c₀(a)·c₁(b)` (`src/neighbour_evidence/phase.rs:33-36`,
`:138-140`). For `c₀(a) = 1` and `c₁(b) = 0.5` the lean is exactly 0.5, so the consensus gives the
het only about a 50% chance.

Two §2.7 rules then link with no further check:

- **"not kept / `|L|` ≥ `T_link` → linked by reads"**: two consecutive candidates the population
  could not link;
- **"a refined heterozygote the consensus does not lean on can be phased by reads"**: this record
  passed no lean test at all.

**Failure.** Take a record the pass wrote heterozygous where the sample is really homozygous. This
can happen even where the sample has reads, since the pass overturns a read-based homozygote (owner
decision 2), and hidden records are re-estimated without the sample's reads (phase output spec §5).
At such a record every fragment shows the same allele. The sign of each fragment then depends only
on which copy it came from at the other record, so two fragments agree half the time. Two agreeing
fragments give 7.8 nats, above `T_link` = 4.6 (§2.5's own arithmetic).

The wrong-link rate on such links is therefore about P(het wrong) × ½. At the lean floor that is up
to a quarter, against a ceiling of 1 in 100. The 64-founder outbred cohort is where these links
matter most, and it is also where the pass's hets are least certain.

**Fix.** Apply §2.6's floor to every link that rests on reads: the record's written `GQ` (the pass's
refined `GQ`) must be at least 20 at both ends. Delete the §6 claim that `LEAN_MIN` protects. Add a
unit test of a pass-called het whose reads all show one allele, and check that it is never
read-linked.

---

## Major

### M1. The ceiling budget does not add up even with the `GQ` floor (§2.6, §1.1)

`T_link` = ln 99 spends the whole 1-in-100 budget on read error. Separately, `GQ` ≥ 20 allows each
end's het to be wrong 1 time in 100, and a wrong het yields a confident link of random sign half the
time once two fragments agree (B1). Worst case per link: 1% + 2 × 1% × ½ ≈ 2%.

There is a second gap. The `GQ` at 3× in a cohort comes largely from the cohort's allele-frequency
prior, and the spec never checks that it is calibrated.

**Fix.** Split the budget. For example, require `GQ` ≥ 30 at both ends and keep `T_link`, or
require `T_link` = ln 199 with `GQ` ≥ 20. Add to §7.1 a calibration check of `GQ` on the hets that
were read-linked, at 3× and 30×.

### M2. The weight `w` reverses sign for `e` > ½ (§2.5)

`w = ln([(1−eᵢ)(1−eⱼ) + eᵢeⱼ] / [eᵢ(1−eⱼ) + eⱼ(1−eᵢ)])` is negative when one `e` exceeds ½. At
`eᵢ = 1` and `eⱼ = 0.01` it is ln(0.01/0.99) = −4.6 nats, a full-strength vote for the opposite
relation.

§2.2 caps the byte at 60 but sets no lower bound. Phred 0 to 3 occurs in real input:

- the walk sets the losing mate of an overlapping pair to base quality 0, which is ε = 1, and that
  read "still counts" (`src/locus_generation/pileup/genome_walk.rs:1712-1717`;
  `open_record.rs:3613-3618` measures 9 read-positions in 1,000 on HG002 carrying ε = 1);
- Illumina Q2 tails give ε = 0.63.

A sequence whose only read is such a read has mean error Phred near 0.

**Fix.** Clamp `e` to at most ½, or drop read alleles below Phred 10, and test the boundary.

### M3. The read cap chooses reads by a chain id that depends on walk history (§2.1, §2.4, §7.5)

§2.4 keeps "the 32 reads whose hashed chain id is smallest". Chain ids are not a function of the
read. The allocator counts up across the whole walk, so the same read gets a different id when
earlier ground was or was not walked. The walker's own oracle records this:

> "The `chr2` read is id 0 when its segment is walked alone and id 4 when two `chr1` segments and a
> satellite came first … this cannot be, by the design of the allocator"

(`src/run/walker.rs:1432-1446`).

**Failure.** Wherever a sample has more than 32 reads at a record, the reads kept, and so `L` and
the written phase, change in three cases:

- between a whole-genome run and a `--regions` run over the same ground;
- between call-from-psps, where ids are the file's from gather time, and call-from-alignments over a
  subset;
- between psps gathered over different ground.

This breaks the segment-independence property that `run_streaming.md` §12 asks of the output. §2.1's
statement that direct and psp mode "walk the same regions in the same order … so they should" give
the same ids holds only when the analysed ground is identical. Psps thinned by `read_cap.rs` already
differ from direct mode (`src/run/psp_source/read_cap.rs:33-37`), which the spec does not mention.

**Fix.** Make the cap's key invariant to walk history. One option is the read's rank of first
appearance within its walked region, which `with_chain_ids_renumbered` shows is invariant. State the
direct-versus-psp guarantee as "same analysed ground, unthinned psps".

### M4. §7.5's direct-versus-psp test cannot see an id mismatch (§7.5)

Ids matter only through the cap, and the cap binds only above 32 reads per sample per record. A
fixture at ordinary depth passes whether or not the ids agree, so "which also checks that both give
the same chain ids" is untrue.

**Fix.** Run that test on a fixture with more than 32 reads per sample at the linked records, or
compare the captured `(read rank, allele)` lists directly.

### M5. Paralog mis-mapping is correlated across reads, and `e_floor` does not cover it (§2.5, §6 trap 8)

§6 says `e_floor` stops "a mis-mapped pair in a duplicated region link[ing] two records with full
confidence". It stops one pair. In a collapsed duplication every read from the other copy carries
that copy's alleles at both records, so `n` such fragments give `n` × 3.9 nats in the same
direction. Unmarked PCR duplicates are a second case: two fragments from one molecule, counted as
independent.

What actually guards these sites is the hidden-duplication filter, because both writers require
`PASS` (`src/neighbour_evidence/line.rs:54`; §2.6). That filter is optional (`--hidden-dup-fdr 0`
turns it off). It also has measured misses at 5× (memory: removes 94 true and 26 false at chr20 5×).

**Fix.** Say plainly that `e_floor` handles independent errors only. Name the filter as the guard
against correlated errors and record what happens when it is off. Consider whether reads-only phase
should be refused when the filter is off.

### M6. The 1-in-100 goal is stated for "every cohort", but the spec's own arithmetic says it fails (§1.1, §2.7)

§1.1 says "Hold every phase writer to the phase output spec's ceiling … on every cohort". §2.7 then
computes that the 64-founder outbred cohort with refinement on stays near 4.3% (4.7% × (1 − 0.085))
"even if reads broke every wrong link they reached". The goal is unachievable as written, and §7.1
would fail by design.

**Fix.** State what the build is held to:

- the ceiling for read-made and read-confirmed links on every cohort;
- for population-only links, no worse than `phase_output_2026-11-02.md`;
- the outbred case left to Q3 and §8.

### M7. The non-leaning attach rule tolerates contradicting reads up to `T_link` (§2.7)

"Its read link must agree with the block's relation or be weaker than `T_link`." So a read link
against the block's relation at −4.5 nats (odds of 90 to 1 against) is accepted. A population link
breaks at −2.2 (`T_break`). The two rules disagree about the same evidence.

The rule also never checks two non-leaning records between one pair of candidates against each
other.

**Fix.** Refuse the attachment when the other side's read link is at or below −`T_break` against
the block's relation. Attach several non-leaning records between one candidate pair as a chain, with
each consecutive link tested.

### M8. Capture runs on every cohort locus, most of which are never written (§2.2, §4.2)

The one-record branch of `alleles_of_sample`, the ordinary SNP locus, does not touch chain ids today
("no read is consulted", `src/run/cohort_merge/build.rs:2049-2052`, `:2092-2108`). The capture adds a
copy and a sort of every id for every sample at every cohort locus. That is 90,000 entries (about
1.4 MB) a locus at 3,000 samples and 30×, and 900,000 at 300×, before any cap. Yet
`a_written_genotype_carries_an_alternative` (`src/run/records.rs:76`) discards most loci a few
steps later, and the spec says so itself ("most loci are never written"). §4.2 calls this "expected
to be small" with no number.

There is a further cost: the uncapped lists ride on `VcfRecord` to the ordered thread, because the
32 cap is applied only in the keeper.

**Fix.** Apply the per-record cap at capture; it needs no state from other records. Either measure
the capture at 200 samples before building it, or defer it until after the written-alternative test
by keeping the per-sample observation slices borrowed until the `finish` closure.

### M9. The spill is sized per record but never in total, and the filter is on by default (§4.1)

Pass one parks every called record of the run and pass two reads them all
(`src/run/paralog_filter/mod.rs:9-15`). The filter is on unless `--hidden-dup-fdr 0`
(`src/cli/call_from_psps.rs:74`). Today each spilled record holds its line plus a fixed row per
sample (`SpilledSample`: a window and two `u32`, `spill.rs:150-161`).

At 30× the spec's own figure is about 18 read alleles × 6 bytes ≈ 110 bytes a sample. That is
several times the row it sits beside, and its size grows with depth up to the cap. Pass two pays to
read fragments it never uses.

**Fix.** Estimate the whole-genome spill at 3,000 samples and 30×, with and without Q1. Put the
fragments in a second spill read only by pass three, or in a separately skippable section of the
entry.

---

## Minor

1. **§2.6 arithmetic.** "At 30×, about 29 pairs in 100 are reached, by about 4 fragments each" is
   wrong. 3.93 is the research note's mean over *all* adjacent pairs (§2 table, `E[spanning]`), so a
   reached pair has about 3.93 / 0.29 ≈ 13.6 fragments. "Most reached pairs have exactly one" at 3×
   should read "about half" (research note §2, point 2). "Close to empty" overstates as well: about
   half of the 20% of pairs reached have two or more fragments, enough for `T_link` once `GQ`
   allows.
2. **§4.1 `k`.** `k` = 0.6 is used at 3,000 samples, where records are dense and `k` is near 1, so
   the 3,000-sample rows are low by about 1.6×. The 50-sample rows use a 39,000-record window. At 50
   samples the record budget does not bind (120 M sample-records, `windows.rs:37`), so the window is
   the full 4 Mb. At tomato's density of about 25 records a kb (199,641 records over 8 Mb) that is
   about 100,000 records, so the 50-sample rows are low by about 2.5×.
3. **Q1's saving at 3×.** The saving is overstated. A sample showing 3 reference reads at a site
   with alternative frequency 0.05 keeps P(het) ≈ 0.095 × ⅛ ≈ 0.012, far above 10⁻³. The floor
   bites at 30× and on rare variants. "A few percent in an inbred cohort" needs this depth caveat.
4. **§4.1 round driver.** "16 regions in flight" cannot happen at 3,000 samples, where
   `round_shape_for` returns one region (`src/cli/calling_run.rs:442-452`). The streaming driver
   runs there, and the ordered thread is the only bottleneck.
5. **§2.2 mapping-quality cap.** The cap is already in `q_sum`: each read's error is
   `max(base error, MAPQ error)` (`open_record.rs:3627-3629`). Say so. Also note that
   `q_sum / num_obs` is the *geometric* mean error, measured at 25 to 44 times below the arithmetic
   mean (`open_record.rs:3600-3611`). That gap is the strongest argument for `e_floor`.
6. **§1.1 versus §3.** §1.1 says "without changing any VCF column calling writes today", but §3
   rewrites `GT` and appends `PS` on every line whenever refinement does not run, which includes
   every one-sample run. Every byte-identity digest (`src/cli/cross_platform_digests.rs`) will need
   re-recording; say so.
7. **Cap ordering is unspecified.** The spec does not say whether the cap of 32 is applied before or
   after the "reaches another record" test. The two orders keep different reads. Use the same hash
   as `read_cap.rs` (`splitmix64`), so that thinning and the cap agree.
8. **Duplicate entries.** §2.2 removes an id carrying two *different* alleles. It should also
   deduplicate an id listed twice with the same allele, so one read can never count twice.
9. **Fragments that span non-consecutive hets are ignored.** A pair whose mates cover hets 1 and 3,
   with het 2 in the unsequenced gap, links 1 and 3, but only consecutive pairs are tested (§2.6,
   §2.7). This loses evidence and can leave a block that contradicts reads it never consulted.
   State it, or test each new record against every earlier block member within `R`.
10. **§2.6 phaseable set and holding.** "Phaseable" is defined by a fragment to "another such
    record", which is circular. Combined with the skip rule, deciding a record's link to the right
    can need records up to 2`R` ahead, not `R`.
11. **§2.7 unspecified cases.** Two cases are left open. Does the copying-length break at `S`
    (`phase.rs:375`) apply to read-made links? What happens when a window has no paths
    (`phase.rs:347-349` returns all `None`, so reads are unused)? It also leaves out the `PS` of a
    block whose first record is an attached non-leaning one: §2.5 of the phase output spec names a
    block after its first record's `POS`, and `join_and_name_core_blocks` keys on that name.
12. **Chain-id uniqueness has a hole.** The allocator pairs mates by `qname` alone
    (`chain_id_allocator.rs:290-291`), so two read groups of one sample with colliding read names
    within the mate window share an id. §2.1's "an id is never given to two reads" should name this
    case. Key pending mates by (read group, `qname`) or say why that is unnecessary.
13. **Citations.**
    - `AlleleRemap::candidate_for` is at `src/calling/allele_candidates/mod.rs:390`; `:352` is the
      struct.
    - The per-block restatement of live ids is in the ⚠ paragraph of `psp_chain_id_encoding.md` §4
      and in `psp_file_format.md` §3.2 (`src/psp/chain_ids.rs:23-29`), not in §4's opening.
    - The phase output spec's 36-path rule is in its §6, not §2.3.
14. **§5 scope.** `founder_grid.py:154` passes `--neighbour-evidence`, which the rename breaks. §5
    says it "already calls with the pass off and on" as if it needed no change.
    `neighbour_evidence_phasing.md` also names the command and is missing from the docs list.
15. **§7.2's `e_floor` fit.** The fit mixes genotype error into read error: at a false het about
    half the fragments disagree with the majority. Restrict it to hets in the GIAB truth. The fit
    also cannot see correlated errors (M5). The GIAB v4.2.1 truth VCF's phasing is "not checked";
    check it before relying on it.
16. **Clarity.**
    - Terms used without definition: *ordered thread*, *frontier claim*, *sink*, *elongate*,
      *pending mates*, *consensus*, *copying length `S`*.
    - "Record" means a VCF record in §1.3 and a per-position observation record in §2.2's branch
      names. Pick one word for the second.
    - Several claims carry no number: "share most of their reads", "expected to be small",
      "negligible", "most of those are linked".
    - The simulator is calibrated by construction. Bases are misread exactly as their quality
      states, and there are no duplicates, chimeras, paralogs or indel errors
      (`examples/ng_str_stutter_benchmark.rs:92-106`). So §7.3 can measure reach and how reads
      combine with population links, but not whether `e_floor` holds. Its STR density is the
      simulator's own, not a genome's. Measure the STR-segment loss on tomato or GIAB psps as well.
    - §7.1 needs a link count per cohort. A single sample at 3× gives tens of read-made links,
      which cannot show a rate below 1 in 100 (Q3 already says about 300 are needed).

---

## Checked and correct

None of these needs a change:

- capture site and `moved_to` (`build.rs:645-791`, `:773-781`);
- the two branches of `alleles_of_sample` and their treatment of overlapping mates (`:2055-2062`,
  `:2082-2213`);
- the allocator's lifetime across regions and chromosomes (`chain_id_allocator.rs:253-262`;
  `genome_walk.rs:265`, `:282`, `:915`, `:985`; `generator.rs:1386-1424`);
- region-crossing reads named twice (`locus_generation/mod.rs:393-400`);
- the monotone hash cutoff in thinning (`read_cap.rs:21-29`, `:104-109`);
- the `finish` closures (`callers.rs:1105`, `:7469`) and `evidence_for_output` (`records.rs:622`,
  `:688-697`);
- the renumbering in `drop_alternatives_no_sample_calls` (`assemble.rs:466-520`);
- the ordered-thread call sites (`psp_caller.rs:790-801`, `callers.rs:837-844`) and the round fold
  (`rounds.rs:86-224`, `claim` at `:227`);
- writer and pass-three routing (`writer.rs:148-158`, `:188`, `:217-231`;
  `pass_three.rs:256-257`);
- the pass's call sites (`stage.rs:245`, `phase.rs:36`, `:338`, `:374-384`, `refine.rs:1024`,
  `line.rs:138`);
- the §5 code citations;
- the §2.5 numbers: 3.90, 7.8 and 2.77 nats;
- the phase output report's 55 of 1,170.

The pass reads `(a, b)` in the line's input numbering (`refine.rs:988-992`), as §2.3 step 4
assumes.

---

## Questions for the owner

1. **The genotype floor for read-made links in the pass (B1, M1).** Should those links require the
   pass's refined `GQ` ≥ 20 (or 30), as reads-only phase does? Recommend yes. Without it the
   ceiling cannot hold.
2. **Reads that confidently contradict a population link.** When `|L|` ≥ `T_link` against a kept
   population link, §2.7 only breaks the block. Should the reads' relation be written instead, as
   it is for a link the population did not keep? Recommend break only, until §7.3 shows how often
   such reads are themselves wrong.
3. **Q1 against decision 2.** At 3× the floor saves little (Minor 3), so its case rests on 30× and
   above. Does the owner accept a posterior floor that departs from "kept whatever the genotype"?
4. **The spill (M9).** Should fragments go to a second spill read only by pass three, or should the
   cost at 3,000 samples be measured first?
5. **Reads-only phase with the hidden-duplication filter off (M5).** Should it still run? Recommend
   that it runs, and that the run log warns that collapsed duplications can then give confident
   wrong links.

---

## Dispositions (2026-11-03)

*Decided by the coordinator for the owner; applied to
[`neighbour_evidence_read_phase.md`](../../ng/spec/neighbour_evidence_read_phase.md) the same day.
Every citation the amendment touched was re-opened.*

| finding | disposition | where in the spec |
|---|---|---|
| B1 | Accepted. Every link that rests on reads (reads-only, read-made in the pass, records attached by reads) needs the written `GQ` ≥ 30 at both ends; the pass's refined `GQ` when it runs. The claim that `LEAN_MIN` protects is deleted. Unit test added: a pass-called het whose reads all show one allele is never read-linked. | §2.5, §2.6, §2.7, §6, §7.4 |
| M1 | Accepted. `T_link` raised to ln 199 = 5.29 nats (still two good fragments at `e_floor` 0.01). Budget shown: 0.5% read error + 2 × 0.1% × ½ ≈ 0.6%. `GQ` calibration on read-linked hets at 3× and 30× added. | §2.5, §7.1 |
| M2 | Accepted. `e` clamped at ½; a Phred-0 overlap loser then weighs exactly zero; boundary tests added. | §2.5, §7.4 |
| M3 | Accepted. Cap key = `splitmix64(id − base)`, `base` the first id issued in the walked region. Direct mode: `chain_id_counters().chain_allocations` at region open. psp mode: not stored today; a per-walked-region `base` (8 bytes each) must be added to the psp, with a format version rise; older psps fall back to the raw id with a log line. Guarantee stated as "same output whenever the walked regions are the same; thinned psps excepted". | §2.4, §4 |
| M4 | Accepted. The direct-versus-psp test uses > 32 reads per sample at linked records, adds a whole-contig against `--regions` run, and compares the captured `(rank, allele)` lists directly. | §7.5 |
| M5 | Accepted. `e_floor` stated as independent errors only; the hidden-duplication filter, through `PASS`, named as the guard against correlated ones; with the filter off read phase still runs and the log warns. The simulator cannot test this; GIAB segmental duplications with the filter on and off, and the tomato read check, can. | §2.5, §6, §7.3, §7.7 |
| M6 | Accepted. Goal restated: the ceiling for read-made and read-confirmed links on every cohort; population-only links no worse than `phase_output_2026-11-02.md`; the outbred case left to Q3 and §8. | §1.1, §7.1 |
| M7 | Accepted. Attachment refused when the other side's read link is at or below −`T_break` against the block; several non-leaning records between one candidate pair attached as a chain, each consecutive link tested. | §2.7, §7.4 |
| M8 | Accepted. The cap is applied at capture, before the reach test. Capture is measured at 200 samples before the build commits to every cohort locus; above 5% of merge time it moves to the `finish` closure, re-running the id walk on the minted records still held there (`rounds.rs:163`; the streaming driver to be checked). | §2.2, §2.4, §4.2, §7.6 |
| M9 | Accepted. Fragments go in a separately skippable section of the spill entry, read only by pass three. Whole-genome spill at 3,000 samples and 30×: about 17 TB today, about 48 TB more without Q1's floor, about 2 TB more with it (arithmetic). | §4, §4.1 |
| Question 2 | Accepted: a confident read contradiction of a kept population link only breaks it, because it is also what a missed collapsed duplication looks like; writing the reads' relation is deferred until §7.3 measures it. | §2.7, §8 |
| Minor 1 | Applied: about half of reached pairs have one fragment at 3×; at most about 1 pair in 10 linkable; 13.6 fragments a reached pair at 30×. | §2.6 |
| Minor 2 | Applied: `k` ≈ 1 at both scales; the 50-sample rows use the full 4-Mb window (about 100,000 records); table recomputed. | §4.1 |
| Minor 3 | Applied: Q1's saving is at 30× and above and at rare variants; the 3× worked example is in §4.1 and Q1. | §4.1, §9 Q1 |
| Minor 4 | Applied: above about a thousand samples the round driver runs one region (`round_shape_for`, `calling_run.rs:456`); the ordered thread is the bottleneck. | §4.1 |
| Minor 5 | Applied: the MAPQ cap is already in `q_sum` (`open_record.rs:3627-3629`); the mean is geometric, 25 to 44 times below the arithmetic mean. | §2.2 |
| Minor 6 | Applied: §1.1 no longer claims calling's columns are unchanged; the digests (`cross_platform_digests.rs`) must be re-recorded. | §1.1, §3 |
| Minor 7 | Applied: cap first (at capture), then the reach test, and why the order matters; `splitmix64` as `read_cap.rs` uses. | §2.2, §2.4 |
| Minor 8 | Applied: an id listed twice with one allele keeps one entry. | §2.2 |
| Minor 9 | Applied: fragments over non-consecutive records stated as lost for linking and used as a check: a new block member is tested against earlier members within `R`, and a contradiction at −`T_break` ends the block. | §2.6, §7.4 |
| Minor 10 | Applied: eligible and phaseable defined without circularity; the stage holds `2R`. | §2.6, §4.1 |
| Minor 11 | Applied: `S` applies to read-made links; a window with no paths phases nothing and uses no reads; `PS` is the first record's `POS` whatever its source, with attachment done before `join_and_name_core_blocks` (`phase.rs:450`). | §2.7 |
| Minor 12 | Applied: the qname collision is named (one sample is one psp, which may hold several read groups); recommended for this build: key pending mates by (read group, read name), at the cost of unpaired mates whose `RG` tags differ. | §2.1, §4 |
| Minor 13 | Applied: `candidate_for` at `mod.rs:390`; restatement cited to `psp_chain_id_encoding.md` §4 ("Random access is preserved by restating"), `psp_file_format.md` §3.2 and `src/psp/chain_ids.rs:24-27`. **Not applied**: the 36-path rule is in the phase output spec's §2.3 (its lines 213-222), so the citation to §2.3 stays; §6 there reports the outbred failure under that rule. | §1.3, §2.1, §2.3 |
| Minor 14 | Applied: `founder_grid.py:154` must pass `--refine-genotypes`. **Not applied** for `neighbour_evidence_phasing.md`: it names neither the command nor the option (checked), so it needs no change and is listed as such. | §5 |
| Minor 15 | Applied: the `e_floor` fit is restricted to GIAB truth heterozygotes and cannot see correlated errors; the GIAB truth phasing is to be checked before use. | §7.2 |
| Minor 16 | Applied: ordered thread, frontier claim, sink, pending mates, consensus and copying length defined in §1.3; *elongate* removed; per-position records renamed *minted records*; unquantified claims given numbers or marked unmeasured; the simulator's calibration limit stated (`examples/ng_str_stutter_benchmark.rs:103-106`) and the STR-segment loss measured on tomato and GIAB too; §7.1 states link counts and the about-300 rule. | §1.3, §2.2, §4.1, §7 |
| Question 1 | Adopted as B1/M1 (`GQ` ≥ 30). | §2.5 |
| Question 3 | Kept as the spec's Q1, with the 3× caveat. | §9 |
| Question 4 | Settled by M9: a skippable section of the entry. | §4.1 |
| Question 5 | Adopted: read phase runs with the filter off, and the log warns. | §2.5 |

The spec's former Q4 (thresholds) is folded into §2.5 as starting values; §9 keeps Q1, Q2, Q3 and Q5.
