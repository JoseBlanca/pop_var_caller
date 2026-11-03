# ng — phase from reads: carrying read evidence from the merge to the phase writers

*Spec, 2026-11-03. No code yet; this settles the design. It realises the item deferred in
[`neighbour_evidence_phase_output.md`](neighbour_evidence_phase_output.md) §8, "Phase from reads,
through calling's chain ids", and uses the method of
[`../research/phasing_imputation/05_read_backed_phasing.md`](../research/phasing_imputation/05_read_backed_phasing.md)
§6 (that note's `src/ng/...` paths are stale: the modules are now at the crate root).*

*The owner settled these points on 2026-11-03, and this spec treats them as fixed:*

1. *Read evidence is worth building now, without first measuring how far reads reach. In the
   tomato cohort's calls, the five most heterozygous accessions have 73% of the gaps between
   consecutive heterozygous SNPs under 350 bases (phase output spec §8).*
2. *Evidence is kept whatever the sample's genotype at the record, because refinement can turn a
   homozygote into a heterozygote or back. The unit is the **fragment** (§1.3), and only reads that
   reach two or more variant records are kept.*
3. *There is **one `PS`**. When refinement runs, it writes one phase set built from population and
   read evidence together. When it does not run, calling writes `|` and `PS` from reads alone.*
4. *Fragments are never written to any VCF. They travel in memory beside each record, and through
   the spill when the hidden-duplication filter is on, either to refinement or to the writer of
   phase from reads alone.*
5. *Refinement runs only inside a calling run. The separate `pvc refine-genotypes` command is
   removed, and `--neighbour-evidence {auto,on,off}` is renamed `--refine-genotypes {auto,on,off}`
   on `call-from-psps` and `call-from-alignments` (§5).*
6. *A downstream phaser reading pvc's fragments is not a requirement.*

*Revised the same day after an adversarial review,
[`neighbour_evidence_read_phase_spec_review_2026-11-03.md`](../../reports/reviews/neighbour_evidence_read_phase_spec_review_2026-11-03.md).
Its dispositions are listed at the end of that file. The main changes: every link that rests on
reads needs `GQ` ≥ 30 at both ends and odds of 199 to 1 (§2.5); an error probability is never above
one half; the read cap's key no longer depends on what was walked earlier (§2.4); the goal says what
the outbred cohorts can and cannot expect (§1.1).*

---

## 1. What this is

A read, or a read pair, comes from one molecule. If it covers two heterozygous sites of a sample,
the two alleles it shows sit on the same chromosome copy, unless one base is a sequencing error or
the read is mapped to the wrong place. This is **read-backed phase**. It is the one source of phase
evidence that does not come from the neighbour-evidence pass's copying model.

That matters because of how the population phase fails. The pass writes phase from the cohort's
haplotypes (phase output spec §2). On outbred cohorts it is confidently wrong: at 50 outbred
samples from 64 founders, 55 of 1,170 links inside phase sets are wrong (4.7%,
`phase_output_2026-11-02.md`), against a ceiling of 1 in 100. Its draws all come from one model, so
no number of draws can see the error (phase output spec §6). Reads are an independent check.

calling already knows which read showed which allele. Every observation the generic path mints
names its reads by **chain id**, and the cohort merge uses the ids to decide which reads covered a
whole locus. The ids are then dropped before calling sees them. **This spec keeps them**: it records,
at every written record, which allele each read showed. It then keeps the reads that reach two or
more records and hands them to whichever step writes phase:

- **Refinement runs** (the neighbour-evidence pass): read evidence enters its link test. Reads
  break a population link they contradict, and join records the population cannot link (§2.7).
- **Refinement does not run**: calling writes `GT` with `|` and `PS` from reads alone (§2.6). This
  covers one-sample runs, outbred cohorts under the default `auto`, and `--refine-genotypes off`.

### 1.1 Goals

- Keep, for each sample, every read that reaches two or more written variant records, with the
  allele it showed at each and an error probability.
- Carry that evidence beside each record from the merge to the phase writer, in memory and through
  the spill. No fragment field appears in any VCF. Where refinement does not run, `GT` and `PS`
  change to carry phase from reads (§3).
- Phase from reads alone when refinement does not run. Combine reads with the population link test
  when it does. Write one `PS` either way.
- **Hold the build to these error rates, measured inside phase sets:**
  - links made or confirmed by reads: at most 1 wrong in 100, on every cohort and at every depth;
  - links the population makes alone: no worse than `phase_output_2026-11-02.md` measured on the
    same cohort.
  
  On outbred cohorts at 3×, most population links are out of a fragment's reach. Reads can break
  the wrong ones they reach, but cannot bring the rest under 1 in 100 (§2.7). That case stays open
  (Q3, §8).
- Produce output that is deterministic and identical at any thread count, in both drivers, and in
  direct and psp mode whenever the walked regions are the same (§2.1).
- Remove the `refine-genotypes` command and rename the option (§5).

### 1.2 Non-goals

- **Writing fragments to a file.** Never in the VCF, never in a side file (owner, decision 4).
- **Phase on the STR path.** STR observations carry no chain ids (`SequenceObservation::chain_ids`,
  `src/locus_generation/mod.rs:398-400`), so STR records get no read evidence (§2.1).
- **Linking a read across two walked regions.** A read that crosses a walked-region boundary has
  two ids, so it counts as two reads here (§2.1). This is deferred (§8).
- **Per-read quality.** The observation keeps only sums. The mean over an observation's reads stands
  in for each read's quality, as the research note recommends (§2.2).
- **Ploidy other than two.** No phase is written, as for the pass today.
- **A phase quality field (`PQ`).** It stays deferred, as in the phase output spec.
- **Assembly-style phasing** (WhatsHap's dynamic programming, HapCUT2's max-cut). With short reads
  a block almost never has three sites that disagree, so the pairwise log-odds is the whole method
  (research note §6.2).

### 1.3 Vocabulary

Terms about reads and records:

- **Chain id**: the number the walk gives one read, or one read pair whose two mates share a single
  id (`src/locus_generation/pileup/chain_id_allocator.rs:1-10`). The walk pairs mates by holding
  each first mate's name until its second mate arrives; those held names are the **pending mates**.
- **Walked region**: the stretch of reference one walk covers. On the generic path that is one
  generic region between two STR segments.
- **Minted record**: one sample's observations over one stretch of a locus, as the walk minted them
  (`SampleLocusObservations`). A cohort locus can hold several minted records of one sample.
- **Variant record**: a record calling writes to the VCF, with at least one alternative allele.
- **Read allele**: the allele one read showed at one variant record, numbered as the record numbers
  its alleles (0 = REF), with that read's error probability there.
- **Fragment**: for one sample, one read (one chain id) that has read alleles at two or more variant
  records (owner's definition). A read with a read allele at a single record is not kept.
- **Reach**, `R`: the longest distance in bases between a fragment's first read allele and a later
  one that is still kept. Default 1,000 (§2.4).

Terms about the run:

- **Ordered thread**: the one thread that receives called records in genome order, in both drivers,
  and passes them on to the VCF or the spill. Workers call loci in parallel; this thread puts the
  results back in order.
- **Frontier claim**: the round driver's rule, on the ordered thread, that a locus whose ground an
  earlier locus already owns is dropped (`claim`, `src/run/cohort_merge/rounds.rs:227`).
- **Sink**: where the ordered thread puts each finished record: the VCF writer, or the spill when
  the hidden-duplication filter is on (`CalledRecordSink`, `src/run/paralog_filter/pass_one.rs:92`).

Terms about phase, shared with the phase output spec:

- **Link**: the relative orientation of two heterozygous records of a sample. A link is **cis** when
  the lower-numbered alleles of the two genotypes lie on one copy, and **trans** otherwise.
- **Link log-odds**, `L`: the natural log of how much more likely the fragments covering both
  records are under cis than under trans (§2.5). Positive values favour cis.
- **Consensus**: the pass's average of each sample's two haplotypes over its rounds (phase output
  spec §2.1). Where it leans by at least `LEAN_MIN` = 0.5 at a heterozygote, the pass can orient it.
- **Population link**: a link the pass keeps because all 36 of its drawn paths agree and lean
  (phase output spec §2.3).
- **Copying length**, `S`: the pass's fitted length of the stretches a sample copies from one
  template, 30 or 300 kb. No link is kept between records farther apart than `S` (phase output spec
  §2.3).

---

## 2. The method

### 2.1 Which reads name the same molecule: the id guarantee

**Read evidence requires that two observations at different loci carry the same chain id only if
they come from the same read or read pair.** The written guarantee is narrower:
`psp_chain_id_encoding.md` §1.1 promises equality "over a window of one locus". **This spec
promotes it.** Within one sample's evidence (one psp file, or one sample's alignments in direct
mode):

- **an id is never given to two reads**, with one exception named below. The allocator lives for
  the whole run. `next_id` is kept across chromosome changes (`chain_id_allocator.rs:253-262`,
  `reset`) and across walked regions: each walk borrows the allocator and hands it back
  (`src/locus_generation/pileup/genome_walk.rs:265`, `adopting_chain_ids`; `:282`,
  `into_chain_ids`; `src/locus_generation/pileup/generator.rs:1386-1424`);
- **within one walked region, a read keeps one id at every locus it covers.** The psp stores the real
  ids and restates the whole set of live ids at the start of every block
  (`psp_chain_id_encoding.md` §4, "Random access is preserved by restating";
  `psp_file_format.md` §3.2; `src/psp/chain_ids.rs:24-27`), so equality survives block boundaries;
- **a read that crosses a walked-region boundary gets a new id in the second region.** Entering a
  region resets the pending mates (`genome_walk.rs:915`, `begin_region`; the reset at `:985`), and
  the field's documentation says such a read is "met twice and named twice"
  (`src/locus_generation/mod.rs:393-397`). The read becomes two shorter fragments. This loses
  evidence but never creates a false link.

**The exception: two reads with one name in one sample.** Pending mates are keyed by read name alone
(`chain_id_allocator.rs:290-291`). A sample is one psp, but one psp can hold several read groups
(157 of 1,707 surveyed tomato samples carry more than one library, `read_groups.md` §1). If two of
its libraries reuse a read name within the mate window, two reads from different libraries share an
id. The merge's own documentation names the case: "a merged file in which two libraries reuse a read
name" (`src/run/cohort_merge/build.rs:2270`). Such a pair would be a false fragment. **Recommended
for this build: key pending mates by (read group, read name).** The cost is that a pair whose two
mates carry different `RG` tags, which is legal input (same documentation), is no longer paired. It
becomes two fragments, a loss rather than a false link. The change alters the ids a walk writes, so
psps gathered before it keep the exception.

`psp_chain_id_encoding.md` §1.1 should be amended to state this guarantee, and to forbid any future
encoding that renumbers ids per block. Apart from the exception, the code already meets it.

**Two things lose reach and are accepted:**

- **STR segments.** On the generic path a walk stops at every STR segment (glossary, *Generic
  path*). A read pair spanning an STR segment is therefore two fragments, and a heterozygote on one
  side cannot be linked by reads to one on the other. STR records carry no read alleles, because STR
  observations name no reads. When refinement runs, an STR record is phased by the population alone,
  as today. When it does not, an STR record is never phased. Generic records that the edge-deletion
  gathering merges into a widened STR record lose their read alleles too (§2.3). *Unmeasured*: the
  share of a sample's consecutive heterozygote pairs within `R` whose two records lie on opposite
  sides of an STR segment (§7.3).
- **Thinning.** A psp written under a looser read cap is thinned per record by a hash of the chain
  id (`src/run/psp_source/read_cap.rs:104-109`). A read kept at one share is kept at every larger
  one, so thinning removes read alleles but never invents one. The walk's own cap,
  `--max-reads-per-position`, discards reads per position without naming them; such a read simply
  has no read allele there.

### 2.2 Capture: in the merge, before the ids are dropped

The ids and the composed alleles coexist in one place only: the emit callback of
`AlleleTable::assemble` (`src/run/cohort_merge/build.rs:645-791`; callback `:673-721`). There each
allele a sample's reads showed is added to a tally keyed by allele and read group, and the ids go no
further. **Capture happens in that callback.** For each covering sample it fills a list of
`(chain id, preliminary allele index, error Phred)`:

- **The sample has one minted record in the locus** (`AlleleBacking::OneSequence`): one entry per
  id in the sequence's `chain_ids`. This branch does not read ids today ("no read is consulted",
  `build.rs:2049-2052`), so capture adds that work here.
- **Several minted records** (`AlleleBacking::OneRead { records, sightings }`): one entry for the
  read, whose id every sighting shares. Reads the merge removed are absent, because they are never
  emitted: missing from one minted record, partial at one, or showing two things at one
  (`alleles_of_sample`, `build.rs:2082-2213`).

**The error Phred.** It is the observation's mean read error, `q_sum / num_obs`, converted from
natural-log error to Phred. Two facts about this number:

- each read's error in `q_sum` is already the larger of its base error and its mapping error
  (`minted_ln_read_error`, `src/locus_generation/pileup/open_record.rs:3627-3629`). So HapCUT2's
  cap of base quality at mapping quality (research note §1.4) is already applied;
- a mean of logarithms is a **geometric** mean. It sits far below the arithmetic mean error:
  44 times below on HG002 at 300×, 25 times on the tomato cohort (`open_record.rs:3609-3611`). The
  number is optimistic for any one read, which is the strongest reason for the error floor of §2.5.

Where a read spans several minted records, the weakest record's mean is used, as
`share_of_one_read` already does for support (`build.rs:1780`). The value is rounded to an integer
and stored in one byte, capped at 60. Rounding at capture means every later step sees the same
integer, so §2.5's arithmetic runs on identical inputs in every driver and mode.

**Duplicates and disagreements.** The list is sorted by id. An id listed twice with the same allele
keeps one entry, so one read never counts twice. An id with two different alleles is removed: a pair
whose two mates overlap and disagree puts one id in two sequences of one minted record, and
`alleles_of_sample` emits both on the one-record branch (its documentation describes the case,
`build.rs:2055-2062`). This follows HapCUT2, where disagreeing overlapping mates contribute nothing
at that site (research note §1.4).

**The cap of 32 is applied here**, per sample and locus (§2.4 gives the key). It needs no state from
other loci, so applying it at capture keeps a deep sample's list at 32 entries for the rest of the
route.

**The preliminary allele index is renamed with the rest of the table**: `sort_canonically` returns
`moved_to` (`build.rs:773`), and the capture's indices are rewritten through it exactly as the
support rows are (`:774-781`). The list is added to `SampleSupport` (`build.rs:1267`) as a new field.

**Cost, and the condition on building it this way.** Capture copies, sorts and caps every id of every
sample at every cohort locus, and most cohort loci never become records. At 3,000 samples and 30×
that is 90,000 entries a locus before the cap. **The build measures the capture at 200 samples
before committing to it** (§7.6). If it costs more than 5% of the merge's time, capture moves to the
`finish` closure (§2.3) and runs only for loci that become records. That is possible because the
sample's minted records are still held while the locus is called: in the round driver the call runs
inside `cache.with_observations` (`rounds.rs:163`). The capture would keep borrowed slices of the
locus's minted records until the closure, and re-run the id walk of `alleles_of_sample` there. Whether
the streaming driver holds them as long must be checked at the build.

### 2.3 Through calling: from merge alleles to the record's alleles

A read allele has to end in the numbering of the line the phase writer reads. That numbering
changes four times, and every change already exists for `AD`:

1. **The merge table to the candidates.** In the `finish` closure, which receives the locus's
   inference and its `AlleleRemap` (streaming driver `src/run/callers.rs:1105`, round driver
   `:7469`), each read allele's merge index goes through `AlleleRemap::candidate_for`
   (`src/calling/allele_candidates/mod.rs:390`). `None` means candidate selection dropped that
   allele, so the read allele is dropped. The candidate id is the record's allele index, as
   `evidence_for_output` uses it for `AD` (`src/run/records.rs:622`; the loop at `:686-694`). This
   happens only for loci that become records, after `a_written_genotype_carries_an_alternative`
   (`records.rs:76`).
2. **Alternatives no sample calls.** `assemble_record` drops them and renumbers the rest
   (`src/vcf/assemble.rs:137`; `drop_alternatives_no_sample_calls`, `:466`). The read alleles are
   renumbered with the same table, and those on a dropped allele are removed.
3. **The edge-deletion gathering** (`src/run/edge_deletion_gathering.rs`). Records that pass through
   keep their read alleles unchanged. A generic record merged into a widened STR record loses its
   read alleles, because the merged record's alleles are haplotypes over a wider span, and the STR
   record has none. The count lost this way goes to the run report beside the gathering's own counts.
4. **The pass's allele subsetting.** `rewrite_line` drops alternatives the pass did not call
   (`src/neighbour_evidence/line.rs:138`). The pass reads read alleles in the line's *input*
   numbering, which is the numbering it reads genotypes in (`refine.rs:988-992`), decides links, and
   carries orientation through the renumbering as it already does for population phase (phase output
   spec §5). A read allele on an allele the pass drops is uninformative, because the written genotype
   cannot contain that allele.

The read alleles ride on the record: `VcfRecord` gains a field holding, per run sample, the sorted
list of read alleles. Everything that moves a record (the gathering, the frontier claim, the sink)
moves them with it. A locus the frontier claim drops takes its read alleles with it.

### 2.4 Keeping fragments only: the fragment keeper

Most read alleles belong to reads that reach one record only. Owner decision 2 keeps only fragments,
so a step must find which reads reach a second record. That needs records in genome order and
per-sample memory across records. **It runs on the ordered thread, right after the edge-deletion
gathering releases each record and before the record is prepared for the VCF or the spill** (psp
mode `src/run/psp_caller.rs:790-801`, direct mode `src/run/callers.rs:837-844`).

**The ordered thread is the same in every driver.** In the round driver, workers build and call
whole building regions, and the results are folded back in genome order on the calling thread
(`merge_and_call_in_rounds`, `rounds.rs:86-226`; `take` is "the sink, on the calling thread, in
genome order"). Workers produce read alleles per locus and keep no state across loci. The per-sample
"recent reads" state lives only on the ordered thread, so the building regions' boundaries and the
worker count cannot change the result.

**The rule.** The keeper holds each record until the newest record's position is more than `R`
bases past it, or until the contig ends. When it releases a record, it keeps, for each sample, a read
allele only if the same read has a read allele at another record no more than `R` bases from the
read's first one. Records leave in the order they arrived; the keeper never reorders. A read allele
beyond `R` from its read's first one is dropped and counted. The cap of §2.2 was applied before this
test, at capture; the order matters, because capping first and testing reach second keeps different
reads than the reverse.

- **`R` = 1,000 bases**, a starting value (Q2). It covers ordinary paired-end inserts: the founder
  simulator draws fragments of 350 ± 50 bases (`examples/ng_str_stutter_benchmark.rs:146-147`),
  and HapCUT2 joins mates up to 1,000 bases by default (research note §1.4). A longer `R` only holds
  more records.
- **The cap: at most 32 read alleles per sample per record**, a starting value (Q2). WhatsHap caps at
  15 fragments a site for a similar reason (research note §1.3). It bounds memory per record whatever
  the depth (§4.1).
- **Chain ids become per-sample fragment numbers**, a `u32` counted from 0 within each contig in order
  of first read allele, which halves what is stored. They are never written anywhere. A contig change
  clears them.

**The cap's key must not depend on what was walked before.** Chain ids count up across a whole walk,
so the same read has a different id when other ground was walked first. The walker's own test records
it: a read is "id 0 when its segment is walked alone and id 4 when two `chr1` segments and a satellite
came first" (`src/run/walker.rs:1437-1440`). A cap keyed on the raw id would change which reads a deep
sample keeps, and so the written phase, between a whole-genome run and a `--regions` run over the same
ground. **The key is the read's rank within its walked region, hashed**: `splitmix64(id − base)`,
where `base` is the first id the allocator issued in that walked region, and `splitmix64` is the hash
`read_cap.rs` already uses (`read_cap.rs:145`). A walked region is walked identically alone or inside
a longer walk (the oracle `walker.rs:1432-1446` describes), so its ids differ between the two only by
`base`, and `id − base` is the same. The 32 reads with the smallest key are kept. Neighbouring records
of similar depth hold mostly the same reads, so they mostly keep the same ones; how much they overlap
is unmeasured.

**Where `base` comes from:**

- **Direct mode**: the generator reads the allocator's counters when a region opens
  (`chain_id_counters`, `genome_walk.rs:294`). `chain_allocations`, the number of ids minted so far
  (`chain_id_allocator.rs:185`), is the next id to be issued, so it is `base`. The generator hands it
  with each minted record of the region.
- **psp mode: not available today and must be added.** No psp records it. The gather writes, for each
  walked region, its `base`: 8 bytes a walked region, a few megabytes for a whole genome walked as
  several hundred thousand generic regions. Where it sits is the psp format spec's decision
  (`psp_file_format.md`), and the format version rises. The reader hands `base` with each minted
  record it decodes. A psp written before the change has none; for it the key falls back to the raw
  id, and the run log says that the phase of samples deeper than 32 reads then depends on the walked
  ground.
- **In memory**: `SampleLocusObservations` gains `chain_id_base: u64`, 8 bytes a minted record.

**The guarantee, stated: the same output whenever the walked regions are the same, in either mode;
thinned psps excepted**, since thinning keys on the raw id (`read_cap.rs:108`) and already makes
psp-mode calls differ from direct-mode calls (`read_cap.rs:33-37`).

The keeper's work is per sample and independent between samples. When the ordered thread is too slow
(§4.2), a batch of released records can be processed in parallel over chunks of samples. The result is
the same, because no sample reads another's state.

### 2.5 A link's log-odds

**Whoever phases counts fragments for whichever two records it treats as consecutive heterozygotes of
the sample** (owner decision 2): calling with its own genotypes, the pass with its refined ones, which
may lie several records apart.

Take sample `s`, records `i` and `j`, and genotypes `{aᵢ, bᵢ}` and `{aⱼ, bⱼ}` with `a < b` at each.
A fragment with read alleles `xᵢ` and `xⱼ` says nothing about this link if either allele is outside
its genotype; it showed a third allele or an allele the genotype lacks. Otherwise it says *cis* when
`(xᵢ = aᵢ)` equals `(xⱼ = aⱼ)`, and *trans* when not.

Each read allele has an error probability

```
e = min( max(10^(−q/10), e_floor), ½ )
```

where `q` is its error Phred (§2.2).

- **The floor, `e_floor`**, stands for errors that base and mapping quality do not see, and that
  strike reads independently: a chimeric fragment, a PCR error, a base the geometric mean of §2.2
  rates too well.
- **The ceiling of ½**: a read allele with `e` above ½ would vote for the opposite relation (the
  formula below reverses sign). Phred 0 to 3 is real input: the walk gives the losing mate of an
  overlapping pair base quality 0, an error of exactly 1, and that read still counts
  (`genome_walk.rs:1712-1717`; 9 read-positions in 1,000 on HG002, `open_record.rs:3613-3618`).
  Clamped at ½, such a read allele contributes a weight of exactly zero.

The fragment's weight is HapCUT2's edge weight (research note §1.4): the log of the chance it shows
what it shows if its claim is true, over the chance if false.

```
w = ln( [(1 − eᵢ)(1 − eⱼ) + eᵢ·eⱼ] / [eᵢ(1 − eⱼ) + eⱼ(1 − eᵢ)] )
L = Σ over fragments of (+w if the fragment says cis, −w if trans)
```

The sum runs in fragment-number order, in `f64`, so it is the same at any thread count.

**Correlated errors are outside this model.** In a collapsed duplication, every read from the other
copy carries that copy's alleles at both records; `n` such fragments add `n` × 3.9 nats in one
direction. Unmarked PCR duplicates are two fragments of one molecule counted as two. `e_floor` does
not cover either. **The guard is the hidden-duplication filter**: both phase writers use `PASS`
records only (`line.rs:54`; §2.6), and the filter is what removes a collapsed duplication from `PASS`.
The filter is on by default and off with `--hidden-dup-fdr 0` (`src/cli/call_from_psps.rs:74`). It
also misses some duplications at low depth. With the filter off, read phase still runs, and the run
log warns that collapsed duplications can then give confident wrong links. Unmarked duplicates are the
alignment preparation's job; nothing here detects them.

**The thresholds**, starting values to be checked by §7:

| threshold | value | meaning |
|---|---|---|
| `T_link` | ln 199 = 5.29 nats | reads make or keep a link: odds of at least 199 to 1, 1 wrong link in 200 if the model holds |
| `T_break` | ln 9 = 2.20 nats | reads contradict a population link enough to break it: odds of 9 to 1 against it |
| `e_floor` | 0.01 | the independent error per read allele that no quality sees; fitted from HG002 first (§7.2) |
| `GQ_link` | 30 | the written `GQ` both ends of a read-resting link must have |

**What one fragment is worth.** With `e_floor` = 0.01, one fragment contributes at most
`ln(0.9802 / 0.0198)` = 3.90 nats, however good its bases. Two agreeing fragments of good quality give
7.8 nats. So a link needs at least two good fragments to reach `T_link`, and one fragment alone can
break a population link: a contradicting read allele of Phred 13 (e = 0.05) beside one at 0.01 gives
2.77 nats.

**Why a genotype floor, and the budget.** A read-resting link has two ways to be wrong:

- **the reads mislead**: at `L` ≥ `T_link`, at most 1 time in 200 (0.5%) if the error model holds;
- **one end is not really heterozygous.** At a true homozygote every fragment shows the same allele
  there, so each fragment's sign depends only on which copy it came from at the other record. Two
  fragments agree half the time, and two agreeing fragments pass `T_link`. A wrong heterozygote at
  either end therefore gives a confident link of random sign: wrong with probability about ½.
  `GQ` ≥ 30 bounds each end's chance of being a wrong heterozygote at 1 in 1,000.

Together: 0.5% + 2 ends × 0.1% × ½ ≈ **0.6%**, inside the 1% ceiling. The `GQ` is the written one: the
pass's refined `GQ` when refinement runs (`rewrite_line` replaces `GQ`, `line.rs:7`), calling's when it
does not. This rests on `GQ` being calibrated. At 3× a cohort's `GQ` comes largely from the allele
frequency prior, and that is checked, not assumed (§7.1). `LEAN_MIN` gives no such protection: a lean
of 0.5 is reached with one copy's allele probability at 1 and the other's at 0.5, about an even chance
that the record is heterozygous (`phase.rs:33-36`, `:138-140`).

`T_break` is lower than `T_link` because breaking a link costs only block length, and the owner's rule
is that a short block costs little and a wrong link a lot.

### 2.6 Phase from reads alone: when refinement does not run

The writer of read-only phase is a stage in `VcfWriter`, beside the neighbour-evidence stage and used
when that stage is absent (`src/vcf/writer.rs:148-158`), at ploidy two. It sits in the writer, not in
calling, because `FILTER` is final only there: with the hidden-duplication filter on, the verdict is
written in pass three, which feeds the spilled lines back through `write_line`
(`src/run/paralog_filter/pass_three.rs:256-257`). Like the pass, it reads `GT`, `GQ` and `FILTER` from
the line and rewrites `GT`.

**Which records of a sample it may phase.** A record is **eligible** when it is `PASS`, its `GT` is
heterozygous and its `GQ` is at least 30 (§2.5). An eligible record is **phaseable** when it has a read
allele of a fragment that also has a read allele at another eligible record. Eligible records with no
such fragment are skipped, not treated as gaps, so two heterozygotes on either side of one without
reads can still be linked.

**Blocks.** Walk the sample's phaseable records left to right. Each consecutive pair is linked when
`|L|` ≥ `T_link`; `L`'s sign gives cis or trans. A pair without such a link ends the block. A block of
one record is written unphased, as in the phase output spec §2.2. `PS` is the `POS` of the block's
first record (phase output spec §2.5). The first record is written `a|b` (copy 0 carries its lower
allele), and each later record's orientation follows from the links.

**Fragments that skip a record.** A pair whose mates cover records 1 and 3, with record 2 in the
unsequenced gap between them, links 1 and 3 but is not counted for the consecutive pairs 1–2 and 2–3.
The evidence is lost for linking. It is still used as a check: when record `h` joins a block, every
earlier record of the block within `R` of `h` is tested against it, and if any such pair has `L` at or
below −`T_break` against the relation the block implies, the block ends before `h`.

**Holding.** A link needs fragments at both records, and fragments reach at most `R` bases. Deciding
whether a record is phaseable needs its fragments' later read alleles, up to `R` ahead, and deciding
its link to the next phaseable record can need another `R`. The stage therefore holds a line until the
stream is more than `2R` past it, and writes lines in arrival order.

**What it gives across the range** (arithmetic from research note §2, not measured). At one
heterozygote a kilobase and 350-base inserts:

- at 3×, about 1 adjacent pair in 5 is reached by any fragment, and about half of those by exactly
  one. One fragment never reaches `T_link`, so at most about 1 pair in 10 can be linked, before the
  `GQ` floor removes more. Sparse, as the owner expected;
- at 30×, about 29 pairs in 100 are reached, by about 13.6 fragments each on average (the note's 3.93
  fragments a pair, over all pairs, divided by 0.29). Nearly all reached pairs can be linked.

Dense heterozygosity (the tomato accessions' 73% of gaps under 350 bases) raises the reachable share.
A single deep sample gets read-backed phase with no cohort at all.

### 2.7 Reads in the pass's link test: when refinement runs

The pass phases in `SamplePhase::phase_heterozygotes` (`src/neighbour_evidence/phase.rs:338`).
Records whose consensus leans by at least `LEAN_MIN` (`phase.rs:36`) are **candidates**; consecutive
candidates are linked when enough of the 36 paths agree and the two are no farther apart than `S`
(the `linked` closure, `phase.rs:374-384`). For a pair of consecutive candidates, write `L_pop` for `L`
signed so that positive supports the relation the consensus gives them. Read evidence changes the
test as follows:

| population link (36 paths) | reads | result |
|---|---|---|
| kept | `L_pop > −T_break` (silent, weak or agreeing) | linked, the population's relation |
| kept | `L_pop ≤ −T_break` | **broken**: the block ends |
| not kept, both ends `GQ` ≥ 30, within `S` | `|L|` ≥ `T_link` | **linked by reads**, with the reads' relation |
| otherwise | | not linked, as today |

**A confident contradiction only breaks, for now.** Where reads contradict a kept population link with
`|L|` ≥ `T_link`, the reads' relation is not written in place of the population's. Two sources
disagree, and a confident read contradiction is also what a collapsed duplication the filter missed
looks like (§2.5). Until §7.3 measures how often such contradictions are themselves wrong, the block
ends there.

**A read-made link joins two blocks, and may flip the second one as a whole.** Each candidate's
orientation comes from its own consensus. If the reads say the relation is the opposite of what the
two consensus orientations give, every record of the following block is swapped. That is exact for
the same reason a window join swaps a whole window (phase output spec §2.4): a block's phase has
meaning only relative to itself. The walk runs left to right, so each swap decision is made once.

**Refined heterozygotes the consensus does not lean on can be attached by reads.** Such a record is
not a candidate today. After the blocks are built, the non-leaning heterozygotes between two
consecutive candidates `p` and `q`, with `GQ` ≥ 30, are taken in order as a chain starting from `p`:

- each is attached when `|L|` ≥ `T_link` links it to the previous record of the chain (`p`, or the
  last record attached), with the reads' relation;
- the chain stops at the first record that fails, and records after it stay unphased;
- if `p` and `q` are in one block, the chain's last record is tested against `q`: when `L` against
  `q` is at or below −`T_break` against the relation the block implies, every record of this chain is
  left unphased.

Attaching never changes an existing link or a population orientation.

**Unchanged and clarified:**

- the copying-length break at `S` applies to read-made links too, though read links are always
  shorter than `R` and `S` is at least 30 kb;
- the break at two records with one `POS`, and single-record blocks written unphased;
- a window whose rounds drew no paths phases nothing today (`phase.rs:347-349`); reads are not used
  there either;
- `PS` is the `POS` of a block's first record, whether that record is a candidate or an attached one.
  Attachment runs before blocks are named, so `join_and_name_core_blocks` (`phase.rs:450`) sees the
  final first record;
- joins across windows (phase output spec §2.4). Read links are shorter than `R` ≤ 1 kb and windows
  overlap by their margin, so every read link falls inside one window. The joins compare written
  orientations, which now include read-made ones.

**What this can do for the outbred cohorts, and what it cannot.** Reads break the wrong population
links they reach. On the founder simulator, 8 to 9% of a sample's consecutive true heterozygotes lie
within 500 bases (phase output spec §8). Even if reads broke every wrong link they reached there, the
64-founder outbred cohort would stay near 4.3% (4.7% × 0.915), above the ceiling. The simulator
places SNPs sparsely (`SNP_RATE` = 10⁻³, `examples/ng_str_stutter_benchmark.rs:189`), and the tomato
data show far more short gaps. So the test of this spec on outbred material needs a denser simulator
(§7.3). Where reads cannot reach, they cannot help; §2.8 and Q3 address that.

### 2.8 The read check of population links

Where reads reach a population link with `|L|` ≥ `T_link` and both ends have `GQ` ≥ 30, they give an
estimate of how often population links are wrong that does not come from the model's draws. **The
stage counts, per sample and in total, the population links reads reached that way, and how many of
those the reads contradicted.** The run log and the run report state the counts. This is the run-level
check the phase output spec §8 deferred "to phase from reads".

It is reported, not acted on, in this stage. The counts are biased: reads reach only links under `R`,
and short links are probably the population model's most accurate, so the rate is a lower bound for
longer links. Whether to use it as a rule that tightens population-only links by itself is Q3.

---

## 3. What the VCF says

- **No fragment field, anywhere** (owner decision 4).
- **Refinement runs:** as the phase output spec §3, with phase now from both sources under one `PS`.
  The `##neighbourEvidence=` header line gains `read_link_min=5.29,read_break_min=2.20,
  read_error_floor=0.01,read_gq_min=30,read_reach=1000`. It keeps its name, which names the pass
  rather than the option.
- **Refinement does not run, ploidy two:** `GT` written `a|b` and a `PS` field appended after the
  existing `FORMAT` keys on read-phased heterozygotes, `.` elsewhere. `PS` is written on every record,
  so the `FORMAT` string is the same on every line. The header declares `GT` and `PS` with the
  declarations the pass already uses (`src/vcf/header.rs:607`, `:610`), and gains
  `##readBackedPhase=<link_min=5.29,error_floor=0.01,gq_min=30,reach=1000>`. The declarations are
  written even when no record is phased, so the header does not depend on the data.
- **This changes calling's output in every run where refinement does not run**, including every
  one-sample run: the `GT` declaration, the `PS` column, and `|` where reads phase. Every recorded
  byte-identity digest (`src/cli/cross_platform_digests.rs`) has to be re-recorded.
- **`REPCN`/`REPCNB`** are not affected: STR records are never read-phased.

---

## 4. Where it sits in the code

| what | existing code | change |
|---|---|---|
| capture | `AlleleTable::assemble`, emit callback (`build.rs:673-721`); `SampleSupport` (`:1267`) | per-sample `(chain id, allele, error Phred)` list; renamed through `moved_to`; duplicates merged, disagreeing ids removed; capped at 32 by key (§2.2, §2.4) |
| the cap's key | `SampleLocusObservations`; the generator (`genome_walk.rs:294`); the psp gather and reader | `chain_id_base` per minted record; a per-walked-region `base` stored in the psp (§2.4) |
| to the record | `finish` closures (`callers.rs:1105`, `:7469`); `evidence_for_output` (`records.rs:622`); `assemble_record` (`assemble.rs:137`, `:466`) | merge index → candidate → record allele; dropped alleles removed; a read-allele field on `VcfRecord` (§2.3) |
| gathering | `EdgeDeletionGathering` | carries read alleles; drops the absorbed generic records' and counts them |
| fragment keeper | new, on the ordered thread after the gathering (`psp_caller.rs:790-801`, `callers.rs:837-844`) | holds `R` of records; keeps fragments; numbers them (§2.4) |
| sink and spill | `CalledRecordSink` (`pass_one.rs:92-247`); `SpillEntry` (`spill.rs:133-161`); pass three (`pass_three.rs:256-257`) | the spill entry gains the fragments in a separately skippable section, read only by pass three (§4.1); pass three hands them to the writer with the line |
| writer | `VcfWriter::write_record` (`writer.rs:188`), `write_line` (`:217-231`) | `write_line` takes the record's fragments beside the line and passes both to whichever stage runs |
| link log-odds | new module, shared | §2.5's weight, clamp and sum; used by both writers |
| reads-only phase | new stage in `VcfWriter` | §2.6 |
| the pass | `NeighbourEvidenceStage::push` (`stage.rs:245`), `HeldLine`; `phase_heterozygotes` (`phase.rs:338`), called from `refine.rs:1024` | `HeldLine` holds the fragments; the window's solve passes them and the refined `GQ` to the link test; §2.7's table and attachment; §2.8's counts |
| the walker | `ChainIdAllocator` pending mates (`chain_id_allocator.rs:290-291`) | keyed by (read group, read name) (§2.1) |
| headers | `header.rs`, `stage.rs` (`as_header_value`) | §3 |
| the command and the option | §5 | §5 |

### 4.1 Memory

One stored read allele is about 6 bytes (a `u32` fragment number, the allele, the error Phred), 8
with overhead. Write `k` for the share of a record's read alleles that reach another record within
`R`. It depends on how dense the records are: roughly `1 − exp(−300 / spacing)` for 300 sequenced
bases a fragment, which is 0.6 at a record every 300 bases and near 1 at a record every 40. The
figures below are arithmetic, not measurements.

**The fragment keeper and the reads-only stage** hold `R` and `2R` of records: at a record every 40
bases, 25 to 50 records. At 3,000 samples and 30× that is 25 to 50 records × 3,000 samples × 30 read
alleles × 8 bytes = 18 to 36 MB.

**The pass** holds a whole window's lines, so it holds their fragments too. Read alleles per sample per
record are min(depth, 32) × `k`. Two window sizes apply. At 3,000 samples the window is the phase
output spec's floor of 39,000 records (§4 there). At 50 samples the record budget
(`budget_sample_records` = 120 million, `src/neighbour_evidence/windows.rs:37`) does not bind, so the
window is its full 4 Mb (2 Mb core and 1 Mb margins); at the tomato cohort's density, about 25 records
a kilobase, that is about 100,000 records. Records at both scales are dense, so `k` ≈ 1:

| cohort | per record | per window | beside the pass's own store at that size |
|---|---|---|---|
| 50 samples, 3× | 1.2 kB | 0.12 GB | small |
| 50 samples, 30× | 12 kB | 1.2 GB | comparable |
| 3,000 samples, 3× | 72 kB | 2.8 GB | the store's two buffers are 1.9 GB |
| 3,000 samples, 30× | 720 kB | 28 GB | too large |
| 3,000 samples, 300× (capped at 32) | 770 kB | 30 GB | too large |

**At the top of the range this does not fit, because reads are kept at every sample and record,
including the many where the sample is certainly homozygous.** There a read allele says nothing about
phase, since both copies carry one allele. Q1 proposes a floor: drop read alleles where calling's
probability that the sample is heterozygous is below 10⁻³. The cost then scales with the number of
plausibly heterozygous genotypes. *Unmeasured*: the share of sample-record pairs above such a floor.
At 30× and above it is close to the share of heterozygous calls: in human data, a sample is
heterozygous at roughly 4 records in 100 of a large cohort's records, which would bring the
3,000-sample, 30× window to about 1.1 GB. **At 3× the floor saves much less**: a sample showing three
reference reads at a record whose alternative frequency is 0.05 keeps a probability of being
heterozygous of about 0.095 × ⅛ ≈ 0.012, twelve times the floor. At 3× the floor bites mainly at rare
variants.

**The capture in the merge** is transient and capped at 32 per sample (§2.2). The streaming driver
holds one locus at a time. The round driver holds a round's called records until the round folds back;
above about a thousand samples it runs one region at a time (`round_shape_for`,
`src/cli/calling_run.rs:456`, with its table at `:440-449`), so at 3,000 samples the round driver does not run and only the
ordered thread is a bottleneck (§4.2).

**The spill.** Pass one parks every called record, and pass two reads them all to score the
duplication filter (`src/run/paralog_filter/mod.rs:9-15`). Each entry today is the line plus, per
sample, a window and two counts (`SpilledSample`, `spill.rs:150-161`), 16 bytes. Fragments are put in
a separately skippable section of the entry, so pass two skips them and only pass three reads them.
Whole-genome arithmetic at 3,000 samples and 30×, before compression, at a record every 35 bases
(about 88 million records; 1000 Genomes has about 85 million variants for 2,504 samples):

- the spill today: a line of about 150 kB plus 48 kB of sample rows, about 200 kB a record, about
  17 TB;
- fragments without Q1's floor: 3,000 × 30 × 6 bytes = 540 kB a record, about 48 TB more;
- fragments with the floor, at 4 heterozygous samples in 100: about 22 kB a record, about 2 TB more.

At that scale the spill is already beyond any ordinary disk, so the comparison that matters is
relative: fragments triple the spill without the floor and add about a tenth with it.

### 4.2 Time

- **Capture**: per read per locus, on the workers. *Must be measured* at 200 samples before the build
  commits to capturing at every cohort locus (§2.2).
- **The fragment keeper**: one hash-map lookup per read allele, on the ordered thread. At 3,000
  samples and 30× without the floor, a record carries 90,000 read alleles; at tens of nanoseconds
  each that is a few milliseconds a record, and at a record every 35 bases a few minutes of one thread
  a megabase. That would be the run's bottleneck. Q1's floor divides it by about 25; the fan-out over
  sample chunks (§2.4) divides it by the thread count. *Must be measured* at 200 samples first.
- **The link sums**: per fragment per tested link, a handful of floating-point operations.

---

## 5. Removing the `refine-genotypes` command

**Owner, 2026-11-03: refinement runs only inside a calling run.** The command read only VCF lines,
and fragments are never in a VCF, so a refinement run from a file could never use read evidence and
would give different output from the same calls.

**What is removed**, as part of this stage's build:

- `src/cli/refine_genotypes.rs` and `src/cli/refine_genotypes/tests.rs`, the `RefineGenotypesArgs`
  and `RefineGenotypesCliError` exports (`src/cli/mod.rs:23`, `:44`), and the subcommand
  (`src/main.rs:29`, `:61`; `src/cli/command_line.rs:13`);
- the refusals that existed only for it: an input already refined, and an input already phased
  (phase output spec §3, last bullet);
- `VcfWriter::create_with_header` as a separate entry point (`writer.rs:128-171`); its only other
  caller is `create` (`:119`), so it can fold into it.

**What is renamed:** `--neighbour-evidence {auto,on,off}` becomes `--refine-genotypes {auto,on,off}`
on both calling commands (`src/cli/call_from_alignments.rs:272`, `src/cli/call_from_psps.rs:252`),
with the same values and default. `--phase-min-path-agreement` stays on both. In prose the step is
still *refinement* (glossary); "refine-genotypes" is only the option's name.

**What else changes:**

- **Tests** in `refine_genotypes/tests.rs` that check the pass's behaviour through the command move
  to the stage's own tests, which feed lines to `NeighbourEvidenceStage` directly, or to a calling
  run on a small fixture. Tests of the refusals and of header rewriting are deleted.
- **Docs**: `doc/user/neighbour_evidence.md`, the glossary's *Refinement* entry (which names
  `pvc refine-genotypes`), `neighbour_evidence_pass.md` (its ruled question on a stand-alone command, line 492)
  and the phase output spec's §3 are updated to the option. `neighbour_evidence_phasing.md` names
  neither the command nor the option and needs no change.
- **Research scripts**: `founder_grid.py` passes `--neighbour-evidence` (`founder_grid.py:154`) and
  must pass `--refine-genotypes`. `compare_templates.py`, `validate_fill_check.py` and
  `tomato_fill_test.sh` in `doc/devel/ng/research/neighbour_evidence/` run the command on an existing
  calling VCF and must instead rerun calling with `--refine-genotypes on`. That costs a calling run per
  measurement instead of a refinement of a file; on the founder simulator's cohorts that is minutes.

---

## 6. Traps

- **A chain id means nothing outside its sample.** Every structure here is per sample; never compare
  ids across samples (`build.rs:2248`).
- **One id can name two sequences at one minted record** (overlapping mates that disagree). Drop it
  there (§2.2); counting it twice would let one molecule vote both ways.
- **Allele numbers change four times** (§2.3): `moved_to`, `AlleleRemap`, the drop of uncalled
  alternatives, and the pass's subsetting. Missing one sends a read's vote to the wrong allele
  silently. Each has a test (§7.4).
- **The list must be sorted by id** before deduplication, the disagreement check and the cap, and
  fragments summed in a fixed order, or two runs can differ in the last bit of `L` near a threshold.
- **The cap's key is the rank in the walked region, never the raw id** (§2.4). A raw-id key passes
  every test whose fixture has fewer than 33 reads a sample at a record.
- **The line carries no fragments.** Anything that rebuilds a record from its line, as pass three
  does, must take the fragments from the spill entry beside it. A line written without them is phased
  by the population alone with no error raised.
- **Same-`POS` records**: the legal tie is a generic record and an STR record at one position
  (`writer.rs:174-186`). The STR record carries no read alleles, so fragments never meet a tie.
- **A wrong heterozygote gives a confident link of random sign.** The `GQ` ≥ 30 floor (§2.5) is what
  keeps it out, for both writers. `LEAN_MIN` does not.
- **`e_floor` covers independent errors only** (§2.5). Do not raise it to cover a duplicated region;
  that is the hidden-duplication filter's job.
- **An error above ½ reverses a vote** unless clamped (§2.5).
- **The keeper's hold is in bases, not records.** A contig end releases everything held; a contig
  change clears every sample's fragment numbers.

---

## 7. How we know it works

1. **The ceiling, and how many links it rests on.** On every cohort, links made or confirmed by reads
   are at most 1 wrong in 100 inside phase sets, and population-only links are no worse than
   `phase_output_2026-11-02.md`. Reported separately for links from the population only, reads only,
   both agreeing, and population links broken by reads (with how many of those were right). Each row
   states its link count: with no wrong link, about 300 links are needed to show a rate under 1 in 100
   with 95% confidence, and a row with fewer is reported as not decided. **`GQ` calibration**: on the
   heterozygotes that read-made links touch, at 3× and at 30×, the observed share of wrong genotypes
   against the share `GQ` states, binned by `GQ`.
2. **Fit `e_floor` on deep data before tuning anything else.** On GIAB HG002 at 30×
   (`benchmarks/giab/`), at records that are heterozygous in the GIAB truth, for pairs reached by three
   or more fragments, the share of fragments that disagree with the majority is the empirical
   per-fragment error. Restricting to truth heterozygotes keeps genotype errors out of the fit: at a
   false heterozygote about half the fragments disagree. `e_floor` is the error that makes §2.5's model
   match it. Then score reads-only phase against a truth phase: the GIAB truth VCF if it carries phased
   genotypes, which has to be checked first, or else the trio HG002/HG003/HG004 called together, whose
   Mendelian phase is an independent truth. Repeat at 3× and 10× by downsampling, and at 300× if those
   data are obtained.
3. **The founder simulator**, which writes paired 2 × 150 reads from 350-base fragments and a phased
   truth VCF (`examples/ng_str_stutter_benchmark.rs:15-30`, `:143-147`). It needs one addition: **a
   setting for the SNP rate**, now fixed at 10⁻³ (`:189`), so that a cohort can have tomato-like dense
   heterozygosity, for example 10⁻². The simulator misreads bases exactly as their quality says and has
   no duplicates, chimeras, paralogs or indel errors (`:103-106`), so it measures reach and how reads
   combine with population links, never whether `e_floor` holds. Then:
   - the nine cohorts of `phase_output_2026-11-02.md` as they are, with refinement on: the error rates
     of item 1, and the change against that report's table;
   - the outbred cohorts again at the dense SNP rate: whether read-made links meet the ceiling there,
     how many wrong population links reads break, and the read check's estimate (§2.8) against the true
     error of the population links;
   - single samples at 3×, 10× and 30× with refinement off: reads-only phase against truth;
   - scoring with `phase_vcf_score.py`, which already reads `PS`.
4. **Unit tests**: capture on the one-record and several-record branches; an id listed twice with one
   allele kept once; an id with two alleles at one minted record dropped; renaming through `moved_to`,
   through a `None` remap, and through `drop_alternatives_no_sample_calls`; a generic record absorbed
   by an STR record losing its read alleles; the cap keeping the 32 smallest keys, and the same reads
   whether the region is walked alone or after other regions; the keeper keeping a read seen at two
   records within `R` and dropping one seen once and one seen beyond `R`; `w` and `L` on hand values at
   the thresholds; `e` clamped at ½, so a Phred-0 read allele weighs exactly zero and one at Phred 2
   does not reverse a link; each row of §2.7's table; **a pass-called heterozygote whose reads all show
   one allele is never read-linked**; a record with `GQ` 29 never read-linked; a read link flipping the
   following block; a chain of two non-leaning records attached in order, one refused because its link
   to `q` is at −`T_break`, and a non-consecutive fragment ending a block; `PS` of a block whose first
   record is attached; the read check's counts; fragments through the spill and pass three, with pass
   two skipping them.
5. **Invariance**: the same bytes at 1, 4 and 30 threads; the round driver writing what the streaming
   driver writes (the rule `psp_caller.rs:780-784` states); the hidden-duplication filter on and off
   writing the same phase on a fixture where it filters nothing. **Direct and psp mode** write the same
   phase from the same walked regions, on a fixture with more than 32 reads per sample at the linked
   records, so the cap binds and a key that depended on walk history would show. The same comparison is
   made between a whole-contig run and a `--regions` run over one walked region. The captured
   `(rank, allele)` lists are also compared directly between the two modes.
6. **Cost**: the capture, the keeper, the pass's held fragments and the spill, at 50 and 200 simulated
   samples and on tomato at 3×, with and without Q1's floor; extrapolated to 3,000 samples with the
   assumption stated. The capture's share of the merge time decides §2.2's condition. The ordered
   thread's time is reported on its own.
7. **Real data**:
   - tomato (`benchmarks/tomato1`, 63 accessions): the read check of §2.8; deep accessions
     read-phased with refinement off, compared with refinement on;
   - the share of reachable heterozygote pairs lost to STR segments (§2.1), on tomato and GIAB psps,
     since the simulator's STR density is its own;
   - correlated errors, which the simulator cannot produce: on GIAB HG002, links inside annotated
     segmental duplications against links outside them, with the hidden-duplication filter on and off.

---

## 8. Deferred, with a recommended home

- **One id per read across walked regions.** The walker now lives for a chromosome
  (`generator.rs:1355-1366`), so pending mates could survive a region boundary and a read spanning an
  STR segment would keep one id. Home: after §7.7 measures what the split costs.
- **A per-read quality.** The mean hides one bad read among three, and at 3× that read may be a link's
  only witness. Home: when a link's weight is shown to matter at 3×.
- **`PQ`**: `L` gives one directly for read-made links. Home: with the phase output spec's `PQ` item.
- **Writing the reads' relation over a population link they confidently contradict** (§2.7). Home:
  after §7.3 measures how often such contradictions are right.
- **The pass's sampler using reads**, as SHAPEIT2 multiplies a read likelihood into its draws
  (research note §3.1), so reads would shape the haplotypes and not only the links. Home: if §7.3
  shows that changing links is not enough on outbred cohorts within reach.

---

## 9. Questions for the owner

The starting values of §2.5 (`T_link` = ln 199, `T_break` = ln 9, `e_floor` = 0.01, `GQ` ≥ 30) are
in the text and are checked by §7; they need no ruling now.

**Q1. Drop read alleles where the sample is almost certainly homozygous?** Owner decision 2 keeps
reads whatever the *called* genotype, because refinement can make a homozygote heterozygous. This
proposes a floor on the *posterior*, not on the call: keep a sample's read alleles at a record only
where calling's probability that the sample is heterozygous there is at least 10⁻³. **Recommend:
adopt it.** Without it, at 3,000 samples and 30×, the pass holds about 28 GB of fragments a window,
the spill triples, and the ordered thread does 90,000 lookups a record (§4); with it, all three follow
the plausible heterozygotes. **The saving is at 30× and above, and at rare variants.** At 3× a sample
with three reference reads at a common variant stays above the floor (§4.1), so the 3,000-sample, 3×
case stays near 2.8 GB a window either way. What the floor gives up: a sample whose reads say
homozygous with odds over 1,000 to 1 that refinement nonetheless turns heterozygous loses read
evidence there. §7.6 measures both sides.

**Q2. `R` = 1,000 bases, and the cap of 32 read alleles per sample per record?** **Recommend: both,
as starting values.** `R` covers standard inserts; a long-insert library loses links beyond it, and
raising `R` costs only held records. The cap keeps 300× samples from costing ten times what 30×
costs; WhatsHap uses 15.

**Q3. Should the read check (§2.8) tighten population-only links by itself?** The phase output spec
asks for "a signal that tightens on its own where the phase is poor, never by excluding the cohort".
Where reads find, say, more than 1 in 100 confident population links contradicted for a sample,
population links that no read confirms could be refused for that sample. **Recommend: report the check
now and decide after §7.3.** The counts per sample are small at 3×: showing a rate under 1 in 100 with
no contradiction takes about 300 reached links. The rate is also biased towards short links (§2.8).

**Q5. Phase from reads alone in a one-sample run?** It follows from decision 3 and costs little.
**Recommend: yes.** It is the one place a single deep sample gets any phase. At 3× it writes little
(§2.6). (Q4, on the thresholds, is folded into §2.5 as starting values.)
