# F1 — `call-from-alignments`: a cohort of CRAMs in, a VCF out

**Date:** 2026-09-01. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone F step F1. **Design:** [`typed_regions_cli.md`](../../ng/spec/typed_regions_cli.md)
(the command surface), [`vcf_output.md`](../../ng/spec/vcf_output.md) (the record format),
[`run_streaming.md`](../../ng/spec/run_streaming.md) §5.1 and §9,
[`arch/run_streaming.md`](../../ng/arch/run_streaming.md) §3.4 and §5.

---

## What landed

**A person can now call a cohort of alignment files from the command line and get a VCF.**

```
pop_var_caller_exp call-from-alignments \
    --reference ref.fa [--catalog ref.fa.repeats.parquet] \
    --alignment a.cram --alignment b.cram ... \
    [--regions analysed.bed] --output calls.vcf.gz \
    (--parameters fitted.parameters.toml | --defaults) [--ploidy 2]
```

Three pieces, in the order they were built.

### 1. The run hands each record over as it finishes it (`src/ng/run/callers.rs`)

`call_cohort_handing_each_record_over(&genotyper, &mut hand_over)`. Same walkers, same merge,
same three calls per locus as `call_cohort` — what differs is that the answer leaves as a
`VcfRecord` the moment it is finished and nothing is kept. It returns a `WrittenCohort`: the
records handed over, the loci called that established no variant, and the two refusal lists and
walk tallies `CalledCohort` already carried.

**`call_cohort` is unchanged and stays the oracle.** It keeps a whole genome of called loci,
which is what an oracle wants and what a command cannot afford; every Milestone D test is written
against it, and the new path is checked against it rather than against itself.

**How the record is built without the locus outliving its evidence.** A record needs what each
sample's reads showed and which of them no *written* allele explains, and both live in the
merge's `CohortObservation` and in candidate selection's leftover — neither of which survives the
call. So `call_one_generic_locus` gained a closure: it is handed the inference, the allele
remapping and the leftover while all three are still in scope, and what it returns is what the
driver keeps. `call_cohort` passes a closure returning the inference unchanged.

### 2. What a record needs that the called locus does not carry (`src/ng/run/records.rs`, new)

The run's half of the seam `ng::vcf::assemble` describes from the other side. Three functions:

- **`evidence_for_output`** — the per-sample `AD` and `DP − ΣAD`, the cohort-pooled mapping
  qualities, the corrected site quality with its two artifact penalties, and the filter.
- **`padding_base_beside`** — the reference base an empty allele is padded with.
- **`a_written_genotype_carries_an_alternative`** — whether the locus reaches the file at all.

**A fourth numbering meets here, and it is not a sample axis.** The merge's allele table and the
record's are different tables, because candidate selection drops alleles between them. Every
per-allele count is keyed by `AlleleRemap` and never by the merge's own index: a read on a
dropped allele belongs in `DP` and in no `AD` slot (spec §7), and using the merge's index would
put it in the wrong one. The fixture that pins this has a **hole** in the middle of its table —
five merge alleles, one dropped, so the survivors take dense ids 1, 2, 3 — because a table
without a hole cannot tell a correct remapping from one that merely counts up.

**The corrected site quality.** `LocusInference::uncorrected_site_quality` is `pub(crate)` and
carried a `dead_code` waiver saying its caller was "the ordered output stage … which is step 11's
plan and not this one's". This is that stage; the waiver is spent and removed, and the compiler
said so on the first build.

### 3. The subcommand (`src/pop_var_caller_exp/call_from_alignments.rs`, new)

One module owning its `Args`, its `run_call_from_alignments` and its `#[non_exhaustive]` thiserror
enum; the variant on `PopVarCallerExpCommand` and the dispatch arm in `src/main_exp.rs`. Shape
copied from `estimate_contamination.rs`.

**`--parameters` and `--defaults` are one clap group, exactly one required.** A group rather than
a pair of conflicting flags, so a run naming neither is told both answers: with
`conflicts_with`/`required_unless_present` the message named only `--parameters`, which the test
caught.

**The run's ploidy is the parameters', not the flag's.** A supplied file states the ploidy its
numbers were fitted at (`parameters_file.md` §3.2), and the records' `GT` has to be enumerated at
the one the model scored at. `--ploidy` therefore applies on the `--defaults` path and is
documented as ignored otherwise.

---

## The padding base — the one thing F1 had to decide

**VCF cannot spell an empty allele**, so an insertion's or a deletion's record is written by
prefixing every allele with the reference base beside its span and moving `POS` one left; at a
span starting at a contig's first base the base to the right is appended instead, `POS` unmoved
(spec §5). **The locus does not carry that base** — the merge gathers the bases of the span and
no more.

**Decision: the run holds one reference accessor of its own for the output, and fetches the base
inside the calling loop at the moment a record with an empty allele is built.** The accessor is
minted from the same `WalkReference` the walkers' accessors come from, so it shares the parsed
`.fai` and the contig table and has its own cursor; it is minted before `walkers()` consumes the
caller. It reads one byte per such record and calls `evict_before` after each fetch, so its
window does not grow into a resident contig over a genome — loci arrive in genome order, so
nothing asks for a base behind the one just fetched. It costs one more open file, inside
`DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES` (32, of which 8 are named).

**A base that cannot be read stops the run** — `RunError::PaddingBaseUnreadable`, naming the
locus. Production's repeat-tract writer invents the letter `N` at a span it cannot read beside
([`vcf_out.rs:405-435`](../../../../src/ssr/cohort/vcf_out.rs)); spec §5 declines to port that,
on the grounds that a base the reference does not contain, at an unshifted position, is a record
that parses and lies.

**The case with no base on either side is unreachable and is refused rather than special-cased.**
It needs a span that starts at a contig's first base *and* ends at its last — a deletion of an
entire contig, which no aligned read can carry, since a read carrying a deletion must match
reference on at least one side of it. Asked for anyway, the fetch itself returns
`RefSeqError::OutOfBounds` and the run stops;
`a_padding_base_the_reference_cannot_serve_is_refused_rather_than_invented` provokes exactly that
state on a two-base contig.

### ⚑ On today's path the base is never fetched, and the first draft of this report implied otherwise

**The generic mint anchors its indels.** `ReadEvent::footprint_span`
([`decompose.rs:53-54`](../../../../src/ng/locus_generation/pileup/decompose.rs)) gives an
insertion a reference span of 1 — its anchor base alone — and a deletion a span of `len + 1`,
the anchor plus the deleted run. So a five-base deletion reaches the record as `REF ACGTTG`
against `ALT A`: the alternative spells a base and is never nothing. **No allele a generic locus
is called over is empty**, so `padding_base_beside` returns `None` at every locus a run produces
today. `cohort_merge.md` §1.3 says the same thing from the merge's side — *"an insertion's span
is 1 (its anchor base); a deletion's span is the deleted run plus its anchor"* — and
`a_deletion_shaped_the_way_the_generic_mint_shapes_one_needs_no_padding_base` pins it.

**The empty allele spec §5 was written for is the repeat-tract path's full-tract deletion**, and
that path is unbuilt. So this machinery is built, unit-tested against an in-memory reference, and
unreached end to end — which is also why no fixture driving `call_cohort_handing_each_record_over`
exercises it.

**It cannot be left out on those grounds.** `VcfRecord::new` asserts a padding base is carried
**exactly** when some allele is empty, so a run that did not compute one would panic at the first
record that needed it rather than quietly write a wrong one. What F1 built is the answer ready
for the step that needs it.

**One consequence worth recording, because a review predicted the opposite.** A left-padded
record's `POS` moves one base left, and `VcfWriter::check_order` admits a tie only
generic-then-tract — so a padded generic record landing on a preceding generic record's `POS`
would abort the run. That state needs an empty generic allele, which the paragraph above rules
out; when the tract path lands, a padded record is a *tract* record and the tie the writer already
admits is exactly the one spec §5 describes. No change is needed.

---

## Which loci reach the file

**Spec §9: a locus no written genotype carries an alternative at is not written.** There is no
gVCF and no reference block, so the record's absence is the file saying *nothing here*. The rule
is on the calls the **file** would write, not on the ones the loop made: a sample whose reads
said nothing is written `./.` (spec §7.1) and carries no allele into the file, so it cannot be
what keeps a locus in it. A rule reading `genotype` without asking `reads_were_uninformative`
would write a record with an `ALT` no written sample carries;
`a_heterozygote_no_sample_showed_reads_for_does_not_keep_a_locus_in_the_file` is what pins that.

**The count of the dropped ones is in the answer** —
`WrittenCohort::loci_called_but_not_written` — because *called* and *written* differing is a fact
about a run rather than an accident of it. The fixture that reaches it is a sample showing two
different sequences at one position, one read each, with four reference reads: the merge builds
the position on the cohort's **pooled** two non-reference reads, and candidate selection then
asks each sequence separately and drops both at its floor of two, leaving the reference alone.
`candidate_alleles.md` §6.2 puts that at more than one built locus in four on both benchmarks, so
it is the ordinary case rather than a corner.

---

## What was measured

All in the development container at the working tree of this step.

| check | before F1 | after F1 |
|---|---|---|
| `cargo test --lib` | 5,818 passed, 13 ignored | **5,857 passed**, 13 ignored |
| `cargo test --lib ng::run` | 378 passed | **406 passed** |
| `cargo fmt --check` | exit 0 | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 | exit 0 |
| `cargo doc --no-deps` | 26 unresolved-link errors, 23 redundant-link warnings | unchanged |

**39 tests added**: 22 in `ng::run::records`, 6 in `ng::run::callers`, 11 on the command surface.

### The end-to-end check is a real VCF, read back

`a_runs_records_become_a_vcf_a_reader_can_parse` drives
`call_cohort_handing_each_record_over` into a real `VcfWriter` over a temporary file, then reads
the file back: the `#CHROM` line ends with the run's two sample names in the run's own order,
there are as many record lines as `records_written`, their first two columns are
`chr1 15` and `chr1 30` — the contig by **name**, not by index — and every line has nine fixed
columns plus two samples with `FORMAT` `GT:GQ:DP:AD`. What this catches that the record
assertions cannot is a record the *writer* refuses: an order it will not accept, a column it
cannot encode.

### The shapes the brief called out as untested

- **A cohort of one sample** — `a_cohort_of_one_sample_writes_its_records`: two records, one
  sample column each. This is the end of the range `design_principles.md` §0 commits to and the
  shape most likely to have a guard written `<= 1` where it meant `== 0`.
- **An output that will not take a record** —
  `an_output_that_refuses_a_record_stops_the_run_naming_the_locus`: the sink refuses the first
  record, the run comes back with `RunError::RecordNotWritten` naming `contig 0:15-15`, and the
  sink is **not called again** (asserted on a counter).

  **⚑ The sink stops; the walk does not, and that is a real cost.**
  `merge_cohort_handing_each_locus_over` takes a sink that returns nothing, so nothing on this
  side can end the merge: it runs to the end of the analysed ground, decoding every remaining
  read of every sample, before the error is returned. On the fixture that is invisible; on a
  cohort whose disk fills at the first chromosome it is the rest of the genome decoded for
  nothing, and `arch/run_streaming.md` §3.4's "iteration ends at the first `Err`" is not what
  happens. Fixing it means the merge's sink saying *stop* — one `ControlFlow` through both
  drivers and the region builder — which changes the merge's interface and is not F1's.

---

## Deviations from the plan, absorbed

Two, both recorded here and in the plan's F1 entry.

1. **The run gained a second entry point.** The plan's F1 reads as "add a subcommand", and a
   subcommand that drove `call_cohort` would hold a genome of called loci in memory. Arch §3.4
   already gives a caller a stream of `VcfRecord`s, so this follows the architecture rather than
   departing from it; `call_cohort` is untouched.
2. **The run applies spec §9's rule about which loci reach the file, and counts what it dropped.**
   The alternative — handing every called locus over and letting the command decide — would build
   the evidence and fetch a padding base for loci that are then thrown away, and would put an
   emission rule in a place a second command would have to repeat.

## What F1 does not do, and where it goes

- **No parameters file is written beside the output**, so the header's `##parametersFile` line is
  left off rather than naming a file that is not there. That is F2, and `beside_the_vcf` and
  `ParametersFile::write_beside_the_vcf` already exist.
- **The summary the command prints is six counts**, not a run report. F3 is the report: what the
  run refused and why, samples with no reads, every parameter defaulted rather than fitted.
- **Every locus goes down the SNP/indel path.** Repeat-tract candidate selection is specified and
  unbuilt, so a tract in the analysed ground is counted as ground this caller has not built yet.
  A run over tract-rich ground is short, not wrong.
- **Single-threaded**, Milestone E being deferred.

## What the review found and what was done about it

Three agents reviewed the uncommitted step in isolated worktrees: correctness with mutation
testing, design fidelity against the specs, and one whose only job was to read what a person
sees. **The third drove the real command end to end** — a 4.6 kb two-contig reference, its
catalog, two 24× samples in BAM and in CRAM — and `bcftools view`, `query`, `stats`,
`norm -f ref.fa -c e` and `index -t` all accept the output with exit 0, with both planted
genotypes right. That is the strongest evidence F1 works, and it is what the fixed fixtures
could not give.

**The correctness pass ran 19 mutations and killed 14.** All five survivors were missing tests
rather than wrong code — the agent wrote throwaway probes to confirm two of them, and the shipped
code answered correctly in both. Five tests were added to close them, and the sixth gap it named
is closed by the insertion fixture below:

- **the leftover's numbering** — `unmatched` runs parallel to the *merge's* covering list where
  the sample rows run parallel to the *run's* order, and no fixture had both a covering list that
  differed from the run's and a sample with a dropped candidate. Indexing it by the run sample
  understates that sample's `DP` with nothing to notice;
- **`AD` summed over a sample's read groups** — every fixture gave each sample one `@RG`, so
  assigning rather than summing changed nothing observable. On a two-library sample it makes `AD`
  disagree with the pooled mapping-quality count and `VcfRecord::new` refuses the record;
  `read_groups.md` §1 puts that at 157 samples in 1,707 surveyed;
- **the contig-start branch is exactly one position wide** — relaxing `== 1` to `<= 2` passed
  everything, because the empty-allele fixtures sat at 1, 3 and 7;
- **which error wins when both the sink and the merge fail** — unpinned, and the run stops either
  way;
- **discarding the fetched padding base passed all 5,845 tests**, because no fixture drove an
  indel through the command. That is now
  `an_insertion_is_written_as_a_record_that_needs_no_padding_base`: three reads carrying a
  two-base insertion, `REF A` against `ALT ACC`, `padding_base()` `None`. **An insertion rather
  than a deletion, and the fixture reference is why** — every module here shares a hundred `A`s,
  so a deletion sits in one homopolymer and left-alignment slides it off the record, which is what
  D1 measured. Inserted `C`s introduce bases the reference does not have and survive.

**One reading of the spec was taken and the owner may overturn it.** `DP` now includes
`SampleSupport::reads_removed_as_evidence` — reads named at some of a sample's records inside the
locus and not at all of them. Spec §7 says `DP` is *"every read observation the sample had at the
locus, whether or not a written allele explains them"*, and those reads were observed there;
leaving them out makes `DP` understate the depth at exactly the loci that span several of a
sample's records. `SampleEvidenceForOutput`'s own doc names two of the three sources, written
before this seam existed and calling its own shape provisional; `UnmatchedSupport`'s doc excludes
them from the *leftover* because they carry no quality sum, which settles the read likelihood and
not the file. The three sources are disjoint by construction — the leftover is computed from
`supported`, which these reads never reach — and
`the_three_kinds_of_unexplained_read_are_added_together` pins the sum.

**Fixed in this step:**

- **`--ploidy` above 16 panicked**, after the whole cohort had been opened, with a source path
  and a `RUST_BACKTRACE` note. `Ploidy::try_new` turns down zero and nothing else; the read
  likelihood's copy-share table asserts at `MAX_PLOIDY_COPIES`. Polyploid crops are ordinary —
  sugarcane runs to about twelve copies — so the number is now judged before anything is read and
  the ceiling is named. Two tests, one of them checking the refusal is not off by one at 16.
- **`--ploidy` was silently discarded when a parameters file was given.** Spec §3.2 puts the
  ploidy in the file precisely so a run cannot be paired with a different one "without saying
  so". The flag is now an `Option`, so only a ploidy somebody *typed* is compared, and a
  disagreement is refused naming both numbers.
- **`--output` naming a directory, or a directory that does not exist**, was discovered after the
  reference was read and every file opened — and in the directory case after the last locus had
  been called, leaving the in-flight `.tmp` beside the person's directory. Both are visible from
  the path alone and are now refused first. Permission is still the writer's to discover, since
  it is only answerable by writing.
- **A contig longer than a `u32` was narrowed silently** by an `as` cast, which would resolve the
  run's ground against a wrong length. Refused, per the rule `typed_regions.rs`'s own
  `ContigTooLong` records.
- **The summary printed six counts and none of them said how much ground the run could not speak
  for.** Measured by the review: a run over a 60-base `AT` tract at 24 reads a sample printed
  every count as zero and exited successfully — indistinguishable from a clean genome, and the
  opposite of what this module's own header promises. The two ground counts are now printed, and
  the two counts that *are* parts of one total are indented under it.
- **Help text.** The subcommand now says, where a person will see it, that repeat tracts are
  analysed and not called and that the run is single-threaded — both were in module prose no user
  reads. The Markdown asterisks and the `doc/devel/...` paths are out of the flag help; `--catalog`
  now leads with *build it first*; `--alignment` says sample names come from `@RG SM` and are the
  VCF's columns.

**Recorded and not acted on**, each because it reaches past this step:

- **The walk does not stop when a record is refused** — see the note above. It needs the merge's
  sink to be able to say *stop*.
- **Three error messages from the read layer are worse than F1's own.** A file that is not a BAM
  reports noodles' *"failed to fill whole buffer"*; the no-index message names only `.csi` when
  `.bai` also works and never mentions `--build-index-if-missing`, the flag in this very command
  that fixes it; and both chains restate one fact four or five times before naming the file.
- **`GenomeRegion`'s `Display` writes `contig 0`** — its coordinates are 1-based and inclusive,
  as everything in ng is, so it is the chromosome alone a reader cannot match — and it is reachable
  through `RunError::RecordNotWritten` — the one message whose job is to say how far a partial
  file got, and the one a person cannot match against the VCF they are holding. Already on the
  standing list; F1 makes it user-facing.
- **The header declares eight repeat-tract fields this mode cannot produce** (`STR`, `RU`,
  `PERIOD`, `REPCN`, and three tract FILTERs). Legal VCF, and a reader who greps for `STR` will
  conclude the genome has no repeats.
- **Nothing in the file says whether the numbers were fitted or defaulted**, except the echoed
  `##commandline`. F2 fills `##parametersFile`; whether a defaults run should say so in the header
  as well is F2's to settle. There is also no `##fileDate`.
- **Two `--alignment` paths that share an `@RG SM` become one sample silently.** Passing one file
  twice printed `samples 1` and nothing else.
- **Nothing is printed while a run is going.** Single-threaded over a cohort of CRAMs, that is a
  long silence with no way to tell a slow run from a hung one.
- **⚑ GQ is worth a look and is nobody's step here.** On 24 reads a sample, a clean
  homozygous-reference call (`AD 24,0`) came back at **GQ 74** beside a heterozygote (`AD 12,12`)
  at **GQ 99**; at `--ploidy 4` the same reads gave GQ 9, 13 and 29, so a routine `GQ>=20` filter
  would discard nearly every tetraploid call. Measured on the review's own fixture, not on a
  benchmark.

## ⚑ Owed to the owner's documents

**`typed_regions_cli.md` still does not record the four subcommand names.** `generate-psps`,
`generate-census`, `call-from-psps`, `call-from-alignments` were agreed on 2026-08-28 and are
written nowhere; F1 built under them, and
`the_subcommand_is_spelled_call_from_alignments` pins the one it added. Recorded in
`PROJECT_STATUS.md`; the spec is the owner's to edit.
