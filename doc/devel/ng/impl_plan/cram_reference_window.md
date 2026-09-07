# The CRAM decode against a window of the reference — implementation plan

**Status:** draft, 2026-09-07. Branch `ng-cram-window`. This turns the settled design —
[`alignment_cursor.md`](../spec/alignment_cursor.md) §10 *"The CRAM reference bases — a window
per slice, from a reader of the decode's own"* (amended 2026-09-07), with the arch in
[`arch/alignment_cursor.md`](../arch/alignment_cursor.md) §1.1–§1.3 and §4 — into build order.
It is **not** a place for new design. It follows the `ng-cram-perf` branch, merged into `main`
at `52b7b787`, whose research note
[`cram_read_path_2026-09-04.md`](../research/cram_read_path_2026-09-04.md) built and measured
the noodles half (§7) and left the ng half unbuilt (§9, item 2).

**What this buys, in one line.** The CRAM decode stops holding a whole chromosome of the
reference — 198 MB for tomato chromosome 1, 493 MB for human chromosome 1, measured — and holds
instead one buffer per open file the size of the largest slice span that file has met, about
9 kb on a real coordinate-sorted CRAM. **What it costs**, priced while building A1 rather than assumed: one more reference reader per
open CRAM — one `open(2)`, one window of about 9 kb, about 18 µs to build — and the run's
`RLIMIT_NOFILE` refusal split by format to match, because a CRAM cursor now needs three
descriptors where a BAM needs two.

---

## Scope

**In:**

- `decode_container_at` decoding each slice against exactly the span its header declares, with
  the no-window call kept for slices that need no external bases (spec §10, points 1–2).
- `AlignmentFile::cursor` taking the reference **factory** in place of one accessor, and the
  CRAM arm minting a second reader from it for the decode to own (spec §10, points 1 and 5).
- `OpenReference` reduced to the description plus the open-time FASTA check; the repository, the
  one-contig bound, `unbounded` and `bases_for_contig` deleted; `CramAlignedReadsReader` trading
  its `repository` field for that reader (spec §10, point 6).
- The measurement that closes the research note's §9 item 2: resident memory and wall time on
  the real whole-genome tomato CRAM, VCF byte-identical.

**Out (later, with a home):**

- **A shared cohort window for the decode.** Not needed at the slice spans measured on a real
  file; returns only if a real whole-genome human CRAM shows megabase spans. Home: spec §12,
  last bullet — a research note measuring such a file first.
- **Multi-slice containers.** `decode_next_container` already guards against decoding a
  container once per slice, and the guard is untested because noodles' writer puts one slice
  per container. This plan does not change that; it is the fixture gap
  `aligned_reads_reader/cram.rs` documents.
- **The three things broken on `main` before this branch** — the contamination assertion at
  `tests/ng_calling_loop_calls_genotypes.rs:1235` and the two examples that no longer compile
  against `ClosedLocusRanges` (`ng_candidate_selection_probe`, `ng_cohort_merge_real_cost`).
  Not touched here; they are what `main`'s C2d work owes.

## Principles (how the order was chosen)

- **The plumbing lands before the change it carries, so the change can be bisected.** A1 threads
  the accessor down to the decode and changes no byte of output. A2 then switches the decode to
  the window and is the *only* commit in which decoded bases can move. An off-by-one in the window
  offset is a wrong read reconstructed silently, not a panic — which is exactly the kind of step
  the authoring skill says to isolate.
- **Verify against ground truth, twice.** A CRAM and a BAM written from the same records must
  yield the same reads (`t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam`), and on the
  real file the harness's digest over every field of 600,000 decoded reads must stay at
  `a127ecf5083f535a002a5f461f150ede`. Self-consistency proves nothing here.
- **Delete only after the replacement is proven.** The repository goes in A3, once A2's oracles
  are green with the repository already unreachable (an empty one is passed from A2 on).
- **Measure at the end, on the real file, not the fixture.** The research note's §2 showed the
  benchmark fixtures are a different workload on the CRAM read path; B1 uses the whole-genome
  CRAM in the main worktree and nothing else.
- **Container builds.** `cargo` through `scripts/dev.sh`; the real CRAM lives in the main
  worktree, so the container needs `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks`.

## Preconditions (already in place — the executor confirms each)

- **The noodles half is merged and tested.** `vendor/noodles-cram` carries
  `Slice::reference_span`, `Slice::records_over_window`, the `ReferenceSequence::Window`
  variant and the too-small-window refusal (`FORK.md` §5); its 298 tests pass; `main` is at
  `52b7b787` or later.
- **A factory already exists one level up.** `SampleReads::cursor<R, F: FnMut() -> R>(contig,
  make_reference)` calls the factory once per file and passes the result to
  `AlignmentFile::cursor<R: RawRefSeq + ContigTable>(contig, reference: R)`, which probes it with
  a zero-length fetch for the cursor's contig and moves it into `AlignmentCursor::over_records`.
  A1 moves the call site of that factory down one level, from `SampleReads` into
  `AlignmentFile`.
- **Every accessor a cursor is given is already `Send + 'static`.** The production sites pass
  `WindowedRefSeq` (`walker.rs`, `generator.rs`, `ssr.rs`, `parity.rs`) and the tests pass
  `InMemoryRefSeq`. The `Rc`-based `SharedReference` in `parity.rs` is `!Send` but is never handed
  to a cursor — it serves the production walk's `MultiChromRefFetcher`.
- **`RawRefSeq` is object-safe** — `&self` methods with concrete types only — and
  `impl<T: RawRefSeq + ?Sized> RawRefSeq for Arc<T>` exists; `WindowedRefSeq` implements it with
  raw (uncanonicalised) bytes.
- **A fetch behind the window re-reads.** `RawChromReader::prepend_backward` extends the buffer
  backward; `evict_before` is a hint (`ref_seq.md` §6).
- **The fixture CRAM exercises the external-reference path.** `indexed_cram` writes through
  `pileup::per_sample::cram_files::build_cram`, a noodles writer with a reference repository set,
  whose default preservation map has `external_reference_sequence_is_required = true` — so the
  window path, not an embedded reference, is what its slices take.
- **The harness exists.** `examples/ng_cram_decode_layers.rs` has a `--only window` pass that
  fetches each slice's span itself and hashes every field of every decoded record.
- **The real file is reachable.** `benchmarks/tomato_big_cram/DRR000741.p1.cram` and its `.crai`
  in the main worktree; the reference at `$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`
  (mounted read-only by `dev.sh`); contig `SL4.0ch01`.
- **The test bar before step 1** is `cargo test --lib --tests` green except the one
  pre-existing failure named under *Out* — 6,275 library tests and every integration binary but
  `ng_calling_loop_calls_genotypes`. `cargo fmt --check` and
  `cargo clippy --all-features -- -D warnings` are **also red on `main`**, at 29 formatting hunks
  and 3 lifetime-elision errors in `run/cohort_merge/`, none in files this plan touches; do not
  read either as this plan's doing.

---

## The steps

### Milestone A — the decode reads its bases through a reader of its own

**A1. The decode is handed a reference reader of its own, and nothing decodes differently yet.** ✅
`AlignmentFile::cursor<R: RawRefSeq + ContigTable + Send + 'static>(contig, make_reference: impl
FnMut() -> R)` takes the factory in place of one accessor; `SampleReads::cursor` forwards its own
factory rather than calling it. The BAM arm mints one reader, for the cursor, exactly as today.
The CRAM arm mints a second and hands it to `CramAlignedReadsReader`, which owns it as
`Box<dyn RawRefSeq + Send>`. **Both go through one `check_reference_reader` closure** — the
contig-table comparison and the servability probe, the two checks that are about a reader rather
than about the argument — so minting a reader and checking it are one operation. `decode_container_at` keeps its
signature: it still decodes with `records_discarding_tags` and the repository, and gains the
reader as a parameter in **A2**, where the body reads it — a parameter threaded a commit early
and named `_reference` would compile just as well after A2 forgot to rename it, and nothing in
this tree would warn (`aligned_reads_reader/mod.rs` allows dead code module-wide).

**Nothing below `AlignmentFile::cursor` is touched**: `AlignedReadsReader` and
`RegionRawAlignedReads` keep their signatures and their freedom from a reference bound, so
`read/reference_free_first_filter.rs` and the tests that drive `read_next` directly are
untouched. *(A first attempt threaded the accessor down `read_next` instead and had to make that
module construct a reference, which is the one thing its docs say it must not do; discarded
2026-09-07, and spec §10 point 5 now records why.)*
*Depends:* —. *Source:* spec §10 points 1, 5 and 6; arch §1.2, §1.3, §4;
`read_filtering_stages.md` §5.

**A2. The decode fetches each slice's span and decodes over it.** ☐ **Own commit — do not
bundle.** Guarded by the two oracles named under *Verification*, green before and after.
`decode_container_at` gains a `reference: &dyn RawRefSeq` parameter, passed
`&*self.reference_reader` from `decode_next_container`. Per slice: `slice.reference_span()`;
when `Some((contig, start, end))`, `reference.fetch_raw_into(ContigId(contig), start,
end − start + 1, &mut scratch.window)` and
`slice.records_over_window(&scratch.window, start, …)`; when `None`, the existing
`records_discarding_tags`. `RecordScratch` gains `window: Vec<u8>`, cleared and refilled
per slice like its other buffers. The repository argument becomes `fasta::Repository::default()`
at the one call site, so from this commit on the accessor is the only source of bases — a slice
that somehow fell through to the repository would fail loudly rather than decode against a
chromosome that is not there. A fetch error is an `io::Error` on the decode, charged like any
other read fault.

**The existing fixture cannot see the failure this step isolates.** Every CRAM fixture in the
tree is written against an all-`A` reference (`cram_files::build_fasta`), and the in-memory
accessor the cursor tests hand out is all-`A` too (`fixture_reference_bases`) — the two are
byte-identical, which is why `t8` stays valid, but a window fetched one base off still yields
`A`s and still matches the slice's stored digest. So A2 adds, in `container.rs`: (i) a fixture
CRAM written against a **non-periodic** reference — a fixed pseudo-random base sequence from a
seeded generator, through a `build_fasta` variant that takes the bases — decoded through a
`WindowedRefSeq` over that FASTA, with every record's bases, qualities and CIGAR compared field
for field against the BAM twin of the same records; a wrong offset moves the bases and fails
it, and a wrong span fails the digest; (ii) a scripted accessor that returns fewer bases than
the span asked for makes the decode return an error, never a container. The real-file digest
under *Verification* is the second offset-sensitive oracle, and the step is not done until both
are green.
*Depends:* A1. *Source:* spec §10 points 2–4; `FORK.md` §5; research note §7.

**A3. The repository goes.** ☐
`CramAlignedReadsReader::new` and its struct lose `repository`, and `decode_container_at` loses
the argument; `AlignmentFile::cursor`'s CRAM arm stops calling `bases_for_contig`; `OpenReference` loses `bases`, `bases_for_contig`,
`resident_contig`, `bound_to_one_contig`, `unbounded`, the `build_fasta_repository` import, and
the three tests of the one-contig bound. The open-time check in `AlignmentFile::open` — today
`reference.bases()` — becomes a check that the reference names a FASTA and that its `.fai` opens
(`WindowedRefSeq::read_index`), so `CramNeedsReferenceFasta` and the missing-`.fai` fault both
still fire at open (`a_cram_against_a_fai_only_reference_is_refused_at_open` stays green). The
module docs of `reference.rs` are rewritten for what the type now is: a description and an
open-time proof, not a cache. `ReferenceBasesError::Build` keeps its meaning — the FASTA could
not be opened — with the `.fai` read as its source.
*Depends:* A2. *Source:* spec §10 point 6; arch §1.1, §4; `alignment_file.md` §5 (the reuse row).

> **Checkpoint A:** the decode holds no chromosome, every oracle under *Verification* is green,
> and `cargo test --lib --tests` matches the pre-existing bar. Pause for review.

### Milestone B — measured on the real file, and the note closed

**B1. The measurement, and the report.** ☐
Two runs, both against the main worktree's whole-genome tomato CRAM:

1. The harness, `--only window` against `--only contig`, 60 containers from `SL4.0ch01`, seven
   repeats, under `/usr/bin/time -l` on the host build — per-layer seconds, peak resident, and
   the digest, which must read `a127ecf5083f535a002a5f461f150ede` in both passes.
2. `pop_var_caller_exp call-from-alignments --reference … --alignment DRR000741.p1.cram
   --regions <a BED of SL4.0ch01:1-10,000,000> --defaults --threads 1 --output …`, on `main` at
   the merge base and on this branch, three repeats alternated: wall time, peak resident, and the
   two VCFs compared line for line except `##commandline` and `##parametersFile`.

Written up as `doc/devel/reports/implementations/ng_cram_reference_window_<date>.md`, with the
numbers stated against the note's §6–§7 table and the slice spans of spec §10 point 3. The
research note's §9 item 2 is marked done with a link, and its §7 gains one line saying the ng
half is built. Numbers are quoted with the file and the range they were measured on, per
`CLAUDE.md`.
*Depends:* A3. *Source:* research note §6–§7, §9; spec §10 point 3.

> **Checkpoint B:** the report exists, the VCF is byte-identical, and the resident figure is
> stated. Pause for review — this is the point at which the owner decides whether the branch
> merges.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | **Ground truth, and two oracles that can see a wrong offset.** (1) `t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam` in `open_bam.rs`: a CRAM and a BAM written from the same records yield the same reads in the same order — green before A2, green after; **blind to a window offset error**, because its reference is all-`A`. (2) A2's new non-periodic-reference test in `container.rs`: the CRAM decoded through a `WindowedRefSeq` over a pseudo-random FASTA equals its BAM twin field for field — a one-base offset fails it. (3) `ng_cram_decode_layers --only window` on the whole-genome tomato CRAM: every field of 600,000 decoded reads hashes to `a127ecf5083f535a002a5f461f150ede`, the value the contig-resident decode produces. Plus A2's short-window refusal test and A3's unchanged open-time refusal tests. |
| B | The report: peak resident and wall time on the real file, VCF byte-identical over 10 Mb of `SL4.0ch01`, one sample, one thread. |

## Out of scope (next plans)

- **Measuring a real whole-genome human CRAM's slice spans** — the one number spec §10 point 3
  says is arithmetic rather than measured. Home: a research note under `doc/devel/ng/research/`
  the first time such a file is to hand; it decides whether spec §12's shared-window bullet is
  ever built.
- **Locus generation's own cost** — research note §9 item 1, the largest item left in a serial
  run at 32 % of busy CPU, none of it in noodles. Home: its own review and plan.
- **CRAM 3.1 codecs** — research note §4's warning that nothing here was measured on rANS Nx16
  or fqzcomp input. Home: re-run the harness on such a file before claiming anything.
