# ng — the alignment cursor: implementation plan

*Draft, 2026-08-02. Turns the settled design —
[spec](../spec/alignment_cursor.md) (what and why) and
[arch](../arch/alignment_cursor.md) (types and interfaces) — into build order. **No design
happens here.** If a step turns out to need a decision, stop and take it back to the spec.*

---

## Scope

**In.** Replacing the per-region read query with a long-lived cursor: a per-file reader that stays
positioned, keeps the reads it has already decoded and filtered, and hands them back when the next
region can use them. The sample-level cursor that merges k files. The walker and generator changes
that follow from there being no per-region stream.

**And moving *every* BAM/CRAM read in ng onto it.** When this plan is done there is one way to read
an alignment file and the old one is gone: `SampleReads::reads_in_region`, `RegionReads`, the
reader pool and `readers_opened` are deleted, and nothing references them. That is Milestone F, and
it is not a tidy-up — it is 13 files beyond the two generators, inventoried there.

**Out.**

- **The per-chromosome reference-base registry** — deferred whole (spec §12). Its trigger is the
  first parallel run over CRAM. A cursor takes its bases once at construction until then.
- **The parallel fan-out.** This plan builds the shape it needs and runs single-threaded.
- **Coalescing regions at the caller.** Refuted — generic regions are not adjacent (spec §13, N2).
- **Removing the per-read copy into the caller's buffer** — measured at +2.8 % and not worth a
  redesign before there is a working baseline (spec §13).

## Principles (how the order was chosen)

- **The algorithmic heart before the plumbing.** The forget rule (spec §6) is one comparison, and
  it is the only part that can silently lose reads. It is built and tested against an in-memory
  reader before any file is involved.
- **Simplest implementation first, as the oracle for the next.** `RecordReader::InMemory` exists
  so the rule can be driven from a scripted list of records. BAM follows, CRAM last.
- **Verify against ground truth, not self-consistency.** Every milestone is proven by output
  identity against the current code on real data, not by its own tests. **1,471 unit tests passed
  while the first attempt lost 3,830 loci** (spec §11).
- **Isolate the step whose failure is silent.** Retention that drops a read it should have kept
  produces a wrong genotype, not a crash. Those steps land as their own commits.
- **Types first, then implementation**, within every milestone.
- **Incremental, with pauses.** One milestone, then stop.

## Preconditions (already in place)

- The design is settled: spec and arch are final, and eleven adversarial-review findings are
  resolved in them.
- **The baseline numbers exist and are reproducible.** `examples/ng_generic_walk_probe.rs` is in
  the working tree; chromosome 21 prints `loci=236081 observations=251786 reads_admitted=54709`
  and chromosome 1 `loci=1541788 observations=1647161`.
- The suite is green at `d95ce8b`: `cargo test --lib ng::` gives 1,471 passed.
- The oracle to extend exists:
  `t5_the_indexed_query_returns_exactly_what_a_linear_scan_returns`.
- **The probe must be committed first** — it is the only instrument that measures the walk rather
  than a tool's output buffer, and every milestone below is verified with it.

---

## The steps

### Milestone A — the instrument, and the types (no behaviour change)

- ✅ **A1.** Commit `examples/ng_generic_walk_probe.rs` with a test module and the fixture builders
  it shares with the bench. *Depends:* —. *Source:* spec §11. — `a400f73`
- ✅ **A2.** `CursorError` (`WrongChromosome`, `Io`), and `AlignmentFile::contigs()`.
  *Depends:* A1. *Source:* arch §1.4, §2.1. — `aaae5db` (`Io` landed as `ReadRecord`)
- ✅ **A3.** `RecordReader` as an enum with only the `InMemory` arm: finds nothing, unpacks
  nothing, yields a scripted list in position order. *Depends:* A2. *Source:* arch §1.3. — `382667b`
- ☐ **A4.** `ReadFilter::source_mut`, and a test that repositioning through it reaches the source.
  One accessor; nothing else in `filtering.rs` moves. *Depends:* A2. *Source:* arch §2.3.

> **Checkpoint A:** the types compile, nothing behaves differently, the probe is committed and
> still prints the baseline numbers. Pause for review.

### Milestone B — the forget rule, against the in-memory reader

- ☐ **B1.** `AlignmentCursor` over `RecordReader::InMemory`: `move_to_region`, `next_read`,
  `contig()`, and the kept reads. *Depends:* A3, A4. *Source:* arch §1.2, §2.2.
- ☐ **B2. The forget rule — its own commit, do not bundle.** Reuse when the new region starts at
  or after the last one served; otherwise drop and reposition. Evict a kept read once it ends
  before the current region's start. **Oracle:** a scripted reader driven through ascending,
  backward, overlapping, adjacent and far-apart regions must return exactly what a linear scan of
  the same scripted list returns. *Depends:* B1. *Source:* spec §6.
- ☐ **B3.** The counters — reads kept, replayed, decoded, repositions — and the test that a
  forward walk decodes each read once. *Depends:* B2. *Source:* spec §11.5.

> **Checkpoint B:** the rule is correct on a reader with no file behind it, and its failure modes
> are exercised before any real input can hide them. Pause for review.

### Milestone C — the BAM arm

- ☐ **C1.** `RegionRecords`: the region narrowing, the sorted early stop, read-group resolution
  and the tally — lifted out of `BamRegionSource`, written once. *Depends:* B2.
  *Source:* spec §5, arch §2.3.
- ☐ **C2.** `BamRecordReader`: index query, positioned reading, and the one-record pushback for
  the record the early stop consumes without yielding. Keeps nothing across regions.
  *Depends:* C1. *Source:* arch §1.3.
- ☐ **C3. Wire the cursor to BAM — its own commit.** **Oracle:** the extended
  `t5_…returns_exactly_what_a_linear_scan_returns`, driving *a run of ascending regions through
  one cursor* rather than a single query. *Depends:* C2. *Source:* spec §11.3.
- ☐ **C4.** `SampleCursor`: k file cursors, the argmin merge, and the `Single` arm kept free of
  dynamic dispatch. *Depends:* C3. *Source:* arch §2.4.

> **Checkpoint C:** BAM reads through the cursor. The probe must print the baseline numbers
> exactly, at 400 bp, 10 kb, 100 kb and whole-contig region sizes. **This is the first point the
> saving can be measured** — record it against 5.18 s on chromosome 21. Pause for review.

### Milestone D — the callers

- ☐ **D1.** The walker holds the cursor for the run instead of taking an iterator at
  construction: `move_to_region` forwarded, and a per-region reset of `WalkerState`, `pending`,
  `done` and `stop_after`. **The largest edit in this plan.** *Depends:* C4.
  *Source:* spec §3.
- ☐ **D2.** The generic generator: cursor per chromosome, minted at the boundary, and
  `make_reference` deleted — which drops a type parameter from `PileupGenerator`.
  *Depends:* D1. *Source:* spec §3, perf review L2.
- ☐ **D3.** The STR generator, same change. *Depends:* D2. *Source:* `ssr.rs:375`.

> **Checkpoint D:** both generators run through cursors; the STR dump and the generic dump are
> byte-identical to their committed baselines. **The old API is still there and still used by the
> callers in Milestone F** — nothing is deleted yet. Pause for review.

### Milestone E — the CRAM arm

- ☐ **E1.** `CramRecordReader`: the `.crai` walk and container decode, keeping nothing across
  regions. The existing single-container cache stays exactly as it is — a within-query
  optimisation. *Depends:* C2. *Source:* spec §5.
- ☐ **E2. Wire CRAM — its own commit.** **Oracle:** the BAM/CRAM parity test (`t8_…`) must still
  show the two formats returning identical reads, now through cursors. *Depends:* E1, D3.
  *Source:* spec §11.
- ☐ **E3.** Measure on a tomato CRAM and record it. CRAM is unmeasured in the perf review; this
  is the first number for it. *Depends:* E2. *Source:* spec §1.

> **Checkpoint E:** both formats read through cursors, with CRAM measured for the first time.
> Pause for review.

### Milestone F — move every remaining reader across, then delete the old path

Until this milestone two ways of reading a file coexist. F ends that. **Every call site below was
found by grep on the clean tree; none is optional, because F4 deletes what they call.**

- ☐ **F1. The acceptance anchors.** `ng_generic_loci_dump` and `ng_ssr_loci_dump` onto cursors.
  Their output is asserted byte-identical, so converting them is itself the proof that the
  migration moved nothing. *Depends:* E2. *Source:* spec §11.
- ☐ **F2. The measurement harnesses.** `ng_generic_walk_probe`, `benches/ng_generic_pileup_perf`,
  `dhat_ng_merge`, `dhat_ng_open_files`. **`dhat_ng_merge` needs care** — it measures the merge's
  own allocation cost by draining one region twice, so it must keep measuring the merge and not
  the cursor's kept reads. *Depends:* F1. *Source:* arch §2.4.
- ☐ **F3. The stage-1 differential.** `parity.rs` calls `reads_in_region` once
  (`:4030-4042`) to feed one read stream to both walkers. It is `#[cfg(test)]`
  (`pileup/mod.rs:131-133`) so it breaks the test build, not the release build — which is why it
  is easy to miss. **First establish whether it still has a job:** the design records that from
  plan 3's A2 the two walkers differ on purpose and this harness "dies by design". If it is
  vestigial, delete it; if not, convert the one call site. **Ask before doing either** — 4,233
  lines is not a decision for a build order. *Depends:* F1. *Source:* `pileup/mod.rs:128-133`.
- ☐ **F4. The research tools.** `ng_ssr_aligner_bakeoff`, `ng_ssr_anchor_firm_validate`,
  `ng_ssr_cohort_stutter`, `ng_ssr_divergent_reads`, `ng_ssr_gain_loss`, `ng_normalizer_screen`.
  Several assert against committed baselines; those must not move. *Depends:* F1.
  *Source:* spec §11.
- ☐ **F5. Delete the old path — its own commit.** `SampleReads::reads_in_region`, `RegionReads`,
  `ReaderHandle`, `BorrowedReader`, the pool, `readers_opened` (**ten read sites across nine
  tests**), and `region_query.rs` itself. Plus the test at `locus_generation/mod.rs:882` and the
  doc link at `ref_seq.rs:611`. **Verification is mechanical:** `cargo build --all-targets` and
  `cargo test` green, and `grep -rn "reads_in_region\|RegionReads\|readers_opened" src/ examples/
  benches/` returning nothing. *Depends:* F2, F3, F4. *Source:* arch §4.

> **Checkpoint F:** there is exactly one way to read a BAM or CRAM in ng. Both dumps byte-identical,
> every example and bench building, the grep clean. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | suite green, probe prints the baseline numbers unchanged |
| B | a scripted reader matches a linear scan of the same script, over ascending, backward, overlapping, adjacent and far-apart regions |
| C | the extended oracle over **a sequence** of regions through one cursor; probe output identical at four region sizes; first saving measured |
| D | both dumps byte-identical to their committed baselines; suite green |
| E | BAM/CRAM parity test green through cursors; first CRAM measurement recorded |
| F | `cargo build --all-targets` green; both dumps and every baseline-asserting tool byte-identical; grep for the old API returns nothing |

**Two things no unit test can prove, so they are checked by hand at every checkpoint:** the probe's
output identity on real data, and the counters showing reads decoded approaching the true count
rather than a multiple of it. All thresholds are absolute counts from one tandem-repeat-targeted
fixture — regression anchors against themselves, not properties of the generator (spec §11).

## Out of scope (next plans)

- **The per-chromosome reference-base registry** — its own change, triggered by the first parallel
  run over CRAM (spec §12).
- **The parallel fan-out** — the plan that uses one cursor per worker.
- **Serving a sorted batch of regions in one sweep** — `alignment_file.md` §3.3 raised it; this
  plan is its prerequisite, not its delivery.
- **One cursor shared by both generators** — halves the kept reads, ties two generators'
  lifetimes together; revisit with the fan-out (spec §12).
