# ng — the chain-id column's encoding: the build order for the experiment

**Status:** plan, 2026-08-19. Nothing built.

This turns the experiment designed in
[`../spec/psp_chain_id_encoding.md`](../spec/psp_chain_id_encoding.md) into build order: what to
write, in what sequence, on which data, and what each step has to produce before the next one
starts. **It is not a place for new design** — §3 names the one question the spec leaves open,
says what must be measured before it can be answered, and hands it back to the spec rather than
deciding it here.

**Most of what this plan builds is throwaway**, and each step says which kind it is. The three
encoders are the exception: if one of them wins, it is the code that ships, so it is written
where `cargo test` can reach it and not inside an `examples/` binary.

---

## 1. What this decides

**Three ways of storing one column, and the question is which one the psp writer should use.**
The column names, for every observation of every record, which reads supported it — one `u64` per
read. Since the ruling of 2026-08-17 it names every read, including the ones that agreed with the
reference, at every position of the genome
([spec §1](../spec/psp_chain_id_encoding.md)); before that ruling it named about 3.4% of them, at
variable positions only. The three candidates are the spec's arms:

| arm | what it stores |
|---|---|
| **A — as today** | a length, then one raw little-endian `u64` per id, the whole block zstd-compressed |
| **A′ — delta + varint** | the first id as a varint, then the gaps between consecutive ids as varints |
| **B — differential** | per record, which reads *started* being named and which *stopped*; every observation but one keeps its own list, and the remaining one's is derived by subtraction |

The output of this plan is a table of measurements and a recommendation. **The owner decides on
the numbers** (spec §7); no threshold is set here.

---

## 2. Scope

**In.**

- A probe that measures the column's shape on real alignments — how many ids per record, how far
  an id reaches, how often it stops and starts again — and answers spec §10 Q4.
- A capture format, so the same real column can be encoded three ways without re-walking the
  alignments each time.
- The three encoders and decoders behind one interface, with a reconstruction oracle.
- The measuring driver, its sweep over block size, and the report.

**Out, and where each goes.**

- **The psp file itself** — byte layout, block cuts, the index, the trailer, versioning — is
  [`run_streaming.md`](../spec/run_streaming.md) §10's, still unwritten. This plan measures a
  column's encoding in isolation and hands the winner over.
- **Everything that needs a real psp file**: write-side peak RSS in a real writer, the cost of
  seeking to a block on disk, and the merge's end-to-end wall time (spec §7, items 2–4). §11
  hands these to the psp writer's own plan, with the arithmetic this plan leaves it.
- **The in-memory shape of `SequenceObservation::chain_ids`** — a narrower id, an inline buffer
  at low depth. Deferred by the spec (§9) to the measurement
  [`cohort_merge.md`](../spec/cohort_merge.md) §8 already owes. Nothing here depends on it, and
  nothing here should be read as answering it: **this plan is about bytes on disk, not about the
  heap the merge holds.**
- **Re-opening the 2026-08-17 ruling**, the local-index alternative (rejected on production's
  experience, owner, 2026-08-18), and a read-major file. All three are closed in spec §8.

### 2.1 Why this runs before the psp writer exists, against spec §7's "when"

Spec §7 says the experiment runs *"after run_streaming.md §10's psp encoding spec lands and ng
has a writer to modify. Before that there is nothing to measure."* **This plan starts now, and
the change of order is deliberate.** Three reasons, and if the owner rejects them the plan waits
without any other change:

1. **Four of the five quantities spec §7 asks for do not need a file.** Compressed bytes are
   exact from an encoder plus zstd; encode and decode cost belong to the codec; the reader's
   resident live set is the reader's. Only the three named in §2's Out list need a writer, and
   §11 keeps them.
2. **Arm B's saving is a function of block size, and the psp encoding spec has to choose one.**
   Arm B restates the whole live read set at the start of every block, so small blocks erase its
   advantage. Run after §10 chooses, this experiment inherits a number it should be informing.
   Milestone D therefore sweeps block size rather than picking one, and the sweep is a result the
   psp encoding spec can use.
3. **If arm A′ wins, the writer gets written once.** Building the writer raw and re-encoding it
   later is the more expensive order.

---

## 3. ⚠ One design question the spec does not answer

**Arm B cannot be coded from the spec as written.** Spec §5 names the gap and does not close it:

> *"Records are not positions, and one sample's records can overlap. […] The live-set stream is
> indexed by position; the id lists belong to records. The mapping between them is the part to
> get right before writing any bytes."*

The unresolved part is what a record's derived list is derived *from*. Arm B stores a per-position
set of live reads and derives one observation's ids as *that set minus the ids the record's other
observations name*. When two of a sample's records cover the same position — a wide deletion
record staying open while shorter records open and close around it
([`open_record.rs:1463-1466`](../../../../src/ng/locus_generation/pileup/open_record.rs)) — a read
can be named by one of them and silent at the other, and subtracting one position's live set from
one record's lists then invents a read the mint never folded. Spec §5's third trap is the same
failure arriving by a second route: a read the depth cap discarded, or one that produced no
observation, is live over the ground and named nowhere.

**This plan does not decide it.** Milestone A measures how big the problem is — how many records
are wider than one base, how often two of a sample's records overlap, and how often a read live
across an overlap is named by one record and not the other. Checkpoint A takes those numbers back
to the spec, which amends §4 before Milestone C3 is written. **C3 is blocked until it does.**

Two facts already in hand bound the answer and belong in that discussion rather than being
re-derived: the generic mint writes a record at **every covered position** — 96,605 records per
sample over 100,000 bases of tomato SL4.0 (commit `bc109406`) — so the great majority of records
are one base wide; and each observation's `chain_ids` are already **sorted ascending and
deduplicated** at the mint
([`open_record.rs:903-904`](../../../../src/ng/locus_generation/pileup/open_record.rs)).

---

## 4. Principles (how the order was chosen)

- **Measure the shape before writing a codec.** Milestone A touches no encoding. It answers spec
  §10 Q4, it produces the numbers §3's design question needs, and it can invalidate arm B's
  arithmetic (spec §4.1) before a line of arm B exists.
- **The simplest arm first, as the oracle for the next two.** Arm A is a port of what production
  already writes; it is built and round-tripped first, and every later arm is correct exactly
  when it reconstructs the lists arm A does.
- **Capture once, encode many times.** The alignments are walked once into a capture file, and
  the three arms are measured on that file. A driver that re-walks the reads per arm would spend
  its time in the generator — 14 to 23 times the merge's cost on the tomato panel (commit
  `bc109406`) — and would make the encoders' differences invisible.
- **Isolate the step whose failure is silent.** Arm B's derived list is the one place a wrong
  answer is not a crash: derive one id too many and the reference observation gains a read that
  does not exist, and the merge composes an allele for it (spec §5). C3 lands as **its own commit,
  not bundled**, with the reconstruction oracle green before and after.
- **Real alignments only.** Both of production's largest wins on this column came from skew that
  synthetic data does not have
  ([`locus_stream_shape_experiments.md`](locus_stream_shape_experiments.md) §5). Every number
  here comes from reads on disk.
- **Both committed corners, at both ends of both axes.** One sample and 63; three reads a
  position and three hundred. A scheme that wins at one depth and loses at another must show it
  (`CLAUDE.md`, *the range, not the example*).
- **Container builds.** All `cargo` through `./scripts/dev.sh` on a machine that has a runtime;
  `cargo` directly on `rick`, where the binary lands in `target/release` (CLAUDE.md).

---

## 5. Preconditions (already in place, confirm before A1)

- **The generic locus generator runs over real alignments and names every read.**
  `PileupGenerator` + `SampleReads` + `LeftAlignPreparer`, driven exactly as
  [`examples/ng_cohort_merge_real_cost.rs`](../../../../examples/ng_cohort_merge_real_cost.rs)
  drives them — the reuse target for A1's driver, so no ingestion code is written here.
- **`SequenceObservation::chain_ids` is populated on the generic path** and empty on the STR path,
  where a locus is one record and `ReadWitness` already answers the question
  ([`locus_generation/mod.rs:289-310`](../../../../src/ng/locus_generation/mod.rs)). The column
  this plan measures is a generic-path column; the report must say so.
- **Production's list codec exists and is the arm-A port target**: `encode_list_column`,
  `encode_list_column_csr`, `decode_list_column`
  ([`src/psp/block.rs:327-380`](../../../../src/psp/block.rs)), together with the LEB128 varint
  helpers beside them.
- **`zstd = "0.13"` is a direct dependency** (`Cargo.toml`), so compressed sizes need no new crate.
- **The fixtures of §6 are on disk.**

---

## 6. Fixtures — the two corners, and the depth ladder that comes free

| fixture | what it is | which axis it moves |
|---|---|---|
| `benchmarks/tomato1/crams/*.bench.cram` (63) with `regions_n160_200kb.bed` | the tomato panel, about three reads a position | cohort size 63, depth 3× |
| the same, first CRAM only | one sample over the same ground | cohort size 1 |
| `benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam` | HG002, one sample, 30× | depth 30× |
| `benchmarks/ssr_hg002/bam/{5x,300x}/…` | the same sample at 5× and 300× | depth, across a sixtyfold range |

The tomato reference is `$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`
(`benchmarks/tomato1/bench.config.sh`), mounted read-only in the dev container (CLAUDE.md).

**⚠ The HG002 fixture is tandem-repeat-targeted, not whole-genome**
([`locus_stream_shape_experiments.md`](locus_stream_shape_experiments.md) §5): its walk types many
more reference bases than it delivers covered ones. That does not distort a *per-record* byte
count, which is what this plan reports, but any per-base or per-second figure taken from it must
say which fixture it came from.

**The depth ladder is not decoration.** Arm B's arithmetic (spec §4.1) is entirely a function of
depth and read length: at depth `d` arm A carries `8d` bytes per record and arm B carries about
`2d/L` entries. Running 5×, 30× and 300× is three more runs of a driver that already exists and
turns that arithmetic into a measured curve.

---

## 7. Instrument and discipline

Carried from the rounds that produced this module's trustworthy numbers
([`locus_stream_shape_experiments.md`](locus_stream_shape_experiments.md) §6):

- **Compressed bytes are the primary result**, and they are exact, deterministic and
  host-independent. Report the chain-id column's zstd-compressed bytes and, beside them, the same
  capture's other columns, so a saving is quoted as a share of the file and not only of itself.
- **Compare after compression, never before** (spec §5). A raw-byte table is inadmissible on its
  own; it may appear only next to the compressed one.
- **Wall-clock comparisons are not admissible on the macOS host** — 6 performance and 12
  low-energy cores. Encode and decode cost is reported as **instructions retired**
  (`/usr/bin/time -l`, floor-subtracted, minimum of three runs a side, arms alternated within one
  script), or as wall time measured on `rick`, and the report says which.
- **The reader's peak resident live set is a first-class result, not a diagnostic** (spec §6): it
  is state arm B adds that arms A and A′ do not have, and at the top of the committed cohort range
  it is `samples × depth` ids held while a block is scanned.
- **Lines of source per arm is a reported column.** The trade is bytes against complexity and the
  complexity side needs a number too (spec §7's *"how much complexity is worth how many bytes"*).

---

## 8. The steps

### Milestone A — the column's shape, before any encoding (throwaway)

**A1. The shape probe.**  ☐
New `examples/ng_chain_id_column_shape.rs`. Drives the generic locus generator over one sample's
reads with the same ingestion, preparation and config as
[`ng_cohort_merge_real_cost.rs`](../../../../examples/ng_cohort_merge_real_cost.rs), retains
nothing, and prints one `key=value` counts header plus the distributions below. `NG_SHAPE_SAMPLES`
and `NG_SHAPE_REGIONS` cap the walk as the sibling probes' knobs do.
*Depends:* preconditions. *Source:* spec §4.1, §5, §10 Q4.

**A2. What it counts.**  ☐
Per sample: records; ids per record (mean and distribution); records wider than one base; pairs of
that sample's records that overlap; **how often a read is named at one record and silent at an
overlapping one** (§3's question); the reach of an id — first and last record naming it — and
**how often an id stops being named and is named again later** (spec §10 Q4, the mate gap); the
per-record share of ids that sit on the reference-matching observation, which is the share arm B's
derivation would remove; and the two unnamed classes, `reads_without_observation` and
`reads_discarded_by_cap`.
*Depends:* A1. *Source:* spec §4.1, §5, §10 Q4.

**A3. Run both corners and write the shape report.**  ☐
Tomato panel at 1 and 63 samples; HG002 at 5×, 30×, 300×. Report to
`doc/devel/ng/reports/chain_id_column_shape_<date>.md`. It must state, in bytes per record, what
each arm is *predicted* to cost from these counts — the spec's §4.1 arithmetic with measured
numbers in it — so Milestone D's measurement can be compared against a prediction rather than
against nothing.
*Depends:* A2. *Source:* spec §4.1, §7.

> **Checkpoint A — pause for review, and the spec is amended here.** Two things leave this
> checkpoint: an answer to spec §10 Q4 (does the stream need a re-entry form, or is every arrival
> alike), and the numbers §3's open design question needs. **The spec's §4 is amended before C3
> is written.** If A3 shows arm B's predicted saving is small at both corners, the plan may stop
> here and recommend A′ without building B — that is a legitimate outcome, not a failure.

### Milestone B — capture the real column once (throwaway)

**B1. The capture format and its writer.**  ☐
Deliberately dumb and uncompressed: a small header naming the fixture, the sample and the walked
regions, then per record its contig, start, end, and per observation a flag for whether its bases
match the reference, its `num_obs`, and its id list. Plus the record's `reads_without_observation`
and `reads_discarded_by_cap`, which arm B needs to not invent reads. Written by A1's probe behind
a `--capture <path>` argument. **No compression and no cleverness** — this file is the experiment's
input, and anything smart in it is a fourth arm nobody asked for.
*Depends:* A1. *Source:* spec §7 (*"the same walk can write all three"*).

**B2. Capture every corner, and commit one small one as a fixture.**  ☐
Captures go to the project-local `tmp/` (CLAUDE.md); one small capture — one tomato sample over a
few kilobases, a few hundred kilobytes — is committed under `tests/data/` as the unit tests'
input, so the round-trip oracle runs without benchmark data on disk.
*Depends:* B1. *Source:* spec §7.

> **Checkpoint B — pause.** Cheap, and the thing to confirm is that a capture read back gives the
> observations the walk produced, record for record.

### Milestone C — the three arms behind one interface (this is the code that may ship)

Home: `src/ng/psp/chain_id_column.rs`, a new `psp` module holding only this. **The home is
provisional**: [`run_streaming.md`](../spec/run_streaming.md) §10 owns the psp module's layout and
may move the file, which costs nothing. What it may not be is an `examples/` binary — the winner
ships, and shipping code is code `cargo test` reaches.

**C1. The interface, and arm A.**  ☐
One trait with an encode side (records in, bytes out) and a decode side (bytes in, the same lists
out), plus the block boundary as an explicit parameter. Arm A calls production's
`encode_list_column_csr` / `decode_list_column` and the LEB128 helpers **as they are**
([`src/psp/block.rs`](../../../../src/psp/block.rs)); no list logic is re-derived. Unit tests on
B2's committed capture.
*Depends:* B2. *Source:* spec §3 (arm A).

**C2. Arm A′ — delta and varint.**  ☐
First id as a varint, then the gaps. **No sorting step**: the mint already emits each list sorted
ascending and deduplicated
([`open_record.rs:903-904`](../../../../src/ng/locus_generation/pileup/open_record.rs)), and a
test pins that so the encoding does not silently depend on an invariant nobody guards.
*Depends:* C1. *Source:* spec §3 (arm A′), §10 Q2.

**C3. Arm B — the differential stream. Own commit, do not bundle.**  ☐
Arrivals and departures per record, one named observation whose list is derived by subtraction,
and the live set restated at every block start so blocks stay independently decodable (spec §4).
**Blocked until the spec's §4 answers §3's question.** Its guard is spec §5's inequality: a derived
list's length must be at most the observation's `num_obs` and at least half of it, since at most
two mates share an id — the same check the walk's differential against production already asserts
([`parity.rs:2245`](../../../../src/ng/locus_generation/pileup/parity.rs)). A violated guard is a
corrupt-input error of the reader (spec §6), never something the merge can see.
*Depends:* C1, Checkpoint A. *Source:* spec §4, §5, §6.

**C4. Arm B, second variant: every observation keeps its list.**  ☐
The same stream with no derivation — the arrivals and departures, and explicit lists everywhere.
This is spec §10 Q3's other half, and having both from one implementation is what lets the report
price the derivation on its own: it is where most of arm B's saving is and where its only silent
failure is.
*Depends:* C3. *Source:* spec §10 Q3.

**C5. The reconstruction oracle.**  ☐
A property test over every committed capture: for each arm, every record's every observation
decodes to **exactly** the `Vec<ChainId>` the capture holds — same ids, same order, since the mint
sorts. Plus a fuzz-style test over generated captures that include the shapes §3 and spec §5 name:
overlapping records, an id that stops and restarts, and a record with reads in both unnamed
classes. **This replaces spec §7's correctness gate** (*"the merge's output must be identical
across the three"*), which needs a writer: identical reconstructed lists is the same guarantee,
one stage earlier and far cheaper to run.
*Depends:* C2, C3, C4. *Source:* spec §7.

> **Checkpoint C — pause for review.** Three arms, one interface, and every arm reconstructing the
> same lists. Nothing is measured yet.

### Milestone D — the measurement (throwaway driver, permanent report)

**D1. The measuring driver.**  ☐
New `examples/ng_chain_id_column_arms.rs`: reads a capture, and for each arm × each block size
reports the chain-id column's zstd-compressed bytes, its raw bytes beside them, the encode and
decode cost by §7's instrument, decode of a single block in isolation (the cost a seek pays once
it has the bytes), the reader's peak resident live set, and the arm's source-line count.
*Depends:* C5. *Source:* spec §7's report list, items 1, 3 and 5.

**D2. The block-size sweep.**  ☐
Sweep the block boundary across a range wide enough to bracket any plausible psp block — the same
discipline as the merge's building-region sweep, which crossed a fiftyfold range (commit
`bc109406`). Arm B restates the live set at each block start, so this sweep is where its saving
either survives or is eaten; arms A and A′ should be nearly flat, and if they are not, that is
itself a result for the psp encoding spec.
*Depends:* D1. *Source:* spec §4 (restating per block), §2.1 of this plan.

**D3. Run every corner and write the report.**  ☐
Both corners, both cohort sizes, three depths, the full sweep. Report to
`doc/devel/ng/reports/chain_id_column_arms_<date>.md`, one table per corner with the same columns
so they can be laid side by side, each measured number beside A3's prediction for it, and a
recommendation with its trade stated in both directions: the bytes saved, and the reader state and
source lines paid.
*Depends:* D2. *Source:* spec §7.

> **Checkpoint D — the owner decides, on the numbers.** The spec deliberately sets no threshold
> (§7). Whatever is chosen, spec §10's questions 1, 2 and 3 are answered in the report and marked
> settled, and the psp encoding spec is told which encoding to write.

---

## 9. Verification summary

| milestone | proven by |
|---|---|
| A — the shape | the counts are printed and add up: every id the walk minted is in exactly one of the classes A2 counts, checked against the generator's own `reads_admitted` / `reads_without_observation` / `reads_discarded_by_cap` counters, which the sibling probes already print |
| B — the capture | a capture written and read back gives the same records and the same lists the walk produced, record for record, on the committed fixture |
| C — the arms | C5's reconstruction oracle: all three arms decode to exactly the mint's `Vec<ChainId>`, on committed captures and on generated ones carrying the awkward shapes; arm B additionally passes spec §5's `num_obs` inequality on every derived list |
| D — the measurement | each measured byte count is compared against A3's prediction; a gap of more than about a factor of two between them means one of the two is wrong and is chased before the report is written |

---

## 10. What could make this exercise worthless

- **Measuring the column against nothing.** A saving of 40% on a column that is 5% of the file is
  not worth a new error class. Every table reports the column's share of the whole capture, not
  only the column.
- **A capture that is not what the writer would see.** If B1's capture drops the two unnamed read
  classes, arm B's derivation is measured against a world where its hardest failure cannot happen.
  B1 carries them for that reason.
- **A block size chosen to flatter arm B.** Hence the sweep, not a number.
- **Synthetic ids.** Chain ids are allocated monotonically from zero across a whole walk
  ([`chain_id_allocator.rs`](../../../../src/ng/locus_generation/pileup/chain_id_allocator.rs)),
  so their magnitudes, their gaps and their clustering within a block are all properties of real
  data. Generated captures are for the correctness oracle only, never for a byte count.
- **Reporting a saving that zstd already delivers.** Arm A′ exists in the experiment precisely so
  arm B is not credited with a saving a few lines also buy (spec §3).

---

## 11. Out of scope — handed to the psp writer's plan

These three of spec §7's reported quantities need a real file and are handed, with this plan's
numbers as their input, to [`run_streaming.md`](../spec/run_streaming.md) §10's encoding spec and
the implementation plan that follows it:

- **Write-side wall time and peak RSS in a real writer.** D1 gives the encoder's own cost; what a
  writer adds around it is the writer's.
- **The cost of seeking to a random block on disk.** D1 measures decoding one block in isolation,
  which is the work; the seek and the I/O around it belong to the file.
- **The merge's end-to-end wall time.** It changes only through decode cost, which D1 measures per
  block — but the composition should be confirmed once, on a real cohort, after the writer exists.

And one item the spec defers and this plan must not be read as answering: **what an observation
costs in memory**, and whether the in-memory id should be narrower — [`cohort_merge.md`](../spec/cohort_merge.md)
§8's owed measurement (spec §9).
