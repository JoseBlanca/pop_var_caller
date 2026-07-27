# ng read groups — implementation plan

**Status:** draft, 2026-07-27. The build order for making the read group a first-class object:
parsed once, identified run-wide, and carried on every read. Design is settled in
[`spec/read_groups.md`](../spec/read_groups.md) (what and why) and
[`arch/read_groups.md`](../arch/read_groups.md) (types and interfaces). This plan turns that into
order; it is **not** a place for new design. The spec's one open question (which tag spells the
experiment, spec §13) does **not** block: the fallback is defined, so the code is writable now and
closing the question later changes one function.

Follows [`read_input.md`](read_input.md), whose `AlignmentFile` and `SampleReads` this modifies.

---

## Scope

**In:** `ReadGroupId` in `types.rs`; a new `src/ng/read/input/read_groups.rs`; a new
`src/ng/read/aligned_read.rs`; changes to `input/{mod.rs, open_bam.rs, region_query.rs}` and
`read/filtering.rs`; the removal of the per-file sample machinery; and the migration of every
`SampleReads::open` call site (eight examples, two in-source tests).

**Out (later plans / owners):**

- **Splitting an observation's counters by read group** — `ObservedSequence`
  ([`locus_generation/mod.rs:113`](../../../../src/ng/locus_generation/mod.rs#L113)) is untouched.
  This plan carries the identifier as far as the read; a later revision of
  [`locus_generation.md`](locus_generation.md) owns the observation side.
- **An ng-owned prepared read**, so the identifier survives the generic path — spec §12; a
  revision of [`read_preparation.md`](read_preparation.md).
- **Per-read allocation reuse** in the new read type — already parked as measure-first in
  [`read_filtering.md`](read_filtering.md). Unblocked by this work, not part of it.
- **Demultiplexing a multi-sample file in one pass**, and **lifting one-sample-per-open** — both
  spec §12, both deliberately not built.
- **Emitting read-group columns from the stutter dump** — the analysis that motivated this work; it
  consumes the result rather than being part of it.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **Build beside, then switch.** Milestone A adds the read-group table without a single reader
  touching it, so nothing observable can change until B. Each later milestone changes one thing.
- **Behaviour-preserving until the capability step.** Through B and C1 the caller still refuses a
  file that declares several `@RG` — the check only *moves*. The new capability arrives in C2, on
  its own, where it can be tested for what it is. The temporary guard is labelled as temporary.
- **The existing suite is the regression oracle.** Every alignment fixture in the repo declares
  exactly one `@RG`, so the `Sole` path must reproduce today's read stream exactly. That is an
  external oracle, already written, and it covers A through C1.
- **Isolate the silent failures.** Three steps here fail *quietly* — a mis-copied field in the new
  decode, a record resolved to the wrong read group, a foreign read counted as a drop. None crashes;
  each moves a number an error model will later read. Marked **own commit** below, with the oracle
  named and green before *and* after.
- **Reuse over rewrite.** `ReadFilter`, the record sources, the reader pool, the merge and
  `compute_adaptor_boundary` are called as-is; only what the design changes is touched.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); a native host
  build at completion.

## Preconditions (already in place — confirm before A1)

- **`read_input` is complete**: `AlignmentFile` ([`open_bam.rs:72`](../../../../src/ng/read/input/open_bam.rs#L72)),
  `SampleReads` ([`mod.rs:345`](../../../../src/ng/read/input/mod.rs#L345)), the region sources
  ([`region_query.rs:64`](../../../../src/ng/read/input/region_query.rs#L64)), the merge, and
  `AlignmentFileError` ([`mod.rs:59`](../../../../src/ng/read/input/mod.rs#L59)).
- **The code being replaced**: `sample_names()` and the `SampleNames` enum
  ([`open_bam.rs:856`, `:835`](../../../../src/ng/read/input/open_bam.rs#L835)),
  `AlignmentFile::sample_name()` (`:417`), `agreed_sample_name`
  ([`mod.rs:550`](../../../../src/ng/read/input/mod.rs#L550)), and the `MultipleSampleNames` /
  `MissingSampleName` variants (`:132`, `:151`).
- **The decode being copied**: `record_buf_to_mapped_read`
  ([`alignment_input.rs:803`](../../../../src/bam/alignment_input.rs#L803)) and `MappedRead` (`:78`).
  `compute_adaptor_boundary` (`:883`) is already `pub(crate)`; **`cigar_to_ops` (`:1105`) is
  module-private** and must be made `pub(crate)` or copied — decide when writing C1.
- **Test fixtures**: `header(sort_order, contigs, read_groups)`
  ([`test_fixtures.rs:36`](../../../../src/ng/read/input/test_fixtures.rs#L36)) already takes several
  read groups as `(id, sample)` pairs; `indexed_bam` (`:148`) and `read_named` (`:119`) build files.
  A2 extends them.

---

## The steps

### Milestone A — the read-group table (nothing else touched)

**A1. `ReadGroupId` in `types.rs`.**  ✅
`pub struct ReadGroupId(pub u32)` with `get()` and the standard derives, doc-commented as an index
into the run's table. A shared-vocabulary addition, so it lands alone; add it to
`ng_step_interfaces.md` §1. *Source:* arch §1.1.

**A2. Extend the header and record fixtures.**  ✅
`header(...)` gains library and platform per read group; a `read_named_with_read_group` helper tags
a record with `RG:Z`. Nothing consumes them yet — they are what A4 onward is tested with.
*Depends:* A1. *Source:* arch §Test & bench shape.

**A3. Scaffold `read_groups.rs`: the types and the error, no logic.**  ✅
`ReadGroup`, `NameWithOrigin`, `NameOrigin`, `ReadGroups`, `SampleReadGroups`,
`ReadGroupResolution`, `RecordOwner`, and `ReadGroupError` with a doc comment per variant saying
when it fires. Declared and compiled; nothing calls them. *Depends:* A1.
*Source:* arch §1.2, §1.3, §2.

**A4. Parse one header's `@RG` records, with the two hard errors.**  ✅
A pure function over a `sam::Header` returning this file's read groups, or `NoReadGroups` /
`MissingSampleName` — two distinct variants with distinct remedies in the message. Unit-tested
against header literals, the way the existing `@RG` tests are
([`open_bam.rs:1026`](../../../../src/ng/read/input/open_bam.rs#L1026) onward). *Depends:* A2, A3.
*Source:* spec §6.

**A5. The synthesized names.**  ✅
A missing `LB` becomes sample + `@RG ID` + the file name without extensions, marked `Synthesized`; a
missing experiment becomes the library. Tests: the composed value, the origin marker, and that a
declared tag is passed through untouched. *Depends:* A4. *Source:* spec §6.

**A6. `build_read_groups(paths)`.**  ✅
Read every header, apply A4 and A5, assign identifiers in input-file order then header order, group
by sample, and raise `DuplicateSynthesizedLibrary` when two files with the same name in different
directories collide. Tests: identifiers stable across a shuffled `paths` order only in the way the
design specifies; the by-sample grouping; the collision error naming both full paths.
*Depends:* A5. *Source:* spec §5, §6; arch §3.1.

> **Checkpoint A:** the table is built and correct, every existing test still passes, and no reader
> knows it exists. Pause for review.

### Milestone B — open from read groups

**B1. Switch `SampleReads` and `AlignmentFile` to the table.**  ☐
`SampleReads::open` takes `(&SampleReadGroups, &ReadGroups)`, enforces one sample per open over the
read groups it was handed, builds a `ReadGroupResolution` per file and passes it to
`AlignmentFile::open_as`, which stores it. **Temporary guard, removed in C2:** a file whose header
declares several `@RG` is still refused, so behaviour is unchanged. One commit, because a signature
change is atomic in Rust — it includes the ten call sites: eight examples (`dhat_ng_merge`,
`ng_normalizer_screen`, `ng_ssr_aligner_bakeoff`, `ng_ssr_loci_dump`, `ng_ssr_gain_loss`,
`ng_ssr_anchor_firm_validate`, `ng_ssr_divergent_reads`, `ng_ssr_cohort_stutter`) and two in-source
tests ([`locus_generation/mod.rs:614`](../../../../src/ng/locus_generation/mod.rs#L614),
[`ssr.rs:1990`](../../../../src/ng/locus_generation/ssr.rs#L1990)). `ng_ssr_cohort_stutter` loses
its hand-rolled grouping ([`:175`](../../../../examples/ng_ssr_cohort_stutter.rs#L175)) — the
pre-pass replaces it. Every tool asserts exactly one sample and errors otherwise, which is what it
did before. *Depends:* A6. *Source:* arch §3.2; spec §4, §5.

**B2. Delete the per-file sample machinery.**  ☐
`sample_names()`, `SampleNames`, `AlignmentFile::sample_name()`, `agreed_sample_name`, and the
`MultipleSampleNames` / `MissingSampleName` variants — all dead after B1. A pure deletion, so it
lands separately and the diff is readable. *Depends:* B1. *Source:* arch §5 (the four removal rows).

> **Checkpoint B:** the whole suite passes with reads unchanged; the sample name now comes from the
> table and from nowhere else. Pause for review.

### Milestone C — the read carries its read group

**C1. `AlignedRead` and the ng decode.**  ☐ **Own commit — do not bundle.**
`src/ng/read/aligned_read.rs`: the read type with `read_group: ReadGroupId` and without
`source_file_index`; `RecordSourceError`; the decode copied from production's, resolving the read
group through the `Sole` arm. `RawRecord::decode` and every ng consumer move to the new type in this
commit — a return-type change is atomic. **Oracle:** decode the same record both ways and assert
field-by-field equality against production's `MappedRead`, plus the whole existing read-input suite,
whose fixtures are all single-`@RG`. A mis-copied field here produces wrong reads, not a crash.
*Depends:* B2. *Source:* arch §1.4, §2, §3.3.

**C2. The `PerRecord` arm, and lifting the several-`@RG` guard.**  ☐ **Own commit — do not bundle.**
Resolve `RG:Z` per record against the file's map; a missing or undeclared tag is fatal and names the
read. Remove B1's temporary guard. **Oracle:** a new fixture declaring several `@RG` with known
per-record tags, asserting each read's identifier; plus the single-`@RG` suite still byte-for-byte
unchanged, which is what proves the fast path was not disturbed. A record resolved to the wrong read
group is silent — it moves a library's error-model input, nothing more. *Depends:* C1.
*Source:* spec §7; arch §1.3, §3.3.

**C3. The other-sample skip and its tally.**  ☐ **Own commit — do not bundle.**
`RecordOwner::OtherSample` records are not yielded and are counted apart from the drop categories.
**Oracle:** one file declaring read groups for two samples, opened twice; each open sees only its
own reads, the drop categories are untouched, and the two opens' reads reunited equal the file.
Silent failure: a foreign read counted as a quality drop corrupts a drop-rate diagnostic without
failing anything. *Depends:* C2. *Source:* spec §9; arch §1.3.

> **Checkpoint C:** every read carries its read group; a multi-`@RG` file reads correctly; a
> multi-sample file serves two opens. Pause for review.

### Milestone D — per-read-group counts

**D1. Key `ReadFilterCounts` on the read group.**  ☐
`AlignmentFile::counts()` and `SampleReads::counts()` return per-read-group tallies rather than
per-file ones, so a drop rate stays attributable when a file holds several read groups. Update the
counts assertions in the existing tests. *Depends:* C3. *Source:* arch §3.2; spec §8.

**D2. End-to-end fixture.**  ☐
One integration test walking the full stack — build the table from two paths, open both samples,
read a region from each — asserting the identifiers on the reads and the per-read-group counts.
This is the test that would have caught any of the three silent failures at the seam rather than in
isolation. *Depends:* D1. *Source:* spec §11.

> **Checkpoint D:** feature complete against the spec's obligation list. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | unit tests over header literals: the three hard errors, the synthesized names and their origin marker, identifier order, the by-sample grouping. The rest of the suite untouched. |
| B | the existing read-input suite, green with reads unchanged — the sample name moved, nothing else did. |
| C | **C1:** field-by-field parity against production's decode of the same record. **C2:** a multi-`@RG` fixture with known per-record tags, plus the single-`@RG` suite unchanged. **C3:** a two-sample file whose two opens partition its reads exactly. |
| D | per-read-group counts asserted in the existing tests, plus the end-to-end fixture. |

The spec's full test obligation list is **spec §11**; every item there lands in one of the steps
above.

## Out of scope (next plans)

Repeated from *Scope* so nothing is dropped silently: the observation-side split of counters by read
group (a revision of [`locus_generation.md`](locus_generation.md)); an ng-owned prepared read for
the generic path (a revision of [`read_preparation.md`](read_preparation.md)); per-read allocation
reuse ([`read_filtering.md`](read_filtering.md)); single-pass demultiplexing and lifting
one-sample-per-open (both spec §12, unowned until an input needs them).
