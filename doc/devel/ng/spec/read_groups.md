# ng — read groups: making the library visible to every read

*Status: design spec, 2026-07-27. Settled with the owner in discussion; **no code yet**. Sits
between [`alignment_file.md`](alignment_file.md) (one file) and [`sample_reads.md`](sample_reads.md)
(one sample), and changes decisions in both — see §8. Under [`ng_proposal.md`](ng_proposal.md) §1
and [`../arch/module_layout.md`](../arch/module_layout.md). Naming: **STR** in prose, `ssr` in
code. Code-facing companion: [`../arch/read_groups.md`](../arch/read_groups.md) (types &
interfaces).*

*Supersedes two shipped decisions, both recorded in §8: a file whose `@RG` records name several
samples is an error (`alignment_file.md` §3.1 check 4), and the sample name is a single string on
each file handle. Both were built and both go.*

---

## 1. What this is — scope and non-goals

Every read was produced by one library preparation on one sequencing run. The SAM header records
that as a **read group** (`@RG`): an identifier, the sample the DNA came from (`SM`), the library
it was prepared into (`LB`), and the platform. ng reads those records today, extracts one string
from them — the sample name — and throws the rest away
([`open_bam.rs:856`](../../../../src/ng/read/input/open_bam.rs#L856)). From that point on nothing
downstream can tell which library a read came from.

**This spec makes the read group a first-class object**: parsed once per file, given a run-wide
identifier, stamped on every read, and available wherever a read or an observation is.

**Why it matters, and why now.** The STR error model fits two parameters — the per-base error `ε`
and the PCR stutter level — **per chemistry group**, and the reason is explicit in
[`../../specs/ssr_cohort_mark2.md`](../../specs/ssr_cohort_mark2.md) §4.4: *"PCR stutter level and
per-base error are properties of the library protocol and DNA preservation, not of the sequencing
depth."* A library prep is not an individual. We surveyed a 2,085-file, 68-project tomato archive
with [`rick_sample_manifest.sh`](../../../../benchmarks/ssr_tomato1/scripts/rick_sample_manifest.sh)
and the two came apart in 157 of 1,707 samples: 133 samples carry two libraries, 20 carry three,
and four carry 7, 16, 16 and 42. Merging those by sample name — which is what the pipeline does
today — destroys exactly the contrast a per-library error model needs. So the information has to
survive the read path before any of that work can start.

**What this design does:**

- parses every `@RG` record once, at a header pre-pass, and gives each one a run-wide identifier
  (§4, §5);
- stamps that identifier on every read (§7);
- fills in a library name when the file does not declare one, and an experiment name when the file
  does not declare that (§6);
- fails loudly, at open, on the input states that cannot be repaired (§6).

**Non-goals — deliberately excluded:**

- **How the information is used.** Grouping libraries into chemistry groups, estimating `ε` or
  stutter per group, per-library likelihood terms — all downstream, all out of scope. This spec
  only guarantees the information arrives.
- **Judging the metadata.** We trust what the user gives us. A `PL:ilumina` typo, a library name
  that is really a run accession, a library that is one lane of another — none of that is ours to
  second-guess. (In the surveyed archive 94% of `LB` values are `<project>_<run>` strings written
  by a re-headering script. We record them as given.)
- **Demultiplexing a multi-sample file in one pass.** Rejected for now, with the cost accepted and
  written down — §9.
- **Grouping samples for the cohort.** A locus's cohort view is assembled above locus generation;
  see §9 and [`locus_generation.md`](locus_generation.md).

**This design does not** change the merge, the region query, the index handling, the filtering
policy, or the order guarantee. It widens a channel that already exists: the per-file tag
`MappedRead.source_file_index` already travels from the record source to the decoded read, along
exactly the path a read-group identifier needs.

---

## 2. Where it sits

```
input paths ─▶ §5 header pre-pass ─▶ ReadGroups (by identifier, and grouped by sample)
                (headers only)              │
                                            ▼
                            per sample: SampleReads::open(its read groups)
                                            │
                     file ─▶ AlignmentFile ─▶ region source ─▶ §7 stamp ─▶ ReadFilter ─▶ merge
                                                                              │
                                                                              ▼
                                                         one ordered stream of ng reads,
                                                         each carrying a ReadGroupId
                                                                              │
                                                                              ▼
                                                            locus generation ─▶ observations
```

The pre-pass is new. Everything to its right exists and changes only in what it carries.

---

## 3. What a read group is, in code

A read group is the record itself plus the two names we may have had to invent. `file`, `id` and
`sample` are read from the input and never modified. Library and experiment names are grouping
keys computed from those atoms when the file does not supply them. Because the atoms survive
alongside them, no naming choice we make here can destroy information: a consumer that dislikes
our grouping can always regroup on `(sample, id)` or on the file.

```rust
/// One `@RG` record, or the one read group a file with a single record has.
pub struct ReadGroup {
    /// The file that declared it. `Arc<Path>` matches what `AlignmentFile` already holds.
    pub file: Arc<Path>,
    /// `@RG ID`, verbatim. A label and never an identity — see §4. SAM makes it unique
    /// within its file and says nothing across files; **this caller requires it unique
    /// across the whole run** and refuses a repeat (§6, the owner's ruling of 2026-09-04).
    pub id: Box<str>,
    /// `@RG SM`. Required; absence is a hard error (§6).
    pub sample: Box<str>,
    /// `@RG LB`, or synthesized (§6).
    pub library: NameWithOrigin,
    /// The sequencing experiment. Falls back to the library (§6).
    pub experiment: NameWithOrigin,
    /// `@RG PL`. Carried for reports only — nothing keys on it (§6).
    pub platform: Option<Box<str>>,
}

/// A name used for grouping, plus whether the file declared it or ng made it up.
pub struct NameWithOrigin {
    pub value: Box<str>,
    pub origin: NameOrigin,
}

pub enum NameOrigin {
    Declared,
    Synthesized,
}
```

`NameOrigin` costs one discriminant and buys the ability to say, in any later report, whether a
grouping came from the file or from us. That distinction is what keeps a chemistry-group diagnostic
honest, and it cannot be recovered afterwards.

---

## 4. Identity — one integer, run-wide

**A read carries one number.** `ReadGroupId` is an index into a flat table of every read group in
the run; one lookup answers both "which library?" and "which file?", because the `ReadGroup` holds
its own path.

```rust
pub struct ReadGroupId(u32);

/// Every read group in the run, in assignment order. `ReadGroupId` indexes it.
pub struct ReadGroups { /* Vec<ReadGroup> */ }
```

**Not the `@RG ID` string.** Two reasons, and the second is the one that would bite at runtime.
`ID` is unique only within its file *by the SAM specification*, so keying on the string would risk
fusing two unrelated read groups. And a string identifier means an `Arc<str>` or a hash lookup per
read, on a loop that carries millions of reads through the region queries the whole read-input cost
model was built around ([`alignment_file.md`](alignment_file.md) §3.3). A generated integer cannot
collide and costs nothing.

**⚠ Since 2026-09-04 the first reason is a rule rather than a hazard**: a run whose files declare
one id twice is refused outright (§6). The generated integer stays, for the second reason and
because it is what the psp records, but nothing now depends on it to keep two same-named lanes
apart — that shape does not reach the caller.

**Scoped to the run, not to the sample.** An earlier draft scoped the identifier space to one
`SampleReads`, with each file receiving a base offset. It fails on the case §9 requires: a file
holding two samples is opened twice, so a per-open numbering would hand the *same physical read
group* a different identifier in each open — exactly what the identity property below forbids. The
run-wide space has one owner and one number per `@RG` record, whatever opens it.

**The identifier space is run-wide; an open is not.** `SampleReads` stays a **single-sample**
object, and that is enforced rather than merely arranged: every read group handed to
`SampleReads::open` must name the same sample, checked before any read flows. A file whose read
groups name several samples is opened once per sample (§9), and each open is given only its own
sample's read groups — the other samples' reads are not filtered out late, they were never that
open's to begin with.

So the contract of the whole read layer is unchanged: **a file, or a set of files, is opened to get
the reads of exactly one sample.** What changes is only that a *file* no longer has to belong to
one sample for that to hold. The check still earns its keep even though the pre-pass (§5) groups
read groups by sample and therefore never trips it: tools and tests assemble an open by hand, and
this is the guard that stops a foreign file being read as part of a sample.

**Not a store that files register themselves with as they open.** That was the other rejected
shape, and it fails on a case this design has to support: the same file opened twice for two
samples (§9) would register its `@RG` records twice and mint two identifiers for one physical read
group. Registration must be keyed on the file and must happen once — which is what §5 is for.

**Identity is per `(file, @RG)`; grouping is by tag value.** Two files declaring `ID:1` are always
two read groups with two identifiers. Whether they should be *treated* as one thing is a question
about their library names, answered later by whoever cares. Keeping the two apart is what stops an
ambiguity in the data becoming an ambiguity in the identifier.

**A tool that ignores the run-wide space is caught downstream, not here.** The pre-pass promises
identifiers unique across *the paths it was given*, and a caller that hands it one file at a time
gets nothing from that promise: every sample's first read group comes back as identifier `0`. The
pre-pass cannot see this — it is one call per file, each correct on its own — so the refusal belongs
where the samples meet. A cohort refuses two samples that claim one identifier, naming both
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §7 check 10); a cohort tool
must call the pre-pass once over every alignment and then open one `SampleReads` per entry of the
by-sample view, which is what the run-wide space is for.

**Assignment order is fixed and does not depend on scheduling**: input-file order, then header
order within a file. Files may still be opened concurrently; the identifiers are assigned in a
serial pass. Two runs over the same input list produce the same identifiers, which is what lets any
identifier reach a report or an output file.

---

## 5. The header pre-pass

Before any sample is opened, read every input file's header — headers only — and build the table:

```rust
pub fn build_read_groups(paths: &[PathBuf]) -> Result<ReadGroups, ReadGroupError>;

impl ReadGroups {
    /// One read group, by its identifier.
    pub fn get(&self, id: ReadGroupId) -> &ReadGroup;
    /// The same read groups, grouped by the sample they name — one entry per sample,
    /// in first-seen order. This is what each open is built from.
    pub fn read_groups_per_sample(&self) -> &[SampleReadGroups];
}

/// One sample and the read groups that name it.
pub struct SampleReadGroups {
    pub sample: Box<str>,
    pub read_groups: Vec<ReadGroupId>,
}
```

The by-identifier and by-sample views are the **same read groups**, not two collections, so they
live on one type.

The pre-pass parses each file's `@RG` records, applies the §6 rules, assigns identifiers, and
groups the result by sample name.

**This is not extra work bolted on for an edge case — it is the step that decides what the samples
are.** You cannot build a per-sample file list without reading the headers, and ng's tooling
already does this by hand: [`ng_ssr_cohort_stutter.rs`](../../../../examples/ng_ssr_cohort_stutter.rs)
opens each file once just to probe its sample name, then groups the paths itself. The pre-pass
replaces that with one shared step.

**Cost.** Each file's header is read twice — once here, once at its real open. A header read is a
file open and one block; the expensive part of opening, parsing the index
([`open_bam.rs:82`](../../../../src/ng/read/input/open_bam.rs#L82)), happens only at the real open.
Not measured; expected to be lost in the noise. If it ever isn't, the pre-pass could hand each
parsed header forward, at the cost of holding every header in memory at once.

---

## 6. What is a hard error, and what we fill in

Three states are hard errors, all detected before any read flows, and all reported with the file
and the remedy. One is lenient by design.

| invariant | where checked | on failure |
|---|---|---|
| the file declares at least one `@RG` | pre-pass | error naming the file, and saying it must be re-headered |
| every `@RG` carries `SM` | pre-pass | error naming the file and the offending `@RG ID` |
| **no two `@RG ID`s in the run are equal** | pre-pass, when the files are opened | error naming the id, both samples, both files, and `samtools addreplacerg` |
| every record names a declared read group — **only in a file that declares several** (§7) | first read of the record | fatal to the run, naming file, read and position |
| `LB` present | — | absent → synthesized, below |

**The uniqueness rule is the owner's, 2026-09-04, and it refuses input the SAM specification
allows.** Two files may legitimately each declare `ID:1` — one sample's two lanes, or two samples
each aligned by a pipeline that names its read group the same way. Within a sample the collision is
silent and costly: one lane's reads counted as another's library. Across samples nothing merges,
since a lane's identity here is the generated integer. It is refused in both shapes anyway, because
in practice a repeated id is a mistake — the same file passed twice, or read groups copied between
files — and because a run whose lanes cannot be told apart by the name they carry is one whose
provenance nobody can follow afterwards: every report, every parameters file and every error
message names a lane by its id.

**The same rule is held over stored files**: a cohort of psps is refused for a repeated id exactly
as a cohort of alignment files is, so psp mode cannot call what direct mode turns down
([`run_streaming.md`](run_streaming.md) §6.2).

Two of the three already exist in some form — `MissingTag` and `NoReadGroups`
([`open_bam.rs:841`](../../../../src/ng/read/input/open_bam.rs#L841)) — but they collapse into
one error variant, `MissingSampleName { read_group: Option<String> }`, whose `None` arm means "no
read groups at all" and whose `Some` arm means "this read group has no `SM`". Those are two
different problems with two different fixes. **Split them**, and put the fix in the message: a user
who hits either has to go and re-header a file, and the error is the only thing that tells them
which.

### The synthesized library name

When a read group carries no `LB`, its library name is built from **sample, `@RG ID`, and the input
file's name without its extensions**. All three are needed:

- without the file name, two files that each declare `ID:1 SM:X` with no `LB` produce the same
  library name and are silently fused;
- with it, the name is unique as long as no two input files share a name, because `@RG ID` is
  already unique within a file.

**If two files do share a name — same file name in different directories, which is exactly how
per-sample directory layouts look — that is a hard error** naming both full paths. The alternative,
appending a disambiguating integer, was rejected for a specific reason: the suffix would depend on
input order, so the same data listed differently would produce different library names, and those
names reach reports and output. An order-dependent identifier is worse than an error.

Using the full path instead of the file name would also be collision-free, and was rejected because
it embeds machine-specific paths in a label that should be stable across machines.

**What the file name costs.** Two lanes of one un-`LB`'d library, split across two files, now get
two library names where the older `(sample, id)` rule would have given them one. Nothing is lost:
`sample`, `id` and `file` are all still on the record (§3), so a consumer that wants those lanes
merged can group on `(sample, id)` itself.

### The synthesized experiment name

`SRX` is not a tag the SAM specification defines, so it will be absent far more often than `LB`.
When it is absent, **the experiment name is the library name** — not `(library, id)`. Several read
groups legitimately share one declared library; that is a library sequenced across several lanes,
and those lanes are one preparation, so one experiment. Falling back to `(library, id)` would split
them.

Falling back to the coarser value is safe here for the same reason as above: the read group's `id`
survives, so anything wanting the finer split can take it.

### What stays soft

`PL` is carried and **nothing keys on it**. The archive contains `ilumina` (misspelled) and a
project that encoded a second value in it until
[`fix_cram_headers.sh`](../../../../benchmarks/tomato1/scripts/fix_cram_headers.sh) rewrote the
headers. It is a report column.

---

## 7. Stamping the read

Which read group a read belongs to is normally the record's own `RG:Z` auxiliary tag — but only a
file that declares several read groups has a question to answer. **How a file's records are resolved
is decided once, at open, from how many `@RG` records its header declares:**

- **Exactly one declared read group → every record belongs to it, and the tag is not read.** Not "we
  skip the check": there is nothing to check. A read in such a file has exactly one group it could
  belong to whatever its tag says.
- **Several declared → resolve `RG:Z` per record** against a small name-to-identifier map built at
  open. A linear scan over the file's handful of read groups is enough. **A record with no `RG:Z`,
  or one naming a group the header did not declare, is a hard error**, fatal to the run — with
  several candidates there is no way to assign it, and guessing would misattribute a read to the
  wrong library.

**Why the first rule is right and not merely cheap** (owner, 2026-07-27). A file whose header gained
an `@RG` from `samtools reheader` — a header edit that does not rewrite the records — has untagged
records and one declared read group. That is a real and ordinary input. Under a rule that demanded a
tag on every record it would fail, and the only remedy would be rewriting every record in the file
to add a tag carrying no information. Every file in the surveyed archive declares exactly one `@RG`,
so this is also the universal path today, and it costs nothing per read: the same stamp-before-read
mechanism that carries `source_file_index` now
([`filtering.rs:438`](../../../../src/ng/read/filtering.rs#L438)).

**The accepted risk, stated once.** In a single-`@RG` file, a record whose tag names some *other*
group is assigned to the declared one without complaint. Detecting it would mean reading the tag on
every read of every file to catch a case that a hand-edited header is about the only way to produce.
Accepted.

**The trap: this keys on the file, not on the open.** A file declaring three read groups, opened for
a sample that owns only one of them (§9), still resolves per record — the other two samples' reads
have to be recognised in order to be left out. The one-declared-group shortcut applies only when the
*file* declares one, in which case the file is single-sample anyway.

**Where resolution happens.** In the record source's `read_next`, which already stamps the per-file
tag on the reused buffer and is the first point where the filled record is available. It writes a
resolved identifier into the buffer; `decode` copies it into the read. A resolution failure is an
`Err` from `read_next`, which the record-source contract already treats as fatal to the run
([`filtering.rs:343`](../../../../src/ng/read/filtering.rs#L343)).

**Shape it as an enum, not an `Option<ReadGroupId>`.** The obvious signature hands the resolver a
default identifier that is `Some` for a single-group file and `None` otherwise — but then `None`
means "several, use the map", which is a sentinel rather than an absence, and this repo's specs pin
`Option<T>` to absence ([`read_filtering.md`](read_filtering.md) §2.4). Two variants say what the
two modes are and make the fast path visible in the type:

```rust
/// How this file's records are assigned to read groups — fixed at open.
pub enum ReadGroupResolution {
    /// The file declares exactly one; every record is this one.
    Sole(ReadGroupId),
    /// The file declares several; each record's `RG` names which.
    PerRecord(/* name → ReadGroupId */),
}
```

**The per-read error needs a real error type.** `RawAlignedRead::decode` returns `io::Result<MappedRead>`
today — a deviation recorded in
[`../impl_plan/read_filtering.md`](../impl_plan/read_filtering.md), forced by reusing production's
decoder. An `io::Error` here would produce a message with no way to find the offending read. The
ng-owned read type (§8) makes an ng error type possible, carrying path, read name and position.

---

## 8. What changes in the modules that exist

| what | today | after |
|---|---|---|
| several `@RG` in one file naming several samples | error at open, `SampleNames::Several` ([`open_bam.rs:839`](../../../../src/ng/read/input/open_bam.rs#L839)) | not an error; several read groups per file is normal |
| the sample-name agreement check | twice — within a file (`sample_names`) and across files ([`agreed_sample_name`, `mod.rs:550`](../../../../src/ng/read/input/mod.rs#L550)) | once, at `SampleReads::open`, over the read groups that open was handed |
| `AlignmentFile::sample_name() -> &str` ([`:417`](../../../../src/ng/read/input/open_bam.rs#L417)) | one string per file | gone; the sample is a property of the read-group table |
| `SampleReads::open(&[PathBuf], …)` ([`mod.rs:363`](../../../../src/ng/read/input/mod.rs#L363)) | takes paths, derives the sample by agreement | takes the sample's read groups, which know their paths |
| per-file `ReadFilterCounts` ([`mod.rs:413`](../../../../src/ng/read/input/mod.rs#L413)) | keyed on the file, "never summed" | keyed on the **read group** |
| the read type | production's `MappedRead` ([`alignment_input.rs:78`](../../../../src/bam/alignment_input.rs#L78)), reused | ng-owned, with `read_group: ReadGroupId` |

Three of these deserve their reason.

**The two sample checks collapse into one.** Within-file and across-file agreement are the same
question — *do these read groups name one sample?* — and it does not care which file each read group
came from. After the change there is one check, at `SampleReads::open`, over the read groups that
open was given (§4).

**One sample per open stays** (owner, 2026-07-27). Note what the old constraint conflated. It
forbade a *file* from naming several samples, and that goes: several read groups per file, naming
several samples, is a normal input (§9). What stays is that a `SampleReads` serves one sample.
Lifting *that* would mean demultiplexing a single pass over a file into several samples' streams —
deferred in §12, and the read-group identifier is what would make it a local change.

**Counts move to the read group** because the reason they exist stops being true otherwise. The
rationale written at [`mod.rs:406`](../../../../src/ng/read/input/mod.rs#L406) is that a bad run
shows up as one file with an anomalous drop rate, so summing erases the signal. Once a file may
hold several read groups, a per-file tally conflates them and erases the same signal. Keying on the
read group restores it.

**ng gets its own read type** (owner, 2026-07-27). Adding a field to production's `MappedRead` was
the alternative; it is a small change — one non-test constructor,
[`record_buf_to_mapped_read`](../../../../src/bam/alignment_input.rs#L803), every other construction
site being `cfg(test)` — but it makes production carry an ng concept. The ng type copies the
constructor, not the logic: `compute_adaptor_boundary` and `cigar_to_ops` are pure and stay
production's. It also drops `source_file_index`, which the read group makes redundant, and it
unblocks the per-read allocation work already parked as measure-first in
[`../impl_plan/read_filtering.md`](../impl_plan/read_filtering.md) — not to be done here.

**One thing this leaves open.** On the generic path a read becomes a `PreparedRead`
([`src/pileup/walker/mod.rs`](../../../../src/pileup/walker/mod.rs)), production's type, which
carries no source tag at all — so the identifier dies there. The STR path is unaffected: its
generator works on the read type directly (`kept: Vec<MappedRead>`,
[`ssr.rs:312`](../../../../src/ng/locus_generation/ssr.rs#L312)). Owning the read type makes owning
a prepared-read type later a second copy rather than a surprise. Deferred (§12).

---

## 9. Several samples in one file

**Decision (owner, 2026-07-27): evidence gathering stays per sample, and an open serves exactly one
sample (§4).** If a file holds several samples, it is opened once per sample — by separate threads
if that is how the work is scheduled — and each open is given, and reads, only its own sample's read
groups. Multi-sample files are unusual, and a demultiplexing path that every read goes through is a
poor trade against a rare input.

What that requires:

- **Filter by read group at the file's edge.** Each `SampleReads` yields only reads whose read
  group belongs to its sample. When every read group in the file belongs to this sample — the
  universal case today — the predicate is settled at open and the per-read path is untouched. Only
  a genuinely shared file pays a test per read, and it has already resolved `RG:Z` by then anyway.
- **A read belonging to another sample is not a drop.** It must not land in `ReadFilterCounts`,
  whose categories answer "how did this read group behave"; a foreign read says nothing about
  quality. Give it its own tally so the accounting still closes.
- **Nothing may assume a path appears once in a run.** Error text, per-file reporting, any future
  cache: key on the read group or the sample, not the path.

**The cost, accepted** (the "reader pool" below became *one cursor per worker* at
[`alignment_cursor.md`](alignment_cursor.md)'s Milestone F — a cursor opens its own descriptor
and holds it, so nothing is pooled and nothing is lent; the accounting is otherwise unchanged):
per multi-sample file, per sample, one index parse, one reader pool, and for
CRAM one reference repository; records shared by several samples are decoded once per sample. All
of it is invisible on single-sample inputs, and every file in the surveyed archive is
single-sample.

**What needs no change.** The duplicate-read check
([`DuplicateReadAcrossFiles`, `mod.rs:622`](../../../../src/ng/read/input/mod.rs#L622)) is scoped to
one sample's files, so opening one path under two samples never trips it while passing the same
path twice within one sample still does. The reader pool is per `AlignmentFile`, so two opens get
two independent pools with nothing to coordinate.

**Downstream, for context, not decided here.** Cohort consumers need every sample's observations at
one locus. That vec is assembled *above* the generator by a driver that walks the segment stream
once and runs the per-sample generators in step — the shape
[`ng_ssr_cohort_stutter.rs`](../../../../examples/ng_ssr_cohort_stutter.rs) already prototypes. The
per-sample generator contract stays as it is. One note for whoever writes that driver:
`SampleLocusObservations` carries `region`, `reference_bases` and `kind`
([`locus_generation/mod.rs:34`](../../../../src/ng/locus_generation/mod.rs#L34)), all properties of
the locus rather than of a sample, so a bare vec of them duplicates both flanks of an STR per
sample.

---

## 10. Cross-cutting concerns

**Performance.** The per-read cost is unchanged on every file in the surveyed archive: one integer
stamped, no auxiliary-tag lookup, no allocation. The new costs are one extra header read per file
and, on multi-sample files only, a duplicated index parse per sample.

**Memory.** The run-wide table holds one `ReadGroup` per `@RG` record — 2,085 of them for the whole
surveyed archive, since every file there declares exactly one. Not a consideration.

**Errors.** Everything checkable from a header is checked in the pre-pass, before a run spends an
hour. The only per-read error is an undeclared read group, and it is fatal, matching the existing
convention that a decode failure ends the run while a filtered read is tallied.

**Concurrency.** Identifier assignment is a serial pass over an ordered list, so it is deterministic
whatever order the opens finish in. Nothing else here is shared: each `AlignmentFile` owns its
reader pool, and two opens of one path share nothing.

---

## 11. Test obligations

- A file with no `@RG` fails at the pre-pass, and the message names the file.
- An `@RG` with no `SM` fails at the pre-pass, and the message names that `@RG ID`. Distinct error
  from the previous case.
- Two files each declaring `ID:1` get two identifiers, and each read resolves to its own.
- A file with no `LB` produces a library name containing sample, `@RG ID` and file name, marked
  `Synthesized`.
- Two input files with the same name in different directories, both without `LB`, fail with both
  full paths in the message.
- A read group with no experiment tag takes its library's name; several read groups sharing one
  declared library share one experiment name.
- Identifiers are unchanged when the input files are opened concurrently, and change with the input
  order only in the way §4 specifies.
- A file declaring one `@RG` whose records carry no `RG` tag reads normally, every read carrying
  that group's identifier — the `samtools reheader` case (§7).
- A file with several `@RG` resolves each record to the right read group; in that file a record with
  no `RG` tag, and one naming an undeclared group, are both fatal errors naming the read.
- `SampleReads::open` handed read groups that name two samples fails before any read flows, naming
  both samples — the single-sample-per-open invariant (§4), including when the two read groups come
  from the same file.
- One path opened under two samples yields, for each, only that sample's reads, with the foreign
  reads counted outside the drop categories.
- Reads carry the same identifiers through the k-way merge, and the existing merge tests still pass
  unchanged.

---

## 12. Deferred, with a home

- **An ng-owned prepared-read type** (§8), so the identifier survives the generic path. Home: a
  revision of [`read_preparation.md`](read_preparation.md).
- **Per-read allocation reuse** in the ng read type (§8) — already parked as measure-first in
  [`../impl_plan/read_filtering.md`](../impl_plan/read_filtering.md). Unblocked by this work, not
  part of it.
- **Dropping the one-sample constraint** (§8) — a deletion plus a grouping step, once there is an
  input that needs it.
- **Demultiplexing a multi-sample file in a single pass** (§9) — worth revisiting only if such
  files become common. The read-group identifier makes it a local change: route by
  `ReadGroup.sample` after decode.

---

## 13. Resolved decisions & open questions

**Resolved (owner, 2026-07-27):**

- **A file with no `@RG` and an `@RG` with no `SM` are hard errors** — no file-name fallback for the
  sample. An earlier draft synthesized a sample name from the file name; rejected because a missing
  `SM` is a broken header, and inventing a name hides it.
- **A record's read group is only required where it is needed** (§7): a file declaring one `@RG`
  assigns every record to it without reading the tag; a file declaring several treats a missing or
  undeclared tag as fatal. An earlier draft required the tag everywhere, which would have rejected
  the ordinary `samtools reheader` output and left the requirement unenforceable on exactly the
  files where it was demanded.
- **Missing `LB` is lenient**, synthesized from sample, `@RG ID` and file name, with a hard error on
  a residual collision. Beat: a disambiguating integer (order-dependent labels), and the full path
  (machine-specific labels).
- **`SRX` follows the same policy and falls back to the library**, not to `(library, id)`, which
  would split the lanes of one preparation.
- **The identifier space is the run**, not the sample. Beat: per-file identifiers with a base offset
  (a file opened for two samples would get two numberings of one read group), and a store files
  register with at open (mints duplicate identifiers when one file is opened twice).
- **Several read groups per file is normal, and they may name several samples; an open still serves
  exactly one sample** (§4), enforced at `SampleReads::open`. The old check forbade the first; only
  the second survives.
- **Evidence gathering stays per sample**; a multi-sample file is opened once per sample (§9).
- **ng owns its read type** rather than adding a field to production's `MappedRead` (§8).

**Open — confirm before code:**

- **Which tag spells the experiment.** `SRX` is not in the SAM specification, so this is a choice:
  read a tag literally named `SRX` from the read group's other fields, read `DS`, or make it
  configurable. `PU` is not a candidate — it is the platform unit, which is the run or lane, not
  the preparation. *Leaning:* look for `SRX`, fall back to the library, and revisit if a real input
  ever carries something else. *What would settle it:* the survey did not look for experiment tags
  at all; re-running it to count which of these tags actually appear in the archive would decide it
  in an afternoon.
