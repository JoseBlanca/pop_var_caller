# ng — the reference's tandem-repeat catalog, built in the pass that already reads the FASTA

*Design spec, 2026-08-10. **No code yet — this settles the design.** **Independent work**: it stands
on its own, is meant for `main`, and is implemented as its own job. Nothing in it depends on the
parameter pre-pass, though
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.3 is one consumer and is why
the question came up.*

*It stands on the tandem-repeat scanner ([`ssr_repeat_scanner.md`](ssr_repeat_scanner.md)), the
reference-info reader ([`reference_info.md`](reference_info.md)), and step 3's classification policy
([`typed_regions.md`](typed_regions.md)). `src/ssr/` is frozen production: everything said about it
here is a record, not a change.*

*Naming: **STR** in prose, `ssr` in code.*

---

## 1. What this is

**Scan the reference for tandem repeats once per genome, and write what was found beside it.**

**The reason it has to exist is the parameter pre-pass's STR sample.** Fitting the stutter model needs
`cap` loci drawn at random from each (period, repeat count) stratum, and "the `cap` lowest hashes in
the stratum" requires **the stratum enumerated first**
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.2) — a genome-wide list of
STR loci with their strata, which nothing in ng can produce today. Not a speed problem: without this
file the selection cannot be made at all.

Today's scan happens inside every walk that needs it, so a fifty-sample cohort detects the same
repeats fifty times and gets the same answer fifty times, and the answer is unavailable to anything
that has not started a walk — so even *how many STR loci does this genome have, and how are they
distributed over repeat counts?* cannot be asked. That is the same gap seen from the other side.

**It is a catalog, and that is all it is.** It records where the reference repeats. What any other
part of the code does with those records — which it calls an STR, which it sends to a stutter model,
which it ignores — is that part's business, and the file imposes none of it.

**What is in it:** every tandem repeat of period 1 to 6 that reaches the copy floor
`[5, 5, 4, 4, 4, 3]` and has at least **15 bp of sequence beside it on each side**. Nothing else is
filtered: no purity floor, no satellite cap, no bundling, no homopolymer rule.

**Why those two thresholds and no others.** The copy floor sits below the count at which any measured
tomato library stutters (§4.1), so it removes rows nothing would ever route while leaving every
calling floor reachable by filtering. The 15 bp is what makes a tract addressable at all: with less
than that beside it, no read can be anchored to it and no caller can use it.

### 1.1 The contract, in one sentence

**The catalog holds enough that a reader can derive the genome's segmentation, and the STR loci under
any copy floor it chooses, without opening the FASTA** — provided its policy is no more permissive
than the one the file was built at (§4).

That is what makes the two things the consumer wants possible from a file: **segment the genome**, and
**enumerate the STR loci by (period, repeat count) so they can be counted per stratum and sampled**.

### 1.2 Goals

1. **One scan per (reference, scan settings)**, not one per sample or per run — done by a command
   whose only job is that (§2.6), and reused by every run afterwards.
2. **A reader derives, it does not re-scan.** Changing a copy floor, a purity floor or a satellite cap
   is a filter over the file.
3. **Attach the work to a pass that already exists** rather than adding a traversal of the reference
   (§2).
4. **Refuse rather than under-report** when a reader asks for something the file cannot answer (§4).

### 1.3 Non-goals, and what it does not do

- **It does not classify.** No purity floor, satellite cap or bundling is applied when the file is
  written, and the copy floor it does apply (§4.1) is below every floor a caller routes on. Which of
  these tracts is an STR is the reader's decision, and §5 is how it is made.
- **It does not replace step 3.** Region typing keeps scanning as it does today; this file is
  **designed so that step 3 could read it later** (§8), and wiring that up is not part of this work.
- **It does not touch production's catalog.** `src/ssr/catalog/` builds one through the external
  `trf-mod` binary; ng has no trf-mod and does not depend on it
  ([`typed_regions.md`](typed_regions.md) §5.2 keeps it as a comparison oracle, never a dependency).
- **It stores no bases** (§3.3).
- **It is not a partition of the genome.** `type-regions` writes one of those
  ([`typed_regions_cli.md`](typed_regions_cli.md)); a partition assigns every base a kind, and this
  file records only where repeats are. §8 says what their relationship should become.

---

## 2. Where the work attaches

### 2.1 The pass that already reads the whole FASTA

**`src/ng/reference_info.rs` streams the reference from byte zero**, reconstructing each contig's
geometry and computing per-contig and whole-reference MD5s "in one buffer" (`reference_info.rs:17-21`).
It is the only place in ng that reads the entire file, and it already has the pattern this work needs:
`read_reference_verifying_or_creating_fai` writes a sibling `.fai` when one is missing "so the next run
starts faster".

**That the bases are already in memory is the whole reason to attach here.** Reading the FASTA is the
expensive part — 800 MB on tomato, 3 GB on human — and while the pass digests a contig, that contig's
sequence is sitting in a buffer. The scan is arithmetic over bytes already paid for. **The catalog
must never cost a second read of the reference**, and every decision in this section follows from
that.

**`RefSeq` is not that place**, despite being the surface every consumer fetches bases through
(`src/ng/ref_seq.rs`). Its three implementations serve `(contig, range)` requests — a whole resident
contig, a sliding window, or a synthetic map. None of them reads the file end to end, and a windowed
walk deliberately never holds more than its buffer.

*Worth knowing rather than assuming: "all reference access goes through one class" is the intent, not
yet the fact.* `RefSeq` and `reference_info` are two doors onto the same file — the first for bases,
the second for the file's shape — and both go through `crate::fasta` beneath. This work touches the
second.

### 2.2 Decision: an observer seam, not a new dependency

**`reference_info` is a deliberate leaf** — its own module doc says it "imports `crate::fasta`
(read-only …) and noodles, and knows nothing else about ng" (`reference_info.rs:12-15`). Putting the
scan inside it would pull the tandem-repeat scanner and step 3's vocabulary into that leaf and invert
the layering.

**So the pass gains a seam, and the catalog builder sits above it.** While streaming, the pass hands
each contig's canonical bases forward, in coordinate order, to an observer; the catalog builder is one
such observer, in its own module, above `tandem_repeat` and beside `region_typing`. `reference_info`
learns nothing about repeats.

**The seam is already half-justified by what the pass does**: it computes MD5s, which is itself a
value derived from the bases as they stream by. The catalog is a second such value. The seam is
"things computed from one forward pass over the reference's bases", and the digest is its first
member.

*Rejected: a separate pass over the FASTA.* It is simpler and costs a second full read — 800 MB on
tomato, 3 GB on human — of bytes the first pass already has in a buffer. Worth reconsidering only if
§2.4's concurrency turns out to be awkward, in which case a standalone builder is a small
rearrangement rather than a redesign.

*Rejected: building it inside `reference_info`.* It reads as the shortest path and costs the leaf
property that makes that module safe to depend on from anywhere.

### 2.3 One contig at a time, whole — which is what removes the chunk-boundary trap

**A tandem repeat does not respect a buffer boundary**, and a scanner fed buffer-sized chunks needs a
margin carried across each one, a rule for which core a straddling detection belongs to, and a cap on
the repeat length it can promise to catch whole — step 3 ties that margin and the satellite cap into
**one field** for exactly this reason (`region_typing/mod.rs`, `TypedRegionConfig::max_str_len`).

**Decision: the observer accumulates a contig's bases as they stream past and scans it whole**, with
`find_tandem_repeats` over one slice. Then there is no boundary to straddle, no margin to carry, and
**no length at which the geometry stops promising a whole tract** — a satellite of any size comes out
as one row.

**What it costs is the contig resident while it is scanned**: 90 MB for tomato's largest chromosome,
250 MB for human chromosome 1, and one such buffer per worker (§2.4). That is this design's memory,
and §10.4 measures it.

**What it buys, beyond simplicity, is two fewer things a reader can be refused on.** With every tract
recorded whole, the satellite cap and the bundle radius stop being properties of the file and become
pure filters over stored spans — a reader may choose any value of either (§4.2). Only the copy floor,
the period range and the minimum flank still bound anything.

*If the memory ever bites*, `scan_windowed` and `WindowCursor` are the drop-in (`tandem_repeat.rs`),
and the margin, the attribution rule and the cap on whole capture all come back with them. It is a
contained swap, not a redesign — but it is a worse file, so it is the fallback and not the plan.

### 2.4 Concurrency: as little as it takes

**This file is built once per reference, by a command whose only job is that** (§2.6), so scan
throughput is not on any path that runs twice. The design spends nothing on it:

- **sequential by default** — each contig is scanned as the stream finishes it, and the pass waits;
- **`--threads N` scans up to N contigs at once**, each on its own buffer, which is the whole of the
  parallelism. Rows are written in the reference's contig order regardless, so the file is
  byte-identical at every thread count (§10.5).

**Nothing else is pipelined**, and the digest keeps its own shape: it is computed from the same bytes
as they stream, in order, exactly as it is today. The catalog is a second value taken from that one
pass (§2.2) — coupling the two is the point, and it is the only coupling there is.

**One consequence to state, because it changes an existing behaviour.** The pass can run in the
background today (`read_fai_verify_in_background`), with the walk proceeding meanwhile. A consumer that
needs the catalog must join it first — not for the detections, which are on disk (§2.5), but for the
digests that say the file describes this reference. Validation is a *prerequisite*, not a background
nicety, for anything that selects loci.

### 2.5 One run builds it; every other run reads it

**The digest pass is not optional and never was** — a run needs the contig geometry and the digests
whether or not anyone wants repeats. The scan is what rides along with it, and **it rides in exactly
one kind of run: the dedicated build (§2.6).** Three states, and only one of them pays for detection:

- **the build command** → the observer scans as the bases stream by, and the file is written beside
  the reference. This is the only case that scans;
- **any other run that wants a catalog** → **no detection is done at all.** It reads the file, and the
  digests the pass computes anyway are what validate the header (§4.3) — so checking the catalog costs
  nothing beyond a comparison. If there is no catalog, or it does not match, the run stops and says
  which command would produce one;
- **a run that wants no catalog** → the pass is what it is today, and the observer is never attached.

**Building only in its own run is what removes the ordering trap.** Whether a catalog on disk
describes *this* reference is not known until the digests are complete — by which time the bases have
streamed past and nothing can ride along any more. A run that could build on demand would therefore
have to choose, before it knows the answer, between scanning defensively every time and being unable
to recover. **Neither choice is needed once building is a separate run**: a calling run only ever
validates and proceeds, or validates and stops.

Never silently rebuild under a name that already means something else.

### 2.6 The command

**`pop_var_caller_exp repeat-catalog`**, a second subcommand of the experimental binary that
[`typed_regions_cli.md`](typed_regions_cli.md) §2 introduces, beside `type-regions`:

```
pop_var_caller_exp repeat-catalog --reference ref.fa [--output ref.fa.repeats.parquet]
                                  [--threads N] [scan knobs]
```

- **`--reference` is the only required argument.** `--output` defaults to a sibling of the FASTA, the
  same shape `write_fai` already uses, so that a later run finds the catalog without being told where
  it is;
- **the scan knobs are the §4.1 axes** — period range, the per-period copy floor table, the minimum
  flank and the two scoring weights — with the defaults of §4.1 rather than step 3's
  calling defaults, since this file is built to outlive any one policy. Every value used is written
  into the header (§3.4), so a later reader checks against what was actually run, not against what the
  defaults say today;
- **`--threads` sets how many contigs are scanned at once** (§2.4). It is a speed knob only: the file
  is byte-identical at every thread count, and §10.5 pins that;
- **an existing output file is not overwritten without `--force`.** A run that is reading the catalog
  while another rebuilds it in place is a failure with no error message, and the dedicated command is
  the one place where "I meant to replace it" can be said out loud;
- **the output is Parquet, one row group per contig** (§3.5), with the header of §3.4 in the footer's
  metadata. `duckdb` opens it directly, which is how it is inspected by hand.

**What it prints when it finishes** is the tally of §5.4: rows written, and the count per period. That
is the measurement §9.1 asks for, so the command produces it as a by-product rather than needing a
separate harness.

---

## 3. What the file holds

### 3.1 One row per repeat found

| field | what it is | why it is here |
|---|---|---|
| contig | the contig id, resolved through the header's table | rows are coordinate-ordered within it |
| detected start, end | the span the scanner reported, **1-based inclusive** | ng's convention (`ng_step_interfaces.md` §1). *Trap: `RepeatInterval` is 0-based half-open in memory (`tandem_repeat.rs:208-219`) and `SsrSegment` is 1-based inclusive — the conversion happens at this edge and is exactly where an off-by-one would live* |
| trimmed start, end | the same tract cut back to whole motif copies at both ends, **or absent** when no clean cut exists | **this is the locus**, and the reason it is stored is §3.2 |
| period | motif length in bp | half of the stratum key, and the axis a reader filters on |
| score | the Ruzzo–Tompa segment total | step 3's `min_score` gate reads it. **Not a TRF score** and not on the same scale — `segment_criteria.rs:509-512` says so, and a threshold carried across from TRF would not mean here what it meant there |
| motif | the repeat unit, verbatim and phase-faithful | **needs the bases** — §3.2 |
| purity | fraction of the **trimmed** tract matching a perfect motif tiling, in `[0, 1]`; absent when the trim is | **needs the bases** — §3.2 |

**Both spans are stored because classification uses both, at different steps.** Bundling and the
pre-screen run on the **detected** span (`segment_criteria.rs:973-977`, `prefilter`); the copy floor,
the purity floor, the flank test and the emitted locus all use the **trimmed** one
(`finish_locus`, `:1030-1063`). A reader holding only one of them would have to read the FASTA to get
the other, which is the cost this file exists to remove.

**Repeat count is not stored**: it is the **trimmed** span over the period, which is the count
`classify` measures and the count a stratum is keyed by. Storing it as well would invite the two to
disagree.

**And the copy floor the builder applies is measured on the detected span**, exactly as `prefilter`
does, not on the trimmed one. That is what keeps bundling identical: a row the builder drops is a row
the reader's own pre-screen would have dropped before bundling, so no tract loses a neighbour it
would have bundled with. The file therefore holds some rows whose *trimmed* count is below the floor,
and that is correct rather than sloppy.

**Rows overlap, and resolving them is the reader's job.** The scanner already drops a period-multiple
re-detection of the same tract (`tandem_repeat.rs:474-476`), so what remains is genuine: two different
tracts whose spans intersect, or a tract detected at two primitive periods. The file keeps both,
because which one survives depends on floors and a purity threshold the file does not apply — and the
reader resolves them with **the same `prefilter` and `classify` that a live scan uses** (§7), which is
the only way §10.1's differential can come out identical rather than merely close. The low-level
scanner interface is documented for exactly this consumer: one "for consumers that resolve overlaps
themselves (the STR catalog)".

**Row order is part of the format**, because §6 asks for a byte-identical file: contigs in the
reference's own order, then by start, then by period, then by end. The last two are tie-breaks that
only overlapping rows reach, and without them a worker-count change could reorder them.

### 3.2 Why motif and purity must be computed here

**This is the finding that decides whether the file is worth anything.** The scanner's output type
carries only four fields — `start`, `end`, `period`, `score` (`tandem_repeat.rs:208-219`). The motif is
the upper-cased, period-length prefix of the tract, sliced from the bases inside `classify`; the purity
is recomputed over the tract, and step 3's floor is "applied after recomputation"
(`segment_criteria.rs:499-502`).

**Both need the sequence, and so does the trim.** Cutting a tract back to whole motif copies
(`minimal_trim`) compares bases against the motif, so a file of bare intervals would send every reader
back to the FASTA to answer *"where does this locus actually start?"* and *"is it pure enough?"* —
the whole cost this work exists to remove. All three are computed during the pass, when the bases are
already in the buffer, so they cost nothing extra to obtain.

**None of the three depends on the criteria**, which is what makes storing them sound: the trim reads
the tract and its motif, the motif is the tract's first `period` bases, and the purity is a count over
the trimmed tract. A different copy floor or purity floor changes which rows survive, never what these
three say — §9.4 is the check that this holds, and §10.1 is what fails loudly if it does not.

### 3.3 What is deliberately not stored: the bases

**Production's catalog embeds the local reference** (`ref_seq` + `ref_seq_start`) because it is "the
*only* reference-bearing input the downstream stages need" (`src/ssr/catalog/mod.rs:5-8`). **ng's
stages have `RefSeq`**, and step 3 already made the matching call: `RegionKind::SsrSegment` is
"**Coordinates, no bases**: the bases are in the reference the caller already has open"
(`region_typing/mod.rs:174-176`).

Embedding would multiply the file by the tract lengths and create a second copy of the truth that can
drift from the FASTA. It is not stored.

### 3.4 The header, which in the file is the footer's metadata (§3.5)

- **the contig table: one row per contig, with its name, its length, and its MD5**, plus the
  whole-reference MD5 over all of them. The pass computes every one of these anyway
  (`reference_info.rs:17-21`), so storing them costs the header a few dozen bytes per contig and
  nothing in compute. Two of the three do work:
  - **the per-contig MD5s are how "is this the right reference?" is answered**, and they answer it
    better than the whole-reference digest alone: that one says only *something changed*, while the
    contig table names **which** contig, which is the difference between a usable error message and a
    puzzle. They also catch a reference that is the same sequences in a different order, which a
    reader's coordinates depend on and a length table would not notice;
  - **the lengths are load-bearing rather than informational**: classification's flank test needs the
    **contig's** length, and `segment_criteria.rs:865-886` documents at length what goes wrong when a
    window's length is passed instead — a locus with a silently truncated flank, or none;
- **the scan settings it was built at** — the period range, the per-period copy floor table, the
  minimum flank and the two scoring weights — the four things the builder actually applied. §4 is what
  these are for. The satellite cap, the purity floor, the score floor and the bundle radius are **not**
  here: the builder applies none of them, and a reader chooses them freely (§4.2);
- **the tool version**, so that a change in the detector invalidates the file even when every setting
  matches.

**The header's digest is what a consumer checks against**, and it must be strong enough to stand in
for "this file describes this reference under these settings" —
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §5.1 leans on exactly that.

### 3.5 Decided: Parquet, one row group per contig

**The file is Parquet**, written and read through the `parquet` crate (the Rust Arrow project), with
**one row group per contig** and the header of §3.4 in the file's key-value metadata.

**Every column here is one Parquet encodes to almost nothing.** `contig` and `motif` are dictionaries
— a genome has tens of contigs, and there are only 5,460 possible primitive motifs of period 1 to 6;
`start` and `end` are ascending within a contig, so they delta-encode; `period` is one byte, `purity`
a float, `score` a small integer. Text spends 8 to 10 bytes writing "12345678" that the column store
writes as a delta of a few bits.

**Four things it gives that the alternative needed extra machinery for:**

- **column pruning**, and the two heaviest methods use it: `count_loci_per_stratum` reads `period`,
  `start`, `end` and nothing else; `sample_loci_per_stratum` adds `contig`. Three or four columns out
  of seven, with `motif` — the widest — untouched by either;
- **the per-contig seek**, from row-group boundaries and their statistics, in place of the byte-offset
  table a TSV needed in its header;
- **truncation caught for free.** A Parquet file's metadata lives in the footer, so a builder that
  died mid-write produces a file that will not open — where a truncated TSV reads as a short but
  well-formed one, which §6 calls out as the failure with no symptom;
- **the header as typed metadata** rather than `##` lines to re-parse: the contig table, the criteria,
  the weights and the tool version go in the footer's key-value map, and reading them is a seek to the
  end of the file.

**What it costs, and why that is acceptable.** It adds the Arrow dependency tree and a format nothing
else in the repo reads, and `bgzip -dc | awk` stops working. Parquet is standard enough now that the
dependency is unremarkable, and the analysis it replaces is better served anyway: `duckdb` reads a
Parquet file directly, so "count loci per period at floor 6" is one SQL line over 20 million rows
rather than an `awk` pipeline.

*Rejected: a bgzip-framed TSV*, the shape production's catalog uses (`src/ssr/catalog/io.rs:1-13`,
`234-285`). It was the front-runner while the dependency looked expensive, and everything it needed —
a byte-offset table per contig, atomic writes to make truncation detectable, text parsing on every
read of every column — Parquet either provides or makes unnecessary. Its one real advantage is that
the repo already knows how to write it.

*Rejected: SQLite (`rusqlite`).* Its strengths here are unused — no updates, no transactions, one
writer, no ad-hoc joins — while B-tree pages plus an index on `(contig, start)` typically land 2 to 4
times the size of the compressed rows, on what is already the largest derived file ng writes.

**Two things the implementation must pin, because Parquet makes them easy to get wrong.**
§6 asks for a byte-identical file, so the **compression codec and level, the row-group boundaries (one
per contig) and the writer-version string** are all fixed by us rather than left to defaults that can
move with a crate upgrade. And **Arrow types stay inside the io module** — the methods of §5.3 hand
back ng's own types, so a format change stays one module deep, as it just did.

---

## 4. The permissiveness rule

**The catalog is built at settings at least as permissive as any reader will ask for, and it records
them. A reader whose policy is more permissive on any bounded axis is refused, not served.** This is
the rule that keeps §1.1's contract honest, and without it the file is a trap: it would answer, and
the answer would be short.

### 4.1 The axes that bound a reader, each for a reason in the code

- **The copy floor, and it is `[5, 5, 4, 4, 4, 3]` for periods 1 to 6.** Step 3 hoists the weakest
  copy floor its policy applies **into the detector**, so that "the intervals that cannot survive it
  are never materialised" (`region_typing/mod.rs`, `TypedRegionConfig::effective_scan`). That is
  output-neutral *for that path*, and it is what a catalog must be careful with: intervals never
  materialised cannot be filtered back in. So the builder hoists the table's **minimum** (3) into the
  scanner and applies the per-period floor as it writes rows — nothing a reader could ask for is lost
  in the detector.

  **Why these numbers.** The calling floors — the copy counts at which a tract starts to stutter,
  measured over the tomato archive — are `[8, 6, 6, 6, 5, 4]` (`segment_criteria.rs`,
  `MinCopies::default`). The catalog sits below every one of them: three repeats below at period 1,
  one below at period 6. A caller can therefore move its routing floor anywhere inside that gap by
  filtering, which is the question the file exists to keep open. A reader asking for less than the
  table at any period is refused (§4.3).
- **The period range.** The scanner "honours its range now", emitting only primitive periods within
  it. A reader asking for a wider range gets silence rather than an error unless the header is
  checked. **Build at the widest** — and that means **including period 1**, since whether homopolymers
  belong on the STR path is a live question: "Widening the floor to 1 puts homopolymers on the STR
  path; narrowing the ceiling takes penta/hexamers off it. Both are questions ng exists to answer"
  (`segment_criteria.rs:483-486`). Homopolymers are in the file from 5 copies up — far enough down to
  ask the routing question at every floor a caller would consider, and short of the 1- to 4-copy runs
  that would dominate the file and settle nothing.
- **The minimum flank, 15 bp on each side.** A tract closer than that to a contig end is not in the
  file, so a reader content with less would come up short at exactly those tracts. **Build at the
  smallest anyone will ask for**, and 15 bp is exactly the flank ng's STR locus generator fetches by
  default — `SsrGeneratorConfig::flank_bp` is `DEFAULT_BUNDLE_THRESHOLD`, 15 (`locus_generation/ssr.rs:132`,
  `segment_criteria.rs:299`) — so a caller at the default is served exactly and one asking for more is
  served by arithmetic on stored coordinates. It bites only within 15 bp of a
  contig end, which is a handful of tracts per contig — but §10.1's differential has to know about it
  (§5.1), or it fails there and nowhere else.

### 4.2 The axes that bound nothing, and the two that must match exactly

**The purity floor, the score floor, the satellite cap and the bundle radius are pure filters over
stored fields**, so any value of any of them is serviceable from any catalog. The last two are here
rather than in §4.1 because §2.3's whole-contig scan records every tract whole: there is no length at
which the file stops being complete, so *where a locus becomes a satellite* and *how near two tracts
must be to bundle* are decisions a reader makes over spans it can already see. They need no bound, and
the header records them only for provenance.

**The two scoring weights are neither bounded nor free.** A different match reward or mismatch
penalty is a different set of tracts, not a subset of the
stored one — the weights decide where each Ruzzo–Tompa segment starts and stops, so a tract can grow,
shrink or vanish, and no filter over the file reproduces that. They are therefore checked like the
tool version: **equal or refuse**, never "at least as permissive". A reader that wants other weights
wants another catalog, and `repeat-catalog` builds it in one run (§2.6).

The same holds for the stored `score` itself: it is comparable across readers only because every row
came out of one weighting, which is why the header carries the weights beside the score floor.

### 4.3 Refusing

A reader states its policy; the builder's header states what was scanned. **More permissive on any
axis of §4.1 → refuse, naming the axis and both values.** Not a warning, and not a silent rebuild
under the same filename: a short answer here is a wrong genome-wide count, and nothing downstream
would notice.

**The FASTA is the source of truth for every digest.** The catalog's header does not state what the
reference is; it states what the builder claims it was, and the claim is checked against MD5s
recomputed from the FASTA in this run. The same holds for the read files' `@SQ M5`. Nothing is trusted
because it was written down — which is why the pass computes the digests unconditionally (§2.5), and
why there is no mode in which a stored digest stands in for the reference.

**The same refusal, and the same specificity, for the wrong reference.** The pass's per-contig MD5s
are compared against the header's contig table (§3.4); a mismatch names the contig and both digests,
and a contig present in one and not the other is named as such. "The catalog does not match the
reference" sends the reader looking through a whole genome; "chromosome 4 differs" does not.

**A mismatch is reported, never repaired.** The run stops and says which contig and which two digests;
it does not rebuild the catalog, fall back to a live scan, or carry on with the rows it has. Somebody
has pointed the run at a reference that is not the one the catalog was built from, and that is a fact
about the inputs, not a condition for this code to work around — silently absorbing it turns a wrong
reference into wrong genotypes with nothing in the log.

**The same digest checks the read files.** The per-contig MD5 stored here *is* the `@SQ M5` a
BAM/CRAM header carries (`reference_info.rs:42-44`), so reference, catalog and read files are checked
against one value rather than three vocabularies, and a read file digesting to something else is the
same class of error and gets the same treatment: named, reported, fatal.

### 4.4 What permissiveness costs

**More detections stored than any one caller keeps, and how many is unmeasured.** The floors of §4.1
admit tracts three repeats below the calling floor at period 1 and one below it at period 6, so the
file holds several times the loci a single policy routes to the STR path; the overlapping-interval
structure also means the row count is not the locus count. **Soft: every number about this file's size
is arithmetic until §10.4 measures it.** If it comes out badly, period 1 is the knob — 5-copy
homopolymer runs will be the densest rows in the file — and raising it costs the ability to ask
whether homopolymers route at 5, 6 or 7 copies.

---

## 5. What a reader derives

### 5.1 The segmentation, with no FASTA open

Given the file and a policy, in order: drop detections outside the period range or below the score
floor; keep those meeting the copy floor for their period and the purity floor; bundle survivors
within `bundle_threshold` of another repeat; call a tract longer than the satellite cap a satellite;
clamp flanks at the contig's own end; everything else is generic.

**Nothing in that list reads a base.** The two steps that look as though they might do not:
`bundle_threshold` "is a distance, and nothing here reads the bases it measures"
(`segment_criteria.rs:530-533`), and the flank clamp needs the contig's length, which is in the
header (§3.4).

**The acceptance test is a differential** (§10.1): the segmentation derived from the catalog must equal
`partition_resident` run over the same reference at the same policy. That function is already
positioned as an independent oracle for the windowed walk (`region_typing/mod.rs`), so it serves here
too.

**With one stated exception, and it is the file's 15 bp flank.** Step 3 drops a locus only when its
flank clamps to *zero* — a tract abutting base 1 or the contig's last base (`segment_criteria.rs`,
`RejectionReason::FlankClamped`) — so a live scan keeps tracts 1 to 14 bases from a contig end that
the catalog never recorded. The differential is therefore run **at a reader flank of 15 bp or more**,
which every calling policy already satisfies, and a reader asking for less is refused by §4.1 rather
than served short.

### 5.2 The type every question is asked in: `StrRepeatCriteria`

**One value says which tracts count as STRs** — how small they may be, how large, how much room they
need beside them — and it is the same value on both sides of the file:

```
StrRepeatCriteria {
    periods:          1..=6,               // motif lengths considered at all
    min_copies:       [5, 5, 4, 4, 4, 3],  // per period 1..6; a wider period falls to the last
    min_flank_bp:     15,                  // sequence required on each side, to the contig's end
    max_str_len_bp:   500,                 // longer than this is a satellite, not a locus
    min_purity:       0.8,                 // fraction matching a perfect tiling
    min_score:        0,                   // the Ruzzo–Tompa segment total
    bundle_threshold: 15,                  // two tracts nearer than this are one bundled locus
}
```

**It is renamed from `StrRepeatFloors` because it now carries a ceiling too** — the satellite cap — and
a type called *floors* holding a maximum is a name that lies. Everything a reader states about which
tracts it wants is in here, so `genome_segments`, `str_loci`, `count_loci_per_stratum` and
`sample_loci_per_stratum` each take exactly one policy argument.

**Three of the seven are what the builder honours** — `periods`, `min_copies`, `min_flank_bp` — and
those are what the header stores and what §4.3 compares. The other four are read-time filters over
stored fields (§4.2): any value of them is serviceable from any catalog, and a reader that changes its
purity floor or its satellite cap needs no new file.

The refusal is then one comparison, `built.serves(asked)`, rather than a check spelled out again at
each call site: **a reader is served when its floors are no lower than the built ones on every period
and on the flank, and its period range no wider**; otherwise it is refused, naming the period, the
flank or the range, and both values.

### 5.3 The methods

Names first, and each says what it returns rather than how it works:

| method | what it gives back |
|---|---|
| `open_checking_against_reference(path, &ReferenceInfo)` | the catalog — or **one of two distinct errors: there is no catalog at that path**, which names the command that builds one (§2.6), or the catalog is not this reference's, which names the contig whose MD5 differs (§4.3). A caller can act on the first and cannot on the second, so they are never the same error. Reads the header only — no repeat row is touched |
| `create_beside_reference(reference, criteria)` | writes the file next to the FASTA during the digest pass (§2.6) and returns where it put it |
| `genome_segments(criteria, region)` | the genome's segments in coordinate order — STR segments and the generic spans between them (§5.1). **Streams**, and `region` narrows it to one contig or one interval, which is what the parallel sharder needs to hand whole segments to workers |
| `str_loci(criteria, region)` | the surviving STR loci alone, without the generic spans. What a sampler picks from, and what a caller enumerating repeats wants |
| `count_loci_per_stratum(criteria, region)` | how many loci in each (period, repeat count) stratum. What the build command prints (§2.6) and what the sampling consumer reweights by (§5.4) |
| `sample_loci_per_stratum(criteria, region, cap, seed)` | up to `cap` loci from each stratum, the ones whose `hash(contig, start, seed)` is lowest. **The seed is the caller's**, because it is a run input in the consumer that needs this ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.4) and two runs at the same seed must keep the identical loci — it is not a property of the file |
| `repeats_in_region(region)` | every tandem repeat the file holds there, exactly as stored and with **no floor applied** — the rows themselves, not a filtered view of them. `region` may be the whole reference. Used by surveys of what the scan emitted (§9.1) and by the comparison against production's catalog (§8) |
| `build_settings()` and `contigs()` | the criteria the file was built under and the weights it was scored with, and the contig table with lengths and MD5s — answerable without reading the body |

**The four methods that take `criteria` — `genome_segments`, `str_loci`, `count_loci_per_stratum`,
`sample_loci_per_stratum` — can all refuse**, and always for the same reason: the caller is asking
about tracts the file does not contain — a copy floor below the built table, a flank below the built
one, a period outside the built range, or scoring weights that differ at all (§4.2). Each names the
axis and both values.

**`repeats_in_region`, `build_settings` and `contigs` take no criteria and never refuse** — they
report what is there. That is what makes the file usable for a question nobody has posed yet, without a re-scan.

**What the sampling costs, since it is the method that sounds expensive and is not.**
`sample_loci_per_stratum` is one forward pass with a bounded heap: each surviving locus gets a
64-bit hash of `(contig, start, seed)` and is offered to its stratum's heap, which holds at most
`cap` values ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.2). **Nothing
is sorted, nothing is randomly accessed, and no locus is visited twice.** Working memory is
strata × `cap` — a few hundred strata at a cap in the thousands is a handful of megabytes — and the
result does not depend on the order rows arrive in, so it survives sharding: merging two shards is
taking the lowest `cap` of the union.

**So the sampler costs the scan, and nothing above it.** Which is why `count_loci_per_stratum` and
`sample_loci_per_stratum` must be answerable **in the same pass** when a consumer wants both, as the
pre-pass does: counting is a tally beside the heaps, not a second traversal.

**Why the flank is a minimum and not an exact match.** The file holds tracts with at least 15 bp on
each side, so a reader whose rule is *at least 15* is served exactly, one asking for *at least 30* is
served exactly too — its extra requirement is arithmetic on coordinates the file already has — and
one asking for *at least 5* is refused. Step 3's own rule today keeps any tract with a non-zero flank
(`RejectionReason::FlankClamped`), which is *less than 15*, so it is refused by this same test until
it adopts a flank floor.

### 5.4 The strata, which is what the sampling consumer needs

For each surviving locus the stratum is **(period, repeat count)**, the count being `span / period`.
Tallying gives, for every stratum, how many loci the genome — or a region subset — holds. That is the
whole input to
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3: enumerate, count per stratum,
keep the `cap` lowest hashes in each — which is `count_loci_per_stratum` followed by
`sample_loci_per_stratum` (§5.3), so that consumer holds no selection logic of its own.

**And it answers the question at a floor that has not been chosen yet**, which is the point: the same
file yields the strata at today's calling floors `[8, 6, 6, 6, 5, 4]`, or one repeat below them, or
two, without a second scan.

---

## 6. Cross-cutting concerns

**Memory.** The observer holds **one contig's bases plus that contig's rows**, and under `--threads N`
one such pair per worker (§2.3, §2.4). The largest tomato chromosome is ~90 Mb and human chromosome 1
is ~250 Mb, so sequential building peaks a little above the largest contig and N workers multiply it.
That is the number to watch, and §10.4 measures it. Nothing accumulates across contigs.

**Compute.** The scan is the cost: a lag-`p` self-comparison plus a Ruzzo–Tompa pass, for each period
in the range, over every base. The digest pass it attaches to is I/O-bound, so **the combined pass runs
at the scan's speed, not the read's** — that is the honest statement, and §2.4's per-contig workers are
what keep it from being a serial multiple.

**Errors.** Three failures, and they are different: a **stale** catalog (a contig's MD5 or the
reference's disagrees with the pass's) is refused, naming the contig; an **under-permissive** catalog
(§4.3) is refused, naming the axis and the two floors; a
**truncated** catalog — the builder died mid-write — must not be readable as a short but valid one,
which Parquet gives directly: no footer, no file (§3.5). Write atomically anyway, as `write_fai`
already does (`reference_info.rs:816`), so a half-written file never sits under the expected name.
**None of the three is a
warning, a fallback or a rebuild**: each stops the run with the contig or the axis named. An error
here is about the inputs the user chose, and the failure it prevents — a genome-wide count taken over
the wrong genome — leaves no other trace.

**Determinism.** Same reference and same settings must give a byte-identical file, independent of the
thread count — the same property step 3 states of `window_bp` and pins with a window-invariance test.
Two things carry it: the row order of §3.1 (contig, start, period, end), which is what makes
"byte-identical" reachable at all once contigs finish out of order, and the fixed Parquet write
settings of §3.5 — codec, level, row-group boundaries and writer version — since any of those left to
a default would change the bytes without changing the content, and would move under a crate upgrade.

---

## 7. Reuse over rewrite

| what | existing code | how it is reused |
|---|---|---|
| detection | `src/ng/tandem_repeat.rs` — `find_tandem_repeats` | used as-is, one call per contig (§2.3). It takes a slice and a period range and is already documented as the interface "for consumers that resolve overlaps themselves (the STR catalog)". `scan_windowed` and `WindowCursor` are the fallback of §2.3, not the plan |
| the whole-FASTA pass | `src/ng/reference_info.rs` | gains the observer seam of §2.2; nothing else changes |
| the sibling-file pattern | `reference_info.rs` — `sibling_fai_path`, `write_fai`, `read_reference_verifying_or_creating_fai` | the same shape: a derived file beside the reference, written atomically, created on request, validated before use |
| motif and purity | `src/ng/region_typing/segment_criteria.rs` — the slicing and recomputation inside `classify` | the two computations are lifted so both `classify` and the builder use one implementation; **neither is re-derived here** |
| classification, at read time | the same module — `classify`, `prefilter`, `SsrSegmentCriteria` | used as-is by §5's reader, fed intervals from the file instead of from a live scan |
| the parity oracle | `region_typing::partition_resident` | §10.1's differential |

**The parity oracle is `partition_resident`.** This is a port in one direction only — the same
classification over the same detections, reached by a file instead of a scan — so identical output is
the bar, not a rough agreement.

---

## 8. Deferred, with a recommended home

- **Step 3 reading the catalog instead of scanning.** Designed for and not wired: §5.1's derivation is
  exactly what the walk does, so the remaining work is the walk taking its detections from a reader
  rather than from `scan_window`. Two things have to be true first, and both are checkable: the walk's
  own machinery is finished (`partition_windowed` and `TypedRegionIterator` are the parts still to
  come, `region_typing/mod.rs:9-11`), and §10.1's differential is green. **Home:** a follow-up to
  [`typed_regions.md`](typed_regions.md).
- **The relationship to `type-regions`' partition file.** That command writes a partition; this writes
  detections. The partition could become a *view* over this file rather than a second scan's output,
  which would also give it a reader it currently lacks — its spec lists "reading the file back into
  ng" among its non-goals precisely because nothing consumed it. **Home:**
  [`typed_regions_cli.md`](typed_regions_cli.md) §8, where that deferral already lives.
- **A region subset.** This file describes a whole reference; a run under a `--regions` BED wants
  counts within the BED. That is an intersection at read time and needs no second file, but nothing
  here specifies it. **Home:** the consumer that first needs it, which is
  [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md).
- **Comparing against production's catalog.** `src/ssr/catalog/` over the same reference is an
  independent yardstick for the detections, in the same way
  [`typed_regions.md`](typed_regions.md) §8 uses it for classification. Not a dependency, ever.
  **Home:** a report, not a spec.

---

## 9. Open questions

1. **How many detections does the scan emit, and how big is the file?** — OPEN, and the first thing to
   measure, because §4.4's whole trade rests on it. *Leaning:* affordable — a Parquet row of these
   seven columns encodes to a handful of bytes (§3.5), so even tens of millions of rows is a file
   measured in tens to low hundreds of megabytes — but 5-copy homopolymer runs are dense and nobody
   has counted them. **Settled by:** a scan of the tomato
   reference at periods 1–6 and the floors of §4.1, tabulated by period.
2. **What copy floor to build at?** — **SETTLED, 2026-08-10: `[5, 5, 4, 4, 4, 3]` for periods 1 to 6**,
   on the owner's call. It is below every calling floor of `[8, 6, 6, 6, 5, 4]`, by three repeats at
   period 1 and one at period 6, so the routing floor stays movable by filtering; and it is high
   enough that the 1- to 4-copy tracts nothing would ever route never reach the file. If question 1
   comes back badly, period 1 is the retreat, and raising it must be recorded as foreclosing the
   low-copy homopolymer routing question rather than as a tuning choice.
3. **Does the combined pass slow the reference read?** — **NOT A QUESTION.** It does, and it does not
   matter: the scan runs once per reference, in a command whose only job is that (§2.6), and the file
   exists because **the per-stratum random selection of STR loci needs the genome enumerated**
   ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.2) — that selection cannot
   be made at all without it, at any speed. Wall clock is still reported by §10.4, as a fact about the
   command, not as a decision input.
4. **Is classification's coordinate derivation independent of the criteria?** — OPEN, and it is the
   assumption under §1.1's contract. `classify` derives a locus's coordinates "by arithmetic on
   detector output" and recomputes purity over the tract; if any of that arithmetic reads a value from
   `SsrSegmentCriteria`, then a purity stored at write time is not the purity of the segment a
   different policy would produce. *Leaning:* independent, from a read of the gates — they filter
   rather than reshape. **Settled by:** §10.1's differential run at several policies, which fails
   loudly if it is not.
5. **What file format?** — **SETTLED, Parquet (§3.5)**, on the owner's call that the Arrow dependency
   is unremarkable today. It was the better file all along — column pruning for the two heaviest
   methods, dictionary and delta encoding on every column, row groups in place of an offset table,
   truncation caught by the missing footer — and only the dependency held it back. §10.4 still reports
   the size and the scan time, now as facts rather than as a trigger.
6. **What satellite cap to build at?** — **DISSOLVED.** It was a build setting only while the scan was
   chunked and the margin capped what could be caught whole. With §2.3's whole-contig scan every tract
   is recorded whole, so the cap is a reader's filter over stored spans (§4.2), carried in
   `StrRepeatCriteria` with a default of 500 bp against step 3's calling default of 100. Changing it
   never needs a new file.

---

## 10. How we know it works

1. **The derived segmentation equals the scanned one.** Run `partition_resident` over a reference at a
   policy; derive the segmentation from the catalog at the same policy; the two must be identical —
   regions, kinds, coordinates, motifs, purities. **Run it at several policies**, including ones that
   differ from the catalog's build settings on every bounded axis, since that is what §9.4 turns on and
   a single-policy check would pass either way.
   **Include overlapping rows deliberately**: a tract detected at two primitive periods, and a pair of
   intersecting tracts. The derived segmentation must make the same choice as `partition_resident`,
   since both run the same `prefilter` and `classify` (§3.1) — and a differential over ordinary
   sequence can miss this case entirely.
2. **A more permissive reader is refused.** Ask for a copy floor below the file's, a flank below it,
   and a period outside its range. Each must produce an error naming the axis and both values — **not
   a short answer**, which is the failure this rule exists to prevent and the only one that is silent.
   **A reader with different scoring weights is refused even though every floor matches** (§4.2), the
   refusal that looks wrong until you see that other weights make other tracts. **And the mirror
   case must pass**: a satellite cap, purity floor, score floor or bundle radius anywhere a reader
   likes is served, since none of them is a property of the file.
3. **The strata tally matches a direct count.** Enumerate the loci at a given floor from the catalog
   and by scanning; the per-stratum counts must agree exactly. This is the number
   [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.5 reweights by, so an error
   here propagates into a diversity estimate rather than into a crash.
4. **Size and time are measured, not assumed.** §6's statements are arithmetic. Report the file size
   and the row count by period, the pass's wall clock with and without the scan, and the observer's
   peak memory, on the tomato reference and on a human one. **This is what settles §9.1 and §9.3.**
5. **The thread count does not change the output.** Build the same reference sequentially and at
   several `--threads` values; the files must be byte-identical, which is the row order of §3.1 doing
   its job. **Include a reference whose contigs differ in size**, so the workers finish out of order
   and the writer has to put them back.
6. **A repeat longer than the satellite cap comes out as one row.** Put a 2 kb tract in a fixture: the
   file must hold it whole, at any thread count, and a reader must be free to call it a locus or a
   satellite (§4.2). This is the property §2.3's whole-contig scan exists for, and the one a chunked
   scan would silently break.
7. **The contig's two ends behave as specified.** A repeat 15 bp or further from either end is in the
   file; one closer than that is not, and that is the only place the catalog holds less than a live
   scan (§5.1). Test the kept case and the dropped case at each end.
8. **Staleness is caught, and the message says where.** Change one base on one contig and re-read the
   catalog: it must refuse **naming that contig**, not the reference as a whole. Then reorder two
   contigs of equal length — the lengths still match and only the per-contig MD5s catch it. Then
   change a scan setting; then the tool version. Each must refuse. A catalog whose header matches but
   whose body was truncated mid-write must not read as a valid short one.
9. **A present, valid catalog costs no detection.** Build once with `repeat-catalog`, then run a
   consumer that wants one: the consumer must do no scanning at all (§2.5), and its wall clock must
   sit with the no-catalog run of §10.4 rather than with the building one. A consumer run with **no**
   catalog present must stop and name the command, not scan and not proceed without one.
10. **Motif and purity match `classify`'s.** For every detection, the stored motif and purity must equal
   what `classify` computes from the bases — the same function, so this is a wiring test, and it is
   the one that fails if the coordinate conversion of §3.1 is off by one.
