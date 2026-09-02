# ng — STR observations through a run: the routing policy, the per-sample walk, and the cohort merge

*2026-09-02. No code yet — this settles the design. Companion documents:
[`calling_loop_ssr.md`](calling_loop_ssr.md) (what consumes the observations this document
produces), [`locus_generation_ssr.md`](locus_generation_ssr.md) (the tract generator, built),
[`cohort_merge.md`](cohort_merge.md) (the merge, built), [`typed_regions.md`](typed_regions.md)
(classification, built), [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) (selection,
settled on paper, whose §2 named the merge gap this document closes).*

*What prompted it:
[`ng_str_path_losses_2026-09-02.md`](../../reports/ng_str_path_losses_2026-09-02.md). At 30×
and 50× on the GIAB benchmark, every truth variant on ground the run classifies as repeat is
missed — recall exactly 0.000 against the production caller's 0.855–0.990 — and the run
classifies about seven times more of the reference as repeat than ng's own calling floors
would.*

---

## 1. What this is

A calling run today produces observations — *what did this sample's reads show here* — for
ordinary sequence only. This document specifies the changes that make a run produce **two
kinds of observation, per sample and then per cohort**: the SNP/indel kind it already
produces, and the repeat-tract kind, so that the calling loop can be handed either and told
which it holds.

Three changes, in dependency order:

1. **The routing policy becomes the caller's, with the catalog as a candidate source only**
   (§2). A run classifies a stretch as an STR locus when it clears thresholds the *user* set,
   defaulting to ng's measured calling floors — never merely because the catalog file holds a
   row there.
2. **The tract slot in the walk's generator set is filled** (§3), so each sample's walk emits
   repeat-tract observations beside the SNP/indel ones. The generator exists; this is wiring
   plus two accounting obligations.
3. **The merge carries the locus kind to its output** (§4), so a cohort observation states
   whether it is a tract — and, when it is, its motif — instead of dropping both.

### 1.1 Goals

- A run's definition of *an STR locus* is a parameter with a measured default, satisfiable
  from the CLI, and recorded in the run's parameters file — beside what any fitted file it
  was handed recorded — so what a run routed with is always on the record.
- Every sample's walk yields tract observations wherever the routing says *tract*, with the
  same per-read-group accounting the SNP/indel path already reports.
- The merge emits cohort observations of both kinds, each carrying its kind, over ground that
  partitions exactly as it does today.
- Nothing on the SNP/indel path changes byte-for-byte where the routing did not change the
  ground it sees.

### 1.2 Non-goals, and what this document does not do

- **It does not call anything.** What the calling loop does with a tract observation —
  selection, likelihood, prior, quality, the record — is
  [`calling_loop_ssr.md`](calling_loop_ssr.md).
- **It does not give bundles a caller.** A bundle — a cluster of repeats none of which has
  clean flanks — keeps its unfilled slot and stays counted as ground the caller cannot yet
  speak for (§8).
- **It does not move the routing frontier.** Which repeats *ought* to go down the STR path —
  the period × length question [`typed_regions.md`](typed_regions.md) §5.2 leaves open, with
  a graft onto the generic path as a live third answer — stays open. This document makes the
  frontier the user's knob and sets its default; it does not measure where it belongs.
- **It does not touch the parameter pre-pass**, whose own selection already asks the catalog
  with its own criteria ([`SelectionTerms`,
  `joint/census.rs:509`](../../../../src/ng/parameter_estimation/joint/census.rs)) — except
  to require that the two sets of criteria agree (§2.3).

### 1.3 Vocabulary

- **The catalog** — the parquet file of every tandem repeat found in the reference, built
  once by `repeat-catalog` at deliberately low floors
  ([`repeat_catalog.md`](repeat_catalog.md)): 5 copies for a homopolymer down to 3 for a
  hexamer, tracts to 500 bp. **A row in it is a candidate, not a finding.**
- **An STR locus** — a candidate that clears the *caller's* thresholds. Only these are routed
  to the tract generator.
- **A locus generator** — the code that turns one stretch of reference plus one sample's
  reads over it into a locus: the distinct sequences the reads showed, each with its
  support. One slot per kind of stretch
  ([`GeneratorSet`, `locus_generation/mod.rs:841`](../../../../src/ng/locus_generation/mod.rs)).
- **Ground** — the stretches of reference sequence a run analyses, measured in bases; the
  run report's own word (*"analysed ground"*).
- **Observation** — per sample, a
  [`SampleLocusObservations`](../../../../src/ng/locus_generation/mod.rs) (what one
  sample's reads showed at one locus); per cohort, a
  [`CohortObservation`](../../../../src/ng/run/cohort_merge/build.rs) — the merge's
  unification of every sample's, and the thing calling consumes.

---

## 2. The routing policy is the caller's; the catalog only supplies candidates

**Decided 2026-09-02 (owner):** *"the catalog should be considered only a source of possible
STR loci. For us an STR locus is only what is above the threshold given to the caller by the
user. […] What's in the catalog, although it is technically a repeat, is not an STR locus for
us — only a candidate to be a possible STR locus."*

This was always the file's design —
[`repeat_catalog.md`](repeat_catalog.md) builds below every calling floor precisely so that
*"a caller can move its routing floor anywhere inside that gap by filtering, which is the
question the file exists to keep open"* — and the run does not honour it:
[`call_from_alignments.rs:845`](../../../../src/pop_var_caller_exp/call_from_alignments.rs)
asks the catalog with `StrRepeatCriteria::default()`, which *is* the file's own storage
floors. Everything the file holds becomes an STR locus of the run.

**What that costs, measured** (GIAB, three samples' confident regions,
[the report](../../reports/ng_str_path_losses_2026-09-02.md)): 32,577 of HG002's 572,037
analysed bases routed as repeat against 4,930 under the calling floors — 6.6×, 7.3× and 7.0×
across the samples; 2,753 of HG002's 2,992 routed tracts are sub-floor repeats, mostly
five-base homopolymer runs; and the extra admissions collide into bundles, 656 of them
against 17, which is where 45 of the benchmark's 55 lost truth indels sit.

### 2.1 The knob

The run's routing criteria are a
[`StrRepeatCriteria`](../../../../src/ng/repeat_catalog/criteria.rs) of its own — the same
type the catalog reader already takes, so nothing downstream changes shape:

- **period range** (default 1–6),
- **per-period copy floors** (default `[8, 6, 6, 6, 5, 4]` —
  [`MinCopies::default`, `segment_criteria.rs:402`](../../../../src/ng/region_typing/segment_criteria.rs),
  **measured** over the tomato archive on 2026-08-10 as the copy counts at which a repeat
  starts to stutter; below them, that type's own documentation says, *"the generic SNP/indel
  caller handles the tract"*),
- **purity floor** (default 0.8 — inherited, never measured),
- **satellite cap** `max_str_len` (default 100 bp — inherited from production, never
  measured),
- **minimum flank** (pinned at the catalog's 15 bp; a reader cannot ask below what the file
  was built at).

Exposed on `call-from-alignments` as flags named exactly as `type-regions` names them
(`--min-copies`, `--min-period`/`--max-period`, `--max-str-len`, `--min-purity` —
[`typed_regions_cli.md`](typed_regions_cli.md); resolved, §9). The defaults are the calling
floors, so a run that types nothing routes as ng means to route.

**The plumbing is one line by construction**: `segments_over` already passes a
`StrRepeatCriteria` to `catalog.genome_segments`, and everything the criteria feed —
admission, bundling, whole-tract emission at BED edges
([`segments_of_contig_in`, `repeat_catalog/segments.rs:72`](../../../../src/ng/repeat_catalog/segments.rs))
— is built and shared with the development tools, which already ask with the calling floors
(`examples/shared/catalog_regions.rs` and every STR dump).

### 2.2 A candidate below the floors is generic sequence, not a hole

Already the classification's rule —
[`typed_regions.md`](typed_regions.md) §2.2, *"a rejected repeat is generic territory, not a
hole"* — restated here because it is what makes §2.1 recover variants rather than merely
shrink a number: a sub-floor candidate's bases fall back to the generic path, whose locus
generator is built and calling at 0.984/0.982 recall (SNPs/indels, 30×). The catalog-reading
path already implements the fallback
([`repeat_features_of_contig`](../../../../src/ng/repeat_catalog/segments.rs), the same
admission the scan runs).

### 2.3 The criteria are recorded, and agreement with a fit is the user's call

Two different obligations, and only the first is enforcement:

- **Against the catalog**: a reader asking for less than the file was built at is refused —
  built, [`StrRepeatCriteria`'s admissibility check,
  `criteria.rs:101`](../../../../src/ng/repeat_catalog/criteria.rs). This one is not policy:
  the rows below the file's floors were never written, so the request cannot be served.
- **Against a fitted parameters file — the user decides, the run records. Decided 2026-09-02
  (owner):** *"The user could supply any priors to the caller; it is their decision to use
  the same routing criteria or not."* The pre-pass already records what it asked the catalog
  for (`SelectionTerms::ssr_criteria`, *"what this run asked the catalog for — a reader
  chooses its floors freely within what the file was built at"*), and the parameters file
  carries that census identity as named, comparable terms
  ([`bindings.rs`](../../../../src/ng/calling/parameters_file/bindings.rs)). The calling
  run's own routing criteria are written into the parameters file it saves beside its VCF,
  so a run whose routing differed from its fit's is a fact on the record — visible, never
  blocking. The trade the user is taking when the two differ is worth one sentence here: a
  tract admitted by the run but outside the fit's selection is scored from strata fitted
  over other loci, or from the stated defaults — which the per-cell warrants already label.

---

## 3. Per-sample tract observations: fill the slot, and pay two accounting debts

**The generator exists and works.**
[`SsrGenerator`](../../../../src/ng/locus_generation/ssr.rs) turns one tract segment into one
locus — reads fetched over the tract, each aligned to read off its repeat length, tallied
into a `SampleLocusObservations` whose `kind` is `LocusKind::Ssr(SsrDetail)` with the
motif and both flanks
([`locus_generation/mod.rs:429`](../../../../src/ng/locus_generation/mod.rs)). Demonstrated
on real data 2026-09-02: at 30×, 18–29 length-pinning reads per covered tract. What a run is
missing is the wiring:
[`generic_path_generators`, `walker.rs:1588`](../../../../src/ng/run/walker.rs) builds
`GeneratorSlot::Unfilled(NotImplemented)` for both repeat slots.

### 3.1 Construction

One `SsrGenerator` per sample's generator set (a generator carries per-segment state and
scratch and is never shared — the run spec's own rule), built beside the pileup generator
from the same `WalkReference` accessors:

- **aligner**: the unit-robust delimiter, algorithm 4u — the bake-off winner and
  [`ng_ssr_loci_dump`](../../../../examples/ng_ssr_loci_dump.rs)'s default. Not a knob here;
  the bake-off is recorded and a run gets its winner.
- **config**: [`SsrGeneratorConfig::default()`](../../../../src/ng/locus_generation/ssr.rs)
  — flank 15 bp (must satisfy `flank_bp <= bundle_threshold`, the constructor's own check),
  read cap 1,000 per locus (production's number, **soft, never measured**).

### 3.2 The two accounting debts, both already on record

- **Read-filter tallies.** The trait's `read_filter_counts` defaults to empty and
  `SsrGenerator` does not override it
  ([`locus_generation/mod.rs:507`](../../../../src/ng/locus_generation/mod.rs);
  the 2026-09-01 run-report work recorded this as the thing that *"will bite whoever fills a
  tract slot"*). The generator keeps its own reader and retired counts; the override is owed
  in the same change that fills the slot, or every run's per-read-group drop rates
  under-report silently — plausible numbers, no error.
- **Generator counts.** `GeneratorSet` already sums a filled slot's `LocusCounts` and exposes
  per-kind counts (`ssr_counts`), so the run report's ground partition moves tract bases from
  *not called — repeat tracts this caller has not built yet* into *called* with no counting
  change. The report's **wording** is the obligation: it must read truthfully in both states —
  a run with the tract slot filled and bundles still unfilled prints a smaller, honest *not
  called* line, not a zero and not the old sentence.

### 3.3 What does not change

The walk's contract is untouched: segments are the reference's partition, no observation
crosses a segment boundary, one kind per cohort locus is asserted where loci chain
([`close.rs:636`](../../../../src/ng/run/cohort_merge/close.rs)), and the tract's exemption
from the span cap is already exhaustive on the kind
([`judge`, `close.rs:110`](../../../../src/ng/run/cohort_merge/close.rs)). Concurrency
invariance (the E2 oracle — byte-identical output at every thread count) must hold with the
slot filled; nothing in the tract generator reads anything another sample wrote, so nothing
new crosses threads.

---

## 4. The merge: cohort observations of both kinds

**The gap, named by [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §2:** the merge
closes tract loci deliberately, but
[`CohortObservation`, `build.rs:974`](../../../../src/ng/run/cohort_merge/build.rs) carries
only the region, the allele table and the per-sample rows. The kind — and with it the motif,
so the period, so the repeat counts — is dropped at assembly.

**The fix:** the kind travels from the per-sample observations to the cohort observation.

- `ClosedLocus` exposes the locus kind. Its members all carry one (`SampleLocusObservations::kind`)
  and the walk already asserts they agree by discriminant; the closed locus states the shared
  kind once instead of leaving the reader to take the head member's.
- `CohortObservation` gains `pub kind: LocusKind`, cloned from the closed locus at
  `CohortObservation::over`. The clone is a motif (inline, ≤ 6 bytes) and two boxed 15-byte
  flanks per *tract* locus only; `LocusKind::Generic` is free.
- Everything already right stays: per-`(allele, read group)` support rows (the stutter model
  is per read group and needs them apart — §2's own second finding), partial observations
  kept apart from alleles rather than padded into them
  ([`build.rs:1050`](../../../../src/ng/run/cohort_merge/build.rs)), the variability filter,
  and the note already booked in the merge for tracts longer than a read
  ([`build.rs:1025`](../../../../src/ng/run/cohort_merge/build.rs)).

**Rejected: a parallel `SsrCohortObservation` type.** Two types would force every holder
between the merge and the caller — the cache, the parallel cover, the driver — to carry two
streams and re-interleave them in genome order, which the one-kind-per-locus assertion makes
pure cost. The kind is a field because a cohort locus *has* exactly one.

---

## 5. What this hands on

After these three changes a run yields, in genome order, cohort observations each stating
`Generic` or `Ssr(detail)` — and the driver still calls only the first:
[`call_one_generic_locus`, `callers.rs:813`](../../../../src/ng/run/callers.rs) is
unconditionally the SNP/indel path. **The branch on the kind is deliberately not here.** It
needs the tract's candidate selection to have something to dispatch *to*, and both belong to
[`calling_loop_ssr.md`](calling_loop_ssr.md). Until that lands, a tract cohort observation is
built and then set aside, counted — which is already an honest improvement over never
building it, and keeps the two documents independently buildable.

**Build order this implies:** §2 alone is shippable first and recovers about four in five of
the currently lost truth variants with no new calling code (routing them to the generic
path). §3 and §4 are shippable together behind it, inert to the VCF until the loop's branch
lands, visible in the run report's counts.

---

## 6. Cross-cutting concerns

- **Memory.** Segments: a tract segment's payload (motif + coordinates) already exists for
  every catalog row admitted; *raising* the floors shrinks the resident segment list
  (~100,000 typed regions on the tomato benchmark's 8 Mb today). Cohort observations: +1
  `LocusKind` per locus; tract loci add two 15-byte flanks each. The merge's peak-heap
  attribution is unchanged in kind.
- **Determinism.** The routing criteria are inputs to the segmentation, which is computed
  once and shared across samples; same inputs, same segments, and the E2 byte-identity
  oracle must pass with the tract slot filled.
- **Errors.** No new failure modes: criteria below the catalog's floors are already a typed
  refusal, and a criteria difference against a fitted file is not an error at all — it is the
  user's choice, on the record (§2.3).

## 7. Reuse map

| what | existing code | reuse |
|---|---|---|
| criteria type, admissibility vs the file | [`StrRepeatCriteria`](../../../../src/ng/repeat_catalog/criteria.rs) | as is; the run builds one from flags instead of `default()` |
| classification at chosen floors, bundling, whole-tract emission | [`repeat_catalog/segments.rs`](../../../../src/ng/repeat_catalog/segments.rs) | as is — already parameterised, already what the dev tools call |
| the tract generator | [`SsrGenerator`](../../../../src/ng/locus_generation/ssr.rs) | as is, plus a `read_filter_counts` override |
| per-sample → cohort collation, partials, kinds asserted | [`cohort_merge/{close,build}.rs`](../../../../src/ng/run/cohort_merge/build.rs) | one field added, nothing reshaped |
| criteria in the run's recorded identity | [`SelectionTerms`](../../../../src/ng/parameter_estimation/joint/census.rs), [`parameters_file/bindings.rs`](../../../../src/ng/calling/parameters_file/bindings.rs) | the calling run's criteria join the comparison |
| oracle for the routing | [`examples/ng_typed_region_dump.rs`](../../../../examples/ng_typed_region_dump.rs) | reproduces a run's routing byte-for-byte; grew the `catalog|calling` switch for exactly this measurement |

## 8. Deferred, with a recommended home

- **A caller for bundles.** Under the calling floors bundles hold 4 lost SNPs and 5 lost
  indels of the GIAB benchmark's 2,391 truth variants, so nothing here blocks on them —
  but they are the one kind with neither a generator nor a design. Home: their own spec,
  after this document and [`calling_loop_ssr.md`](calling_loop_ssr.md) land; the payload
  `LocusKind::SsrBundle` deliberately carries nothing until then.
- **Where the routing frontier belongs** — the period × length measurement, including the
  DRAGstr-style graft option. Home: [`typed_regions.md`](typed_regions.md) §5.2's question,
  answerable cheaply once both paths run (call the same ground both ways, score both against
  truth).
- **Satellites** stay a permanent, counted refusal; unchanged.

## 9. Resolved decisions and open questions

- **The catalog is a candidate source; the caller's thresholds define an STR locus — decided**
  (owner, 2026-09-02; §2). Supersedes the run's routing on the file's floors, and retires
  [`typed_regions.md`](typed_regions.md) §5.2's *"v1 = the catalog's set"* scaffolding as any
  run's behaviour.
- **The kind is a field on the one cohort-observation type, not a second type — decided**
  (§4; confirmed by the owner, 2026-09-02).
- **Routing criteria that differ from a fitted file's are the user's decision, not a
  refusal — decided** (owner, 2026-09-02; §2.3). The run records both sets and calls on.
- **OPEN: should period-1 tracts route to the STR path at all at the default floors?** The
  calling floor says 8+ copies stutter; the generic path demonstrably handles ≤ 7. Leaning:
  keep period 1 in the default range and let the frontier measurement (§8) rule; confirm
  before changing the default range.
- **The CLI spelling is flat flags, named as `type-regions` names them — resolved
  2026-09-02** (`--min-copies`, `--min-period`, `--max-period`, `--max-str-len`,
  `--min-purity`): one vocabulary for the same five knobs across the tool that previews a
  classification and the run that calls on it, and no criteria file to version. Settled here
  rather than left to the plan, which may not decide design.

## 10. How we know it works

- **Routing parity:** `ng_typed_region_dump` at the run's criteria matches the run's own
  report bases exactly (it does today: 539,460/32,577 on HG002 — keep that as the oracle as
  the criteria become flags).
- **Recovery:** on GIAB at 30×, the truth variants on ground that moves generic must be
  found at the generic path's own rate — the report's upper bound predicts overall recall
  ≈ 0.97 SNPs / ≈ 0.94 indels from routing alone; measure against it.
- **No regression:** where routing is unchanged (ground generic under both criteria), the
  VCF is byte-identical.
- **Slot filled:** the run report's ground partition keeps summing to 100%; per-read-group
  drop tallies are non-zero on a run with tract ground (the §3.2 debt's test); E2
  byte-identity at 1–16 threads still passes.
- **Merge kind:** a tract cohort observation reaching a probe carries the same motif the
  generator minted, on real data (`ng_ssr_loci_dump` ground truth beside the merge's
  output).
