# ng generic locus generator — Milestone D (D1)

**Plan:** [`locus_generation_pileup_generator.md`](../../ng/impl_plan/locus_generation_pileup_generator.md)
· **Spec:** [`locus_generation_pileup.md`](../../ng/spec/locus_generation_pileup.md) §3, §12, §13
· **Branch:** `ng-generic` · **Date:** 2026-07-29

D1 turns the differential around. Until now it compared `PileupRecord`s: ng's locus was laid
back out as production's older, smaller type by `to_pileup_record`, and everything B2 added
was merged or dropped on the way out. D1 projects **forward** — production's record said in
ng's terms — so every field ng carries that production cannot say becomes a difference with a
name, a count, and a test.

## 1. What landed

- **`project(&PileupRecord) -> SampleLocusObservations`** in `parity.rs`, total: `PileupRecord`
  and `AlleleSupportStats` are both destructured exhaustively, so a field added to production's
  types stops the file compiling rather than being quietly left out of the oracle. Every drop is
  a named divergence class, and the two drops that are *not* classes (`windowed_gc` /
  `windowed_coverage`, filled at the `.psp` seam; `placed_start`, which ng stops computing at
  B2) are bound and dropped by name at one site each.
- **`WalkOutcome` now holds ng's type on both sides.** `drive_production` maps each record
  through the projection; `drive_ng` projects nothing. Every comparison helper moved with it:
  `comparable`, `comparable_exact_q_sum`, `float_only_divergences`, `RecordEvidence` →
  `LocusEvidence`, `record_chain_ids` → `locus_chain_ids`, `classify_record` →
  `classify_locus`.
- **The permanent anchor**, `ng_agrees_with_production_where_production_fabricated_nothing`,
  on a new fixture (§3 below).
- **The census**, `every_divergence_from_production_is_one_of_the_six_named_classes`: every
  class counted, every class required to fire, and any difference outside them a panic.
- **Two unit tests the projection needed of its own**
  (`the_projection_says_everything_a_record_says`,
  `the_projection_orders_rows_as_the_walk_does`) — a differential cannot tell a right
  projection from a symmetrically wrong one.
- **`open_record::coverage_order` lifted to `pub(super)`** so the projection sorts rows with
  the walk's own comparator instead of a second spelling of it.
- **`to_pileup_record` is off the differential** (and off `mock_reference.rs`). It survives for
  the 44 inherited tests only; §6 carries that decision.

## 2. The finding that reshaped the step: the anchor's predicate was insufficient

Spec §3 and the plan both say the permanent anchor is *"loci where every folded read witnessed
the whole footprint must agree with production forever"*. **It does not hold**, and the anchor
found the counter-example on its first run:

```
seed 0x5eed0001 case 11 (complete), locus 30, region 42..49, reference AGTGTTTC
  ng    …  "AGTGTTTCC" Complete n=1   (7 rows, each n=1)
  prod  …  "AGTGTTTC"  Complete n=2   (6 rows; the REF allele holds two reads)
```

One read witnessed all eight positions **and an insertion inside them**. ng reports it
`Complete` with nine bases; production folded it into the REF bucket with eight. Production is
wrong about a read that witnessed everything.

**The cause is that `refold_live_reads` is ng's, not production's.** `grep` settles it:
`src/pileup/walker/open_record.rs` has `widen` and no re-fold; ng's copy has both (A3, option
(b)). Production's `widen` appends reference bases to every bucket and revisits no read — so a
record that widens *after* a read folded into it leaves production holding that read's
haplotype computed against a stale, narrower footprint. Both directions occur: reference bases
where the read's own belong, and events inside the widened region missed entirely.

This is **not** a special case of divergence class 1. Class 1 is "a read did not witness the
whole footprint"; here the read did. Filing it under class 1 would put reads production
*mis-folded* into the count of reads production *credited with bases they never sequenced* —
which is the deliverable — and spec §13.2 asks for those two numbers separately.

**So the census names six classes, not five.** Class 6 is production's stale widen. This is a
correction to spec §3's table, carried to Checkpoint D (§7) rather than edited into the spec.

## 3. The anchor, rebuilt on a fixture that contains no fabrication

With the predicate insufficient, the anchor needs a fixture where the *cause* is absent rather
than a filter over one where it is not. `generate_uniform_events` gives every read on a contig
**one shared event set** (one CIGAR, one start; bases, qualities, MAPQ, strand and pairing all
still vary). Then every read in a record carries its own copy of every widening event, so every
read is re-folded by that event and none is left stale. Production's append-to-every-bucket
still runs; it is simply overwritten, per read, by the read's own subtract-then-add.

It does **not** prevent widening, and an earlier draft of this fixture claimed it did — the
assertion `record_widen_events == 0` failed at 7 widens on the first case. A pileup opens a
record at every covered position, so a deletion anchored at one widens the record already
standing there. The fixture's real property is *who is left stale*, and what is asserted is:

- every read on a contig carries the same CIGAR and start (the fixture's defining property);
- widened records were reached at all (`widens > 0` — the class this anchor protects is a
  claim *about* widened records, so a fixture that stopped producing them would test the easy
  half);
- at every qualifying locus, ng equals the projection under `comparable`.

Measured, and recorded in the test so nobody mistakes which half is working: on this fixture
the every-read-`Complete` filter excludes **nothing** — 216,203 of 216,203 loci qualify —
because uniform events leave no way to be blind over part of a footprint. What is asserted is
therefore that the two walkers agree at *every* locus of a fabrication-free fixture.

## 4. The numbers

Default run (400 cases × 4 seeds, debug):

```
the anchor: 216,203 of 216,203 loci (100.0%) had every read witness the whole footprint,
  11,644 of them spanning more than one base, across 11,644 widened records; all identical
  to the projection, field for field. 103 agree only after q_sum is rounded to 1e-9.

the divergence census over 256,974 loci (one read group): 255,149 identical, 1,825 differing,
  863 of those agreeing once q_sum is rounded.
  class 1 partial witness 1,484   class 3 counters 34,106   class 4 unsupported bucket 27,891
  class 5 row order 34,683        class 6 stale widen 264
  at two read groups: class 2 on 190,542 loci, of which 105,324 carry one allele in two rows.
  the deliverable: production credited 2,787 reads over 1,484 loci with 8,239 reference bases
  they never sequenced (5.55 fabricated bases per fabricating locus).
```

Host-native soak, `PVC_PARITY_CASES=5000 cargo test --profile soak --lib …parity`, green:
3,262,582 loci censused, 2,700,954 anchored, class 6 at 3,074, deliverable 34,549 reads /
18,598 loci / 101,216 reference bases (5.44 per locus). Class 6's 264 at the default count is
the same population the retired `EvidenceIntact` class reported as "264 records (0.15 %)" —
which is corroboration that the sixth class is that class, now named instead of tolerated.

## 5. What the census asserts, beyond counting

- **Region, reference bases and kind are identical** at every locus. A2 moves rows; it does not
  move a record's existence or footprint.
- **The read accounting balances exactly**: ng's observations plus `reads_without_observation`
  equal production's observations. A read folds into exactly one bucket per record on
  production's side, so ng either emits it or counts it out — there is no third state, and
  anything else means evidence was created or lost rather than moved.
- **The per-`bases` evidence reconciles** wherever classes 1, 3 and 6 are all absent — spec §3's
  own wording for class 2, and what stops a split hiding a lost read.
- **Chain ids are compared by equality**, not by subset, wherever the reads landed in buckets
  with the same bases. The old record-level comparison asserted "ng's ids are a superset of
  production's"; **this census disproved that** at `seed 0x5eed0001 case 22, locus 37`. ng's
  rule is per read (did *this read* agree with the reference across everything it witnessed);
  production's is positional (`allele_index == 0`). The two coincide exactly when the read
  witnessed the whole footprint and both walkers bucketed it the same way — so that is where
  equality is asserted, which is a stronger claim than the subset it replaces.
- **Two passes, one and two read groups.** At one group class 2 cannot fire, so the
  unlisted-divergence panic has its full force; at two it fires almost everywhere, which
  exercises the class and the split but leaves that panic nothing to catch. Each pass asserts
  what it is in a position to assert — including that class 2 is **silent** at one group, so a
  bug tagging rows with arbitrary groups cannot look like the class working.

## 6. The eleventh test on this branch that could not fail

`ng_emits_no_allele_bucket_without_support` claimed A3's eviction. Its own doc named the
mutation it was written to catch — moving `evict_unsupported_alleles` above the contributor
fold loop, which strands every bucket that loop empties. D1 applied that mutation:

**199 tests passed. 198 in the module, plus the rest of the 2,720-test suite. Nothing caught
it.**

The reason is structural, and B2 introduced it: ng's rows are derived from `folded_reads`, per
read, so a bucket no read is folded into produces no row and leaves **no trace in the emitted
locus at all**. The test was meaningful while the walk emitted `PileupRecord`s and stopped
being so the moment it did not. The projection it ran through could not have saved it either —
`to_pileup_record` also derives rows from `observed_sequences`.

Fixed where the property lives: a `debug_assert!` at the top of `finalise`, over the allele
table itself. Every walk in the suite now checks it — the census alone runs it over ~257,000
loci, and the `soak` profile keeps it armed at release speed. With the assert in place the same
mutation fails immediately. The parity test keeps a narrower, real claim under a name that says
it (`every_emitted_row_carries_a_read`: the re-derivation cannot mint a row for nothing, which
is the B1 mistake it was written to avoid).

## 7. Mutation table

Every new test was mutated. Reverted after each.

| mutation | caught by |
|---|---|
| `evict_unsupported_alleles` moved above the fold loop | **nothing, before D1** — now `finalise`'s `debug_assert!`, via the census |
| `project` maps `fwd` onto `num_fwd`→`placed_left` | `the_projection_says_everything_a_record_says`, the anchor, the census, and the error differential (4 tests) |
| `sort_rows` drops the `ReadCoverage` tie-break | `the_projection_orders_rows_as_the_walk_does` |
| `coverage_of` always returns `Complete` | 9 tests, including the census (class 1 stops firing) |
| `finalise` reports `reads_without_observation` as 0 | 4 tests, including the census's read accounting |
| the anchor fixture stops sharing event sets (a fresh CIGAR per read) | the anchor's fixture-property assertion, by name |
| `finalise` reverses every row's `bases` — a corruption that is *not* the stale-widen shape | the census, as an **unlisted** divergence: *"production's row "CAT" is not any ng row's bases plus a reference tail"* |

The last one is the one that matters for class 6, which is the only class recognised from the
*difference* rather than read off ng's rows — an inherent liability, treated as one. Its shape
is required to be: the footprint spans more than one position; the bases genuinely fail to
reconcile; no evidence moved; and every production row is some ng row's prefix followed by a
suffix of the reference bases. A base corruption that does not fit is reported, not absorbed.

## 8. Recorded deviations

- **Six divergence classes, not five** (§2). The plan's D1 line says "every divergence falls in
  one of the **five** named classes"; the sixth was found by building the anchor and cannot be
  folded into any of the five without contaminating the deliverable. Raised at Checkpoint D.
- **The anchor's fixture is new** (`generate_uniform_events`) rather than `generate_complete`,
  because the spec's predicate alone does not select an agreeing class (§2, §3).
- **`open_record::coverage_order` became `pub(super)`** — a one-line visibility widening inside
  ng's own module, to avoid a second spelling of a comparator that must not drift.
- **`ng_walk_in_groups(case, groups)`** deals the fixture's reads round-robin into read groups.
  The read group reaches nothing but the row key in `open_record.rs`, so the walk is unchanged;
  without it class 2 is a branch nothing takes.
- **`to_pileup_record` is not deleted.** D1 retires the *hazard* — it is no longer the only view
  of the emitted type — and removes its two differential callers. Its remaining callers are the
  44 inherited tests and two `open_record.rs` fixtures, whose subject is the walk rather than
  the shape of the emitted type. Converting them is 49 assertion sites of hand translation
  against a B2 decision that was recorded with a reason ("67 hand-translated assertions are 67
  chances to re-express a test slightly weaker than it was"), so it is a decision and not a
  cleanup. **Carried to Checkpoint D with a recommendation: do it after D3**, when the dump tool
  and the throughput number are in and the inherited suite is the last consumer left.

## 9. Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib`: **2,722 passed**, 0 failed, 5 ignored (2,720 before D1; +12 net in
  `parity.rs`, minus the tests folded together).
- `cargo doc --no-deps`: still **12** pre-existing unresolved intra-doc links; none new.
- Host-native soak at 5,000 cases per seed, `--profile soak`: green (§4).
- `cargo test --all-targets --all-features` is not run: it panics in
  `benches/psp_writer_perf.rs:386` for reasons predating this branch.

---

# Milestone D — D2: the dump tool, and the two counters it forced

**Date:** 2026-07-30 · same plan, spec and branch.

D2's deliverable is `examples/ng_generic_loci_dump.rs` and the six fixtures spec §12 owes,
"because no inherited test exercises the defect". Building it forced two counters into
existence, released a fourth copied file, and found two of its own fixtures unable to fail.

## 1. The tool

Following `ng_ssr_loci_dump.rs`: a `#`-prefixed `key=value` counts header, a bare TSV column
line, one row per observed sequence. The pipeline is the real one end to end — `ReferenceInfo`
→ `SampleReads` → `TypedRegionIterator` → `GeneratorSet` with the **`Generic` slot filled** →
`SampleLocusObservationsIterator`. It is driven through the *dispatcher*, not the generator
directly, so the tool exercises the routing as well as the walk.

```text
# generic_loci=40 generic_regions=1 generic_region_bp=120 records_outside_region=0
# reads_admitted=2 reads_declined_by_preparer=0 reads_silent_over_footprint=0 reads_without_observation=0 reads_discarded_by_cap=0
# rows_complete=40 rows_observed=0 reads_complete=60 reads_observed=0
# record_widen_events=0 column_depth_truncations=0 regions_in=1 regions_handled=1 loci_emitted=40
contig  start  end  ref_bases  depth  read_coverage   read_group  observed  reads  chain_ids
chr1    19     24   CCGTTA     2      observed:0+5    0           CCGTT     1
```

The three surfaces waiting for a first caller all have one now:

- **`WindowedRefSeq::with_shared_index`** — the generator needs two accessors and calls the
  factory at every region; a fresh `new` parses the whole `.fai` on its first fetch (189 µs on
  a GRCh38-shaped reference against a ~120 µs per-region walk constant). The tool parses the
  index once and shares it, which is the 19 µs path.
- **`GeneratorSet::generic_counts()`** — the boxed generator's nine counters, read through the
  dispatcher. The tool errors out if the slot reports `None`, so "the generic slot is the
  filled one" is checked rather than assumed.
- **The `#`-prefixed counts header** — four lines, and every one of them is what makes spec
  §13's read accounting checkable.

`read_coverage` prints the **run** (`observed:<offset>+<positions>`) rather than a side label:
on the generic path a read can be blind in the *middle* of a footprint, where `partial:left`
would be a lie.

One knob, `PVC_GENERIC_REGION_CHUNK_BP=n`, splits each `Generic` region into *n*-bp pieces
before generating. It is a real knob — in the pipeline the regions are handed in and their size
is the caller's — and it is what makes the region-boundary fixture a comparison of two
configurations over the *same* reads rather than two different fixtures.

## 2. The two counters spec §13's read accounting forced

Both were named in Milestone C's report as structurally zero. They are live now.

- **`reads_silent_over_footprint`.** `ActiveRead` gains `ever_contributed: Cell<bool>`, set in
  the contributor loop; `ActiveReads` tallies the reads that leave without it, at both exits;
  `RunSummary` carries the total and `fold_region_walk` sums it. A `Cell` because the
  contributor loop walks the set through `iter()` and queries each read's cursor inside that
  borrow. Set **before** the mate-overlap collapse and before the depth cap, so a read either of
  those removed is counted where it belongs and not here as well.

  This is the read *neither* per-locus counter can see: every base `N` or adaptor-masked, so it
  produced no observation and never reached the fold that records
  `reads_without_observation`. Before the flag, the only honest thing to say about it was
  nothing.

- **`reads_declined_by_preparer`.** A counter on the shared `ReadPreparation` cell, taken by
  `end_walk` like the shed error beside it. It reads **zero on every real run** — no v1 preparer
  declines anything, the only step that could was BAQ, deferred — and that is why it exists: a
  read the preparer declines never reaches the walk, so `reads_admitted` cannot account for it,
  and "no preparer declines today" is a fact a counter can state and a missing field cannot.

`active_read_set.rs` is therefore **released from `copy_fidelity.rs`** (the fourth of eight), in
this commit, with the release recorded in its own header and in the guard's table — the pattern
the plan set at A0. The guard caught the divergence on the first run, which is what it is for.

## 3. The six fixtures

Each is written so a **span-derived** or **fill-preserving** implementation fails it.

| # | fixture | what only the events know |
|---|---|---|
| 1 | a read adaptor-masked over part of a five-position record | its alignment spans the whole footprint, so a span-derived coverage says `Complete` |
| 1b | the same record carries `Complete` and `Observed` rows **with the same bases** | one read witnessed five positions and deleted three; the other witnessed two and was silenced — production cannot tell them apart |
| 2 | an interior `N`, and a ref-skip, inside a footprint | the witness is non-contiguous, so there is no observation at all and the read is counted |
| 3 | a record widened past an **expired** read | production appends the reference bases to a bucket whose read is gone; ng's row still says five |
| 4 | an indel whose own anchor base was masked | the one recorded residual — exactly one borrowed reference base |
| 5 | two read groups on one allele | two rows, summing back to the one-group total per locus |
| 6 | a deletion across a region boundary | the halo, checked against the same fixture walked as one region |

Plus the four whole-fixture properties: one locus per covered position and none for an uncovered
one; the read accounting; the per-read chain-id rule; and byte-identity across runs. `push_locus`
asserts the global consistency check on **every** run of the tool, not only in tests — no row
claims more locus positions than its events account for, which is deliberately *not* an equality
(an insertion adds bases without positions, a deletion the reverse).

## 4. Three things the fixtures got wrong first, and how

Worth recording because each is a way a fixture can look right and test nothing.

**The real read filter drops reads under 30 bp.** The first draft's indel reads carried 26
sequenced bases, so `reads_admitted` was 3 of 4 and two fixtures failed for a reason that had
nothing to do with the walk. A read *aligned* over 30 positions cannot end inside a 16-position
footprint either, which is why fixture 3's short read carries a 17-base soft clip: the minimum
counts sequenced bases.

**ng's real read preparation left-aligns the fixtures' deletions.** A deletion at 21..25 whose
preceding base equals its last deleted base shifts to 20..24 — so the record opened one position
earlier than the test asserted and there was nothing at all where it looked. Every deletion in
these fixtures is now chosen so that `ref[anchor] != ref[anchor + len]`, which is the condition
under which no left-shift is available; fixture 3 uses one deliberate shift and says so.

**A hole outside every footprint costs a read nothing.** Fixture 2's first `N` sat at a position
no record spanned, so `reads_without_observation` stayed zero while the test asserted it did
not.

## 5. Two of the fixtures could not fail, and mutation is what said so

| mutation | before | after |
|---|---|---|
| `coverage_of` always returns `Complete` (a span-derived implementation) | fixtures 1 and 3 fail — as spec §13 requires. *(Milestone D's review re-ran this and counted **four** failing dump fixtures, not two: the region-boundary and chain-id ones fail as well. Wider coverage than claimed, not narrower.)* | unchanged |
| the chain-id rule becomes production's positional `allele_index == 0` | **all ten passed** | fixture rewritten; now fails |
| the halo removed from the region query | **all ten passed** | fixture rewritten; now fails |
| the read group dropped from the row key | fixture 5 fails | unchanged |
| `ever_contributed` never set | the silent-read test fails | — |
| the declined tally never incremented | the declined-read test fails | — |

**The chain-id fixture** had a matching read and a deletion read — both whole-footprint
witnesses, and the per-read rule and the positional rule *coincide* on those. The case that
separates them is a **partial** row whose bases are a prefix of the footprint's reference bytes:
ng gives it no id (it agreed with the reference across everything it witnessed) while the
positional rule tags it, its row not being `alleles[0]`. The test now runs on the masking fixture
and asserts that row directly, plus a whole-footprint match and a genuine departure.

**The boundary fixture** chunked a read set that all started at 11 — inside the first piece —
and a read overlapping the unwidened query span is fetched either way. The halo only matters for
a read that begins **after** the region's end and still folds into a record anchored inside it,
so the fixture now adds one starting one base past the boundary, and asserts its row is present
before comparing the two configurations.

That makes twelve tests on this branch that could not fail, three of them in this milestone.

## 6. Recorded deviations

- **The tool is generic over the preparer**; `main` uses ng's real `LeftAlignPreparer` and the
  tests wrap it, setting an adaptor boundary on the incoming `AlignedRead` by qname. Deriving the
  boundary from mate geometry *is* production's code, tested there; a BAM whose geometry happened
  to land the boundary where a fixture needs it would make the fixture a test of that derivation
  and would move whenever it changed. Everything after the boundary is the real preparation.
- **One region stream for the whole run**, with the per-contig typed-region walks chained,
  because the locus iterator owns the stream, the reads *and* the generator set — and the
  generator must not be rebuilt per contig, its chain-id allocator being run-lifetime.
- **`#[allow(clippy::arc_with_non_send_sync)]`** at the one site that wraps a `WindowedRefSeq`:
  the generator's `new` takes `Arc<R>`, chosen at C1 when the only accessor anyone had was an
  in-memory one that is `Send + Sync`. A file-backed accessor holds a reader behind a `RefCell`
  and is neither, and the walk is single-threaded anyway (arch §9), so **`Rc` would say what is
  true** — carried to Checkpoint D as a finding rather than changed from an example.

## 7. Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib`: **2,724 passed** (2,722 after D1; +2 for the counters).
- `cargo test --example ng_generic_loci_dump`: **10 passed**.
- `cargo doc --no-deps`: still 12 pre-existing unresolved links.
- Host-native soak, `PVC_PARITY_CASES=3000 --profile soak`: green — the walk's new `RunSummary`
  field is dropped by name in `parity.rs`'s counter comparison, production having no counterpart.

---

# Milestone D — D3: the two measurements

**Date:** 2026-07-30 · same plan, spec and branch. Both numbers are deliverables, not
by-products (plan D3), and both are taken **host-native in release** — the parity harness's
environment variables do not reach the container, and production's reachable `debug_assert!` on
overlapping mates carrying deletions would abort a debug walk almost immediately.

## 1. The size of production's defect, on real reads

Run through the parity module's `#[ignore]`d real-data differential, which since D1 prints the
six-class census and the three-number deliverable over the same code path the synthetic census
uses.

| data | region | loci | reads | production credited … |
|---|---|---|---|---|
| GIAB HG002 **300×** BAM | chr1:1,000,000–6,000,000 | 48,905 | 55,054 prepared | **871 reads over 162 loci (0.33 %) with 1,550 reference bases they never sequenced** |
| GIAB HG002 **10×** BAM | chr1:1,000,000–1,200,000 | 1,390 | 75 prepared | 0 — the window carries no partial witness at all |
| tomato **CRAM** (SRR7279481) | SL4.0ch01:1,000,000–6,000,000 | 96,260 | 4,616 prepared | 6 reads over 4 loci (0.00 %) with 58 reference bases |

Class counts on the 300× run: partial witness 162, group split 0 (one read group), counters 35,
unsupported bucket 60, row order 9,561, **stale widen 0**.

**What the number says.** The defect is real and it is *small* on this data: about one locus in
300 at 300×, and 9.6 reference bases per affected locus. It is also concentrated exactly where
the port predicted — the loci with long deletions — and the 10× row is the important control:
a shallow window over the same coordinates produces **zero**, because a partial witness needs a
read that stops inside a multi-position footprint, and at 10× over a tandem-repeat benchmark
there are barely any multi-position footprints to stop inside.

**On the indel-deficit hypothesis, this is not yet a verdict.** The measurement the hypothesis
needs is over the loci where production's indel calls go missing, and this is a 5 Mb window of
one chromosome of a tandem-repeat benchmark. What D3 establishes is the *method* and an order of
magnitude: the fabrication exists, it is countable, and at 0.33 % of loci it is a candidate
contributor rather than an obvious cause. Class 6 being **zero** on real data is worth as much:
production's stale widen, which the synthetic census fires 264 times in 257,000 loci, needs a
record to widen *after* a read folded into it, and real 300× data over 5 Mb produced none.

## 2. Throughput against production's `pileup`, one human chromosome

Same reference, same BAM (HG002 TR 30×, so both walk the same reads), chr1, **single-threaded on
both sides** (`--threads 1`; ng's parallelism is deferred whole, so the comparison would
otherwise be against four workers).

```
ng   ng_generic_loci_dump … chr1 > file      34.70 s wall   461 MB peak RSS   72 MB of TSV
ng   the same, output to /dev/null           32.59 s wall   461 MB peak RSS
prod pop_var_caller pileup --threads 1       10.20 s wall   559 MB peak RSS   27 MB .psp
```

**ng is 3.4× slower, or 3.2× discounting the TSV it writes and production does not.** Peak RSS is
**lower** than production's, which is the memory property spec §7 asked for holding at
chromosome scale.

### Where the time goes, from the counters and one throwaway probe

The dump's own header decomposes the run:

```
generic_loci=1541788  generic_regions=613682  generic_region_bp=240227974  records_outside_region=883083
reads_admitted=374437  record_widen_events=2514  column_depth_truncations=0
```

- **Region typing is 5.92 s of the 32.59 s (18 %) — a pass production does not make at all.**
  Measured with a throwaway probe that drives `TypedRegionIterator` over chr1 and counts,
  nothing else; deleted, per this branch's practice. It types 1,227,363 regions, 613,682 of them
  `Generic`.
- **The regions average 391 bp**, and Milestone C's review measured a region at ~0.12 ms to set
  up — *290 loci of walking*, with the conclusion that **regions under ~300 bp cost more to open
  than to walk.** This run is squarely in that regime. (The 0.12 ms constant was measured before
  `with_shared_index`; at 26.7 s of non-typing time over 613,682 regions the real figure here is
  ~43 µs per region *including* its walk, so the fix landed at C's review is already paying.)
- **883,083 records — 36 % of the 2.42 M the walk finalised — are discarded at the region
  clamp.** That is halo work: with a 391 bp region and a 5,000 bp halo, every region walks far
  past its own end, and the stop rule keeps that bounded rather than free.
- **The allele list is not the problem on this data**, which is worth saying because spec §7
  named it as the first place to look: `column_depth_truncations` is 0 and `record_widen_events`
  is 2,514 over 1.5 M loci, so records are essentially never deep and essentially never widen.

**So the first lever is the region grain, not the fold.** Spec §7's rule applies — a bad number
is a performance problem to solve, not a design to reconsider — and the two candidate fixes are
both outside this plan: hand the generator **fewer, larger** `Generic` regions (their size is the
caller's choice, and the dump's `PVC_GENERIC_REGION_CHUNK_BP` already shows the output is
invariant to it), and stop re-entering a region for work the previous one's halo already did.
Parallelism, deferred whole, is a third: production's default is four workers, where this
comparison used one.

## 3. What the measurement found in the differential itself

**The `q_sum` comparison's tolerance was a fixed absolute grain, and 300× data broke it.** The
first deep run failed as an *unlisted divergence* — the census doing its job — on a locus with
414 observations:

```
locus 919: ng's rows do not sum back to production's per-allele totals, and none of a
partial witness, a counted-out read or a stale widen explains it
  ng   {[84]: … q_sum_rounded: -3360392684715 }
  prod {[84]: … q_sum_rounded: -3360392684716 }
```

One grain apart. The comparison rounded both sides to 1e-9 **absolute**, on an argument that was
really a statement about depth: "the grain is nine decimal places on values of order −3 to −50,
where the smallest real difference is a whole read's `ln` contribution — order 1". At 300× a
locus's `q_sum` is order −3,400, where 1e-9 absolute is 3.4 × 10¹² grains and a reordered
accumulation lands one apart.

Fixed twice over:

- **Relative, not absolute.** `Q_SUM_TOLERANCE = 1e-9` is now a *relative* allowance,
  `|a − b| ≤ ε · max(1, |a|, |b|)`, so the headroom against a real one-read difference stays six
  orders of magnitude at any depth instead of only at low ones.
- **A tolerance, not a rounding.** Rounding decides equality by which side of a grain boundary
  each value falls on, so two sums one ulp apart can round *apart* however fine the grain — at
  millions of loci that happens. `ComparableLocus` and `LocusEvidence` now carry hand-written
  `PartialEq`s: every integer field exactly, `q_sum` within the tolerance, and both destructure
  exhaustively so a field added to either type stops the comparison compiling.

Both loci types' comparisons therefore no longer use `SampleLocusObservations`' derived
`PartialEq`, which compares `f64`s bit for bit — right for the shared type, wrong for a
differential whose subject is two implementations summing the same addends in different orders.

## 4. Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib`: **2,724 passed**, unchanged by the tolerance fix.
- The three real-data runs above are the measurement; the perf artefacts and the typing probe
  are deleted, the numbers living here.
