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
