# ng — what a locus observation says a read witnessed

**Status:** draft, 2026-07-30. **No code yet — this settles the design.** Code-facing companion:
[`../arch/locus_witness_representation.md`](../arch/locus_witness_representation.md) (the types and
signatures); build order in
[`../impl_plan/locus_witness_representation.md`](../impl_plan/locus_witness_representation.md).
Changes the **shared** locus type
([`locus_generation.md`](locus_generation.md)), so it moves both generators: the generic pileup one
([`locus_generation_pileup.md`](locus_generation_pileup.md)) and the STR one
([`locus_generation_ssr.md`](locus_generation_ssr.md)).

Raised by Milestone D of the generic generator (owner, 2026-07-30) after the two gaps below were
found while measuring what that generator gets right.

---

## 1. What we are building, and why now

A locus's evidence is a table of **sequence observations** — one entry per distinct thing the reads
showed there. An entry holds the bases, how much of the locus the reads spanned, which read group
they came from, and the pooled support of every read that agreed; reads matching on all three of
those key fields share one entry. Today an entry can say **"this read witnessed one contiguous run
of the locus"** and nothing else. Two things follow, and both are live today:

- **A base the read did not witness can appear in `bases`, presented as though it did.** When an
  indel's anchor position has no `Match` event, the fold takes that base from the reference so the
  indel has an anchor.
- **A read whose witness has a hole in it is dropped entirely** — no observation, just a tally in
  `reads_without_observation`. It witnessed two genuine runs; neither is recorded.

Neither is a policy choice anyone made. Both are the type being unable to say what happened, and
the fold doing the least-bad thing available. **The type's own documentation already defers this
decision to now** — `mod.rs:227-231` says of the `Observed` variant:

> Revisit when the **generic** path mints its first run: it needs runs flush with neither border (a
> read blind in the middle of a footprint), which neither constructor expresses, so the full
> constructor set — and with it the case for sealing the variant — is only knowable then. Building
> it now would be designing against one producer and guessing at the second.

The generic path now mints those runs. Both producers exist, so the constructor set is knowable.

**The type carrying an entry is renamed `ObservedSequence` → `SequenceObservation` — decided
(owner, 2026-07-30).** The value is an observation and what was observed is a sequence; the old
name said the reverse, and stopped being true when the identity widened to
`(bases, read_witness, read_group)` — from then on two entries could carry identical bases, so
"an allele from two read groups is two observed sequences" reads as a contradiction. Production
names the same shape the same way (`AlleleObservation`, `pileup_record.rs:138`), and the module's
prose has called these observations all along, including the constructor whose doc reads *"An
observation of `bases`"* while its return type says `ObservedSequence` (`mod.rs:935`). **This spec
therefore says "observation" and never "row" or "cell"** — two nicknames the old name forced into
the prose, neither ever defined, and one of them ("row") also naming two different types in the
crate and a matrix row in the aligner.

### Goals

1. **A read with a hole is an observation, not a tally** — its witnessed positions are recorded
   rather than discarded, so step 7 has something to censor. **This is the goal RNA-seq makes
   load-bearing** (§8): a spliced read whose junction falls inside a record widened across it is
   discarded whole today.
2. **`reads_without_observation` narrows to what its name says**: reads that witnessed *nothing*
   inside the footprint.
3. **The STR path's output does not move.** It mints only `Complete` and one flush-left or
   flush-right run (`ssr.rs:770`, `:821-822`, `:889`, `:989`), all of which the new
   representation must express identically.

*A fourth goal was here — "an observation says, per base, whether the read witnessed it" — and the
measurement demoted it to deferred work (§3.2, §6). It is not abandoned; it is unpaid-for.*

### Non-goals

- **Consuming partial observations.** Step 7's censored likelihood — freebayes' `1/k` prefix/suffix
  scheme — stays where `locus_generation_pileup.md` §10 put it. This spec makes the evidence
  available; what a model does with it is the model's decision.
- **Changing the no-fabrication rule.** The rule is right and is not reopened. This is about
  *saying* what it already does.
- **Per-base quality.** Observations carry quality as summed moments (`q_sum`, `mapq_sum`,
  `mapq_sum_sq`); nothing here proposes per-base arrays.
- **Sealing `ReadWitness`'s fields.** Still open, still for the same reason (§6).

### It does not

- change which loci exist, or their footprints;
- change `num_obs`' meaning — it stays *how many reads showed this*, one read counted once;
- add a second shape holding the same information (a side-channel beside the observations).

---

## 2. What the code does today — the part worth reading before designing

### The borrowed base is per **indel event**, not per observation

The residual is documented as *"one base inside an event the read genuinely witnessed"*
(`open_record.rs:1100-1111`), which is true of each occurrence and easy to read as "at most one per
observation". **It is not.** `apply_events_into` has **two** borrow sites, each guarded per event:

```rust
// open_record.rs:1217-1224 — Insertion
if offset < ref_len && offset >= consumed_until {
    allele_seq.push(ref_seq[offset as usize]);
}
// open_record.rs:1228-1238 — Deletion
let anchor_inside = *anchor_ref_pos >= record_pos;
if anchor_inside && offset < ref_len && offset >= consumed_until {
    allele_seq.push(ref_seq[offset as usize]);
}
```

`consumed_until` tracks how far emitted bases have reached, so the borrow fires exactly when **no
`Match` event emitted that position's base**. A widened footprint can hold several indel events
from one read, each with its own anchor, so **an observation can borrow one base per such event**. Any
design that carries "the anchor was borrowed" as a single flag or a single byte under-describes the
phenomenon.

The condition is also broader than adaptor masking: it is *no `Match` at that position*, which is
masking, an `N` base, or a position no `Match` covered.

**It has now been measured, and it is a corner in fact as well as by construction: 8 occurrences in
225 million event-folds, never two in one observation** (§8). That is why the design for it is
deferred (§3.2, §6) rather than built. The shape above still matters — it is what rules out the
single-flag designs, the day a count arrives that asks for one.

### The drop path

`fold_read_into_record` (`open_record.rs:1335-1345`) takes `apply_events_into`'s `None` — a
non-contiguous witness — records the read id in a per-record set, and subtracts any prior
contribution. The comment there is worth carrying into the new design: the path is reached **at
every position the record is affected at**, which is why the tally is a set of read ids rather than
a counter, and why a read that folded contiguously can *become* non-contiguous when the window
widens across an interior gap.

### What the run currently costs

An observation already heap-allocates **twice** — `bases: Box<[u8]>` (`mod.rs:147`) and
`chain_ids: Vec<ChainId>` (`mod.rs:191`) — so a mask is not a change of kind. The scale it has to
be cheap at: **1,647,161 observations over 1,541,788 loci** on chr1 of a 30× HG002 BAM, i.e.
observations ≈ loci on real data. Milestone B measured its own observation work at **+15.1 % wall /
+24.5 % allocations**, with the per-observation `bases` clone alone 2.2 % — so a **third allocation
per observation is the thing to avoid**, and an inline representation for common sizes is the
requirement, not a nicety.

*Where the observation counts come from:* Milestone D's chr1 throughput run, from the dump's own
`rows_complete` / `rows_observed` header fields (`rows_observed` was renamed `rows_partial` at
plan D4, with the variant name it was left over from). The loci count is in that run's report
(`ng_locus_generation_pileup_generator_d_2026-07-29.md`); the two observation counts are not, so
re-running the dump is the only way to reproduce them.

---

## 3. The design: a witnessed set

**An earlier draft of this section opened "the two gaps are one gap", and the measurement in §8
retired that claim.** They are two gaps. They sit on different axes — one in locus positions, one
in read coordinates — they arise from different mechanisms, and they occur at rates four orders of
magnitude apart. Bundling them made "is this worth building?" unanswerable, because the answer
differs for each. So this section settles **one** of them: a read's witness is a set of locus
positions, not one run. The other, which bases the read supplied, is deferred with its measurement
(§6).

### 3.1 `ReadWitness::Observed` carries a witnessed set, not a run

```rust
pub enum ReadWitness {
    /// The read reached both borders and witnessed every position between them.
    Complete,
    /// The positions the read witnessed, in locus coordinates — one run, or several.
    Partial { positions: WitnessedLocusPositions },
}
```

**The variant is `Partial`, not `Observed`.** Next to `Complete`, "observed" is not a contrast — a
complete witness was observed too — and once the enum itself says *witness*, the word adds nothing.
`Partial` says the one thing that distinguishes it. The name was avoided while `PartialLeft` /
`PartialRight` were recent, since those were removed for being side-tagged; the payload here is
visibly a set of positions and carries no side, so the confusion does not arise.

**`Complete` stays a variant.** It is the overwhelmingly common case — 1,646,289 of 1,647,161
observations on the chr1 run — and keeping it makes `complete_observations` (`mod.rs:123`) a cheap
equality rather than "a set equal to the whole footprint", and keeps the STR path's call sites
unmoved.

**Chosen over two alternatives, both of which were live:**

| option | why it lost |
|---|---|
| **one observation per run** — a holed read becomes two observations | breaks `num_obs` as a read count: one read would contribute 2. Deduplicating needs a read identity on the observation, and `chain_ids` cannot serve — a read agreeing with the reference carries none (`open_record.rs:473-483`), which is the common case |
| **`Observed { runs: SmallVec<[(u16, u16); 2]> }`** | keeps `num_obs` honest, but leaves `bases` as the concatenation of runs — a string the read never showed as one sequence, which is production's fill in a new costume |

A set subsumes both: one run, two runs, N runs are one shape, and `is_flush_left` /
`is_flush_right` / `positions_covered` become derivations from it, as they are already derivations
from the run (`mod.rs:329-347`).

### 3.2 What `bases` still cannot say — deferred, not solved

`SequenceObservation` gains no field here. A borrowed anchor base stays indistinguishable from a
sequenced one, exactly as today, because the measurement says the case barely exists: **8 borrowed
bases in 225 million event-folds, and never two in one observation** (§8). The design for it — a set
of indices into `bases`, on the read axis — is recorded in §6 so it can be built the day a
measurement asks for it.

Recording why it is a *separate* design and not part of the set above, since that is the thing this
spec got wrong first: the two live on different axes and the type says so today. `read_witness` is
in **locus positions**, `bases` is *"allele content, in read coordinates"* (`mod.rs:216-217`). An
insertion adds bases without positions and a deletion positions without bases, so one index cannot
address both.

### 3.3 The representation must be canonical — the trap

An observation's identity is `(bases, read_witness, read_group)` (`open_record.rs:536-540`), and
observations are merged by comparing it. **If two equal witnessed sets can have two
representations, two reads with the same witness stop sharing an observation** — silently, as extra
observations rather than as a wrong number. So:

- `WitnessedLocusPositions` is a **validated newtype with a private field**, one representation per
  set — runs sorted, non-empty, non-adjacent, non-overlapping;
- `Eq` and `Hash` are the type's own, not derived over a representation that admits duplicates;
- `open_record::witness_order` (`:261`, lifted to `pub(super)` at D1) must extend to a **total**
  order over sets, because `finalise` sorts observations with it and the output has to be deterministic
  (`locus_generation_pileup.md` §7).

This is the one place the design can go wrong quietly, which is why the newtypes are part of the
spec and not left to the arch doc.

### 3.4 Inline below a bound, heap above

A locus is at most `max_record_span` positions, and C1 capped that knob at `u16::MAX` = 65,535
(`MAX_RECORD_SPAN_CEILING`, `generator.rs:43`; the spec's own §7 still types the knob `u32`), so a
fixed-width inline mask cannot cover the worst case. It can cover every real one: C1's reasoning
for the cap is that "a locus is at most ~100 bp of reference, and a 5,000 bp record is already
unreachable with Illumina reads" (`locus_generation_pileup_generator.md`, C1). **The distribution
of real footprint lengths has not been measured**, which is why the bound is §8's question and not
a number here.

**So: inline up to a bound, heap beyond it.** The bound is a number to pick with a measurement
(§8), not by taste. `LocusLen::from_positions` **saturates** at `u16::MAX` (`mod.rs:263-265`) — the
hazard C1 capped the knob for — and the same saturation must not silently truncate a witnessed set;
whatever the encoding, the out-of-range case is an error or an assertion, not a clamp.

---

## 4. What else moves

| what | where | how it changes |
|---|---|---|
| the shared type and the field holding it | `mod.rs:145`, `mod.rs:44` | `ObservedSequence` → `SequenceObservation` (§1), and `observed_sequences` → `observations` — decided (owner, 2026-07-30). `SampleLocusObservations::observations` says the whole thing at the definition, so the field need not repeat "sequence"; `sequence_observations` was the alternative, exactly the plural of the element type and proof against drift, and it lost on length at 100 sites. Mechanical and wide: **39 uses of the type across 7 files, 102 of the field** |
| the fold's own accumulator | `open_record.rs:237` | `ObservationRow` holds the same information as the public type in the fold's layout, so it needs a name that is not "row" once the public one is an observation. **26 uses with `ObservationKey` across 4 files.** Settled in the [arch doc](../arch/locus_witness_representation.md) §3: `KeyedObservation`, with `ObservationKey` keeping its name |
| `ReadWitness` + constructors | `mod.rs:213-347` | `Observed` becomes `Partial` and its payload becomes a set; `from_left`/`from_right` keep their signatures and build single-run sets; a third constructor for an interior run, which is what the deferred note asked for — landed as `from_witnessed_runs(runs, locus_len)`, which subsumes it and is the only constructor that may answer `Complete` (plan D3, owner 2026-07-31; [arch](../arch/locus_witness_representation.md) §1.1 carries the rule). **The type is `ReadCoverage` today and is renamed here** — "coverage" reads as depth, which the type's own doc already has to correct; 193 uses of the type and 91 of the field across 12 files. Decision record in the [arch doc](../arch/locus_witness_representation.md) §3 |
| `is_flush_left` / `is_flush_right` / `positions_covered` | `mod.rs:329-347` | derivations from the set; same signatures |
| `num_obs_along_locus` | `mod.rs:69-104` | iterates the set instead of one run. **Its clamp stays** — the comment there explains why the bound is not expressible on the type, and a set does not change that |
| `complete_observations` | `mod.rs:123` | unchanged (`Complete` is still a variant) |
| the fold | `open_record.rs::apply_events_into` | returns the witnessed set instead of a `RefSpan`, and no longer returns `None` for a hole. **The two borrow sites are untouched** — that half is deferred (§3.2) |
| the drop path | `open_record.rs:1335-1345` | narrows to "witnessed nothing"; the set-of-read-ids mechanism and its reason survive |
| `witness_of` (`coverage_of` today) | `open_record.rs:184-212` | resolves a set against the final footprint rather than a span |
| the STR generator | `ssr.rs:770, 821-822, 889, 989` | call sites unchanged if the constructors keep their signatures — which is the point of keeping them |
| both dump tools | `examples/ng_ssr_loci_dump.rs`, `examples/ng_generic_loci_dump.rs` | must **print** the set, or this becomes another surface nobody reads |
| the generic census | `pileup/parity.rs` | must **count** holed witnesses and the positions inside them — the counters the §8 measurement used, kept rather than thrown away |

**The parity oracles.** The STR dump's byte-identity is the oracle for "nothing moved on the STR
path" — the same oracle that caught the `PartialLeft`/`PartialRight` reshape at plan 1's Milestone
B. For the generic path, the anchor
(`parity::ng_agrees_with_production_where_production_fabricated_nothing`) must stay green, and the
divergence census must gain a class for "ng now emits an observation where it used to count the read out"
rather than absorbing it into an existing one.

---

## 5. Cross-cutting concerns

**Memory and speed.** The whole cost is per observation, and observations ≈ loci on real data (§2).
The requirement is: **no third allocation on an observation that witnessed one run.** §8 measured
that as every witness in 225 million event-folds of DNA-seq, and two runs in the RNA-seq case, so
an encoding with two runs inline pays nothing on the common path. Anything that allocates per
observation should be rejected at review, not measured afterwards.

**Errors.** No new error type and no new fallible path. The one hazard is a set that cannot be
represented (a locus past the inline bound, or past `u16::MAX`), and that must fail loudly rather
than clamp — see §3.4.

**Determinism.** Covered by §3.3: canonical representation, total order, and the existing sort in
`finalise`. Spec §7's byte-identity-across-runs claim is unchanged and its test
(`parity::ng_emits_the_same_bytes_in_a_second_process`) still applies.

**Concurrency.** None. The walk is single-threaded and parallelism is deferred whole.

---

## 6. Deferred, with a recommended home

- **Saying which bases the read supplied — deferred on measurement, design settled.** A borrowed
  anchor base stays indistinguishable from a sequenced one. The design when it is wanted: a
  validated newtype of **ascending indices into `bases`**, on the read axis, empty in the common
  case, named for what it holds rather than for provenance (the arch doc's `UnwitnessedBases`). Not
  an `Option` — an empty set already costs no allocation once the encoding is inline, so the
  `Option` would only add a second spelling of "nothing".

  **Why it is not built:** §8 measured **8 borrowed bases in 225 million event-folds**, across
  human and tomato, and **never two in one observation**. The reachability argument stands — the
  borrow is guarded per indel event, so a widened footprint holding two indels from one read can
  borrow twice — but nothing in the data has ever asked for it. **Home: this spec.** Build it when
  a measurement, or a data type not yet tried, produces a non-trivial count.
- **The rename in the code's own doc comments.** The three sibling specs are done — 
  [`locus_generation.md`](locus_generation.md),
  [`locus_generation_pileup.md`](locus_generation_pileup.md) and
  [`locus_generation_ssr.md`](locus_generation_ssr.md) now say "observation" throughout — but the
  code still carries both old nicknames: `mod.rs:138` calls it "a table of cells" and `ssr.rs:1202`
  says "one row per (allele, read group) cell" in a single assertion message. They move with the
  code rename, not before it, since the identifiers move at the same time. **What must not move:
  the aligner's ~340 "row"s are dynamic-programming matrix rows and are correct**, as is the dump
  tools' TSV line — the STR dump's own prose about "rows" of output stays.
- **Consuming the newly-available partial evidence** — step 7's censored likelihood.
  `locus_generation_pileup.md` §10 already owns it; this spec only stops the evidence being
  discarded.
- **Sealing `ReadWitness::Observed`'s fields.** Still open for the reason `mod.rs:219-231` gives —
  a run clamped against *some* `LocusLen` proves nothing about the locus it ends up on. §3.3 makes
  the *set* a validated newtype, which is a different question from sealing the variant. Revisit
  with the arch doc.
- **Whether `reads_without_observation` survives at all.** After this change it counts reads that
  witnessed nothing in a footprint they overlapped. That is a real class (every base masked or
  `N`), and `reads_silent_over_footprint` already counts the run-level version, so the two should
  be compared before either is kept. Home: the arch doc's counts section.

  **The question sharpened at C3, and the answer is now more likely "no" (Milestone C review).**
  On the *generic* path the counter looks structurally unreachable: `process_position` skips a
  contributor whose event window is empty, `refold_live_reads` visits only reads already folded —
  whose window over a widening record can only grow — and every event surviving
  `events_overlapping` clips to a non-empty run, so `apply_events_into` cannot answer "nothing
  witnessed". The class is real and the route to it is not. **On the STR path it is the opposite:**
  the same counter carries four reasons, and C0 added the largest of them
  (`OutsideTract`, 6,704 reads on one tomato chromosome). So the two generators disagree about
  what the field is for, which is the thing to settle — not whether the number is currently zero.

---

## 7. How we know it works

1. **The STR dump is byte-identical** across the change, on the committed fixture and on a tomato
   CRAM. The STR path mints only `Complete` and one flush run, so any movement is a defect.

   **One exception, and it is a rebaseline rather than a hole in the oracle (owner, 2026-07-30).**
   The STR path also minted a flush run covering **zero** positions — a read that clips the locus
   *window* and never enters the tract, 6,704 times against 7,085 genuine partials on chr01 of
   tomato `SRR7279503`. `WitnessedLocusPositions` cannot express it, which is how the
   representation work found it. Those reads are not in the locus and the SNP/indel path owns their
   bases, so the STR path now discards and counts them
   ([`locus_generation_ssr.md`](locus_generation_ssr.md) §3). The dump moved **once**, before any
   step of this spec's own change, and only by deleting rows: 3,180 gone, every one a partial
   witness with empty bases; `obs_partial` 13,789 → 7,085, `reads_without_observation`
   2,561 → 9,265; `obs_complete`, the locus count and every non-empty observation unchanged. Every
   step from C1 on is byte-identical against **that** baseline.

   **A second deliberate move, at the end of Milestone D (owner, 2026-07-31), and this one is one
   column of one label.** The dumps gained a fourth label, `partial:both`, for a partial witness
   that touches both borders of the tract — a read that anchored one flank and whose repeat, laid
   down from that flank, covers the tract end to end without measuring the allele. It was
   previously reported as `partial:left`, *including for reads anchored on the right*. Verified
   line by line against the previous baseline: **2,530 of 8,135 observation rows differ, every one
   of them in the `read_witness` column alone, and every transition is `partial:left` →
   `partial:both`.** The two header lines, the `depth` column, the row order and every other field
   are unchanged. New baseline `tmp/witness_baseline/ssr_dump_partial_both.tsv`.
2. **The generic anchor stays green** and the census gains its new class, counted and floored — no
   divergence may be absorbed into an existing class.
3. **A fixture per new capability, each written to fail the old representation:** a read blind in
   the middle whose two runs are both recorded; a holed read that is *not* counted in
   `reads_without_observation` any more. **The spliced fixture from §8 is the third and the most
   valuable** — it is the only one drawn from a real failure rather than constructed to exercise a
   branch, and under the old representation the read vanishes from the record entirely.
4. **The observation-identity property:** two reads with the same witness share one observation.
   Mutation: perturb the set's canonical form and watch the observation count inflate — if nothing
   fails, §3.3's trap is unguarded.
5. **The memory requirement, measured, not asserted:** allocations per observation on the chr1 run,
   against Milestone D's baseline of 1,647,161 observations at 461 MB peak RSS.

---

## 8. Resolved decisions & open questions

Measured 2026-07-30 with a throwaway probe in `apply_events_into`, counting borrows and holes
directly. It changed no behaviour — the 275 locus-generation tests, including the parity anchor and
the byte-identity check, stayed green — and it was validated on the two shapes it exists to count
before any number was trusted, because zero is also what a miswired probe reports.

| run | event-folds | borrowed bases | holed witnesses | witnesses over one run |
|---|---|---|---|---|
| HG002 chr1, 30×, tandem-repeat tiers | 43,084,914 | 0 | 0 | 0 |
| tomato SRR7279503 | 89,557,864 | 8 | 0 | 0 |
| tomato SRR7279510 | 92,268,192 | 0 | 0 | 0 |

- **How often does the borrowed base fire — resolved: 8 times in 225 million event-folds, never
  twice in one observation.** Adaptor masking was live in all three runs (210 to 2,871 bases
  silenced), so this is not "the trigger never happened": it happened and almost never coincided
  with an indel anchor. **Consequence: that half is deferred** (§3.2, §6). The per-indel-event
  finding in §2 stands as a description of the code; it just has no population.
- **How much evidence does the current drop discard — resolved, and the answer differs by data
  type.** On DNA-seq, nothing: zero holed witnesses in 225 million folds. The earlier leaning
  ("an interior `N` or a ref-skip is ordinary") was wrong, and structurally so — a ref-skip is an
  RNA-seq CIGAR op, and modern Illumina puts `N`s at read ends where they cannot make a hole.

  **On RNA-seq it discards whole reads, and that is why §3.1 is being built.** Demonstrated on a
  fixture: a spliced read with a 15 bp intron, plus a read whose 20 bp deletion widens the record
  across it to positions 28–48. The spliced read witnessed **6 of those 21 positions** — three in
  each exon, two runs — and is **absent from the record entirely**; the drop path fired 6 times for
  that one read, once per position the record was affected at. Shorten the deletion and the
  footprint stops before exon 2, and the same read is recorded normally. **One base of footprint
  width decides it**, and the boundary is 17 against 18: at 17 the footprint ends at 45, one
  position short of exon 2 at 46. *(This paragraph said 16 until 2026-07-31. That number came from
  the throwaway probe; the permanent fixture — `pileup/tests.rs`,
  `one_more_deleted_base_is_what_turns_the_spliced_read_into_a_holed_one`, plan D6 — asserts both
  sides of the real boundary, and a Milestone D review found 16 to be two positions short rather
  than one.)*

  **What gates the rate, and is not known:** an intron never widens a record on its own, because a
  `Skip` emits no event (`cigar_cursor.rs:333`) and footprints grow from events
  (`open_record.rs:1790`). The hole needs an indel allele spanning the intron. Whether that is
  common depends on how the aligner treats short introns — some emit a deletion where the truth is
  a splice. **Untested: no RNA-seq alignment was available.** It does not gate the design, since
  the failure is demonstrated and the fix's shape does not change with the rate.
- **What inline bound, and what encoding — resolved: runs, two inline.** Every witness in 225
  million DNA-seq folds is one run; the RNA-seq case is two. A bitmask is out — it would pay 16
  bytes to say "positions 0..40" for a case that is one run in every observation measured.
- **Does every consumer of `positions_covered` survive the move to a set — resolved: yes, and the
  audit is small.** Grepping `positions_covered` and `offset_in_locus` across `src/` and
  `examples/` finds them only inside `src/ng/locus_generation/` and `examples/ng_generic_loci_dump.rs`.
  The paralog filter, named in an earlier draft of this question, does not read them: it is
  production code and ng's depth derivation is `num_obs_along_locus`, which iterates runs instead
  of one range and keeps its clamp. The dump's own invariant check
  (`offset_in_locus + positions_covered <= footprint`) becomes a check per run.
- **How often does the hole fire on real RNA-seq?** *No leaning: the fixture proves reachability,
  not frequency, and the aligner's junction behaviour decides it.* Settled by running the same
  probe over one spliced BAM. Worth doing before the implementation plan orders the work, and not
  worth blocking the design on.

  **The counter that answers it now exists and needs no probe (owner, 2026-07-31).**
  `reads_with_holed_witness` and `hole_positions` are on the walk's own `RunSummary` and on
  `PileupGeneratorCounts`, and `ng_generic_loci_dump` prints them in its header — so pointing that
  tool at a spliced BAM answers this question directly. They were first put only on the divergence
  census, which cannot answer it: that census is `#[cfg(test)]` and compares against production's
  walker, so it only measures loci where **production** also produced a record.

- **How much STR evidence is a partial witness that touches both edges of the tract?** *Measured
  2026-07-31, on chr01 of tomato `SRR7279503`:* **2,530 of 6,216 partial rows, 41%.** These are
  reads anchored at **one** border whose repeat, counted in read bases, reached or passed the
  reference tract's length — so the trimmed run covers the tract end to end. They are honest as
  *witnesses* (the read did see every reference position) and they are correctly **not**
  `Complete` (the read ran out, so it did not measure the allele). What collapses is the *label*:
  every one of them prints as a left-edge partial, including the reads anchored on the right. The
  open decision is whether the dumps should spell this case separately.
