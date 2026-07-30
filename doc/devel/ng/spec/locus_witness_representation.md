# ng — what a locus observation says a read witnessed

**Status:** draft, 2026-07-30. **No code yet — this settles the design.** Code-facing companion:
[`../arch/locus_witness_representation.md`](../arch/locus_witness_representation.md) (the types and
signatures); implementation plan not yet written. Changes the **shared** locus type
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

1. **An observation says, per base, whether the read witnessed it** — so a borrowed anchor base is
   distinguishable from a sequenced one.
2. **A read with a hole is an observation, not a tally** — its witnessed positions are recorded
   rather than discarded, so step 7 has something to censor.
3. **`reads_without_observation` narrows to what its name says**: reads that witnessed *nothing*
   inside the footprint.
4. **The STR path's output does not move.** It mints only `Complete` and one flush-left or
   flush-right run (`ssr.rs:770`, `:821-822`, `:889`, `:989`), all of which the new
   representation must express identically.

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

**Nobody has measured how often this fires**, on any data. It is a corner by construction and may
be a rare one; see §8.

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
`rows_complete` / `rows_observed` header fields. The loci count is in that run's report
(`ng_locus_generation_pileup_generator_d_2026-07-29.md`); the two observation counts are not, so
re-running the dump is the only way to reproduce them.

---

## 3. The design: one witnessed set, two places it is read

The two gaps are one gap. A read's witness is **a set of locus positions**, and its bases carry
**which of them the read supplied**. Today the first is forced into one run and the second is
assumed to be "all of them".

### 3.1 `ReadWitness::Observed` carries a witnessed set, not a run

```rust
pub enum ReadWitness {
    /// The read reached both borders and witnessed every position between them.
    Complete,
    /// The positions the read witnessed, in locus coordinates.
    Observed { witnessed: WitnessedPositions },
}
```

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

### 3.2 `bases` carries per-base provenance

```rust
pub struct SequenceObservation {
    pub bases: Box<[u8]>,
    /// Which bases the read supplied. `None` = all of them, the common case.
    pub base_provenance: Option<BaseProvenance>,
    // … unchanged: read_witness, read_group, num_obs, num_fwd, q_sum,
    //   mapq_sum, mapq_sum_sq, placed_left, chain_ids
}
```

`Option` so the common observation costs 8 bytes and no allocation, and so "nothing was borrowed" is said
by saying nothing. The same reasoning is already written on the generator's
`chain_ids: Option<ChainIdAllocator>` (`generator.rs:539-546`): an `Option` rather than a
placeholder, because "a placeholder starting at zero is the state this whole arrangement exists to
avoid, and it would fail silently. Absent, it fails loudly."

**Why not fold this into the coverage set.** The two live on different axes and the type says so
today: `read_witness` is in **locus positions**, `bases` is *"allele content, in read
coordinates"* (`mod.rs:216-217`). An insertion adds bases without positions and a deletion
positions without bases, so one index cannot address both. Merging them is the mistake
`locus_generation_pileup.md` §8 warns about in the other direction.

### 3.3 The representation must be canonical — the trap

An observation's identity is `(bases, read_witness, read_group)` (`open_record.rs:536-540`), and
observations are merged by comparing it. **If two equal witnessed sets can have two
representations, two reads with the same witness stop sharing an observation** — silently, as extra
observations rather than as a wrong number. So:

- `WitnessedPositions` and `BaseProvenance` are **validated newtypes with private fields**, one
  representation per set (no trailing empty word, no unsorted runs if runs are the encoding);
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
| the shared type and the field holding it | `mod.rs:145`, `mod.rs:44` | `ObservedSequence` → `SequenceObservation` (§1), and `observed_sequences` → `observations` — decided (owner, 2026-07-30). `SampleLocusObservations::observations` says the whole thing at the definition, so the field need not repeat "sequence"; `sequence_observations` was the alternative, exactly the plural of the element type and proof against drift, and it lost on length at 100 sites. Mechanical and wide: **39 uses of the type across 7 files, 100 of the field** |
| the fold's own accumulator | `open_record.rs:237` | `ObservationRow` holds the same information as the public type in the fold's layout, so it needs a name that is not "row" once the public one is an observation. **26 uses with `ObservationKey` across 4 files.** The arch doc's call; leaning `KeyedObservation` |
| `ReadWitness` + constructors | `mod.rs:213-347` | `Observed`'s payload; `from_left`/`from_right` keep their signatures and build single-run sets; a third constructor for an interior run, which is what the deferred note asked for. **The type is `ReadCoverage` today and is renamed here** — "coverage" reads as depth, which the type's own doc already has to correct; 193 uses of the type and 91 of the field across 12 files. Decision record in the [arch doc](../arch/locus_witness_representation.md) §3 |
| `is_flush_left` / `is_flush_right` / `positions_covered` | `mod.rs:329-347` | derivations from the set; same signatures |
| `num_obs_along_locus` | `mod.rs:69-104` | iterates the set instead of one run. **Its clamp stays** — the comment there explains why the bound is not expressible on the type, and a set does not change that |
| `complete_observations` | `mod.rs:123` | unchanged (`Complete` is still a variant) |
| the fold | `open_record.rs::apply_events_into` | returns the witnessed set instead of a `RefSpan`, and no longer returns `None` for a hole; the two borrow sites record provenance |
| the drop path | `open_record.rs:1335-1345` | narrows to "witnessed nothing"; the set-of-read-ids mechanism and its reason survive |
| `witness_of` (`coverage_of` today) | `open_record.rs:184-212` | resolves a set against the final footprint rather than a span |
| the STR generator | `ssr.rs:770, 821-822, 889, 989` | call sites unchanged if the constructors keep their signatures — which is the point of keeping them |
| both dump tools | `examples/ng_ssr_loci_dump.rs`, `examples/ng_generic_loci_dump.rs` | must **print** the set and the provenance, or this becomes another surface nobody reads |
| the generic census | `pileup/parity.rs` | must **count** borrowed bases and holed witnesses, for the same reason |

**The parity oracles.** The STR dump's byte-identity is the oracle for "nothing moved on the STR
path" — the same oracle that caught the `PartialLeft`/`PartialRight` reshape at plan 1's Milestone
B. For the generic path, the anchor
(`parity::ng_agrees_with_production_where_production_fabricated_nothing`) must stay green, and the
divergence census must gain a class for "ng now emits an observation where it used to count the read out"
rather than absorbing it into an existing one.

---

## 5. Cross-cutting concerns

**Memory and speed.** The whole cost is per observation, and observations ≈ loci on real data (§2).
The requirement is: **no third allocation on an observation that borrowed nothing and witnessed one
run.** `Option<BaseProvenance>` gives that for provenance; the coverage set needs an inline encoding
for the common sizes. Anything that allocates per observation should be rejected at review, not
measured afterwards.

**Errors.** No new error type and no new fallible path. The one hazard is a set that cannot be
represented (a locus past the inline bound, or past `u16::MAX`), and that must fail loudly rather
than clamp — see §3.4.

**Determinism.** Covered by §3.3: canonical representation, total order, and the existing sort in
`finalise`. Spec §7's byte-identity-across-runs claim is unchanged and its test
(`parity::ng_emits_the_same_bytes_in_a_second_process`) still applies.

**Concurrency.** None. The walk is single-threaded and parallelism is deferred whole.

---

## 6. Deferred, with a recommended home

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

---

## 7. How we know it works

1. **The STR dump is byte-identical** across the change, on the committed fixture and on a tomato
   CRAM. The STR path mints only `Complete` and one flush run, so any movement is a defect.
2. **The generic anchor stays green** and the census gains its new class, counted and floored — no
   divergence may be absorbed into an existing class.
3. **A fixture per new capability, each written to fail the old representation:** a read blind in
   the middle whose two runs are both recorded; an observation that borrowed **two** anchor bases (the case
   §2 shows is reachable and which no current test covers); a holed read that is *not* counted in
   `reads_without_observation` any more.
4. **The observation-identity property:** two reads with the same witness share one observation.
   Mutation: perturb the set's canonical form and watch the observation count inflate — if nothing
   fails, §3.3's trap is unguarded.
5. **The memory requirement, measured, not asserted:** allocations per observation on the chr1 run,
   against Milestone D's baseline of 1,647,161 observations at 461 MB peak RSS.

---

## 8. Open questions

- **How often does the borrowed base actually fire, and how many per observation?** *Leaning: rare, but
  unmeasured — and the design does not depend on the answer, only its priority does.* Settled by
  counting borrowed bases per observation over the parity soak and over the D3 real-data runs. **Do this
  before writing the arch doc**: if the answer is zero on real data, this half of the change may
  deserve a `debug_assert` and a note rather than a field.
- **How much evidence does the current drop discard?** *Leaning: more than the borrow, because an
  interior `N` or a ref-skip is ordinary.* Settled by the same three numbers D3 produced for the
  fabrication — loci, reads, and witnessed positions inside the dropped runs — which turns "ng
  forgets" into a size. **This is the number that says whether the change is worth making.**
- **What inline bound, and what encoding?** Options: a bitmask (`u128` inline covers 128
  positions), or a `SmallVec` of runs (two inline covers every case seen so far). *Leaning: runs,
  because the common case is one run and the STR path only ever produces one — a bitmask pays 16
  bytes to say "positions 0..40".* Settled by the distribution of runs per witness on real data,
  from the same instrumentation.
- **Does every consumer of `positions_covered` survive the move to a set?** Audit list:
  `num_obs_along_locus`, the STR path's flush predicates, the paralog filter's depth derivation.
  *No leaning — this is a read-the-code task for the arch doc, and if one consumer genuinely needs
  a single run, that is a finding about the consumer.*
- **Should the coverage set and the base provenance be one type?** §3.2 says no, on the axes
  argument. *Leaning: keep them separate.* Confirm when the arch doc writes the two newtypes; if
  they end up with identical shapes and one is always derivable from the other, revisit.
