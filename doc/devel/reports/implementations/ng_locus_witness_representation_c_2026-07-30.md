# ng — the locus witness representation, Milestone C (the fold)

*Implementation report, 2026-07-30. Plan:
[locus_witness_representation.md](../../ng/impl_plan/locus_witness_representation.md) Milestone C.
Design: [spec](../../ng/spec/locus_witness_representation.md) §1, §3.1, §4, §7;
[arch](../../ng/arch/locus_witness_representation.md) §2. Branch `ng-pileup-generator`,
worktree `pop_var_caller-ng-pileup`.*

**Status: C1 and C0 complete. C2–C4 not started.** This report is extended per step and
committed with each of them. C0 was added mid-milestone by owner decision — see its section,
which is at the end because that is when it was found, not where it sits in the order.

---

## C1 — `apply_events_into` hands over a set, and still discards a hole

### What the step is

The fold stops speaking in `RefSpan` — one run, inclusive of both ends — and speaks in
`WitnessedRefPositions`, the canonical set B3 built. **It still throws away a witness with a
hole in it**, which is what makes the step byte-identical by construction: every set the fold
mints holds exactly one run, so every number downstream is the number it was.

### Departures from the plan, recorded

Two, both about a signature the arch sketched before the Milestone B review changed the
constructor set underneath it.

- **`apply_events_into` returns `bool` and fills a buffer, where arch §2 sketches
  `-> Option<WitnessedRefPositions>`.** The same section's *contract* is what forced it: "the
  runs are accumulated into a buffer the caller owns and the callee clears, like
  `allele_seq`, so a fold allocates nothing per read". Returning the set means constructing
  it inside, which is `take_from` — and `take_from` moves the buffer's storage out, leaving
  the next read to allocate again. The whole point of the B review's M2 fix is that
  `refold_live_reads` refills **in place**, swapping, so the buffer inherits the old
  witness's storage. A function that returns the set cannot offer that to its caller. So the
  buffer is a parameter, the `bool` says whether anything was witnessed, and the caller picks
  `take_from` or `refill_from`. The arch's signature predates the constructors it would have
  to use.

  The contract is stated on the function in both directions, because only one half is new:
  `true` leaves at least one non-empty run in the buffer; **`false` leaves it empty**, so a
  read that got no observation cannot leak runs into the next read's witness.

- **`FoldedReadState` is `Clone`, no longer `Copy`.** A set is not `Copy`.
  `refold_live_reads` used to lift the whole state out with a `*` and rebuild it; it now
  copies the four scalar facts out of a `&`-pattern and refills the witness through the
  `&mut`. **The exhaustive-destructure discipline the old code's comment argues for is kept,
  and now applies at both ends** — the read side names all five fields, and the write side
  destructures the destination with no `..`, so a field added to `FoldedReadState` still
  cannot be silently carried stale.

`fold_read_into_record` uses `take_from` rather than `refill_from`, as the plan says. Worth
recording that this function *also* runs on a re-fold at a later walker position, where the
`remove` drops the read's previous witness and its storage with it — which costs nothing
while a witness is one or two runs, since those are inline, and is noted in the code as the
place to copy `refold_live_reads`' shape if a multi-junction witness ever makes it allocate.

### Two conventions became one

`RefSpan` was inclusive; the set is half-open. Three `± 1`s existed to bridge the two and all
three are gone in the same commit:

| site | was | is |
|---|---|---|
| `apply_events_into`'s return | `end: end - 1` | the half-open run, pushed as-is |
| `witness_of`'s `past_last` | `witnessed.end.saturating_add(1)` | `last_end` |
| `witness_of`'s intersection assert | `witnessed.end >= record_pos` | `last_end > record_pos` |
| `read_agreed_with_reference` | `(end - pos).saturating_add(1)` | `end - pos` |

The last two are the traps the B review named ([review
M3](../reviews/ng_locus_witness_representation_b_2026-07-30.md)), and both compile either
way. What catches a mistake in them is that **the derived expectations in the existing tests
were not edited** — `positions_covered: 3`, `offset_in_locus: 13`, `ReadWitness::Complete`,
the chain-id assertions — only the literal extents, which gained a `+ 1` on their end.

`witness_of` and `read_agreed_with_reference` both still read the witness's **enclosing
extent** rather than walking its runs, which is exact while a set holds one run and is
labelled in both places as what C2 and C3 respectively have to change. Clamping the enclosing
extent is precisely how the hole gets swallowed, so the note says so where the code is.

### The one new test, and the mutation it fails under

`apply_events_into_clears_the_witness_buffer_it_was_handed`. C1 is what introduces a
caller-owned runs buffer, shared across every read of every record, so "the callee clears it"
becomes a contract that can be forgotten. Deleting `witnessed_runs.clear()`:

```
assertion `left == right` failed: the previous read's run must not survive into this one's
witness
  left: [(100, 105), (102, 103)]
 right: [(102, 103)]
test result: FAILED. 292 passed; 1 failed
```

The leak is silent twice: nothing panics, and the witness is a third of an observation's
identity, so the read either stops sharing an observation with the reads it agrees with or
shares one claiming positions neither saw. The test asserts both halves of the contract.

The `#[cfg_attr(not(test), expect(dead_code, …))]` came off `take_from`, `refill_from` and
`runs`, which C1 wires. It stays on `from_half_open_runs` and `positions_covered` with
corrected reasons: the fold owns a buffer rather than building sets from literals, and it
measures against the *footprint* rather than against the runs' own total. Both attributes are
`expect` rather than `allow`, so wiring either later fails the build.

### Validation

| check | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --examples --all-features` | **2,824** passed, 0 failed (2,823 before; +1 is the new test) |
| `ng::locus_generation` | **293** passed (292 before) |
| STR dump, tomato `SRR7279503` chr01, 11,318 lines | **byte-identical** to `ssr_dump_a2.tsv` |
| `parity::ng_agrees_with_production_where_production_fabricated_nothing` | green |
| `parity::ng_emits_the_same_bytes_in_a_second_process` | green |
| `parity::every_divergence_from_production_is_one_of_the_six_named_classes` | green |

No test expectation changed except the literal extents described above; no census counter
moved.

---

## C0 — a partial witness must have witnessed something (added mid-milestone)

### How it was found

C2 swaps `ReadWitness::Partial`'s payload for a `WitnessedLocusPositions`, and that type has no
empty representation. Writing it turned up a shape the spec, the arch and the plan had all
assumed away: the STR path mints `Partial { offset_in_locus, positions_covered: 0 }`. Worse than
unrepresentable — `Partial { 0, 0 }` and `Partial { len, 0 }` are *different values* today
(flush-left against flush-right, different sort key, different bucket), and a set of positions
cannot tell them apart, because both are the empty set.

### What those reads are, measured before deciding anything

A probe on the mint ([ssr.rs:816](../../../../src/ng/locus_generation/ssr.rs#L816)), over chr01 of
tomato `SRR7279503` — 13,789 partial mints:

| | `reach == 0` | `reach > 0` |
|---|---|---|
| count | **6,704** | 7,085 |
| read's extracted region over the locus window | min 1, **median 16**, max 39 bases | min 2, **median 35**, max 128 |
| region shorter than the 30-base flank | 6,440 of 6,704 (96 %) | 2,184 of 7,085 (31 %) |

and their offsets are degenerate in exactly one way each — 3,446 `border=Left` with
`tract_start == tract_end == region_len`, 3,258 `border=Right` with `tract_start == tract_end == 0`.

**They are reads that overlap the locus window by less than the flank and never reach the tract.**
A read whose slice over the window is 16 bases, against a 30-base left flank, cannot have entered
the repeat; the aligner crosses the left junction at the region's own end, `anchor_firm`'s
`if !right_anchored { tract_end = read_len }`
([ssr_anchor_firm.rs:710-715](../../../../src/ng/alignment/ssr_anchor_firm.rs#L710)) puts the far end
at the same offset, and the span comes out empty. The `Right` cases are the mirror image, entirely
inside the right flank. The lower bound they supply is "the tract is at least **0** long".

**A first framing of this was wrong and is recorded so it is not repeated.** It read the count as
"half the reads are anchored on the flank but cover zero tract bases", which invited a decision
about *representation*. The owner rejected the premise, the geometry above is what checking it
produced, and the real question was never about the witness type: it is whether the STR generator
should mint these observations at all.

### The decision

**Owner, 2026-07-30: they are not in the locus.** The locus is the tract; these reads are in the
flank, and the SNP/indel path analyses those bases. The STR path discards them and counts them.

`NoObservationReason::OutsideTract` and `SsrGeneratorCounts::outside_tract`, with the guard in
`classify_read`: a `FromLeft`/`FromRight` span that is empty is no observation. `Between` with an
empty tract is untouched — a fully deleted tract is a real measurement. Recorded in
[`locus_generation_ssr.md`](../../ng/spec/locus_generation_ssr.md) §3, which is the design's home;
the dump's `reads_without_observation` line now sums all four reasons, which is why the third
smallest of them no longer stands for the whole.

### The test, and the mutation it fails under

`a_read_covering_only_a_flank_is_outside_the_tract`, asserting both borders because they arrive by
mirror-image routes and a fix catching one would halve the population and look like it worked.
Deleting the guard:

```
assertion `left == right` failed: a read that stops at the tract's first base witnessed no
tract position
  left: Observed { bases: [], read_witness: Partial { offset_in_locus: 0, positions_covered: 0 },
                   q_sum: -0.0 }
 right: NoObservation(OutsideTract)
```

which is the defect in its own words.

### The oracle moved, once, and only by deletion

This is the one step of the plan that changes the STR dump. On chr01 of tomato `SRR7279503`:

| | before | after |
|---|---|---|
| dump rows | 11,318 | 8,138 |
| `obs_partial` | 13,789 | 7,085 |
| `reads_without_observation` | 2,561 | 9,265 |
| `obs_complete` | 15,404 | 15,404 |
| `ssr_loci` / `zero_coverage` | 213,344 / 211,277 | 213,344 / 211,277 |

The arithmetic closes exactly: −6,704 partial observations, +6,704 reads without one, and
**3,180 removed rows of which 3,180 are a partial witness with empty bases** — no complete
observation, no locus and no non-empty observation moved. New baseline:
`tmp/witness_baseline/ssr_dump_outside_tract.tsv`; every step from C1 on is byte-identical against
it.

### Validation

| check | result |
|---|---|
| `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --examples --all-features` | **2,825** passed, 0 failed |
| `ng::locus_generation` | **294** passed |
| STR dump | moved as tabulated above, and only so |

### What is left open

The fetch queries the tract **plus its margin**, so the delimiter is still handed regions that
cannot span the tract — 6,440 of the 6,704 were shorter than one flank. Classifying them at the
output is correct and cheap; not handing them to the aligner at all would be cheaper. Not done
here, because it is a change to the read fetch and this milestone is about the witness. **Home:**
[`locus_generation_ssr.md`](../../ng/spec/locus_generation_ssr.md) §2, when someone measures what
the aligner spends on them.

---

## C2 — `witness_of` resolves a set, and `Partial` carries one

### Scope: C2 absorbed C4, D1 and D2 (owner, 2026-07-30)

C2's own contract — "clamp **each run** into the footprint, never the enclosing span" — only
means anything if `ReadWitness::Partial` carries a `WitnessedLocusPositions`. That swap does
not compile alone. In one commit it forces:

| plan step | what the swap forces | why it cannot wait |
|---|---|---|
| **C4** | `sort_key` borrows: `(u8, &[(u16, u16)])` | `(u8, u16, u16)` cannot be a total order over a set, and `ReadWitness` stops being `Copy` |
| **D1** | `num_obs_along_locus` iterates runs | it read the two `u16` fields directly |
| **D2** | `is_flush_left` / `is_flush_right` delegate | one line each; `WitnessedLocusPositions` already implements both |
| part of **D4** | the generic dump's and the census's per-`Partial` invariants | they destructure the fields |

13 files carry `ReadWitness`, ~200 sites. Raised before writing any of it; the owner chose to
absorb rather than split. **C3 is untouched by this** — it still lands alone and bisectable.
D4 keeps its real content (the three dumps' label drift, sharing `witness_label`'s
derivation); D1 and D2 are marked landed-at-C2 in the plan.

### The constructors return `Option`, which arch §1.1 did not ask for

`from_left` / `from_right` are `-> Option<Self>` where the arch keeps them infallible and §5
says "call sites unchanged". C0 is why: a clamped run covering **no** position cannot be a
`WitnessedLocusPositions`, and the three alternatives were all worse — panicking on a public
constructor, fabricating a position, or leaving the empty case expressible on a type whose
producer had just been fixed for minting 6,704 of them.

The result is better than the signature it replaced: the STR path's only mint now reads

```rust
let Some(read_witness) = (match border { … }) else {
    return Classified::NoObservation(NoObservationReason::OutsideTract);
};
```

so C0's guard and the type give the **same** answer instead of being two decisions that can
drift. Two production call sites and ~20 test sites moved; the plan's D3 note ("`from_left` /
`from_right` keep their signatures") is superseded.

**One arch contract deliberately not implemented.** §1.1 says all constructors "return
`Complete` when the result covers the whole locus". They do not, and must not here: on the STR
path an expanded allele gives a reach past the reference tract, so `from_left(n ≥ len, len)`
is reachable and would flip observations from `Partial` to `Complete` — moving the STR output
that spec §1 goal 3 protects. The existing full-locus merge (`from_left(len, len)` and
`from_right(len, len)` are the same value) is preserved: both build the same one-run set.
**Raise at D3**, which owns the constructor set.

### What `witness_of` does now

Each run is intersected with the footprint and rebased onto the locus; runs falling wholly
outside are dropped rather than becoming runs of no positions; the survivors go through
`WitnessedLocusPositions::from_half_open_runs`. **`Complete` is decided on the set's total
coverage**, not on its outermost edges — a witness flush at *both* borders can still have a
hole, and calling that complete is the fabrication C2 exists to prevent.

The empty-after-clamp case is an `expect`, not an `Option` return (the plan flagged the
choice). The argument is the fold's: every run is clipped into `[record_pos, record_end)` when
it folds, a record's anchor never moves, and a widen only extends the right — so no run can
lie outside the *final* footprint. Returning `Option` would have made every caller invent a
policy for a state the fold cannot produce.

### Four new tests, each failing under the defect it names

C2 moves no byte — the fold still discards a holed witness until C3 — so its change is only
visible on sets built directly. All four mutations were run and the output quoted:

| test | mutation | result |
|---|---|---|
| `witness_of_a_witness_with_a_hole_is_not_complete_however_far_its_ends_reach` | `witness_of` clamps the enclosing extent | **3 failed**, incl. this one |
| `witness_of_clamps_each_run_rather_than_the_extent_enclosing_them` | *(same)* | fails |
| `witness_of_drops_a_run_that_falls_outside_the_footprint_entirely` | *(same)* | fails |
| `depth_over_a_witness_with_a_hole_leaves_the_hole_at_zero` | `num_obs_along_locus` sums the enclosing extent | **1 failed** |
| `witness_order_is_total_over_witnesses_of_several_runs` | `sort_key` returns only the first run | **1 failed** |
| `a_constructor_asked_for_no_positions_answers_none` | `from_left` saturates a zero-length run to one position | **1 failed** |

The first is the milestone's point in one fixture: a read touching both borders of a 5..21
footprint and blind between them is credited with all sixteen positions under the enclosing
reading, and with the six it saw under the set.

### One test rewritten rather than kept

`a_run_whose_end_would_overflow_is_still_flush_right` named a hazard C2 removed: `is_flush_right`
had to *add* offset and length, and the Milestone B review showed a `wrapping_add` there
survived the whole suite. A run now carries its end, so there is no addition; and the
offset-and-length spelling that could overflow lives in one constructor that rejects the sum,
pinned by `an_empty_input_or_an_empty_run_is_rejected_rather_than_dropped`. What remains
testable — the predicate at the last representable position — is
`a_run_reaching_the_last_representable_position_is_flush_right`, which fails if the predicate
reads the wrong end of the set.

### Validation

| check | result |
|---|---|
| `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --examples --all-features` | **2,831** passed, 0 failed (2,825 before) |
| `ng::locus_generation` | **300** passed (294 before) |
| STR dump vs the C0 baseline | **byte-identical** |
| the three parity anchors | green |

The generic dump's `read_witness` column renders one `<offset>+<positions>` per run, comma
separated — identical for the one-run witnesses the generic path mints today. **D4 still owns
what that column finally says**, and the label drift across the three STR dumps with it.
