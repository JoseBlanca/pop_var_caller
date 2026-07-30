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
