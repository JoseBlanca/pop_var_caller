# Milestone C — behaviour-change safety review

Scope: did each of `ebe3685` (C1), `6805e42` (C0), `761d53e` (C2), `82b13a0` (C3),
`fc7d839`, `3fecdf6` change exactly what it claimed, and nothing else?

Method: a differential harness rather than the edited tests. A `review_differential_render`
test was appended to `examples/ng_generic_loci_dump.rs`'s test module at each commit,
rendering eleven fixtures (plain coverage, the masking fixture masked / unmasked / chunked,
an interior `N`, a ref-skip, the accounting fixture, two read groups, an insertion, a
widen-across-a-hole, masked+holed) through the tool's own `report.render()`. 507 lines of
TSV per commit, diffed byte for byte. Plus mutation runs, and hand re-derivation of every
flipped expectation from the fixture geometry.

Baseline reproduced: `cargo test --lib --all-features ng::locus_generation` →
`303 passed; 0 failed; 1 ignored`. Full lib suite: `2756 passed; 0 failed; 5 ignored`.

---

## Claims verified as sound

### 1. C1 is byte-identical on the generic path (VERIFIED)

```
$ diff diff_preC1.txt diff_C1.txt && echo "IDENTICAL: pre-C1 vs C1"
IDENTICAL: pre-C1 vs C1
```

Independent of the edited test expectations. The three `± 1`s are also algebraically
equivalent, checked by hand:

- `apply_events_into`'s `end - 1`: the internal `witnessed` tuple was already half-open
  (`(event_start, event_end)`), so dropping the `- 1` and pushing the pair is the identity.
- `witness_of`: old `witnessed.end.saturating_add(1)` on an inclusive end ≡ new `last_end`
  on a half-open one, for every value below `u32::MAX`.
- The assertion flip: `end_incl >= record_pos` ⟺ `last_end - 1 >= record_pos` ⟺
  `last_end > record_pos`. Exact, and a `debug_assert` either way.
- `read_agreed_with_reference`: old
  `(end_incl.saturating_sub(pos)).saturating_add(1)` vs new `end_excl.saturating_sub(pos)`.
  These agree whenever `end_excl > pos` and **differ** when `end_excl <= pos` (old → 1,
  new → 0). Unreachable: `apply_events_into` only records a run when
  `event_start < event_end` with `event_start = max(anchor, record_pos)`, so every stored
  run has `end > start >= record_pos`; and a record's anchor never moves (`grep -n "\.pos = "
  src/ng/locus_generation/pileup/open_record.rs` → no hits; `widen` only extends right).

`canonicalise_runs` running over the fold's single run is a no-op (`sort_unstable` on one
element, merge loop `1..1` empty, `truncate(1)`).

The `FoldedReadState` `Copy` → `Clone` change is behaviour-preserving: the fields the
re-place used to copy out and write back (`contribution`, `chain_id`, `read_group`) are now
simply left in place, and nothing between the read and the write touches that map entry —
`subtract_contribution` / `add_contribution` mutate `alleles`, not `folded_reads`.

### 2. C2 moves no byte (VERIFIED)

```
$ diff diff_C1.txt diff_C2.txt && echo "IDENTICAL: C1 vs C2"
IDENTICAL: C1 vs C2
```

Every derived quantity checked for one-run equivalence, not just the dump:

| site | pre-C2 | post-C2 | equal over one run? |
|---|---|---|---|
| `mod.rs:74` `num_obs_along_locus` | `from = min(off,len)`, `to = min(from+cov, len)` | `from = min(start,len)`, `to = min(end,len).max(from)` | yes — `.max(from)` is a no-op when `off <= len`, and both give an empty range when `off > len` |
| `witness.rs:341` `is_flush_right` | `off.saturating_add(cov) >= locus_len` | `last.end >= locus_len` | yes — the constructor rejects the overflow the `saturating_add` existed for |
| `witness.rs:329` `is_flush_left` | `off == 0` | `first.start == 0` | yes |
| `parity.rs` `fabricated_ref_bases` | `footprint - positions_covered` | `footprint - positions.positions_covered()` | yes |
| `parity.rs` `ObservationIdentityWithoutGroup` | `(&[u8], (u8,u16,u16))` | `(&[u8], (u8,&[(u16,u16)]))` | yes — injective relabelling, same dedup |
| `ssr.rs:1141` bucket key `HashMap<(bases, ReadWitness, group)>` | `(off,cov)` | `[(off, off+cov)]` | yes — bijection at fixed `off` |

### 3. The `sort_key` order claim (VERIFIED — proven, not assumed)

Old `(1, offset, covered)` vs new `(1, &[(start, end)])` with `start = offset`,
`end = offset + covered`.

*Proof.* The tag decides first in both. At equal tags, the first component is `offset` in
both. At equal offsets, `end1 < end2 ⟺ offset+cov1 < offset+cov2 ⟺ cov1 < cov2`. Slice
length is 1 in both, so no length tiebreak arises. The map is therefore order-isomorphic —
including the two cases the review asked about: **two runs sharing a start** (decided by
`end` ≡ decided by `covered`, same direction) and **a run whose end orders differently from
its length** (only possible at *different* starts, where the first component already
decided identically in both keys).

Confirmed exhaustively over a 132 × 132 grid (offsets 0..12 × lengths 1..12, both
directions) by a temporary probe:

```
test ng::locus_generation::witness::tests::review_sort_key_order_is_unchanged_over_one_run_witnesses ... ok
test result: ok. 304 passed; 0 failed; 1 ignored
```

`sort_key` is also injective (`Complete` → tag 0 with an empty slice; every `Partial` has a
non-empty canonical slice), so the borrow does not make any sort order-dependent. All four
call sites use it inside a `.cmp()` in a `sort_by` comparator or a `HashSet` key
(`ssr.rs:1205`, `parity.rs:1364`, `parity.rs:2310`, `open_record.rs:753`) — no
`sort_by_key`, so the lifetime change introduces no behavioural difference.

### 4. C3 is the only commit that moves generic-path output, and it only *adds* rows (VERIFIED)

```
$ diff diff_C2.txt diff_C3.txt
178,179c178,179
< # locus_sum_reads_without_observation=1 locus_sum_reads_discarded_by_cap=0
< # rows_complete=31 rows_observed=0 reads_complete=84 reads_observed=0
---
> # locus_sum_reads_without_observation=0 locus_sum_reads_discarded_by_cap=0
> # rows_complete=31 rows_observed=1 reads_complete=84 reads_observed=1
192a193
> chr1	20	24	CGTTA	2	observed:0+2,3+2	0	CGTA	1	2
216,217c217,218
...
230a232
> chr1	20	24	CGTTA	2	observed:0+2,4+1	0	CGA	1	2
...
267a270
> chr1	20	24	CGTTA	2	observed:0+2,3+2	0	CGTA	1	2
...
433a437
> chr1	20	24	CGTTA	1	observed:0+2,3+2	0	CGTA	1	1
```

Four of eleven sections changed. Every change is either a **new** row or a header counter
following from it. **No existing row was modified or deleted, and no `chain_ids` field on
any surviving row moved.** That is a stronger statement than the commit makes.

### 5. `fc7d839` and `3fecdf6` change no generic-path byte (VERIFIED)

```
$ diff diff_C3.txt diff_3fecdf6.txt; echo "exit=$?"
exit=0
```

`fc7d839` is comments only, mechanically confirmed:

```
$ git show fc7d839 -- examples/ src/ | grep -E "^[-+]" | grep -vE "^[-+][-+][-+]" \
    | grep -vE "^[-+]\s*(///|//|$)"
(no output)
```

`3fecdf6`'s non-comment changes are exactly: the bakeoff's four-reason sum, the
`classify_read` guard removal, moving `*allele_index = new_index;` after the refill, and
four added tests plus one strengthened assertion. Nothing else.

### 6. The guard removal in `3fecdf6` is equivalent on every reachable input (VERIFIED)

`ssr.rs:752` (was the `RepeatSpan::FromLeft(t) | FromRight(t) if t.is_empty()` arm) →
`ssr.rs:865-870` (`partial()`'s `let Some(read_witness) = … else { return
NoObservation(OutsideTract) }`).

With `tract.start == tract.end`:
- `reach = 0`;
- `from_left(0, len)` → `covered = 0` → `one_run_from_offset_and_length(0, 0)` →
  `from_half_open_runs([(0,0)])` → `canonicalise_runs` rejects `start >= end` → `None`;
- `from_right(0, len)` → `covered = 0` → `(len, len)` → same rejection → `None`;
- both return `Classified::NoObservation(NoObservationReason::OutsideTract)`.

**Same reason, same counters.** `tally` (`ssr.rs:1169-1176`) increments the per-locus
`reads_without_observation` and `counts.outside_tract` from the enum value alone, so the two
routes are indistinguishable downstream. **No side effect in between**: `partial` takes no
`qual_buffer`, borrows two slices without mutating, and never reaches `ln_p_err_sum`.

`locus_len == 0` is not a hazard: `SsrSegment::tract_len()` is `end - start + 1`
(`segment_criteria.rs:223`), so it is at least 1. And a *non-empty* tract can never be
turned into `None` by the constructors: `covered = min(reach, locus_len) >= 1`. So C2's
`Option` fires **exactly** where C0's guard fired, and nowhere else.

(One residual difference is recorded as F2 below.)

### 7. C3's seven flipped tests are the intended consequence, and five are stronger (VERIFIED)

Every expectation re-derived from the fixture geometry, independently of the test:

- **`apply_events_a_hole_in_the_middle_is_recorded_as_two_runs`** — `ref_seq = b"ACGTA"` at
  100, events `Match@100`, `Match@103`. Run opens `(100,101)`; `103 > 101` so a hole is
  pushed and `(103,104)` opens; the tail push closes it. Runs `[(100,101),(103,104)]`,
  `positions_covered = 2`. Bases: offsets 0 and 3 are both `< ref_len`, so `b"AT"`. ✓
- **`a_read_blind_inside_a_footprint_is_recorded_with_a_hole_in_its_witness`**, the hardest
  pair:
  - *interior `N`*: `blind` is `M30` from 11 with position 22 = `N`, so the cursor emits no
    event there. Over the record `20..=24` (`record_pos = 20`, `end = 25`) the events are
    20, 21, 23, 24. `(20,21)` merges 21 → `(20,22)`; `23 > 22` → push, open `(23,24)`,
    merge 24 → `(23,25)`. Rebased: `(0,2),(3,5)` → label `observed:0+2,3+2`, and
    `positions_covered = 4 ≠ 5` so `Partial`. Bases = ref[20],ref[21],ref[23],ref[24] =
    `ref_span(20,21) ++ ref_span(23,24)`. **Both expectations reproduce exactly.** ✓
  - *ref-skip*: `M11` (11..21), `Skip 2` (22,23), `M19` (24..). Events at 20, 21, 24.
    `(20,22)`, then `24 > 22` → `(24,25)`. Rebased `(0,2),(4,5)` → `observed:0+2,4+1`;
    bases `ref_span(20,21) ++ ref_span(24,24)`. ✓
- **`a_read_whose_witness_splits_when_the_record_widens_stays_in_it`** — `wide` is `M25`
  from 1 with reference position 11 = `N`; the record anchors at 5 and widens to `5..=20`
  (16 positions, `end_excl = 21`). Events 5..10 then 12..20. Runs `(5,11)` and `(12,21)`;
  rebased at 5 → **`[(0,6),(7,16)]`**, `positions_covered = 6 + 9 = 15` of 16. ✓ Both
  literal expectations reproduce. `folded == 3` is opener + widener + wide, each `num_obs`
  1.
- **`a_read_with_a_hole_is_counted_neither_as_capped_nor_as_witnessing_nothing`** —
  `holey` = `ACNTA` over `5..=9`; runs `(5,7)`,`(8,10)` → offsets `(0,2),(3,5)`, 4 of 5
  positions. `reads_partial == 1` is a **new** assertion; the pre-C3 test had none.
- **`a_read_folding_at_four_positions_of_one_record_is_one_observation`** — gained
  `reads_complete == 1`, `reads_partial == 1` and `num_obs == 1` from four folds, replacing
  a single `reads_without_observation == 1`. Strictly stronger.
- **`a_read_whose_witness_splits_at_a_widen_keeps_its_observation`** — gained
  `holed_state.witnessed.runs().len() == 2` and `num_obs sum == 2`.
- **`every_read_the_walk_saw_is_accounted_for`** — `> 0` → `== 0` is a tightening, and the
  differential's `accounting` section confirms it end to end.

**None of the seven was weakened to pass.** Five gained assertions.

### 8. C0's arithmetic (VERIFIED from the code path)

"Every removed row is a partial witness with empty bases" is structural, not observational.
The removed arm fires only when `tract.is_empty()`, and pre-C0 `partial()` built
`bases: region_seq[tract].into()` — an empty slice by definition of the guard — with
`reach = 0`, hence `Partial { offset_in_locus: 0 | locus_len, positions_covered: 0 }` and
`q_sum = ln_p_err_sum(&[]) = -0.0`. That is the exact value the commit quotes.

**Could a legitimate observation have been dropped?** No.
- The arm names `FromLeft | FromRight` only. `RepeatSpan::Between` with an empty tract — a
  fully deleted tract, a real measurement — still reaches `complete_or_low_quality`
  (`ssr.rs:791`) and still yields `ReadWitness::Complete` with empty bases. Untouched.
- A dropped row contributed nothing to `num_obs_along_locus` (an empty range), nothing to
  allele evidence (empty bases) and `-0.0` to `q_sum`. Only a phantom count is lost.
- The `zero_coverage` figure being unchanged is right for the right reason: the STR dump's
  test (`ng_ssr_loci_dump.rs:96-101`) requires `observations.is_empty() &&
  reads_without_observation == 0`, and the reclassified reads land in the second term, so a
  locus whose only reads were outside the tract does not become a zero-coverage locus.
- `reads_without_observation` 2,561 → 9,265 is exactly +6,704, matching `obs_partial`
  13,789 → 7,085. The 3,180 rows are the *distinct observations* those 6,704 reads merged
  into (at most two per locus: one left-flush and one right-flush key, same empty bases).
  Internally consistent.

Both live consumers of the reason set now sum all four
(`ng_ssr_loci_dump.rs:254-257`, `ng_ssr_aligner_bakeoff.rs:376-379`); a repo-wide grep for
`no_border_anchored` finds no third.

### 9. C1/C3 cannot reach the STR path at all (VERIFIED)

```
$ grep -rn "apply_events_into\|apply_events(" src/ | grep -v open_record.rs
(only doc-comment references in witnessed_ref.rs)
```

The fold is confined to `open_record.rs`. The STR generator mints witnesses only through
`from_left` / `from_right`, which are one run by construction, so C3's hole branch is
structurally unreachable from it.

### 10. `3fecdf6`'s M1 fix is genuinely pinned (VERIFIED by mutation)

Restoring the pre-C3 enclosing-extent body of `read_agreed_with_reference`:

```
thread '…::a_witness_with_a_hole_never_counts_as_agreeing_with_the_reference' panicked at
src/ng/locus_generation/pileup/open_record.rs:3260:9:
the read said nothing about 101, so it cannot have agreed with the reference there — and
`ACG` equalling the enclosing slice is not evidence that it did
test result: FAILED. 302 passed; 1 failed; 1 ignored
```

---

## Findings

### F1 — Minor. `finalise`'s cap exclusion is now a guard that cannot fail, and nothing says so

`src/ng/locus_generation/pileup/open_record.rs:683-690`

```rust
reads_discarded_by_cap: self
    .reads_discarded_by_cap
    .iter()
    .filter(|read_id| {
        !self.folded_reads.contains_key(read_id)
            && !self.reads_without_observation.contains(read_id)   // ← this clause
    })
    .count() as u32,
```

C3 removed the *cause* the second clause defended against, so no input can now distinguish
the two-clause filter from the one-clause one. Deleting the clause:

```
test result: ok. 303 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out
```

This is exactly the class the `3fecdf6` sweep was built to find — it found five — and this
one survived it, even though the commit's own comment at line 677-682 identifies the state
as unreachable ("the `!contains` below is the statement of that rather than a live filter").
The comment argues the clause is worth keeping; the review's own standard is that a kept
guard needs a test that fails when it is deleted.

**Fix.** Either (a) add a unit test that reaches `finalise` with a read id in *both*
`reads_discarded_by_cap` and `reads_without_observation` — the record can be built directly,
as `a_witness_with_a_hole_never_counts_as_agreeing_with_the_reference` does — or (b) if the
state is to be treated as genuinely unreachable, replace the clause with a `debug_assert!`
that the two sets are disjoint, which states the invariant instead of silently absorbing its
violation.

### F2 — Nit. The removed `tract.is_empty()` guard also covered an inverted span; `partial()` does not

`src/ng/locus_generation/ssr.rs:752` (the removed arm) and `ssr.rs:850`.

`Range::is_empty()` is `start >= end`, so the C0 guard swallowed `start > end` as well as
`start == end`. `partial()` does not:

```rust
let reach = (tract.end - tract.start).min(u16::MAX as usize) as u16;   // ssr.rs:850
```

On an inverted range this underflows — a panic in debug; in release it wraps to a huge
`usize`, clamps to `u16::MAX`, and then `region_seq[tract]` (`ssr.rs:872`) panics with a
slice-index message that names neither the cause nor the read. Where the guard answered
`OutsideTract`.

It is unreachable in practice: all three delimiters funnel through
`TractReadout::classify` (`ssr_best_path_flat_gap.rs:335`), and `ssr_anchor_firm.rs:707-712`
sets the unanchored side to `0` or `read_len`, so `FromLeft` always has `end = read_len` and
`FromRight` always `start = 0`. But `classify`'s own ordering check is a **`debug_assert`**
(line 340) on a struct whose fields are public precisely so the B3 parity harness can build
one raw — and this repo's own comment (`mod.rs:87-90`) records that a debug-only guard
compiles out of the build it actually runs.

**Fix.** Add `debug_assert!(tract.start <= tract.end)` at the top of `partial()`, or note in
its doc block that the inverted case is answered by a panic rather than by `OutsideTract`,
so the commit's "changes no byte" claim is scoped to the reachable input set it is true of.

### F3 — Nit. A wrong number in an assertion message added by C3

`src/ng/locus_generation/pileup/open_record.rs:2895`

```
"and it is not counted out either — the tally means 'witnessed nothing' now, \
 and `wide` witnessed nineteen of the record's sixteen... positions on both \
 sides of its N. Record: {widened:?}"
```

`wide` witnessed **fifteen** of the record's sixteen positions — which the very next
assertion in the same test states (`positions_covered() == 15`). The message would mislead
whoever reads it on a failure.

**Fix.** "fifteen of the record's sixteen positions".

### F4 — Nit. The accounting fixture's headline identity cannot fail

`examples/ng_generic_loci_dump.rs:1028-1035`

```rust
assert_eq!(
    report.reads_complete + report.reads_observed,
    report.rows.iter().map(|row| u64::from(row.reads)).sum::<u64>(),
    "the per-class read totals are the rows' reads, split by coverage"
);
```

`push_locus` (`ng_generic_loci_dump.rs:181-221`) increments `reads_complete` / `reads_observed`
by `obs.num_obs` and then pushes a row with `reads: obs.num_obs`, in the same loop iteration.
The two sides are the same additions, so the assertion is a tautology over any input.

Pre-existing rather than introduced by C3 — but C3 edited this test and its message calls it
"the identity still holds", which reads as evidence it is not.

**Fix.** Either drop it, or make it an identity that can fail — e.g. assert against
`reads_admitted` minus the counted-out classes, which is the statement the surrounding
comment actually makes.

### F5 — Nit (confirmation, not a defect). `reads_without_observation` on the generic path

Confirmed independently: the differential shows it at 0 in every one of the eleven fixtures
after C3, including the two built specifically to populate it. `a_read_folding_at_four_
positions_of_one_record_is_one_observation`'s own doc concedes its original property is now
"vacuous by construction". `3fecdf6` already records this ("structurally unreachable … the
class is real and the route to it is not") and spec §6 owns the question, so this is a
confirmation of the author's own note rather than a new finding.

---

## Bisectability

`git bisect` over the generic dump gives a clean single culprit:

| commit | generic dump render |
|---|---|
| `ebe3685^` (pre-C1) | baseline |
| `ebe3685` (C1) | **identical** |
| `761d53e` (C2) | **identical** |
| `82b13a0` (C3) | 4 sections gain one holed row each; 4 counters move |
| `3fecdf6` (HEAD) | **identical to C3** |

`6805e42` (C0) touches only `ssr.rs` and the STR dump; the fold it would have to go through
(`apply_events_into`) is confined to `open_record.rs`, so it cannot move the generic path
either.
