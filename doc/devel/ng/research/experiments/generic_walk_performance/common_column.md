# The ordinary column: how common it is, and what specialising it is worth

**Two numbers first.**

**The ordinary column — one covered base where every read simply matches and there is
nothing to reconcile — is 7,898 in 10,000 columns at ~130× and 7,789 in 10,000 at 30×.**
Measured, not estimated, by counting every column of both fixtures.

**Specialising it is worth −34.3 % of the walk at ~130× and −23.9 % at 30×.** Instructions
retired, start-up subtracted, min of 3 runs a side, alternating, ranges disjoint by an order
of magnitude. All four acceptance dumps `cmp`-identical; probe counters exact.

| fixture | depth | baseline walk | with the fast lane | change |
|---|---|---:|---:|---:|
| tomato `SL4.0ch01`, 1 M loci | ~130× | 215.813 G | **141.703 G** | **−34.34 %** |
| HG002 chr1 | 30× | 110.152 G | **83.858 G** | **−23.87 %** |

Ranges: tomato 215.813–216.358 G against 141.703–141.896 G; chr1 110.152–110.220 G against
83.858–83.895 G.

The **cursor hint** — the second, independent question — is a **regression** and is reverted:
**+1.4 % at ~130×, +0.9 % at 30×**. Section 6.

---

## 1. How the frequency was measured

A census on the walk itself, behind `PVC_COLUMN_CENSUS=1`, counting every column that reaches
the fold and classifying it against the strict predicate the task named: a one-base record
opened at the walker position (no open record overlaps it), no contributor carrying an
`Insertion` or `Deletion`, no contributor whose read has a `D` op anywhere (so
`events_overlapping` over the one-base window cannot differ from `events_at`), no
`bq_zero_in_window` or `bq_override_at_walker_pos`, no depth-cap truncation, one read group.

Verbatim, `PVC_COLUMN_CENSUS=1` on the two fixtures (`ng_generic_walk_probe`, stderr):

```
=== 130x tomato 1M loci ===
census_columns=1445179
census_columns_ordinary=1141408
census_contributors=123613443
census_contributors_ordinary=95441475
census_reject_record_already_open=2776
census_reject_indel_event=3431
census_reject_read_has_deletion=167404
census_reject_mate_overlap=151946
census_reject_depth_cap=0
census_reject_multi_read_group=0
census_reject_read_has_indel=247396
census_columns_simple=1070919
census_contributors_simple=90011927
census_columns_simple_with_mate=1196245
census_contributors_simple_with_mate=103049404
=== 30x HG002 chr1 ===
census_columns=2430138
census_columns_ordinary=1892834
census_contributors=43084509
census_contributors_ordinary=30147176
census_reject_record_already_open=5267
census_reject_indel_event=3885
census_reject_read_has_deletion=310923
census_reject_mate_overlap=262962
census_reject_depth_cap=0
census_reject_multi_read_group=0
census_reject_read_has_indel=478644
census_columns_simple=1747437
census_contributors_simple=27553275
census_columns_simple_with_mate=1949151
census_contributors_simple_with_mate=33214831
```

Read as natural frequencies:

| | ~130× tomato | 30× HG002 chr1 |
|---|---:|---:|
| columns | 1,445,179 | 2,430,138 |
| reads folded (contributors), mean per column | 123,613,443 (85.5) | 43,084,509 (17.7) |
| **columns that are ordinary** | **7,898 in 10,000** | **7,789 in 10,000** |
| **reads folded in an ordinary column** | **7,721 in 10,000** | **6,997 in 10,000** |

The read-weighted number is the one that matters for cost, and it is lower than the column
number at 30× because ordinary columns are shallower than average there.

**What disqualifies the other two columns in ten**, as a share of all columns:

| reason | ~130× | 30× |
|---|---:|---:|
| some read in the column has a `D` op anywhere in its CIGAR | 1,158 in 10,000 | 1,279 in 10,000 |
| a mate overlap fires (two contributors share a chain id) | 1,051 in 10,000 | 1,082 in 10,000 |
| a record is already open over this base | 19 in 10,000 | 22 in 10,000 |
| an `Insertion` or `Deletion` is anchored here | 24 in 10,000 | 16 in 10,000 |
| depth cap, or more than one read group | 0 | 0 |

The two reasons that matter are each about one column in ten, and they overlap, which is why
they do not sum to the 21 % shortfall.

---

## 2. What was built

`src/ng/locus_generation/pileup/fast_column.rs`, a new 322-line module. The general path is
untouched and takes every column that does not qualify.

### The predicate, decided before any per-read work

Three tests, in this order, all before the loop:

1. `active_reads.len() <= max_snp_column_depth` and the set is non-empty. The contributor
   count is at most the active-read count, so this bounds the column cap from above without
   knowing which reads contribute.
2. `open_records.find_overlapping(p, p + 1).is_none()`. Then the record this column would
   open is one base wide, anchored at the walker, and can never be found again — the next
   position's events start at `p + 1`, and half-open intervals that touch do not overlap. No
   widen, no re-fold, no read folded twice.
3. Per read, `cursor.matches_only()` — a new one-bit field on `CigarCursor`, true when the
   CIGAR has no `I` and no `D` op. The first read that fails it aborts the pass.

A fourth is decided *after* the pass, from a sort the general path runs anyway: no two
contributors share a chain id. A column that fails it is handed back having cost only the
pass.

**Handing back is exact.** The pass writes only into scratch this module owns and sets
`ever_contributed`, which the general path sets for the same reads a moment later. Nothing
else is touched.

`matches_only` is *conservative* where the strict census predicate is exact — a read with an
insertion 40 bases away disqualifies the column. That is the whole gap between the census's
7,898 in 10,000 and what actually fires:

| | ~130× | 30× |
|---|---:|---:|
| strict predicate (census) | 1,141,408 | 1,892,834 |
| what the built predicate admits (census `columns_simple`) | 1,070,919 | 1,747,437 |
| **columns the fast lane actually took** (`fast_columns`) | **1,069,716** | **1,746,198** |

74.0 % and 71.9 % of all columns. The 1,203 and 1,239 columns of difference are where the
gate's use of the *active set* rather than the *contributor list* for the indel test is
conservative.

### What the fast lane does per read

A 40-byte `PlainContribution`: `read_id`, `chain_id`, `read_group`, `base`, `ln_q`, forward
bit, placed-left bit, `mapq`. `cursor.match_at(p, read)` returns `Option<(base, bq)>` —
`events_at`'s answer in scalars, since for a `matches_only` cursor the other two arms of
`events_at` are unreachable. Then one sort by `read_id`, one fetch of the reference base, one
grouping pass over at most a handful of `(base, read_group)` buckets, and the locus.

It replaces, per read per column: an `EventsAt` `SmallVec` of 40-byte enums; a
`ReadContribution` push; `find_overlapping`; `open_new` and its reference fetch;
`apply_events_into` (clear, reserve, loop, push, run-accumulate); `WitnessedRefPositions::take_from`
and `canonicalise_runs`; an `AHashMap` entry, hash and probe; a `FoldedReadState` store; and
at close a `Vec` collect, a `witness_of`, a `read_agreed_with_reference`, and a linear search
over the observation list.

### The three places byte-identity was at risk

**`q_sum` order.** The general path adds one `f64` term per read per observation in ascending
`read_id`, established by `keyed_observations_counting`'s sort. The fast lane sorts its own
compact buffer by `read_id` and accumulates in that order — the same sequence of additions.
This was checked first, and it held: the tomato generic dump (1,718,914 lines, `q_sum` printed
per observation) is `cmp`-identical.

**Chain ids.** `read_agreed_with_reference` reduces, for a one-base record with a complete
witness, to `base == reference base` — confirmed against the code, not assumed: bucket 0
holds the record's reference bytes, a read matching them lands there, and the function's
`allele_index == 0` branch then answers `first == 0 && past_last == reference.len()`, which a
one-base complete witness satisfies. So the fast lane pushes a chain id exactly when the
read's base differs from the reference, then sorts and dedups as `finalise` does.

**Observation order.** `finalise` sorts on `(bases, read_witness, read_group)`. Every witness
here is `Complete` and every `bases` is one byte, so that is `(base, read_group)`.

### The one step that had to be reproduced, and was not obvious

The fast lane finishes the locus *at* the position it walks; the general path leaves a
one-base record open and drains it one step later. **That one step is observable**, and two
tests found it — a fully-consumed walk cannot show it, and the four dumps were already green
when both failed:

- `parity::both_walkers_report_the_same_error_on_the_same_malformed_input` — a walk that
  **aborts** (a reference fetch past the contig end) loses whatever is still open. Emitting a
  step early handed the consumer one locus more before the error: *"ng emitted 6 stream items,
  production 5"*.
- `generator::tests::an_abandoned_walk_does_not_leak_its_active_reads_into_the_next_region` —
  a walk **abandoned** part-way stops where its consumer stopped pulling, and a locus offered
  a step early moves that point back a position, so fewer reads are admitted and fewer chain
  ids allocated: *"left: 3, right: 4"*.

The fix is a one-slot `WalkerState::sealed` holding the locus for exactly one step. One slot
is always enough, and for the same reason the predicate gives: the fast lane fires only when
no record covers its base, so every record still open ends at or before it and the table is
empty by the next position — which is also why emitting the held locus first is always
coordinate order. It is also read by `reached_stop`, which asks the record table for its
lowest anchor; without that a region ending on an ordinary column stops one position early
and loses the out-of-region record the next column would have produced (91,572 of them across
chr21's 102,938 regions, in the first draft).

---

## 3. Instruction counts

Instrument: `instructions retired` from `/usr/bin/time -l`, floor-subtracted, min of 3 runs a
side, **alternating baseline/fast within one script** so a drift on this shared host hits both
sides. Two separately-built release binaries — the baseline is a pristine build of
`6fbbd093764662ed2496acde39424c8ee234ea1c`, not a switched-off code path.
`PVC_TRUST_REFERENCE_INDEX=1` throughout (`reference_check=trusted_unverified` on every run).

Raw, verbatim (`tmp/fast_column/final_ab.txt` in the worktree):

```
== floors (baseline binary / fast binary) ==
          1307993164  instructions retired [floor tom base]
          1310781790  instructions retired [floor tom fast]
           350275601  instructions retired [floor chr1 base]
           349647510  instructions retired [floor chr1 fast]
          1306114293  instructions retired [floor tom base]
          1309525654  instructions retired [floor tom fast]
           349323741  instructions retired [floor chr1 base]
           349580883  instructions retired [floor chr1 fast]
          1308811856  instructions retired [floor tom base]
          1310382851  instructions retired [floor tom fast]
           349768135  instructions retired [floor chr1 base]
           349779794  instructions retired [floor chr1 fast]
== tomato ~130x, 1M loci ==
loci=1000000 seconds=9.350            386908160  maximum resident set size         217663789865  instructions retired [tom base]
loci=1000000 seconds=6.541 fast_columns=1069716            391495680  maximum resident set size         143096515559  instructions retired [tom fast]
loci=1000000 seconds=9.267            386809856  maximum resident set size         217196103029  instructions retired [tom base]
loci=1000000 seconds=6.565 fast_columns=1069716            386662400  maximum resident set size         143012953126  instructions retired [tom fast]
loci=1000000 seconds=9.312            386596864  maximum resident set size         217119415429  instructions retired [tom base]
loci=1000000 seconds=6.501 fast_columns=1069716            390348800  maximum resident set size         143205163345  instructions retired [tom fast]
== HG002 chr1 30x ==
loci=1541788 seconds=6.226             21299200  maximum resident set size         110568634183  instructions retired [chr1 base]
loci=1541788 seconds=5.254 fast_columns=1746198             19775488  maximum resident set size          84209849368  instructions retired [chr1 fast]
loci=1541788 seconds=6.271             21741568  maximum resident set size         110566355651  instructions retired [chr1 base]
loci=1541788 seconds=5.242 fast_columns=1746198             21102592  maximum resident set size          84207401354  instructions retired [chr1 fast]
loci=1541788 seconds=6.208             21479424  maximum resident set size         110501387445  instructions retired [chr1 base]
loci=1541788 seconds=5.208 fast_columns=1746198             18956288  maximum resident set size          84244512409  instructions retired [chr1 fast]
```

Floor-subtracted walk instructions:

| fixture | side | min | range | floor | walk (min) |
|---|---|---:|---:|---:|---:|
| tomato ~130×, 1 M loci | baseline | 217.119 G | 217.119–217.664 | 1.306 G | **215.813 G** |
| | fast lane | 143.013 G | 143.013–143.205 | 1.310 G | **141.703 G** |
| HG002 chr1 30× | baseline | 110.501 G | 110.501–110.569 | 0.3493 G | **110.152 G** |
| | fast lane | 84.207 G | 84.207–84.244 | 0.3496 G | **83.858 G** |

**−34.34 % at ~130×, −23.87 % at 30×.** Disjoint on both, with the two sides' ranges more than
a thousand times their own widths apart.

**Wall time is not admissible** on this host (6 performance and 12 low-energy cores, three
other agents measuring) and is quoted only because it agrees: tomato 9.267–9.350 s →
6.501–6.565 s; chr1 6.208–6.271 s → 5.208–5.254 s.

**Peak RSS is neutral.** Tomato 386.6–386.9 MB against 386.7–391.5 MB; chr1 21.3–21.7 MB
against 19.0–21.1 MB. Both inside the run-to-run spread already recorded for these fixtures.

An earlier reading of the same A/B, taken with the same binary switched off by
`PVC_FAST_COLUMN=0` rather than against a separate baseline build, gave −34.26 % and
−23.68 %. The two agree to within a tenth of a point, so the field and the parameter the
change adds to the general path cost nothing measurable.

---

## 4. Gates

All four acceptance dumps **`cmp`-identical** to the stored copies in
`tmp/perf_review_2026-08-04_ng-generic-walk/`, verbatim:

```
OK generic chr21
OK ssr chr21
OK generic tomato
OK ssr tomato
  251792 …/generic_chr21.txt
    4406 …/ssr_chr21.txt
 1718914 …/generic_tom.txt
   11945 …/ssr_tom.txt
```

Probe counters on chr21, baseline then fast lane:

```
loci=236081
observations=251786
reads_admitted=54709
loci=236081
observations=251786
reads_admitted=54709
fast_columns=262498
```

Exact, and **262,498 columns took the fast lane** on that walk — the counter is unconditional
(one relaxed atomic increment per ordinary column, against the ~18 reads that column would
otherwise have folded) precisely so a fast lane that silently never fires cannot look like one
that works.

Validation, in debug:

- `cargo test --lib` — **2,882 passed; 0 failed; 5 ignored**. Same as the clean tree.
- `cargo test --examples` — 33 targets, all ok.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo doc --no-deps` — 12 unresolved links, the recorded baseline; no new ones.

**The cross-check that makes a second path defensible.** In debug builds the fast lane asserts,
per read per column, that `match_at` and `events_at` give the same answer:

```rust
debug_assert_eq!(
    active.cursor.match_at(walker_pos, &active.read),
    match active.cursor.events_at(walker_pos, &active.read).first() {
        Some(ReadEvent::Match { base, bq_baq, .. }) => Some((*base, *bq_baq)),
        Some(_) | None => None,
    },
    …
);
```

It is placed where the two paths genuinely duplicate a computation, and only there —
everything after that loop is a sum over what it produced. It is armed in every debug walk in
the suite; the parity census alone runs it over ~257,000 loci.

---

## 5. The honest complexity cost

**Two code paths that must agree, through the hottest code in the walk**, and the gate only
catches a disagreement on inputs that exercise both. Specifically:

- **The fast lane re-derives `read_agreed_with_reference`, `witness_of`, the observation key
  and the observation sort as closed-form special cases.** Each is a paragraph of reasoning
  in `fast_column.rs` rather than a call. If any of the four general-path rules changes —
  what makes an observation's identity, how a witness resolves to `Complete`, how chain ids
  are attributed — the fast lane silently keeps the old rule. Nothing in the type system says
  so. The debug cross-check does not cover this: it pins `match_at` against `events_at`, not
  the closed forms against `finalise`.
- **The one-step emission delay is a shared invariant with no single owner.** `sealed`,
  `close_aged_records_into`, `flush_chromosome_into`, `begin_region` and `reached_stop` all
  have to agree that a held locus is a record. Two of the five were wrong in the first draft
  and only two tests in 2,882 found it.
- **`matches_only` is a second, coarser statement of the same thing `spans_only_its_anchors`
  says.** Two predicates about indels now live on `CigarCursor`, one exact for the fold's
  reuse rule and one conservative for the fast lane's gate.
- **The `fold_capacity` heuristic no longer sees ordinary columns.** It learns the fold
  table's size from the last record to close, and the fast lane closes no records — so it
  now learns only from general-path records. Nothing observable depends on it, and the
  measurement above is with it in that state, but it is a coupling nobody wrote down.
- **`PVC_FAST_COLUMN=0` and the `PVC_COLUMN_CENSUS=1` block are measurement scaffolding left
  in the diff.** The census costs the general path one `OnceLock` load per column when off;
  that cost is inside the measured numbers on both sides. Both should be decided on before
  this lands.

What the change does **not** cost: no `unsafe`, no new dependency, no lifetime, no build flag,
no change to any emitted byte, and no change to the general path's own code.

---

## 6. The cursor hint, measured separately — a regression, reverted

`CigarCursor` is stateless and rescans the CIGAR from op 0 at every position, and the walk
queries each read at monotonically increasing positions. A per-read `Cell<u32>` holding the op
index that answered last time turns the scan into a check whenever
`offsets[hint].ref_pos <= walker_pos` still holds, falling back to a scan from op 0 when it
does not. Applied to `events_at_linear` and to the fast lane's `match_at`; `events_overlapping`
left alone, since it queries windows that start behind the walker. Answers are unchanged, which
is what `cursor_queries_are_stateless` asserts, so that test needed no edit, and all four dumps
stayed `cmp`-identical.

Measured on top of the fast lane, three runs a side, alternating:

| fixture | depth | fast lane | fast lane + hint | change |
|---|---|---:|---:|---:|
| tomato `SL4.0ch01`, 1 M loci | ~130× | 141.658 G | 143.692 G | **+1.44 %** |
| HG002 chr1 | 30× | 83.854 G | 84.565 G | **+0.85 %** |

Ranges disjoint on both. **The prize was small and the bookkeeping is bigger than it.** An
Illumina CIGAR under this cursor is one to three ops, so the scan the hint skips is one or two
integer compares, while the hint costs a `Cell` load, a bounds check, a running
"last op at or before" variable through the loop body, and a `Cell` store. The task's
prediction — *"for a 1–3-op Illumina CIGAR the prize is small"* — holds, and on these two
fixtures the read population never reaches the many-op regime where it would not.

**Reverted.** Not carried in the diff below. Recorded so it is not re-run.

---

## 7. The next lever, sized

**Mate overlap is the largest single thing the fast lane hands back**, and it is the one that
could plausibly be taken. It disqualifies 1,051 in 10,000 columns at ~130× and 1,082 in 10,000
at 30×, and admitting it would raise the read-weighted coverage by 9 to 13 points:

| | ~130× | 30× |
|---|---:|---:|
| reads folded in a column the fast lane takes | 7,282 in 10,000 | 6,395 in 10,000 |
| the same, if mate overlap were handled | 8,336 in 10,000 | 7,709 in 10,000 |

For a one-base record with only `Match` events, everything `resolve_mate_overlap_at_pos` does
reduces to rewriting one `bq` scalar per read — `bq_zero_in_window` is "this read's quality is
zero here" and `bq_override_at_walker_pos` is "this read's quality here is *this* instead". So
it is expressible. It is not free: it needs the pair-finding sort's *result* rather than just
its yes/no, and it reintroduces the tie-break rules the fast lane currently gets to not know
about — a second copy of the reconciliation, which is exactly the complexity cost §5 is about.
**Not attempted, and not recommended without a measurement of its own.**

---

## 8. Where the code is

Worktree: `/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-af6f29682c6996432`,
detached at `6fbbd093764662ed2496acde39424c8ee234ea1c`, left as measured.

- Diff of tracked files: `tmp/fast_column/fastlane_final.diff` (521 lines)
- The new module: `src/ng/locus_generation/pileup/fast_column.rs` (322 lines, untracked;
  copy at `tmp/fast_column/fast_column.rs.final`)
- Raw measurement output: `tmp/fast_column/final_ab.txt`
- Scripts: `tmp/fast_column/gate.sh` (four dumps), `tmp/fast_column/ab2.sh` (the A/B),
  `tmp/fast_column/ab3.sh` (the cursor hint)
- The two binaries the A/B compared: `tmp/fast_column/probe_baseline`,
  `tmp/fast_column/probe_fast`

`git diff --stat`:

```
 examples/ng_generic_walk_probe.rs              |  14 ++
 src/ng/locus_generation/pileup/cigar_cursor.rs |  84 ++++++++++
 src/ng/locus_generation/pileup/genome_walk.rs  | 214 ++++++++++++++++++++++++-
 src/ng/locus_generation/pileup/mod.rs          |  55 +++++++
 src/ng/locus_generation/pileup/open_record.rs  |   2 +-
 5 files changed, 361 insertions(+), 8 deletions(-)
```

`open_record.rs`'s one changed line makes `phred_to_ln_perr` `pub(super)`.
`copy_fidelity.rs`'s two pinned files, `decompose.rs` and `chain_id_allocator.rs`, are
untouched.

Copies alongside this report, so it is self-contained:

- `common_column_fastlane.diff` — the tracked-file diff
- `common_column_fast_column.rs` — the new module
- `common_column_final_ab.txt` — the raw A/B output quoted in §3
