# Fixes applied — ng psp E4 (the review of `75eb7d58`)

*2026-08-28. Answers [the review](ng_psp_e4_2026-08-28.md). Step E4 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. The Blocker: a change that removed a guard

**`gap_from` was made saturating in this very commit, and that turned a refusal into silence.**

An observation's reads go on the wire as ascending gaps. Before the change, a list that was not
ascending produced a wrapped gap and the reader refused the file by name — *"a gap of
18446744073709551611 past id 7, which is past every id there is"*. After it, the same list is
accepted and names **different reads**: measured, `[3, 3]` reads back as `[3, 4]` and `[7, 3]` as
`[7, 8]`. An observation gains a read nothing folded, with neither side reporting anything, which
is spec §5's failure reached from the writer's side.

All three review agents found it, and two added the part that makes it worse: the `debug_assert`
offered as the safeguard is removed by this repo's release profile, and **the commit's own diff
weakened the one test that wrote an unsorted list** — `vec![7, 3]` became `vec![3, 7]`, a test
edited to satisfy a new precondition on a `pub` path that documents none.

**Fixed by removing the precondition rather than documenting it.** `encode_read_list` makes the
list a set itself: one pass to check, and a copy only in the branch nothing takes, since both of
ng's pileup paths already sort and deduplicate. There is nothing left for a caller to get wrong,
and `a_list_that_is_not_a_set_is_written_as_the_set_it_names` fails if the normalisation goes.

`gap_from` keeps wrapping rather than saturating, and its comment now says plainly that the arm is
**unreachable** — swapping it fails no test, because no caller can reach it — rather than implying
the line is load-bearing.

## 2. The other Blocker, and the reason both were invisible

**The residual index bound had no test.** Deleting `residual_at > declared_observations` left all
241 `ng::psp` tests green while a body claiming observation 200 of a one-observation record flipped
from refused to accepted. `a_residual_index_past_the_records_observations_is_damage` holds both
sides now — the refusal, and the sentinel that must still be accepted.

**And the multi-record fixtures never derived a residual at all.** `twelve_records_in_order` gives
chain ids to two observations per record but left `num_obs` at `a_rich_record`'s 137 — so two
identifiers against 137 reads fails the writer's own inequality, and **not one of the twelve
records ever derived anything**. Two agents measured the consequence: a probe that panics on the
derive path fired in **8 tests of 241**, and replacing `the_live_reads_of`'s body with an empty set
left `237 passed; 4 failed` — so the `&LiveSet` threaded through twelve call sites was inert. The
whole record-level oracle — the skip patterns, the truncation sweeps, *a body stands on its own* —
ran past the step it had just been extended for.

`num_obs` moves with the list now, in `record.rs`'s fixture and in `block.rs`'s two. **The same
probe fires in 20 tests.**

## 3. The guard had slack exactly where paired-end data lives, and it is closed

Spec §5 proposes bounding the derived count by the observation's read count: at most `num_obs`, at
least half of it. **That window is `num_obs / 2` wide, and its slack is exactly the number of read
pairs whose two mates both cover this record** — which is the shape paired-end data has. Measured
by the review: a residual naming two reads with `num_obs = 4`, against a live set carrying two
identifiers no observation named, derives a list of four, passes `4 <= 4` and `8 >= 4`, and the
reference allele silently gains two reads.

**The record now carries the residual's length** — one varint, `residual-read-count`, written only
when something is derived — so the check is an **equality**. The inequality stays beside it as a
second, independent statement: that the declared count could describe those reads at all.
`a_live_set_carrying_reads_no_observation_named_is_refused` is the test, in the regime where the
inequality alone would pass.

## 4. The other Majors

| what | what was done |
|---|---|
| **`observation-reads` was declared under `ChainIdChanges`**, an encoding of a different shape — two counted runs where a list is one, so a later reader stepping over it would consume seven bytes where the field is three | `FieldEncoding::ChainIdList` is its own scheme, and `skip_unknown_field` measures it as one run. Two tests: a later writer's list field stepped over correctly, and an assertion that `record_fields()` declares each of the two chain-id fields under the scheme that describes its bytes — nothing round-trips a manifest here, so that assertion is what holds it |
| **`the_chain_ids_round_trip_exactly` passed when nothing was derived at all** — its "the saving is real" assertions compared two direct `encode_read_list` calls and never touched the writer's choice | it asserts the choice, and measures the body against the same record made underivable: about the residual's own list shorter |
| **The `live_reads` parameter carried the format's most dangerous precondition, undocumented on two `pub` functions** — and the doc ninety lines above one of them still said *"These two functions take no running state and are given none"* | both say what the set must be and what a wrong one costs; the body-codec doc says which one thing in a body is not in the body |
| **The tie-break was asserted in prose and pinned by nothing** — reversing it left all 4,778 green | `a_tie_for_the_largest_list_goes_to_the_lower_index` |
| **`residual_reads`' advance loop survived deletion** | the fixtures that now derive reach it: three tests fail without it |

## 5. The Minors worth naming

- **The parity citation was half true.** Both the report and two doc comments said the guard is
  *"the same inequality the walk's own differential against production asserts"*; `parity.rs`
  asserts the **lower** bound only, and the upper is this reader's own. Corrected where it appears.
- **Five field counts in the docs were stale** — the head has five fields and the body
  twenty-two, so a record declares twenty-seven.

## 6. What the review confirmed rather than found

- **14,720,000 fuzzed inputs across seven seeds**, no panic, no hang, peak heap **10,110 bytes**,
  and nothing allocated from a declared length.
- **The harness was shown able to fail three ways**, including an injected `reserve(count)` caught
  by a heap ceiling at input 5 — *"34359739906 bytes asked for"* — which an RSS ceiling missed
  because the operating system handed out the 34 GB lazily.
- **Two oracles that matter for Milestone D held**: `Malformed` never heals when bytes are
  appended, and `Truncated` never becomes `Malformed` when bytes are removed.
- **One mutation was proved equivalent rather than reported as a survivor** —
  `decode_read_list`'s `into.clear()` — by 3,000,000 fuzz inputs giving byte-identical counters.
- **Every number in the report and commit message re-derived and holds.**

## 7. One arm reported as unreachable rather than counted

`gap_from`'s wrapping subtraction cannot be reached: the arrivals and the departure positions are
ascending by construction, and an observation's reads are made so by `encode_read_list` before they
get there. Swapping wrapping for saturating therefore fails nothing, and that is recorded at the
line — the same treatment `begin_next_block`'s decoder reset and `apply_the_changes_just_parsed`'s
idempotence guard already have.

## 8. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp` | 246 passed; 0 failed |
| `cargo test --lib` | 4,783 passed; 0 failed; 14 ignored |

**Ten defects re-injected against the strengthened tests:**

| defect | tests failed |
|---|---|
| **an unsorted list is written as it stands** (Blocker 1) | 1 |
| **the residual index bound is dropped** (Blocker 2) | 1 |
| the declared residual count is not checked | 1 |
| the residual's read count is not written | 15 |
| the tie-break is reversed | 1 |
| `residual_reads` never advances past a stored id | 3 |
| the list is declared under the changes' encoding | 1 |
| the skip arm for a list reads two counted runs | 1 |
| the derivation subtracts nothing | 2 |
| `gap_from` saturates instead of wrapping | **0 — unreachable**, §7 |

Nine for ten, with the tenth named rather than counted.
