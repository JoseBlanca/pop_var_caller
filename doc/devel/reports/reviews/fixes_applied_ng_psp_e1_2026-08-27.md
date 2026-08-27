# Fixes applied — ng psp E1 (the review of `b9ef37e3`)

*2026-08-27. Answers [the review](ng_psp_e1_2026-08-27.md). Step E1 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. The Blocker, and it is the decision E1's report calls deliberate

**A `Truncated` fault left the live set half-advanced.** `read_changes` applied a record's
departures before it read the arrival count, so any buffer that stopped in the arrival half
returned `Truncated` — the class whose whole contract is *fetch more bytes and re-parse this record
from its first byte*. The retry then resolved the same departure positions against a set that had
already shrunk.

Measured over a record that departs one read and gains one, at each of its six cut points:

- **five of six retried to `Ok` with a read silently gone** from the live set, and gone for every
  later record of the block, because the set is carried forward;
- on a second fixture, **seven of nine cuts turned a good record into `Malformed`**.

Which of the two you got depended only on whether the shortened set still had a byte at that
position. This is verbatim what spec §8 names — *"a parse that half-advances that state before
failing corrupts every record after it, plausibly"* — and it is the sentence `block.rs`'s own retry
arm quotes. **Two of the three review agents found it independently**, one of them with a property
test that failed on its first generated case.

**And E1's report §2.3 argued for it.** The rationale was that applying the departures early kept
positions resolvable against "the set as it stands at that moment". That was not true of the code:
the departure loop had already finished before the apply, and `encode_changes` binary-searches the
*pre-departure* set on the writing side too. The eager apply bought nothing and cost the retry.

**Fixed** by decoding the whole record before applying any of it. The wire format is unchanged —
both sides already resolved a position against the set the previous record left. Two tests now hold
it: the truncation sweep retries on the *same* reader and compares against an uninterrupted read,
and a property test does the same over generated records at 600 cases. **Re-injecting the eager
apply fails three tests.**

## 2. The Majors

| what | what was done |
|---|---|
| **No property-based test for a codec**, in a module the checklist names, in a crate that already has `proptest` and uses it next door in `record.rs`. The Blocker was found by one | three property tests at 600 cases each: a truncated read retried against the uninterrupted answer; any run of records round-tripped through **one** buffer, advancing on the returned byte count; and arbitrary bytes either refused or leaving the live set a sorted set |
| **The truncation sweep could not fail on the property its docstring claimed** — it built a fresh reader inside the cut loop, so it checked the fault's *class* and never the state the fault left behind | the reader is hoisted, the retry happens on it, and its set is compared against an uninterrupted read at every cut. Filed separately from the Blocker for a reason: without this the defect is re-introducible with a green suite |
| **No fixture made either merge interleave.** Replacing `derive_changes`' `Ordering::Greater` arm and `apply_arrivals`' `else` arm with `panic!` left all 214 `ng::psp` tests green — every fixture allocated ids in walk order, so an arrival was always the largest | `an_id_that_goes_live_again_below_the_ids_allocated_since_reads_back`: a returning id sorts **below** everything allocated since it left, which is the only shape that reaches those arms, and spec §6 puts that at 83 % of ids on the human sample and 91 % on tomato. Both arms now carry a comment saying so |
| **`read_changes` could return the buffer's length instead of the bytes it used** and pass every test, because every fixture handed it a `Vec` sized to exactly one record's stream | a stream with eight bytes behind it, and a property test that reads a whole run out of **one** buffer by advancing on the count. That count is how Milestone E3 finds the next field |
| **`expect` on the writer's hot path with no `// PANIC-FREE:` marker** | marked, naming where the guarantee comes from — `LiveSetChanges`'s fields are private and its only two producers are `derive_changes` and `read_changes` |
| **`#[derive(Default)]` was a second constructor for the per-block state**, and `new()` was the path that used it — so the comment justifying the one-field struct held on only one of two paths | `Default` routes through `at_block_start`, and the unused `Clone` is gone. There is genuinely one constructor now |

## 3. The Minors worth naming

- **A departure count of three billion against a live set of 300** was bounded only by a record
  *body*'s byte ceiling, borrowed from a different container — this stream lives in the head. A
  departure is a position in the live set and the positions are strictly ascending, so
  `departures > live.len()` is exact, and it names the fault at the count rather than at whichever
  position first runs off the end.
- **A record that departs an id and arrives the same id was accepted**, with nothing saying whether
  it should be. **Decided: it is damage.** `derive_changes` puts an id in at most one of the two
  lists, so no writer produces it; honouring it would report one departure and one arrival for a
  record where nothing changed, and those two counts are what Milestone E4's residual arithmetic is
  specified against. The check now asks the set as the previous record left it, which is simpler as
  well as stricter.
- **A symmetric change to both halves of the gap arithmetic** — which keeps the codec
  self-consistent while changing what every previously written file means — was killed by exactly
  one test, incidentally, through a byte literal in a *damage* test whose docstring is about
  something else. `a_run_of_arrivals_costs_its_gaps_biased_by_one` pins the exact bytes of a
  multi-entry run on both halves; every other byte-exact test in the module holds a run of one,
  where the bias never applies.
- **Both `is_empty` accessors could be inverted** without failing any of the 214 tests, and neither
  had a caller. `an_empty_set_and_an_empty_change_say_so` gives them one and pins the claim
  `LiveSetChanges::is_empty`'s doc makes about the common case at depth.
- **`new()` on both collaborators had no doc at all**, where the question a caller has is *must I
  call `start_block` before the first record?* — the one whose wrong answer produces a file that
  parses perfectly and is wrong from a block's first record.
- **`read_changes` was 68 lines with eleven decision points**, and both kinds of `u64` — identifiers
  and positions — lived in one scope. Split into three, with the shared ascending-gap step in one
  helper so the two error strings cannot drift. The ordering that makes the encoding correct now
  reads in four lines.
- **`live()` and `changes()` after a fault** now say what they hold: `live()` is exact, which is
  what makes the retry safe, and `changes()` is not.
- The edges nothing covered — an empty slice, `u64::MAX`, a live set of a thousand — have a test.

## 4. ⚠ Two numbers of mine were wrong, and both were in prose rather than in a fixture

The review re-derived every figure in the commit message, the report and the status entry. **The
measured ones all reproduced exactly**: 3,257 bytes against 106,166 is 32.60×, "about 300 live at
once" is exactly 300, the test counts add up, and every percentage quoted from the specs matches
them verbatim. **Wrong:**

1. **`43.78` was given as the raw-identifier baseline for both corners.** It is the deep corner's
   alone: spec §6's table gives 1.020 bytes a position at 11.4 reads and 43.78 at 293. The sentence
   named the depth for both differential figures and then one baseline for both, so it reads as a
   hundredfold saving where the real pairings are **2.4× at the shallow corner and 6.8× at the
   deep one**. The number was real and the corner it belongs to was missing — and the missing
   corner is the one that flatters the design.
2. **The restatement's 12 % was quoted "at the settled block size".** Spec §6 measured it with
   blocks cut every **1,500 positions**; the settled genomic block size is 100 kb, sixty-seven
   times larger, over which the same per-block restatement amortises to a fraction of a percent.
   The figure overstated the cost of the very decision the paragraph was defending.

Both corrected in the module doc and in [the E1 report](../implementations/ng_psp_e1_2026-08-27.md).

## 5. What the review confirmed rather than found

- **No input makes this reader panic, hang or grow.** 1,500,000 fuzzed inputs — a corpus of real
  streams damaged by byte overwrites, bit flips, truncations and appends, plus uniform random bytes,
  each fed to a reader primed with intact records so the live set it resolves against is a real one.
  0 panics, 0 hangs (slowest single input 630 µs against a 200 ms deadline), and **nothing sized
  from a declared count**: the reader never held more than 976 bytes across its four vectors on
  inputs up to 96 bytes.
- **And the harness that says so can fail**, shown three ways — including catching the removed
  already-live check at case 2,053, on an unsorted live set.
- **`module_structure` filed no findings above Nit**, and the `pub(super)` widening is exactly what
  `chain_ids.rs` uses with nothing surplus.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::chain_ids` | 23 passed; 0 failed |
| `cargo test --lib ng::psp` | 224 passed; 0 failed |
| `cargo test --lib` | 4,761 passed; 0 failed; 14 ignored |

**Fifteen defects re-injected against the strengthened tests, one at a time, on a clean copy:**

| defect | tests failed |
|---|---|
| **the departures are applied before the arrivals are read** (the Blocker) | 3 |
| the writer does not restart at a block | 2 |
| the reader does not restart at a block | 1 |
| a departure is written as its identifier, not its position | 13 |
| an arrival for a read already live is accepted | 2 |
| a departure position is not bounded by the live set | 1 |
| the departure count is not bounded by the live set | 1 |
| the writer's gap loses its bias | 12 |
| **the gap loses its bias on both halves at once** — a self-consistent codec that changes what every existing file means | 4 |
| the interleaving arm of `derive_changes` drops the arrival | 2 |
| the interleaving arm of `apply_arrivals` takes the live id instead | 2 |
| `read_changes` reports the buffer's length, not what it used | 2 |
| a record naming the same id twice keeps both | 2 |
| `LiveSet::is_empty` inverted | 1 |
| `LiveSetChanges::is_empty` inverted | 1 |

Fifteen for fifteen.
