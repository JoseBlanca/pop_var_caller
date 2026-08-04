# Milestone C — module structure & API fitness review

Commit reviewed: `4e0e02f` (head of Milestone C, both prior reviews applied).
Baseline reproduced: `ng::locus_generation` **304 passed**, and restored green after every experiment.
Everything below that says "verified" was compiled; the compiler's own words are quoted.

Method note: all experiments ran in two container starts, applying a variant, compiling, and
restoring pristine files between each. The tree is back at `4e0e02f` with no changes.

---

## Verdicts on the four recorded departures

### 1. `apply_events_into` returns `bool` and fills a caller-owned buffer — **right call, wrong residue**

The departure itself is correct and the justification is the arch's own. Arch §2's *Contract*
paragraph says "the runs are accumulated into a buffer the caller owns and the callee clears,
like `allele_seq`, so a fold allocates nothing per read", and `refill_from` must swap so a
re-fold of a spilled witness allocates nothing per (read × widen). A function that *returns*
the set takes the buffer's storage with it. The signature sketch and the contract in the same
section are inconsistent; the implementer resolved in favour of the contract, which is the
deliverable. Sound.

What the `bool` costs is finding **F1** below: the fact now lives in two places and three sites
reconcile them by hand.

### 2. `from_left` / `from_right` return `Option<Self>` — **right call**

Verified from the code rather than the report. Arch §5's "call sites unchanged" is a promise
about **STR output** (spec §1 goal 3), and that promise is kept — the dump is byte-identical.
The `Option` bought something concrete: `3fecdf6` was able to delete the duplicate
`tract.is_empty()` guard in `classify_read` because the constructor already answered first
(`ssr.rs:874-883`), collapsing two decisions that could drift into one. A `-> Self` signature
would have had to panic, fabricate a position, or leave `Partial { 0, 0 }` expressible on a
type whose producer had just been fixed for minting 6,704 of them.

### 3. Arch §1.1's "constructors return `Complete` when the result covers the whole locus" — **claim VERIFIED, and understated**

The claim is reachable and the consequence is real, but the report justifies it on *output
movement* when the stronger argument is *correctness*.

Reachability, from the code: `ssr.rs:863` computes `reach` in **read** bases
(`tract.end - tract.start`), `ssr.rs:864` computes `locus_len` from the **reference** tract
length, and `ssr.rs:878-880` hands one to the other. An expanded allele reads further than the
reference tract has positions, so `from_left(reach ≥ len, len)` is an ordinary input, not a
corner. It is already pinned by a fixture:
`ssr::tally::tests::an_expanded_allele_merges_the_two_sides_into_one_observation`
(`ssr.rs:1372-1401`) builds `from_left(9, LocusLen(6))` and asserts
`counts.observations_partial == 2`.

Output movement, confirmed: under the arch contract both of those become `Complete`, so
(a) `observations_partial` → 0 and `observations_complete` → 2, which is the STR dump's
`# obs_complete= obs_partial=` header line; and (b) they enter `complete_observations()`
(`mod.rs:147`), which is the STR dump's `depth` column on **every row of the locus**. Two
columns of the byte-identity oracle move. Claim verified.

**The stronger reason, which the report does not give.** A read anchored at *one* border whose
read-coordinate reach happens to equal or exceed the reference tract length has **not reached
the second border** — it ran out of read. Its evidence is a lower bound. `Complete` is defined
as "the read reached both borders and witnessed every position between them" (`witness.rs:198`)
and is the gate on `complete_observations()`, i.e. on what a likelihood may score as an exact
allele length. Implementing arch §1.1's contract would therefore score a lower bound as an exact
length — a defect, not merely a moved byte. **Recommendation for D3:** record the semantic
reason alongside the output one, because the output argument alone invites "rebaseline the
oracle and implement it".

### 4. C2 absorbed C4, D1 and D2 — **right call, forced**

Spot-checked and true: `sort_key`'s `(u8, u16, u16)` cannot be a total order over a set (C4),
`num_obs_along_locus` read the two `u16` fields directly (D1, now `mod.rs:110-125`), and the
flush predicates delegate in one line each (D2). None compiles without the payload swap. The
thing that mattered — C3, the only behaviour change — still lands alone. The residue is that
D4 lost part of its content to C2 and what remains is label work only, which the plan records.

---

## Can the API support Milestone D? — yes, verified by compiling D's next steps

I wrote each of D's remaining steps against the current types and compiled them.

### `ReadWitness::from_run` (arch §1.1) — composes, no missing piece

```rust
pub fn from_run(offset_in_locus: u16, positions_covered: u16, locus_len: LocusLen) -> Option<Self> {
    let start = offset_in_locus.min(locus_len.get());
    let end = offset_in_locus.saturating_add(positions_covered).min(locus_len.get());
    WitnessedLocusPositions::from_half_open_runs([(start, end)]).map(|positions| Self::Partial { positions })
}
```

Compiles; `from_run(3,4,len) == partial_run(3,4)`, `from_run(8,40,len) == partial_run(8,2)`,
`from_run(12,4,len) == None`, `from_run(3,0,len) == None`. Note it must *not* borrow
`one_run_from_offset_and_length`, which rejects rather than clamps — the interior constructor
needs the clamp its two siblings have, so the shape is `from_half_open_runs` over pre-clamped
ends. **Test result: 19 passed, including the two probes.**

A multi-run `from_runs(runs, locus_len)` also compiles and is the only constructor that can
answer `Complete`; it is the shape `witness_of` currently hand-rolls (`open_record.rs:208-238`).
Worth considering at D3 so the "clamp each run, then decide `Complete` on total coverage" rule
has one home instead of two.

### D4 — both dump tools print the set: **compiles, no API gap**

Rewriting `ng_ssr_loci_dump::witness_label` to `-> String` emitting
`partial:{side}:{off}+{len},...` compiles against the public API (`ReadWitness::Partial`'s field
is public, `WitnessedLocusPositions::runs()` is public). The only churn is the `ObservationRow`
field type and two test expectations that assert the exact label
(`ng_ssr_loci_dump.rs:659`, `:768-774`) — D4's own work, not a missing accessor.

### D5 — census counts holed witnesses **and the positions inside them**: compiles, but see F7

Added to `DivergenceCensus`, next to `fabricated_ref_bases` (`parity.rs:1738`):

```rust
if positions.runs().len() > 1 {
    self.holed_witness_reads += u64::from(observation.num_obs);
    let first_start = positions.runs().next().expect("a witnessed set is never empty").0;
    let last_end    = positions.runs().last().expect("a witnessed set is never empty").1;
    self.hole_positions +=
        (u64::from(last_end - first_start) - u64::from(positions.positions_covered()))
        * u64::from(observation.num_obs);
}
```

Compiles and measures correctly (probe asserts 4 hole positions over `[(0,3),(7,10)]`). But it
needs two `expect`s on an invariant the constructor already enforces — see **F7**.

---

## Visibility and placement — the boundary is real (compiler-verified)

Three bypass attempts, one `cargo check --lib --all-features`:

```
error[E0603]: module `witnessed_ref` is private
   --> src/ng/locus_generation/mod.rs:33:45
    | fn probe_v1c_name_the_ref_type() -> pileup::witnessed_ref::WitnessedRefPositions {
    |                                             ^^^^^^^^^^^^^  struct `WitnessedRefPositions` is not publicly re-exported

error[E0423]: cannot initialize a tuple struct which contains private fields
   --> src/ng/locus_generation/pileup/open_record.rs:197:5
    |     WitnessedRefPositions(SmallVec::from_slice(&[(9u32, 4u32), (0, 5)]))
note: constructor is not visible here due to private fields

error[E0423]: cannot initialize a tuple struct which contains private fields
  --> src/ng/locus_generation/mod.rs:29:5
   |     WitnessedLocusPositions(smallvec::SmallVec::from_slice(&[(9u16, 4u16), (0, 5)]))
note: constructor is not visible here due to private fields
```

The canonicalising constructors **cannot** be bypassed from anywhere in the crate. Both fields
are private to their own file, so the `pub(super)`/`pub` type visibility never widens
constructibility. `WitnessedRefPositions` is in fact tighter than arch §2 promises: because
`mod witnessed_ref;` is private in `pileup/mod.rs:109`, the type is not merely un-nameable
outside `pileup`, its *module path* is unusable there. Splitting it out of `open_record.rs`
(commit `5790914`) was the right move and it did what it claimed.

The only remaining route to a non-canonical witness is the public
`WitnessedLocusPositions::from_half_open_runs`, which canonicalises — so what an external caller
can build is an *unclamped* set, never an unsorted or overlapping one. That is the documented
position (`mod.rs:91-103`, the clamp in `num_obs_along_locus` is the guard) and it is correct.

---

## The two `expect(dead_code)` attributes — reasons accurate, guard works, shape improvable

Both stated reasons are true: `from_half_open_runs` and `positions_covered` on
`WitnessedRefPositions` are reached only from `#[cfg(test)]` code, and the fold really does go
through `take_from`/`refill_from` and really does measure against the footprint.

Wiring either into non-test code (`witness_of`) does fail the build **as intended**:

```
warning: this lint expectation is unfulfilled     # cargo check --lib
  --> src/ng/locus_generation/pileup/witnessed_ref.rs:67:13
error: this lint expectation is unfulfilled       # cargo clippy --lib -- -D warnings
   = note: `-D unfulfilled-lint-expectations` implied by `-D warnings`
error: could not compile `pop_var_caller` (lib) due to 2 previous errors
```

Two caveats, and one is a cheap improvement — see **F9**.

---

## Findings

### F1 — **Major.** `apply_events_into`'s `bool` puts one fact in two places, and it can be silently ignored

*File:* `src/ng/locus_generation/pileup/open_record.rs:1278-1284` (signature), `:1442-1448`,
`:1513-1565`, `:1691-1770`.

The function's contract is "`true` leaves at least one non-empty run in `witnessed_runs`;
`false` leaves it empty" (`:1248-1254`). That is one fact stored twice, and the three call sites
reconcile the copies by hand: two `expect("apply_events_into reported a witnessed run")`
(`:1447-1448`, `:1564-1565`) and one `debug_assert!(refilled, ...)` (`:1767-1770`).

**The sharper problem is that the answer is ignorable.** The `bool` carries no `#[must_use]`, so
`apply_events_into(buf, runs, pos, seq, events);` compiles and the caller then folds a read that
witnessed nothing. That is the *exact* hazard this milestone reshaped `canonicalise_runs` to
close — its own doc says so (`witness.rs:27-32`): "**The buffer flows through the return value**
so that ignoring the answer cannot compile… no lint says a word." The rule was applied to the
helper and not to the hot function it exists for.

*Fix, compiled and green:* delete the return type. `WitnessedRefPositions::take_from` already
answers `Option`, and it is the type that owns the invariant.

```rust
pub(super) fn apply_events_into(..) { /* no return; push the last run if there is one */ }

// fold_read_into_record
apply_events_into(allele_seq_buf, witnessed_runs_buf, rec_pos, &alleles[0].seq, window_events);
let Some(witnessed) = WitnessedRefPositions::take_from(witnessed_runs_buf) else { /* drop path */ };

// refold_live_reads — the refill must come after the bucket work, so it reads the buffer
if witnessed_runs_buf.is_empty() { /* drop path */ continue; }

// apply_events (test wrapper)
let witnessed = WitnessedRefPositions::take_from(&mut runs)?;
```

*Evidence:* this variant is clean under `cargo clippy --lib --all-features -- -D warnings` and
`ng::locus_generation` is **304 passed, 0 failed**. It removes both `expect`s and the `bool`.
(If the shape is kept instead, `#[must_use]` on the return is the one-line floor.)

### F2 — **Minor.** `sort_key` leaks the encoding that `WitnessedLocusPositions` documents as private

*File:* `src/ng/locus_generation/witness.rs:384-389`, against `:92`.

`:92` states the design promise: "The encoding is private behind [`runs`] so it can still move."
`:384` breaks it in the public API of a publicly re-exported type:

```rust
pub fn sort_key(&self) -> (u8, &[(u16, u16)]) {
    ...
    Self::Partial { positions } => (1, &positions.0),
}
```

`&[(u16, u16)]` *is* the encoding. Any of the moves the doc reserves — a `Run` newtype (arch §4
lists it as an impl-time confirmation), a bitmask, three inline runs — changes this signature and
every caller. Arch §4's "the encoding stays private" is not true today.

*Fix:* derive `PartialOrd, Ord` on `WitnessedLocusPositions` — well-defined precisely because the
representation is canonical, which is this type's whole reason for existing — and return
`(u8, Option<&WitnessedLocusPositions>)`. `Complete` is `(0, None)`, `Partial` is `(1, Some(&p))`;
`None < Some` and the tag already separates them, so the order is unchanged. `ReadWitness` still
gets no `Ord`, which was the stated objection (`:381-383`).

### F3 — **Minor.** `witness_of` narrows run *coordinates* through `LocusLen`, saturating where spec §3.4 forbids a clamp

*File:* `src/ng/locus_generation/pileup/open_record.rs:215-216` (`:222-223` with W3 applied).

```rust
LocusLen::from_positions(u64::from(first - record_pos)).get(),
LocusLen::from_positions(u64::from(past_last - record_pos)).get(),
```

Neither value is a locus length — they are a run's start offset and its end offset. `LocusLen` is
being used as a `u64 → u16` **saturating cast helper**, which is (a) exactly the
"two same-shaped quantities with no way to tell them apart" confusion the newtype's own doc
(`witness.rs:233-245`) says it was minted to prevent, now committed from the inside; and (b) a
silent clamp where spec §3.4 says "the out-of-range case is an error or an assertion, not a
clamp" and arch §2 says the requirement is met "by a gate that already exists and a
`debug_assert`". The gate (`PileupGeneratorConfig::check`) does exist, which is precisely why the
truthful spelling costs nothing.

*Fix:* a file-local narrowing that says what it is and cannot lie —
`fn locus_offset(delta: u32) -> u16 { u16::try_from(delta).expect("the footprint is capped at MAX_RECORD_SPAN_CEILING by PileupGeneratorConfig::check") }`.

### F4 — **Minor.** `witness_of`'s `expect` is the right *shape*; its message is not, and a second failure mode in the same function is unnamed

*File:* `src/ng/locus_generation/pileup/open_record.rs:228-229`, `:199-203`.

On the shape: `expect` is right. `witness_of` is `pub(super)` with two production call sites,
both inside `finalise` (`:599`, `:719`), and the precondition is established by the fold, not by
the caller — an `Option` would make both call sites invent a policy for a state the walk cannot
produce, and the reliability review already judged the unreachability argument sound.

On the message, probed by forcing it:

```
PROBE V7 empty-after-clamp message:
  a folded read witnessed at least one position inside the record it folded into
```

That is the invariant that was *supposed* to hold, stated as a fact. A maintainer gets no runs,
no `record_pos`, no `record_end_exclusive`, and no hint of where the invariant is established.
The panic says what should have been true, not what was.

There is also a second, unnamed way out of this function. An inverted footprint never reaches the
`expect` — it dies inside `u32::clamp`:

```
PROBE V7b inverted-footprint message: min > max. min = 21, max = 5
```

and the width guard above it waves the case through, because `record_end_exclusive.saturating_sub(record_pos)`
reads `0` on an inversion. std's message names neither ng nor the record.

*Fix, compiled and green (`clippy -D warnings` clean, 304 passed):*

```rust
debug_assert!(record_pos <= record_end_exclusive,
    "inverted record footprint {record_pos}..{record_end_exclusive} reached `witness_of`");
...
let positions = WitnessedLocusPositions::from_half_open_runs(clamped).unwrap_or_else(|| panic!(
    "every run of {witnessed:?} fell outside the record footprint {record_pos}..{record_end_exclusive}; \
     `apply_events_into` clips a read's runs into the record when it folds, a record's anchor never \
     moves and a widen only extends the right, so a folded read's witness cannot lie outside the \
     final footprint"));
```

### F5 — **Major.** The four-term reason sum in two dump tools will drift, and nothing stops it

*Files:* `examples/ng_ssr_loci_dump.rs:254-257`, `examples/ng_ssr_aligner_bakeoff.rs:376-379`,
against `src/ng/locus_generation/ssr.rs:177-199`.

Both tools hand-write `no_border_anchored + low_quality + window_truncated + outside_tract`. The
comment at each site records the last failure ("summing three of them would report a fraction…
on tomato chr01 it under-reported by 6,704 of ~9,265 reads") — a comment is not a mechanism, and
the same slip is one field away.

*Evidence it is silent (probe V5):* adding a fifth reason to `SsrGeneratorCounts` and changing
nothing else —

```
cargo clippy --lib --examples --all-features -- -D warnings
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 5.94s
```

Clean. Both dumps now under-report and no test, lint or type says a word.

*Evidence the fix works (probe V6):* moving the sum onto the counts type, written over an
exhaustive destructuring, turns the same edit into a compile error:

```
error[E0027]: pattern does not mention field `probe_fifth_reason`
   --> src/ng/locus_generation/ssr.rs:208:13
    | 208 |  let Self {
    | ...   |
    | 217 |  } = self;
    | |_________^ missing field `probe_fifth_reason`
error: could not compile `pop_var_caller` (lib) due to 1 previous error
```

*Fix:*

```rust
impl SsrGeneratorCounts {
    /// Every reason, and a new one is a compile error rather than a silent under-report.
    pub fn reads_without_observation(&self) -> u64 {
        let Self { reads_fetched: _, reads_discarded_by_cap: _, observations_complete: _,
                   observations_partial: _, no_border_anchored, low_quality,
                   window_truncated, outside_tract } = self;
        no_border_anchored + low_quality + window_truncated + outside_tract
    }
}
```

Both tools call it. (The per-locus scalar in `ssr::tally` is already safe — its `match` on
`NoObservationReason` is exhaustive, so a new variant forces an arm there.)

### F6 — **Minor.** Four `witness_label`s, three spellings, and two of them are already internally inconsistent

| file | line | strings emitted |
|---|---|---|
| `examples/ng_ssr_loci_dump.rs` | 163 | `complete` / `partial:left` / `partial:right` / `partial:interior` |
| `examples/ng_ssr_cohort_stutter.rs` | 154 | `complete` / `partial_left` / `partial_right` / **`partial:interior`** |
| `examples/ng_ssr_aligner_bakeoff.rs` | 198 | `complete` / `partial_left` / `partial_right` / **`partial:interior`** |
| `examples/ng_generic_loci_dump.rs` | 304 | `complete` / `observed:<off>+<len>,…` |

The three STR copies carry a **byte-identical seven-line comment** and differ only in string
literals — which is what drift looks like when it is copy-paste. Two of them already mix
separators *inside one function*: `partial_left`, `partial_right`, `partial:interior`. That is
not a hypothetical future drift; it is present today, and it means a consumer grepping
`partial:` gets one tool's left/right rows and another tool's interior rows only.

Second half: the generic dump still says **`observed`** — the variant name spec §3.1 explicitly
retired ("'observed' is not a contrast with `Complete`… once the enum says *witness* the word
carries nothing"), as do the report fields `rows_observed` / `reads_observed`
(`ng_generic_loci_dump.rs:119, :121`) and the header they print. D4 already owns "the label
drift across the three dumps"; this adds that the fourth spelling uses a word this very spec
removed.

*What would stop it:* one body. Either `#[path = "shared/witness_label.rs"] mod witness_label;`
in each example (three declarations, one implementation — the compiler then cannot let them
disagree), or a method beside the type, which is defensible because the label is a pure
derivation from `ReadWitness` + `LocusLen` and four binaries want it. The `#[path]` form keeps
presentation out of the library and is the smaller change.

### F7 — **Minor.** `WitnessedLocusPositions` guarantees non-empty and offers no accessor that says so

*File:* `src/ng/locus_generation/witness.rs:136-169`.

The type's central promise is "never empty" (`:69-70`, `:104`, and `an_empty_input_or_an_empty_run_is_rejected_rather_than_dropped`).
Every accessor nevertheless hands back an `Option`-shaped answer:

- `runs()` is an iterator, so `first`/`last` are `Option` — D5 needed both and paid two
  `expect("a witnessed set is never empty")` per call site (verified above);
- `is_flush_left` / `is_flush_right` hide the same thing with `is_some_and` (`:156-169`), which
  *silently answers `false`* for a set that cannot exist. A set that could be empty would be
  reported "not flush left and not flush right" rather than caught.

This is an invariant with no home in the read API: it is enforced at construction and then
re-asserted, or quietly discarded, at every consumer.

*Fix:* total accessors on the type that owns the invariant —

```rust
pub fn first_run(&self) -> (u16, u16) { self.0[0] }                      // never empty
pub fn last_run(&self)  -> (u16, u16) { self.0[self.0.len() - 1] }
pub fn span(&self) -> u32 { u32::from(self.last_run().1 - self.first_run().0) }
```

The flush predicates become one-liners over them and stop carrying a branch for the impossible
case; D5's hole count becomes `span() - positions_covered()` with no `expect` at all.

### F8 — **Nit.** `witness.rs`'s module doc names the wrong importer

`src/ng/locus_generation/witness.rs:17` — "`canonicalise_runs` is crate-internal and is imported
by path, by `pileup::open_record`". Since `5790914` it is imported by `pileup::witnessed_ref`
(`witnessed_ref.rs:15`); `open_record` does not import it at all.

### F9 — **Nit.** `#[cfg(test)]` is a shorter and *stronger* guard than the two `cfg_attr(not(test), expect(dead_code, …))`

*File:* `src/ng/locus_generation/pileup/witnessed_ref.rs:64-70`, `:142-148`.

Both reasons are accurate and the guard fires (quoted above). Two caveats:

1. Under `cargo test --lib` the attribute is `cfg`-ed away entirely, so the guard exists only in
   the clippy step. A maintainer who wires `positions_covered` into the fold and runs the tests
   sees nothing.
2. Replacing both with plain `#[cfg(test)]` is strictly better. Verified: `cargo clippy --lib
   --all-features -- -D warnings` and `cargo clippy --all-targets --all-features -- -D warnings`
   are both clean, `ng::locus_generation` is **304 passed**, 12 lines of attribute go away, and
   the guard is upgraded from a lint to `error[E0599]: no method named positions_covered` — a
   hard error in every profile, including `cargo test`.

The `expect(dead_code)` form buys a written reason at the definition, which is worth something;
a one-line doc comment buys the same and the `#[cfg(test)]` keeps the enforcement.

### F10 — **Minor (flag for D3).** "Flush at both borders" does not mean "pinned", and nothing in the type says so

Two distinct witnesses satisfy `is_flush_left() && is_flush_right(len)` without the read having
witnessed the whole locus from both sides:

- **The STR saturating case.** `from_left(reach ≥ len, len)` builds `(0, len)` — flush both ways,
  from a read that anchored **one** border and ran out of read (`ssr.rs:863-883`; pinned by
  `ssr.rs:1372`). `is_flush_right`'s `>=` is deliberate and correctly defended
  (`witness.rs:160-169`, `a_run_reaching_past_the_locus_is_flush_right`), but the consequence is
  that flushness cannot distinguish "reached the border" from "reach exceeded the locus".
- **The holed case.** `witness.rs:578-583` asserts, correctly, that a read blind only in the
  middle "is constrained from both sides".

Spec §3.1 makes the flush predicates "the entire surviving representation of prefix-versus-suffix",
and step 7's censored likelihood is their consumer. A consumer reading them as *fully
constrained* is wrong in the first case and right in the second, and the type offers no way to
tell them apart. `Complete` is the only honest "pinned" test, and `witness_of` correctly decides
it on `positions_covered()` rather than on the outermost edges (`open_record.rs:231-237`) —
which is exactly the distinction the predicates cannot carry.

*Fix (documentation, and it is where departure 3 belongs):* say on `ReadWitness` that flushness
is a statement about **run placement** and never about completeness, that `Complete` is the only
"the length is pinned" test, and record the STR saturating case as the reason arch §1.1's
"return `Complete`" contract stays unimplemented. That is the sentence D3 needs in order to
decide rather than re-derive.

---

## What else I checked and found sound

- **The `RecordWitness` → `RecordWitnessCounts` boy-scout rename** the arch asked for is done
  (`open_record.rs:112`), with a doc that says why it is not a witness.
- **`ObservationRow` → `KeyedObservation`** is done (`open_record.rs:264`), and the dump tools'
  own `ObservationRow` correctly kept the word — spec §6 exempts the TSV line explicitly.
- **`canonicalise_runs`'s visibility is minimal.** `pub(in crate::ng::locus_generation)` is the
  narrowest scope that reaches `pileup::witnessed_ref`, which is a descendant. It cannot be
  tightened, and the "one normaliser, two axes" argument (`witness.rs:34-38`) is the right call —
  the alternative is two copies of a merge rule that can drift while both compile, which is the
  Milestone A finding it cites.
- **`witness_of` clamps per run, not per enclosing extent** (`open_record.rs:208-220`), and
  decides `Complete` on total coverage rather than on the outermost edges (`:231-237`). Both are
  the milestone's point and both are pinned by tests that fail under the enclosing reading.
- **`num_obs_along_locus` iterates runs and keeps its clamp** (`mod.rs:110-125`), with the
  clamp's reason stated correctly — the bound genuinely is not expressible on the type, because
  a witness cannot know its locus.
- **The invariant tests in `witness.rs` are where arch §6 asks for them**, the canonicality
  property is a proptest rather than a fixture, and it asserts the *position set* and not merely
  the canonical shape — which is what makes it catch a normaliser that loses positions.
- **The generic dump's per-run assertions are per run** (`ng_generic_loci_dump.rs:193-207`), so a
  witness whose enclosing extent fits but whose interior run overruns is caught.
- **`ReadWitness` stopping being `Copy`** was absorbed without a `clone()` on the hot path: the
  fold's per-read state holds a `WitnessedRefPositions` (a different type) and the locus-axis
  witness is built once, at `finalise`.

---

## Suggested order

1. **F5** (Major, mechanical, and the one that has already bitten once) — the counts method.
2. **F1** (Major) — drop `apply_events_into`'s `bool`, or at minimum `#[must_use]` it.
3. **F7** then **F2** — the two API-surface fixes D4/D5 will otherwise work around.
4. **F3**, **F4** — the two `witness_of` fixes; both are small and both are already compiled.
5. **F6** — with D4, which owns it.
6. **F9**, **F8**, **F10** — boy-scout, and F10's sentence should land before D3 decides.
