# ng cohort merge — E3: applying the review's fixes

*Fix-application report, 2026-08-18, against
[the E3 review](../reviews/ng_cohort_merge_e3_2026-08-18.md). Every finding was applied.*

## 1. The two that mattered, and both are about the round being invisible

**The round's size was pinned by an inequality a doubled round also passes.** The memory test
asserted only that eight builders hold more than one; the reviewer doubled the round and the
whole suite stayed green, with the held-record table moving from 2/4/4/12/29 to 4/4/12/29/62.
The test now asserts the exact table. Three round-size mutations — forced to one region, to
twice the count in flight, and to the whole analysed region — all die on it, and it is the only
test any of them fails.

**Nothing could tell a parallel merge from a sequential one.** `.par_iter()` → `.iter()` needed
no other change and left the suite green. The round's builders now go through `in_region_order`,
whose bound is `IndexedParallelIterator` — so `.iter()` is a compile error (verified: *the trait
bound `Map<Iter<'_, GenomeRegion>, …>: IndexedParallelIterator` is not satisfied*), and so is an
adaptor that gives up the index. **Whether the builders genuinely occupy several threads is
still untested**, and the module doc says so rather than leaving the next reader to look for the
test.

## 2. The rest

- **`CohortLocusBuilderRegionsInFlight`** replaces the bare `NonZeroUsize`, joining the module's
  three other run parameters. Its default is a rule rather than a constant —
  `one_per_worker_thread()` — because what it should be depends on the machine's cores and on
  how much memory the cohort's width leaves, neither knowable where the other three defaults are
  written. A test pins the rule, since a default ignoring rayon's pool would pass every other
  test in the module.
- **`refuse_malformed_analysed_regions` moved to `mod.rs`**, where all three drivers can find
  the contract they share, instead of `parallel.rs` importing it from the file named for the
  other arrangement.
- **The shared fixtures grew**: `source_of`, `width`, `in_flight`, `render`,
  `refuse_any_difference` and the 600-base three-sample layout the byte-identity claims rest on.
  The two files whose outputs are compared now read the same fixture — the argument that moved
  `member` in the first draft, applied to the rest.
- **`RegionOutcome` is destructured** where the driver consumes it, as `organise.rs` does at its
  own two consumers.
- **`S: Sync` is documented** as what sharing the cache costs and as a constraint the run's
  future `ObservationSource` has to honour — a source built on `Rc` or `RefCell` cannot be
  merged in parallel.
- **A panicking builder is documented and tested.** Two overlapping records of one sample pass
  the cache, which checks only that starts do not go backwards, and trip `build_region`'s
  disjointness assertion on a rayon worker. Rayon joins every builder before re-raising, so the
  cache is unborrowed when the caller regains it; it is not rewound, and `Organiser::finish`
  never runs.
- **Four tests added**: two intervals on the same contig with a gap and an interval shorter than
  one region; a randomised sweep against the oracle with the count in flight drawn too; the
  panicking builder; and the overlap half of the shared guard.
- **Nits**: `Vec::with_capacity` on a caller-supplied count, `as usize` narrowing a `u64`, three
  names, and both module docs still calling this step "still to come".

## 3. Documentation corrections

**One number was wrong.** "at most 15 idle builders once" at a 600-base interval on 20-base
regions is **2**: 600 ÷ 20 is 30 regions, so at 16 in flight the rounds are 16 and 14. Fifteen
is the bound `regions_in_flight − 1` over all interval lengths, quoted against an interval that
determines the answer exactly. The doc now gives the rule, the worked case and the bound as
three separate statements.

Also corrected: `n` was used undefined in the spec amendment (now `builders × k`, with `k` named
as the cohort's sample count); "three of them want it" pointed the pronoun at test modules, of
which there are two, where three *drivers* want the fixture; "no test reaches" the driver's two
end assertions became "no test makes either fire", since every merge evaluates them; and
`organise.rs`'s note said E3 gives the eviction choice to the organiser, where it is the
driver's.

## 4. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `224 passed; 0 failed` (219 before these fixes).
- `cargo test --lib` — `3847 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out;
  finished in 743.74s` (3,832 before this step).
- **`tmp/mutate_e3.sh` — 13 mutations, 11 killed.** The two that survive are the driver's own
  end assertions, documented as safety nets: every merge evaluates them and no fixture makes
  either fire, because both can only trip if the organiser displaced something and nothing a
  builder produces can. A fourteenth, `.par_iter()` → `.iter()`, is now a compile error rather
  than a mutation.

## 5. Follow-ups

- **The round's tail is unmeasured**, and spec §6.2 says the two alternatives — an `RwLock`, or
  windows handed out as owned copies — should not be reached for before it is.
- **That the builders occupy several threads is untested.** It needs a hook the driver does not
  have; the type-level guard is what stands in for it.
- **Arch §4's `Organiser::cache()`**, owed until the run's own source and error types land.
