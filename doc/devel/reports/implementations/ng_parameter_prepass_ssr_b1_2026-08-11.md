# ng step 4, the STR path — B1: one locus's reads, laid out across the offset buckets

*Implementation report, 2026-08-11. Step B1 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — five agents, 33 mutations, 6 real survivors. Design authority:
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.2 and
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1, §5.*

## What the step is

`LocusShape`: how one locus's reads fell across the nine whole-repeat offset buckets, plus how
many showed a length that is not a whole number of motif copies. It is what the sparse
per-stratum table of Milestone B2 is keyed on, so two loci that looked alike become one entry
with a count of two.

Nothing stores the depth — it is the buckets plus the guard — and it is exact rather than binned
because `MAX_LOCUS_READS` bounds it at twelve values, so a ladder over it would save nothing.

**The constructor takes a wider integer than the shape stores.** Counts arrive as `u32`, the
width the cap itself is written in and the width a caller tallying a deep locus holds; the
storage is `u8`. A constructor as narrow as the storage would make `try_new(tally as u8)` the
natural call, under which 260 reads arrive as 4 and validate as an ordinary shallow locus — a
depth wrong 65-fold with nothing to notice it. It also refuses an empty shape: that would enter
the table as a locus, counting towards `MIN_LOCI_TO_FIT` and towards every "loci behind this
fit" number while contributing a likelihood of exactly one to every candidate.

## Recorded deviations from the architecture

1. **A checked `try_new` where arch §2.2 sketches public fields.** The plan's own acceptance test
   ("a shape whose counts exceed the cap cannot be built") needs a constructor that can refuse, so
   the fields are private and the rejection is a `DomainError` — one new variant,
   `SsrLocusShapeReads`, following the enum's convention of naming the quantity that was wrong.
2. **The empty shape is refused as well as the over-cap one.** The plan names only the cap. The
   reason is above; recorded because it puts a second condition behind one error variant, exactly
   as `SsrPeriod` covers zero and over-range in one.
3. **`whole_repeat_depth` for what the design calls the *scored depth*** (arch §2.4). Named as the
   exact complement of the guard, in vocabulary the module already owns; the design's term is
   given in the doc comment so the link survives.
4. **The cap travels inside the error rather than being interpolated into its message.** The
   message lives in `types.rs`, the shared vocabulary, and `MAX_LOCUS_READS` is the STR path's own
   constant — naming it there would point the wrong way up the module tree. It is also a value
   arch §7 expects to move, so an error carrying the number in force says what a run rejected
   against.

## What the review changed

**Blocker — no fixture put a read into either saturating end bucket's arithmetic.** Four mutants
survived all 43 tests: `whole_repeat_depth` dropping either end bucket, `reads_in` folding the
top end into its neighbour, and `try_new` discarding the top end. Each silently under-counts
precisely the loci carrying a long or short allele — the ends hold 0.89% of reads across GRCh38's
typed tracts. Every offset in every fixture was one of `0, ±1, ±2, -4`, and the one fixture that
reached bucket 0 was used only for sorting. Two new tests close it, and I re-ran one of the
mutants against them: dropping the top bucket from `whole_repeat_depth` now fails
`every_bucket_round_trips_through_both_accessors` and
`reads_that_saturated_into_an_end_bucket_are_stored_read_back_and_counted`.

**Blocker — the key-identity test exercised 2 of the 10 places a shape can differ.** The reviewer
replaced the derives with `Eq`, `Ord` and `Hash` keyed on the first eight buckets — ignoring the
top one entirely — and the whole module stayed green. This type is the `BTreeMap` key B2 builds
on, so such a key merges two loci that observed different things into one entry, changing both
the per-shape likelihood and every loci-behind-the-fit number, with nothing going red.
`shapes_differing_in_any_single_bucket_are_distinct_keys` sweeps all ten places through both the
hash and the `BTreeMap`; I reproduced the ignore-the-top-bucket key and confirmed it now fails.

**Major — the accepted range's lower edge was untested.** The shallowest fixture in the file held
three reads, so widening the guard to reject one- and two-read shapes left everything green. A
locus contributing one read is the common shallow case, and rejecting it is silent: the caller's
answer to a rejected shape is to enter nothing, so those loci would vanish from the stratum's
count. Reproduced and killed by `a_shape_of_a_single_read_is_accepted`.

**Minor — one assertion could not fail.** `counts().iter().sum() == scored_depth()` was
character-for-character the body of `scored_depth`, so both sides moved together under every
mutation; both end-bucket mutants passed it. It now compares against the literal array the
fixture was built from.

**Naming, five Minors applied.** `not_whole_repeat()` → `reads_not_whole_repeat()` (a modifier
with the noun missing, where the crate attaches it everywhere else); `counts()` →
`reads_by_bucket()` (two names for one quantity, and a bare plural naming the container rather
than what is in it); `scored_depth()` → `whole_repeat_depth()`; the error message no longer names
one object two ways in one sentence; and the module doc now separates *locus shape* from *entry*
rather than using both for the key.

**Idiomatic, two Minors applied.** All six accessors take `self` rather than `&self` — every
sibling value type in the module does, and the receiver mode was forced only by `counts()`
returning a reference. And every count read singly comes back as `u32` rather than `u8`, which
removes the widening the type was doing inside its own body; the array alone stays in its stored
width, because that is what the fit converts element by element.

**Four wrong numbers, all mine, all about my own code.** The const assertion's doc claimed a cap
above 255 would let a bucket "wrap silently, since this repo's release profile leaves
`overflow-checks` off". The reviewer mirrored `try_new` with the cap at 300 and ran it under both
settings: it **panics** at the narrowing, identically in both builds — `overflow-checks` governs
arithmetic operators, not conversions. The silent-truncation story is real but belongs to the
*signature* argument, where the cast sits on the caller's side. So the assertion buys a compile
error instead of a crash, not instead of a wrong entry, and the comment now says that. Also
"a hundred-fold wrong depth" for 260 narrowing to 4 (it is 65-fold), and "nine `u32::MAX`s
overflow a `u32` sum by a factor of nine" where the fixture offers ten values.

**Declined, with reasons.** Replacing the two `expect`s with `map_err(…)?`: the reviewer proved
both unreachable, so the `Err` arm would be a branch no test can reach, and it would report a
caller error for what could only be an internal bug. They carry the repo's `// PANIC-FREE:`
marker instead, in the house form. Moving `SsrLocusShapeReads` out of `DomainError` into a
module-local error: keeping it there is what lets `SsrEstimationError::Domain` wrap it with the
sample, read group, stratum and fit for free.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib parameter_estimation::ssr` | **47 passed** (33 before this step) |
| `cargo test --lib --bins --tests --all-features` | **3,424 passed, 0 failed, 10 ignored** |

Counted rather than recalled: `grep -c '#\[test\]'` on `stratum_table.rs` gives **14**, and the
suite moved 3,410 → 3,424 with this step.

**Two gates are red on this branch and neither is this step's**: `cargo clippy --all-targets`
fails in four `examples/` files, and `cargo doc` reports 13 unresolved intra-doc links. Both
reproduce with this branch's changes reverted, which is why validation runs `--lib --bins
--tests`.

## A process note on the review fan-out

Three of the five agents reported that they were **not** given an isolated worktree — their
working directory was the main checkout — and each declined to run the prescribed
`git checkout --detach`, which would have destroyed the live tree. Two reviewed read-only and one
built its own worktree by hand. One agent nonetheless edited the main checkout's
`stratum_table.rs` (a cosmetic refactor of the two narrowing calls); the tree was restored from
the reviewed patch before any fix was applied, and `diff <(git diff) step.patch` confirmed it
byte for byte. The `rust-code-review` skill's step-0 instruction assumes the harness supplies the
worktree; when it does not, the instruction is destructive.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-b1/` — five per-category files (reliability, errors, naming,
idiomatic+smells, numbers), the reviewed patch, and the reviewers' probes.
