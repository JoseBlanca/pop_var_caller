# ng read likelihoods — F2: Model A, and where a slip is allowed to land

*Implementation report, 2026-08-25. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step F2 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone F, on
top of `3ca95b35`.*

## 1. What it is

`StutterSubstitutionEmission` — the default STR emission, and the first thing on this path that
produces a number. Two factors:

- **how likely a length change is**, from the distribution Milestone E built, for this
  candidate's stratum;
- **how likely the letters are**, once the candidate has been resized to the read's length —
  composed from the alignment module's `FlatEmission` under the fitted per-stratum rate, never
  re-implemented.

## 2. The correction that had to come first

F1's `SsrCandidate.bases` was documented as **the whole locus as a carrier has it, flanks
included** — carried over from the generic path, where that is exactly right. **It is wrong here**,
and building F2 on it would have compared flanks against flanks and put the tract at the wrong
offset for every resize.

On the STR path the locus *is* the tract. The generator slices exactly the tract into an
observation's bases (`region_seq[tract]`,
[`locus_generation/ssr.rs`](../../../../src/ng/locus_generation/ssr.rs)), the flanks are shared by
every candidate by construction, and comparing two copies of the same context would add nothing.
The doc now says so, and says why the two paths differ.

## 3. Where a slip may land, and the rule that had to be made consistent

In a pure tract, adding a repeat anywhere gives the same bytes. In an **interrupted** tract the
placements differ, so the model enumerates them and averages with equal weight — **whole-repeat
slips only**; a part-repeat change is resized at the tract's end in a single placement. That is
production's split, stated rather than left to be guessed.

**E2 and F2 disagreed about the deepest contraction, and the disagreement was real.** E2's
`unreachable_mass` reports a total contraction as *unreachable* — a read must still show a repeat,
which is the reading that reproduces the two sizes spec §4.2 states. Production's placement
enumeration allows a run to give up everything it holds. So on a pure tract the report would have
called a contraction unreachable while the scoring placed mass on it.

**Made consistent in the scoring**: a placement that would leave the tract with no repeats at all
is rejected. On a pure tract the two halves now agree exactly.

**On an interrupted tract they can still differ, in a stated direction.** A slip lands in *one*
run, so a tract of two runs of two repeats cannot lose three even though it holds four — while
`unreachable_mass` sees only the total and counts that contraction as reachable. **The report
therefore understates the loss on interrupted tracts.** Closing it needs the run structure to
reach the distribution, which takes only a period and a repeat count today. Recorded as an open
question rather than papered over (§6).

## 4. What the tests pin

The three spec §12 property tests, and four more:

| test | what it pins |
|---|---|
| §12.1 — a read identical to its candidate | scores at least `same_length_share × (1 − ε)^length` and within 5% of it |
| §12.2 — direction and size ordered as fitted | a repeat short outscores a repeat long; one repeat outscores two |
| §12.3 — a whole repeat beats a stray base | **under the corrected condition** — a comparison of *products*, asserted to hold for the fixture before the consequence is asserted, so it survives the one-step shares being untied |
| placements | an interrupted tract gives two distinct sequences for a one-repeat expansion; a pure one gives a single placement |
| unreachable | a contraction no single run can absorb scores **exactly zero**, and a pure tract may not vanish |
| part-repeat resize | lengthening continues the motif's tiling from the phase the tract ended on |
| the two factors | the score is the length factor times the letter factor, exactly, across three rates — and one mismatching base costs exactly one `ε/3` in place of one `1 − ε` |

**One fact worth recording, because a test was written against the wrong version of it twice.** On
a **pure** tract the unreachable branch cannot be reached through `emission` at all: the deepest
contraction an observation can ask for is the one that leaves it empty, and a single run can give
up everything it holds. It takes an *interrupted* tract to reach that branch — which is what the
test now uses.

## 5. Validation

| command | result |
|---|---|
| `./scripts/dev.sh cargo test` — library target | **4,390 passed, 0 failed, 14 ignored** (4,383 at F1) |
| `clippy --lib --all-features --tests -- -D warnings` | exit 0 |
| `cargo fmt --check` | exit 0 |

**11 tests** in `ng::calling::likelihood::ssr_emission`, against 4 at F1's close.

`censored_emission` is `unimplemented!()` — deliberately, and not a placeholder that scores. A read
that ran out is not a shorter complete observation, and a stub returning something would let a
truncated read out-discriminate a whole one. G1 builds it.

## 6. Deferred for the owner

Carried forward, plus one new and one sharpened:

1. **`unreachable_mass` understates the loss on interrupted tracts** (§3). It sees a period and a
   total repeat count; reachability is per run. New, and the sharpest of these because it is a
   number the row uses to keep candidates comparable.
2. **Spec §4.2's prose and its figures disagree** about whether a tract can lose its last repeat.
   This branch follows the figures, and F2 now makes the scoring follow them too — so the choice
   is load-bearing in two places rather than one.
3. **Spec §4.2 calls the part-repeat cutoff "10 base pairs"** where it counts a compressed rank.
4. **`arch/read_likelihoods.md` §4.2** sketches `stutter_rates_for` beside the distribution and
   names the context's field `truncated_mass_lost`; both moved.
5. **Spec §12's fourth test asks periods 1 to 6; the tripwire runs 2 to 6.**
6. **`Provenance`'s ladder puts *supplied* below *borrowed***, and `weaker_of` follows it.
