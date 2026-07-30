# Code review — ng generic locus generator, Milestone D (D1–D3)

**Date:** 2026-07-30 · **Scope:** `e626850..31be8ee` on `ng-generic` (D1 `6993704`, D2 `7bfcd8a`,
D3 `31be8ee`) · **Fixes:** this commit ·
**Impl report:** [ng_locus_generation_pileup_generator_d_2026-07-29.md](../implementations/ng_locus_generation_pileup_generator_d_2026-07-29.md)

## The fan-out did not happen, and that is the first thing to say about this review

Five category agents were dispatched, each meant to run in its own git worktree, as the previous
three milestones did. **Seven of eight launches were killed by API 529 overload** — the five
originals and two retries, over about forty minutes. One agent completed: **reliability**, plus
the test-challenge pass, re-launched on a different model tier as a test of whether the overload
was tier-specific. It was.

So Part 1 below is a **one-category review**. Four checklists did not run: `errors`/`defaults`,
`idiomatic`/`smells`/`unsafe_concurrency`, `module_structure`/`naming`/`refactor_safety`, and
`extras` (spec conformance and the measurements). Milestone D should not be treated as fully
reviewed on the strength of Part 1 alone.

**They were run on 2026-07-30 and are Part 2 of this document**, which also diagnoses the
isolation failure below. Part 2 found **11 Majors**, 6 of them a mutation or an input nothing
catches — so the assumption that the four missing categories held findings was correct, and the
margin was not small.

**One process finding, from the agent that did run: the worktree isolation did not hold.** It
worked in `/Users/jose/devel/pop_var_caller-ng-generic` — the author's own worktree — applying and
reverting ten mutations there, and wrote its report into that tree's `tmp/`. No damage resulted
(the tree is clean at `31be8ee` and the suite re-verified at 2,724), but a mutation-heavy agent
sharing the author's tree is exactly the collision the per-agent worktree rule exists to prevent,
and with several such agents it would have been unrecoverable. Worth diagnosing before the next
fan-out.

## Verdict: 2 Majors, both "a mutation nothing catches", one of them Milestone D's own

The reliability agent re-ran all seven mutations the impl report claims are caught and confirmed
every one, then tried three of its own. Two of those three fail **no test** in either suite.

### Major 1 — `flush_all`'s `ever_contributed` guard was untested. **Fixed.**

`src/ng/locus_generation/pileup/active_read_set.rs`

`reads_silent_over_footprint` is fed by the active set's **two** exits and D2 pinned only one.
`expire_passed` — a read the walker has passed — is covered by
`a_read_silent_at_every_position_is_counted_rather_than_lost`. `flush_all` is the other, and on
the generic path it is **not an edge case**: a region walk stops at `region.end` while the reads
reaching into the halo are still active, so *every bounded walk ends by flushing reads that never
expired*. Deleting the guard there — counting every flushed read as silent — left the **whole
2,724-test suite green**, which the author reproduced before fixing.

The counter would then have over-reported on every region of every real run. That is the failure
mode this milestone exists to eliminate: a number nobody can see being wrong.

**Fixed** by `a_read_still_active_when_the_walk_stops_is_counted_by_what_it_contributed`
(`generator.rs`), and the fix needed two goes:

- **The first draft could not fail either.** With one silent read and one contributing read, the
  correct guard and a guard *inverted* to count the contributors both total 1. Two contributing
  reads make the two answers 1 and 2.
- Mutation-verified in **both** directions: guard deleted → 3 against 1; guard inverted → 2
  against 1; guard restored → green.

That is the thirteenth test on this branch found unable to fail, and the second in two milestones
to be a test written *for a review finding*.

### Major 2 — `refold_live_reads`' contributor-skip has no regression test. **Carried, with the reviewer's sketch.**

`src/ng/locus_generation/pileup/open_record.rs`

Deleting `if contributors.iter().any(|c| c.read_id == read_id) { continue; }` changes nothing
observable in either suite (202 lib + 10 dump), and makes the `contributors` parameter unused —
which confirms the skip is that parameter's only remaining use.

**Carried rather than fixed, for three reasons.** It is **A3's** code, not Milestone D's, and the
function's own doc comment already records the gap deliberately: *"Unpinned, and deliberately so…
Mutating the skip away leaves the whole suite green… Do not read the absence of a failing test
here as the absence of a reason."* The reviewer looked for a correctness counterexample and did
not find one — the carried `contribution` makes a double re-place idempotent — so the exposure is
a future edit breaking that invariant silently, plus a possible change in allele **creation
order**, which is observable in the output. And the test wants the record's internal allele table
at a chosen walker position, which is a fixture shape none of the existing ones have.

The reviewer's sketch, for whoever writes it:
`refold_live_reads_skips_a_read_that_is_also_the_widening_contributor` — build a record where read
`r` folds at one position, then at a later position both `r` has an event (so it is a contributor
there) and a second read anchors a deletion that widens the record; assert the allele order
matches what the fold loop alone would produce.

### Not filed: `apply_events_into`'s `run_end.max(event_end)`

The reviewer also mutated `witnessed = Some((run_start, run_end.max(event_end)))` to drop the
`.max()`; no test fails. Filed as a note rather than a finding, with the reasoning: within one
read's CIGAR-ordered event stream `event_end` is non-decreasing, so the two forms agree on every
input the debug-asserted precondition admits, and the `.max()` reads as defence for a case that
precondition already rules out. **Carried as a note to check against `decompose.rs`'s
events-overlapping construction** before anyone relies on it, since the reviewer could not build a
failing input either way.

## What the review confirmed

All seven prescribed mutations fail a test, and **one is caught more widely than the impl report
claimed**: `coverage_of` always returning `Complete` fails 9 lib tests plus **4** dump fixtures,
where the report said two. Corrected in the impl report.

The agent also audited the permanent anchor and its fixture, and could not break either:

- **The floor assertions bound something real.** `anchored > total * 5` and
  `anchored_multi_base > total` stop the fixture degenerating into single-base loci where the
  equality is vacuous; `widens > 0` stops it losing the one property that makes
  `generate_uniform_events` worth more than `generate_complete`. All three read accumulated
  per-run counts, not constants.
- **It could not construct an input where the fixture leaves a read stale**, and traced why by
  hand: every read on a contig is built from one shared template, so wherever one read's event
  triggers a widen every other live read has an event at that position and is a contributor there
  — not the "live but silent" read a stale fold needs. Both paths that can make a folded read stop
  being a contributor mid-record (the column cap, the mate-overlap collapse) operate on the same
  post-truncation list `refold_live_reads` is driven from, so the read is re-folded from its own
  cursor. The argument depends on no two reads on a contig having different starts or CIGARs,
  which is what the fixture-property assertion checks.

## One documentation contradiction. **Fixed.**

`parity.rs` — the anchor's doc said `generate_uniform_events` gives every read one event set "so
**no record widens at all**", and two paragraphs later said the opposite, correctly: "not by
stopping records widening, which it does not". The first was a leftover from the draft whose
`record_widen_events == 0` assertion had already failed at 7 widens; the test itself agrees with
the second (it asserts `widens > 0`). Corrected to say what is true — the shared event set removes
the *staleness*, not the widen — because the contradiction made the property impossible to verify
by reading.

## Validation after fixes

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` **2,725 passed** (2,724 + the new test); the example's 10 tests green;
`cargo doc --no-deps` still 12 pre-existing unresolved links.

---

# Checkpoint D's answers, applied

**Date:** 2026-07-30. The owner answered the three carried items: the `to_pileup_record`
conversion was left to my judgement, and the read-group grain is **settled — per-`@RG` stays**,
because the statistical models for both SNP and STR calling need the evidence per library. No
measurement was needed for that; the requirement is the reason.

## `to_pileup_record` is deleted

The 44 inherited tests now assert on `SampleLocusObservations` directly. Production's positional
idiom is carried across by four accessors on the locus type, test-only:

| production | ng |
|---|---|
| `pos` | `anchor()` |
| `ref_span()` | `footprint_len()` — the region's length, not a row's |
| `alleles[0].support.*`, `alleles[0].chain_ids` | `reference_row().*` |
| `alleles[0].seq` | `reference_bases` — a field, not a row |
| `alleles[1].*` | `first_alt_row().*` |
| `alleles.len()` | `observed_sequences.len()` |

**`reference_row` and `first_alt_row` panic rather than return an `Option`, and that is the point.**
The two ways ng's rows differ from production's buckets are exactly what a silent `None` would
hide: a locus no read matched the reference at has **no** reference row, where production created
the bucket regardless of support; and a locus whose reference-matching reads split by coverage or
read group has **several**, where production had one merged bucket. A test landing on either is
asking a question ng's type does not answer, so it fails loudly and names which case it hit.

**Three inherited tests changed meaning, and each was checking less than it appeared to.** All
three are the same input class: every read carries a non-reference base, so production's
`alleles[0]` was an empty bucket.

- `column_depth_cap_keeps_first_n_of_admission_order` asserted `alleles[0].support.num_obs == 0`
  ("no surviving read matched REF") and `alleles.len() == 3` ("REF + 2 surviving SNP buckets").
  ng emits **two** rows and no reference row. The new assertions say the same fact about the
  evidence without counting a bucket that held none.
- `paired_mates_with_overlapping_positions_share_chain_id` asserted
  `alleles[0].chain_ids.is_empty()` — "REF chain ids are dropped". **That held for free:** both
  mates carry `C` over an `A` reference, so `alleles[0]` was empty, and an empty bucket's id list
  is empty whatever the rule does. The chain-id rule is checked where it can fail — in
  `only_the_reads_that_departed_from_the_reference_carry_a_chain_id` and the dump fixture of the
  same name, both of which put a genuinely reference-matching read beside a departing one.

That is a fourth assertion this milestone found unable to fail, and it was found by the
conversion rather than by mutation: moving to a type that cannot represent an empty bucket made
the vacuity visible.

**Validation:** fmt clean; clippy `--all-targets --all-features -D warnings` clean;
`cargo test --lib` **2,725 passed**; the example's 10 tests green; `cargo doc --no-deps` still 12
pre-existing unresolved links. The `open_record.rs` fixture that used the projection now reads
ng's type; the two doc comments that cited `to_pileup_record`'s losses as live facts now say it is
deleted.

---

# Part 2 — the four categories that did not run, run

**Date:** 2026-07-30 · **Scope:** `e626850..af967d5` (the seven commits above, plus `dad9baf`
and the two doc commits) · **Findings:**
[`tmp/review_2026-07-30_ng-generic-milestone-d/`](../../../../tmp/review_2026-07-30_ng-generic-milestone-d/)
— `_context.md`, `errors_defaults.md`, `idiomatic_smells.md`, `structure_naming.md`, `extras.md`

Four agents, one worktree each, all four completed. **11 Majors**, of which 6 are a mutation or
an input nothing catches. Every Major below was re-verified in the author's tree against the
source before being filed; the two that turned on a claim I could not confirm by reading are
marked where that matters.

## Why the isolation failed last time, diagnosed before dispatching

`isolation: "worktree"` does hand each agent its own tree. It checks it out at **`main`**.
Probed: HEAD came back `2c510a8`, which is exactly `git rev-parse main`.

**That tree does not contain the code under review.**
`src/ng/locus_generation/pileup/` holds 11 files there against 13 on the branch, and the two
missing ones are `generator.rs` and `mock_reference.rs` — the generator itself. Every A/B/C/D
change to the other eleven is absent too.

So the reliability agent was not ignoring its worktree. Its report's first line reads
`Scope: … (worktree /Users/jose/devel/pop_var_caller-ng-generic…)` because the prompt handed it
that path, and that was the only checkout on the machine where the files it was told to mutate
existed. Five such agents would have converged on the author's tree for the same reason, and
that would have been unrecoverable rather than merely untidy.

A second mechanism stacks on top: an **unchanged** worktree is auto-cleaned, and an agent
resumed after that lands in a real checkout. Observed live — the probe agent, resumed, reported
its cwd as `pop_var_caller-ng-generic` and stopped rather than work there.

**The fix, verified before the fan-out**, is a mandatory step 0 in every prompt:
`git checkout --detach af967d5` in the agent's own tree. Tested on a worktree created off `main`
exactly as the harness makes them: `generator.rs` absent before, present after, tree clean, no
fetch needed since it is one object store. All four agents confirmed step 0 in their reports.

**And one recorded claim is wrong.** Milestone C's `intent_and_cost.md` filed a cross-category
note saying `scripts/dev.sh` "cannot build a worktree under `.claude/worktrees/` — it silently
builds the main checkout instead, then fails with `Exec format error`". That is not a property of
`dev.sh`: `PROJECT_DIR` is derived from the script's own location, so an agent's own copy mounts
and builds the agent's own tree. Measured — from a worktree under `.claude/worktrees/`, its own
`dev.sh` compiles `pop_var_caller v0.1.0 (…/worktrees/probe-build)` in 31 s, runs 15 tests green,
and puts its artefacts in that worktree's `target-container` (1.3 GB). What the C note describes
is invoking `/Users/jose/devel/pop_var_caller/scripts/dev.sh`, which is what C's own context file
told agents to do in the words "by absolute path". Both were corrected in this fan-out's prompts.

## The measurements reproduce exactly, and the harness is deterministic

The decisive question for a differential shipped as a deliverable, answered: **yes.**

Every synthetic number in the D1 report came back identical — the anchor's 216,203 of 216,203,
all six class counts, the deliverable's 2,787 reads / 1,484 loci / 8,239 bases, and the
5,000-case soak's 3,262,582 loci with class 6 at 3,074. Five runs across three profiles
(`soak` ×2, `soak` at 5,000 ×2, plain `debug`), with the report lines byte-identical between
repeats.

The seeding turns out to be **total** rather than mostly. The four seeds are a `const` and every
generator takes a `&mut SplitMix64` with no external state; `ng_walk_in_groups` deals read groups
by `index % groups` over the case's read order, so there is no RNG in that path at all; and
`AHashMap` order — the one non-seeded input — is neutralised upstream by taking reads in
`read_id` order.

D3's real-data numbers were **not** re-derived: no BAM exists in an agent's worktree.

## A. Three counters, and four assertions that cannot fail

This is the branch's defining defect and it accounts for 3 of the 11 Majors.

### Major — `generator.rs:678`: the declined-read tally's `std::mem::take` has no test at all

The line above it already explains the trap: *"Taken, not read: the cell outlives every walk, so
a tally left in it would be folded again at the next region's end — the same shape as the shed
error below."* The guard is right, the reason is written down, and **nothing tests it.**
Replacing the `take` with a plain read leaves all 213 tests green; the reviewer's two-region
probe then reports **5** declined reads where 3 were declined.

`a_read_the_preparer_declines_is_counted_and_never_admitted` walks **one** region, which is
exactly the shape that cannot see it. Its sibling `reads_silent_over_footprint` *is* pinned,
because that one is folded in `fold_region_walk` and this one is not.

That makes **three** counters on this branch whose guard the whole suite could not see: the
allocator's `reset`/`summary` pair at Checkpoint C, `flush_all`'s `ever_contributed` in Part 1
above, and now this. All three are spec §8's trap, and in all three the author wrote the guard
*and* the comment explaining it. **The pattern is not that the trap is unknown — it is that a
guard against triangular summing is invisible to a single-region test, and nearly every test
here walks one region.**

### Major — `tests.rs:378`: the fourth assertion `dad9baf` lost, and the subscript is why

`deletion_record_does_not_double_count_ref_reads` binds
`let ref_allele = &anchor.observed_sequences[0];` and asserts `num_obs == 1` on it under the
message *"REF: 1 obs from r1 only"*. But `finalise` sorts rows by `bases`
(`open_record.rs:601-610`), and this locus holds `"C"` (r2's deletion) and `"CGTA"` (r1's
reference match). `"C"` sorts first. So `ref_allele` **is the deletion row** — the same row the
test then finds again by `bases == b"C"` and asserts `num_obs == 1` on a second time. The
reference row is never inspected, and re-creating the exact historical bug the test is named for
(r1's REF row reporting `num_obs = 4`) leaves it **green**.

**The conversion's own table has the right accessor for this site** — `alleles[0].support.*` →
`reference_row().*` — and this site did not use it. Production's `alleles[0]` was positional and
*was* the REF bucket; ng's `observed_sequences[0]` is the lexicographically first row. Same
subscript, different meaning, no compiler error. That is precisely the risk the B2 decision
recorded — *"67 hand-translated assertions are 67 chances to re-express a test slightly weaker
than it was"* — materialising at one of the 49 sites, and it is the **fourth** such assertion in
this milestone after the three the conversion itself reported.

### Major — `examples/ng_generic_loci_dump.rs:164-172`: the global §13 assertion never reads `bases`

`push_locus`' doc restates spec §13 verbatim — *"No row claims more locus positions than its
events account for — the consistency check between `bases` and `read_coverage`"* — and what is
asserted is `offset_in_locus + positions_covered <= footprint`, a bound of `read_coverage`
against the *region* that never mentions `bases`, and which skips `Complete` rows. Weakened to
`reach <= footprint + 1000`, **all ten dump fixtures stayed green.**

The spec clause as literally worded may be unimplementable — `positions_covered` is *derived
from* the events, so checking it against them is tautological, and §8's own trap says no
inequality relates `bases.len()` to the footprint in general. That makes this a **spec edit plus
a doc correction**, not a code fix: say what is checkable, and stop restating a stronger claim
above weaker code.

### Major — `parity.rs:1419, 1497-1529`: `float_only` is printed and unasserted

`float_only` is the count of loci that agree only within the `q_sum` tolerance, and its stated
purpose is that *"the tolerance is shown to be doing work rather than quietly matching nothing"*.
It has a floor of zero and no ceiling. Measured, by injecting a uniform **relative** `q_sum`
error into every emitted row:

| injected relative error | result |
|---|---|
| 5 × 10⁻¹⁰ | **all 12 tests green**; anchor `float_only` 103 → **215,659** of 216,203; census 863 → **253,143** of 256,974 |
| 2 × 10⁻⁹ | 3 tests fail |

So the discrimination threshold sits where it was designed to, and the tolerance is not vacuous —
but the census keeps its headline invariant while the quantity that *explains* the headline moves
by three orders of magnitude. The census already has this exact guard on the class-1 deliverable
("a ceiling as well as a floor"); the same reasoning applies verbatim.

**The tolerance itself is sound, and D3 already fixed the part that was wrong.**
`Q_SUM_TOLERANCE = 1e-9` is applied **relatively** — `|a − b| ≤ 1e-9 · max(1, |a|, |b|)` — so the
carried note asking for "a relative tolerance rather than an absolute one" is **already
discharged**. Loose enough at large `q_sum`: summing N `f64` terms accumulates ≈ N · 2.22 × 10⁻¹⁶,
which at the 414-observation locus D3 hit is 9.2 × 10⁻¹⁴, four orders inside the tolerance, and
does not reach it until N ≈ 4.5 × 10⁶ observations at one locus. Tight enough at small `q_sum`:
the `max(1.0)` floor makes it absolute at 1e-9 below |q_sum| = 1, and the smallest non-zero
per-read contribution is BQ = 1, i.e. ≈ 0.2303 — eight orders above the floor. A read *can*
contribute exactly 0.0, but moving it also moves `num_obs`, `num_fwd`, `mapq_sum` and
`mapq_sum_sq`, all compared exactly.

### On the running tally: it is unreconcilable, so this report stops keeping one

The reports disagree, and they cannot all be right: D1 §6 calls its find *"the eleventh"*, D2 §5
then says *"twelve … three of them in this milestone"* — which cannot follow from an eleventh
plus D2's own two — Part 1 above says *"the thirteenth"*, and the Checkpoint D section says
*"a fourth … in this milestone"* on a different basis again.

**Milestone D's members, named instead of numbered** — nine, of which four are new here:

1. `ng_emits_no_allele_bucket_without_support` — eviction order (D1 §6).
2. D2's chain-id dump fixture — the positional and per-read rules coincided on it.
3. D2's region-boundary dump fixture — the halo was not exercised.
4. `paired_mates_with_overlapping_positions_share_chain_id`'s empty-bucket assertion — found by
   the type conversion, not by mutation.
5. `flush_all`'s `ever_contributed` guard (Part 1) — **and its first fix draft, which also could
   not fail.**
6. `deletion_record_does_not_double_count_ref_reads` — this review, §A above.
7. `generator.rs:678`'s `std::mem::take` — this review; no test at all, not a weak one.
8. `push_locus`' global assertion — this review.
9. `float_only` — this review.

A branch-wide total is not worth reconstructing and citing a wrong one is worse than citing none.
What is worth carrying is the **failure shapes**, which now number four: a fixture that masks the
mutation it was written for; an A-vs-B equivalence test with a bug in the shared arm; a
single-region test blind to a cross-region fold; and a subscript whose meaning changed under a
type conversion with no compiler error.

## B. The public surface says things that are not true

### Major — `mod.rs:155-161`: `run` being `pub` leaves half the config bypass open

The same file states the rule at `:80-84`: *"The vocabulary is bound `pub(crate)`, not `pub`: it
is an internal aid for the copies, not an ng-flavoured public alias for production's types…
which matters from plan 3 on, when ng's walker starts to diverge and two live paths to one name
would stop being harmless."* Then `:155-161` exports `WalkerError`, `PileupWalker`, `RunSummary`,
`run` and `DEFAULT_MAX_ACTIVE_READS` — production's walker vocabulary, under identical names, in
ng's **public** API.

`run` is the one with teeth. `to_walker_config` was made `pub(super)` to force every caller
through `PileupGeneratorConfig::check()` and its `MAX_RECORD_SPAN_CEILING`; with `run` still
`pub`, an external caller builds a `WalkerConfig` from production's `pub` `Default` and runs ng's
walker directly. The reviewer did exactly that from an external example and reached
`max_record_span = 983,025` — **15× the ceiling `check()` enforces.**

### Major — the `Arc<R>` cannot buy thread-safety, and the reason on file is the weaker one

Carried since D2 as *"the walk is single-threaded (arch §9), so `Rc` would say what is true"* —
an argument from intent. The structural fact is stronger: `PileupGenerator` holds
`preparation: Rc<RefCell<ReadPreparation<P>>>` (`generator.rs:538`), so the generator is
**already `!Send` and `!Sync`** whatever the reference accessor is. No consumer can require
`Send + Sync`, so nothing can break. Made in the reviewer's tree: clippy clean under
`--all-targets --all-features -D warnings`, 2,725 lib + 10 dump tests green.

**And the citation the old argument leaned on does not exist.**
`doc/devel/ng/arch/locus_generation_pileup.md` has §1–§5 and no §9. Four code sites cite
"arch §9" for single-threadedness — `locus_generation/mod.rs:632`,
`pileup/generator.rs:334`, `pileup/active_read_set.rs:55`,
`examples/ng_generic_loci_dump.rs:352`. The claim itself is real and lives in
`spec/locus_generation.md:676` ("Parallelism — deferred whole (§9)"); the document named is the
wrong one, in all four.

### The re-export audit, settled by four whole-tree builds

Your carried item, answered per name. **A bare-name grep cannot do this** — production has
identically-named `WalkerError`, `PileupWalker`, `RunSummary`, `run` and
`DEFAULT_MAX_ACTIVE_READS` — so each name was demoted in turn and the tree rebuilt. One method
note worth keeping: demoting `PileupGenerator` makes the *lib* target treat the module as dead,
and under `-D warnings` that cascade stops cargo ever reaching the example, so the answering
build is `cargo check --all-targets --all-features` **without** `-D warnings`.

| name | what broke when demoted | verdict |
|---|---|---|
| `PileupGenerator` | `E0603`, `examples/ng_generic_loci_dump.rs:50` | keep `pub` |
| `PileupGeneratorConfig` | `E0603`, `examples/ng_generic_loci_dump.rs:50` | keep `pub` |
| `PileupGeneratorCounts` | `E0603`, `locus_generation/mod.rs:442` | keep `pub` — `pub(super)` compiles but leaves a `pub` enum variant's payload unnameable |
| `WalkerError` | `E0603`, `locus_generation/mod.rs:555` | keep `pub`, same reason |
| `MAX_RECORD_SPAN_CEILING` | nothing — `unused import` | keep `pub`: no build consumer, but an external caller needs it to assert against `ceiling` |
| `PileupGeneratorConfigError` | nothing — `unused import` | keep `pub`: it is the `Err` of a `pub fn new`; `?` does not need it nameable, matching on it does |
| `PileupWalker` | nothing, non-test **and** `cfg(test)` | **delete the re-export** |
| `RunSummary` | nothing outside `pileup` | **demote to a private `use`** |
| `run` | nothing outside `pileup` | **demote to a private `use`** — this is the bypass above |
| `DEFAULT_MAX_ACTIVE_READS` | nothing outside `pileup` | **demote to a private `use`** |

So Checkpoint C's "nine of ten have no consumer" was wrong in both directions: four are load-bearing,
two have no build consumer but a real external reason, and **four should be demoted or deleted**.
Demote rather than delete anything intra-doc-linked (`mod.rs:74`, `generator.rs:95`) or
`cargo doc` gains a thirteenth error.

## C. The measurements, and what §13 asked for

### Major — §13.2 asks for the fabrication triple twice; the census ships it once

Spec §13.2 asks for *"how many loci, how many reads, and how many reference bases production
credits to reads that never sequenced them — and, **separately, the same three numbers** for the
reads `widen` extended after they had already left the active set"*. `DivergenceCensus` carries
`fabricating_loci` / `fabricated_reads` / `fabricated_ref_bases` for class 1, and for class 6 a
bare `stale_widen: usize` incremented once per locus. `parity.rs:1665-1672` states the
requirement — *"§13.2 asks for the two numbers **separately**"* — and then meets half of it.

So the census can say production mis-folds reads at 264 loci (3,074 in the soak) and cannot say
whether that is 3,074 reads or 300,000, or how many reference bases moved. Class 6 is the class
this milestone **discovered**, so it is the one a reader has least intuition for. `stale_widen_shape`
already walks exactly the rows the count needs.

### Major — the dump's counts header cannot balance, and its doc says it does

Three of the five read counters are once-per-read run-level values; `reads_without_observation`
and `reads_discarded_by_cap` are **per-locus** quantities the tool sums over loci. They print on
one `#` line under a module doc reading *"Every read the walk saw has to be somewhere… Printing
all of them is what makes 'every fetched read is accounted for' a thing a reader can check rather
than a claim (spec §13)."* Measured on the accounting fixture with both depth caps forced to 1:

```
# reads_admitted=4 reads_declined_by_preparer=0 reads_silent_over_footprint=1 \
  reads_without_observation=0 reads_discarded_by_cap=55
```

`reads_discarded_by_cap = 55` against `reads_admitted = 4`. The only identity actually asserted,
`reads_complete + reads_observed == Σ row.reads`, says nothing about admission — and the code
comment beside it knows so, while the module doc a user reads does not.

**Spec §13's own bullet has the same defect**, asking for the identity in a form that mixes the
two units, so this is a spec edit as much as a code one. The balancing quantity that *can* exist
is "distinct read ids appearing in at least one row", which needs a run-level flag beside
`ever_contributed`.

Both new counters are correctly **outside** the parity module's per-locus read balance, and
`reads_declined_by_preparer` has **no oracle against production at all** by construction — the
parity module drives `genome_walk::run` directly rather than `PileupGenerator`, so no differential
in it can see that counter. Worth recording; not a defect.

### Major — D3's 461 MB peak RSS is the dump tool's buffer, not the generator's

`DumpReport.rows` is a `Vec<ObservationRow>` that grows for the **whole run**, and `render()`
then materialises the entire TSV as one `String`, so both are live at peak. `ObservationRow` is
**152 bytes** with five heap allocations behind it — including a fresh `contig: String` and a copy
of `ref_bases` **per row**. chr1 emitted `generic_loci=1541788` and rows ≥ loci, so the vector's
spine alone is **≥ 234 MB of the 461 MB**, before the per-row allocations and before the 72 MB
string.

The D3 report reads that number as *"Peak RSS is **lower** than production's, which is the memory
property spec §7 asked for holding at chromosome scale"*, and
`impl_plan/locus_generation_pileup_generator.md:303` repeats it.

**The number does not support the conclusion, and the shape is the one §7 forbids.** §7's clause
is that everything the generator holds is bounded by depth, not by region length — and this
buffer is bounded by region length. So the measurement is not merely uninformative about §7; it
measures the forbidden shape and reports the total as a pass. The direction happened to be
favourable, which is what made it easy to leave.

Confidence: **high** on the arithmetic and the mechanism, both re-checked against the struct in
the author's tree; **medium** on the word "dominates", since chr1 was not re-run.

## D. The dump tool's input handling, four Majors

New code, 1,440 lines, and the first thing in this plan that reads files a user chose.

| site | behaviour |
|---|---|
| `:514` | `eprintln!("error: {error}")` drops the whole `#[source]` chain — a missing BAM reports `error: reading the input files' read groups failed` with neither the path nor `No such file or directory`, both of which are in the chain |
| `:477-479` | `PVC_GENERIC_REGION_CHUNK_BP` swallows every parse failure via `.ok().and_then(…ok())`, so `0`, `-1`, `abc` and out-of-range all **silently disable chunking** (`generic_regions=1`, exit 0) — turning off the boundary check the knob exists for |
| `:308` | plain `at + chunk - 1` on the env value: `18446744073709551615` panics with `attempt to add with overflow` in debug, wraps silently in release |
| `:485-486` | the `#[must_use] VerificationHandle` is bound to `_verify` with no comment, so the library's `Drop` warning prints on **every successful run**, and the FASTA is read and verified twice from two separate caches |

Minor, same area: an unmatched contig filter yields an all-zeros report and exit 0; the read
filter's own drops appear in no counter, which the file's comment admits;
`build_index_if_missing: true` writes a `.bai` beside the input, undocumented; `begin_segment`
absorbs inverted, zero-width and off-contig regions as `Ok(None)` while only an unknown contig
errors.

**Verified and not filed:** Milestone C's error invariant holds and is *structural* — a bare `?`
on all four source types is compiler-rejected, and adding `#[from]` to a region-carrying variant
fails with `deriving From requires no fields other than source and backtrace`. Neither new
counter overflows. `MAX_RECORD_SPAN_CEILING` is pinned by two tests. A too-low column depth cap
is fully counted, not silent. `PileupGeneratorConfig::default()` passes the defaults rules.

## E. Documentation drift, and the completeness that matters

**Spec §3's owed edits are wider than two sentences.** The "five classes" count survives in
**six** more places, and the false anchor predicate in **five** more. Two of those repeats are
worse than the originals:

- `parity.rs:2039` is the **panic message an unlisted seventh class will print** — it will say
  "none of spec §3's five classes is present" while `DivergenceClasses` has six fields. The
  census's whole value is that an unnamed difference panics rather than being absorbed; the
  message should not undercount the named set.
- `arch/locus_generation_pileup.md:455-456` states the anchor as *"every folded read **spans** the
  final footprint"* — the false predicate **and** the span-derived vocabulary that §6's boxed
  warning exists to forbid ("a coverage tag derived from the alignment span is therefore blind to
  all four"). This is the arch doc restating the anchor in the words the whole §4/§6 change
  removed.
- `spec/locus_generation_pileup.md:1111` is a **second spec site**, in §13, outside the owed §3
  table.
- `impl_plan/locus_generation_pileup_port.md:131-132` — **plan 2, nowhere corrected**: the
  sentence that told plan 3 what to build.

Also: `parity.rs:1164` cites `two_read_groups_split_rows_without_moving_evidence`, which **exists
nowhere in the repository** — the guarantee is genuinely discharged elsewhere (`evidence_by_bases`
at every census locus, plus the dump's own two-group fixture), so it is a stale pointer, but on
this branch a doc citing a test that does not exist is the exact shape a reviewer must not trust.
Four sites still say `q_sum` is "rounded to 1e-9" after D3 replaced rounding with a tolerance, and
`comparable`'s doc links `[Q_SUM_GRAIN]` and `[round_q_sum]`, both deleted — invisible to
`cargo doc` only because `parity.rs` is `#[cfg(test)]`.

**The arch doc's inventory:** eleven structure divergences, plus the missing §9 above. `mod.rs` is
still credited with `PileupGenerator` and with a shim deleted at A0; `generator.rs` and
`mock_reference.rs` are absent; "eight copies" is now three per `copy_fidelity.rs:102`;
`tests.rs` is described as "verbatim" and as "dies with plan 3", both false; the generic
parameters, the field list and the `PileupGeneratorCounts` listing are all stale.

## Carried, with a decision attached

**`refold_live_reads`' contributor-skip** (Part 1, Major 2) — still carried; no category found a
counterexample either.

**`apply_events_into`'s `run_end.max(event_end)`** — checked against the events construction, as
Part 1 asked. Every call passes `events_overlapping(rec_pos, rec_end, &active.read)`: one read's
cursor. For a single read's CIGAR the reach is non-decreasing — a deletion of *L* at anchor *a*
reaches *a+L+1* and the next match sits at *a+L+1*, reaching *a+L+2*; at a tied anchor the
`Match → Insertion → Deletion` rank gives reaches 1, 1, *L+1*. So the `.max()` is a no-op on every
real input, which is why no test fails.

**But Part 1's reason was not quite right.** It said the two forms agree on every input the
debug-asserted precondition admits. That precondition is *non-strict* — sorted by
`(anchor, kind rank)`, non-decreasing — and it admits two `Deletion`s at one anchor with
`deleted_len` 5 then 2, giving reaches *a+6* then *a+3*, where the forms differ and the `.max()`
is what stops the witnessed run shrinking below what the read witnessed. What rules that input
out is the **cursor's CIGAR walk**, not the assertion.

**So: neither a fix nor a note.** Add the reach-monotonicity to the existing `debug_assert!`,
keep the `.max()` as the release-build defence, and say in the doc that it defends the *stated*
precondition rather than the cursor's output. That converts an untestable `.max()` into an
asserted precondition with a test that **can** fail — a two-deletions-at-one-anchor case trips
the new assertion.

## Process notes worth keeping

- **`rust-code-review/SKILL.md` needs the step-0 re-point in its prompt template.**
  `isolation: "worktree"` will start every future fan-out on `main`, so any review of a feature
  branch hits this. The same file's Milestone-C-era `dev.sh` advice ("by absolute path") is what
  misdirected the last fan-out and should name the agent's own copy instead.
- **~1.3 GB of build artefacts per agent**, not the ~2.4 GB assumed; 12 GB total was reclaimed
  after the fan-out.
- **Four probes and one candidate fix were lifted out of the agent worktrees before deletion.**
  A dirty worktree is evidence, but only until it is cleaned — extract before pruning.
