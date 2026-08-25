# ng read likelihoods — handoff at Checkpoint H, the plan complete

*2026-08-25. Branch `ng-calling-likelihoods`, worktree `../pop_var_caller-calling-likelihoods`,
plan [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). Written for
whoever picks this up next. It follows
[the Checkpoint F handoff](ng_calling_likelihood_handoff_checkpoint_f_2026-08-25.md), whose nine
items still hold except where item 5 below supersedes them.*

## Where the plan stands

**All fifteen steps are done — A1 through H2, Checkpoint H.** Step 7 exists as a set of pure
functions and the calling loop can consume them.

| milestone | steps | state |
|---|---|---|
| A–D — the generic path | A1…D1 | ✅ |
| E — the stutter distribution | E1–E3 | ✅ |
| F — the STR emission seam | F1–F3 | ✅ |
| G — the censored term | G1 | ✅ |
| H — the STR row | H1, H2 | ✅ **Checkpoint H** |

Suite: **4,532 passed, 0 failed, 14 ignored**, in the container and natively on the host; **219 in
`ng::calling::likelihood`**. Build with
`/Users/jose/devel/pop_var_caller-calling-likelihoods/scripts/dev.sh cargo test`; validate with
`cargo clippy --lib --all-features --tests -- -D warnings` — **not** `--all-targets`, still red on
`main` for unrelated reasons.

**The branch is seven commits ahead of `main` and not merged.** `main` was merged *in* on
2026-08-25 (`677b8fa4`), which brought the compiler pin to **1.98**. A second session is building
the calling loop on `ng-calling-loop`; it consumes this code and does not edit it, and it is on a
`main` that does not yet contain G1, H1 or H2. **Merging this branch is the first thing to decide.**

## The first thing to do

**Nothing in the plan.** Decide whether to merge to `main`, and if so do it — the calling loop is
blocked on the row's real signature, and until the merge it is coding against a `main` where
`censored_emission` is `unimplemented!()` and `ssr.rs` does not exist.

## The nine things a successor must not rediscover

**1. Four claims in the design documents were measured false and corrected this session.** Do not
re-derive them, and do not trust an exactness claim in `read_likelihoods.md` without running it.

- **§5.2 and §12's thirteenth test** said a read that ran out is *always* less discriminating than
  a whole read of the same bases. It is not: where one candidate is shorter than the stretch the
  read got through and the other longer, the truncated read separates them by **5.661 nats against
  1.586**, because under the longer candidate it needs no slippage at all. Corrected in four
  places, each with a note saying what moved.
- **§12's sixth test** asked that a read nothing explains leave every pair of genotypes exactly as
  far apart, *bit for bit*. No implementation of the formula can: `(a + k) − (b + k)` is not
  `a − b`. Corrected to a relative bound.

**2. Measure floating-point agreement against the entries, not against their difference — this
cost two wrong corrections in one day.** Units in the last place *of a difference between two
genotypes* measure the true rounding error scaled by `|entry| / |separation|`, and that ratio is
set by something unrelated to the code: the same fixture reports **16 units at 3 junk reads and
3,072 at 300**. My first correction to §12.6 quoted 64 and was measuring exactly that artefact.
Relative to the entries' own magnitude the answer is one `f64::EPSILON` and stays there.

**3. The STR mixture's three terms must be probabilities of the same event, and two of them are
deliberately different expressions.** Spec §2.1 writes the junk term `λ · U(o)` — a function of the
observation. A whole read allows one length; a truncated one allows every length at or above what
it witnessed. Before H2 fixed it, the junk term was a point mass while the other two were tails,
and a truncated read was preferentially explained as somebody else's DNA by **two orders of
magnitude**, widening with the tract. **But the junk count is a normaliser, not a membership test**
— it is floored at one length, so a read at a length no candidate reaches still has somewhere to
go, which is §4.5's whole purpose. **The seed is a real distribution and does test membership.**
Making the two one expression breaks one or the other.

**4. `emission` scores 22 length changes the support calls unreachable, and that is left alone on
purpose.** All 22 are part-repeat contractions that would leave the tract below one repeat — a
two-base `CA` tract scores a one-base read at 5.390e-4 while `unreachable_mass` counts that mass as
unplaced. It is milestone F's code; repairing it moves scores at one- and two-repeat tracts, so it
waits on spec §4.2's open contraction question. Pinned by
`the_scoring_and_the_support_disagree_only_where_the_open_question_says_they_do`, which asserts the
count and the direction.

**5. The reachability boundary is in three places and item 2 of the Checkpoint F handoff is
superseded.** `unreachable_mass` and `reachable_length_changes` now share one `contractable_repeats`
so they cannot drift; `enumerate_placements` states the same rule at **run** grain and differs
deliberately on interrupted tracts. Tests tie all three —
`the_placements_and_the_support_agree_on_a_pure_tract` fails on the byte-level mutation and moving
`contractable_repeats` alone fails 11 tests across both files. Before this session, moving either
one alone left the other green.

**6. What the calling loop owes, and the half that is still open.** The loop must convert the
genotype prior's seed — one entry per **candidate**,
`fill_seed_share_per_candidate` — into a distribution over the locus's **reachable lengths**,
parallel to `SsrLocusParameters::reachable_lengths`. **The open half**: whether a reachable length
that no candidate spells gets zero or the prior's geometric evaluated there. Both satisfy
everything the row writes down, and they differ exactly at the contaminant §4.5.1 exists for. Most
entries are zero either way — about 39 reachable lengths against 5 candidates at a dinucleotide
locus.

**7. The review fan-out found eight blocking defects across three steps and seven were tests that
could not fail.** Not wrong code — code with no evidence. The recurring shape is **a fixture in
which everything is identical along the axis under test**. Three traps worth carrying forward:

- `Fixture::of_groups` gives each read group its own slippage row, so a test that varies something
  *else* across groups passes under its own mutant. Use `Fixture::sharing_parameters`.
- A "junk" read must be past the cutoff of the **longest** candidate, not the shortest — one past
  a short candidate's cutoff is still inside a longer one's, so the read is not junk.
- Calling a pure function twice with identical arguments proves nothing. The first version of the
  per-sample test did exactly that and killed **0 of 17** mutations; what has teeth is pinning the
  junk floor's own size.

**8. Two operational things.** The shell's working directory drifts back to the **main checkout**
on every background-agent event — drive everything by the worktree's absolute path and check
`pwd` if a file looks wrong. And an edit anchored on a `fn` signature can match the **trait
declaration** instead of the `impl`; one such match deleted two struct definitions before it was
caught. Anchor on a body line.

**9. `fill_reachable_lengths` lives in `ssr_emission.rs`, not `ssr.rs`** — beside `SsrCandidate`
and `period_of`, which it needs. The tell that it belonged there was having to widen `period_of`
to `pub(super)` to reach it from the row.

## Open questions, none blocking

1. **Spec §4.2's contraction boundary** — prose and figures disagree; the branch follows the
   figures. Item 4 above waits on it.
2. **Which value a reachable length no candidate spells should carry in the contamination seed** —
   item 6.
3. **How two candidates spelling one length share that length's seed mass** — spec's own open
   question 3.
4. **`unreachable_mass` understates the loss on interrupted tracts** (Checkpoint F handoff item 3),
   unchanged.
5. **`arch/read_likelihoods.md` §4.1 still sketches the row's signature with nine loose
   arguments** and a scalar contamination fraction. The built shape groups the locus's five
   parameters into `SsrLocusParameters` (clippy's argument limit at 1.98) and takes
   `&[ContaminationView]` (spec §4.5.1 asks per read group). Worth repointing.
6. **Spec §12's items 14–19** are measurement runs needing genotypes end to end —
   [`calling_bakeoffs.md`](../../ng/impl_plan/calling_bakeoffs.md).

## Loose in the tree

`proptest-regressions/ng/alignment/stutter.txt` is untracked. It records a seed from a deliberately
injected mutant rather than a real failure, so it was kept out of the commits; delete it or keep
it, but it is not evidence of anything.

## The standing habit, and this session's evidence for it

**Every number about your own work is measured before it is written.** Wrong-first-time this
session: a censored-tail loss quoted at 2 parts in a million to 2 in a thousand that measures 0 to
2.06 in a hundred; a sweep size given as 12,960 that is 8,364; a junk-cancellation bound of 8 units
that is 64, then 64 that is the wrong *unit* entirely; a naive emission cost of 18 under a rule
yielding 27; a cache saving called "a factor of ten" that is 3.5; and a claim that the length
support removes the cohort, where it removes the *reads* and leaves a 2.2–2.6 fold swing. **Every
one was caught by a reviewer running it or by asserting it in a test.** Quotations from the
specification, by contrast, were correct throughout — it is the claims about one's own fixture that
go wrong.
