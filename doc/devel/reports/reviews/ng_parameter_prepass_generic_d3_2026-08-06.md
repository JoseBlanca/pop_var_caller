# ng step 4, D3 — review of `fit_by_profile_scan`

**Date:** 2026-08-06. **Reviewed:** commit `753bcd1f`. **Fixes:** the commit after it.
**Agents:** three, each in its own worktree detached at `753bcd1f`, covering ten categories.

| agent | categories | outcome |
|---|---|---|
| reliability | `reliability`, `errors`, `refactor_safety` | 10 new mutations, **1 survived**, 2 unlocalisable panic paths; 4 Majors, 1 Minor |
| structure | `module_structure`, `naming`, `idiomatic`, `smells` | 8 Minors and 4 nits, every one applied and rebuilt before filing |
| numbers | contract against the design, and every quantitative claim | 2 Majors, 1 more Major on the seam's generality, 6 Minors; **10 of 10 mutations reproduced exactly** |

## Verdict

**The scan does what the design says**, checked clause by clause by the numbers agent against arch
§4.2, spec §3.1 and the plan: one `for` over the ladder with no `break` and no second pass, so
nothing is a refinement stage or an early exit; summing log-likelihoods across ploidy groups is the
right composition, since the groups partition the cells and share the noise parameter; the tie
rule's two halves are genuinely bound, with reversing either alone failing the test. What the review
found was in what the code *reports when it fails*, in what the tests could not see, and — again —
in the prose.

**⛦ The zero-sites guard is written per ploidy and only a single-ploidy scan tested it.** With one
ploidy in the scan, "this ploidy's cells" and "every cell in the slice" are the same set, so
replacing the per-group sum with a whole-slice sum left the suite green. A read group covering a
haploid contig that contributed no sites is the ordinary case this scan exists for — one error rate
across every ploidy the group covered — and the guard is the only thing that names *which* ploidy
was empty; two frames down the message is "every cell carries zero weight", which names neither.

**⛦ Two panic paths reached a message nobody could localise.** The scan checked the *width* the
model appended but not the *values*. A model that goes wrong at one rung of 161 surfaced as
"cell 3, genotype 0: NaN is not a log-likelihood" — no rung, no ploidy. **The rung is the part no
later frame can recover**, because by the time the climb refuses the table the noise parameters are
gone. Both are now checked where the rung is in scope, with the weight distinguishing a fault from a
shape: a cell carrying no sites may legally say no genotype produced it, and a test pins that so the
check cannot be tightened into one that refuses a legal table.

**⛦ Convergent finding, both agents independently: the scan drops the climb's `converged` flag.**
`mixture_weights.rs` keeps it `pub(crate)` for exactly two readers — "the profile scan … and the
tests" — and neither used it. A climb that exhausted its passes is scored **below** its own summit,
so a neighbouring rung can win on that alone and the argmax stops being the argmax; there is no
channel to report it, since `FitTermination` covers the outer alternation of Milestone E. The
measurement behind this is D1's: on a four-cell fixture the slowest truth took 1,234 passes, against
a real table of up to 583 cells climbed 161 times. Now a `debug_assert!` naming the rung and the
pass count — debug only, because a slow climb is a data condition rather than a bug
(`spec/parameter_prepass.md` §3.1) and the release path must not abort on one. It fires on no
existing fixture.

**⚠ A wrong number in the author's own prose — the sixth round in ten — and this time a wrong
*story*, which is worse.** D3's report explained that the rail fixture's unnormalised columns "handed
the win to the flattest rung, on total column mass rather than on fit". The numbers agent
re-introduced the fault and printed the five rung scores: the flattest rung scored **worst of the
five**, about 160,000 nats below the winner. What the extra mass actually did was move the argmax
**one rung inward**, from 4 to 3. The first telling would send the next reader looking for a symptom
that does not occur; the real one — an argmax one short of the end — is far easier to mistake for
rounding. Two smaller faults in the same fixture: `flat` was one over the *genotypes* where a column
is a distribution over the *cells*, and the comment's general claim that "blending a distribution
with a flat vector does not preserve the sum" is false, since a convex combination of two vectors
that each sum to one sums to one. With `flat` correct the renormalising loop is a no-op, and is kept
so that it is a checked property rather than an assumption.

**⚠ "Two tests hold it" was half true.** Of the two tests for "every rung is scored", only the
two-humped curve kills the natural early-exit mutation. The rung-recording fixture **survives** it,
because the model is asked before the score comparison that would break — it can only catch an exit
taken before the model call. Said so now.

**⚠ The tie rule and the rail flag both presuppose an ordered one-dimensional ladder, and nothing in
the signature says so.** The generic binding is sound and tested. But spec §3.1 says the STR path
replaces the scan with a search from several starting points over three parameters: a set of starts
has no direction, and a three-parameter grid flattened into a slice has no single edge, so a point
interior in one axis and on the boundary of another would read `argmax_at_ladder_end == false` on a
railed fit. The field is `pub` and its own doc calls it "the only thing standing between a railed fit
and a plausible-looking number" — on the path with thousands of unexamined fits. **Recorded for the
STR plan; nothing to change here.**

## Applied

From the structure agent, every item built and re-verified before filing:

- **`WeightedCell` moved from a bound on the function to a bound on `NoiseModel::Cell`**, so the
  seam is stated in one place. With it on the function, a model whose cell knew neither its ploidy
  nor its site count compiled and failed only at the one call site that scanned it.
- **`ScanResult` moved to `profile_scan.rs`**, beside the one function that builds it; its rail-flag
  doc had been repeating the paragraph on `fit_by_profile_scan` almost word for word, and a contract
  stated in two files drifts. `fitting/mod.rs` now holds the seam and `FitTermination`, and dropped
  three imports with it.
- **`ScanResult::frequencies` and `MixtureWeightsFit::genotype_weights` are both
  `genotype_frequencies`** — the same numbers under two names, in two types that meet in one
  assignment. The *procedure* keeps its statistical name (`mixture_weights`, `fit_mixture_weights`);
  the *values* are named for what they are.
- **`Winner` deleted.** It was a second, private copy of `ScanResult` — three of four fields under
  the same names, with `log_likelihood` a bare `f64` in one and a `LogProb` in the other, exactly
  the wrapper `LogProb` exists to stop going missing. The result is now built inside the loop where
  `rung` is in scope; both decisions the refactor touched were re-mutated and both still die.
- **`ln_likelihood` → `ln_likelihood_row_major`, `log_likelihood` → `rung_log_likelihood`**: one
  letter apart, in the same loop body, meaning different things.
- **`plan`/`group` → `plans`/`plan`**: in this module "group" already means *read group*, the grain
  E1 will wrap a loop around, so `group.ploidy` read as if a read group had a ploidy.
- `weights_under` generalised so the rail test stopped carrying a twelve-line copy of it;
  `# Panics` completed; `let better` → `scores_at_least_as_well`.

## Mutation record

**Twenty-eight in total, twenty-eight killed.** Ten before review; ten more from the reliability
agent, of which one survived and two more reached unlocalisable messages; eight run here after the
merge, of which two survived and are closed.

The two that survived the merge are worth naming, because both were tests that could not fail:

- **`HashMap` in place of the `BTreeMap`** the determinism claim rests on. The order-stability test
  as first written used two ploidies with equal group scores — and a sum of equal numbers is
  order-independent whatever the container does, while at two keys a `HashMap` iterates them in the
  same order often enough to pass anyway. It now uses **six** ploidies with site counts spread over
  four orders of magnitude, so the partial sums differ in the last bits and the chance of a matching
  order is one in 720. The mutation now dies.
- **The uniform start's width.** My first two attempts at mutating it were no-ops
  (`genotypes.max(2)` is `genotypes` for every fixture here), which is worth recording as a lesson
  about mutation testing rather than about the code: a mutation that does not change behaviour is
  not evidence of coverage. `genotypes + 1` and a non-distribution start are killed by 13 tests each.

## Verification

Container throughout. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → **3,108** in the library binary
(from 3,100); `cargo test --doc ng::parameter_estimation` → 1 passed; `cargo doc --no-deps --lib` at
12 unresolved links, none in `parameter_estimation`. `ng::parameter_estimation` 198 → **206** tests,
`profile_scan` 11 → **20**.

## Left for the owner

- **⛦ Arch §5.2's coupled loop cannot be implemented as written, and E2 is where it lands.** The
  scan's signature takes no incoming frequencies, so §5.2's step 1 — "each read group's error rate,
  from the read-group table, **at the genotype frequencies the previous iteration produced**" —
  has nothing to accept them. The scan is a pure function of its read-group table, so step 1 returns
  the same rung every iteration, the stopping condition "every read group's winning rung is
  unchanged" is met at iteration 2 for **every** sample, and `MAX_COUPLED_FIT_ITERATIONS = 20`, the
  oscillation argument and "only the 157-in-1,707 multi-library samples iterate" all describe a loop
  that cannot run more than twice. Arch §5.2 half-concedes this in its own correction paragraph
  ("an earlier version of this paragraph said the frequencies were held fixed…"), but the
  consequence is sharper than it records. **A design question for E2.**
- **The architecture's module table still puts the scan in `fitting/mod.rs`** and has no row for
  `profile_scan.rs`. Left for the owner, as A4's `depth_bins.rs` row was.
- **Spec §3.1 point 1 still says the error rate is scanned "coarsely at first"**, which its own
  point 3 and arch §4.2 contradict. The implementation follows point 3; point 1's wording is the
  stale half.
- **The prefactor cost D2's review measured is still unclaimed** — three `lgamma` per cell per rung,
  ~98% of the inner loop, identical at all 161 rungs. `profile_scan.rs` is where a per-cell cache
  would live and does not have one.
- **`ng::types::GenotypeFrequency` exists and the scan reports bare `f64`.** Bare is right here;
  E2 is the only place the newtype boundary can be crossed, and is worth a check in its review.
- **Plan D3's oracles landed as toy-model analogues** — a rung index rather than an error rate,
  which is all that is available while nothing has read a locus. The literal form ("a table
  generated outside Phred 10–50") arrives with E1, and neither the plan nor the report said so.
