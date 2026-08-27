# Can the new tests fail? — a mutation review of `ng-prepass-handover`

Worktree `/Users/jose/devel/pop_var_caller-review-tests` at `346abaaf`. Baseline: `cargo test --lib`
is **4,928 passing, 0 failing, 11 ignored**, which matches the report's own table.

24 single-edit defects were planted in the production code the new tests cover, one at a time, each
restored and byte-compared against a backup before the next. **Eight survived.** The tree was clean
(`git status --short` empty) at the end.

---

## The survivors

### Blocker — none.

### Major 1. Three of the four ownership checks are never exercised

`src/ng/calling/run_parameters.rs:229, 244, 253` — the `its_own_read_group` calls in the
minted-error loop, the contamination loop and the substitution-rate loop.

`from_prepass`'s doc names four quantities that are checked against the sample's own read groups,
and describes what missing one costs: *"the value would land under a real identifier and score
another library's reads under this sample's chemistry"*. Only the **error-rate** call has a test.
Deleting each of the other three leaves all 36 tests in the file green:

| deletion | fate |
| --- | --- |
| `its_own_read_group(..., "a minted-error total")` removed | **lived** (36/36 pass) |
| `its_own_read_group(..., "a contamination estimate")` removed | **lived** (36/36 pass) |
| `its_own_read_group(..., "a substitution rate")` removed | **lived** (36/36 pass) |
| `its_own.contains(&group)` replaced by `its_own.contains(&group) || true` | died — `a_rate_under_another_samples_library_is_refused` |

`a_rate_under_another_samples_library_is_refused` mis-files an **error rate** (`generic[2].error_rate`),
so it is the only one of the four routes it can see. The other three checks fire in the suite only
as collateral of two *other* mutations (reversed sample order, contamination read under the first
sample's name), which is why they look covered and are not.

Fix shape: the refusal test wants three siblings, or one test parameterised over the four
quantities. The `what` label is also unchecked — see Note 3.

### Major 2. The repeat-tract length assertion cannot fire in any test

`src/ng/calling/run_parameters.rs:203-211`.

Mutation: the second `assert_eq!` compares `repeat_tract_by_sample.len()` against
`generic_by_sample.len()` instead of `of_each_sample.len()` — the exact defect the brief names,
two per-sample lengths checked against each other rather than against the run's sample table.
**Lived**, 36/36.

The only wrong-length test, `per_sample_results_that_do_not_cover_the_run_are_refused`, shortens the
**generic** list, so the first assertion always fires first and the second is dead. Two consequences:

- a run whose repeat-tract list is short would index past the end of it and panic with
  `index out of bounds` instead of the written message;
- the test's `#[should_panic(expected = "cover 2 samples and the run's read-group table names 3")]`
  is a substring of **both** messages — `"the SNP/indel results cover …"` and
  `"the repeat-tract results cover …"` — so even a test that did shorten the repeat-tract list could
  not tell which assertion fired. Making the expectation `"the repeat-tract results cover"` and
  adding the mirror case fixes both halves.

### Major 3. The seed's diversity moment is algebraically invisible to the tracing test

`src/ng/calling/run_parameters.rs:262-265` — the `joint.fitted_diversity()` argument.

Mutation: `joint.fitted_diversity()` → `ExpectedHeterozygosity::try_new(joint.expected_heterozygosity / 2.0).ok()`
(the diversity halved, the frequency untouched). **Lived**, 36/36.

The reason is arithmetic, not fixture: `every_field_of_the_assembled_run_comes_from_its_own_input`
checks only

```
seed.alpha_alt_total() / (seed.alpha_ref() + seed.alpha_alt_total()) == expected_alternative_frequency
```

and the builder sets `alpha_ref = A(1−f)`, `alpha_alt = A f`, so that ratio is **exactly `f` for any
total `A`** — the diversity cancels. The seam's second moment is therefore checked only where the
pair is *inconsistent* enough to panic inside `total_for_diversity` (which is what kills the swap in
the report's row 5, below). Any diversity that stays consistent with the frequency travels
unchecked. Asserting `seed.alpha_ref() + seed.alpha_alt_total()` against the total the two fitted
moments imply would close it.

### Major 4. A whole ploidy can be dropped from the substitution projection with the library green

`src/ng/parameter_estimation/ssr/mod.rs:768-772` — `substitution_rate_by_stratum`.

Mutation: `.filter(|(key, _)| key.ploidy == Ploidy::try_new(2).expect("two genome copies"))` inserted
before the `map`. **Lived — all 4,928 library tests pass.**

Both new fixtures are diploid-only: `three_strata_at_unlike_substitution_rates` builds
`SsrAccumulators::new(diploid())` and `key_at` hard-codes `Ploidy(2)`; the seam fixture's
`repeat_tract_parameters_of` sets `ploidy: diploid()` on all four keys. So the `StratumKey`'s third
axis never varies anywhere in the new work, and a haploid stratum's substitution rate silently
vanishing is invisible. This matters at the project's own stated range — a genome with a haploid sex
chromosome is exactly the case `GenericSampleParameters::rates` documents as "two entries".

The pre-existing `the_substitution_rate_is_keyed_by_the_runs_ploidy` (line 1406) does vary ploidy,
but it calls `assemble` directly and never touches the projection.

### Minor 5. `insert` → `entry().or_insert()` is invisible (and benign here)

`run_parameters.rs:232`, `error_rate_by_read_group.insert(...)` →
`.entry(group).or_insert(rate.clone())`. **Lived**, 36/36 — and it is *not* a defect: the two differ
only when one key arrives twice, and `its_own_read_group` refuses that first. Recorded so the
coverage claim is not overstated: the "first writer wins vs last writer wins" choice is untested
because it is unreachable, which is the right reason.

### Minor 6. The one-sample test cannot tell "not identified" from "never joined"

`a_run_of_one_sample_assembles_and_is_uncontaminated` asserts `view.contamination_is_absent()`, which
is also what dropping the contamination join entirely produces. Demonstrated in passing: under the
`take(len − 1)` mutation (row 6 of the table below) that sample's single estimate is dropped
outright, and the test stayed green — only the three-sample tracing test failed. The test's other
three assertions (sample count, read-group count, scale) are real; this one is satisfied by absence
of work as well as by absence of contamination. Asserting `read_group_count() == 1` alongside a
`was_measured`-style check, or asserting the fit's estimate reached the run before the flattening,
would separate them.

### Note 7. The `what` label on the ownership check is unchecked

`run_parameters.rs:229` — passing `"a minted-error total"` where the error-rate loop should pass
`"a fitted error rate"`. **Lived**, 36/36. The `should_panic` expectation is
`"which is not one of its own"`, which every one of the four shares. Since that string is the whole
point of `what` (the doc says *"`what` names the quantity so the message says which of the four went
astray"*), the expectation could name the quantity instead.

### Note 8. `its_own.is_empty()` as an escape hatch survives

`run_parameters.rs:597` — `assert!(its_own.is_empty() || its_own.contains(&group), …)`. **Lived**,
36/36. Low impact: a sample in `ReadGroups` always has at least one library, so the arm is not
reachable from real input. Listed for completeness of the battery.

### The five together

`R1 + R2 + R3 + R4 + R11` applied **simultaneously** — the repeat-tract length check pointed at the
wrong list, three of four ownership checks gone, and the seed's diversity halved — and
`cargo test --lib` reports **4,928 passed, 0 failed**.

---

## Fixture permutation check (brief item 2)

| fixture | two things alike? |
| --- | --- |
| `A_RUNS_LIBRARIES` (4 libraries × 6 columns) | **No.** Rates, minted means, fractions, tract rates all distinct, and the four calibration scales come out 0.5 / 0.25 / 0.75 / 0.125 — also distinct. Sample column `[0,0,1,2]` keeps library index ≠ sample index. This one is built the way the report claims. |
| `A_RUNS_INBREEDING` `[0.10, 0.20, 0.30]` | No. |
| `the_shared_stratum` (one stratum for all four libraries) | Deliberately shared, and the doc comment gives the right reason — a swapped rate answers with a neighbour's number rather than *absent*. Correct as written. |
| **`repeat_tract_parameters_of` / `key_at` — ploidy** | **Yes — every key is diploid.** Survivor 4. |
| **`three_strata_at_unlike_substitution_rates` — read group** | **Yes — all three strata are `ReadGroupId(0)`** (`tract_at_a_known_rate(…, 0, …)` and `key_at`'s hard-coded id). Forcing `read_group: ReadGroupId(0)` in the projection leaves the ssr test green; only the seam's tracing test, in another file, catches it. The test's doc claims "a key dropped, added or paired with a neighbour's rate all show" — true on the stratum axis, false on the read-group axis. |
| **`two_declared_batches` — samples 1 and 2 share a batch** | **Yes.** The test asserts `batch_of_sample(1) == BatchId(1)` and `batch_of_sample(2) == BatchId(1)`: two inputs with the same value, so the batching axis cannot see samples 1 and 2 exchanged. Batches `{rg0,rg1}` / `{rg2,rg3}` are forced by the library-to-sample map; three batches, or asserting the read-group→batch map instead, would restore the permutation. |
| `UNLIKE_LIBRARIES` (estimate.rs) | No — 8 reads at −7 nats against 12 at −9. Both halves of each library's total differ from the other's, and the arithmetic half of the assertion is genuinely independent of the tally. |

---

## The report's mutation claims (brief item 3)

All five rows of the "Five deliberate defects, all killed" table reproduce, with the attributions as
written:

| the report's row | reproduced? |
| --- | --- |
| samples read in reverse order → tracing test + refusal test | **Yes.** `generic_by_sample[len − 1 − index]`: 2 failed, exactly those two. |
| contamination looked up under the first sample's name → both | **Yes.** `joint.contamination.values().next()`: 2 failed, exactly those two. |
| tract rates taken from the first sample's → tracing test | **Yes.** `repeat_tract_by_sample[0]`: 1 failed, the tracing test. |
| ownership check disabled → the refusal test | **Yes.** `contains(&group) \|\| true`: 1 failed, `a_rate_under_another_samples_library_is_refused`. |
| the seed's two moments swapped → tracing test + one-sample test | **Yes, but by a mechanism worth stating.** A literal swap does not compile (the two moments are different newtypes); wrapping each raw number into the other's type gives 2 failures, those two. The one-sample test fails because an inconsistent (f, d) pair **panics inside `total_for_diversity`** — that test asserts nothing about the seed. So the pair is guarded only against inconsistency, not against a wrong-but-consistent diversity (survivor 3). |

Counts:

- **"the permutation left all 201 other tests in that file green"** — the number holds, the phrase
  "that file" does not. Pairing each key with its neighbour's rate leaves
  `cargo test --lib ng::parameter_estimation::ssr` at 201 passed / 1 failed / 1 ignored — but that
  filter spans the whole `ssr::` tree (203 tests over five files, one of them ignored:
  `ssr::slippage::tests::the_search_recovers_a_known_truth_and_no_start_beats_it`). `ssr/mod.rs`
  itself holds 91. The mutation also killed `every_field_of_the_assembled_run_comes_from_its_own_input`
  in `run_parameters.rs`, which the report does not mention (it was written before the seam existed).
- **"Both left the other ten tests in that file green"** (estimate.rs) — **holds, twice.** Emptying
  the map: 10 passed, 1 failed. Reversing the values against the keys: 10 passed, 1 failed. The file
  has exactly 11 tests.
- **"now: 4,928 passing, 0 failing, 11 ignored"** — **holds** on the current tree.

---

## Every mutation tried

| # | file:line | mutation | fate |
| --- | --- | --- | --- |
| 1 | run_parameters.rs:203 | 2nd length assert vs `generic_by_sample.len()` | **lived** |
| 2 | run_parameters.rs:229 | minted-error ownership check deleted | **lived** |
| 3 | run_parameters.rs:244 | contamination ownership check deleted | **lived** |
| 4 | run_parameters.rs:253 | substitution-rate ownership check deleted | **lived** |
| 5 | run_parameters.rs:232 | `insert` → `entry().or_insert()` | **lived** |
| 6 | run_parameters.rs:243 | contamination loop `.take(len − 1)` | died — tracing |
| 7 | run_parameters.rs:216 | sample order reversed on both lists | died — tracing + refusal |
| 8 | run_parameters.rs:236 | contamination from `values().next()` | died — tracing + refusal |
| 9 | run_parameters.rs:217 | `repeat_tract_by_sample[0]` | died — tracing |
| 10 | run_parameters.rs:597 | ownership assert `\|\| true` | died — `a_rate_under_another_samples_library_is_refused` |
| 11 | run_parameters.rs:264 | fitted diversity halved | **lived** |
| 12 | run_parameters.rs:227 | coefficient taken from sample 0 | died — tracing |
| 13 | run_parameters.rs:229 | `what` label crossed | **lived** |
| 14 | run_parameters.rs:616 | `.chain(minted.keys())` dropped from the id union | died — `a_minted_total_without_its_fitted_rate_is_refused` |
| 15 | run_parameters.rs:597 | ownership assert `its_own.is_empty() \|\|` escape | **lived** |
| 16 | run_parameters.rs:262 | seed's two moments swapped (re-wrapped) | died — tracing + one-sample |
| 17 | estimate.rs:201 | `minted_errors: BTreeMap::new()` | died — new estimate test (10 others green) |
| 18 | estimate.rs:201 | keys zipped against `values().rev()` (the two totals swapped) | died — new estimate test (10 others green) |
| 19 | ssr/mod.rs:768 | keys zipped against `values().cycle().skip(1)` | died — new ssr test + seam tracing test |
| 20 | ssr/mod.rs:768 | `.filter(ploidy == diploid)` | **lived — whole library** |
| 21 | ssr/mod.rs:768 | key rebuilt with `read_group: ReadGroupId(0)` | died — **seam test only**; ssr test green |
| 22 | joint/fit.rs:405 | `fitted_alternative_frequency` returns `expected_heterozygosity()` | died — **seam test only** |
| 23 | parameter_estimation/mod.rs:281 | "no default" dropped from the message | died — new message test |
| 24 | — | mutations 1+2+3+4+11 applied together | **lived — 4,928 passed, 0 failed** |

Two intended mutations turned out not to be expressible, which is a point in the design's favour:

- **crossing the error-rate union with the minted union** — the two maps hold
  `Estimate<ErrorRate>` and `MintedReadErrors`, so every crossing of them, at the loops or at the
  `assemble` call, is a compile error;
- **swapping the seed's two moments literally** — `ExpectedAlternativeFrequency` and
  `ExpectedHeterozygosity` are distinct newtypes, so the swap only exists once each raw `f64` is
  re-wrapped by hand (mutation 16).

Also noted: `JointFit::fitted_alternative_frequency`, added on this branch, has no test of its own —
mutation 22 is caught only by the seam's tracing test one module away.

## Safety

Every file was backed up under `tmp/*.bak` before its first edit and restored with `cp`, then
verified with `diff <backup> <file>` after **every** mutation (all reported `RESTORED-CLEAN`). Final
`git status --short` is empty; `tmp/` is gitignored. Nothing was fixed — report only.
