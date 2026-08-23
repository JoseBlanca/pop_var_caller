# Reliability checklist

**Purpose.** Behavioural correctness verified by tests — coverage, robustness, regression protection, test naming.

**Triggers.** Look at every public item, every non-trivial private function, every previously-fixed bug, every `unsafe` safety condition, every doc-comment invariant. Read the existing test files for the scope and identify gaps.

**Skip when.** Never skipped — even snippets without tests get findings (the missing tests *are* the findings).

## Rules

- **Tests exist** for every public item and every non-trivial private function. Each missing test is a finding with a proposed test name in `function_returns_expected_on_condition` form.
- **All tests pass** under `cargo test --all-targets --all-features`. If the orchestrator already ran this and quoted output, do not re-run; use that output. A failing test is at least Major; a failure on a correctness-critical path is Blocker.
- **Coverage classes.** Every function under review is covered against: happy path, every error variant it can return, boundary values (`0`, `1`, `MAX`, empty, single-element, very large), Unicode and non-ASCII strings, malformed input, and concurrent access where shared state exists. Each missing class is a finding.
- **Property-based or fuzz tests** are required for: parsers, serializers/deserializers, any pure function over a structured input domain, and any function whose correctness depends on an algebraic law (associativity, idempotence, round-tripping). Use `proptest` / `quickcheck` or a `cargo-fuzz` target.
- **Concurrency tests** are required when the code uses `unsafe`, `Arc`, `Mutex`/`RwLock`, atomics, channels, or shared `async` state.
- **No flaky tests.** Time-dependent tests use an injected clock; network-dependent tests use a mocked transport or a feature flag excluded from default CI; randomness uses a seeded RNG with the seed printed on failure; ordering-dependent tests use explicit synchronization, not `sleep`. Each violation is at least Major.
- **Doc examples** compile and run as doc tests. Use of `no_run` or `ignore` requires a comment explaining why.
- **Regression tests** exist for every doc-comment invariant, every previously-fixed bug, and every `unsafe` safety condition. The test would fail if the invariant were violated.
- **Test names** describe the behavior under test and the expected outcome (`parse_returns_error_on_empty_input`, not `test_parse_2`). The name alone communicates the bug on a CI failure.

## Restoring a mutation is part of running it

**A mutation that is not provably reverted is a defect you introduced.** Verify the revert by
**content** — `git diff HEAD -- <paths>` and read the `+`/`-` lines — never by a summary count: a
mutation usually edits one line *inside* code the diff already counts as added, so `--stat` is
identical before and after. Scripts that mutate in a loop are where this bites, because the log
they print is written before the last restore has landed.

Measured in this repo: two mutations from a verification script survived into a commit whose
message quoted a full-suite pass. The suite had genuinely passed on the clean tree minutes
earlier; the restores reached disk afterwards, and `--stat` could not see them.

## Mutation testing: prove the mutant differs before recording a survivor

You are in your own worktree so that you can change the code and re-run, which finds what reading does not. **But a mutation that does not change behaviour is not evidence of coverage.** Before recording a survivor, show the mutated code takes a different path on at least one fixture — print from both arms, assert the two answers differ, or diff the outputs. Then say which.

A no-op mutant is not a finding, and reporting one as a survivor sends the author to write a test against a hazard that does not exist. Three shapes account for most of them:

- **A guard on a condition no fixture reaches.** `genotypes.max(2)` is `genotypes` for every fixture in the file; a `lender != group` filter on a branch that group cannot enter is dead code.
- **A widening that never narrows.** `Vec::resize` truncates as well as extends, so deleting a preceding `clear()` changes nothing when every slot is overwritten.
- **A reordering with the same fixed point.** Two blocks of an alternating fit can be swapped without moving where it converges; only the intermediate reports differ.

Say plainly in your report how many mutations you ran, how many survived, and how many changed no behaviour — the three numbers are different and only the first two are findings.

## Challenge the tests that exist, not only the ones that do not

The pass below asks what input class is *missing*. This one asks the inverse, and it is where the Blockers are: **on one plan in this repo, ten of sixteen review rounds had a Blocker that was a test unable to fail** — not wrong code.

**For every existing test in the changed code, name what in its fixture makes the asserted failure reachable.** If nothing does, that is the finding, and the fix is the fixture rather than the assertion. Two causes recur:

- **A fixture where several wrong implementations give the same number.** Four different ways of averaging a per-window posterior coincide exactly when every window carries the same weight; three separate defects then leave the suite green. Ask of each assertion: *which wrong implementations also satisfy it?* If the answer is "several", the fixture is not discriminating and the test is decoration.
- **A fixture in a regime where the interaction under test does not occur.** A test of a *coupling* between two quantities, run where they do not overlap, passes whatever the code does — measured, a coupled fit tested at a depth where a real variant and a sequencing error are never confusable returned the right answer from **every** starting point, including one a hundred times off, and four mutations survived it. **Choose the fixture against the regime where the thing being tested stops being distinguishable**, not against the regime that is convenient.

A test that genuinely cannot separate two rules is not always a defect — sometimes no fixture can. When that is the case, say so on the test itself rather than leaving the next reader to rediscover it, and file the limitation rather than the test.

## Challenge tests (additional pass)

For each non-trivial function in the changed code, identify at least one input class **not** covered by existing tests, name the specific bug that input class would expose, and provide the test as code. Add these under a `## Missing tests` heading at the bottom of your output file (in addition to per-rule findings). Each entry:

- Proposed test name in `function_returns_expected_on_condition` form
- Input class covered (e.g. "empty input", "header longer than 4 KiB")
- Specific bug it would catch
- Test body as code (or an unambiguous specification)
