# Code review — ng calling foundations, step B3 (`LocusInference`, `SampleGenotypeCall`)

*2026-08-21. Branch `ng-calling-foundations`, reviewed at `56b65ae0` (steps B1+B2) with the step's
uncommitted working-tree diff applied. Two category sub-agents, each in its own git worktree.
Per-category audit trail in the gitignored `tmp/review_2026-08-21_ng-calling-b3/`.*

## 1. Scope

**What was reviewed:** the working-tree diff of step B3 of
[`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md) — 287 insertions and 2
deletions, all in `src/ng/calling/mod.rs`: `SampleGenotypeCall`, `LocusInference`, its constructor,
five tests and two import lines.

**Deliberately out of scope:** `CandidateAlleles` and `ExpectedAlleleCopies` in the same file,
reviewed and committed as steps B1+B2.

**Categories dispatched.** Two agents for 287 lines adding two record types:

| category | reason |
|---|---|
| `reliability` | always; the only category that mutation-tests. Asked directly about the public fields versus the constructor, the `passes` bound, whether sample order is checkable, and whether the STR seed marker could be enforced off the SNP/indel path |
| `naming` + `idiomatic` + `smells` + `defaults` + `errors` | **one agent for five**: two plain record types, one constructor, no new error type |

Not dispatched: `module_structure` (no module moved; the folder was reviewed at B1+B2),
`unsafe_concurrency`, `tooling`, `extras`, `refactor_safety`.

## 2. Verdict

**Approve-with-changes.** No Blockers. Five Majors and six Minors, all applied except one
recommendation this report disputes on the evidence (§6, "The one disputed recommendation"). The
mutation pass ran **15 mutations, 2 survived, 0 changed no behaviour**, and both survivors were
weakenings of the constructor's own checks.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt` | 0 | clean |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 3.83s` |
| `cargo test --lib ng::calling` | 0 | `17 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |

Both agents built in their own worktrees. Findings labelled "Needs verification": **zero**.

## 4. Open questions and assumptions

1. **Is a locus that converged on its first pass representable?** Raised by `reliability` as an
   author question rather than a finding, and answered in §6.
2. **Where should `Provenance`'s ordering live?** The field name `weakest_provenance` presumes an
   ordering that `Provenance` does not define. Out of this step's files; carried as a follow-up.

## 5. Top 3 priorities

1. **M4 — the constructor's invariant held only at the instant it returned.** `alleles` was a
   public field and `CandidateAlleles::admit` takes `&mut self`, so a value that passed through
   `new` could be widened against unchanged copies afterwards, from a call site that looks like
   ordinary use of the discovery round.
2. **M3 — the suite's only coverage of the STR seed marker sat on a SNP/indel fixture**, the one
   kind the field's own documentation says can never raise it.
3. **M1 — the copies/alleles check was tested in one direction only**, and the untested direction
   is the one the pipeline produces.

## 6. Findings

### Major

**M1: the copies/alleles pairing was tested only where copies are *narrower*.** **Category:**
reliability. **Confidence:** High. Weakening `assert_eq!` to `assert!(copies.len() >=
alleles.len())` survived all 17 tests; a probe showed pristine panicking on 3 copies over a
2-allele table where the mutant returned a built `LocusInference`. **Wider is the direction the
pipeline produces**: the final prune *shrinks* the allele table, so a copies vector not re-cut
alongside it is left long, and since every consumer indexes by `AlleleId` the trailing entries
would ride along unread — a pruned-away allele's expected copies travelling in the output with
nothing to notice.

**M2: `passes == 1`, the documented lower bound, was never built successfully.** **Category:**
reliability. **Confidence:** High. Tightening the bound to `passes > 1` survived all 17 tests: the
one fixture passing `passes: 1` panicked on the *earlier* copies assert before the `passes` check
was reached, so no test ever observed a locus accepted at the boundary.

**M3: the only test setting `seed_diversity_unreachable` used a `LocusKind::Generic` fixture.**
**Category:** reliability. **Confidence:** High. The field's own doc says it is never set on the
SNP/indel path. The agent proved this was not a misreading: adding a guard rejecting that
combination compiles and fails exactly that one test. So the suite's only coverage of the flag
being `true` sat in the regime where the quantity cannot arise, and would have stayed green under
an implementation that wired the STR seed marker onto the generic path.

**M4: `new`'s invariants were `&mut`-reachable after construction.** **Category:** reliability.
**Confidence:** High. Two probes:

```
PROBE-LITERAL: built with passes=0 and 3 copies over 2 alleles, no panic
PROBE-MUTATE: after new(), admit() widened alleles to 3 against 2 copies; passes now 0
```

The second is the one that matters — a struct literal is visible in a diff, but
`inference.alleles.admit(…)` looks like ordinary use of the discovery round. The agent's remedy is
the reason this was cheap: **only two of the eight fields carry a cross-field invariant**, so
closing it costs two one-line accessors, not eight, and the other six stay public with the house
style intact. `#[non_exhaustive]` would not help, since the loop that builds these is inside the
crate.

**M5: the constructor pays a positional constructor's full cost for a check a struct literal
skips.** **Category:** idiomatic. **Confidence:** High. Eight positional arguments with `bool` at
positions 5 and 8. The agent answered the lint question empirically rather than from memory: it
added 2-, 3- and 4-bool probe functions and ran clippy, and only the 4-bool probe warned — so
`fn_params_excessive_bools` staying silent here is a threshold artifact (the default is 3), not
evidence the shape is safe. Its recommended remedy was a struct literal plus a `checked(self)`
method; M4's private-field remedy was taken instead, which subsumes the bypass half.

### Minor

- **Mi1 (errors, reliability — convergent): an empty `per_sample` was neither rejected nor pinned**,
  and two of the five tests built a locus naming no sample's call at all. Adding a rejection left
  all 17 tests green, which is what makes it a gap rather than a decision.
- **Mi2 (naming): the module header's inventory was not updated** for the two new types, though they
  meet its stated criterion.
- **Mi3 (naming, reliability — convergent): `weakest_provenance` cited a section that does not
  exist.** `arch/read_likelihoods.md` runs §0 to §7; provenance propagation is its §1.4, and the §8
  is in the *spec*. Inherited from a stale reference in `arch/calling_priors.md`.
- **Mi4 (idiomatic): `weakest_provenance` names an ordering `Provenance` does not define.**
  `FittedHere`, `Borrowed`, `Defaulted` and `Supplied` are four names, not a scale, and where a
  *supplied* value sits against a *fitted* one is genuinely open.
- **Mi5 (reliability): a test doc claimed a case a length check cannot reach.**
- **Mi6 (naming, on the diff's factual claims): the seed-marker doc stated as settled what the
  architecture marks provisional.** The arch says the ceiling behaviour stands "until Q2 settles the
  policy"; the field's doc stated it flatly.

### Nits

`new`'s doc opened directly with `# Panics`, leaving the item without a summary line; nothing
checked `region.start <= region.end`, though `GenomeRegion`'s own doc says a caller requiring
well-formedness must say so and a called locus is such a caller; the pass count `50` is a literal
that will need to become `DEFAULT_MAX_PASSES` when that constant lands; `passes` reads as a
collection at a use site where `pass_count` would not (the arch fixes the name); and
`seed_diversity_unreachable` parses most naturally as *the seed's diversity is unreachable* when the
quantity that could not be reached is the measured gene diversity (the arch fixes this name too).

Also recorded as accepted with reasons rather than filed: the two-`bool`-fields smell (all four
truth-table rows are meaningful), the eight-field threshold, `assert!` over `Result` (both
quantities internally computed, matching the crate's stated line), and the absence of any `Default`
impl — a derived one would produce `passes: 0`, the value `new` rejects.

### The one disputed recommendation

`reliability` recommended asserting that `converged == true` requires `passes >= 2`, reasoning that
convergence is a comparison between two passes so there is nothing to compare against after one,
and that a locus reporting `converged: true, passes: 1` would mean the loop compared against a
reused scratch buffer's stale contents. It filed this as a recommendation plus an author question
rather than a finding, which was the right call — **and the answer is that a first-pass convergence
is legitimate.**

`spec/calling_em_loop.md` §2's loop begins with an *initialisation* E-step outside the frequency
loop — "one E-step with NO prior — reads only — then sum: the cohort's expected copies" — and §3
explains why: the prior needs the cohort's expected copies, which do not exist until a pass has
produced them, so the first pass runs on the reads alone. The frequency loop's first pass therefore
compares its own M-step output against that reads-only estimate, which is freshly computed and not
stale scratch. A locus with `converged: true, passes: 1` means the prior barely moved the copies
away from what the reads alone said — a real outcome, and one this caller's hardest committed case
(a cohort of one at three reads a position) is where to expect it.

**Not applied**, and the field's doc now says so, so the next reader does not re-derive it. A test
pins the accepting side of the bound instead (M2).

## 7. Out of scope observations

- **`Provenance` has no ordering** (`src/ng/parameter_estimation/mod.rs`). The comparison this
  step's field name depends on has no home yet, and it blocks the loop step rather than this one.
- **`arch/read_likelihoods.md` §8 is cited and does not exist**, both here and in
  `arch/calling_priors.md` §5. A design-doc fix, which this loop does not make.
- **The three aggregate gates red on `main`** are unchanged; recorded in
  [A1's review](ng_calling_a1_2026-08-21.md) §7.

## 8. Missing tests to add now

All applied: the wider direction of the copies/alleles check; a locus accepted at `passes: 1`; the
STR fixture for the seed marker plus a SNP/indel locus rejected for carrying it; an empty
`per_sample` rejected; and a backwards region rejected.

## 9. What's good

- **Two agents, both measuring rather than asserting.** One proved the lint's silence was a
  threshold artifact by adding 2-, 3- and 4-bool probes and running clippy; the other proved every
  survivor changed behaviour before writing it up, and proved the seed-marker guard lands clean by
  adding it and counting which tests fail.
- **The remedy was scoped to what actually carries an invariant.** "Two accessors, not eight" is
  what made M4 cheap enough to take, and it came from counting the fields that participate.
- **The transposition worry I raised was answered with evidence and partly dismissed.** The fixtures
  do set the two flags to opposite values, so a swap in the constructor's body fails two tests — the
  agent said so rather than agreeing with the premise.
- **The disputed recommendation was filed as a question, not a finding**, with the spec sentence it
  rested on quoted — which is what made it checkable against the sentence that overrules it.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::calling
./scripts/dev.sh cargo test --all-targets --all-features   # expect the pre-existing bench panic
```
