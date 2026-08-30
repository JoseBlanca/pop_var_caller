# Code Review: ng calling loop — A1, the three shared types

**Date:** 2026-08-25
**Reviewer:** rust-code-review skill (orchestrator), four category sub-agents in isolated worktrees
**Scope:** the working-tree diff of step A1 of [calling_loop.md](../../ng/impl_plan/calling_loop.md), captured as commit `5843f60a` over branch point `bbcf2165`
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** a diff — `git diff bbcf2165 5843f60a`, 1,273 insertions across two files.
- **Reviewed against:** commit `5843f60a45e4e1de5f9565475059eb74821debf8` on branch `ng-calling-loop`, kept alive as `refs/review/ng-calling-loop-a1`.
- **In-scope files:**
  - [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs)
  - [doc/devel/reports/implementations/ng_calling_loop_a1_2026-08-25.md](../implementations/ng_calling_loop_a1_2026-08-25.md)
- **Deliberately out of scope:** `src/ng/calling/likelihood/` and `src/ng/calling/allele_candidates/`, which two other sessions own on their own branches — consumed here, never edited. Everything outside `src/ng/calling/`.
- **Categories dispatched, and why:**
  - **reliability** — always; and this diff is almost entirely invariants and their tests.
  - **errors** — always; the module deliberately has no `Result`, so whether each check is held at the right level is the whole question.
  - **naming** — always; the diff introduces about thirty new public names.
  - **defaults** — the diff adds a default type parameter, a derived `Default`, a poison value and an in-band "absent" signal.
  - **Not dispatched:** `module_structure` (one source file, and the type placement is settled by the architecture), `unsafe_concurrency` (no `unsafe`, no shared mutable state in the diff), `tooling` (a diff, not a crate), `smells` and `idiomatic` and `refactor_safety` and `extras` — a **deliberate trim**, recorded here rather than left to be inferred: four agents each hold a 3–5 GB worktree and two other sessions were building on the same machine's container runtime at the time. The three trimmed always-on categories are the gap in this review's coverage.

## 2. Verdict

**Request-changes.** One Blocker and thirteen Majors, of which twelve are fixable inside this step. The code as written is correct on every path the tests reach; what the review found is a set of paths the tests do not reach, three checks the spec asks for and the diff does not make, and one API shape that reproduces a rule the sibling module already ships a named spelling for.

## 3. Execution status

Run by the orchestrator in the container from the reviewed tree:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features --tests -- -D warnings` | 0 | no warnings |
| `cargo test --lib` | 0 | `4502 passed; 0 failed; 14 ignored` (branch point: `4488 passed; 0 failed; 14 ignored`) |

**`cargo clippy --all-targets` was not run, and the reason expired while this review ran.** It was excluded because it was red at this branch's compiler pin, 1.97.1, from pre-existing lints in examples and benches with none under `src/`. On 2026-08-25 `main` moved the pin to 1.98 (`54a0fd96`) and made that command exit 0 (`f3c8c797`), which this branch does not yet carry. **So the wider scope is available once this branch takes main's compiler pin, and running it is on the fix-application checklist rather than settled here.**

Run by sub-agents in their own worktrees, quoted from their files:

- `cargo test --lib ng::calling` on the pristine tree → `474 passed; 0 failed; 3 ignored`.
- **All 16 release-held assertions in the implementation section downgraded to `debug_assert` and the module's tests run under `--release`** → `19 passed; 16 failed`. The failures map onto 15 of the 16 assertions; exactly one survives untested (finding **M9**).
- `cargo test --release --lib ng::calling` on the pristine tree → `461 passed; 4 failed`, the four failures all pre-existing in `likelihood/` (see §7).

Findings labelled "Needs verification": **1** (Mi10, whose fix is proposed but was not compile-verified — the sub-agent's verifying run was killed by a session interruption and its file says so).

**Mutation testing, reliability category: 6 mutations run, 2 survived, 0 changed no behaviour.** Both survivors are proved rather than assumed: a tree carrying the shipped code plus the six proposed tests runs `41 passed; 0 failed`, and re-running the first survivor against that tree fails on `ssr_evidence_against_a_bundle_allele_table_is_accepted`. The written file still carries an unfilled `RESULTS_SURVIVED` placeholder on its summary line — the numbers here come from the sub-agent's closing report, and its per-finding mutations each name their own outcome in place.

**Five queued mutations never ran, and the cause was the orchestrator's.** Part-way through that batch the reliability worktree's `Cargo.toml` and git registration were removed by this session's worktree cleanup, ending all builds in it. The findings that would have rested on those five are stated from grep-level facts and the passing baseline instead, and each says so where it stands. `src/ng/calling/mod.rs` in that worktree was restored and verified by content against a pristine copy of `5843f60a` before it was lost. **The gap is real: five mutation results this review would otherwise have are missing**, on top of the three always-on categories §1 records as untrimmed.

## 4. Open questions and assumptions

1. **Is `SampleGenotypeCall::Missing` the right name, given that a zero-coverage sample is `Called`?** Affects **Mi6**. The spec's own word for the output is *missing* (§9), and the variant qualifies the *call*, not the data — but a geneticist may read "missing data". Resolved below by keeping the name and making the zero-coverage case explicit in the doc, rather than by a rename that would diverge from the spec's vocabulary.
2. **Should the module's release-held assertions be pinned by CI?** Affects **M14**. They are release-held by the source text and by nothing else. Adding the gate is blocked behind four pre-existing `--release` failures in another branch's module, so it cannot be done here.
3. **Where does the genotype-table-versus-allele-count check belong — A1's `prepare_for`, or A2's `call_locus`?** Affects **M2**. Both sub-agents that raised it note the seam is A2's. Resolved below in favour of A1, because `prepare_for` is the one point where the shape is fixed.

## 5. Top 3 priorities

1. **B1** — two cells of the path-agreement matrix are untested, and one of the surviving mutations routes a repeat bundle to the SNP/indel row: a plausible genotype at every sample, nothing failing.
2. **M1** — a `CallingScratch` that was never sized answers six of its eight accessors with an empty slice, so the cohort's expected copies come back as a plausible `0.0` — the exact outcome the `NaN` poisoning exists to prevent, through the door the poison cannot cover.
3. **M2** — the spec's *first* named caller bug, a genotype table that disagrees with the allele count, is asserted nowhere, because no function in the diff holds both objects.

## 6. Findings

### Blocker

**B1: src/ng/calling/mod.rs:466 — the repeat-bundle cells of the path-agreement matrix are untested, and two mutations survive**
**Categories:** reliability. **Confidence:** High.
`LocusKind` has three variants, so the (evidence variant × kind) matrix has six cells; three are exercised and none of them involves `LocusKind::SsrBundle`. Measured: deleting `| LocusKind::SsrBundle` from the `Ssr` arm, and separately widening the `Generic` arm to `LocusKind::Generic | LocusKind::SsrBundle`, each leave `474 passed; 0 failed`. `SsrBundle` is a live kind — `run/cohort_merge/close.rs:120`, `parameter_estimation/joint/census.rs:2091` and `parameter_estimation/ssr/locus_offsets.rs:92` all branch on it, and `close.rs:120` groups it with `Ssr(_)` exactly as this method does. The second mutation is the one that bites: a bundle whose alleles say `SsrBundle`, handed SNP/indel evidence, would pass the check and be scored by the generic row.
**Fix:** add the accept cell and the refuse cell as tests.

### Major

**M1: src/ng/calling/mod.rs:659 — six of eight scratch accessors answer an unprepared instance with an empty slice**
**Categories:** defaults. **Confidence:** High.
`#[derive(Default)]` is the only constructor and produces a value unusable until `prepare_for`. `lg_row` and `sample_copies` catch the omission through `row_range`'s assertion; `expected_copies`, `expected_copies_mut`, `expected_copies_prev`, `advance_expected_copies`, `seed_concentration` and `seed_concentration_mut` return the field directly, which on an unprepared scratch is an empty `Vec`. Measured by the sub-agent: the sum over `expected_copies()` is `0.0` and a fold over `expected_copies_mut()` runs zero iterations and writes nothing. Step B2's oracle is a **bitwise** comparison of summed expected copies, which an all-zero sum passes against another all-zero sum.
**Fix:** one `assert_prepared()` helper, called first in each of the six.

**M2: src/ng/calling/mod.rs:711 — the spec's first named caller bug is asserted nowhere**
**Categories:** errors; reliability (cross-category). **Confidence:** High.
Spec §8 names three caller bugs to assert: a genotype table that disagrees with the allele count, a non-finite likelihood reaching the loop, a sample count that disagrees between evidence and parameters. The third is checked; the second belongs to the loop A1 does not build; the **first is checkable here and is not**. `prepare_for` sees the `GenotypeTableView` and never the `CandidateAlleles`; `assert_agrees_with` sees the alleles and never the table. `GenotypeTable::build` takes the allele count as a bare `usize`, and a discovery round admitting an allele is exactly what makes a correct table stale. With a table one allele narrow, every per-allele buffer is sized for the old count and the first thing to notice is `ExpectedAlleleCopies::new` at the far end of the locus.
**Fix:** `prepare_for` takes the `CandidateAlleles` alongside the view and asserts the two agree.

**M3: src/ng/calling/mod.rs:839 — `advance_expected_copies` hands back the pass-before-last's real values, un-poisoned**
**Categories:** errors. **Confidence:** Medium (the loop that could fail to write it does not exist yet; the gap in the type is certain).
Every other buffer is `NaN`-filled by `prepare_for`, and `fill_poisoned`'s own comment calls stale reuse "the bug this exists to avoid". This method reintroduces exactly that state and states the contract in prose — "must be written in full" — in a file whose `assert_agrees_with` exists to turn prose contracts into runtime facts. The cost reason given does not survive: the fill is `alleles` writes against a pass costing `samples × genotypes` evaluations each carrying an `lgamma`. A pass that skips an allele leaves it at its pass-*n−2* value, and the convergence test is a per-allele maximum change between two passes — so a stale entry near its neighbour lets the locus report `converged = true` on a number no pass wrote.
**Fix:** fill with the unwritten value on the way out; state a fact instead of a requirement.

**M4: src/ng/calling/mod.rs:553 — the empty contamination slice is a second, weaker spelling of a rule this file already imports a named one for**
**Categories:** defaults. **Confidence:** High.
`ContaminationMixture::new` (`likelihood/mod.rs:865`) **refuses** the empty spelling — *"one named way to say it, so that a caller reaches the decision rather than the shortest thing that compiles"* — and ships `uncontaminated()` and `is_absent()` for it, with a test pinning the refusal. `FrozenParameters::new` reproduces the emptiness test those three exist to replace and offers neither name. A consumer writing `contamination.get(rg).map(|v| v.fraction).unwrap_or(0.0)` silently turns *not estimable* into *estimated and clean*, which spec §3.6 of the read likelihoods refuses. This is the failure `arch/calling_em_loop.md` §2.1 names about the allele cap, on a different field. Nothing outside this file names `FrozenParameters` yet.
**Fix:** refuse the empty slice in `new`; add `uncontaminated()` and `contamination_is_absent()`.

**M5: src/ng/calling/mod.rs:679, :684, :796, :815 — the cohort's expected copies and one sample's share a name stem, a type and a length**
**Categories:** naming. **Confidence:** High.
`expected_copies()` is the cohort's; `sample_copies(i)` is one sample's. Both are `&[f64]`, both exactly `allele_count` long, and only one name carries its owner. The file's own `ExpectedAlleleCopies` doc says both quantities exist, that the prior's leave-one-out term is the difference of the two, and that "confusing them is a wrong prior rather than a compile error". Renaming is free today and costs a call-site audit once Milestones B and C land.
**Fix:** `cohort_` / `sample_` in every name, reusing the word `LocusInference::cohort_expected_copies` already uses.

**M6: src/ng/calling/mod.rs:935 — "this happens on the SNP/indel path only" is a doc-comment invariant with neither an assertion nor a test**
**Categories:** reliability. **Confidence:** High.
The spec states it twice (§5.0.1, §9). `LocusInference::new` enforces the mirror-image invariant one line up — `seed_diversity_unreachable` refused on a `Generic` locus, pinned by a test — but nothing refuses a `Missing` call at an `Ssr` or `SsrBundle` locus. Spec §9 says such a sample "has no posteriors and no expected copies at all", so a `Missing` wrongly emitted at a tract also removes that sample from the cohort's expected-copies denominator: wrong allele frequencies with a well-formed VCF record beside them.
**Fix:** the assertion beside its sibling, plus the test.

**M7: src/ng/calling/mod.rs:553 — the constructor refuses a run with no samples and accepts one with no read groups**
**Categories:** errors; reliability — **convergent**. **Confidence:** High.
`!inbreeding.is_empty()` is asserted on the reasoning that "a run has at least one sample"; `calibration` gets no such check, and an empty one *also* satisfies the contamination check, because `contamination.is_empty() || …` is true when both are empty. `read_group_count()` then returns 0 and every read-group lookup indexes an empty slice — an out-of-range index at whichever locus first carries a read, naming neither the read group, the locus, nor the axis.
**Fix:** mirror the sample-axis assertion.

**M8: src/ng/calling/mod.rs:389, :405 — the two empty-evidence messages name neither the locus nor the path, and are byte-identical**
**Categories:** errors. **Confidence:** High.
Both constructors carry the same text verbatim and neither interpolates the `region` it was just handed. Both facts are free — `region` is an argument and the constructor's name is the path — and both are the standard the rest of the file sets: `assert_agrees_with` interpolates the region, the path, the kind and both counts. On a caller committed to whole-genome loci and thousands of samples, a panic identifying no locus turns a one-line log into a re-run under a debugger.
**Fix:** interpolate the region and separate the two paths.

**M9: src/ng/calling/mod.rs:410 — `LocusEvidence::ssr`'s empty-list guard has no test in any build profile**
**Categories:** errors; reliability — **convergent**. **Confidence:** High.
The one test exercises `generic` only. Measured: of the 16 release-held assertions downgraded to `debug_assert!` and run under `--release`, this is **the only one** no test notices. The two constructors are a copy-paste pair, which is the shape in which one half drifts from the other — and the implementation report's test table lists the case once, as though it covered both.
**Fix:** the mirror test.

**M10: src/ng/calling/mod.rs:272 — `ExpectedAlleleCopies::new` documents refusing an infinite count; only `NaN` and negative are tested**
**Categories:** reliability. **Confidence:** High.
Weakening `is_finite()` to `!is_nan()` admits `+∞` and leaves the suite green. `+∞` is the arithmetic-gone-wrong shape a log-domain sum produces, it survives the constructor, rides into `LocusInference`, and the convergence comparison of two infinities is `NaN` — the failure the `NaN` check exists to stop, arriving through the untested door.
**Fix:** the third test.

**M11: src/ng/calling/mod.rs:580, :591, :599, :607 — the only positive `FrozenParameters` fixture has one read group and one sample, and two accessors are never read at all**
**Categories:** reliability (two findings, merged here). **Confidence:** High.
With both counts equal to 1, every wrong pairing of the two accessor bodies passes: `read_group_count()` returning `inbreeding.len()`, or `sample_count()` returning `calibration.len()`. `sample_count()` is what `assert_agrees_with` compares the evidence against, so a swap turns the run-order guard into a read-group-count guard — invisible on the tomato cohort, where each of 63 accessions has one library, and wrong on any run where a sample is sequenced twice. Separately, `calibration()` and `inbreeding()` are read by no test: making `inbreeding()` return `&self.inbreeding[..0]` leaves the suite green.
**Fix:** one fixture with one library and three samples, asserting all four.

**M12: src/ng/calling/mod.rs:667, :669, :675, :677 — four of the nine sized buffers have no accessor, so their sizing is unobservable**
**Categories:** reliability; errors and defaults (cross-category) — **convergent**. **Confidence:** High.
`prior_row`, `posterior_row`, `sample_concentration` and `prior_allele_scratch` are sized and poisoned and reachable by nothing. Sizing `posterior_row` with `allele_count` instead of `genotype_count` is therefore not merely untested but unobservable from outside the module — and at a diploid biallelic locus those are 3 and 2, so a swap is a legal length. That is precisely the hazard `prepare_for` takes a `GenotypeTableView` rather than two integers to prevent; the prevention is applied to the two buffers a test can see and not to the four it cannot. `#[derive(Debug)]` counts as a read of each field, which is why `dead_code` stays quiet.
**Fix:** give the four the accessor pairs the others have, and assert their lengths.

**M13: src/ng/calling/mod.rs:660, :858 — the three sub-scratches and the default type parameter are never instantiated**
**Categories:** reliability. **Confidence:** High.
`ssr_row_mut`, `generic_row_mut` and `selection_mut` have no call site anywhere in `src/`. All five `CallingScratch` values in the file are `CallingScratch::<()>::default()`, so `SsrRowScratch<StutterSubstitutionScratch>` — the whole subject of the implementation report's §2.3 — is never built. `#[derive(Default)]` on a generic struct is bounded by `SsrEmissionScratch: Default`, so the shipped configuration compiling is not implied by the tests compiling.
**Fix:** a test that builds the shipped configuration and touches all three.

**M14: no release gate holds any of this module's assertions — DEFERRED**
**Categories:** errors. **Confidence:** High.
The module's no-`Result` design rests on structural checks holding in release (spec §8), and the implementation report states it as an achieved property. The only test command in `.github/workflows/ci.yml:47` is a debug run, in which `assert!` and `debug_assert!` are indistinguishable; no `--release` invocation exists anywhere in `.github/workflows/`. So the property is held by the source text and by nothing else, and a later edit can downgrade any of these with a green CI.
**Why it is deferred rather than fixed:** the gate cannot be added today. `cargo test --release --lib ng::calling` at this commit is `461 passed; 4 failed`, and all four failures are pre-existing `#[should_panic]` tests in `src/ng/calling/likelihood/` — another branch's file, which this session must not edit. See §7.

### Minor

- **Mi1 — `lg_table` / `lg_row` put spec notation in the API, and the doc comment beside them names the wrong quantity.** (naming, High.) `Lg` appears nowhere else in `src/` and is expanded nowhere in code. Worse than the abbreviation: the comments call the buffer "the read likelihood", where `spec/read_likelihoods.md:29` reserves that for `Lr`, one read against one allele. What the buffer holds is the **genotype** likelihood. *Fix:* rename to the spec's own words and correct the three comments.
- **Mi2 — `FrozenParameters`' five fields name their topic, not the axis they are indexed by.** (naming, High.) The type's own doc says the two axes "are not interchangeable"; no name carries which. `parameters.inbreeding()[i]` gives the reader no way to tell which `i` is legal. *Fix:* `…_by_read_group` / `…_by_sample`, `seed` → `prior_seed`, `ssr_strata` → `ssr_slippage_fits`.
- **Mi3 — three words carry the scratch's doc comments and none is defined there:** *concentration* (defined one module down), *stratum*, *fold* (the functional-programming accumulate — borrowed jargon). (naming, High.) *Fix:* one clause each at first use.
- **Mi4 — `fill_poisoned` borrows a term of art, and the `NaN` it fills with is an unnamed literal in nine places.** (naming; defaults — **convergent**, High.) "Poison" in Rust means a lock left unusable by a panicking thread. The five accessors that hand a poisoned buffer out say nothing about its state; `advance_expected_copies` is the exception and shows what they lack. *Fix:* a named constant plus a plainly-named helper, and one doc line on each accessor.
- **Mi5 — `selection` / `selection_mut` is half a name.** (naming, High.) The crate's term is *candidate selection*. *Fix:* `candidate_selection_mut`.
- **Mi6 — `SampleGenotypeCall::Missing` names the VCF symbol, not the condition, and the case a geneticist would call "missing" is the other variant.** (naming, **Medium**.) A zero-coverage sample gets `GenericSampleEvidence::empty()`, is scored, and comes back `Called`. *Verification step attached by the sub-agent:* confirm that zero-coverage samples really are `Called`. They are. See §4 question 1 for the resolution.
- **Mi7 — the default type parameter picks one arm of the emission seam invisibly.** (defaults, High.) `SsrRowScratch<ModelScratch>`, the neighbour it wraps, gives the same parameter no default, so the convention in the seam is the opposite of the one adopted here. Softest of the set: a mismatch is a type error, not a wrong genotype. *Fix (sub-agent offered two):* keep the default and document it, or drop it and match the sibling.
- **Mi8 — `ssr_strata`'s "supply a gather over no outcomes" is prose with no expression a caller can copy.** (defaults, Medium.) The test module knows the incantation and hides it in a helper. *Fix:* name the expression in the doc bullet.
- **Mi9 — the two `checked_mul` `expect` messages name neither operand**, though both are locals one line up, and these are the only two panic sites in `prepare_for` that pathological real data can reach: `genotype_count` grows as `C(A + P − 1, P)`. (errors, High.)
- **Mi10 — `LocusInference::new` checks the copies vector against the allele table but never the called genotypes' allele ids.** (errors, Medium, **needs verification** — the fix was not compile-verified.) `Genotype::new` checks only non-emptiness, and `CandidateAlleles::admit`'s own comment records that the final prune "renumbers every id above it — so … every `AlleleId` minted before it goes stale". An out-of-range id reaches the VCF's `GT` as an index past the ALT list. The check costs `samples × ploidy` comparisons per locus.
- **Mi11 — the implementation report uses three plan step codes bare** (`E3`, `E2a`, `B1`) where it glosses three others (`A1`, `E1`, `C3`), plus "the three joins" and "the two fan-out plans". (naming, High.)
- **Mi12 — the implementation report's test table credits one test with "the positive path on both paths",** which reads as complete coverage of the discriminant; measured, it covers two of the three `LocusKind` variants. (reliability, High.) The third is **B1**.
- **Mi13 — three small coverage gaps** (reliability, merged): `region()` is asserted on the `Generic` arm only; `CandidateAlleles::is_empty` and `ExpectedAlleleCopies::is_empty` are documented "never true" and never called; and `CallingScratch` is documented "allocated once per worker", which requires `Send`, with nothing checking it.

### Nits

Grouped, not enumerated: `prepare_for` → `prepare_for_locus` (the preposition dangles at the call site); `assert_agrees_with` does not say what agreement is checked; `prior_allele_scratch` reads as "a prior allele" and `_scratch` inside `CallingScratch` says nothing the type has not; `row_range(sample, width)` takes two bare `usize`s positionally in the type whose `prepare_for` deliberately refuses to, and its message does not say which of the two flat tables was indexed; `assert_agrees_with` formats `LocusKind` with `{:?}`, which renders a tract's flanks as decimal byte arrays; test fixtures named as bare adjectives (`diploid()`, `outbred()`, `frozen()`); `advancing_makes_this_passs_copies_the_previous_passs` reads as two dropped apostrophes; *seam* is architecture-doc metaphor doing load-bearing work in a module doc; and `GenericLocusSample` names the pair by its coordinates where the value is evidence **plus** a ruling.

## 7. Out of scope observations

- **`src/ng/calling/likelihood/` holds its checks in debug, not release.** Four `#[should_panic]` tests there fail under `--release` at this commit: `a_sample_carrying_more_than_its_whole_batch_is_a_caller_bug`, `one_pair_appearing_twice_is_a_caller_bug`, `rows_out_of_pair_order_are_a_caller_bug`, `rows_out_of_read_group_order_within_one_allele_are_a_caller_bug`. Not introduced here and not this session's to fix — but **any release gate added for M14 has to deal with them first**. Follow-up: raise with the `ng-calling-likelihoods` branch.
- **`doc/devel/ng/arch/read_likelihoods.md` §2.1 still puts `genotype_must_be_missing` on `GenericSampleEvidence`,** which the shipped type does not carry. The implementation report's §2.2 and §6 record the divergence and say the amendment is the owner's call.
- **`doc/devel/ng/arch/calling_em_loop.md` §2** still sketches `CallingScratch` with public fields and a single `concentration: Vec<f64>`, and `SampleGenotypeCall` as a struct. Same disposition.
- **`SelectionScratch`'s doc comment** in `allele_candidates/mod.rs` says it "becomes a field of `CallingScratch` when that type exists". It now is one. Another branch's file.
- **`StratumFits` has no named empty constructor.** The right home for Mi8's rule is a `StratumFits::none_fitted()` in `parameter_estimation/joint/`, the way `ContaminationMixture::uncontaminated` works. Outside `src/ng/calling/`; routed as an API note.
- **The most dangerous join the spec names has no owner until E1.** Spec §5.0 warns that `LocusSelection::unmatched` is parallel to the merge's *covering* samples while `LocusEvidence::per_sample` is one entry per *run* sample. `LocusEvidence`'s doc correctly assigns the conversion to the input edge, which means it is checked nowhere in A1. Worth an explicit line on E1's checklist.

## 8. Missing tests to add now

| test | input class | bug it catches |
|---|---|---|
| `ssr_evidence_against_a_bundle_allele_table_is_accepted` | repeat evidence × `SsrBundle` alleles | a `paths_agree` arm that lost `SsrBundle`, refusing every bundle in the run |
| `generic_evidence_against_a_bundle_allele_table_is_refused` | SNP/indel evidence × `SsrBundle` alleles | a `Generic` arm widened to accept `SsrBundle` — a bundle scored by the generic row |
| `ssr_evidence_naming_no_sample_at_all_is_refused` | empty per-sample list, repeat path | the `ssr` constructor's guard being dropped |
| `expected_copies_reject_an_infinite_count` | `+∞` | `is_finite()` weakened to `!is_nan()` |
| `read_group_count_and_sample_count_are_two_different_axes` | one library, three samples | the two accessors reading each other's field; a truncated `inbreeding()` |
| `preparing_a_locus_poisons_the_previous_passs_copies_as_well` | two loci of one shape, previous-pass buffer written between | locus *n* declared converged on locus *n−1*'s copies |
| `a_repeat_tract_locus_cannot_carry_a_missing_call` | `Missing` at a tract | the SNP/indel ruling wired onto the wrong path (needs M6's assertion) |
| `frozen_parameters_refuse_a_run_with_no_read_groups` | empty calibration list | the read-group axis going missing (needs M7's assertion) |
| `the_shipped_emission_scratch_is_the_default_type_parameter` | the configuration the run actually uses | the shipped `Scratch` losing `Default`; the three dead sub-scratch accessors |
| `an_unprepared_scratch_is_refused` | `default()` without `prepare_for` | M1's silent empty slices |
| accessors + length assertions for the four unreachable buffers | — | `posterior_row` sized with `allele_count` (3 against 2 in the fixture) |

## 9. What's good

- **`prepare_for` taking a `GenotypeTableView` rather than two integers** is the right shape for the hazard: at a diploid biallelic locus the allele and genotype counts are 2 and 3, and passed as bare integers a swap leaves every buffer a legal length. (Two sub-agents cited it approvingly; M12 is that it was not applied to all nine buffers.)
- **The scratch fixture writes each sample's own index into that sample's row of both flat tables**, so a wrong window names the sample it came from. Verified by the reliability sub-agent: slicing `lg_row` with `allele_count` is killed by that test, and so are the `lg_row_mut` and `sample_copies` equivalents.
- **`preparing_a_locus_overwrites_the_previous_locus_of_the_same_shape` really does catch the `Vec::resize` failure** — verified by replacing `fill_poisoned`'s body with a bare `resize`.
- **`GenericLocusSample` pairing the evidence with the ruling rather than carrying two parallel slices** is the shape spec §5.0's closing paragraph asks for.
- **The diff honours arch §2.1's allele-cap decision**: no cap appears anywhere in `calling/mod.rs`, and `CandidateAlleles::admit` says so in its own doc rather than introducing one.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features --tests -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
- **New, introduced by this review:** `./scripts/dev.sh cargo test --release --lib ng::calling::tests` — the only run that can tell `assert!` from `debug_assert!`, and the one that established M9. It is green for `ng::calling::tests` at this commit; the four `likelihood/` failures of §7 appear only under the wider `ng::calling` filter.
