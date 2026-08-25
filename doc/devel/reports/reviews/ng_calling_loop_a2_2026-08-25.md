# Code Review: ng calling loop — A2, the calling seam and its configuration

**Date:** 2026-08-25
**Reviewer:** rust-code-review skill (orchestrator), five category sub-agents in isolated worktrees
**Scope:** the diff of step A2 of [calling_loop.md](../../ng/impl_plan/calling_loop.md), captured as commit `a2d28e62`
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** a diff — `git diff a2d28e62~1 a2d28e62`, 931 insertions across three files.
- **Reviewed against:** `a2d28e62cad4c33ebf69f7a7c5f8fd130b0084d9` on branch `ng-calling-loop`, kept alive as `refs/review/ng-calling-loop-a2`.
- **In-scope files:** [src/ng/calling/inference/mod.rs](../../../../src/ng/calling/inference/mod.rs) (new), two hunks in [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs), and [the step's implementation report](../implementations/ng_calling_loop_a2_2026-08-25.md).
- **Out of scope:** `src/ng/calling/likelihood/` and `src/ng/calling/allele_candidates/`, owned by two other branches — consumed, never edited. The previous step's types in `calling/mod.rs`. Everything outside `src/ng/calling/`.
- **Categories dispatched:** `defaults` (the diff is nine constants, five `Default` impls and a validator), `errors` (it introduces the folder's only `Result`), `reliability` (all of it is validation logic), `naming` (about thirty new public names), and **`idiomatic`** — dispatched deliberately because A1's review had to trim it and recorded that as a coverage gap, and because this diff is the branch's first public API surface.
- **Not dispatched:** `module_structure`, `unsafe_concurrency`, `tooling`, `smells`, `refactor_safety`, `extras`.

## 2. Verdict

**Request-changes.** Three Majors, of which two were demonstrated by building the failing case rather than argued. The code was correct on every path its tests reached; what the review found is a value that passes validation and makes a whole run lie, a check nothing forces anyone to call, and a set of tests that could not fail for a third of the values they existed to pin.

## 3. Execution status

Run by the orchestrator in the container from the reviewed tree, at main's 1.98 compiler pin:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | no warnings — **the wide scope**, available since this branch took main's pin |
| `cargo test --lib` | 0 | `4523 passed; 0 failed; 14 ignored` (A1 left 4,517) |
| `cargo test --release --lib ng::calling::inference` | 0 | `6 passed; 0 failed` |

**Mutation testing, reliability category: 16 run, 8 survived, 1 of those changed no behaviour** — so seven real gaps. Each survivor's behaviour change was demonstrated by diffing a 72-case `validate()` fingerprint against the clean tree, not asserted. That worktree's source was restored and **verified by content** against `git show a2d28e62:…` before the findings were written.

**Two sub-agents also built failing cases rather than reasoning about them**, and both are quoted in §6: a refused configuration handed to `call_locus` through a `Box<dyn LocusGenotyper<()>>`, and a convergence threshold of `1.01` and `50.0` returning `Ok(())`.

Findings labelled "Needs verification": **0**.

**One process defect, recorded rather than quietly fixed.** The `defaults` sub-agent's first probe ran a `cd` that landed in the **orchestrator's** checkout rather than its own worktree, so that probe was compiled and run in the tree under review. The agent reverted it and re-checked by md5 and `git status`; the orchestrator independently confirmed the working tree was byte-identical to `a2d28e62` (`git diff --stat refs/review/ng-calling-loop-a2` empty) before any fix was applied. No result in this review rests on that tree's state at the time.

## 4. Open questions and assumptions

1. **Should the discovery bar reuse the merge's `MinAltReads`, or state its own two primitives?** Affects **Mi3**. Resolved in favour of reuse: the type is denominator-agnostic by construction — its own doc says *"the numerator is the caller's, and the callers count different reads"* — so sharing it costs nothing and makes a negative share unrepresentable. The **numbers** stay separate, because spec §4.1's third open question sweeps discovery's independently.
2. **Validated newtype, or a release-held assertion at the seam?** Affects **M2**. Two sub-agents proposed different fixes for the same finding. Resolved in favour of the newtype, because `idiomatic` built it and measured the cost: 53 insertions / 12 deletions, one call site.
3. **Does the seam need an error channel?** Not filed as a finding. `arch/calling_em_loop.md` §3.2 requires arms B and C to *reject* a non-zero `InbreedingF` rather than ignore it, and `call_locus` returns no `Result`. Those arms are the bake-offs plan's; recorded in §7.

## 5. Top 3 priorities

1. **M1** — a convergence threshold of 1.0 passes validation, stops every locus in the run after one pass, and flags each as settled.
2. **M2** — a configuration `validate()` refuses reaches `call_locus` unimpeded; nothing forces the check that is the entire point of the step.
3. **M3** — three of the nine shipped values are not pinned by the test whose job is pinning them; two are compared against the constants they were built from.

## 6. Findings

### Major

**M1: inference/mod.rs:339 — the convergence threshold is bounded below and not above, and 1.0 is a silent lie**
**Categories:** errors, defaults, reliability — **convergent, all three**. **Confidence:** High.
The check is `is_finite() && > 0.0`. The quantity compared against it is a per-allele change **already divided by the cohort's chromosomes** (spec §6), so it lies within `[0, 1]` and cannot exceed 1 — which makes any threshold at or above 1 satisfied by a loop that has done nothing. Measured on the unmodified tree: `conv=1.01 -> Ok(())`, `conv=50.0 -> Ok(())`. Production validates the same field within `(0, 0.1]` and its constant says why in as many words — *"a `1.0` threshold would always exit after one EM iteration"* — and spec §6 quotes production's reason for the division as being that the threshold *"and its validation range"* carry over unchanged in meaning. **The lower half was inherited and the upper half was dropped.** The tell is internal: the discovery share, twenty lines away in the same function, is bounded at both ends.
**Consequence:** every genotype in the run is the initial pass's, emitted with `converged = true` — and that flag exists precisely so that a genotype from a loop that did not settle is a different claim from one that did.
**Fix:** restore `CONVERGENCE_THRESHOLD_RANGE_MAX = 0.1` with the reason in place, and name the range in the message.

**M2: inference/mod.rs:338, :486 — a configuration `validate()` refuses reaches `call_locus` unimpeded**
**Categories:** errors, idiomatic — **convergent**. **Confidence:** High.
`CallingLoopConfig` has four public fields, a `Default`, and a `validate()` nobody is obliged to call; the seam takes `&CallingLoopConfig`. `validate()` has no call site outside the file's own tests. **Demonstrated:** a configuration with `convergence_threshold = NaN`, `slippage_refit.max_rounds = 3` and `discovery.mode = AgainstFullConvergence` is refused by `validate()` and then handed to `call_locus` through `Box<dyn LocusGenotyper<()>>` — it compiles and runs, with no refusal anywhere.
The argument that settles it is in-folder: the sibling `CandidateSelectionConfig` needs no validation **at all**, because both its fields are newtypes that cannot hold an illegal value; and `SsrSegment`, the house's error-shape precedent, keeps its fields private behind a fallible `new` precisely so *"the accessors are infallible by construction"*.
**Consequence:** the plan's rule for this step — an unbuilt setting *"rejected loudly at config validation, never silently ignored"* — rests on a caller remembering a method call, in a folder whose whole policy is that wrong-answer-without-a-crash failures get held checks rather than conventions.
**Fix:** `validate(self) -> Result<RunnableCallingLoopConfig, _>`; the seam takes the token. Measured by the `idiomatic` sub-agent on its own build: 53 insertions / 12 deletions in one file, exactly one call site touched, tests green.

**M3: inference/mod.rs:521 — three of nine shipped values are not pinned by the test that exists to pin them**
**Categories:** reliability. **Confidence:** High (three surviving mutations).
`the_shipped_configuration_is_one_this_caller_will_run` asserts five constants against **literals**, and those are discriminating — mutating `DEFAULT_MAX_PASSES` or the shape pull-back kills the test. But `slippage_refit.round_convergence_threshold` is never read by the test at all, and the discovery share and round cap are compared against the very constants `Default` read a moment earlier — **an identity that holds for whatever value the constant is edited to**. Measured: setting the round threshold to `0.5`, the share to `0.9` and the discovery round cap to `1` each leaves the whole suite green.
**Corroborating, and it points at the same constant:** the implementation report says "eight `DEFAULT_*` constants" where the file has nine, and the one that dropped out of that count is exactly the one with no assertion.
**Consequence:** these are the numbers the design's open questions Q2–Q4 will move deliberately, and a deliberate move should have to edit a test that states the old value.
**Fix:** compare against literals throughout, and add the missing field.

### Minor

- **Mi1 — `PullBackOutOfRange { which: &'static str }` discriminates two closed cases with a string.** (naming **Major**, errors, idiomatic — **convergent, three categories**, High.) A public variant of a `#[non_exhaustive]` enum, whose two values are matched back by string literal in the tests: a typo compiles and the `matches!` simply stops matching. `DiscoveryMode` two lines below is carried as a type. *Fix:* a two-variant enum with `Display`; the message text is unchanged.
- **Mi2 — three counts are typed wider than their domains.** (idiomatic, High.) `max_passes`, the discovery bar's `min_reads` and `DiscoveryConfig::max_rounds` are all documented as at least 1. Two are defended at run time by an error variant each; **the third is not defended at all** — a switched-on discovery loop with a zero cap runs no rounds and reports finding nothing, which is the failure `DiscoveryConfig::default()` is written out to prevent. `NonZeroU32` deletes the two variants and closes the third hole with no check. The house already does this next door (`MinAltObs(pub NonZeroU32)`).
- **Mi3 — the discovery bar restates a rule the merge already ships.** (defaults Medium, idiomatic cross-category — convergent.) A read floor plus a share of one sample's own reads, both of which must clear, is `MinAltReads`. Holding it as a bare `u32` and `f64` makes `validate()` re-derive both range checks by hand and costs two error variants. This is the argument the module doc already makes against holding the allele cap here, applied to the bar.
- **Mi4 — `shape` is the wrong word and it is taken.** (naming, High.) `shape_pull_back_pseudocounts` pulls back two numbers the crate already names `Slippage::shorter_share` and `Slippage::fall_off` — in the very struct `StratumFits` hands the calling path — so the file invents a third spelling for one of them while using the crate's word for the other. And `LocusShape` already means a read histogram in `parameter_estimation/ssr/`. **`level` survives**: it is exactly `Slippage::level`, and its constant glosses it in the first clause.
- **Mi5 — `DiscoveryBar` and `min_spanning_read_share` do not say whose reads.** (naming, High.) Spec §4.1's rule is a share of *one sample's* spanning reads; HipSTR's other admission route in the same table is cohort-wide, so both readings are live and only the doc comment separates them.
- **Mi6 — the check-ordering rule is pinned by one fixture, and the order between the two not-built refusals is unpinned.** (reliability, mutation-proven.) Swapping them changes what a run is told while leaving every test green — and the both-on case is the one a measurement harness hits first, since a bake-off arm sets both.
- **Mi7 — the value an error carries and the message it prints are never asserted.** (reliability, High.) Every test reaching an `f64`-carrying variant matches it with `{ .. }`, which reads no payload, so an implementation answering every refusal with `ConvergenceThresholdOutOfRange { threshold: 0.0 }` passes the suite. No test renders a message at all.
- **Mi8 — the finiteness half of three of the four `f64` checks is untested.** (reliability, mutation-proven.) `is_finite()` weakened to `!is_nan()` still refuses a not-a-number and admits an infinity, so a not-a-number case alone leaves the weakening alive. It survived on both pull-backs and on the round threshold.
- **Mi9 — seven of nine fields do not name the constant that fills them.** (defaults, High.) Exactly one does. `DiscoveryConfig::max_rounds` matters most: `4` appears nowhere a caller reading that struct would meet it.
- **Mi10 — `DEFAULT_ROUND_CONVERGENCE_THRESHOLD` states the provenance of its rule and not of its value.** (defaults, High.) `1e-3` is production's, used twice — one threshold on the shape coefficients and one on the level multiplier. That ng collapses the two into one number is unremarked.
- **Mi11 — the stand-in implementation hard-codes ploidy 2 and a two-entry copies vector**, against a ploidy-generic sibling line. (reliability.)
- **Mi12 — `name()` is asserted only to be non-empty.** (reliability.)
- **Mi13 — the implementation report's counts are wrong in two places.** (naming, defaults — convergent.) It says "eight `DEFAULT_*` constants" (there are nine) and the error enum's doc says "three variants carry the offending `f64`" (four do). Its §2.1 also says an associated type "would have worked equally" for the emission scratch; it does not compile — `error[E0207]`.

### Nits

Grouped: `DiscoveryNotBuilt` renders the mode with `{mode:?}`, putting a Rust identifier in a log line; none of the eight messages names the field an operator actually set; `max_passes` has no upper bound where production caps its analogue at 500; `share.is_finite()` is redundant beside a `[0, 1]` range test; the four configuration structs are not `#[non_exhaustive]` where the error enum is; `CallingLoopConfig` has a `Default` but no `pub const DEFAULT: Self`, where the sibling in the same folder has both; and the doc's *"report it as the re-fitted arm"* uses **arm** — a word appearing on eleven lines of the file and defined on none — in the two places that reach an operator rather than a reader of the source.

## 7. Out of scope observations

- **`call_locus` has no error channel, and `arch/calling_em_loop.md` §3.2 needs one.** That document requires arms B and C to *reject* a non-zero `InbreedingF` rather than ignore it — *"a silently dropped `F` on a selfing panel is the failure this is guarding"* — and the seam as built cannot express a rejection. Those arms belong to the bake-offs plan; raise it there before either is written.
- **`StratumFits` still has no named empty constructor**, carried over from A1's review.
- The **naming** sub-agent notes that `LocusEvidence::Generic` inherits `Generic` from `LocusKind`, which every doc comment in both steps has to gloss as "SNP/indel" — a crate-wide decision, not this step's.

## 8. Missing tests to add now

| test | input class | bug it catches |
|---|---|---|
| literal assertions for the three unpinned shipped values | the shipped configuration | a constant edited without a test noticing |
| infinity cases on both pull-backs and the round threshold | `+∞`, `−∞` | `is_finite()` weakened to `!is_nan()` |
| `both_unbuilt_settings_at_once_report_the_slippage_refit` | both outer loops on | the two refusals swapped, changing what a run is told |
| `a_refused_value_travels_in_the_error_and_into_its_message` | a refused value | every refusal answered with the same payload |
| the ordering generalised over every range check | each range failure × both unbuilt | a not-built check hoisted above one particular range check |
| `discovery_and_candidate_selection_share_a_rule_and_not_its_numbers` | the two bars | a sweep of one silently moving the other |

## 9. What's good

- **The two not-built refusals say what accepting them would cost**, not merely "unsupported" — the `errors` sub-agent's verdict was that this earns its length, because it tells an operator why they cannot shrug and proceed.
- **The range-before-not-built ordering is right and is discoverable**, with a comment at the branch and a test of its own. Checked exhaustively: all seven range checks are unconditional and every one precedes both refusals.
- **`candidates` by value is the right signature.** `LocusInference` stores the table by value, so `&mut` would force a per-locus clone of a `Vec<Box<[u8]>>` at a seam whose scratch exists because allocation was measured at about one cycle in six.
- **The type parameter on the trait is forced, not preferred** — an associated type is `error[E0207]` here, and every implementor would need its own parameter plus `PhantomData`.
- **`DEFAULT_DISCOVERY_MAX_ROUNDS`'s doc is the strongest of the nine on the question it needs to answer**: it says the number is this step's own, that it is soft, what it is twice, and which plan should set it.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --release --lib ng::calling::inference`
