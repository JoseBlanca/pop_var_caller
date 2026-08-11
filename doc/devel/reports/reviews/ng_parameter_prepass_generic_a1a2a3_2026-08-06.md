# Code Review: ng_parameter_prepass_generic_a1a2a3
**Date:** 2026-08-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** commit `5b3a646` — step 4's parameter pre-pass, plan steps A1+A2+A3 (module tree, four constrained scalars, the error-rate ladder)
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** one commit's diff — `5b3a6468a5e84d787e5ba79d5ec9609a294b91c1` on `ng-parameter-estimation`.
- **Reviewed against:** that commit, checked out detached in six isolated worktrees.
- **In-scope files:** [src/ng/parameter_estimation/mod.rs](../../../../src/ng/parameter_estimation/mod.rs), [fitting/mod.rs](../../../../src/ng/parameter_estimation/fitting/mod.rs), [fitting/mixture_weights.rs](../../../../src/ng/parameter_estimation/fitting/mixture_weights.rs), [generic/mod.rs](../../../../src/ng/parameter_estimation/generic/mod.rs), [generic/depth_and_alt_reads.rs](../../../../src/ng/parameter_estimation/generic/depth_and_alt_reads.rs), [generic/histogram.rs](../../../../src/ng/parameter_estimation/generic/histogram.rs), [generic/runs.rs](../../../../src/ng/parameter_estimation/generic/runs.rs), [src/ng/types.rs](../../../../src/ng/types.rs), [src/ng/mod.rs](../../../../src/ng/mod.rs), and the step's implementation report.
- **Deliberately out of scope:** the preceding housekeeping commit `ce3f0b4`; Milestones B–G, unwritten by design.
- **Categories dispatched:** `reliability` (always), `errors` (always), `naming` (always — and the category this project's history says is highest-yield on a docs-heavy commit), `module_structure` (the commit creates a module tree), `idiomatic` + `defaults` (constants and a public builder), `smells` + `refactor_safety` (six near-identical constructors). `unsafe_concurrency` skipped — no `unsafe`, `Arc`, lock, atomic, channel or thread. `tooling` skipped — `Cargo.toml` untouched.

### 2. Verdict

**Approve-with-changes.** No Blocker. One Major that is an untested rejection direction rather than a live defect, three Majors of documentation accuracy where a stated reason is checkably false, and a set of Minors concentrated on one theme: **the ladder builder coerces where the newtypes beside it reject.**

### 3. Execution status

Run in the container before dispatch, and quoted verbatim into every sub-agent prompt:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::types::` | 0 | 26 passed |
| `cargo test --lib ng::parameter_estimation` | 0 | 3 passed |
| `cargo test --all-targets --all-features` | 101 | 2901 passed, 1 failed, 5 ignored |
| `cargo doc --no-deps --lib` | 101 | 15 unresolved links — 12 pre-existing, 3 introduced here |

**`cargo doc` was not in the author's validation set and should have been.** `Cargo.toml` sets `broken_intra_doc_links = "deny"`, and `fmt`/`clippy`/`test` do not cover rustdoc. Two categories found the three new errors independently.

The single test failure is
`ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes`
(seed `0x5eed0001` case 18, `record_widen_events` 4 against production's 3), confirmed pre-existing at `HEAD` by stashing every uncommitted change. Out of scope; reported to the owner.

Findings labelled "Needs verification": **zero.** Every behavioural claim below was produced by mutating the code in an isolated worktree and re-running, with the output quoted in the per-category file.

### 4. Open questions and assumptions

1. **Does `ErrorRate` sharing `DomainError::ErrorRate` with the two aligner constructors need fixing?** (affects M3, Mi7) — three constructors in two modules render one identical message. Splitting the variant touches `alignment/`, outside this step. Deferred as an owner decision.
2. **Is `DomainError` still the right shape now that it carries six variants?** (affects the errors nits) — one enum per crate against "one error type per fallible operation". A deliberate ng-wide convention; worth a decision now the enum has doubled, not a change in this commit.

### 5. Top 3 priorities

1. **M1** — two of the three rate constructors have a rejection direction no test reaches, and widening it leaves the suite green. `InbreedingF` accepting 1.5 is the live one: its own doc says a user supplies it on the command line.
2. **Mi3** — the ladder's rung count comes from a saturating `as u32` cast with no ordering assertion. Swapped constants give a silent **one-rung** ladder, which sets the endpoint-argmax flag for every read group and destroys the one bit that distinguishes a railed fit from a plausible number.
3. **M4** — `generic/histogram.rs`, the first doc a Milestone-B implementer reads, describes a two-number key and says "nothing is lost". That is the key the design rejects for a multi-library sample, and it contradicts its own sibling file.

### 6. Findings

#### Major

**M1: src/ng/types.rs:313,338 — each `[0, 1]` constructor has one rejection direction no test reaches**
**Categories:** reliability
**Confidence:** High.
`each_constrained_rate_rejects_out_of_range_with_its_own_variant` passes `-0.01` and `1.01` to `ErrorRate` but only `1.5` to `GenotypeFrequency` and only `-0.5` to `InbreedingF`. Widening `GenotypeFrequency` to `(-1.0..=1.0)` and `InbreedingF` to `(0.0..=2.0)` in a worktree left `26 passed; 0 failed`. These four constructors are the only enforcement; nothing downstream re-checks, so a widened bound is an `F = 1.5` reaching step 8's genotype prior and producing a wrong genotype rather than an error.
**Fix:** a proptest over the whole `f64` line with a dense `-2.0..3.0` arm — `f64::ANY` alone never lands in `(1, 2]` and passed against the mutant — plus a `to_bits()` round-trip, which also closes the silent-clamp class.

**M2: src/ng/types.rs:357 — `Ploidy`'s doc contrasts it with "the unchecked newtypes above", which are all checked**
**Categories:** naming
**Confidence:** High.
The four newtypes directly above `Ploidy` — `MismatchFraction`, `ErrorRate`, `GenotypeFrequency`, `InbreedingF` — all have private fields and checked constructors. The arch doc says "elsewhere in `types.rs`", which is true; "above" was introduced in this commit and never was. A reader checking the claim stops trusting the paragraph, which is where the reason `Ploidy` rejects zero lives.

**M3: src/ng/types.rs:646-648 — the stated reason for the `is_finite` guard is false, and the guard is dead**
**Categories:** naming, errors, idiomatic, smells, reliability — **five of six, convergent**
**Confidence:** High.
`RangeInclusive::contains` is `start <= x && x <= end`, so it already rejects `NaN`, `+∞` **and** `-∞`. Four agents independently deleted `!x.is_finite() ||` from the three constructors and re-ran: `26 passed; 0 failed`, including the test that exists to pin the behaviour. The test's doc comment and the implementation report both claim `is_finite` is what catches `-∞`. The same file states the correct reason 100 lines earlier on `MismatchFraction`, so two comments in one file give contradictory accounts of the same predicate.

**M4: src/ng/parameter_estimation/generic/histogram.rs:5-9 — "nothing is lost" describes the key the design rejects**
**Categories:** naming
**Confidence:** High.
The doc says a site "reduces to two numbers … and nothing is lost". The arch says the opposite for the case the table exists to handle: with the library forgotten, the likelihood is exactly flat along every combination holding the share-weighted mean rate fixed, and no amount of genome separates the libraries — which is why `SiteKey` has an `Attributed` arm at all. The sibling `generic/mod.rs` states it correctly, so the two placeholder docs disagree.

#### Minor

**Mi1: src/ng/parameter_estimation/mod.rs:30-102 — `WindowIndex` and the ladder sit one level above their only consumer**
**Categories:** module_structure. **Confidence:** High.
The recorded justification — that `fitting/` reads the ladder — is contradicted by the design: `fit_by_profile_scan` takes the ladder as a parameter, and `fitting/mod.rs`'s own doc says it knows nothing about markers or windows. `WindowIndex`'s two consumers are both under `generic/`. `parameter_estimation/mod.rs` is the level the STR sub-unit will share, so an STR fit would inherit a window size and a ladder of per-base error rates it has no use for. The agent moved both blocks to `generic/mod.rs` in its worktree: no other `use` line in the crate changed, 3 tests green.

**Mi2: src/ng/types.rs:262 — the banner states consumers that do not exist, in the present tense, and miscounts what step 4 fits**
**Categories:** module_structure, naming. **Confidence:** High.
"the four constrained scalars step 4 fits and steps 7 and 8 consume" — steps 7 and 8 are not in the tree, and `Ploidy` is an *input* to the fits, as the same file says 90 lines later. The placement is right and should stay; the argument for it is stronger than the one written (defining them in `parameter_estimation/` would make the later steps import from a sibling stage module).

**Mi3: src/ng/parameter_estimation/mod.rs:90 — the rung count comes from a saturating cast, and two ways of mis-stating the ladder fail silently**
**Categories:** idiomatic, errors, refactor_safety. **Confidence:** High.
Measured: `MAX = 50.1` leaves the count at 161 with the top rung still at Phred 50, so the constant named `MAX` is not the ladder's top. Swapping `MIN` and `MAX` gives a negative `f32` that saturates to `0` on the `as u32` cast, producing a **one-rung** ladder; a step of `0.0` gives 4,294,967,295 rungs. Only `assert_eq!(ladder.len(), 161)` catches the first two — and `the_error_rate_ladder_rungs_are_a_constant_ratio_apart` **passed on the collapsed ladder**, because `windows(2)` over one element iterates zero times. The test that checks the shape is silent exactly when the shape collapses.

**Mi4: src/ng/parameter_estimation/mod.rs:70-74 — the two ladder-edge constants restate their value instead of stating the choice**
**Categories:** defaults. **Confidence:** High.
`INBREEDING_WINDOW_BP` is the model of how to do it — units in the type, "fixed, not a knob", the reason, the spec citation — and `ERROR_RATE_LADDER_STEP_PHRED` matches. The two edges get one line each that converts Phred to a probability and stops. Both sources exist and neither is cited: Phred 10–50 is DRAGstr's own grid range, and the architecture says what happens to a read group outside it.

**Mi5: src/ng/types.rs:288-341 — five copies of one `[0, 1]` predicate; extract a helper**
**Categories:** smells. **Confidence:** High.
The commit defends the repetition on the grounds that each constructor names its own `DomainError` variant. That fails for the type the commit is centred on: `ErrorRate::try_new` returns the *same* variant `FlatEmission::try_new` and `SsrSequenceMarginal::try_new` already return, so the argument covers two of four. And naming a variant never required duplicating the predicate — a tuple-variant constructor is a `fn(f64) -> DomainError` and can be passed as a value. The agent wrote and ran the replacement: 26/26 pass, clippy clean, no macro.

**Mi6: src/ng/types.rs:385 — `DomainError` derives `PartialEq` over `f64`, so a `NaN` rejection is not equal to itself**
**Categories:** refactor_safety, errors, reliability. **Confidence:** High.
Demonstrated: `assert_eq!(ErrorRate::try_new(f64::NAN), Err(DomainError::ErrorRate(f64::NAN)))` fails printing `left: Err(ErrorRate(NaN))` / `right: Err(ErrorRate(NaN))` — two sides that render identically. The existing tests dodge it with `is_err()`, but nothing records why. A `NaN` input is precisely what these constructors exist to reject, so this is not a corner case.

**Mi7: src/ng/types.rs:390-393 — `DomainError::ErrorRate` has three producers and one message**
**Categories:** errors, smells. **Confidence:** High.
`ErrorRate::try_new`, `FlatEmission::try_new` and `SsrSequenceMarginal::try_new` all render `per-base error rate {0} is not a finite probability in [0, 1]`. A bad rate in the step-4 pre-pass and one in the aligner's setup are indistinguishable in a log.

**Mi8: src/ng/types.rs:366 — `Ploidy::try_new` is never offered a large copy number**
**Categories:** reliability. **Confidence:** High.
The test is named "accepts every real copy number" and checks three of 255. Mutating the guard to `copies == 0 || copies > 8` left the suite green. Polyploids are in scope — the type's own doc says ploidy varies by region — so a ceiling slipped in later would reject a legitimate hexaploid region.

**Mi9: src/ng/parameter_estimation/mod.rs:90 — the ladder is pinned by literals, not by the constants it is documented to be derived from**
**Categories:** reliability. **Confidence:** High.
Both tests assert `161`, `0.1`, `1e-5` and a written-out `10^0.025`; nothing checks that the top rung *is* `ERROR_RATE_LADDER_MAX_PHRED`. Replacing `.round()` with bare truncation leaves all three green, so the rounding is untested.

**Mi10: the implementation report's validation table is wrong in two places**
**Categories:** reliability, naming. **Confidence:** High.
The full-suite row reads "2,904 passed, 1 failed" where the run says `2901 passed; 1 failed; 5 ignored` — 2904 is the "filtered out" figure from a different row, and the ignored count is omitted. Separately, "the first thing that reads a number is Milestone D" misses Milestone B: `add_site` derives a bin and sums depths, and B3's `mean_depth_in_cell` is a mean, isolated in its own commit precisely because getting it wrong lands the fit 5.2 rungs off.

**Mi11–Mi16: documentation claims that do not survive checking.** Each **Category:** naming, **Confidence:** High.
- "**Phred appears here and nowhere else**" is refutable with one `grep`: `types.rs` has `MapQual` and `BaseQual`, both Phred-scaled. True of step 4, not of the crate.
- "**coarsest**"/"**finest**" rung for the ladder's endpoints: the rungs are evenly spaced, so no rung is coarser than another. What varies is the error rate — Phred 10 is the *noisiest* read group. The test carried the metaphor into bindings named `coarser`/`finer` over two probabilities.
- "**Finer than a caller can feel**" is asserted where the cited spec section explicitly marks that number soft and unmeasured — and it is the one number the spec singles out as worth checking first.
- `INBREEDING_WINDOW_BP`'s "chosen against the runs a genome actually carries, **not against a run's data**" — the contrast has no referent anywhere in the file or the spec.
- The module doc attributes the pure-reference skip to **production's caller**; the spec says the record is still staged and written and the skip is the heterozygosity accumulator's alone. The spec put that sentence there to stop exactly this misreading.
- `mixture_weights.rs` never uses the words "mixture" or "weight"; `runs.rs` never uses the word "run". Both are the design docs' terms, so a reader arriving from a design doc and one arriving from the file tree get two vocabularies with nothing tying them together.

#### Nits

Three unresolved intra-doc links at `parameter_estimation/mod.rs:12-14`, where the same references are spelled as plain code spans in the six sibling files. No `#[must_use]` on `error_rate_ladder()`, which ng uses elsewhere for exactly this shape. `WindowIndex`'s doc compares it to `ContigId` to justify *not* being in the shared vocabulary — but `ContigId` is in the shared vocabulary, so the comparison supports the opposite of the clause attached to it. "unchecked" at `types.rs:357` against "unconstrained" everywhere else, for one concept. `top_rung` holds a rung *index*. `histogram.rs`'s worked example uses numbers that do not match the design-doc table it cites. `fitting/mod.rs` is the only placeholder that does not name the milestone that will fill it. `.expect()` renders `Debug`, so the carefully worded `Display` message never appears in the panic. `WindowIndex::get()` is public and untested. The `ng/mod.rs` re-export list has no membership criterion — 10 of 14 public items, omitting `ReadGroupId` and `GenomeRegion`, which the same design doc has step 4 reading (pre-existing, extended here).

### 7. Out of scope observations

- **`ng::locus_generation::pileup::parity`** fails at `HEAD` independently of this work (§3). A bug report against locus generation.
- **`FlatEmission::try_new`** (`alignment/emission.rs:288`) still takes a raw `f64` and re-implements the invariant `ErrorRate` now owns. A natural first adoption of the new newtype, in a later change.
- **`cargo doc` has 12 pre-existing unresolved links** across `ssr.rs`, `region_typing/`, `tandem_repeat.rs`, `types.rs` and `var_calling/`, so the doc build is red independently of this commit.

### 8. Missing tests to add now

1. `the_constrained_rates_accept_exactly_the_probabilities_and_round_trip` — the whole `f64` line with a dense `-2.0..3.0` arm plus `f64::ANY`. Catches a widened or flipped bound in either untested direction (M1) and any silent clamp of an accepted value. Verified failing at `x = 1.012771917376405` against the widened `InbreedingF`.
2. `ploidy_accepts_every_non_zero_copy_number_and_round_trips` — all 256 `u8` values. Catches an upper bound added later (Mi8); verified failing at `copies = 9` against the `copies > 8` mutant.
3. `the_error_rate_ladder_ends_at_the_phred_constants_it_is_built_from` — first rung, last rung and count read through the constants rather than literals. Catches a ladder whose top rung is not `MAX` (Mi9); verified reporting `last rung 1.0232929922807538e-5 vs 0.00001` at a step of 0.3.
4. `the_ladder_constants_divide_into_a_whole_number_of_rungs` — the precondition `.round()` silently relies on. Verified reporting `got 133.3333280351429` at a step of 0.3.
5. A non-vacuity assertion in `the_error_rate_ladder_rungs_are_a_constant_ratio_apart`, so a collapsed ladder cannot pass by iterating nothing (Mi3).

### 9. What's good

- **`INBREEDING_WINDOW_BP` is the reference implementation of a documented constant** in this codebase: units in the type, the "fixed, not a knob" decision stated, the spec section cited, and the value pinned by a test whose doc says why the number is load-bearing. Two categories cited it as the standard the ladder constants should meet.
- **The four-types-not-one-`Probability` decision holds up**: `module_structure` confirmed the placement in `types.rs` independently and supplied a better argument for it than the commit's own.
- **The deviation from A1 — no empty `#[cfg(test)]` blocks in the five documentation-only files — was explicitly endorsed** by the reliability agent as the right call rather than a gap.
- **`Ploidy` deriving `Ord` where the three `f64` rates derive only `PartialOrd`** is the correct distinction, and `ploidy_orders_by_copy_number` pins it.
- **Every `expect`/`unwrap` audit came back clean**: exactly one `expect` in the whole in-scope surface, carrying its `// PANIC-FREE:` comment, and proven unreachable for every rung the ladder produces.

### 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `./scripts/dev.sh cargo doc --no-deps --lib` — **new to this feature's validation set**; it is what surfaced the intra-doc links, and `Cargo.toml` denies them.

Per-category files kept as an audit trail in `tmp/review_2026-08-06_ng-param-A1A2A3/`.
