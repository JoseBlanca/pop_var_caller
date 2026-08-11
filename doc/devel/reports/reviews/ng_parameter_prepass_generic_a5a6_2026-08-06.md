# Code Review: ng_parameter_prepass_generic_a5a6
**Date:** 2026-08-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** commit `b9ef1e8` — step 4's parameter pre-pass, plan steps A5 and A6 (the output types, the step's error, the fit floors)
**Status:** Request-changes → resolved

---

### 1. Scope

- **What was reviewed:** one commit's diff across four files in `src/ng/parameter_estimation/` plus an `impl Display for Ploidy` in `src/ng/types.rs`. Types and constants only; nothing computes.
- **Reviewed against:** `b9ef1e8ec1949aec9326c0c281b8bfed1991c553`, checked out detached in three isolated worktrees.
- **Categories dispatched:** `reliability` + `errors`; `naming` + `defaults` + `module_structure`; `idiomatic` + `smells` + `refactor_safety`. `unsafe_concurrency` and `tooling` skipped — no `unsafe`, no concurrency primitive, `Cargo.toml` untouched.
- **Out of scope:** `depth_bins.rs` and the error-rate ladder (reviewed at A4); Milestones B–G.

### 2. Verdict

**Request-changes.** Six Majors. One is a wrong number a caller would read; two are wrong numbers in doc comments, one of which traces upstream into the spec and the architecture; one is in the implementation report; and two are gaps where a mutation left the suite green. All resolved.

### 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::parameter_estimation` | 0 | 34 passed |
| `cargo doc --no-deps --lib` | 101 | 12 unresolved links, all pre-existing, none in these files |
| `cargo test --all-targets --all-features` | 101 | one pre-existing failure in `ng::locus_generation::pileup::parity` |

Findings labelled "Needs verification": **zero**.

### 4. Open questions and assumptions

1. **Does `ParameterEstimationError::Domain` need splitting?** (affects M6) `#[error(transparent)]` names neither the sample nor the operation, `#[from] DomainError` collapses five origin constructors into one variant, and `Domain` is a mechanism name rather than a domain one. The shape is inherited verbatim from `arch` §5.4, so changing it amends the architecture. Recorded as an owner item.
2. **Do the spec and architecture need correcting?** (affects M2) `spec/parameter_prepass_generic.md:785` and `arch/parameter_prepass_generic.md:1128` both describe the runs-model genome as "29% covered by runs" where the research note's realised `F` for that row is 0.2629. The code is fixed; the design documents are not this workflow's to edit.

### 5. Top 3 priorities

1. **M1** — `SampleRates`' accessors return a *different quantity*, not `None`, when the vector is shorter than the ploidy implies. A ploidy-2 set holding one entry answers `homozygous_non_reference_rate()` with the homozygous-**reference** rate: near 1.0 where the truth is near 0.001.
2. **M2** — two documented measurements are wrong: "29% covered by runs" against a realised 0.2629, and a precision claim about `MIN_SITES_TO_FIT` that runs backwards.
3. **M4/M5** — two mutations leave the whole suite green: replacing `Display for Ploidy`'s body with a constant, and swapping `RunsModelStarts::default()`'s two fields.

### 6. Findings

#### Major

**M1: generic/mod.rs — `SampleRates`' accessors read by dosage, and nothing guarantees the dosages are there**
**Categories:** reliability, idiomatic — convergent
**Confidence:** High. Measured:
```
PROBE diploid len=1: het=None homalt=Some(1.0)
PROBE tetraploid len=2: het=None homalt=Some(0.5)
PROBE diploid len=7 (sum 2.8): het=Some(0.2) homalt=Some(0.7)
```
`homozygous_non_reference_rate()` returns the *last* entry, so a diploid holding only entry 0 hands back the homozygous-reference rate under the homozygous-non-reference name — the exact inversion, near 1.0 where the truth is near 0.001. No test built a wrong-length set.

**The commit's recorded assumption was wrong on the half that matters.** It reasoned that entries are already `GenotypeFrequency` so none can leave `[0, 1]`, and that the length-and-sum check belongs with the Milestone E fit. The *sum* half holds — no accessor reads the sum, so an unchecked sum cannot produce a wrong number until that fit exists. The *length* half does not: `homozygous_non_reference_rate()` depends on the length today. The `Option` return added to guard the empty case guards only the empty case, and pushes an unactionable `None` into every call site while failing to stop the confidently wrong read.

**M2: two documented measurements do not match their source**
**Category:** naming (documentation accuracy)
**Confidence:** High.
- **"29% covered by runs"** (`runs.rs`). The genome that produced `F` = 0.2634 against a collapsed 0.0000 is the 3-per-kb row of research §3.4, whose **realised `F` is 0.2629** — 26%. The nearest 29% in the note (0.2886) is a different draw in §3.5. The code copied this from `spec` §6.5 and `arch` §5.3, **which carry the same error**, so the fix here does not close it upstream.
- **`MIN_SITES_TO_FIT`'s precision clause.** "six million read observations pin an error rate to one part in eighty, so a fit is precise long before this" runs backwards: six million observations is two million sites at three reads — **200 times** the 10,000-site floor. The measurement says nothing about precision below the floor; the same arithmetic at 10,000 sites gives ~30,000 observations and about one part in five. The clause is new here; the architecture stops before it.

Eleven other numeric claims in these files were checked against the research note, both specs and the architecture, and **all eleven hold exactly**.

**M3: the implementation report states a wrong reason for one of its own tests**
**Categories:** reliability, naming
**Confidence:** High. The report and the commit message both say a `>=` slipped into `observed_heterozygosity()` "would let a haploid answer, and only this test would notice". Neither half holds: `1 >= 2` is false. Measured — `>= 2` is caught by the **tetraploid** test, and `<= 2` is the mutation that lets a haploid answer. The test is genuinely valuable; the stated reason for it was wrong.

**M4: `impl Display for Ploidy` has no test**
**Category:** reliability
**Confidence:** High. Replacing its body with `write!(f, "MUTANT")` leaves the whole `ng::` lib suite green. The A6 oracle asserts the sample, the site count and the floor — never the ploidy.

**M5: `RunsModelStarts::default()`'s two fields can be swapped with every test passing**
**Categories:** reliability, defaults
**Confidence:** High. Both are `SmallVec<[f64; 3]>`; the tests assert only length 3, ascending, and inside `(0, 1)` — properties both triples share. A swap puts the spread on the wrong axis, which is exactly the silent `F` = 0.0000 the type exists to prevent. The impl plan pins the values by name; nothing executable did.

**M6: `ParameterEstimationError::Domain` breaks three error rules at once**
**Category:** errors
**Confidence:** High. `#[error(transparent)]` names neither sample nor operation, against the standard the module's own test doc states; `#[from] DomainError` collapses five origin constructors into one variant; and `Domain` is a mechanism name. Inherited verbatim from `arch` §5.4 — see open question 1.

#### Minor

**Mi1: `RunsModelStarts::separations` is inverted** — 0.05 means the states are guessed *far apart*, 0.75 close together, on a `pub` unvalidated field whose misconfiguration is the one path to a silent `F` = 0.0000. **Category:** naming.

**Mi2: `RunsModelStarts::default()` was undocumented** — behaviourally significant, no doc on the impl, and its values written nowhere the compiler or a reader could check, where the architecture states them on the fields. **Category:** defaults.

**Mi3: the floors lack the citations and soft/fixed labels `INBREEDING_WINDOW_BP` sets as the standard**, and `MAX_COUPLED_FIT_ITERATIONS` justifies "generous" with a fact about how many samples iterate rather than how many iterations they need. **Category:** defaults.

**Mi4: `mod.rs` hard-wires the SNP/indel floors into a `#[non_exhaustive]` error the STR path is documented to extend** — an STR failure would print "need 10000" at the wrong grain. Carrying `floor` on the variant fixes it and deletes the parent→child import. **Category:** module_structure.

**Mi5: `MIN_WINDOWS_TO_FIT_INBREEDING` is the only runs-model item outside `runs.rs`**, where the architecture declares it inside §5.3 — the section that maps to that file. **Category:** module_structure.

**Mi6: eleven raw `f64`s** across `RunsModelFit`, `StartOutcome` and `ScanResult` where `GenotypeFrequency`, `InbreedingF` and `LogProb` already exist one file away. `LogProb` is the free one — unconstrained, `PartialOrd`, and comparison is the only thing `log_likelihood` is for. Applying it cost two lines. **Category:** smells.

**Mi7: `starts_tried`'s "sorted best first" is prose on a `pub` field** with no accessor, on the field documented as the only thing separating a real `F` = 0 from a failed search. **Category:** smells.

**Mi8: the floor assertions pass vacuously** — `message.contains("10000")` is also true of a floor of 100,000. **Category:** reliability.

**Mi9: `homozygous_non_reference_rate`'s doc conflates its own value with a sum.** The quantity belonging to the (individual, reference) pair is *heterozygosity + this rate*, per spec §5. **Category:** naming.

**Mi10: `inside_het_floor` is named a floor and documented as a fitted rate.** **Category:** naming.

#### Nits

`SmallVec<[GenotypeFrequency; 3]>` contradicts its own doc comment, which justifies the vector by "at `P = 4` there are five entries" — measured, five entries spill; widening to 5 costs 16 bytes per `SampleRates`. `[StartOutcome; 9]` is right at 288 bytes inline. `#[non_exhaustive]` on `ParameterEstimationError` is justified by a comment naming in-crate additions, where the attribute has no same-crate effect. `implied_f` reads as a formula symbol.

### 7. Out of scope observations

- The pre-existing `ng::locus_generation::pileup::parity` failure, unchanged.
- The `29%` slip in `spec/parameter_prepass_generic.md:785` and `arch/parameter_prepass_generic.md:1128`.

### 8. Missing tests to add now

1. `a_set_with_the_wrong_number_of_entries_for_its_ploidy_is_rejected` — the inversion (M1).
2. `a_set_that_does_not_sum_to_one_is_rejected_within_a_rounding_tolerance`.
3. `ploidy_displays_as_the_bare_copy_number` (M4).
4. `the_default_starts_are_the_values_the_design_specifies_on_each_axis` (M5).
5. `best_start_reads_the_head_of_the_best_first_ordering` (Mi7).
6. The floor assertions rewritten against `MIN_SITES_TO_FIT.to_string()` (Mi8).

### 9. What's good

- **Eleven of thirteen numeric claims checked out exactly**, including the harder ones — 0.2634 against 0.0000, nine starts against five, 0.23 and 0.84 at 1,200 windows, 1,550 and 157 of 1,707, 125 and 10³⁵.
- **The four-file split is sound and reproduces the architecture's own module tree line for line** — both directories have real siblings, no `super::super::`, no back-references.
- **No derive here is the `DepthBinEdges` trap**; `Eq` is correctly absent everywhere it must be, and `SampleRates` compares reflexively because `GenotypeFrequency::try_new` rejects `NaN`.
- **Refactor safety came back clean where it counts**: adding a field to `GenericSampleParameters` is a compile error at every site; zero hits for `Default::default`, `derive(Default)` or `match` across the four files; renaming an error-variant field errors through thiserror's format-string resolution.
- **`CoupledFit`/`GenericSampleParameters` field overlap was checked and cleared** as two different things that happen to share a shape, not duplication.

### 10. Commands to re-verify

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --all-targets --all-features`, `cargo doc --no-deps --lib` — all via `./scripts/dev.sh`.

Per-category files kept as an audit trail in `tmp/review_2026-08-06_ng-param-A5A6/`.
