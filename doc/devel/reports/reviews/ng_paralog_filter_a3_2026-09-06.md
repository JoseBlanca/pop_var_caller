# Code Review: ng_paralog_filter_a3
**Date:** 2026-09-06
**Reviewer:** rust-code-review skill (orchestrator) over two category sub-agents, each in its own worktree
**Scope:** commit `2f6a6c32` on `ng-paralog-filter` — the prior, the false-discovery curve and the cut copied into `src/ng/paralog/`, and the guard extended for a copy that takes a span of a larger file
**Status:** Request-changes (one Blocker), all resolved before the commit was amended

---

### 1. Scope

- **Reviewed against:** `2f6a6c32`, branch `ng-paralog-filter`, parent `4b3a8287`.
- **In-scope files:** [copy_fidelity.rs](../../../../src/ng/paralog/copy_fidelity.rs), [production_parity.rs](../../../../src/ng/paralog/production_parity.rs), [mod.rs](../../../../src/ng/paralog/mod.rs) (ng's own); [prior.rs](../../../../src/ng/paralog/prior.rs) and [calibration.rs](../../../../src/ng/paralog/calibration.rs) (copies — their `use` lines, placement, visibility and ng's headers only); [src/ng/mod.rs](../../../../src/ng/mod.rs).
- **Out of scope:** the copied logic in all five copies. Frozen; observations routed to §7.
- **Categories dispatched:** `reliability` (the guards are the only thing asserting the port), `module_structure` + `errors` + `naming` (the span copy, the new sanctioned kind of line change, and the prose that states the guard's contract).
- **Skipped:** `defaults` (every `DEFAULT_*` in scope is production's), `unsafe_concurrency` (none), `tooling` (no manifest change), `smells`/`idiomatic` (folded into the second agent's naming pass on this small a diff), `extras` (no parser, no untrusted input).

### 2. Verdict

**Request-changes.** The copies are production's and the numbers agree — a reviewer re-derived the generator in Python from scratch and matched all six pinned shape counts and all four pinned locus counts exactly. What the review found is one **Blocker** in the guard's own sanction check, four Majors where the differential or ng's own tests used the implementation as their oracle, and a set of doc claims that A3's own changes had made false. All are fixed.

### 3. Execution status

| command | result |
|---|---|
| `cargo test --all-features --lib "ng::paralog::"` | `ok. 77 passed; 0 failed` |
| `cargo test --all-features --lib --bins --tests` | 6,352 lib tests passed; one integration test failed, pre-existing |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 3 errors in `cohort_merge/`, pre-existing |
| `cargo fmt --check` | dirty on 9 files, none this commit's, pre-existing |

Findings labelled "Needs verification": **0**.

### 4. Open questions and assumptions

1. **`lr_threshold` is written into the VCF header for provenance, and it is not what `flags` decides on.** See **M4**. The relation is now pinned as "agree to within one bin"; whether the *header* should carry the bin's edge rather than its centre is a question for step C4 and, behind it, for production's own filter.

### 5. Top 3 priorities

1. **B1** — two of the guard's three sanction checks tested only for substring presence, so a declared substitution could carry an arbitrary edit into a frozen copy with the suite green.
2. **M4** — the equivalence `calibration.rs` documents between the flag and the recorded cut is false by half a bin, and the recorded cut is what goes in the header.
3. **M2** — the differential stopped one step short of the verdict: `flags` and `posterior` were never compared with production's.

### 6. Findings

#### Blocker

**B1: [copy_fidelity.rs:331](../../../../src/ng/paralog/copy_fidelity.rs#L331) — two of the three sanction checks were containment tests, and the check had no rejecting-direction test at all**
**Categories:** reliability (filed), errors (convergent) · **Confidence:** High
The module header claimed the rule was exact and asserted. Only the visibility check was written that way; the other two tested that the two lines *contain* certain substrings and left everything else on the line unconstrained. Demonstrated twice on the real files: declaring `locus_score.rs`'s import repoint with an ng-side line of `use crate::ng::paralog::{GridSpec, SfsPriorSpec}; const SMUGGLED_THRESHOLD: f64 = 0.05;`, and putting that line into ng's copy, left `ng::paralog::copy_fidelity` at `6 passed; 0 failed`; likewise a doc-link repoint rewritten to `/// Tuning knobs. The fallback prior is 0.05, not 0.03.` **The declaration is the one thing the content comparison cannot cross-check** — the comparison rewrites production's file *using* the declaration — so it has to be checked directly. And the `assert!` had no test in its failing direction for any kind. Separately, when it did fire it blamed ng's copy, where the guard uses `THE GUARD COULD NOT RUN` everywhere else.

#### Major

**M1: [production_parity.rs:702](../../../../src/ng/paralog/production_parity.rs#L702) — the ratio streams never left the histogram's range, never ran dry and never failed to converge**
**Category:** reliability · **Confidence:** High
Four `prior.rs` mutations survived the differential: both of `bin_index`'s clamps, the empty-histogram seed, and the `converged` flag forced true. Each is a real behaviour change (each fails one of production's transcribed tests). And unlike the scorer's stream, nothing pinned what the ratio streams contained: narrowing `-20.0 * next_unit() - 1.0` to `-2.0 * …` would shrink what the calibration half compares with no test saying so.

**M2: [production_parity.rs:738](../../../../src/ng/paralog/production_parity.rs#L738) — the differential stopped one step short of the verdict**
**Category:** reliability · **Confidence:** High
π, `converged`, the curve at fifteen probes and the cut at five targets were compared; `flags` and `posterior` — *is this record dropped* and *what number goes in `PARALOG_POST`* — were compared with nothing. Not a reachability limit: `crate::var_calling::paralog_filter::calibrate` is `pub(crate)` and nameable from ng's test module.

**M3: [mod.rs:115](../../../../src/ng/paralog/mod.rs#L115) — `flags`'s oracle was `flags`'s own body, at one target**
**Category:** reliability · **Confidence:** High
Two wrong implementations passed: `<` in place of `<=` (the six probes take only three distinct q-values, none equal to `0.01`), and a literal `0.01` in place of `self.target_fdr` — an implementation that ignores the operator's knob entirely.

**M4: [mod.rs:110](../../../../src/ng/paralog/mod.rs#L110) — `lr_threshold` was never read, and the equivalence its doc states is false**
**Category:** reliability · **Confidence:** High
`calibration.rs` documents `flags` as *"equivalent to `lr >= lr_threshold` by the curve's monotonicity"*. Written as a test it fails on the module's own fixture: the recorded cut is `-7.85` and `flags(-7.90)` is `true`, because `lr_threshold_for_fdr` returns the crossing bin's **centre** while `flags` decides on the bin. Every ratio in the lower half of that bin — `0.05` wide at the shipped resolution — is flagged while sitting below the recorded cut. **`lr_threshold` is what the run writes into the VCF header for provenance**, and both ng fixtures set it to `None`, so nothing read it at all.

**M5: [copy_fidelity.rs:220](../../../../src/ng/paralog/copy_fidelity.rs#L220) — the recipe `src/ng/mod.rs` publishes for finding reaches into production cannot see this one**
**Category:** module_structure · **Confidence:** High
A3 adds `include_str!("../../var_calling/paralog_filter/calibrate.rs")`, the first reference from `src/ng/paralog/` to a pipeline-stage module. `src/ng/mod.rs` says the way to find such reaches is `grep -rn 'use crate::' src/ng | grep -v 'use crate::ng'` "rather than a list here that goes stale" — and an `include_str!` is not a `use`, so the stated method has a blind spot at exactly the dependency that points at a stage module.

**M6: [copy_fidelity.rs:37-59](../../../../src/ng/paralog/copy_fidelity.rs#L37) — the header still stated the pre-A3 rules, and asserted two guarantees the code did not give**
**Category:** naming · **Confidence:** High
The heading read "paths into `src/paralog/`" and the text said "a repoint may only turn a `crate::paralog` path into a `crate::ng::paralog` one — asserted", when five of A3's seven declared lines are the new visibility kind. And "line for line, and then byte for byte" had gained an exception — the span-end trim — documented 270 lines below, inside the extractor's body.

#### Minor

**Mi1: the span-end normalisation also swallowed a CRLF and trailing spaces on the last *content* line**, so "byte for byte" had a one-line hole and nothing tested it.
**Mi2: `ends_before` had no "occurs exactly once" check**, so a marker that came to match earlier would truncate production's side and the failure would blame ng's copy.
**Mi3: three doc comments A3 wrote still counted two of what are now three, or two of what are now seven** — `WhyRepointed`'s own doc, `Repoint`'s, and `calibration.rs`'s header.
**Mi4: `src/ng/mod.rs`'s crate-level inventory of `paralog` was not updated**, and under-reported the module by half.
**Mi5: `src/ng/mod.rs` states one visibility rule (widen a *production* item so a parity test can see it) and A3 introduced its mirror image** (widen ng's *own copy*), written down only inside a `#[cfg(test)]` module. The two share a verb and differ in subject — the pair a reader conflates.
**Mi6: leaving `cohort_inbreeding` out of the span was recorded in the commit, the report and `calibration.rs`'s header, but not in `copy_fidelity.rs`** — the module that carries the guarded-copy and release tables and exists so that a dropped entry cannot look like a deliberate one. And nothing linked the end marker to the item: renaming it while keeping its doc sentence left the guard green.
**Mi7: `CalibrationConfig`, its `Default` and `DEFAULT_FALLBACK_PARALOG_PRIOR` are public re-exports with no test in ng**, and the coupling their doc says only a comment holds — that the fallback matches `DEFAULT_EM_START` — is unenforced.

#### Nits

- `Repoint` no longer names what it holds, now that half its instances point nothing anywhere.
- `posterior` returns `Some(NaN)` for a `NaN` prior, since `NaN <= 0.0` and `NaN >= 1.0` are both false. Production's, frozen.

### 7. Out of scope observations

The five equivalent mutants in production's scorer from A2's review still stand. `posterior`'s `NaN`-prior hole above is production's. `cargo fmt` in the container still rewrites nine files this branch never touched.

### 8. Missing tests to add now

All applied — see the fix-application report. In summary: the sanction's rejecting direction for each of the three kinds; a span whose end marker production no longer has, and one that occurs twice; the verdict compared with production's across four stream shapes and five targets; an empty histogram and an EM that cannot converge; a **narrow histogram** where the end bins are not saturated; a target sweep and the `<=` boundary for `flags`; and the true relation between the recorded cut and the flag.

### 9. What's good

- **The span's end boundary is the only place either side is normalised**, and the commit says so, with the measurement that forced it (`cargo fmt` strips the trailing blank).
- **The visibility check was exact from the start**, and is the model the other two were brought up to.
- **`ends_before`'s failure is reported as the guard breaking, not as the copy drifting** — the distinction the whole module turns on.
- **The four ratio-stream shapes are named for what they are** (nothing duplicated, a tenth, everything, unscorable values) rather than drawn uniformly, because the EM's answer turns on the distribution's shape.

### 10. Commands to re-verify

```
scripts/dev.sh cargo test --all-features --lib "ng::paralog::"
scripts/dev.sh cargo test --all-features --lib --bins --tests
scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
```

Per-category audit trail: `tmp/review_2026-09-06_ng_paralog_a3/{reliability,structure_errors_naming}.md`.
