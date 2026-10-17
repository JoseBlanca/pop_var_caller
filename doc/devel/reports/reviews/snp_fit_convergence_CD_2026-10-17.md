# Correctness review: SNP/indel fit convergence, plan steps C1–D2 (2026-10-17)

*Written by the reviewing agent; filed by the lead (the harness refused the agent's file write).*

**Scope.** Commits `6dc9e123..dda9e929` (C1 the regularised incomplete Beta moved to `genetics.rs`; C2
`cohort_covariance`; C3 `fit/determined.rs` and its reports on `JointFit`; D1 the run-log line; D2
parameters file format 4 with two informational tables), reviewed as committed (an export of `dda9e929`).

## Verdict

**No Blocker and no Major.** The trade-cost maths, the covariance, the delta-method errors, the
incomplete Beta, and format 4's reading, writing and validation are correct; no VCF can move. The
59-sample verdicts were confirmed with exact (mpmath) derivatives: the density's trade cost 1.3467282879
(code 1.346728287), the duplicated pair's 0.1157063099 (code 0.1157063100), the error of S 2.03984889e-4
(code 2.0398488887e-4). Four Minor findings, two Nits.

**Checks (export):** `cargo fmt --check` and `clippy --all-targets --all-features -D warnings` exit 0;
the lib in release: 5,195 passed, the 3 debug-only should-panic tests fail; the incomplete-Beta,
`determined`, verdict, digest, `parameters_file` and `calling_run` tests pass. The script and `sed -i`
edits: `git diff --check` clean, no CR line endings, no mojibake; every removed line is the C1 move, the
argument reorder, or an import or signature change; the golden TOMLs differ only in the version line and
the two tables.

## The questions asked

- **The 4-sample trade cost of 4.66e16:** the verdict is right and rounding cannot flip it. At a = 15.4,
  b = 50, P(f < 0.01) is 2.5e-17, so 1 − I is exactly 1.0 in floating point; S's gradient equals σ's to
  1e-10 and V_trade = Var(σ) − Cov²/Var(S) is pure cancellation (NumPy: exactly 0; Rust: one ulp). For
  rounding to give a cost below 2, σ's squared relative error would need to exceed about 5e14.
- **Delta-method errors:** correct; they reproduce the spec's tables, except near a share's lower bound
  (Minor 1).
- **Incomplete Beta near the bounds:** 80 points at x = 0.01 against mpmath at 60 digits, shapes 0.02 to
  50: I accurate to 8.7e-14 relative, 1 − I to 2.5e-12; the symmetry branch is never taken at x = 0.01;
  a shape's central difference is accurate to about 1e-7. The move is arithmetic-identical; the argument
  reorder is correct.
- **Format 4:** version-3 files read with both tables absent; version 5 refused; a version-3 build
  refuses a version-4 file; migrated version-1/2 files are written back as 4; a calling run writes both
  tables back (`src/cli/calling_run.rs:640-654`); a fit cannot write a file it would refuse.
- **Calling reads none of the new tables:** it reads only `reference_concentration`,
  `alternative_concentration_total` and `rung` (`to_run_parameters.rs:246-248`); `bindings.rs:708-710`
  ignores both tables; the digests' VCF checksums are unchanged.

## Findings

### Minor 1 — `determined.rs:264-273` (`Group::gradient`): central difference wrong near a share's lower bound

The step is `1e-5 × |value|`; `p_fixed_alt`'s lower bound is 1e-12 (`fit.rs:3332`), giving a step below
one ulp of σ. Measured: ∂S/∂p_fixed_alt −0.0867 against an exact −0.0533 (at 0.84, 1e-12, 0.02, 4.28);
−0.694 against −0.893 (at 0.9, 1e-12, 0.6, 2.2); ∂E[f]/∂p_fixed_alt 0.694 against 0.786. At the
59-sample fit (p_fixed_alt 2.96e-5) the error is about 1e-8. These gradients feed the written errors of
S, Hexp and E[f] and V_trade; a cohort sequenced against its own reference plausibly sits on that bound.
**Fix:** S, σ, E[f] and Hexp are linear in p_invariant, p_fixed_alt and δ — use exact gradients there
(∂S/∂p_invariant = ∂S/∂p_fixed_alt = −(1 − I); ∂S/∂δ = 1 − I; ∂E[f]/∂p_fixed_alt = 1 − a/(a+b),
∂E[f]/∂p_invariant = −a/(a+b); likewise Hexp) and central differences for the shapes only.

### Minor 2 — `census_fit.rs:559`: no test catches a wrong condition for writing the density's table

Mutating `rung == SeedRung::FittedCurve` to `true` fails nothing new; the writers
(`census_fit.rs:571-618`) are covered only by the macOS md5 digests (a swap of written Hexp and E[f] was
caught only there). **Fix:** a unit test on a `JointFit` with known reports: on `FittedCurve` every
written field equals its report field, `separately_determined` included; on every other rung
`fitted_density` is `None` while `duplicated_positions` is still written.

### Minor 3 — `determined.rs:307-319, 338-346`, `standard_errors.rs` (`a_cost`): "no trade" gets two wordings

Cancellation landing on exactly 0 gives +∞ ("cannot change without changing it"); one ulp above gives
4.66e16 ("over a million"): the same state. `a_cost`'s doc says such a cost "is rounding", but drawn
cohorts measured genuine costs of 1.4e6 and 4.3e6. **Fix:** treat `trade_variance <= 1e-12 * variance` as
no trade (+∞); reword the doc: a cost that large is rounding or a trade far past the line.

### Minor 4 — `standard_errors.rs` (`in_three_figures`): a cost just under 2 prints "2.00" beside a mark

Decimals are chosen before rounding: 1.996 prints "2.00" next to a verdict saying it is below 2 (also
0.9996 → "1.000", 999.6 → "1000"). **Fix:** choose decimals after rounding; never let a marked pair's
printed cost reach 2.

### Nits

- `validate.rs:128-180` accepts the new tables in a file stamped version 3, `fitted_density` on a rung
  other than `fitted_curve`, and a pair whose halves disagree on `separately_determined`; within the spec
  ("nothing else is checked"), but refusing format-4 tables in a version-3 file would match the migration's
  refusal of format-3 keys in version-1/2 files (`from_toml.rs:295`).
- The test `a_version_three_file_reads_without_the_tables_and_a_version_five_file_is_refused` injects its
  unknown key with `.replace`, editing both tables at once.

### Spec-level observation

At one sample `fixed_non_reference_share` keeps its error (2.08e-5), computed with the unidentified
shapes held fixed — the reason the spec gives for withholding p_invariant's error, though p_fixed_alt is
probably near-independent of the shapes. The code does what the spec says.

## The implementer's decisions

| step | decision | judgement |
|---|---|---|
| C1 | move the existing function, not a second copy | correct; arithmetic identical |
| C2 | `[[Option<f64>; 8]; 8]`, a dropped parameter `None` across its row and column | correct; tested against the dense inverse |
| C3 | `trade_cost` `None` when Var(S) is not positive finite | consistent; cannot occur inside the bounds (1 − I ≥ 0.011) |
| C3 | derived numbers keep their errors when the pair is marked | matches the spec's example tables |
| D1 | "over a million" | safe for the verdict; wording needs Minor 3 |
| D2 | `parameters_file_of` takes the `JointFit` | fine |
| D2 | the density's table written on `fitted_curve` only | matches §2.1.3; untested (Minor 2) |

Housekeeping: an accidental build of the live tree (read-only) warned that `HarshCase` and `a_harsh_case`
in the other agent's in-progress `standard_errors.rs` are unused. Scratch: `tmp/review_CD/` (`py/vtrade.py`,
`py/pfa_bound.py`, `py/t59.py`, `py/ib_grid.py`, logs, `src_copy/`).
