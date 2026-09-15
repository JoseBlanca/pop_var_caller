# Code Review: portable_float_B2
**Date:** 2026-09-14
**Reviewer:** rust-code-review skill (orchestrator, single-pass — see §1)
**Scope:** uncommitted diff of plan step B2 (`doc/devel/implementation_plans/portable_float.md`, Milestone B): std transcendental calls in `src/alignment/` and `src/locus_generation/` replaced by `crate::float`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the working-tree diff `git diff HEAD -- src/` in the worktree
  `pop_var_caller-portable-float`, on top of `342ce9ec` (B1, `src/float.rs`). 12 files, +186 / −152.
- **Reviewed against:** `HEAD` = `342ce9ec` plus uncommitted changes.
- **In-scope files:**
  `src/alignment/emission.rs`, `src/alignment/ssr_anchor_firm.rs`, `src/alignment/ssr_anchor_robust.rs`,
  `src/alignment/ssr_best_path_flat_gap.rs`, `src/alignment/ssr_best_path_unit_slip.rs`,
  `src/alignment/ssr_marginal_sequence.rs`, `src/alignment/ssr_marginal_whole_read.rs`,
  `src/alignment/ssr_noise_robust.rs`, `src/alignment/ssr_robust_indel.rs`,
  `src/alignment/ssr_unit_robust.rs`, `src/alignment/stutter.rs`,
  `src/locus_generation/pileup/minted_error_census.rs`; direct callees `src/float.rs`; the
  recorded-answer oracles `src/alignment/delimit_parity.rs` and
  `src/alignment/testdata/delimit_parity_production.tsv`; the conversion script `tmp/b2_convert.py`.
- **Deliberately out of scope:** `calling`, `genetics.rs`, `types.rs`, `paralog`,
  `parameter_estimation` (steps B3–B5); `examples/` (not converted, by plan); `src/float.rs` itself
  (reviewed as B1).
- **Categories:** the skill's per-category fan-out with one mutation-testing agent per worktree was
  **not** run. The dispatcher forbade running `cargo` and editing source, which removes the thing
  the fan-out's isolation exists for; the diff is a mechanical rewrite whose review questions (a)–(e)
  were set by the dispatcher. They were answered in one pass that covers `reliability`
  (semantic equivalence, tests), `refactor_safety` (completeness), `naming`/doc accuracy, and
  `extras` (stable output, diff matches intent). `unsafe_concurrency`, `module_structure`,
  `tooling`, `defaults`: no triggers in this diff.

## 2. Verdict

**Approve-with-changes.** Every converted expression means what it meant before, apart from the
library that rounds it; no call is left behind. Two things need fixing before commit: the new
constant's doc comment gets the rounding backwards (Mi1), and the new test was inserted inside an
existing test's doc comment, fusing the two (Mi2).

## 3. Execution status

- **Commands run by the reviewer:**
  - `rg` completeness sweeps over `src/alignment` and `src/locus_generation` (method calls,
    line-start chains, `f64::`/`f32::` paths, `libm::`, `mul_add`) — exit 0, results in §6 and
    finding "Completeness" below.
  - `uv run --no-project python tmp/b2_review_probe.py` — exact `ln` of the `f64` nearest 0.01 by
    60-digit `Decimal`, and macOS's `math.log`/`math.exp` bits. Output quoted verbatim:
    ```
    apple math.log(0.01) 0xc0126bb1bbb55515
    dist to ..515 4.3321049331195365730183607901833744555462888E-16 dist to ..516 4.5496792638817157503706925545431880444537112E-16
    pos in ulp from 515: 0.48775165406317581014898509352445702647298434557990888538112
    ```
- **Gates reported by the dispatcher (not re-run):** `cargo fmt`, `clippy -D warnings`, `cargo doc`,
  full container suite 4,794 passed / 0 failed, macOS `alignment::` 326 passed, identity-oracle
  checksums unchanged.
- **Commands not run:** all `cargo` commands, by instruction.
- **Needs verification:** 2 (Q2, Q3).

## 4. Open questions and assumptions

1. **Does the plan record the exception?** Milestone B's checkpoint revision says "Constants are
   not written out as bits" (`portable_float.md` line 122). B2 writes one out. The reason is sound,
   but the plan should say so, and say what happens when B3–B5 hit the next recorded answer that a
   libm rounding flips: write out another constant, or re-record the answer? Affects Mi1, Mi4.
2. **Needs verification — libm's bits for `ln(0.01)`.** The doc says libm gives `…516`. The
   reviewer could not run the `libm` crate; the claim is consistent with the reported flip of
   case 28. Affects Mi1.
3. **Needs verification — did any other `TransitionCosts` field change bits?** `GAP_EXTEND_PROB =
   float::exp(-1.0)`, `ln(1 − 2·2.9e-5)`, `ln(2.9e-5)`, `ln(1 − e⁻¹)` and `ln(e⁻¹)` now come from libm.
   macOS's library rounds all five correctly (0.22 to 0.42 of a step from exact, probe above), so
   they move only if libm rounds one of them wrongly. The parity suite passing says no tie flipped;
   it does not say no bit moved. Affects §6 "Output coverage".

## 5. Top 3 priorities

1. **Mi1** — the doc comment on `LN_GAP_OPEN_PROB_TRACT` says libm is correctly rounded and the
   platform libraries are not. The opposite is true, and the constant's justification should not
   rest on it.
2. **Mi2** — the new test's doc comment is glued onto the tail of
   `the_tract_transitions_leaving_a_match_sum_past_one`'s, so `rustdoc`/readers see one test with
   two unrelated explanations and the other test with none.
3. **Mi3** — the per-quality table lost its only check that its literal bits are the recorded
   ones; a 1-to-4-step edit to one of its 512 values would now pass the suite.

## 6. Findings

### Blocker

None.

### Major

None.

### Minor

#### Mi1: src/alignment/ssr_best_path_flat_gap.rs:206-213 — The constant's doc comment inverts which library rounds `ln 0.01` correctly

**Confidence:** High (for the rounding fact); the libm bits are Needs verification (Q2).
**Categories:** reliability (diff's own quantitative claims, skill step 8a)

[src/alignment/ssr_best_path_flat_gap.rs](../../../../src/alignment/ssr_best_path_flat_gap.rs#L206-L214)
says: *"glibc and Apple's libm both give `0xc0126bb1bbb55515`, and libm … gives the correctly
rounded `…516`"*.

Checked by exact arithmetic on the `f64` nearest 0.01 (§3): the true value of its logarithm lies
4.33e-16 from `…515` and 4.55e-16 from `…516`, i.e. 0.488 of a step from `…515`. **`…515` is the
correctly rounded result**; a library returning `…516` is 0.512 of a step off. Apple's `log`
(checked, Python's `math.log` on this host) returns `…515`.

Claims in the comment, each checked:

| claim | result |
|---|---|
| exact value lies almost halfway between two `f64`s | CHECKED-CORRECT (0.012 of a step from the midpoint) |
| Apple's libm gives `…515` | CHECKED-CORRECT |
| glibc gives `…515` | not checkable on this host |
| libm gives `…516` | Needs verification (Q2); consistent with the case-28 flip |
| `…516` is "correctly rounded" | **WRONG** — `…515` is |
| one unit apart | CHECKED-CORRECT |

This matters beyond pedantry: as written, the comment reads as "we keep the platforms' slightly
wrong value to protect old answers", which invites a later reader to "fix" it back to libm's value.
The truth makes the constant easier to defend: it is the exactly rounded `ln`, and libm is the one
off by half a step.

**Fix** (wording):

```rust
/// `ln(GAP_OPEN_PROB_TRACT)`, **written out as bits**, not computed. The exact value of `ln 0.01`
/// lies within 0.012 of a step of halfway between two adjacent `f64`s. glibc and Apple's library
/// round it correctly, to `0xc0126bb1bbb55515`; libm — which [`crate::float`] uses — returns `…516`,
/// half a step off. That one unit breaks a tie in the repeat-tract delimiter: … A test keeps it
/// within one unit of `float::ln(GAP_OPEN_PROB_TRACT)`.
```

The commit message, if it repeats the claim, needs the same correction.

#### Mi2: src/alignment/ssr_best_path_flat_gap.rs:973-997 — New test inserted inside another test's doc comment

**Confidence:** High
**Categories:** naming / documentation

The diff adds the new doc lines directly after the last line of an existing `///` block, with no
blank line:

```rust
    /// The match→match cost is built from the **flank** gap-open and never recomputed under
    /// the tract regime, so inside the tract the three transitions leaving a match sum to
    /// about 1.02. Reproduced from production deliberately; pinned so that fixing it is a
    /// visible choice rather than an accident (spec §4.2).
    /// The written-out tract gap-open cost is `ln(GAP_OPEN_PROB_TRACT)` to within one unit in the
    /// …
    #[test]
    fn the_written_out_tract_gap_open_cost_is_the_log_of_its_probability() {
    …
    #[test]
    fn the_tract_transitions_leaving_a_match_sum_past_one() {
```

[src/alignment/ssr_best_path_flat_gap.rs](../../../../src/alignment/ssr_best_path_flat_gap.rs#L973-L997).
The 1.02-sum explanation now documents the wrong test, and the test it belongs to has no doc.

**Fix:** move the new test (with only its own three doc lines) above the `/// The match→match cost
is built…` block, so that block sits directly on `the_tract_transitions_leaving_a_match_sum_past_one`
again.

#### Mi3: src/alignment/emission.rs:612-630 — The table's recorded bits are no longer checked by anything

**Confidence:** High
**Categories:** reliability (test challenge), refactor_safety

Before B2, `per_quality_table_matches_the_dindel_model` asserted bit equality on Linux against
glibc's `powf`/`ln`, which verified the doc's claim that the 512 literal values are "the bits glibc
produced". After B2 it uses a tolerance of `4·ε·max(|v|,1)` everywhere. The tolerance change is
correct — libm and glibc round some entries differently — and the doc now says exact bits "are fixed
by the literal itself". But the literal is only fixed until someone edits it: a typo that moves one
entry by up to four steps (a transposed low hex digit, say `…d27` for `…d29`) passes every test,
and the module doc says such a move "can move a measured repeat by a byte". `delimit_parity` would
catch it only if it happens to flip a tie among its 11,992 cases.

Which wrong tables pass the new assertion: any table within 4 steps per entry of libm's — including
libm's own table, the one the plan decided *not* to adopt.

**Fix:** add a pin on the bits, independent of any maths library — see §8,
`per_quality_table_bits_are_the_recorded_glibc_bits`.

#### Mi4: doc/devel/implementation_plans/portable_float.md:122 — The step contradicts the plan's checkpoint revision without amending it

**Confidence:** High
**Categories:** extras (diff matches stated intent)

The Checkpoint A revision says constants are not written out as bits once calls go through libm.
B2 introduces `LN_GAP_OPEN_PROB_TRACT` as bits, for a different reason than the one the revision
rules out (a recorded answer from a deleted caller, not compiler folding). The plan should record
the exception and the rule for the next one (Q1) — B3–B5 convert 325 more calls, several of which
feed recorded-answer tests (`calling/quality_parity.rs`, `calling/loop_parity.rs`,
`paralog/production_parity.rs`).

**Fix:** a bullet under "Revised at the checkpoint" naming `LN_GAP_OPEN_PROB_TRACT`, why, and the
policy for further flips.

### Nits

- [src/alignment/ssr_marginal_sequence.rs](../../../../src/alignment/ssr_marginal_sequence.rs#L269-L277)
  and [#L795](../../../../src/alignment/ssr_marginal_sequence.rs#L795): prose still says "that single
  `.ln()` is the whole of the conversion", "`0f64.ln()` is `-∞`", "`.ln()` already maps `0` to `-∞`".
  The code now calls `float::ln`; the statements remain true of libm's `log(0)`, but once B6 bans
  `f64::ln` the docs name a banned method. Suggest "`float::ln`" / "`ln 0 = −∞`".
- [src/alignment/emission.rs:674](../../../../src/alignment/emission.rs#L674): `assert_eq!(match_ln,
  float::ln(2.225_073_858_507_201_4e-308))` compares a glibc-recorded table entry to libm's `ln` for
  exact equality — the "two libraries in one assertion" pattern the plan wants gone. It passes
  (the two agree there), but a libm upgrade could fail it with the table still correct. Compare to
  the bits `0xc086232bdd7abcd2` (the Q0 match entry) instead, which is what the test means.
- [src/alignment/emission.rs:867](../../../../src/alignment/emission.rs#L867)
  `uniform_base_ln_is_ln_of_a_quarter`: same pattern, lower risk (`ln 0.25 = −2 ln 2`, and every
  library tested rounds it alike). Leave, or pin bits.
- The `cfg_attr(not(test), expect(dead_code, …))` on `GAP_OPEN_PROB_TRACT` is correct, but its
  `reason` string is the only place that says the constant is now documentation-plus-test-anchor;
  consider one sentence in its doc ("Not read by shipped code; the shipped cost is
  [`LN_GAP_OPEN_PROB_TRACT`]").

### Answers to the dispatcher's questions

**(a) Semantic equivalence — every hunk.** No hunk changes meaning beyond std → libm/`float::powi`.
Checked specifically:

- *Receiver precedence:* every `(expr).ln()` became `float::ln(expr)` with the parenthesised
  receiver as the argument, including the three-line `SlipCosts` chains in all six aligners
  (`(a * b)\n.ln()\n - s [- margin]` → `float::ln(a * b) - s [- margin]`). No unary minus applied to
  a method call anywhere in the diff. `(lo - hi).exp().ln_1p()` → `float::ln_1p(float::exp(lo - hi))`,
  inner/outer order preserved. `(x.max(FLOOR)).ln()` → `float::ln(x.max(FLOOR))`, `max` still
  inside.
- `float::ln(1e-4 / 3.0).abs()` (`ssr_robust_indel.rs`): `.abs()` still applies to the logarithm,
  as `(1e-4f64 / 3.0).ln().abs()` did.
- *Literal suffixes:* `0.5f64`-style receivers became unsuffixed literals passed to `fn(f64, …)`
  parameters, so they infer `f64`; `1.0 / 3.0`, `0.02 / 3.0`, `0.000_1`, `2.225…e-308` likewise.
- *Integer exponents:* every `powi` exponent was already `i32` (`as i32`, `i32::try_from(..).unwrap_or(i32::MAX)`);
  `float::powi` takes `i32`. The `usize as i32` wrap for `read.len()` is unchanged and behaves
  identically under B1's `powi` (B1's edge test covers negative and `i32::MIN`/`MAX` exponents).
- *Closures / generics:* the `reached` closure (`stutter.rs:487`) and the `LazyLock` initialiser
  keep their signatures; no generic contexts or iterator adaptors (`.map(f64::ln)`) were involved.
- *One intentional non-mechanical change* (`ln_gap_open_tract`) — see Mi1; the test tying it to
  `float::ln` accepts at most one step and also asserts `TransitionCosts::new()` uses it.

**(b) Completeness.** No transcendental call remains in `src/alignment/` or
`src/locus_generation/`:
- `rg '\.(ln|exp|powf|powi|log10|log2|log|ln_1p|exp_m1|sin|cos|tan|…|mul_add)\s*\('` → only doc-comment
  prose in `ssr_marginal_sequence.rs` (lines 269, 270, 274, 277, 760, 795).
- `rg -U '^\s*\.(ln|exp|…)\b'` (line-start chains) → none.
- `f64::ln`-style paths, `f32::`, `libm::`, `mul_add` → none (the two `f64::powf`/`f64::ln` hits at
  `emission.rs:162-163` are historical prose, accurate).
- `locus_generation`'s per-read error comes from `open_record.rs:3228` `phred_to_ln_perr`, which
  multiplies by `LN_10` — arithmetic only, already portable.
- Indirect dependence remains through code B3 converts: `crate::types` has two std `ln` calls, both
  in tests (`types.rs:2325`, `2341`); `calling/likelihood/ssr_emission.rs:617` builds a `FlatEmission`
  with a rate that comes from the fit (B5).

**(c) Docs touched by the conversion.** The script edited one doc comment deliberately
(`emission.rs` "pins the literal against `float::ln(0.25)`" — correct). No string literal or code
example was altered. Stale prose: Nits above. Wrong prose: Mi1. Misplaced prose: Mi2.

**(d) Tautologies versus deliberate std comparisons.**
- *Deliberate comparison lost:* `per_quality_table_matches_the_dindel_model` on Linux compared the
  table to glibc exactly; now a tolerance (Mi3).
- *Still independent (not tautological):* the six `SlipCosts` reconstruction tests compare
  `float::ln` of the costs against `float::ln(model.probability(..))`, which goes through
  `float::powi` and the stutter arithmetic — a different path, 1e-12 tolerance. The `stutter.rs`
  decay tests build the expectation with `float::powi` against `Regime::probability`'s own
  `float::powi` at 1e-15: same function on both sides, but the shares and step counts are spelled
  independently, so a wrong exponent or share still fails. `ssr_marginal_sequence.rs:782` compares
  against a test-local port of production's `align_subst` (both sides now libm, both were std) —
  character unchanged.
- *Tautological by design, unchanged:* `marginal_probability_is_the_single_logarithm_of_the_linear_probability`
  (`to_bits` of `float::ln(linear)`) — it pins the one-log boundary and should share the function.
- *Mixed-library exact equalities:* `emission.rs:674`, `emission.rs:867` (Nits).

**(e) Output coverage.** Shipped calls converted, and what observes them:

| shipped call | reaches output via | recorded-answer test |
|---|---|---|
| `TransitionCosts::default` (5 costs, `GAP_EXTEND_PROB`) | every tract aligner | `delimit_parity` (algorithm 3 only, 11,992 cases, 12,000-line TSV) |
| `SlipCosts::from_model` in `ssr_unit_robust.rs` | **the shipped SSR generator** (`locus_generation/ssr.rs:2297`, algorithm 4u) | **none** — only the identity oracle |
| `SlipCosts` in anchor_firm / anchor_robust / noise_robust / robust_indel / unit_slip | research aligners, tests | none (measured, not verified — `delimit_parity.rs` module doc) |
| `FlatEmission::try_new` | `calling/likelihood/ssr_emission.rs:617` | calling parity tests (`quality_parity`, `loop_parity`) |
| `StutterModel` `reached` / `Regime::probability` `powi` | calling's SSR likelihoods | calling parity tests; `powi` is bit-identical to std at run time (B1 test) |
| `ssr_marginal_whole_read::ln_sum2`, `SsrSequenceMarginal` | research marginals | none |
| `minted_error_census` `exp` | only when `PVC_MINTED_ERROR_CENSUS=1`; measurement output | none needed |

So the gap is the shipped generator's `SlipCosts` (`ln` of fitted stutter shares): a rounding flip
there is seen by no recorded-answer test, only by the identity oracle's five checksums, which the
dispatcher reports unchanged. That covers the oracle's inputs, not the full input range (Q3; the
parity oracle's VCF at B7 is the stronger check).

## 7. Out of scope observations

- `src/types.rs:2325` and `:2341` — test-only `p.ln()` / `1e-300f64.ln()`; B3's scope, listed so the
  ban (B6) does not surprise.
- `delimit_parity.rs` covers algorithm 3 only, while the shipped generator uses algorithm 4u. Not
  new, but every B-step conversion that touches `SlipCosts` inherits it; the B7 parity VCF is the
  place to close it.

## 8. Missing tests to add now

**`PER_QUALITY_LN` (`src/alignment/emission.rs`)**

- `per_quality_table_bits_are_the_recorded_glibc_bits` — input class: the whole table; catches a
  one-to-four-step edit to any entry, which the tolerance test (Mi3) cannot. Spec:

  ```rust
  #[test]
  fn per_quality_table_bits_are_the_recorded_glibc_bits() {
      // FNV-1a over the 512 bit patterns, in table order. Recorded from the table as committed at
      // 75722336 (glibc, aarch64 Linux). A change here is an edit to the table, not rounding.
      let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
      for q in 0..=u8::MAX {
          let s = PerQualityEmission::new().scores_for(BaseQual(q));
          for bits in [s.match_ln.to_bits(), s.mismatch_ln.to_bits()] {
              for byte in bits.to_le_bytes() {
                  hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
              }
          }
      }
      assert_eq!(hash, RECORDED_TABLE_HASH);
  }
  ```

  (Record `RECORDED_TABLE_HASH` by running once on the committed table.)

**`TransitionCosts` (`src/alignment/ssr_best_path_flat_gap.rs`)**

- `transition_costs_are_pinned_to_the_bit` — input class: the five shipped costs; answers Q3 once
  and for all by recording each field's bits, so a future libm change that moves one is reported by
  name rather than only if it flips a parity tie.

## 9. What's good

- `tmp/b2_convert.py` makes every replacement an exact string with an expected count, so a pattern
  that matched zero or two times fails loudly instead of converting partially.
- `the_written_out_tract_gap_open_cost_is_the_log_of_its_probability`
  ([ssr_best_path_flat_gap.rs:981](../../../../src/alignment/ssr_best_path_flat_gap.rs#L981)) counts
  steps between bit patterns rather than using a relative tolerance, and also asserts the constant is
  the one `TransitionCosts` actually uses — it fails both on a stale constant and on a constant
  nobody reads.
- The `cfg_attr(not(test), expect(dead_code, reason = …))` turns "the probability constant is now
  only documentation" into something the compiler re-checks: if shipped code starts reading it
  again, the build says so.

## 10. Commands to re-verify

- Reviewer ran: `uv run --no-project python tmp/b2_review_probe.py` (rounding of `ln 0.01`);
  the `rg` sweeps listed under (b).
- New, for the author:
  - `./scripts/dev.sh cargo test --lib alignment::ssr_best_path_flat_gap` after moving the test (Mi2).
  - `./scripts/dev.sh cargo test --lib alignment::emission` after adding the table pin (Mi3).
  - `./scripts/dev.sh cargo doc --no-deps` and inspect the two flat-gap test docs.

**PROJECT_STATUS.md** was not updated: the plan schedules its entry for step C3, the feature has no
block yet, and the dispatcher asked for the report only.

## Fixes applied (2026-09-14)

| finding | what was done |
|---|---|
| Mi1 rounding described backwards | Confirmed: the implementer's own check took the logarithm of exact 0.01, not of the `f64` nearest it. The constant's doc comment and the plan's §8 now say glibc and Apple return the correctly rounded `…515`, and libm returns `…516`, half a step off |
| Mi2 new test inside another test's doc comment | the match→match doc block moved back above `the_tract_transitions_leaving_a_match_sum_past_one` |
| Mi3 table's exact bits no longer checked | `per_quality_table_keeps_its_written_bits`: an FNV-1a digest over all 512 bit patterns, `0x0454_3206_1128_dbc0`, recorded from the table unchanged since `75722336` |
| Mi4 plan does not cover the exception | plan §8 records it and states the rule: a constant is written out as its old bits only where converting it moves a recorded answer |
| nits | the `.ln()` prose mentions stay: the ban (B6) covers calls, not prose. The two exact comparisons of a recorded value against libm's `ln` pass and are left as they are |
| (e) algorithm 4u's slip costs are observed only by the identity oracle | noted; the oracle's five checksums are unchanged after this step |

Validation after the fixes: fmt, clippy `-D warnings` and doc clean in the container; container suite 4,795 passed, 0 failed, 4 ignored; `alignment::` 327 passed on macOS and in the container. The identity oracle ran on the release build before the fixes, which touched only a comment, a test and a moved doc block, and its five checksums were unchanged.
