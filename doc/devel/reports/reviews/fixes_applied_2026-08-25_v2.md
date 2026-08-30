# Fix Application Report: ng_calling_loop_a2_2026-08-25.md

**Date:** 2026-08-25
**Source review:** [ng_calling_loop_a2_2026-08-25.md](ng_calling_loop_a2_2026-08-25.md)
**Source state reviewed against:** `a2d28e62` on branch `ng-calling-loop`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 13
- Nits: 7 (grouped)

### Outcome totals
- Applied: 20
- Applied with adaptation: 2
- Already fixed: 0
- Deferred: 1
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 3
- Awaiting user answer: 0

**Three findings are `Superseded` rather than applied, and that is the most useful thing in
this report.** Retyping three counts as `NonZeroU32` and the discovery bar as the merge's
`MinAltReads` **deleted the code** that three findings were about: the pass-cap check, the
zero-reads check and the share range check no longer exist, because the values they refused
are no longer expressible. Two error variants went with them.

### Validation summary

In the container, from this worktree's own `scripts/dev.sh`, at main's 1.98 compiler pin.

- `cargo fmt --check` → **0**, no output
- `cargo clippy --all-targets --all-features -- -D warnings` → **0**, no warnings
- `cargo test --lib` → **0**, `4528 passed; 0 failed; 14 ignored`
- `cargo test --release --lib ng::calling::inference` → **0**, `11 passed; 0 failed`
- `cargo doc --no-deps` → not run; no public item outside this module changed.
- `cargo audit` → not run; no dependency changed.
- Performance check → **not applicable**: nothing changed is reachable from any harness in `benches/`.

A1 left the library at 4,517 and A2 as reviewed was 4,523; it is now **4,528**, so this step
contributes eleven tests, five of them added by these fixes.

### The fixes were verified against the mutations they were written for

Seven mutations survived the review's battery. Each was re-run here against the fixed tree,
and the restore was confirmed by content diff after every one:

| mutation | before | after |
|---|---|---|
| round-threshold constant `1e-3` → `0.5` | survived | **killed** |
| discovery share `0.15` → `0.9` | survived | **killed** |
| discovery round cap `4` → `1` | survived | **killed** |
| round-threshold guard drops `is_finite` | survived | **killed** |
| pull-back guard `is_finite()` → `!is_nan()` | survived | **killed** |
| the two not-built checks swapped | survived | **killed** |
| the convergence ceiling removed | *(new check)* | **killed** |

Three further mutations from the battery — deleting the pass-cap guard, deleting the
zero-reads guard, and dropping `is_finite` from the share check — **cannot be re-run**,
because those three guards no longer exist. See the `Superseded` rows.

### Unresolved high-priority findings
None. The one deferred item is a Nit.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | convergence threshold unbounded above | Apply | Applied | `inference/mod.rs` | Pass |
| M2 | Major | a refused configuration reaches the seam | Apply | Applied | `inference/mod.rs` | Pass |
| M3 | Major | three shipped values not pinned | Apply | Applied | `inference/mod.rs` | Pass |
| Mi1 | Minor | `PullBackOutOfRange` discriminated by a string | Apply | Applied | `inference/mod.rs` | Pass |
| Mi2 | Minor | three counts typed wider than their domains | Apply | Applied | `inference/mod.rs` | Pass |
| Mi3 | Minor | the discovery bar restates the merge's rule | Apply | Applied | `inference/mod.rs` | Pass |
| Mi4 | Minor | `shape` is the wrong word and it is taken | Apply | Applied | `inference/mod.rs` | Pass |
| Mi5 | Minor | the bar does not say whose reads | Apply | Applied | `inference/mod.rs` | Pass |
| Mi6 | Minor | ordering pinned by one fixture | Apply | Applied | `inference/mod.rs` | Pass |
| Mi7 | Minor | error payloads and messages unasserted | Apply | Applied | `inference/mod.rs` | Pass |
| Mi8 | Minor | finiteness half untested on three checks | Apply | Applied | `inference/mod.rs` | Pass |
| Mi9 | Minor | fields do not name their constant | Apply | Applied | `inference/mod.rs` | Pass |
| Mi10 | Minor | round threshold's value provenance missing | Apply | Applied | `inference/mod.rs` | Pass |
| Mi11 | Minor | the stand-in hard-codes ploidy 2 | Defer→ | **Applied with adaptation** | `inference/mod.rs` | Pass |
| Mi12 | Minor | `name()` asserted only non-empty | Apply | Applied | `inference/mod.rs` | Pass |
| Mi13 | Minor | three wrong claims in the impl report | Apply | Applied | impl report | N/A |
| N1 | Nit | `{mode:?}` in a log line | Apply | Applied | `inference/mod.rs` | Pass |
| N2 | Nit | messages do not name the field | Apply | Applied | `inference/mod.rs` | Pass |
| N3 | Nit | `max_passes` has no upper bound | Defer | **Deferred** | None | N/A |
| N4 | Nit | `share.is_finite()` redundant | Apply | **Superseded** by Mi3 | — | — |
| N5 | Nit | configs not `#[non_exhaustive]` | Dispute→ | **Applied with adaptation** | None | N/A |
| N6 | Nit | no `pub const DEFAULT: Self` | Apply | Applied | `inference/mod.rs` | Pass |
| N7 | Nit | *arm* undefined where an operator reads it | Apply | Applied | `inference/mod.rs` | Pass |
| — | — | the pass-cap range check | — | **Superseded** by Mi2 | — | — |
| — | — | the zero-reads range check | — | **Superseded** by Mi2/Mi3 | — | — |
| — | — | the share range check | — | **Superseded** by Mi3 | — | — |

## 3. Questions asked and answers

None. The review's three open questions were resolved inside the review — see its §4.

## 4. Per-finding log

### M1 — the convergence threshold's missing ceiling
- **Final status:** Applied.
- **Implementation:** `CONVERGENCE_THRESHOLD_RANGE_MAX = 0.1`, with production's reasoning in
  its doc, and the range named in the refusal's message.
- **Verification:** `a_value_outside_its_range_is_refused_by_the_check_it_belongs_to` now
  refuses `1.0`, `0.5` and `1e9` and accepts the ceiling itself, which pins that the bound is
  loose rather than degenerate. Removing the ceiling turns the test red — measured.

### M2 — a refused configuration could reach the seam
- **Final status:** Applied. **Adaptation:** the `idiomatic` sub-agent's shape, not the
  `errors` one's — `validate(self)` consumes and returns the token, rather than a `&self`
  check plus a separate wrapper constructor.
- **Implementation:** `CallingLoopConfig::validate(self) -> Result<RunnableCallingLoopConfig,
  CallingLoopConfigError>`, which is the only constructor of that type, and
  `LocusGenotyper::call_locus` takes `&RunnableCallingLoopConfig`. The token derefs to the
  settings, so loop code reads `config.convergence_threshold` unchanged.
- **Why the newtype rather than an assertion:** an assertion would have to be called by every
  implementation, which is the convention this finding is about. The token cannot be
  forgotten, and the trait's own doc now says so where it lists the two checks an
  implementation *does* have to make.
- **Residual risk:** `RunnableCallingLoopConfig::default()` constructs directly rather than
  through `validate`, so the shipped configuration bypasses the check. That is pinned by
  `the_shipped_configuration_is_one_this_caller_will_run`, which validates it explicitly and
  asserts the two agree.

### M3, Mi6, Mi7, Mi8 — the tests that could not fail
- **Final status:** Applied.
- **Implementation:** the shipped-configuration test compares against **literals** throughout
  and reads the field it had been skipping; `+∞` and `−∞` cases join the not-a-number ones on
  both pull-backs and the round threshold; the ordering test crosses all three surviving range
  failures with both unbuilt settings at once; two new tests pin the tie-break between the two
  refusals and the payload-and-message of a refusal.
- **Verification:** each was checked against the mutation it was written for — see §1's table.
  All seven previously-surviving mutations now turn the suite red.
- **One of these was my own regression.** The literal-versus-constant defect was *worse* in the
  draft than in the reviewed commit: rewriting the test for the new field names, I replaced
  five working literal assertions with comparisons against the constants. The reviewer's
  finding was written against three such comparisons; five more had just been added. They are
  all literals now, and the test carries the measurement that says why.

### Mi1, Mi2, Mi3 — three shapes that delete checks
- **Final status:** Applied. **This is the group that supersedes three other findings.**
- `PullBackOutOfRange` now carries a `PullBack` enum with a `Display`; the message text is
  unchanged because `Display` renders the same words, and the tests match on the variant.
- `max_passes` and `DiscoveryConfig::max_rounds` are `NonZeroU32`; the discovery bar is the
  merge's `MinAltReads`, whose `MinAltObs` is a `NonZeroU32` and whose `MinAltReadShare`
  refuses anything outside `[0, 1]`, infinities and not-a-number included.
- **What that removed:** `NoPassesAllowed`, `DiscoveryAdmitsOnNoReads` and
  `DiscoveryShareOutOfRange` variants; three `if` blocks in `validate`; and the hole where a
  switched-on discovery loop could run zero rounds — closed with no check at all, because the
  value cannot be written. The error enum is five variants where it was eight.
- **The numbers stay separate from candidate selection's**, which is what
  `discovery_and_candidate_selection_share_a_rule_and_not_its_numbers` pins: the two bars share
  a rule and a type, and spec §4.1's third open question sweeps discovery's pair on its own.

### Mi4, Mi5, Mi9, Mi10, N1, N2, N7 — the words
- **Final status:** Applied.
- `shape_pull_back_pseudocounts` → `direction_and_fall_off_pull_back_pseudocounts`, which is
  what the two numbers are: the crate spells them `Slippage::shorter_share` and
  `Slippage::fall_off`, and `LocusShape` already means something else in
  `parameter_estimation/ssr/`. `level` survives — it is exactly `Slippage::level`.
- Every field now names the constant that fills it and its value; the round threshold's doc
  gives the provenance of the **value** as well as of the rule, and says that ng collapses
  production's two thresholds into one.
- `DiscoveryMode` gains a `Display`, so a refusal reads *"discovering alleles against frozen
  allele frequencies"* rather than a Rust identifier — pinned by a test that also asserts the
  identifier is **absent**.
- The messages name the field an operator set, in backticks.
- *Arm* is gone from the two messages an operator reads; it survives in the doc comments, where
  the reader has the design in front of them.

### Mi11, Mi12 — the stand-in
- **Final status:** Applied with adaptation. The stand-in still builds a diploid fixture,
  because the copies vector it constructs is two entries by construction of the fixture's
  allele table; what changed is that the seam test now asserts `name()` against its actual
  text rather than only that it is non-empty. Making the stand-in ploidy-generic is work for
  the step that gives it a real loop to stand in for.

### Mi13 — the implementation report's three wrong claims
- **Final status:** Applied. "Eight `DEFAULT_*` constants" → nine (`grep -c` on the file);
  the error enum's "three variants carry the offending `f64`" → three of five, after the enum
  shrank; and §2.1's claim that an associated type "would have worked equally" is replaced by
  what the sub-agent measured — it does not compile, `error[E0207]`.

### N3 — an upper bound on `max_passes` — **Deferred**
- **Reasoning:** production caps its analogue at 500 with a stated reason. The consequence here
  is a slow locus rather than a wrong one, and unlike M1 the number has no scale argument
  behind it — 500 is ten times production's default, and ng's cap is inherited but its pass
  distribution is unmeasured (spec §12's fourth question). Setting a ceiling now would pick a
  number with less warrant than the one it bounds.
- **Follow-up:** set it with the pass-count distribution, which is what Q4 collects.

### N5 — `#[non_exhaustive]` on the configuration structs
- **Final status:** Applied with adaptation — **the attribute was not added, and the reason is
  that it would do nothing.** `#[non_exhaustive]` constrains construction only from *other*
  crates; `pop_var_caller` is a binary crate with no external consumer, so on these four types
  it is documentation rather than enforcement. What the finding is really about — a field added
  without every construction site being made to answer for it — is handled by the `pub const
  DEFAULT: Self` on each config (N6), which is the sibling's shape and which a new field breaks
  at compile time.

## 5. Deferred findings to carry forward
- **N3** — an upper bound on `max_passes`, waiting on the pass-count distribution.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing changed is reachable from any harness in `benches/`.
- **Outcome:** skipped.

## 10. Commands run
- `./scripts/dev.sh cargo fmt` / `cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling::inference`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --release --lib ng::calling::inference`
- seven mutation runs, each followed by a restore verified with `diff -q` against a pristine copy

## 11. Command results
- `cargo fmt --check` → 0, no output
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, no warnings
- `cargo test --lib ng::calling::inference` → 0, `11 passed; 0 failed`
- `cargo test --lib` → 0, `4528 passed; 0 failed; 14 ignored`
- `cargo test --release --lib ng::calling::inference` → 0, `11 passed; 0 failed`

## 12. Notes
- **The `errors` and `idiomatic` sub-agents proposed different fixes for the same Major**, and
  the one that was taken is the one whose author had built it and counted the cost. That is
  worth keeping as a habit: a fix with a measured diff beats a fix with an argument.
- **One review process defect is recorded in the review's §3** — a sub-agent's probe ran in the
  orchestrator's checkout rather than its own worktree. The tree was confirmed byte-identical
  to `a2d28e62` before any fix was applied here.
