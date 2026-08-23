# Code review — ng calling foundations, step A1 (`AlleleId`, `Phred`)

*2026-08-21. Branch `ng-calling-foundations`, reviewed at base commit `ee62a518` with the step's
uncommitted working-tree diff applied. Five category sub-agents, each in its own git worktree.
Per-category audit trail in the gitignored `tmp/review_2026-08-21_ng-calling-a1/`.*

## 1. Scope

**What was reviewed:** the working-tree diff of step A1 of
[`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md) — a single-file, purely
additive change to [`src/ng/types.rs`](../../../../src/ng/types.rs), 184 insertions and 0
deletions at the time of review.

**In-scope regions:** the `AlleleId` block, the `Phred` block, the `PHRED_PER_NAT` constant, the
`DomainError::Phred` variant, and the test additions in `mod tests`.

**Deliberately out of scope:** the rest of the repository, and pre-existing code in
`src/ng/types.rs` outside the added regions except where the additions had to match its
conventions.

**Categories dispatched, and why the other six were not.** The diff is 184 additive lines in one
file, defines two newtypes with no callers, and changes no existing behaviour, so a full
eleven-category fan-out would have cost more than it could return. Dispatched:

| category | reason |
|---|---|
| `reliability` | always; and the only category that mutation-tests |
| `errors` | always; the change adds a `DomainError` variant and two fallible constructors |
| `naming` | always; the diff *is* vocabulary, and its doc comments are reader-facing prose |
| `idiomatic` + `smells` | always; **run by one agent**, since a 184-line diff does not warrant two fan-out slots. Findings kept under separate headings |
| `defaults` | the change adds public API and a behaviourally significant reject-don't-clamp contract |

Not dispatched, with reasons: `module_structure` (one file, nothing moved), `unsafe_concurrency`
(no `unsafe`, `Arc`, `Mutex`, atomics, channels or `async` — the crate `forbid`s `unsafe_code`),
`tooling` (`Cargo.toml` untouched), `extras` (no parser, no untrusted input, no hot path, not a
published crate), `refactor_safety` (purely additive, zero callers, no existing behaviour altered).

## 2. Verdict

**Approve-with-changes.** No Blockers. Four Majors, of which three are applied in full and one is
applied in part with the remainder deferred to the step that first fills a `GQ` column. The
mutation pass is the reason this is not an Approve: two mutations survived the suite as submitted,
and one of them changed a value that would reach a VCF column.

## 3. Execution status

Commands run by the orchestrator in the dev container, and passed verbatim into every sub-agent's
prompt so none re-ran them:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 8.16s` |
| `cargo clippy --all-targets --all-features -- -D warnings` | 101 | 18 errors, **all pre-existing** (§7) |
| `cargo test --all-targets --all-features` | 101 | lib `3941 passed; 0 failed; 11 ignored`, all integration binaries ok, then **one pre-existing bench panic** (§7) |

Two sub-agents built and ran the library tests in their own worktrees (`37 passed; 0 failed` on
the `ng::types` filter); `reliability` additionally ran a probe test printing every result as value
**and bit pattern**, and removed it afterwards.

Findings labelled "Needs verification": **zero**. Every finding below is backed by a quoted
measurement or by code the agent read.

## 4. Open questions and assumptions

1. **Is `Result` the right return for a value ng computes itself?** `arch/ng_step_interfaces.md` §1
   nominally points the other way for internally-computed values (a `new` that `debug_assert!`s and
   clamps a float-epsilon overrun). Resolved by the `errors` agent in favour of `Result`, and the
   reasoning is now in the code: a `debug_assert` would let a release build write `NaN` or `-4` into
   a `QUAL` column, and `ln p = -∞` is not an overrun of anything — it is a representable input with
   no image on the output scale, which the policy's own "a *gross* out-of-range value stays a loud
   bug" clause covers. Affects **M3**.
2. **Who owns the `GQ` ceiling?** Not this step. Affects **M4**, whose second half is deferred.

## 5. Top 3 priorities

1. **M1 — `from_log_prob` returns negative zero for a certain call.** Measured, and the test
   written for exactly that input could not see it. Would print as `-0` in a `QUAL` column.
2. **M3 — one error variant for three conditions, only two of which are bugs.** `ln p = -∞` is
   inside `LogProb`'s documented domain; a consumer had to do float forensics to tell a routine cap
   from a broken sum.
3. **M2 — `Phred` was the only constrained newtype in the file without a property test**, and the
   bit-for-bit half of that test is precisely what M1 slipped through.

## 6. Findings

### Major

**M1: `src/ng/types.rs:364` — `from_log_prob` stores negative zero at `p = 1`, and the test for that
input cannot see it.**
**Categories:** reliability, idiomatic (convergent — both agents found it independently, both by
building). **Confidence:** High.

`-PHRED_PER_NAT * 0.0` is `-0.0` under IEEE, and `try_new` accepted it because `-0.0 >= 0.0` is
true. Measured by the `reliability` agent's probe:

```
PROBE from_log_prob(ln 1) = Ok(-0e0) bits=0x80000000
PROBE try_new(-0e0)       = Ok(-0e0) bits=0x80000000
```

The assertion written for the case was `assert_eq!(…get(), 0.0)`, and `-0.0 == 0.0`, so it passed
either way. Confirmed by mutation: replacing `.then_some(Self(q))` with `.then_some(Self(q.abs()))`
flips the stored bits from `0x80000000` to `0x0` and **all 38 tests still passed**.

Why it matters: `Phred`'s stated purpose is VCF's `QUAL` and `GQ`, `format!("{}", -0.0f32)` is
`"-0"`, and `QUAL` is a non-negative field — so the one input the doc comment calls out as *the one
that must not be rejected* is the one that would print with a minus sign. It also made two `Phred`s
that compare equal behave differently under `Display` and `is_sign_negative`.

**M2: `src/ng/types.rs:333-370` — `Phred` was the fourth constrained newtype in this file and the
only one without a property test.** **Category:** reliability. **Confidence:** High.

`ErrorRate`, `GenotypeFrequency` and `InbreedingF` share one; `Ploidy` has one; `SsrPeriod` has
one. Each asserts *accepted exactly when in range* **and** *the accepted value comes back bit for
bit*. `Phred` had three point tests. The missing bit-for-bit half is exactly what M1 slipped
through, and the rates' own proptest doc comment already states the argument that applies here
verbatim: a widened bound on either side is a value reaching the genotype prior instead of an
error.

**M3: `src/ng/types.rs:697` — one `DomainError::Phred` variant for three conditions, and only two of
them are bugs.** **Category:** errors. **Confidence:** High on the collapse and the doc
contradiction; Medium on the remedy.

Three unrelated inputs landed on one variant: `ln p > 0` (probability above one — broken
arithmetic), `NaN` (likewise), and `ln p = -∞` — which is **not** broken. `LogProb`'s own doc says
so: "`f64::NEG_INFINITY` is a legal value, not an error … Every finite `f64` and `-∞` is a valid
log-probability". Production's answer to it is not to fail either: `DEFAULT_MAX_GQ_PHRED`'s doc
says the cap "also prevents `+∞` GQ when EM yields `P(best) = 1` exactly".

So a consumer could only separate *cap and carry on* from *abort* by testing the payload's sign and
finiteness — which is what the submitted tests were reduced to doing. The doc also contradicted
itself: `try_new` said "each of them says the caller's arithmetic went wrong" while
`from_log_prob`'s own test three hundred lines later called the same `-∞` deliberate.

**M4: `src/ng/types.rs:336-372` — `Phred` defers two decisions to its consumers and names neither.**
**Category:** defaults. **Confidence:** Medium. **Assumption:** that the first consumer is the one
`arch/calling_em_loop.md` names, `SampleGenotypeCall.genotype_quality`, a posterior-derived GQ.
Nothing consumes `Phred` in this diff, so the impact is a forecast about the next build order, not
an observed defect.

Two halves. The first: the type doc pointed at `var_calling::posterior_engine` for the cap without
naming the constant, which is the drift the category's magic-number rule forbids — a later change
to `DEFAULT_MAX_GQ_PHRED` would leave this doc describing a ceiling that had moved. The second:
`arch/ng_step_interfaces.md` §1's *internally computed* branch has no constructor here at all, so a
consumer holding a posterior-derived quality must supply production's three numbers — the
`1.0 - f64::EPSILON` pin, the `0.0` floor and the `99.0` ceiling — out of its own head, and `QUAL`
and `GQ` can end up capped at two different values with nothing in the type to compare them
against.

### Minor

- **Mi1 (reliability): the `f64`-width invariant had no test, and the mutation that violates it
  survived.** The doc comment promises scaling at `f64` and narrowing once at the end. Replacing
  the body with `-(PHRED_PER_NAT as f32) * (log_p.get() as f32)` left **all 38 tests passing**,
  though the arithmetic genuinely differs: `ln 1e-300` gives `3000.0` (`0x453b8000`) against the
  mutant's `2999.9998` (`0x453b7fff`). The three existing pairs sit at qualities 30, 20 and 3.0103,
  where the difference is far below the `1e-4` tolerance.
- **Mi2 (reliability, errors — convergent): a third rejection cause, undocumented and untested.**
  `as` from `f64` to `f32` saturates, so a finite log probability below about `-7.8e37` becomes
  `+∞` and is refused. Measured: `PROBE from_log_prob(-1e300) = Err(…inf…)`. Unreachable from real
  data — `ln(f64::MIN_POSITIVE)` is about `-745` — so a documentation and coverage gap, not a
  defect. Naming it stops a future reader "simplifying" the guard to `!q.is_nan()`.
- **Mi3 (reliability): no test for a `NaN` log probability.** `LogProb`'s field is public and
  unconstrained, so a caller's `0.0 / 0.0` reaches `from_log_prob` directly. The behaviour was
  already right; nothing pinned it.
- **Mi4 (reliability): the top of `Phred`'s range was uncovered.** The type doc argues at length
  that capping is the consumer's decision and must not happen inside the type — a contract with no
  test behind it, so a ceiling added later would break nothing.
- **Mi5 (errors): `from_log_prob`'s `Err` names the quality it computed, never the `LogProb` it was
  handed.** An operator reading `-0.4342944` has to divide by `10/ln(10)` to recover the `0.1` that
  was passed.
- **Mi6 (naming): `AlleleId`'s doc called the genotype "the owned allele multiset".** The term
  appears nowhere in `src/` or in either design authority; the thing has a fixed name, `Genotype`,
  and the arch sentence this comment paraphrases uses it. That sentence is the one carrying the
  argument for why a bare `u16` id is safe, so the reader who wants to check it should be able to
  grep the name.
- **Mi7 (defaults): `AlleleId`'s index `0` meant "the reference allele" in prose only.** The
  evidence was in the diff itself — `assert_eq!(AlleleId(0).get(), 0); // the reference allele's
  index`, a trailing comment doing a named constant's job in the type's first use. The public field
  is **not** a finding: it matches `ContigId`, `Position` and `ReadGroupId` and §1's rule for
  unconstrained newtypes.
- **Mi8 (smells): the type doc said crossings are "never an `as` cast", and 42 lines below, inside
  the named crossing, sat `as f32`.** The code is right — the rule governs the *scale* crossing,
  which the multiply performs and names, and std offers no non-`as` narrowing of an `f64` — but a
  reader checking doc against code met an apparent contradiction in the one place the doc was
  trying to make a rule stick.

### Nits

`try_new`'s parameter was `q`, following the file's one terse precedent (`MismatchFraction`'s `x`)
rather than its five explicit ones (`rate`, `frequency`, `coefficient`, `copies`, `bases`), and it
leaked into the doc prose and the tests. "The only constructor." was copied from five types where
it is literally true and is not true here. *Nat* was used before being defined. "Truncated" is the
wrong operation for `4.343`, which is `4.342944819…` **rounded** to four figures — inherited from
`baq/probaln.rs`'s own wording, so fixable here only by not using the verb. The error message was
the only one of the enum's fifteen starting with a capital. `try_new`'s doc said `NaN` "needs no
test of its own" where it meant *no check of its own*, and the test module then wrote one. The
`0.5` pair in the scale test discriminates nothing the two decade pairs do not — measured, it
passes under `PHRED_PER_NAT = 4.343` where both decade pairs fail — while the test's comment
credited "the tolerance" generally. The infinity assertions bound their payload with `_`, so they
would have passed on the wrong infinity.

Not filed, and deliberately: `PHRED_PER_NAT` as an associated const (optional, no behaviour
change); a hand-written `Ord` on `Phred` (its invariant rules out `NaN`, so a total order is
available — flagged for the calling loop if it starts sorting qualities); a `Display` impl (no
consumer); missing `///` on the two `get()` accessors (all ten existing accessors in the file are
the same).

## 7. Out of scope observations

- **The aggregate clippy gate is red on `main`.** 18 errors in `benches/cohort_var_calling_perf.rs`,
  `benches/ng_joint_fit_perf.rs` and `examples/ng_joint_contamination_harness.rs`. None is in a file
  this branch touches. Follow-up: a separate lint-only branch.
- **The aggregate test gate is red on `main`.** `cargo test --all-targets` runs each bench's
  harness and `benches/psp_writer_perf.rs:386` panics with `index out of bounds: the len is 3300000
  but the index is 3300000`. Verified pre-existing by stashing this diff out and re-running that
  bench alone: identical panic, identical line. Its setup loop primes records until a projected-byte
  target is reached, never reaches it, exhausts all 3.3 M, and the body reads one past the end.
  Follow-up: the psp writer's own branch.
- **`cargo doc --no-deps --lib` is red on `main`** for 17 pre-existing unresolved intra-doc links,
  one of them in `src/ng/types.rs` outside this diff. None of the new doc links is among them.
- **`DomainError`'s doc comment says "IEEE equality on the `f64` payloads"** and two variants carry
  `f32` (`MismatchFraction` predates this change). One word, fixed here since this change adds the
  second.
- **`AlleleId`'s doc promises that "an out-of-range id is caught when the table is read"** and
  nothing reads a table yet. Carry into step B2: `CandidateAlleles`' accessor should return an
  `Option`/`Result` rather than indexing, or that sentence becomes a slice-index panic.

## 8. Missing tests to add now

Supplied as complete code by the `reliability` agent, all six applied (three as new tests, three
folded into existing ones — see the fix-application report):
`phred_zero_is_positive_zero_whichever_constructor_made_it`,
`phred_accepts_exactly_the_finite_non_negative_values_and_round_trips`,
`phred_from_log_prob_keeps_full_f64_width_before_narrowing`,
`phred_from_log_prob_rejects_a_log_probability_too_large_for_f32`,
`phred_from_log_prob_rejects_a_nan_log_probability`,
`phred_try_new_accepts_the_largest_finite_quality`.

## 9. What's good

- **The mutation pass is what earned this review its findings.** 15 mutations run, 2 survived, 1
  changed no behaviour — and the three numbers are reported separately, so the no-op is not
  miscounted as a survivor. Both survivors became findings with a measured before/after.
- **The no-op mutation was kept and used as evidence.** Reordering `q >= 0.0` and `q.is_finite()`
  produced byte-identical probe output, which mechanically proves the doc comment's claim that
  `q >= 0.0` — not `is_finite()` — is what rejects `NaN`.
- **Every quantitative claim in the diff's own prose was re-derived rather than re-read**, including
  the two the test comment makes about wrong implementations: `log2` in place of `log10` returns
  99.657845 against the claimed "99.7", and dropping the factor of ten returns exactly 3.
- **Three agents independently checked the same two cross-file claims** (`baq`'s `4.343`,
  production's allele caps) and all three reached the same values from the same lines.
- **`naming` and `idiomatic` each argued a rule *down* where the code was right** — the
  `PHRED_PER_NAT` collision is disarmed by scoping, the `as f32` narrowing has no alternative in
  std — instead of filing a finding for the appearance of a violation.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --all-targets --all-features        # expect the pre-existing bench panic
./scripts/dev.sh cargo test --lib ng::types                     # the step's own surface
```
