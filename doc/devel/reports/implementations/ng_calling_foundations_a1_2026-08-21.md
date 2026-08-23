# ng calling foundations — A1: `AlleleId` and `Phred`

*Implementation report, 2026-08-21. Branch `ng-calling-foundations`, worktree
`../pop_var_caller-calling-foundations`. Step A1 of
[`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md), Milestone A.*

## 1. Plan

Add the first two scalars the calling step needs to ng's shared vocabulary,
[`src/ng/types.rs`](../../../../src/ng/types.rs), following the conventions that file already
keeps and appending at the end of the sections they belong to so the parallel
`ng-calling-prerequisites` branch (which edits only the `InbreedingF` block and inserts a
`DomainError` variant beside `InbreedingF`'s) cannot conflict:

- `AlleleId(pub u16)` — an index into one locus's candidate-allele table. Unconstrained, so a
  public field and no checked constructor, like `ContigId` and `ReadGroupId`.
- `Phred(f32)` — a quality on the Phred scale. Constrained (finite, `>= 0`), so a private field
  and a checked `try_new`, like `MismatchFraction`; plus `from_log_prob`, the named crossing from
  ng's internal natural-log currency, and a `DomainError::Phred` variant appended at the end of
  that enum.

Design authority: [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §Module home and
§2 (which own both names), and [`arch/ng_step_interfaces.md`](../../ng/arch/ng_step_interfaces.md)
§1 (the newtype and validation conventions).

## 2. Assumptions — what the design left open and what was chosen

**`from_log_prob` returns `Result`, and the two values the scale cannot hold are rejected, not
clamped.** The architecture fixes the name and says "validated `>= 0`, conversions named", and it
does not say what the conversion does with `ln p = -∞` — a legal `LogProb` in this codebase, the
documented score of an impossible read line-up — whose Phred is infinite. Nor with `ln p > 0`, a
probability above one, whose Phred is negative.

The choice made: both go through `try_new` and come back as `Err`. Nothing new is invented — the
type states once what a legal quality is, and the conversion inherits it. The alternative, capping
inside the type, was rejected because the cap is a run-level configuration everywhere it already
exists: production caps `GQ` at `max_gq_phred` at the point it fills the column
([`posterior_engine.rs:977`](../../../../src/var_calling/posterior_engine.rs)), not inside a
scalar. A clamp here would silently pick a ceiling for every future consumer.

`ng_step_interfaces.md` §1's validation policy makes `try_new -> Result` the rule for *untrusted*
input and allows a `debug_assert!`-plus-epsilon-clamp `new` for *internally computed* values. A
`LogProb` reaching `from_log_prob` is internally computed, so the policy would permit the softer
door; the harder one was taken because neither rejected case is a float-epsilon overrun — `-∞` and
a positive logarithm are both categorically different numbers, and the policy's own words are that
a gross out-of-range value stays a loud bug.

**No reverse conversion.** `Phred` → `LogProb` is not built, because no consumer exists yet and
the plan names only `from_log_prob`. Recorded here so a later step adds it deliberately rather
than finding the gap under time pressure.

**Where the two types were placed in the file.** Both were appended at the end of the "Scalar
newtypes" section, in the plan's own order (`AlleleId`, then `Phred`), rather than each being
inserted beside its nearest relative (`ReadGroupId` for the one, `LogProb` for the other). The
plan's region discipline asks for appends; an insert in the middle of the section would work
today but is exactly the shape that conflicts when two branches both do it.

## 3. Changes made

One file, purely additive: `src/ng/types.rs`, **+184 / −0** (`git diff --stat`).

- **`AlleleId(pub u16)`** with the file's ergonomic derives and `#[inline] get()`. The doc comment
  carries two contracts: index 0 is the reference allele, and an id is meaningless away from the
  locus it was minted at — the same relationship `Position` has to `ContigId`, which the file
  already explains. The width choice is justified against production's caps: 6 candidate alleles
  per record by default, refused above 16
  ([`per_group_merger.rs:57`](../../../../src/var_calling/per_group_merger.rs),
  [`:139`](../../../../src/var_calling/per_group_merger.rs), enforced at
  [`:162`](../../../../src/var_calling/per_group_merger.rs)).
- **`Phred(f32)`**, private field, with:
  - `try_new(f32) -> Result<Self, DomainError>` — `q >= 0.0 && q.is_finite()`. `NaN` is rejected
    by the first conjunct alone, since no comparison with `NaN` is true; `+∞` by the second.
  - `from_log_prob(LogProb) -> Result<Self, DomainError>` — scales in `f64`, the width `LogProb`
    holds, and narrows to `f32` once at the end, then routes through `try_new`.
  - `#[inline] get() -> f32`.
- **`PHRED_PER_NAT: f64 = 10.0 / LN_10`**, private to the module.
- **`DomainError::Phred(f32)`**, appended at the end of the enum.

### The one deviation worth recording: the constant's sign and its namesake

`PHRED_PER_NAT` was first written as `-10.0 / LN_10`, with the negation folded into the constant.
It was changed to the positive `10.0 / LN_10` with the negation at the call site
(`-PHRED_PER_NAT * log_p.get()`) after finding that **a constant of exactly this name already
exists**: `baq::probaln::PHRED_PER_NAT`, `pub(super)`, holding `4.343`
([`probaln.rs:33`](../../../../src/baq/probaln.rs)). That one is the four-digit truncation htslib
compiled, deliberately kept truncated for byte-parity with htslib and guarded by its own test
([`baq/tests.rs:118`](../../../../src/baq/tests.rs)); its call sites spell the conversion
`-PHRED_PER_NAT * p.ln()`. Two constants sharing a name while differing in *sign* would have been
a trap; sharing a name while differing only in *precision*, with the reason stated in both doc
comments, is not. The new constant's doc comment names the collision and says not to unify them.

## 4. Tests added

Three new tests in `src/ng/types.rs`'s `mod tests`, plus two assertions added to the existing
`unconstrained_newtypes_expose_their_value`:

| test | what it pins |
|---|---|
| `phred_accepts_zero_and_rejects_everything_below_it` | The boundary in both directions: `0.0` is legal (`p = 1`, a call that cannot be wrong) and `-f32::EPSILON` is not; `-1.0`, `NaN` and both infinities are rejected. `NaN` is asserted with `matches!` and `q.is_nan()`, never `assert_eq!`, because `DomainError`'s `PartialEq` is IEEE equality on the float payload — the enum's own doc comment states this. |
| `phred_from_log_prob_matches_the_hand_computed_scale` | The conversion against numbers worked out by hand rather than by the same formula: 30 at one wrong call in a thousand, 20 at one in a hundred, and `10 log10 2 = 3.0103` at a half. Tolerance `1e-4`. |
| `phred_from_log_prob_rejects_what_the_scale_cannot_hold` | `ln 1 = 0` gives quality zero and must **not** be an error; `ln 0 = -∞` errors with a `+∞` payload; a positive logarithm errors with a negative payload. |
| `unconstrained_newtypes_expose_their_value` (extended) | `AlleleId(0)` and `AlleleId(u16::MAX)` round-trip through `get()`. |

## 5. Validation

Run in the dev container (`./scripts/dev.sh`, Apple `container` on macOS), verbatim:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 8.16s` |
| `cargo clippy --all-targets --all-features -- -D warnings` | 101 | **18 errors, every one pre-existing** — see below |
| `cargo test --all-targets --all-features` | 101 | lib `test result: ok. 3941 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 676.01s`; every integration-test binary ok; **one pre-existing bench panic** — see below |

The four in-scope tests pass:

```
test ng::types::tests::phred_accepts_zero_and_rejects_everything_below_it ... ok
test ng::types::tests::phred_from_log_prob_matches_the_hand_computed_scale ... ok
test ng::types::tests::phred_from_log_prob_rejects_what_the_scale_cannot_hold ... ok
test ng::types::tests::unconstrained_newtypes_expose_their_value ... ok
```

### The two aggregate gates that were already red on `main`

Neither is this branch's doing, and both were checked rather than assumed.

**Clippy.** All 18 errors cite three files, none of which this branch touches — `git diff --stat`
reports exactly one changed file, `src/ng/types.rs`:

```
benches/cohort_var_calling_perf.rs:110, :170
benches/ng_joint_fit_perf.rs:335
examples/ng_joint_contamination_harness.rs:292, :462, :463, :489, :537, :575, :612, :613, :629, :634, :639, :678
```

The library alone is clippy-clean at `-D warnings`.

**Tests.** `cargo test --all-targets` runs each bench's harness, and one of them panics:

```
thread 'main' (9093) panicked at benches/psp_writer_perf.rs:386:60:
index out of bounds: the len is 3300000 but the index is 3300000
error: test failed, to rerun pass `--bench psp_writer_perf`
```

Verified pre-existing rather than argued: the change was stashed out
(`git stash push src/ng/types.rs`) and `cargo test --bench psp_writer_perf` re-run, giving the
identical panic at the identical line and exit 101. The bench's setup loop primes records until a
projected-byte target is reached, exhausts all 3,300,000 of them without reaching it, and the body
then indexes one past the end; the file was last touched in commit `21dac5bb`, long before this
branch.

## 6. Trade-offs and follow-ups

- **No reverse conversion** (`Phred` → `LogProb`) — see §2. Add it when a consumer needs it.
- **Nothing consumes either type yet.** The calling module arrives at step B1; until then these
  are exercised only by their own tests.
- **`AlleleId`'s reference-at-zero contract lives in prose, not in the type.** Enforcing it is
  `CandidateAlleles`' job at step B2, which is where the table itself is built.
- **The two red aggregate gates** above are unowned by this plan. Fixing three unrelated
  benches/examples is not in its scope; raised for the owner at Checkpoint A.
