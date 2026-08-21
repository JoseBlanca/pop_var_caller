# Applying the A1 review — ng calling foundations, step A1

*2026-08-21. Branch `ng-calling-foundations`. Input:
[`ng_calling_a1_2026-08-21.md`](../reviews/ng_calling_a1_2026-08-21.md). Every finding in that
report is accounted for below; nothing is silently dropped.*

## Findings table

| id | severity | decision | status |
|---|---|---|---|
| M1 — `from_log_prob` stores negative zero | Major | Apply | **Applied** |
| M2 — no property test, alone among the file's constrained newtypes | Major | Apply | **Applied** |
| M3 — one error variant for three conditions | Major | Apply | **Applied** |
| M4 — the two deferred decisions have no names | Major | Apply in part | **Applied with adaptation** (doc half applied; constructor half deferred) |
| Mi1 — the `f64`-width invariant has no test | Minor | Apply | **Applied** |
| Mi2 — a third rejection cause, undocumented and untested | Minor | Apply | **Applied** |
| Mi3 — no test for a `NaN` log probability | Minor | Apply | **Applied** |
| Mi4 — the top of the range is uncovered | Minor | Apply | **Applied** |
| Mi5 — the `Err` names the quality, not the log probability | Minor | Apply | **Applied** (documentation) |
| Mi6 — "the owned allele multiset" is not a name the reader can look up | Minor | Apply | **Applied** |
| Mi7 — index `0` means the reference allele in prose only | Minor | Apply | **Applied** |
| Mi8 — the doc forbids an `as` cast 42 lines above one | Minor | Apply | **Applied** |
| Nit — `try_new(q)` too terse | Nit | Apply | **Applied** (`quality`) |
| Nit — "The only constructor." is not true here | Nit | Apply | **Applied** |
| Nit — *nat* used before defined | Nit | Apply | **Applied** |
| Nit — "truncated" is the wrong operation for `4.343` | Nit | Apply with adaptation | **Applied with adaptation** |
| Nit — the error message is the only capitalised one of fifteen | Nit | Apply | **Applied** |
| Nit — "needs no test of its own" collides with `#[test]` | Nit | Apply | **Applied** |
| Nit — the scale test credits the tolerance, not the pairs that do the work | Nit | Apply | **Applied** |
| Nit — the infinity assertions bind their payload with `_` | Nit | Apply | **Applied** |
| Nit — `PHRED_PER_NAT` as an associated const | Nit | Dispute | **Won't fix** |
| Nit — a hand-written `Ord` on `Phred` | Nit | Defer | **Deferred** |
| Nit — a `Display` impl for `Phred` | Nit | Defer | **Deferred** |
| Nit — no `///` on the two `get()` accessors | Nit | Dispute | **Won't fix** |
| Cross-cat — `PhredQual` would match `MapQual`/`BaseQual` | — | Dispute | **Won't fix** |
| Cross-cat — no reverse `Phred` → `LogProb` conversion | — | Defer | **Deferred** |
| Out of scope — three red gates on `main`, the `AlleleId` table-read promise | — | Defer | **Deferred** (§4) |

## What changed, and how each was verified

### M1 — the negative zero

`try_new` now normalises the sign of zero, so every route in agrees:

```rust
.then_some(Self(if quality == 0.0 { 0.0 } else { quality }))
```

Pinned by a new test, `phred_zero_is_positive_zero_whichever_constructor_made_it`, which asserts
`is_sign_positive()` and that the value formats as `"0"` — an `assert_eq!(.., 0.0)` cannot see this,
because `-0.0 == 0.0`.

**Verified by re-running the mutation that survived the review.** Removing the normalisation
(`.then_some(Self(quality))`) and running `./scripts/dev.sh cargo test --lib ng::types`:

```
test ng::types::tests::phred_zero_is_positive_zero_whichever_constructor_made_it ... FAILED
test ng::types::tests::phred_accepts_exactly_the_finite_non_negative_values_and_round_trips ... FAILED
test result: FAILED. 39 passed; 2 failed; 0 ignored; 0 measured; 3915 filtered out
```

Two tests kill it where the submitted suite killed none. The file was restored afterwards.

### M2 — the property test

`phred_accepts_exactly_the_finite_non_negative_values_and_round_trips` was added inside the file's
existing `proptest::proptest! { … }` block, modelled on
`the_constrained_rates_accept_exactly_the_probabilities_and_round_trip` beside it: accepted
**exactly when** finite and at or above zero, and the accepted value back **bit for bit**, with
`-0.0 → +0.0` as the one stated normalisation. Same dense-arm construction as the rates', for the
same reason — `f32::ANY` alone essentially never samples next to zero.

It is one of the two tests that kill M1's mutation, above.

### M3 — the error variant split

`DomainError::PhredInfinite` was appended after `DomainError::Phred`, both at the end of the enum,
so the parallel `ng-calling-prerequisites` branch's insertion beside `InbreedingF` stays disjoint
and no discriminant moves. `try_new` returns it for `+∞` before the range test; everything else
keeps `DomainError::Phred(quality)`.

The doc contradiction is gone: `try_new` now says a negative or a `NaN` is broken arithmetic while
an infinity is a probability of exactly zero, "a different event"; `from_log_prob` lists its two
causes against their two variants.

**Verified by mutation.** Deleting the `PhredInfinite` early return so `+∞` falls back to
`DomainError::Phred`:

```
test ng::types::tests::phred_accepts_zero_and_rejects_everything_below_it ... FAILED
test ng::types::tests::phred_from_log_prob_rejects_what_the_scale_cannot_hold ... FAILED
test result: FAILED. 39 passed; 2 failed; 0 ignored; 0 measured; 3915 filtered out
```

**This is one variant more than the plan's step A1 wrote** ("a new `DomainError` variant",
singular). Recorded as a review-driven deviation: it adds no policy, changes no existing behaviour,
and nothing consumes either variant yet — which the `errors` agent named as the cheapest moment for
the split, since after the calling module lands the type would be changing under its consumers.

### M4 — applied in part

**Applied:** the type doc now names the constants instead of the module —
`DEFAULT_MAX_GQ_PHRED` (99, the GATK and bcftools convention) and `GQ_PHRED_RANGE_MAX` (200),
both as intra-doc links, so the doc and the values cannot drift. It also records, in the type where
a consumer will read it, that a quality ng computed itself has no constructor here yet and that the
step which first fills a `GQ` column is where one belongs.

**Deferred:** the `DEFAULT_MAX_QUALITY` constant and the `from_log_prob_capped(log_p, cap)`
constructor. Both decide *what ceiling ng writes*, which is run policy this plan does not own —
`arch/calling_em_loop.md` §2.1 lists the calling loop's configuration and no quality ceiling is in
it, and `calling_loop.md` owns `CallingLoopConfig`. Building the constructor now would fix a
ceiling before the step that needs one exists. Carried as a follow-up on that plan.

### Mi1 to Mi4 — the four coverage gaps

- `phred_from_log_prob_keeps_full_f64_width_before_narrowing` pins the width invariant at a quality
  of 3000, where the ordering shows. **Verified by mutation:** replacing the body with
  `-(PHRED_PER_NAT as f32) * (log_p.get() as f32)` gives `test result: FAILED. 40 passed; 1 failed`
  — that one test, where the submitted suite killed nothing.
- The `f32`-overflow rejection (Mi2) is now documented on `from_log_prob` and asserted inside
  `phred_from_log_prob_rejects_what_the_scale_cannot_hold` (`LogProb(-1e300)` → `PhredInfinite`).
- The `NaN` log probability (Mi3) is asserted in the same test, with the reason it is reachable:
  `LogProb`'s field is public and unconstrained.
- `try_new(f32::MAX)` (Mi4) is asserted inside `phred_accepts_zero_and_rejects_everything_below_it`,
  whose doc now says why the top of the range is deliberately open — a ceiling added inside the type
  later must break a test.

The last three were folded into existing tests rather than added as three more `#[test]`s: each is
one assertion about the same function, and the file's style is a test per behaviour rather than per
input.

### Mi5 to Mi8 and the nits

Documentation and naming, no behaviour change: the `Err` payload and the saturating narrowing are
now stated on `from_log_prob`; `AlleleId`'s doc names `Genotype`; the "never an `as` cast" sentence
now says "never a bare `as` cast between the two scales" and explains what the narrowing inside the
function is; `q` became `quality` in the signature, the prose and the tests; "The only constructor."
became "The one check, and every constructor goes through it."; *nat* is glossed as "per unit of
natural logarithm"; the message is lowercase like its fourteen siblings; "needs no test of its own"
became "needs no check of its own"; and the scale test's comment now names the two decade pairs as
the discriminating ones and quotes the three measured wrong-implementation values (99.657845, 3,
30.000381) rather than crediting the tolerance in general.

**Mi7** added `AlleleId::REFERENCE` and `AlleleId::is_reference()`, with
`allele_id_zero_is_the_reference_allele` pinning both. The public field was kept — the review
agreed it should be. No `Default` impl: a no-argument `AlleleId::default()` silently meaning "the
reference allele" is exactly the hidden significant default the category forbids.

**The "truncated" nit was applied with adaptation.** Two agents noted that `4.343` is
`4.342944819…` *rounded* to four figures, not truncated, and one warned that fixing the word here
alone would make this comment disagree with `baq/probaln.rs`'s own ("4-digit truncation"). The
sentence now avoids the verb entirely — "the same quantity as the four-digit literal htslib
compiled, `4.343`" — so it is accurate without contradicting the file it points at. Correcting
`probaln.rs` is a separate, one-word follow-up on a file this branch does not touch.

**Also fixed, though it is pre-existing:** `DomainError`'s own doc said "IEEE equality on the `f64`
payloads" while two variants carry `f32`. `MismatchFraction(f32)` predates this change, but
`Phred(f32)` is the second, so the word is now "float".

### The four not applied

- **`PHRED_PER_NAT` as an associated const on `Phred`** — disputed. The reviewer's own case for it
  was that the path would then distinguish it from `baq::probaln::PHRED_PER_NAT`; the same reviewer
  established that the two cannot collide, because one is a private module const and the other
  `pub(super)` inside `baq`, so no scope can name both. No behaviour change, no reader benefit that
  the "do not unify them" note does not already give.
- **A hand-written `Ord` on `Phred`** — deferred. The constructor's invariant does rule out `NaN`,
  so a total order is available, but the file's rule is that float newtypes stop at `PartialOrd` and
  no consumer sorts qualities yet. Flagged for the calling loop.
- **A `Display` impl** — deferred for the same reason: `Ploidy` and `SsrPeriod` have one because a
  message renders them today; nothing renders a `Phred` yet.
- **Renaming `Phred` to `PhredQual`** — won't fix. The reviewer flagged rather than recommended it,
  and correctly: the name is fixed by `arch/ng_step_interfaces.md` §1 and
  `arch/calling_em_loop.md` §0, so changing it is a decision about those documents, not a local
  edit.
- **No `///` on the two `get()` accessors** — disputed. All ten existing accessors in this file are
  the same; adding two would be noise and would leave the file inconsistent.

## Validation

Run in the dev container after every fix above was in place:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | clean |
| `cargo doc --no-deps --lib` | 101 | 17 unresolved links, **all pre-existing**; the only one in `src/ng/types.rs` is `SsrSegment` at line 808, in the Motif section this branch does not touch. Every new intra-doc link resolves. |
| `cargo test --lib ng::types` | 0 | `41 passed; 0 failed; 0 ignored; 0 measured; 3915 filtered out` |
| `cargo test --all-targets --all-features` | 101 | lib `test result: ok. 3945 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 666.94s` — four more than the 3,941 the review saw; every integration-test binary ok; the run then hits the same **pre-existing** `benches/psp_writer_perf.rs:386` panic §7 of the review records |

The step's whole surface, after the fixes — seven tests, all passing:

```
test ng::types::tests::allele_id_zero_is_the_reference_allele ... ok
test ng::types::tests::phred_accepts_exactly_the_finite_non_negative_values_and_round_trips ... ok
test ng::types::tests::phred_accepts_zero_and_rejects_everything_below_it ... ok
test ng::types::tests::phred_from_log_prob_keeps_full_f64_width_before_narrowing ... ok
test ng::types::tests::phred_from_log_prob_matches_the_hand_computed_scale ... ok
test ng::types::tests::phred_from_log_prob_rejects_what_the_scale_cannot_hold ... ok
test ng::types::tests::phred_zero_is_positive_zero_whichever_constructor_made_it ... ok
```

The change is now **+357 / −1** in `src/ng/types.rs` (`git diff --stat`), against the **+184 / −0**
the review was given. The one deletion is the `f64` → `float` word in `DomainError`'s doc comment.

The three mutation runs quoted above were each followed by restoring the file from a copy; the one
side effect they left — a new seed line in `proptest-regressions/ng/types.txt` recording the
negative-zero counterexample from the mutant run — was reverted with `git checkout`, so the working
tree holds only the intended change.

## Follow-ups this run created

1. **A named quality ceiling and a capping constructor** (`DEFAULT_MAX_QUALITY`,
   `from_log_prob_capped`) — for the step that first fills a `GQ` column;
   [`calling_loop.md`](../../ng/impl_plan/calling_loop.md).
2. **`CandidateAlleles`' accessor must not be a bare index** — `AlleleId`'s doc promises an
   out-of-range id is caught when the table is read. Step B2.
3. **A reverse `Phred` → `LogProb` conversion** — when a consumer needs one.
4. **`baq/probaln.rs:30`'s "4-digit truncation"** should read "rounding". One word, on a file this
   branch does not touch.
