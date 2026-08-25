# ng candidate alleles — A1: review and the fixes applied

*Review report, 2026-08-24. Branch `ng-candidate-alleles`. Scope: the working-tree diff of step
A1 of [`candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md) — three files, 197 added
lines, **no logic**. Reviewed at `3edab4cd` + the step patch, five agents each in its own
worktree. Fixes applied in the same pass; this file carries both.*

**Verdict: approve with changes, all applied.** 1 Blocker, 7 Major and 23 Minor findings as
filed — **three distinct defects above Minor** once the convergent filings are merged, because
three agents found the same two independently. All three are the same shape: *a wrong number that
produces plausible output instead of a crash*. **Every one was found by mutating the code, none
by reading it.** That is the failure class the plan's own principle ("isolate the silent steps")
names, arriving one milestone earlier than the plan expected it.

---

## 1. Which categories ran, and which did not

Five agents: **defaults**, **reliability**, **naming**, **idiomatic + errors**,
**module_structure**. Not dispatched: `refactor_safety`, `smells`, `tooling`,
`unsafe_concurrency`, `extras` — A1 adds no logic, no concurrency, no `Cargo.toml` change and no
parser, so each would have had nothing to bind to. The skill lists four of those as "always";
this is a recorded departure, not an oversight.

Findings by agent, as filed:

| agent | Blocker | Major | Minor | Nits |
|---|---|---|---|---|
| reliability | 1 | 2 | 2 | 1 |
| defaults | — | 3 | 5 | 4 |
| idiomatic + errors | — | 2 | 6 | 2 |
| naming | — | — | 6 | 6 |
| module_structure | — | — | 4 | 1 |

**Eighteen mutations were run across the three agents that mutated** (defaults 5, reliability 7,
idiomatic+errors 6). Twelve survived the six shipped tests; one further mutation was shown to
change no behaviour and was correctly not counted as a survivor. The overlap is large — the three
converged on the same two defects independently, which is why they are stated once below.

## 2. What was actually wrong

### 2.1 The config's support rule was unpinned — Blocker

Nothing compared `CandidateSelectionConfig::DEFAULT.support` with `DEFAULT_ALLELE_SUPPORT`. All
six tests read the constant directly; none read the config, **which is the value Milestone B/C
will actually consume**. Writing `support: MinAltReads::DEFAULT` instead — a one-token slip, both
names in scope in the same file, same type, same floor — left every test green while the share
dropped from 5 in 100 to 2 in 100.

**What that costs, measured by spec §3.3 rather than argued:** on the GIAB trio at 300× the two
bars keep 2,308 and 5,596 alternatives. The genotype prior divides its concentration by the
alternative count (spec §3.1), so the substitution more than halves the concentration a real
allele starts with, at every high-depth locus, with nothing in the output saying so. At tomato
depth it is invisible, because there the two shares are provably the same rule.

Found by all three mutating agents. Filed Blocker by one and Major by two, the difference being
whether "no consumer exists yet" mitigates it; the Blocker framing is the one this report takes,
because the gap is in the test set and the consumer is two steps away.

**Fixed** by `the_default_config_is_the_two_announced_constants_and_not_the_merges_rule`, whose
last assertion is `assert_ne!(config.support, MinAltReads::DEFAULT)` — the one that makes the
test about the *number* rather than the type. Re-running the mutation now fails that test and
only that test.

### 2.2 `new_const`'s lower bound was untested, and a negative share fails silently — Major

The two `#[should_panic]` tests passed `1.5` and `NaN`. **Both fail the upper comparison alone**,
so deleting `share >= 0.0 &&` survived all six tests, `cohort_merge`'s own 253, and everything
else in the crate.

**A negative share does not crash — it deletes half the rule.** `required_of` computes
`(share × reads).ceil()` and casts to `u32`; a negative product saturates to 0, so
`max(floor, 0)` is the floor at every depth and the high-depth half of the bar stops existing.
`required_of(300)` answers 2 where it should answer 15.

**Fixed three ways.** The range test is now written once, in a private `const fn
is_a_fraction_of_one`, and **both constructors call it** — so the fallible one and the `const` one
cannot come to disagree about what a legal share is, which is the same principle the module
already applies to `MinAltReads` itself. Its tests moved from the consumer's module to
`cohort_merge`'s own, beside the type they guard, and gained the three cases that were missing: a
negative share, both infinities, and both ends of the closed range as *accepted*. And
`the_const_share_refuses_exactly_what_the_fallible_one_refuses` walks ten values through both
constructors under `catch_unwind` and requires the same answer from each.

### 2.3 A cap of 0 or 1 is representable — Major, deferred to Checkpoint A

`max_candidate_alleles: u16` is `pub`, and `CandidateSelectionConfig { …, max_candidate_alleles:
0 }` compiles — one agent built exactly that `const` and `cargo check --lib` accepted it. At 0 or
1 the reference is the only survivor and every alternative becomes a truncation, which is refusal
under another name and is what spec §4.1 rules out.

**Not fixed here, because the fix changes a shape arch §2.1 declares** (`pub u16`) and would
ripple into the repeat-tract plan that inherits the same config. Raised at Checkpoint A. What
landed instead is the obligation, written into the field's own doc comment: `select_generic`
asserts a cap of at least 2 when it lands at step C2. That keeps the hazard from being
rediscovered rather than from existing.

### 2.4 Four measured figures were quoted under conditions they were not measured under — Minor

The naming agent checked **19 figures against the spec and found 18 correct**; the failures were
all about the *conditions*, not the values:

| what the doc said | what the source says |
|---|---|
| "4 loci in 53,935 differ" between 5 in 100 and 2 in 100 on tomato | spec §3.3 measured **0% against 2%**. No 5-against-2 tomato figure exists |
| "removes three quarters of the table" at 300× | true only against the count-only bar (2,308 kept of 10,793). Against the whole table it is 85%, against 2 in 100 it is 59% |
| floor 2→3 loses five, share→10% loses two | both are **30×** figures. At 300× the 10-in-100 bar loses four |
| the cap binds at 23 tomato loci and none of the trio's | §4.2's own header: "*with the bar at 2 reads or 2%*" — not the 5 in 100 this module ships |

All four rewritten to name their baseline and their depth. The "inert below about 40 compared
reads" claim was also replaced by the arithmetic that produces it — `ceil(0.05 × 40) = 2` is the
floor and 41 is the first count where the share asks for more — which is exact where "about 40"
was approximate, and is now what the crossover test pins.

### 2.5 The crossover test did not hold the claim it was written for — Minor

`the_allele_share_is_stricter_than_the_merges_only_once_depth_makes_it_so` was the one test
claimed to be discriminating. **Half of that survived scrutiny.** It does kill the most plausible
wrong implementation — writing `MinAltReadShare::DEFAULT` for the share — but its equality arm
stopped at 20 compared reads, where the floor still decides for every share up to 10 in 100. So
**doubling the share to 0.10 passed it**, one read short of the fixture that would have caught it.

Renamed `the_allele_share_binds_only_above_forty_compared_reads` and carried to 40 and 41, which
pins the share to more than 2/41 and at most 0.05. Re-running the doubled-share mutation now fails
three tests where it previously failed none of this one's assertions.

### 2.6 The docs build's deny-level lint was never being run — Minor

`Cargo.toml` sets `broken_intra_doc_links = "deny"`, and the module doc linked `` [`generic`] ``
— a sibling that does not arrive until Milestone B. `cargo doc --lib --no-deps` reported it.
**The step's verification set never ran `cargo doc` at all**, which is the finding behind the
finding. Fixed to plain code text until the file exists; `cargo doc` is now part of this step's
gate, and the count went from 24 unresolved links to 23 — the 23 are pre-existing on `main`, in
files this step does not touch.

### 2.7 Smaller things, all applied

- `#[allow(clippy::manual_range_contains)]` **suppressed a lint that never fires** — clippy already
  exempts `const` contexts. Measured both ways by one agent: deleting the attribute leaves clippy
  clean, and a temporary non-`const` twin makes it fire. Deleted; the explanation survives as the
  doc sentence it already was.
- **"The only callers are `const` declarations" was untrue in the same commit** — two of the three
  call sites were runtime calls from the step's own tests, and `pub const fn` is an ordinary
  function too. The claim is now scoped, and it names the release profile's `panic = "abort"` and
  points a runtime caller at `new`.
- **The equivalence paragraph omitted the infinities**, which are exactly why `new` carries the
  `is_finite()` call the `const` version drops. Now named, along with `-0.0`, which both accept.
- `PartialOrd` on `CandidateSelectionConfig` had no meaning and no caller. Dropped.
- `src/ng/mod.rs` still called `genotype_prior` "the first of the four"; `calling/mod.rs` described
  step 6 in the present tense where every neighbouring ng module doc says "so far". Both corrected.
- Test and binding renames: `the_cap_default_is_six_and_the_config_carries_it`, and the floor's
  coupling to `MinAltObs::DEFAULT` asserted rather than only its value.
- `required_of(301) == 16` added, because `0.05 × 300` is exactly 15 and the two original
  assertions could not see the rounding direction at all.

## 3. Disputed

**"`Option::unwrap` is const-stable since Rust 1.83, so `MinAltReadShare::new(0.05).unwrap()` may
make `new_const` unnecessary."** Refuted by compiling it:

```
error[E0015]: cannot call non-const associated function `MinAltReadShare::new` in constants
```

`new` is not a `const fn`, so what `unwrap` can do in a `const` never comes into play. `new_const`
is necessary.

## 4. Raised at Checkpoint A, not taken here

Four items, all of which would edit [`../../ng/arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md)
§2.1, which declares these shapes verbatim and which the repeat-tract plan inherits:

1. **`support` and `DEFAULT_ALLELE_SUPPORT` name observed evidence but hold a threshold.** In this
   crate "support" already means what a sample's reads showed (`AlleleSupportStats::num_obs`). The
   crate now carries two `MinAltReads` values differing only in the share, and the field's own doc
   says the point is that they not be confused. `min_allele_support` /
   `DEFAULT_MIN_ALLELE_SUPPORT` keeps the "minimum" at every use site.
2. **`max_candidate_alleles: u16` should be a newtype refusing anything below 2** — §2.3 above.
3. **`new_const` names the mechanism, not the panic.** `pub const fn` is already visible in the
   signature; the `assert!` is not. No precedent for the suffix anywhere in `src/`.
   `new_or_panic` was the suggestion.
4. **Six alleles: production's constant counts the reference, GATK's `--max-alternate-alleles`
   counts alternates.** Production's own comment equates them
   (`per_group_merger.rs:53-57`) and spec §4 inherits the equation. If GATK's six is six
   alternates, ng's cap is the tighter of the two by one allele and the spec's sentence needs a
   clause. **Unverified here** — the vendored `gatk/` tree is in neither checkout.

## 5. Validation after the fixes

All in the container, on the fixed tree:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo doc --lib --no-deps` — 23 unresolved intra-doc links, **all pre-existing on `main` and
  none in this step's files** (24 before the fix).
- `cargo test --lib allele_candidates` — 5 passed.
- `cargo test --lib const_share` — 5 passed.
- **Mutation re-check, four of the survivors re-applied to the fixed tree:** the config's support
  swap now fails 1 test; deleting the lower bound fails 2; making the upper bound exclusive fails
  1; doubling the share to 0.10 fails 3. Each was applied from a file backup and reverted by
  restoring that backup, with the restore verified by `diff -q`.

`cargo clippy --all-targets --all-features` is **red on this tree with 14 errors, none in `src/`**
— `benches/cohort_var_calling_perf.rs`, `benches/ng_joint_fit_perf.rs`,
`examples/ng_joint_contamination_harness.rs`, `examples/ng_joint_duplication_probe.rs` and
`examples/profile_posterior_engine.rs`, all inherited from `main` and none touched by this plan.
The gate this step is held to is therefore `--lib --tests`, stated rather than quietly widened.

## 6. What the fan-out cost, and one process note

Five agents, each with its own worktree, mutation-testing rather than reading. Every one of the
six findings that mattered came from a mutation; **not one came from reading the code**, and the
figure check that found four misattributed measurements came from reading the *spec* against the
diff rather than the diff on its own.

**One thing went wrong on the author's side and is recorded because it nearly cost the step.**
Probing whether `new` could be called in a `const`, the probe was appended to a file holding
uncommitted work and reverted with `git checkout <file>` — which discarded the step's own edit to
that file along with the probe. It was caught immediately by diffing the tree against the review
patch, and restored from that patch. The lesson is the one the plan-driven skill already states
for a different reason: **`git diff` the tree against a known-good copy, not `git diff --stat`,
and never revert a probe with a command that cannot tell the probe from the work.** Every
subsequent mutation in this pass used a file backup and a verified `diff -q` restore.
