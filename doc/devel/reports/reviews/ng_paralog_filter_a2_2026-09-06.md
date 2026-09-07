# Code Review: ng_paralog_filter_a2
**Date:** 2026-09-06
**Reviewer:** rust-code-review skill (orchestrator) over three category sub-agents, each in its own worktree
**Scope:** commit `a8af3866` on `ng-paralog-filter` — the per-locus paralogy score copied into `src/ng/paralog/`, and the differential that asserts it computes what production's does
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** one commit's diff.
- **Reviewed against:** `a8af3866`, branch `ng-paralog-filter`, parent `bffc74bb`.
- **In-scope files:** [production_parity.rs](../../../../src/ng/paralog/production_parity.rs) (ng's own, the review's centre of gravity), [copy_fidelity.rs](../../../../src/ng/paralog/copy_fidelity.rs) (ng's own, gains the `Repoint` mechanism), [mod.rs](../../../../src/ng/paralog/mod.rs), [locus_score.rs](../../../../src/ng/paralog/locus_score.rs) (a copy; its `use` lines, its placement and ng's appended header note only), [src/ng/mod.rs](../../../../src/ng/mod.rs).
- **Deliberately out of scope:** the *logic* of the three copied files. They are production's byte for byte and may not be edited here; sub-agents routed defects in them to `Cross-category` as informational, and five such observations came back (dead code in production's own scorer — see §7).
- **Categories dispatched:** `reliability` (the differential is this step's only real oracle; if it is weak the step has none), `module_structure` + `errors` (the copy now imports production by path and the guard sanctions it — the one structural question the commit raises), `naming` + `smells` + `idiomatic` (the guard and the differential are the only new Rust, and their prose is what a human reads to decide whether to trust the port).
- **Categories skipped, with reason:** `defaults` — every `DEFAULT_*` in scope is production's, frozen. `unsafe_concurrency` — no `unsafe`, no threads. `tooling` — no `Cargo.toml` change. `extras` — no parser, no untrusted input; its "diff matches stated intent" item ran as the orchestrator's step 8a.

### 2. Verdict

**Approve-with-changes.** The copy is production's and the differential agrees with it — both confirmed independently, the second by a reviewer re-deriving the splitmix64 stream outside Rust and matching all six pinned counts exactly. What the review found is that **the differential's reach was narrower than its own documentation claimed**, in a way four surviving mutations proved, and that the guard's new `Repoint` mechanism reached one line short of where it needed to. All are fixed; see the fix-application report.

### 3. Execution status

| command | result |
|---|---|
| `cargo test --all-features --lib "ng::paralog::"` | `ok. 54 passed; 0 failed` |
| `cargo test --all-features --lib --bins --tests` | 6,329 lib tests passed; **one integration test failed**, pre-existing |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 3 errors, `needless_lifetimes` in `cohort_merge/`, pre-existing |
| `cargo fmt --check` | dirty on 9 files, none this commit's, pre-existing |

`cargo test --all-targets` still does not compile on `main`; `cargo doc --no-deps` is red there on 11 unresolved links, though one reviewer ran it scoped and confirmed `paralog` emits none. Findings labelled "Needs verification": **0** — every finding was reproduced by its agent in its own worktree, and each fix was reproduced again in the working tree before and after.

### 4. Open questions and assumptions

1. **How does step A3's copy fit `GuardedCopy`?** A3 puts *two* originals into one ng file — `src/paralog/prior.rs` plus `calibrate.rs:64-118`, a line range in a different module. The type holds one source and no end marker. Affects **Mi1**; recorded for A3 rather than pre-built here.
2. **Nothing constrains what a `Repoint` may say**, only where it applies. Partly closed (a repoint must turn a `crate::paralog` path into a `crate::ng::paralog` one — asserted); the general form is A3's, when files the differential does not cover arrive.

### 5. Top 3 priorities

1. **M1** — four input classes the differential never drew, each with a mutation to ng's copy that survived it. Two of the four violate a spec trap outright.
2. **M2** — ng's copy documents production's type in two rustdoc links, and one is inside the header the guard forbids a copy to touch, so it could not be repointed at all.
3. **M3** — the test whose stated job is to notice a narrowed generator asserted only `count > 0`; a 99-in-100 narrowing passes it.

### 6. Findings

#### Major

**M1: [production_parity.rs:96-140](../../../../src/ng/paralog/production_parity.rs#L96) — four input classes never drawn, four surviving mutations**
**Category:** reliability · **Confidence:** High
Copy numbers were drawn from `[0, 9)` and both the σ₀ slice and the per-pass tables were always built from the locus itself, so four classes never occurred. Each left a branch of the scorer un-compared, and each of these mutations to ng's copy passed all five differential tests when `production_parity` was run alone: narrowing `clamp(0.0, cmax)` to `min(cmax)` (a negative copy number); deleting the `precompute.cohort_size != cohort_size` clause; weakening the σ₀ length check from `!=` to `<` (a slice *longer* than the cohort); deleting the `configs.is_empty()` clause. Decoded, the first left ng scoring one locus at `81.374` against production's `39.028` — more than twice as high — and the fourth returned `0x7ff8000000000000`, a **`NaN` ratio from a scored locus**, which spec §6 trap 4 forbids from ever leaving the scorer. The second and third are spec §6 trap 3 verbatim. **For all four the text guard caught what the differential missed** — the reverse of the ordering `production_parity.rs`'s own header argued for, and the ordering that matters once `coverage_model.rs` is released as already scheduled.

**M2: [locus_score.rs:33](../../../../src/ng/paralog/locus_score.rs#L33), [:65](../../../../src/ng/paralog/locus_score.rs#L65) — ng's `pub` module documented production's type, and one link was unreachable by the mechanism**
**Category:** module_structure · **Confidence:** High
Both lines are rustdoc link definitions reading `` [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel ``. ng re-exports its own, so `cargo doc` on ng's public module sent a reader into the frozen tree for the type ng's scorer actually consumes. `copy_fidelity.rs` claimed the module had exactly one such line; there were three. Worse, the first sits inside production's 33-line module header, which is compared as a byte-verbatim prefix while repoints were applied only *past* the header — **so that line could not be repointed at all as the mechanism was drawn.** A `use crate::` sweep never finds a rustdoc link, which is part of why it survived.

**M3: [production_parity.rs:418-450](../../../../src/ng/paralog/production_parity.rs#L418) — the guard against a narrowed generator asserted only `count > 0`**
**Category:** reliability · **Confidence:** High
Its own doc comment says it exists because "a generator that quietly stopped drawing zero-read samples … would leave the differential passing over a narrower input class than its own documentation says — and nothing else would notice". Replaying the stream with one rate narrowed at a time: the malformed-record rate cut from 1 in 10 to 1 in 1,000 takes those draws from 1,256 to **17** in 12,600 and the test still passes; likewise zero-read (to 61) and absent (to 19); and the very σ₀ bug the author had already found and fixed (4 in 7 degenerate) also passes.

**M4: [production_parity.rs:412-436](../../../../src/ng/paralog/production_parity.rs#L412) — the test that certifies the stream rebuilt it from two hand-typed literals**
**Category:** naming/smells · **Confidence:** High
It opened `Splitmix64(0x5eed_0000_0000_0000 ^ 63u64)` and drew `draw_locus(&mut rng, 63)`; the sweep opens `… ^ cohort_size` over `COHORT_SIZES = [1, 2, 10, 63]`. The two agreed only because `63` and the seed base were typed twice. Change either and this test passes while describing a stream the differential never draws.

**M5: [copy_fidelity.rs:188-195](../../../../src/ng/paralog/copy_fidelity.rs#L188) — the guard's structural-failure branch had no test**
**Category:** errors · **Confidence:** High · **Convergent** with reliability, which found the same gap from the mutation side (deleting the "exactly once" loop, and loosening the exact match to a prefix, both survived).
Both synthetic fixtures passed `repoints: &[]`, so the whole `Repoint` path ran only in its accepting direction, on the one real pair. Two of its three failure modes are caught incidentally by the content compare; the third — a sanctioned line coming to occur twice, rewritten on both sides — is caught by nothing else. The file's own doc says the synthetic tests exist precisely to prevent this.

**M6: [production_parity.rs:383-386](../../../../src/ng/paralog/production_parity.rs#L383) — a doc comment claimed a spec trap was tested that the scorer cannot even reach**
**Category:** naming · **Confidence:** High
The all-absent test's header said it pins spec §6 trap 4. Trap 4 is about `NaN`; `score_locus_for_paralogy` never returns one — `ParalogScore::neutral()` sets every float to `0.0`, and the same file says so twice elsewhere, once in warning terms. As written, a reader concludes a Blocker-class spec trap has a test. It does not: the `NaN` sentinel belongs to the run wiring, which arrives at step B1.

#### Minor

**Mi1: [copy_fidelity.rs:85-137](../../../../src/ng/paralog/copy_fidelity.rs#L85) — `GuardedCopy` pairs one ng file with one whole production file, and A3's next copy is neither.** The plan puts `prior.rs` *and* `calibrate.rs:64-118` into one ng file; the type holds a single source and a `CopyBegins` whose variants both mean "to end of file". `include_str!` and three failure messages also hard-code `src/paralog/`. Deciding this now costs a paragraph. **Category:** module_structure.

**Mi2: [copy_fidelity.rs:404,427-436](../../../../src/ng/paralog/copy_fidelity.rs#L404) — `GUARDED` restated `guarded_copies()`, and the assertion that they agree compared only lengths.** A name swapped in one list and not the other keeps the lengths equal. The surrounding test argues that a hand-maintained list needs an outside authority, then adds a second hand-maintained list beside it. **Category:** module_structure.

**Mi3: [copy_fidelity.rs:243-268](../../../../src/ng/paralog/copy_fidelity.rs#L243) — the divergence message names two causes and there are three.** A deleted `Repoint` declaration lands in that message, and its advice — "whichever side moved, revert it" — would put `crate::paralog::` back into ng's copy, the outcome the whole mechanism exists to prevent. Folding the structural `Err` into the same `Option` **is** right at the sink; the defect is one level up, in the wording. **Category:** errors.

**Mi4: [copy_fidelity.rs:191](../../../../src/ng/paralog/copy_fidelity.rs#L191) — eighteen literal spaces in a failure message, and it does not say the guard is off.** A `\` continuation that was lost; `cargo fmt` cannot see inside a string. When it fires, *nothing is checking ng's copy*, and the reader was not told. **Category:** errors, smells (convergent).

**Mi5: [copy_fidelity.rs:145-158](../../../../src/ng/paralog/copy_fidelity.rs#L145) — `AtTheFirstItemDoc` can extract nothing, and nothing against nothing passes.** If production's file loses every `/// ` line — an item doc rewritten as `/**`, or reordered behind a `#[doc = …]` — both sides extract `""` and every comparison holds. `model_params.rs` is guarded in exactly this mode. A guard failing *open* is the worst of the three outcomes. **Category:** errors.

**Mi6: [production_parity.rs:209-250](../../../../src/ng/paralog/production_parity.rs#L209) — `ParalogScore` read field by field, so a sixth value would go uncompared.** The header claims all five are asserted; nothing keeps that true. This plan's own B1 already requires exhaustive destructuring for the same reason. **Category:** errors.

**Mi7: [src/ng/mod.rs:19-21,31-36](../../../../src/ng/mod.rs#L19) — ng's register of what has landed and of which tests read production is stale.** It said the filter is "so far the per-sample coverage model" after A2 landed the score, and "a test may read production as an oracle, and two do" when a grep finds at least eight. A count wrong by six is not consulted. **Category:** module_structure.

**Mi8: [production_parity.rs:117-118](../../../../src/ng/paralog/production_parity.rs#L117) — "scores nothing 57 times in 100" omits the absent draw.** With σ₀ degenerate 4 in 7 and absent 1 in 8, the share is `1 − (7/8)(3/7) = 5 in 8`. **Category:** smells.

**Mi9: [production_parity.rs:38-39](../../../../src/ng/paralog/production_parity.rs#L38) — one of the two cited precedents for the hand-rolled generator is a xorshift, not a splitmix64.** **Category:** smells.

**Mi10: [copy_fidelity.rs:34-43](../../../../src/ng/paralog/copy_fidelity.rs#L34) — "Each such line is declared as a `Repoint`" was not true**, and the mechanism was described as preventing a reach "at run time" when a rustdoc link reaches back at *read* time. **Category:** smells.

**Mi11: [production_parity.rs:433,331,350-352](../../../../src/ng/paralog/production_parity.rs#L433) — `4.0` written out where `DEFAULT_MAX_RELATIVE_COPY_NUMBER` is the value, and `200` where `LOCI_PER_COHORT_SIZE` is.** **Category:** naming.

**Mi12: the implementation report's quoted mutation values did not reproduce.** The pair `26.61449637978953` / `26.614488944807874` came from a run *before* the σ₀ draw rate was fixed, so the stream had changed under it; on the committed tree the same mutation fails at `cohort size 1, locus 0` with `85.12238217765987` against `85.12233325249034`. Its "7 parts in 10⁹" was also wrong on its own quoted values — about 3 parts in 10⁷. The claim the numbers support holds; only the numbers were stale. **Category:** the orchestrator's step 8a, found independently by reliability.

#### Nits

- `guarded_copies()` allocates a `Vec` where a `const` slice would do.
- `Repoint`'s doc calls it "one line a copy is allowed to change, because leaving it would point ng's code back into production", but the guard constrains only *where* a substitution happens, never what it says: a repoint rewriting a constant would pass and read as sanctioned.
- `next_below`'s modulo is biased for bounds that do not divide `2^64` — fine for drawing test shapes, worth not claiming otherwise; and its non-zero precondition is a doc comment rather than a type.
- The two saturated `200 of 200` pins in `LOCI_THAT_SHOULD_SCORE` carry little signal on their own.
- The differential runs at one parameter setting only, and at no cohort size near spec §4's three thousand.
- `mod.rs` described `model_params.rs` under the rule that exempts `mod.rs` from guarding, when `model_params.rs` is in fact guarded.

### 7. Out of scope observations

Five equivalent mutants pointing at **dead code inside production's frozen scorer**, reported informationally: the `total_reads > 0` clause in the hom-alt guard never discriminates (`0/0` is already `NaN`, so the VAF comparison is false at any `homalt_min_depth`); `ms.sort_unstable()` never reorders; the `usable.is_empty()` early return is a performance guard only, since both marginals give exactly `0.0` on an empty set; and `LogSumExp`'s `>` against `>=` is exactly equivalent. These belong to `src/paralog/`'s backlog. The `GridSpec` public-fields-versus-validated-constructor trade filed in A1's review still stands.

`cargo fmt` in the container rewrites nine files this branch never touched. One reviewer could not establish the cause and rated a rustfmt version difference as likely as drift; either way the tree is not `fmt --check` clean on `main`.

### 8. Missing tests to add now

All applied — see the fix-application report.

- `a_negative_copy_number_is_winsorised_at_zero_on_both_sides`, `tables_built_for_another_cohort_are_neutral_on_both_sides`, `an_empty_carrier_set_is_neutral_on_both_sides`, and a σ₀ slice *longer* than the cohort folded into `a_mismatched_sigma_slice_is_neutral_on_both_sides` — M1's four classes, each verified to kill its mutation with `production_parity` run alone.
- `a_repoint_is_applied_where_declared_and_nowhere_else` — M5, driving every rejection of the `Repoint` path on synthetic strings including the structural one.
- `a_guard_that_can_find_no_copied_content_fails_rather_than_passing` — Mi5.
- The shape counts pinned exactly, and the whole stream derived from `COHORT_SIZES` and one `stream_for` — M3 and M4.

### 9. What's good

- **The commit does not trust its own tests to prove parity, and says why** — production's transcribed tests pass on both trees whatever either computes, which is the argument that makes the differential worth having.
- **Comparison by bit pattern rather than tolerance**, with the reason stated: a port agreeing to twelve decimal places is a port that has diverged.
- **`N = 1` was put in the sweep before anything asked for it**, and it is the case production has never run.
- **The pinned per-size scored counts** turn a later narrowing of the generator into a failing assertion rather than a silent loss of reach — a reviewer re-derived all six counts outside Rust and matched them exactly.
- **`score_through_both` returns two distinct types**, so the two verdicts cannot be transposed the way a same-typed pair could.

### 10. Commands to re-verify

```
scripts/dev.sh cargo test --all-features --lib "ng::paralog::"
scripts/dev.sh cargo test --all-features --lib --bins --tests
scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
scripts/dev.sh cargo fmt --check
```

The last two are expected to fail on `main`'s three lints and nine files; the gate is that neither set grows.

Per-category audit trail: `tmp/review_2026-09-06_ng_paralog_a2/{reliability,module_structure_errors,naming_smells_idiomatic}.md`.
