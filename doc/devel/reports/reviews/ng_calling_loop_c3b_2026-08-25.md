# Code review — ng calling loop C3b: the final pass

**Scope:** the working-tree diff of step C3b of
[`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — 1,612 insertions across
`src/ng/calling/inference/summarise_condition.rs`, `src/ng/calling/mod.rs` and
`src/ng/calling/inference/mod.rs`, on top of `054ad2dc`.
**Date:** 2026-08-25. **Verdict: request changes** — **2 Blockers, 4 Majors, 5 Minors, 1 Nit**, and
**10 of 42 quantitative or mechanism claims wrong**. All applied; see
[the fix report](fixes_applied_2026-08-25_c3b.md).

**Three agents, each in its own git worktree**, each re-pointed at `054ad2dc` and handed the diff
as a patch — the working tree was uncommitted, so a `git checkout --detach` alone would have
reviewed the wrong code:

| agent | brief | outcome |
|---|---|---|
| reliability | tests and mutation testing | **12 mutations run, 6 survived, 3 killed, 3 changed no behaviour** — 2 Blockers, 4 Majors, 1 Minor |
| craft | naming, errors, idiomatic, smells, refactor-safety | 4 Minors, 1 Nit; four fixes compiled in the worktree before being proposed |
| numbers | re-derive every claim the diff makes about its own fixtures | **42 claims checked, 32 correct, 10 wrong** — 3 Major mechanisms, 4 Minor miscounts, 3 Nits |

**Three agents rather than the skill's one-per-category fan-out**, because the machine had 94 GiB
free and a three-agent fan-out on this repository has been measured at 44 GB. The craft agent
carried five checklists; its output keeps them apart.

---

## Blockers — both are tests that could not fail

**B1. No test pinned a per-sample genotype quality in a cohort of more than one.**
`summarise_condition.rs`, the final pass's `score_best_genotype` call. Every fixture asserting a
`genotype_quality` had exactly one called sample; the multi-sample fixtures asserted genotypes,
counts and copies but never a quality. **Measured: overwriting every called sample's quality with
the first called sample's leaves the suite green.** That is precisely the defect
[`calling_quality.md`](../../ng/spec/calling_quality.md) §3.1 says this pass exists to prevent —
the posterior row is one reused buffer, so a quality taken after the walk is the last sample's —
and the number becomes the per-sample `GQ` column, where a wrong-but-legal Phred panics nothing.
**Fixed** by `two_samples_take_their_confidence_from_their_own_posterior_rows`, which pins both
samples' qualities (11.553 and 17.697 Phred) against measured values.

**B2. The expectation's copy count was not pinned to the primary alternative.**
`copies_of_each_allele[usize::from(primary.get())]` is the only place the ninth artifact number is
joined to the allele the other eight were pooled for. Every fixture asserting
`genotype_expected_alternative_reads` was biallelic, where the primary alternative is necessarily
allele 1 — and the one three-allele fixture, built to defeat exactly that hard-coding in the
*choice* of alternative, asserted eight of the nine counts and not this one. **Measured: replacing
the index with the literal `1` leaves the suite green.** At a multi-allelic locus the allele-balance
test would then weigh reads of allele 2 against an expectation computed from copies of allele 1 — on
that fixture, 0 expected against 9 observed, a maximal apparent deficit, and only a deficit is
penalised. **Fixed** by two added assertions in the existing fixture.

## Majors

**M1. The reference row's forward and placed-left counts were never told apart.** Both fixtures
asserting them had those two counts equal (1 and 1), so swapping the two `+=` lines survived. The
two numbers feed *different* binomial tests — strand bias and read-position bias. **Fixed** by
changing the hand-computed fixture's reference row to 3 reads, 2 forward, 1 placed left, and
updating the five assertions and the doc comment that quote it.

**M2. No fixture proved a set-aside sample's reads are excluded from the *choice* of primary
alternative.** The existing fixture pins exclusion from the counts and the expectation, but it is
biallelic, so the set-aside sample's reads cannot change *which* allele is chosen. Measured:
deleting the skip in the pooling walk leaves the suite green. **Fixed** by
`the_primary_alternative_ignores_a_set_aside_samples_reads`, where allele 1 draws 10 reads across
the called samples and the set-aside sample's 20 reads of allele 2 would reverse the choice.

**M3. The repeat-tract test asserted no site quality**, though its name and doc comment both claim
the tract carries one — and the site quality is the *only* number that arm has to get right, since
it skips the artifact summary by design. **Fixed** by comparing against an independently prepared
fold, as the SNP/indel test already did.

**M4. Four of the ten values handed to `LocusInference::new` were not pinned as travelling onto the
record** — the region, `converged`, `passes` and the provenance. Measured: replacing all four with
constants leaves the suite green. `LocusInference`'s own tests cover the fields but construct the
type directly, so they say nothing about this pass's wiring. **Fixed** by
`the_pass_carries_the_region_the_outcome_and_the_provenance_onto_the_record`.

**M5–M7 (numbers agent), all mechanism claims, all wrong as written:**

- *"the M-step refuses the locus loudly on the scratch's `NaN` sentinel"* — it is the **E-step**, in
  the loop's prior-free first pass, before any M-step runs. The conclusion the sentence supports
  ("nothing wrong reaches a run today") stands; the named guard was the wrong one, in both the doc
  comment and the implementation report.
- *"each sample keeps a fifth of its probability on the homozygous reference"* — the fixture's
  converged posterior is `[0.34381, 0.65068, 0.00551]`: **a third**. This is the sentence explaining
  why the copies are two thirds of a copy apart, and at a fifth the arithmetic does not close.
- *"nothing in ng can score a tract yet, since the repeat-tract read-likelihood row is
  `unimplemented!()`"* — **it is implemented**, in `likelihood/ssr.rs` over the shipped
  `StutterSubstitutionEmission`; the only `unimplemented!()` under `src/ng/calling/` is a
  `#[cfg(test)]` oracle's `censored_emission`. The claim was inherited from the plan's own blocker
  note, which this commit corrects: what still blocks a tract end to end is the repeat-tract
  **candidate** path, which is unwritten.

## Minors

- **Naming.** `pooled_primary_alternative` was the module's only noun-phrase free function, and it
  resisted a verb because it did four jobs — including carrying the pass's only two input
  validations. Renamed `pool_reads_and_pick_primary_alternative`.
- **Errors.** Three `expect()` calls without the repo's `// PANIC-FREE:` marker, one of whose
  invariants is already written out at `genotype_table.rs:562`. Markers added, naming the enforcing
  constant rather than restating the claim.
- **Idiomatic.** `add_called_sample` took `ploidy: f64` where `Ploidy` exists and refuses zero — and
  the one call site *unwrapped* a `Ploidy` to hand over the primitive. A zero would have put `inf`
  or `NaN` into the expected-read count. Now takes the newtype.
- **Smells.** `generic_samples` and `primary_alternative` were co-dependent `Option`s re-paired at
  the use site; the combination that cannot arise became a silent no-op. Now one value.
- **The ceiling check had no boundary test**, and the boundary is a value
  `score_uncorrected_site_quality` actually returns. Added.
- **Four miscounts of the author's own code** (numbers agent): a biallelic locus called
  two-alternative, `PooledArtifactCounts` described as pooling nine numbers where it pools eight, a
  "tenth" scratch buffer that is the fourteenth, and an `allow` reason claiming `call_locus` shares
  five of the arguments where it shares four. All corrected.

## Nit

- `LocusEvidence::Generic { per_sample, .. }` hid exactly one field; naming it makes a third field
  on the variant a compile error here. Applied.

## Three judgement calls, all upheld

- **The artifact summary as an `Option`** — right, and for the reason the doc comment gives.
- **`site_quality` private with a `pub(crate)` reader** — keep. The reader carries
  `expect(dead_code)` rather than `allow`, so when the correction stage starts calling it the
  compiler *errors* on the unfulfilled expectation and forces the attribute out.
- **Nine arguments** — keep. All nine types are distinct, so a transposition is a compile error;
  the four pass-through arguments are `LocusInference::new`'s own, and bundling them would move
  that constructor's checks away from the code producing the values.

## Out of scope observations

- **Release-profile `+=` on read counts.** `overflow-checks` is off in `[profile.release]`, so a
  release overflow would wrap silently. The accumulators are `u64` over `u32` addends, so this needs
  on the order of 4×10⁹ observations at one locus; judged unreachable and not changed.
- **The site quality's fold is quadratic in cohort size and runs unconditionally per locus.** A
  performance question for the bakeoffs plan, not a correctness one
  ([`calling_quality.md`](../../ng/spec/calling_quality.md) §13's Q3 carries it).
- **`ArtifactTestCounts` derives `Copy` at about 72 bytes** — the previous commit's, unchanged here.

## Verification

Run in the container from the main checkout, after the fixes:

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings` —
  exit 0.
- `cargo test --lib` — `4672 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `626 passed; 0 failed; 3 ignored`.
- The seven release-held checks downgraded together: `618 passed; 8 failed`, every check reached.
