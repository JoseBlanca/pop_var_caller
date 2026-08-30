# ng calling loop — E3a: ng calls genotypes on the SNP/indel path

**Step:** E3a of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the generic half of the
integration fixture. **E3 split into E3a and E3b in the course of this work**; §5 says why.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §1, §5.0, §9;
[`arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md) §5.1.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

**ng calls genotypes from real evidence on the SNP/indel path.** A new test binary,
[`tests/ng_calling_loop_calls_genotypes.rs`](../../../../tests/ng_calling_loop_calls_genotypes.rs),
hands the merge one cohort locus and takes `LocusInference` out:

```text
per-sample observations → ClosedLocus → CohortObservation::over → select_generic
                       → shape_generic_locus → call_locus → LocusInference
```

No source behaviour changed; the only `src/` edit is one refusal message, which now names E3b.
**This step is a test**, and that is the milestone: every earlier test of this loop handed it a
likelihood table, an evidence view or a candidate set built for the occasion.

## 2. What it supplies and what it runs

**Supplied**, each for a reason:

- **the per-sample observations** — what each sample's reads showed, as the SNP/indel locus
  generator emits them. Turning aligned reads into observations is step 5's and outside this
  plan's Scope;
- **the `ClosedLocus`** the merge is handed. `LocusCloser`'s chaining walk — which groups
  overlapping observations into loci and judges each — is a different subsystem and is not run;
- the run's frozen parameters and the loop's configuration, which are a run's inputs rather than
  a locus's.

**Run rather than supplied:** the merge's allele unification and read attribution
(`CohortObservation::over`), candidate selection (`select_generic`), the input edge
(`shape_generic_locus`), and the loop (`call_locus`).

**What the fixture reproduces rather than assumes: the walk's keep rule.** `merge` asks
`MinAltReads::DEFAULT` of each sample separately — at least 2 non-reference reads and at least 2%
of that sample's own compared reads — and refuses a fixture no sample reaches it with, because
such a locus the real walk would call too quiet and never build. **That check caught a fixture on
its first run** (§4).

**It lives in `tests/` rather than in the library suite** so that the seams this path names are
`pub`. That is the twenty-odd items it imports, not every seam in the middle — `evidence_shaping`
has a `pub(crate)` accessor on this path and the test compiles regardless.

## 3. The ten tests, and what each pins

| test | what it pins |
|---|---|
| three samples at a SNP | `1/1`, `0/1`, `0/0` from 20, 10+10 and 20 reads; **1 pass** |
| the cohort's expected copies | **3.0000019 reference against 2.9999981 alternative** out of six, to 1e-6 |
| one sample | the cohort end of the range: no panel, so the seed alone is the concentration |
| a sample that covered nothing | **the merge lists 2 samples and the loop 3** — the join no type enforces |
| the merge builds it, selection calls it over the reference alone | the 2% bar against the 10% bar |
| two alternatives | three alleles, and **`C` is merge index 3 and candidate id 2** |
| three reads a sample | the same three calls, and **4 passes** — the only fixture where the loop iterates |
| the weakest warrant of two libraries | `Defaulted` beside `FittedHere`, **asserted in both orders** |
| all libraries fitted | the same locus claims `FittedHere`, so the field is not hard-wired |
| the locus's own region and allele table | shape, not arithmetic |

## 4. What the reviews found, and what changed

Two agents in worktrees detached at `8c6664a3` with the diff applied: one hunting tests that
cannot fail, one checking claims and design conformance. **Between them: 2 Blockers, 8 Majors,
and 24 of 64 claims wrong.** Every finding was applied. Every counted figure in the first draft
was correct; every failure was an explanation or a fixture.

**Blocker 1 — the join the fixture existed to test was the identity.** `merge` pushed a
`SampleMembers` entry for every sample, including one with no observations — a shape the real
merge cannot produce, since `SampleMembers`' own contract is *"a sample with nothing here has no
`SampleMembers` at all"*. So the merge's covering-sample list and the run's sample list were the
same list in all seven tests, and **an input edge that used the merge's index as the run's passed
every one of them**. `merge` now skips a sample with no observations, and the silent-sample test
asserts the merge's list before the loop runs. The mutation now fails it.

**Blocker 2 — allele ids could not be told apart from merge positions.** Nothing was ever dropped,
so the remapping was the identity, and the three-allele test resolved each id by looking its bases
up in the output — invariant under any permutation. **Deleting the remapping entirely left it
green.** The fixture now carries a fourth sample whose `G` clears the merge's 2% bar and fails
selection's 10%, so the merge interns `[A, T, G, C]` and the candidates are `[A, T, C]`; the ids
are asserted literally. Both the dropped remapping and a reversed admission order now fail it.

**The keep-rule check caught a third fixture that could not arise.** A test asserting that one
stray read makes no candidate used 1 read of 21 — which fails the *merge's* bar as well, so the
real walk would have discarded the locus. It is now 3 reads of 100, which clears 2% and fails 10%,
and that is the only gap in which such a locus exists.

**The loop never iterated.** Every fixture was at twenty reads a sample, where the emission decides
every genotype and the frequency loop settles in one pass — so **a pass cap of 1 and a convergence
threshold of 0.1 both left the suite green**, and two `assert!(converged)` could not fail. There
is now a fixture at **three reads a sample**, which is where this project's tomato cohort sits, and
it takes **4 passes**; the pass count is asserted exactly in three tests.

**The weakest-warrant fold was tested nowhere in the repository.** One read group everywhere meant
`weaker_of` never combined two provenances: replacing it with *last one wins* passed this binary
and all 4,815 library tests. There are now two libraries with different provenances, asserted in
both orders.

**Six mutations, all previously surviving, now fail a test in this binary**: the input edge's
index, the dropped remapping, a reversed admission order, the pass cap, the convergence threshold,
and the warrant fold.

## 5. The silent sample, and a mechanism the first draft got wrong

**A sample with no reads at the locus comes out `0/1`.** Its reads score every genotype alike, so
what decides it is its own prior: the seed plus what the **other** samples showed, which is what
the leave-one-out subtraction leaves. Its two neighbours are called `1/1` and `0/0`, so they
contribute exactly `[2.0, 2.0]` copies, and against a seed of `[1.0, 0.001]` its concentration is
`[3.0000019, 2.000998]`. **The heterozygote wins by two parts in ten thousand** — posterior
`[0.39985, 0.40005, 0.20009]`, about 2 Phred of genotype quality.

**Two things the first draft said about this were wrong.** It quoted the cohort's total, 3.20
against 2.80, as what the sample was scored against; that total *includes* the sample's own
posterior contribution, which the leave-one-out subtraction removes, so it was scored against an
even 2:2 split. And it attributed the result to the owner's open question recorded against C1 — a
silent sample's flat-pass vote inflating the cohort estimate **other** samples see. That is a
different mechanism: here the subtraction removes its own vote from its own prior, so what decides
it really is its two neighbours.

**It also quoted a sentence from `spec/calling_em_loop.md` §7 that is not in that document.** The
real sentence is §5.0's: *"A sample with no reads at the locus already contributes zero for every
genotype and is decided by the prior alone (§7)."* The misquotation was inherited from an earlier
step's record and restated here in a new source file; it is corrected in the test's doc comment.

**What is worth the owner's attention is smaller than the first draft claimed and still real**: a
silent sample is called, not set aside, and at a polymorphic locus what it is called is whatever
its neighbours make most probable — at 2 Phred of confidence.

## 6. Why E3 split

**E3a is a test; E3b is implementation and one open question.** Five places in
`inference/summarise_condition.rs` take the generic path's per-sample evidence unconditionally and
each needs a tract branch, with `call_locus` refusing a tract in front of them — so E3b changes
library code where E3a changed none.

**The open question is where a tract's two run-level prior inputs live.** `fill_ssr_seed` needs the
cohort's repeat gene diversity and a decay per repeat. Both are checked types in `ng/types.rs`;
**neither is emitted by anything** — `src/ng/parameter_estimation/` produces neither, and the
pre-pass cohort gather that would is unbuilt. What is undesigned is only the **calling-side
carrier**: the quantities themselves have a designed source (`arch/calling_priors.md` §5 names the
pre-pass, and `spec/parameter_prepass_ssr.md` gives the diversity a **Home:**).

**The first draft argued this badly and the correction matters**, because it makes the
recommendation stronger rather than weaker: it said `FrozenParameters` "carries the SNP/indel
spectrum seed and nothing else". It does not — E2 put two run-level STR parameters on it,
`ssr_slippage_fits` and the substitution-rate map. So putting the tract's seed inputs there follows
a precedent rather than setting one. The plan's E3b entry carries the recommendation, the
alternative and what the alternative costs.

## 7. Validation

Run with `./scripts/dev.sh` from the `ng-calling-loop` worktree.

| gate | result |
|---|---|
| `cargo test --test ng_calling_loop_calls_genotypes` | **10 passed / 0 failed** |
| `cargo test --lib` | 4,815 passed / 0 failed / 14 ignored — unchanged |
| `cargo test --release --lib ng::calling --all-features` | 754 passed — unchanged |
| `cargo test --test ng_calling_loop_allocation --features dhat-heap` | 1 passed |
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo doc --no-deps --lib` | 28 unresolved links, **exit 101** — unchanged, none added |

`cargo doc` exits 101 because the crate denies broken intra-doc links, so the pre-existing 28 are
hard errors; the row is here because the check that matters is *this step added none*.

## 8. What the fixture still does not cover

- **The high end of the cohort axis.** One, two, three, four and five samples are here; nothing
  near the thousands this caller commits to, and nothing here says what changes there.
- **The high end of the depth axis.** Three reads and twenty; nothing at several hundred.
- **Partial reads.** Every observation is `ReadWitness::Complete`, so the partial arm of
  `GenericSampleEvidence` is empty in every fixture.
- **Indels.** Every allele is one base, at a one-base locus, in a file whose subject is the
  SNP/indel path.
- **Contamination.** Every run is `FrozenParameters::uncontaminated`, so the loop's per-pass
  reassembly branch is not exercised here — E2a's own fixtures are where it is.

## 9. Banked for the owner

1. **The open question E3b carries** — where a tract's two run-level prior inputs live. §6, and the
   plan's E3b entry, which states the recommendation and what the alternative costs.
2. **A sample with no reads is called `0/1` at a polymorphic locus**, at about 2 Phred. §5. Whether
   that is right is close to, and not the same as, the question already open against C1.
