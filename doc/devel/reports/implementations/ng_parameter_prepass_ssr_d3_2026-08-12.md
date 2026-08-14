# ng step 4, the STR path — D3: a genotype set wider than three

*Implementation report, 2026-08-12. Step D3 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — one agent, 6 mutations, 1 survivor that named a real hole. Design
authority: [`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §3,
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.2.*

## What the step is

The shared climb over the genotype frequencies declared its answer as `SmallVec<[f64; 3]>` — three
being how many genotypes a diploid has on the SNP/indel path, where a genotype is a count of
alternative copies. On the STR path a genotype is an unordered tuple of allele *lengths*, and a
diploid stratum has between 66 and 91 of them. The declared type widens to `Vec<f64>`.

**This is the first of the two changes the sibling plan anticipated `fitting/` would need**, and the
step where the seam is tested rather than argued about.

## Recorded deviations from the plan

**The plan calls it "one line"; it is fifteen, across four files.** `fit_mixture_weights`'s return
type is one of them, but the same width is restated by `MixtureWeightsFit::genotype_frequencies`,
by `ScanResult::genotype_frequencies`, by `fit_read_group_error_rates`' parameter and by four
signatures in the generic path's coupled fit that pass the map through. Nothing about the change is
harder than the plan says — every one is mechanical and the generic path's tests are the proof —
but "one line" understates where the width was written down.

**The climb itself is unchanged**, as the plan requires: the body of `climb_with_cap` differs by
exactly three lines, all of them a type and its constructor. The expectation step, the maximization
step, the monotone-ascent assertion, the stopping rule and the final scoring pass are byte-identical.

## What the review changed

**Major — nothing pinned that the climb reads its start at all.** The new test's own doc comment
claimed the skewed start tested start-independence. It did not, and neither did the pre-existing
`every_interior_start_reaches_the_same_summit`: on a concave surface *"arrived from the skewed
start"* and *"never looked at the skewed start"* produce identical output, so a climb that
validated its argument and then discarded it left the whole module green. Measured — the reviewer
replaced `start.to_vec()` with the uniform point and got 30 passed, 0 failed. The test now compares
where **one** pass lands from each of two starts, before either has had time to converge. I
reproduced the mutant against the fix: it now fails.

**Two wrong claims of mine.** The comment cited `spec/parameter_prepass_ssr.md` §4.2 for the figure
of 91 genotypes; that section gives the `A(A+1)/2` formula and works it at **nine** lengths, for 45,
and never mentions thirteen or 91. The 91 is `arch/parameter_prepass_ssr.md` §3's, now cited there.
And "up to 91" invited the reader to wonder whether a narrow stratum fits inline: it cannot, because
`allele_support` clips only at the low end, so the narrowest stratum the copy floors admit still
spans eleven lengths and 66 genotypes — **22 times the inline three**. Every STR fit spills.

**One rename, because my change made an existing name looser.**
`a_rung_loop_refills_one_buffer_and_allocates_nothing_per_rung` counts no allocations and never did;
it asserts numerical recovery across rungs. With three more `Vec`s inside the climb it read as a
guarantee it was not making, so it is now
`a_rung_loop_refills_one_likelihood_buffer_rather_than_rebuilding_it`.

**The blast radius was checked and is clean.** `generic/runs.rs`'s two `SmallVec<[f64; 3]>` fields
are the runs model's grid of starting guesses, not genotype frequencies, and are correctly
untouched — they are the only ones left in the crate. `SmallVec<[GenotypeFrequency; 5]>` in
`SampleRates` is a genotype-frequency vector but is the SNP/indel dosage representation the STR path
does not travel through; `coupled_fit.rs:689` is the line a later step would need if it ever does.

**The allocation cost was priced rather than asserted.** At diploid width the climb makes ten
transcendental calls per cell per pass, so 5,830 per pass over the generic path's 583-cell table,
against three 24-byte allocations per call — and `fit_mixture_weights` already built a `Vec` per
call before this change. The figures are op counts read off the code, not a benchmark.

## ⛦ Raised for Checkpoint D, not fixed here

Four lines of the architecture documents now describe the pre-change state, including
`arch/parameter_prepass_ssr.md` §7's *"a change this unit forces on the shared module, and it is not
optional"* — which is now done rather than pending. Editing the design documents is the owner's
call.

## Tests

One new, in the shared module. A 45-genotype table — nine allele lengths' unordered pairs, spec
§4.2's own worked example — filled with the **exact expected counts** under a known truth rather
than drawn ones, so the maximum-likelihood point *is* the truth and the tolerance is decided rather
than chosen. Recovered to 1e-9 from the uniform start and from an interior start with 99.9% of its
mass on one genotype. The generic path's three- and five-genotype tests still pass unchanged, which
is the other half of the plan's contract.

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` and
`cargo test --lib --bins --tests --all-features` in the container. Suite 3,517 → 3,518.
