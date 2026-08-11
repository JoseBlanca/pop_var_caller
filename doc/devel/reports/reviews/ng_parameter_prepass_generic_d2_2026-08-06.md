# ng step 4, D2 — review of the generic `NoiseModel`

**Date:** 2026-08-06. **Reviewed:** commit `fcb5906f`. **Fixes:** the commit after it.
**Agents:** three, each in its own worktree detached at `fcb5906f`, covering ten categories.

| agent | categories | outcome |
|---|---|---|
| structure | `module_structure`, `naming`, `idiomatic`, `smells` | 1 Major, 4 Minors, 5 nits — every fix applied and rebuilt before it was filed |
| numbers | `defaults`, `tooling`, `extras`, and every quantitative claim | 0 Majors, 4 Minors, 4 nits; the mathematics derived independently and found correct |
| reliability | `reliability`, `errors`, `refactor_safety` | **stalled after 30 mutations**; its harness and two rounds of proposed tests were recovered from its worktree and run here |

## Verdict

**One Major, eight Minors, and the mathematics is right.** The numbers agent derived the closed form
from spec §1 and arch §5.1 independently and checked the code term for term in both arms; it also
proved the pooled arm is the attributed arm summed over the forgotten splits, on paper for `G = 2`,
`n = 8`, `k = 3` and numerically to 5.6 × 10⁻¹⁷. Nothing in the scoring rule was wrong. What the
review found was in the seam's contract, in the tests' reach, and in the prose.

**⛦ The seam's row-width contract was written in a number only the SNP/indel path has.** The trait's
one method documented itself as appending `ploidy + 1` entries. That is a fact about *dosage*
genotypes, not about noise models: on the STR path a genotype is an unordered tuple of allele
lengths, so a diploid stratum spanning nine lengths has `A(A+1)/2` = **45**, and above diploidy
`C(A+P−1, P)` = 495 at nine lengths and four copies
(`spec/parameter_prepass_ssr.md` §4.2, verified). D3 is the step that turns the appended buffer into
a table, and `GenotypeLikelihoodTable::from_natural_logs` takes the row width as its **only** shape
argument — so the only width D3 could read off the trait was right for the first implementor and
wrong for the second. `from_natural_logs` cannot catch it either: 45 columns read as 3 still
divides, so the mis-shaped table is *accepted* and the climb runs on transposed rows. The trait now
carries `fn genotypes(&self, ploidy) -> usize`, asked of the model. It also repairs the one thing
appending gives up against clearing — under a clearing contract the row width is `out.len()` for
free, under an appending one it is a difference nobody records — and a new test pins the declared
width against the written one at every ploidy from 1 to 6.

**⛦ Two same-typed values transposed, and only one test in twenty-nine can see it.** The structure
agent swapped the two libraries' shares between their read groups in the fixture and ran the module:
**none of the three sum-to-one identities noticed**, and they cannot, because a rule with the shares
swapped is still a probability over the cell space. Only the harness fixture caught it, and only
because that fixture happens to live in a 90/10 world. Re-run here after the fixes, both
transpositions — shares swapped, and rates swapped — are killed by exactly one test, the harness
oracle. **This is the finding to carry into E1**, which computes those shares from read counts and
will have no identity oracle behind the pairing at all. It is also what makes `LibraryNoise`
load-bearing rather than tidy: holding the read group, the share and the rate in one struct is the
only thing preventing the pairing, and the recorded deviation that put the shares in `NoiseParams`
is justified more strongly than the report first claimed.

**⚠ My own justification for the fourth oracle was measurably backwards.** The test compared
genotype *differences* rather than absolute values, on the reasoning that the two `ln Γ`
implementations disagree in the last bits and those prefactors carry no genotype, so cancelling them
sharpens the comparison. The premise is true; the conclusion is not — cancelling a term makes a
check **weaker**. Measured over the fixture's 54 values, `libm` and the harness's Lanczos series
agree to **1.42 × 10⁻¹⁴** absolute, which is under six units in the last place and seventy times
inside the test's own 1e-12 tolerance, so the absolute comparison was available all along. And it
cost a kill: a rule that drops the library share `w_g` adds `−Σ_g k_g·ln w_g`, which also carries no
genotype and slips straight through a difference-only check. Both are now asserted.

**⚠ A wrong number in the author's own prose again — the fifth round in nine, and again a claim
about the reach of the author's own test.** "2,900 cells" for the depth-100-to-124 sweep is
**2,825**: `Σ_{d=100}^{124} (d+1)`. Three smaller ones went with it: "at depths 1 to 9" describes a
test that runs depths 1, 3, 5 and 9; the deviations section is headed "from the architecture" and
lists five items of which one is an addition to a research harness; and "the rest carry a handful"
is false of the four tomato samples carrying 7, 16, 16 and 42 libraries — the conclusion (1,683 of
1,707 fit two inline slots) is untouched.

**⚠ The sweep named for the negative-reference-read guard does not exercise it.** Its `alt_reads`
runs `0..=depth`, so the reference count is never negative and the assertion is never reached; what
those 2,825 cells assert is that every score is a log-probability. Deleting the guard fails only the
one-line refusal test. The report said the sweep was what stood between the model and the 5.2-rung
failure; it is the refusal test. Two further tests now bracket the tolerance from both sides.

**⚠ The sum-to-one identity closes at a whole depth, and a binned cell's mean is not one.** At a
fractional depth the cell space is a binomial series truncated at `⌊n⌋`: at depth 6.5 on two
libraries it comes to 0.998 and 0.061 at the two non-reference genotypes, against 1.0000000000000002
at every whole depth tested. **Not a defect in the rule** — the real table is not truncated, because
each cell carries its own mean rather than sharing one — but the doc invited a reader to carry the
identity over to the case cells are actually scored in, where it has not been established. What
stands there instead is the measured binning bias of the adopted ladder, 0.054 rungs and 0.3%
(research note §4.3). Said so now.

## The mutation record

The reliability agent stalled after 30 mutations, having reported two survivors. Its harness,
pristine copy and two rounds of proposed tests were recovered from its worktree; the tests were
applied and **the whole battery re-run here**, plus three of my own.

**Thirty-three mutations, thirty-two killed.** The two survivors it found are both closed:

- `REFERENCE_READ_TOLERANCE` widened from 1e-9 to **1.0** — now killed by
  `a_cell_a_millionth_below_its_alternative_count_is_refused`, which brackets the constant from the
  fault side while `a_cell_a_whisker_below_its_alternative_count_is_clamped_rather_than_charged`
  brackets it from the rounding side. Without the pair, widening the constant to swallow the
  per-bin-mean bug — the one that lands the fit 5.2 rungs low — would have passed.
- `single()` building its library under read group 0 whatever it was handed — now killed by
  `a_single_library_sample_carries_the_group_it_was_given_and_the_whole_share`. `single()` was also
  the one constructor that bypassed `new()`'s four invariant assertions; it now goes through them.

Nine tests were added from the recovered work: the two tolerance brackets, the two endpoint tests
(`ε` may legally be 0 or 1, where some category's probability is exactly zero and `0 · ln 0` lives),
the `−∞` test (a genotype that cannot have produced the cell is charged `−∞`, not merely a poor
score), the single-library sum-to-one at ploidy 1, 2 and 4, the share-tolerance refusal, the
cell-space-lists-each-cell-exactly-once check, and the declared-versus-written genotype count.

That last one deserves its own line: **the sum-to-one identity would be worth nothing if
`whole_cell_space` were not the cell space** — a cell listed twice and one missing cancel, and the
sum still comes to one. The numbers agent checked it independently by counting (a site of depth 5
over two libraries produces 16 cells: 1 pooled at zero alternative reads, 2+3+4+5 attributed, 1
pooled at five) and confirmed the helper builds exactly those; there is now a test asserting the
keys are distinct and each takes the arm the accumulator would give it.

## Everything else applied

- `Ploidy` printed with `Display` rather than `{:?}`, so the panic reads "a cell of ploidy 4" and
  not "Ploidy(4)". It was the only message in `src/ng/` diverging from a convention the crate wrote
  a `Display` impl and a test to establish.
- `share_weighted_rates` returns a named `ShareWeightedRates` rather than an unnamed `(f64, f64)`,
  whose two callers read it two different ways. All three transpositions of it were caught by the
  sum-to-one identity, so this is readability rather than a latent bug — but D3 adds a third caller.
- `single()` goes through `new()`.
- Bare adjectives named; the two `#[allow(clippy::cast_*)]` in the test helpers justified.

## Checked and correct

The closed form matches spec §1 and arch §5.1 term for term in both arms, with the reference factor
summed over every library in both — the term the plug-in got wrong. The real-`n` binomial
coefficient is the right Gamma extension. `alternative_read_probability` is exactly `ε` at zero
copies and exactly `1 − ε/3` at full dosage, swept over every ploidy 1 to 255 and five error rates,
with **zero** error at the ends. The 18-row `HARNESS` fixture is bit-identical to a fresh
regeneration — all 54 literals — so oracle 4 is not the author checked against the author's
transcription. `--only=oracle` is an exact-equality guard, not a suffix match, so the harness's
measured output could not have moved; it was re-run to completion and is unchanged. Verified
figures: 68%/78% (research note, about the average-share plug-in); 1,550 of 1,707; 35 splits;
5.2 rungs and 0.3%; 31 worlds; ~75,000 cell-scorings per fit; D2 is the fourth of the six
silent-failure steps.

## Verification

Container throughout. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → **3,089** in the library binary
(from 3,080), 69 across the nine other binaries unchanged; `cargo test --doc
ng::parameter_estimation` → 1 passed; `cargo doc --no-deps --lib` at 12 unresolved links, none in
`parameter_estimation`. `ng::parameter_estimation` 178 → **187** tests, `noise_model` 20 → **29**.

## Left for the owner, and for D3

- **`fit_mixture_weights` should become `pub(crate)` in D3** (carried from D1), and D3 should assert
  the model's declared `genotypes()` against what it appends, per cell.
- **The three factorial prefactors are ~98% of the inner loop's arithmetic and never change.**
  Measured by the structure agent: three `lgamma` at 26.1 ns against 0.5 ns for the rest, so a
  diploid cell spends ~78 ns per rung on numbers that are identical at all 161 rungs — roughly 6 ms
  of prefactor against 0.1 ms of everything else per scan, before E2's up-to-20 alternations
  multiply it. They cannot be hoisted through the current trait, because the model is called once
  per (cell, rung) with nowhere to keep a per-cell scratch. **D3 is the last cheap moment.** Three
  options, in the review file; dropping them is explicitly *not* recommended, because the
  sum-to-one tests would then be testing a different function from the one the scan runs.
- **`dump_attributed_oracle` hard-codes its world's rates and shares** rather than looking the world
  up in `worlds()` by name, as the harness's other sections do. If the world list is
  re-parameterised the dump keeps printing the old world under the old name.
