# ng genotype prior — D3: the run's pair spread over one locus's alleles

*Implementation report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step D3 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone D, on top of D2
(`0c024bdb`).*

## 1. What it is

`fill_locus_concentration` takes the two numbers step D2 read off the panel's frequency spectrum —
the reference allele's concentration and the total belonging to the alternatives — and writes one
locus's concentration, one entry per allele:

```text
  out[0]    = α_ref,       floored at MIN_ALT_CONCENTRATION
  out[1..]  = α_alt_total / (number of alternative alleles),
                           floored at MIN_ALT_CONCENTRATION
```

**Sharing the total out rather than repeating it keeps a site's total polymorphism independent of
how many alleles it happens to carry** ([`spec/calling_priors.md`](../../ng/spec/calling_priors.md)
§4). A triallelic site is not twice as polymorphic as a biallelic one merely for holding a third
allele — that would be a statement about the candidate generator rather than about the genome.

It is the port of production's `alpha_from_diversity`
([`genetics.rs`](../../../../src/genetics.rs)) with the fitted pair as input instead of `α_ref = 1`
and `α_alt = θ` hard-coded.

## 2. Three differences from the thing it ports

| | `alpha_from_diversity` | `fill_locus_concentration` |
|---|---|---|
| where the pair comes from | `α_ref` fixed at 1, `α_alt` the diversity | both fitted, from D2's projection |
| what it returns | a fresh `Vec` per call | fills a caller-owned slice — nothing allocates per locus per pass (spec §8) |
| the reference's floor | none: the constant 1 is above it by construction | floored, because `SpectrumSeed::new` admits any strictly positive value and `Concentration`'s invariant is that **every** entry clears `MIN_ALT_CONCENTRATION` |

The floor on the alternatives, its value, and the monomorphic-locus behaviour are carried across
unchanged.

## 3. Two departures from arch §4's signature, and one I made and the review took back

Arch §4 sketches `seed_for_locus(seed: &SpectrumSeed, allele_count: usize, out: &mut [f64])`.

- **Named `fill_locus_concentration`**, matching the module's two other buffer-fillers,
  `fill_expected_spectrum` and `fill_sample_concentration` — the second of which is step C1's own
  rename for the same reason, already ratified and already in arch §3.1. `seed_for_locus` reads
  as a getter returning a seed; it returns nothing and fills a buffer.
- **The seed by value, not by reference.** `SpectrumSeed` is `Copy` at 24 bytes and its accessors
  already take `self`.
- **`class: VariantClass` added**, which the plan's D3 line mandates and arch §4's sketch predates.

**`allele_count` was dropped and then put back, and the reasoning is worth keeping.** I dropped it
on the argument step D2 used for the panel size: `out.len()` is the allele count, so a second
argument is a second place for it to disagree. **That analogy is false here**, and one reviewer
caught it while another endorsed the drop. D2's class weights *are* the data — their length is the
panel size and nothing else. `out` here is a slice of scratch the calling loop owns and reuses, so
its length is a slicing decision. Measured with the check absent: a locus of three alleles handed
an eight-wide buffer at its full length gives every alternative `θ/7` instead of `θ/2` — 3.5 times
too little prior mass, about 5.4 phred off every non-reference genotype, every value still looking
like a value. A buffer *shorter* than the locus is worse, because the entries past its end keep the
previous locus's concentrations and every length check downstream passes.

The check is only load-bearing when the two arguments come from different expressions — the
locus's own allele count against the buffer the worker sliced — and nothing in the type system can
enforce that. It is what arch §4 specified, it costs one integer compare per locus, and dropping it
was mine to get wrong.

## 4. The class argument fixes the shape of the call and settles nothing else

Two of the reviewers pushed on this from different directions and converged. Both are recorded
here because the next step to touch it is the calling loop's, not this plan's.

**It cannot absorb the split, because the diversity arrives inside the seed.** Spec Q1 asks that
*the concentration function* take the class so that splitting later does not touch every call
site. But a run that fits two diversities fits **two seeds**, and whoever calls this has already
chosen between them. The seam that could absorb a two-class pre-pass without a signature change is
`project_spectrum_seed`, which carries the same argument one step earlier. **Which of the two ends
owns the split is not settled by any document**, and applying it at both would apply production's
8:1 ratio twice. Worth settling before the loop freezes its parameters: `FrozenParameters`
([`calling_loop.md`](../../ng/impl_plan/calling_loop.md)) carries **one** `SpectrumSeed`, which
structurally commits the run to the split happening here.

**A locus can carry both classes at once, and no code in ng can tell them apart.**
`LocusKind::Generic` is documented as *"a SNP/indel candidate site"* — one variant covering both —
so the class has to come from the sequences: comparing the reference's length against each
alternative's, which is production's `is_indel` test written on alleles instead of CIGARs, and
which does not exist in ng. A reviewer built a locus with reference `AT` and alternatives `GT` and
`A`; nothing refuses it. So the calling loop, today, has no rule to follow and no helper to call.

**None of this is live.** The argument is unread, and
`both_variant_classes_get_the_same_concentration_today` pins that both values give byte-identical
output, so the loop may pass either with no effect on any measurement.

## 5. What the reviews found

Three agents in isolated worktrees — conformance with prose, the calling loop's seam, and
correctness with mutations. The last was cut off by a rate limit and re-run; it ran twelve
mutations, eleven killed, one proved to change no behaviour, no behaviour-changing survivors.

**The one test of the reference floor was checking nothing under the gate.** Flooring `α_ref` is
this step's one departure from spec §4's written formula, and its only test was
`let _ = Concentration::new(&out);` — but `Concentration`'s per-entry check is a `debug_assert!`
and `Cargo.toml`'s `[profile.release]` sets no `debug-assertions`, so under
`cargo test --release` the only thing that construction still rejects is emptiness. The assertion
is now written out, and the type still constructed beside it so the two stay in step.

**Three doc claims were wrong in the same direction — they justified a guard by a reason that
cannot happen.** The reference floor was explained by the value being *fitted*; a fit cannot make
it bind, because D2's search box bottoms out near `α_ref = 1e-5`, ten million times the floor. What
the floor is really for is that `SpectrumSeed::new` admits any strictly positive value, so a
hand-built seed can sit below it. Likewise "nine orders below the least diverse panel anyone would
call" quoted a human figure as a lower bound on diversity, which is the move `CLAUDE.md` warns
about and the wrong direction besides.

**The `§8` citation against a per-allele class array was a rationalisation.** Spec §8 forbids
allocation *inside the per-sample loop* and explicitly mandates caller-owned scratch sized by
allele count — which is what such an array would be. The honest argument needs no §8: the class per
allele is derivable from `CandidateAlleles`, which hold the bases, on the day Q1 needs it.

**Two doc claims about precision named the wrong case.** The test's tolerance was called "four
orders inside the tightest case (a `θ` of 1 in 10,000 over 5 alternatives)" — but the bound was
*absolute*, so headroom shrinks as `θ` grows and the tightest cell is the largest `θ`, where it is
two orders rather than four. And "the split is exact only where the number of alternatives divides
the total in binary" is not true of this test's own inputs: measured, the error is exactly 0.0 in
all twelve cells, 3 and 5 alternatives included. The bound is now relative — `8 · ε · θ` — which
holds at every total `SpectrumSeed` admits, where the absolute one would have failed on a legal
fitted seed (the sum's worst-case rounding over 9 alternatives at a total of 1e3 is 1.1e-13), and
is four orders tighter at the diversities the test actually runs.

**The monomorphic early return is documentation rather than a guard.** Removing it changes no
output for any input — `out[1..]` is empty there, so the fill writes nothing whatever the division
produced. It stays, with a comment saying which it is.

### Raised, not taken: the seed should be a `Concentration`

A reviewer wrote the calling loop's inner body against the shipped API and found that **skipping
this call entirely is caught by nothing.** A zeroed buffer passes `fill_sample_concentration`
(which takes `&[f64]`), passes `Concentration::new` in release, and reaches the prior's row as
`lgamma(0)` — the row comes back `[NaN, −inf, NaN]`, which is exactly what the
`GenotypePriorModel` contract forbids by name. Downstream the locus runs to the pass cap and is
emitted as unconverged with nothing saying why.

Every *shape* mistake on this path is caught in release; the *omission* is caught nowhere. The fix
is to have this function return a `Concentration<'_>` over what it filled and
`fill_sample_concentration` take that instead of a bare slice — the seed satisfies every one of
`Concentration`'s invariants already. **It changes step C1's committed signature and arch §3.1**,
so it is raised for the milestone review rather than taken here.

## 6. Tests

Nine, all cheap — the function is a division and a fill.

| test | what it pins |
|---|---|
| `a_locus_carries_the_same_total_polymorphism_however_many_alleles_it_has` | spec §4's reason for sharing the total out: the alternatives sum to the run's total to within `8 · ε · θ` at 2, 3, 4 and 6 alleles and three diversities — measured, exactly 0.0 in all twelve cells |
| `the_reference_allele_takes_the_first_entry` | the ordering the concentration is read in — reversing it would tell every genotype row the reference is the rare allele |
| `a_locus_with_no_alternative_allele_gets_only_the_reference` | a monomorphic locus is an answer, matching the ported function |
| `a_cohort_with_no_polymorphism_still_clears_the_floor` | the floor bites at zero diversity, and is 50 million times below a real one — so it cannot perturb an estimate |
| `every_entry_clears_the_floor_the_concentration_type_requires` | the output really is a legal `Concentration`, checked by an assertion that survives `--release` rather than by a construction that does not |
| `both_variant_classes_get_the_same_concentration_today` | Q1's seam: the argument exists, the behaviour does not yet split |
| `a_locus_with_no_alleles_is_refused` | the assertion is what makes the panic name the caller's mistake rather than an index |
| `a_buffer_longer_than_the_locus_is_refused` / `…shorter…` | the mistake the calling loop is most likely to make, in both directions |

## 7. Corrections owed to the design documents — raised, not applied

1. **Arch §4's `seed_for_locus` signature** — the name, the seed by value, and the `VariantClass`
   argument the plan mandates. `allele_count` stays as §4 wrote it.
2. **Spec §4's formula floors only `α_alt`.** The code floors `α_ref` too, because
   `SpectrumSeed::new` admits any strictly positive value while `Concentration` requires every
   entry to clear the floor. One clause records it.
3. **Neither spec Q1 nor arch §4 says which end owns a class split** — the projection or the
   per-locus expansion. §4 above.
4. **`CallingScratch` (`arch/calling_em_loop.md` §2) is two buffers short and one type wrong** for
   the chain this module now exposes: it needs a per-locus seed buffer (which cannot be the
   per-sample `concentration`, since the seed is re-read for every sample), the per-allele scratch
   `PriorRow::new` requires, and `posterior_row` retyped from `Vec<f64>` to `Vec<LogProb>` — the
   prior writes `LogProb` and there is no cast without `unsafe`, which the crate forbids.
5. **`calling_loop.md`'s zero-allocation test cannot be written as specified.** `Cargo.toml`
   forbids `unsafe_code` crate-wide, so no counting `GlobalAlloc` can be installed; the check needs
   a separate dev-dependency crate or a different oracle.
