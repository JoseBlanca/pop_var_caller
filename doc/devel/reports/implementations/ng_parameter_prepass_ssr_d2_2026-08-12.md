# ng step 4, the STR path — D2: how likely each genotype makes a locus's shape

*Implementation report, 2026-08-12. Step D2 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md) — **own commit, do not
bundle**, because its failure is silent. Four review agents; the reliability agent ran 24 mutations
of which 3 survived, and I ran 3 more after the fixes to confirm each bites. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1, §4.2, §10;
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §3. Ported from
[`examples/shared/stutter_model.rs`](../../../../examples/shared/stutter_model.rs)'s
`read_bucket_probs` and `genotype_bucket_probs`.*

## What the step is

`SsrNoiseModel` and its `NoiseModel` implementation: given a locus's shape — how its reads fell
across the nine offset buckets — how likely does each genotype make it? A genotype is an unordered
tuple of allele *lengths*, so a diploid stratum whose alleles span thirteen lengths has 91 of them.
Each read picks one of the locus's copies with equal chance and then slips, so a bucket's
probability is the average of the copies' own bucket distributions, and the shape's probability is
one multinomial over the buckets.

**An end bucket's probability is the sum over every offset it absorbs, never the probability of
sitting exactly on the edge.** That is the whole of D2, and it is why the step lands alone.

## Why this one commits alone

The plug-in is the tempting shortcut and it is wrong twice over, neither of them visible from
outside. Scoring at the edge fails the sums-to-one identity — the buckets come to 0.9488 at a
recorded range of ±1 — and on a stratum whose alleles reach three copies either side it returns the
slippage level **52% low**. Rescaling it so the buckets do sum to one repairs the identity and not
the bias: it then runs **33% high** where 30 in 100 slipped reads take a second step, which is the
regime long tracts sit in. Both figures are the harness's own, re-run at the top of this step and
recorded in D1's report.

## Recorded deviations from the architecture

1. **`type Cell = StratumCell`, where arch §3 sketches `type Cell = LocusShape`.** A `LocusShape`
   cannot be a `fitting::WeightedCell`: it is the `BTreeMap`'s *key*, so it knows neither how many
   loci showed it (the value) nor what ploidy they sit at. `StratumCell` pairs a `StratumEntry`
   with a `Ploidy`. **The sibling path made the identical move** — `generic/histogram.rs` stores
   `CellCounts` and hands out `Cell` through `cells(ploidy)` — and `fitting/mod.rs` already records
   that `WeightedCell` replaced the `(Cell, Ploidy, u64)` tuple the generic architecture sketched.
   The arch sketch predates the trait. **One line of arch §3 is stale and wants the owner's
   approval to correct** (raised at Checkpoint D).
2. **`append_genotype_likelihoods`, not `genotype_likelihoods`** — the trait's real method name,
   and its contract is to append rather than to fill.
3. **The slip kernel is evaluated once per cell rather than once per allele length.** The kernel
   renormalises its truncation over eight powers of the fall-off on every call, and a cell's score
   walks it once per allele length, so hoisting it turns 1,872 `powi` calls per cell into 144.
   Bit-identical — same distances, same values, same accumulation order — and every test passed
   unchanged across the change.

## What the review changed

**Blocker — `WeightedCell::sites()` was exercised by nothing, and the error it guards against is
the design's central one.** Two agents found it independently. Replacing `sites()` with the shape's
read *depth* left all 567 tests of this path green — and that substitution is exactly what turns a
table of loci back into a tally of reads, which is the keying spec §4.1 rejects because the fitted
slippage level then moves **333-fold depending only on where the search starts**. Nothing inside
D2 calls either trait method: the ploidy check inside the scoring resolves to the *inherent*
`StratumCell::ploidy`. A test now goes through the trait, on a fixture whose locus count, whole-repeat
depth and guard count are three different numbers, so every wrong answer differs.

**Major — the heterozygote test could not tell an average from a sum.** It asserted
`L(a,b) = ½(L(a,a) + L(b,b))` at a depth of one read, where the multinomial's arrangements term is
zero — so the identity holds for *any* constant multiple of the truth, and a rule that summed the
two copies instead of averaging them passed it. Raising the depth cannot fix it: the identity is
false at two reads even for correct code. It now also anchors one homozygote against the
single-copy distribution it is built from, which fixes the scale. I reproduced the sum mutant: it
now fails five tests including this one.

**Major — the rescaled plug-in's only reliable catcher was one row that looked like padding.** The
sums-to-one gate is the wrong instrument for a *proper* rule, and it catches this one only through
its zero-level row, where the rescaling divides a distribution that is one at the allele's own
bucket. That row now carries a comment saying so, measured.

**Major — the harness-agreement test called the library's own kernel on both sides**, so it pinned
the bucket-clamping loop and nothing else. The reference's `Slip::p` is now transcribed too, with
its closed-form renormaliser, and the comparison runs over the *genotype* bucket probabilities the
plan names as the oracle rather than only the per-read ones.

**Minor — `MAX_ALLELE_LENGTHS` was pinned in one direction only.** Undersizing panics on an index;
oversizing was silent. Now an equality.

**Three renames**, all from the design review and all because the name did not say what the value
is: `MOST_ALLELE_LENGTHS` → `MAX_ALLELE_LENGTHS` ("most" reads as "the majority of"),
`multiset_count` → `genotype_count` (combinatorics jargon for a value that *is* the genotype count),
and `allele_offsets` → `allele_support` — the important one, because "offsets" reads as the offsets
*observed* when the point is that these are the ones the fit is *permitted* to place mass on.
`for_each_genotype` became public: a fit returns one frequency per genotype in that order, so
whatever reports those frequencies has to walk the same one.

**The mathematics itself was found correct**, verified numerically rather than by reading: agreement
with the reference over 320 parameter settings × 5 strata × every allele offset, worst disagreement
3.3e-16; `genotype_count` exact against a `u128` Pascal triangle over 182 cases; `for_each_genotype`
matching an independent generator over 52 cases; no `NaN`, panic or out-of-range index from any
parameter corner.

## ⚠ Five wrong claims of mine, and four imprecise ones

Every figure quoted from the design documents was correct — 0.9488, 52%, 33%, and all of the
combinatorics. Every wrong claim was my own arithmetic or my own account of my own code.

1. **"a tetraploid's 495 genotypes"** in the ploidy-mismatch test — that test's stratum has thirteen
   allele lengths, so `C(16,4) = 1,820`. 495 is the count at *nine* lengths, which is the sibling
   test's stratum.
2. **"the loop would otherwise run forever or index out of range"** — measured, without its guards
   it emits one bogus genotype and returns. What the guards prevent is a *width disagreement*:
   `genotypes()` would promise zero entries while the loop pushed one, and the scan sizes its table
   from that promise.
3. **"the four checks the design puts before anything is fitted"** — the design says three, the
   harness prints three, and the same patch said three sixty lines earlier. The harness's fourth
   block is the control *fit*, not a pre-fit check.
4. **"keeping the ploidy on the cell is what would let a two-ploidy fit"** — necessary, not
   sufficient, and the sentence hid the real gap. A table is keyed by `(read group, stratum)` with
   no ploidy in it, so two loci of different ploidy that showed the same shape already collapse
   into one entry and no later pairing can split them.
5. **"declares 99.85% of this genotype's reads impossible"** — the arithmetic held (0.0015 of the
   reads survive) but the mechanism did not. The plug-in charges the bucket those reads land in 670
   times too little rather than nothing; only two of the nine buckets are literally impossible.

Imprecise and corrected: "by a factor of thirty" holds only at the widest stratum (fifteen at the
narrowest); "an empty product: zero" reads as the empty product *being* zero, where the value is a
likelihood of one and `ln L = 0`; "the plug-in asks how often a read sits exactly six copies from
where it was" — it asks that of one bucket and slips of two to ten across the nine; and "the fit is
what holds the ploidy map" ignored that the accumulator visibly holds one too.

## ⛦ Raised for Checkpoint D, not fixed here

**A stratum's table can pool loci of different ploidy with no way to tell them apart at fit time.**
`SsrAccumulators` keys by `(ReadGroupId, Stratum)` and never calls `ploidy_at`, where the sibling
accumulator keys by `(read group, ploidy)` and does. Nothing is wrong today — `ConstantPloidy` is
the only `PloidyMap` there is — but the sibling documents this exact failure at
`generic/histogram.rs`: *"haploid sites, which can never be heterozygous, entering the
heterozygosity fit as diploid ones. A wrong fitted rate with nothing to show for it."* The fix is
one key, `(ReadGroupId, Stratum, Ploidy)`, and it reopens a committed step (C5), so it is the
owner's call.

## Tests

Twelve new, twenty-six in the file, 138 across the STR path.

| test | what it pins |
|---|---|
| `the_scoring_rule_sums_to_one_over_the_entry_space` | gate one, per genotype, over four depths and four parameter settings |
| `the_scoring_rule_sums_to_one_at_a_tetraploid_stratum_too` | the same at ploidy four, where the average and the enumeration generalise |
| `scoring_an_end_bucket_at_its_edge_fails_the_sums_to_one_gate` | the shortcut, rejected by the mechanism rather than by a bound |
| `no_bucket_is_charged_a_negative_share_of_a_locus_reads` | gate two |
| `a_silent_kernel_puts_every_read_on_its_own_allele` | gate three, at all thirteen alleles |
| `the_bucket_probabilities_agree_with_the_harness` | the reference's own `Slip::p`, `read_bucket_probs` and `genotype_bucket_probs`, transcribed |
| `a_cell_counts_the_loci_that_showed_its_shape_and_not_their_reads` | loci, not reads |
| `a_heterozygote_is_the_average_of_its_two_copies` | the average, and its scale |
| `a_stratum_has_one_genotype_per_unordered_tuple_of_its_own_allele_lengths` | 45 / 66 / 91 / 13, and the scratch array's width |
| `every_genotype_is_emitted_once_in_non_decreasing_order` | one spelling per genotype, ploidy 1 to 3 |
| `the_guard_reads_do_not_enter_the_length_likelihood` | the factorisation |
| `a_cell_cannot_be_scored_against_another_ploidy_genotypes` | the assertion nothing else reaches |

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` and
`cargo test --lib --bins --tests --all-features` in the container. Suite 3,505 → 3,517.
