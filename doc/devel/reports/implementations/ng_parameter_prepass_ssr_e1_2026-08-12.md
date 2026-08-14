# ng step 4, the STR path — E1: the substitution rate, per stratum

*Implementation report, 2026-08-12. Step E1 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1 and §4.2,
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §4.1.*

## What the step is

The first of the four fits, and the only one that is not a search: mismatched bases over compared
bases, one division per stratum. `substitution_rate_of(&StratumTable)` turns one table into an
`Estimate<ErrorRate>`; `substitution_rates(&SsrAccumulators)` walks every stratum and returns the
map. It goes first because it needs none of the other three — a read's mismatch count is binomial
at this rate whatever tract length the read showed, so the two channels factorise and this one
closes in a division (spec §4.1).

The arithmetic was already there from B3 (`StratumTable::substitution_rate`), proven against a
grid of 100,000 rates. What E1 adds is the per-stratum plumbing, the warrant beside each rate, and
the oracle that the plumbing keeps each stratum's bases to itself.

## The step it was preceded by: the tables gained a ploidy

*Committed separately (`98430bcf`), because it changes an object Milestone D built.*

`SsrAccumulators` keyed by `(read group, stratum)` and held a `PloidyMap` it never asked anything,
so one table could pool a haploid locus with a diploid one. Invisible while counting — an entry is
a count of reads at each offset and knows nothing about chromosomes — and wrong at the fit, which
scores every entry of a table against the genotypes of one ploidy. The key is now a named
`StratumKey { read_group, stratum, ploidy }`, which the four fits, their provenance lists and
`SsrSampleParameters::by_stratum` all share. Owner's call, 2026-08-12: do it before the fits are
built on it. `ConstantPloidy` remains the only map, so every key today carries the same value and
today's tables are unchanged.

## Recorded implementation choices

1. **Two functions, not one.** The table-shaped `substitution_rate_of` is what E4's merge needs — a
   merged pair of strata is refitted from a table that is in no accumulator — and the map-shaped
   `substitution_rates` is what a driver walking the accumulator needs. Neither is a wrapper worth
   deleting.
2. **The count beside the rate is bases compared, not loci.** `Estimate::observations` is
   documented as "reads for a per-read rate, sites for a per-site one"; this rate is per base. A
   locus count would let a hundred shallow loci outrank one deep one on a quantity neither of them
   measures. `StratumTable::bases_compared()` was added to expose it.
3. **A stratum with nothing compared is absent from the map**, not present at zero — the
   `Option` rule `substitution_rate` already states, carried up. Filling such a gap is the
   borrowing step's business (E3) and it marks what it does.

## What the review changed

Three agents: correctness, numbers, and a mutation run over the new code.

**The count beside the rate is bases, and the field it will be reported in said loci.**
`StratumFitSummary::low_slippage_substitution` — and `arch` §4.3's sketch of it — promised "how
many loci stood behind it", where the code emits compared bases. The deciding argument is the
comparison the number exists for: this rate is to be read against the SNP/indel path's rate for
the same library, and that path counts an error rate's observations as reads times the sites they
covered, which is one base per read per site. Two warrants on different scales cannot be compared.
The code is right and the docs were fixed; **`arch` §4.3 is now behind the code and is on the
checkpoint list**.

**The doc promised a fallback the design does not define.** It said that where a stratum has no
compared bases, filling the gap is the borrowing step's business — but the borrowing the design
defines is over the slippage model, not over this rate, and `StratumFit::substitution` is not an
`Option`. Nothing observed can reach the case (a read reaches a table only through a complete
witness, which compares its bases), so this is a contract gap rather than a live defect; the doc
now says so instead of promising a rung.

**Six wrong claims of mine, all in prose, none in the arithmetic.** Every figure about the
fixtures held — 10,000 and 30 compared and mismatched, 4,000 and 120, 10,000 and 300, the pooled
0.0107, the exact-float reasoning. What was wrong: "forty times too small" for 41.7; "two
libraries sequenced the same tracts here" where the fixture has two *different* loci of one
stratum, one library each; "`ConstantPloidy` is the only map that exists", contradicted 760 lines
below in the same file by the test's own two-ploidy map (production's only map is what was meant);
`generic/histogram.rs` cited for the sibling path's ploidy keying, which is in
`generic/accumulators.rs` — histogram.rs says the opposite of its own type, "a table carries no
ploidy of its own"; a claim that a locus count is something "nothing else downstream could" tell
apart, where two public accessors on the table do; and the merge-refusal test's rewritten doc,
which described two shards *disagreeing* about the genome where the fixture hands it two equal
constant maps. That last one mattered most: the guard is pointer identity and the test proves the
identity check, not the disagreement. Its message, its name and its doc now say which.

**And the fixture helper's guarantee was stated against the wrong reference.** It asserted one
differing base against the locus's reference bases, while the accumulator compares each read
against the motif tiled from its first base. The two agree only because every call site passes a
perfect tiling; the assertion is now over the tiling, so a fixture that broke the assumption fails
there rather than reporting a rate no test intended.

**Seven of twenty mutations of the new code passed the whole module**, and three fixtures closed
all seven. They are one pattern: *the fixtures were all alike*.

- **The test named for absent-rather-than-zero could not see it.** Its accumulator held **no
  strata at all**, and an empty accumulator is an empty map under every implementation ever
  written. Two mutants lived there: one turning a stratum with nothing compared into a measured
  zero, one deciding on the locus count instead of the base count — indistinguishable, because
  `StratumTable::default()` has neither. The new fixture is a tract every read shows as entirely
  deleted: the reads witness it completely, so the locus files a shape, and they show no bases,
  so nothing is compared. The table then exists with one locus and zero compared bases, which is
  what the two mutants need to be visible.
- **No fixture had a stratum whose rate is zero**, so a floor of 0.0001 under the rate survived,
  and so did a filter dropping the strata that measured zero — the very distinction the `Option`
  exists to protect.
- **Every fixture compared at least 4,000 bases and every ploidy was 2**, so a minimum-evidence
  gate below 4,000 was invisible, a provenance decided by "at least 10,000 bases" was invisible,
  and dropping the ploidy from the output key changed nothing. The clean-stratum fixture compares
  100 bases, which is a hundredth of the others, and the new two-ploidy fixture puts the same
  stratum of one library at ten-fold different rates on one and on two genome copies.

**I re-ran all seven mutants after applying the fixes** rather than taking the agent's word: the
zero-instead-of-absent and locus-count guard fail the deleted-tract test; the dropped ploidy fails
the two-ploidy test at `rates.len()`; and the evidence gate, the floor, the zero-filter and the
evidence-decided provenance each fail the clean-stratum test, at four different assertions.

## Tests

Seven new.

| test | what it pins |
|---|---|
| `a_stratums_substitution_rate_is_the_rate_its_reads_carried` | 30 mismatched bases in 10,000 compared come back as 0.0030 exactly, with 10,000 as the warrant and `FittedHere` as the provenance |
| `two_strata_of_one_library_keep_their_own_rates` | 0.0030 and 0.0300 at two strata of one library; a pooled fit reports 0.0107 for both |
| `two_libraries_at_one_stratum_keep_their_own_rates` | 0.0030 and 0.0300 at two libraries of one stratum |
| `two_ploidies_of_one_stratum_keep_their_own_rates` | and at one stratum of one library on one and on two genome copies |
| `a_stratum_whose_reads_all_matched_reports_a_measured_zero` | zero is a measurement: not floored, not filtered, and not withheld for thinness at 100 compared bases |
| `a_stratum_with_loci_and_no_compared_bases_is_absent_rather_than_zero` | the case the title of the next test claims and cannot reach |
| `a_table_with_nothing_compared_yields_no_rate_at_all` | absent rather than zero, at the table and at an empty walk |

Two more, from the ploidy key that preceded the step: `two_ploidies_of_one_stratum_are_two_tables`,
and the merge guard's fixture, whose name and message now say what it actually proves — that two
shards built on **separate** map objects are refused even where the two maps agree, because the
check is pointer identity.

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` and
`cargo test --lib --bins --tests --all-features` in the container: 3,441 → **3,448** lib tests, 0
failed, 11 ignored. Also run **natively on the host** — the first host run on this branch — 3,445
passing at the ploidy-key commit, 210 s.
