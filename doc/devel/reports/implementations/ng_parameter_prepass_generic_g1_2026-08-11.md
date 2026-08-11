# G1 — the model-free anchors, and the first evidence that is not our own model

**Step:** G1 of `impl_plan/parameter_prepass_generic.md`, Milestone G. Held through step 4a
because the numbers it asserts are the ones that milestone changed.
**Date:** 2026-08-11. **Code:** `generic/truth_anchors.rs`.

## Why this file exists

Every recovery test in this module generates its data from the model it then fits, so a shared
misspecification cancels and the test passes. Those tests catch gross bugs and **cannot catch
bias**, which is the failure step 4 exists to remove (`arch` §9). G1 is the first evidence that
owes the model nothing: three numbers counted off the GIAB benchmark and the reads.

## What it counts, and the one thing it asserts

- **The error rate** — disagreeing reads over total reads at loci the benchmark calls
  homozygous reference. Reads and not bases: a read carrying a twelve-base deletion disagrees
  once.
- **Heterozygosity** — loci the benchmark calls heterozygous over the loci walked.
- **The every-copy-non-reference rate** — the same over loci where no copy is the reference's,
  which is `1/1` and `1/2` together, because the fit counts non-reference *copies*.

**One assertion: the fitted error rate is not below the counted one.** The confident regions are
the easy regions, so a rate counted there is a floor for a rate fitted over the same loci — the
fit sees every read the count saw and cannot see fewer disagreements than were there. Every
other number is printed, because heterozygosity from easy regions is depleted of the sequence
where variation concentrates and a fitted value above it is expected rather than wrong.

## The measurement

| | model-free | fitted | apart |
|---|---|---|---|
| HG002 30x, error rate | 2.3300 × 10⁻³ | 2.3327 × 10⁻³ | **+0.1%** |
| HG002 300x, error rate | 2.4963 × 10⁻³ | 2.5039 × 10⁻³ | **+0.3%** |

38,665 disagreeing reads of 16,594,596 at 30x; 169,437 of 67,875,294 at 300x. A tenth of a
percent is a fortieth of one rung of the error-rate ladder.

| | truth | fitted | ratio |
|---|---|---|---|
| heterozygosity, 30x | 9.9666 × 10⁻⁴ | 1.0813 × 10⁻³ | 1.085 |
| heterozygosity, 300x | 1.0035 × 10⁻³ | 1.1465 × 10⁻³ | 1.142 |
| every-copy-non-reference, 30x | 5.7444 × 10⁻⁴ | 5.3804 × 10⁻⁴ | 0.937 |
| every-copy-non-reference, 300x | 5.7631 × 10⁻⁴ | 5.4564 × 10⁻⁴ | 0.947 |

**Two of these were already known and one was not.** 9.9666 × 10⁻⁴ and 5.7444 × 10⁻⁴ are, to
every digit, what `research/noise_model_overdispersion_2026-08-10.md` reached independently in
Python through `bcftools`. This code shares nothing with it. The 300x truth heterozygosity had
never been measured — step 4a's reports derived roughly 1.0039 × 10⁻³ from two of that note's
ratios, and the direct count says **1.0035 × 10⁻³**.

**The error-rate agreement is new and is step 4a's doing.** Before that milestone the fit
returned 2.239 × 10⁻³ against a model-free 2.263 × 10⁻³ measured through `samtools mpileup`.
Both numbers moved: the fit because the second class of site now takes its share, and the count
because it is now taken from ng's own reads under ng's own filters rather than from a second
tool that had to be kept in step. They now agree to a tenth of a percent.

## Two decisions worth stating

**Per locus, not per base** (owner, 2026-08-10). `arch` §9 words the denominator as *the
confident regions' length*. The anchor counts loci instead, classing a locus heterozygous when
any of its reference positions carries a heterozygous record, because the fitted heterozygosity
it bounds is a per-locus rate — a per-base truth and a per-locus fit would differ in their unit
before they differed in anything interesting. It is also the rule the research note used, so the
anchor and the ratios step 4a already reports are one measurement rather than two that agree.

**The reads are ng's own.** The research note reached its count through `samtools mpileup`, which
then had to be talked into matching the walk's read filters. The walk already reports, per locus,
how many reads there were and how many disagreed — so the count comes from exactly the reads the
estimator is fitted on, with no second tool to keep in step. Above the depth cap that is fewer
reads than the alignment holds (67.9 M over 549,180 loci at 300x is 123.6 apiece, essentially
every locus drawn down to the cap of 124), and it is still the right comparison, because the fit
sees the identical draw from the identical seed.

## What the review found, and it was the important half

Three agents' worth of habit says reading finds nothing here, and mutation did the work again.

**The anchor could pass while measuring nothing, by two separate routes.**

1. **The single assertion sat inside a `for` loop**, and a loop is a comparison that can run
   zero times. Making it iterate nothing left the anchor green with three console lines of four
   and no assertion run. **My first fix was wrong**: I asserted the map was non-empty, re-ran the
   mutation, and it still passed — the map was never the problem. The fix is to not loop. One
   read group is asserted (the count is pooled over the whole site, so comparing it against
   several libraries' rates would not be like for like), so the rate is taken directly.
2. **The two truth classes were interchangeable.** Reading `0/1` as every-copy-non-reference and
   `1/1` as heterozygous left the error rate byte-identical — both classes leave the
   homozygous-reference denominator alike — and moved only the two printed ratios. Nothing
   asserted the classification at all. Now caught by biology rather than a tolerance: HG002 is
   outbred, so heterozygous loci outnumber those with no reference copy, 550 against 317, and
   the swap inverts it.

**And the counting is blind in one direction**, which is sharper than anything the module doc
had said. Over-classifying pulls loci out of the homozygous-reference denominator and *lowers*
the model-free rate, making the assertion easier; padding every truth record by ten bases
inflates heterozygosity twenty-fold and still passes. The headroom is 0.12% at 30x and 0.3% at
300x. What guards that direction is the classification checks above, not the assertion.

**Three smaller things, all mine.** The printed record count was 4% high — 947 against the 909
the benchmark holds over these spans — because a record straddling two of the 3,142 typed
regions came back from both queries; it now counts distinct records. The claim that a record
occupies its whole reference allele turns out to change nothing, since the pileup already widens
a locus to an indel's span, so it is documented as belt-and-braces rather than load-bearing. And
the module said *bases* where it counts reads.

**Mutations that were already caught**, for the record: classification always returning nothing
(model-free 3.2946 × 10⁻³, red), reading `POS` as 0-based (3.2197 × 10⁻³, red), stripping `chr`
from contig names (panics), and — with no code change at all — running HG002's reads against
HG003's benchmark (2.4881 × 10⁻³, red).

## Not covered

`FILTER` is read by nothing. The whole HG002 benchmark holds 34,211 non-`PASS` records against
4,013,574 `PASS`, about 8 in 1,000, and none of them land in these 100 spans — checked. A
different BED would admit them silently.

The contested-position rule — two records touching one position, resolved to heterozygous —
never fires on either arm. It is reachable and was reached by a mutation; the benchmark simply
does not produce it here.

## Validation

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --lib --bins --tests --all-features` (3,236 passed, 0 failed, 10 ignored — the tenth
is this anchor) and `cargo doc --no-deps --lib` at the 12-unresolved-link pre-existing baseline,
none in this module.
