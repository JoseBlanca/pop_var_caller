# ng step 4, the STR path — C1: which stratum a locus belongs to

*Implementation report, 2026-08-11. Step C1 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — two agents, 19 mutations, 4 behaviour-changing survivors. Design
authority: [`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.3, §4 and
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1.*

## What the step is

`stratum_of`: the reference tract's length over its motif's period. **A pure function of the
reference**, which is the property the whole design rests on — every sample files a tract under
the same stratum whatever its own reads showed, so one sample's stutter can be compared with
another's and a cohort can pool them. A sample whose alleles differ from the reference does not
move between strata; its reads land at an offset instead.

## Recorded deviations from the architecture

1. **Three answers where arch §2.3 sketches `Option<Stratum>`.** The two ways a locus can have no
   stratum mean opposite things to the accumulator that asks: a locus that is not one repeat tract
   is passed over in silence, while a tract whose reference length is not a whole number of motif
   copies is a delimiting fault upstream — counted, reported, never rounded. An `Option` collapses
   them, and C5's counter would have to re-derive from the locus what this function had already
   decided. The reviewer checked this against the arch's own `add_locus` contract ("ignored" versus
   "counted and skipped") and agreed; it also confirmed that an error type would be wrong here,
   since a `?` at the C5 call site would abort a cohort run on one mis-delimited locus.
2. **A panic rather than a saturating cast** on a repeat count past `u32::MAX`. Saturating would
   file the locus under a stratum it does not belong to and say nothing, which is what the type
   exists to prevent; the panic is unreachable twice over (the tract would be a 4.29-billion-base
   allocation, and the catalog refuses to serve a tract past 500 bases).

## What the review changed

**Blocker — the divisor was unpinned across two-thirds of the STR scope.** Every test asserting a
successful stratification used period 1 or 2; periods 3 and 4 appeared only where the assertion
was the *fault* arm, whose reported period comes from the motif and never from the divisor.
Clamping the divisor at 2 or at 3 left all seven tests green — and a sixteen-base tetranucleotide
tract then files at eight copies instead of four, or is refused outright. A table over all six
periods and a property test now cover it; I reproduced the clamp and it fails four tests.

**Major — every non-divisible fixture left a remainder of exactly one.** Narrowing the guard to
"one base past a whole number of copies" survived the suite, silently flooring a fourteen-base
trinucleotide tract to four copies — exactly the rounding the third variant exists to prevent,
with the counter that should catch it reading zero. The table now carries remainders 2, 4 and 5,
and the property test carries all of them.

**Major — the purity test only ever showed reads *shorter* than the reference.** It is the one
test that can see the function read the observations at all, and it uniquely kills a rule taking
the length from the first observation — but a rule taking the *longest* of the reference and the
observations passed everything. A sample carrying an expansion would then migrate to a different
stratum from every other sample at that tract, which is precisely what the test is named after. It
now carries a two-copy expansion.

**Minor — a fourth `LocusKind` would have been passed over in silence.** The `let … else` is now a
kind-by-kind `match`, so a locus kind added later is a compile error rather than a default.

**Two wrong claims of mine, both about wording rather than arithmetic.** The doc quoted the
architecture as saying the counter "should read near zero" — the arch says "should be near zero",
and the string as quoted is from a sibling source file. And the panic note's bound named "a contig
longer than any genome carries", which is about seven orders of magnitude too loose and holds only
for human anyway: a tract never approaches contig length, because the catalog refuses to serve one
past `CATALOG_MAX_STR_LEN_BP` — 500 bases, against step 3's calling default of 100. Both corrected.
Every other figure held, including the copy floors (8 at period 1, 6 at period 2) against spec
§5.1.1 and `segment_criteria.rs`, and the catalog's own looser floors quoted where the sentence was
about what the catalog emits.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib parameter_estimation::ssr` | **80 passed** (70 before this step) |
| `cargo test --lib --bins --tests --all-features` | **3,457 passed, 0 failed, 10 ignored** |

Counted rather than recalled: `grep -c '#\[test\]'` on `locus_offsets.rs` gives **10**, and the
suite moved 3,447 → 3,457. I re-ran three of the reviewer's mutants against the fixed suite — the
clamped divisor, the remainder-one guard and the longest-observation length rule — and each now
fails between one and four tests.

**Two gates are red on this branch and neither is this step's**: `cargo clippy --all-targets` fails
in four `examples/` files, and `cargo doc` reports 13 unresolved intra-doc links.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-c1/` — two per-category files (reliability; errors, naming
and the numbers check) and the reviewed patch.
