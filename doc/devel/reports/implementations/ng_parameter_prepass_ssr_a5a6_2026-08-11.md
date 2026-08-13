# ng step 4, the STR path — A5+A6: what a fit returns, and how it fails

*Implementation report, 2026-08-11. Steps A5 and A6 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — three agents, 19 mutations, 7 survivors of which 4 changed
behaviour. Design authority: [`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md)
§2.4, §4.1–§4.3 and [`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md)
§4.3–§4.6.*

## What the step is

Types only: what a stutter fit returns, the summary a person reads instead of several hundred of
them, the accumulator's counters, this path's error enum and its four constants. Nothing produces
any of it until Milestone E.

**The output has two halves because this path runs one fit per (read group, stratum)** — several
hundred a sample, against the SNP/indel path's four fits in total. A per-stratum record exists so a
wrong fit is traceable; the summary exists because several hundred records are a file nobody opens,
and a flag nobody reads is how a badly-fitted parameter reaches a caller.

**`StratumFit` carries two provenance lists**, which is the least obvious thing here. The level and
the two shares that describe slippage are measured from different populations, and at the bottom of
the repeat range those differ by four orders of magnitude: a stratum of 100,000 loci at five reads
each holds half a million reads, of which about 455 slipped, 77 of those gained a repeat, and 5 of
those gained two. So the same stratum measures its level to about 5% of itself and its fall-off to
about 45% — and the honest answer is to keep the level it measured and borrow the two shares it did
not, with each half saying where it came from.

## Recorded deviations from the architecture

1. **`AllelePairFrequency` instead of `Vec<(WholeRepeatOffset, WholeRepeatOffset, f64)>`** — named
   fields, for the reason `fitting/`'s `WeightedCell` replaced a three-member tuple. The review
   asked for more: its constructor now **sorts the two alleles**, so one genotype has one spelling
   and a fit cannot emit `(-1, 0)` and `(0, -1)` as two entries whose frequencies then fail to sum
   to one.
2. **`SlippageStart::from`, not `start`** — `start.start` asks a reader which of the two is the
   noun; the sibling `StartOutcome` in `fitting/` already uses `from`.
3. **Two named fields instead of `loci_behind_fits: (u64, u64)`** — nothing in a pair of `u64`s says
   which end is which, and this is the tuple `AllelePairFrequency` was created to refuse.
4. **`Domain` is neither `#[from]` nor transparent**, where the architecture sketches
   `#[error(transparent)] Domain(#[from] DomainError)`. It matches the generic path's twin instead:
   a transparent variant forwards the inner message and drops the sample, and `#[from]` would let
   `?` mint the variant at every constructor, so the context would have to be *remembered* rather
   than demanded.

## What the review changed

**Major — no failure named the read group.** The fits are keyed on `(ReadGroupId, Stratum)` — this
module's own opening point — and all three variants carried only the sample. On a four-library
sample, a search that did not settle named one of four fits without saying which. All three now
carry the read group, and `Domain` carries the stratum too; `ReadGroupId` gained a `Display` so a
message can render it as the bare index it is.

**Major — `NoFittableStratumAtPeriod` named the floor but never the shortfall.** A reader could not
tell a period whose thickest stratum held 812 loci from one that held 3, which is the difference
between widening the run and dropping the period. It now names both. Its old closing advice, "supply
the parameters", pointed at a door that does not exist — no config field takes them — and is gone.

**Major — `SsrSampleParameters` derived `Default`.** On this unit's *output* type, where the sibling
`GenericSampleParameters` deliberately does not: `.unwrap_or_default()` would report a never-fitted
sample as an empty one. Removed.

**Four mutation survivors that changed behaviour and no test noticed:** the two floors could move
(the messages interpolate them, so a floor raised tenfold changes what a caller is told, silently);
`{starts}` could be dropped from the unsettled-search message, so "across 4 starting points" became
"across its starting points" — two starts landing 333-fold apart is a different claim from four
doing so; and a default summary's `loci_behind_*` could be seeded at `u64::MAX`, which a `min`-fold
can never leave. Four new tests close them.

**Three tests asserted only over struct literals the test itself wrote**, so no implementation could
fail them. Kept, renamed and **labelled as shape tests** — what they can catch is the *type* losing
its ability to carry two different answers; what they cannot catch is a producer filling both lists
alike, and that assertion is recorded for the milestone where a producer exists.

**Two wrong numbers, both mine, both inherited from the design.** "A stratum can clear
`MIN_LOCI_TO_FIT` by four orders of magnitude and still put five reads behind the fall-off's gaining
arm" — the worked example is 100,000 loci against a 1,000 floor, which is **two** orders; at four
the gaining arm holds about 503 reads, so the sentence contradicted itself. And "the level and the
two shares starve at rates 20,000 apart" restates a ratio whose point is a **change of unit** —
100,000 loci against 5 gaining-arm reads — which the restatement dropped. Both corrected in the
code; **both are still wrong in `spec/parameter_prepass_ssr.md` §4.5 and
`arch/parameter_prepass_ssr.md` §2.4**, which is the owner's to fix.

**Declined:** moving the two floors into a config (they are unmeasured, and a knob invites tuning
what has not been measured — recorded in `MIN_LOCI_TO_FIT`'s doc, which now says plainly that it
carries no derivation and none exists upstream); renaming `fitted_over`.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --all-features` | **3,410 passed, 0 failed, 10 ignored** |
| `cargo test --lib parameter_estimation::ssr` | 33 passed |

Counted rather than recalled: `grep -c '#\[test\]'` on `ssr/mod.rs` gives **19**, and the suite moved
3,399 → 3,410 with this step.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-a5a6/` — three per-category files and the reviewed patch.
