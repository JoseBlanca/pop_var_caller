# ng read likelihoods — C1: the contamination mixture, and a ceiling the model cannot carry

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step C1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone C, on
top of `17a910b5`'s parent.*

## 1. What it is

Some of a sample's reads came from somebody else's DNA. The row now says so:

```text
log Lg(g)  =   Σ  n_o · log[ (1 − c_{r(o)}) · own(o | g)  +  c_{r(o)} · q(o) ]   +   q_sum_other
                o
```

`own(o | g)` is what the row computed before — `k_a/P` for a read the genotype explains, `ε̄/m`
for one it does not. Beside it sits the chance the read is not this individual's at all: `c` is
the share of this read group's reads that came from somebody else, and `q(o)` is how often that
somebody else shows the allele this observation shows. `q(o)` arrives as a parameter; computing
it from the batch the sample was sequenced in is C2.

**Spec §3.6's mixture replaces §3.3's closed form as the row's one shipped path.** There is no
second code path, and no `c == 0` branch.

## 2. The claim this step had to keep rather than assert

**With no mixture, the row computes spec §3.3.** That is what lets contamination default on: a
clean cohort is untouched by the default. It is also the kind of claim that is easy to write and
easy to be wrong about, so it is measured.

- **The explained side is bit for bit.** The mixture multiplies the copy share by `1 − 0` and
  adds `0 · q`, and both are exact in floating point, so `n·ln(k/P)` is the same bits it was.
- **The error side is not**, and cannot be. §3.3 charges `q_sum + n·(log scale − log m)` in log
  space; §8 requires the mixture be evaluated in probability space and logged once, so the two
  forms are separated by an `exp`/`log` round trip.

Measured across **3,552 comparisons** — three loci of two, three and four alleles, four quality
profiles, read counts from 1 to 300, three read-group scales, ploidies 2 and 4, every genotype —
the worst relative disagreement is **2.9 × 10⁻¹⁶**, at three alleles and 300 reads: −1592.6682015614529
from the mixture against −1592.6682015614524 from an independent log-space oracle. In the units a
genotype is decided in that is **2 × 10⁻¹² Phred**. The test's named tolerance is 3.5 times the
measurement.

**The oracle is written out again rather than being the row under a flag** — a shared
implementation would agree with itself whatever either of them did — and it reads its spreads
from `fill_log_error_spreads`, so deleting that filler is a compile error rather than a quiet
subtraction.

### What the sweep actually catches

The plan asks that this test fail if anyone reintroduces production's extra `(1 − ε)` factor or
its allele-count divisor into `own`. Both were injected and measured:

| injected defect | comparisons disagreeing | worst relative |
|---|---|---|
| production's `(1 − ε)` factor on explained reads | 3,172 of 3,552 | 0.014 |
| production's allele-count divisor on the error mass | 2,336 of 3,552 | 0.19 |

The first also failed eleven other tests in the module, the second two. Both files were restored
from a byte-identical copy and the checksum re-checked afterwards.

## 3. The finding: A2 built the wrong floor pair, and this step is where it shows

`ReadGroupCalibration::calibrated_error` was built at A2 for this consumer and clamps into
`[MIN_BASE_ERROR, MAX_BASE_ERROR]`. **The ceiling cannot be applied to what the row charges.**

The reason is not a preference. Spec §2.3 requires that no term be a non-linear function of a
per-read quality, because the merge hands the row a *fold* of reads and `q_sum` recovers only
their geometric mean. `min(x, ½)` is exactly such a function, and the row's own aggregation
fixture is where it bites: a read at Phred 1 is minted at an error of 0.794, so the ceiling
charges it 0.46 nats less — and it binds on every such read taken singly and on **none** of the
folds, whose geometric mean over the 93/1 alternation is 2 × 10⁻⁵. At that fixture's 300-read
end, 150 such reads, capping would move the answer by **69 nats**, against a property pinned to a
relative 2 × 10⁻¹⁴.

So the row charges through a new `charged_error`, which floors and does not cap.

**The floor is different in kind, which is why it stays.** It is not reached by any quality the
read preparation admits: Phred 93 is an error of 5.0 × 10⁻¹⁰, and at 0.37 — the smallest scale
the fixtures use — 1.9 × 10⁻¹⁰, **185 times the floor**. So it changes no answer that can occur
and exists only so a logarithm never sees a zero (spec §8).

**And this is the conservative reading, not a new decision.** §3.3's log-space form charges
`q_sum + n·log scale`, which exceeds zero under exactly the same conditions a capped error would
have clipped. Introducing the cap at C1 would have been the change; not introducing it preserves
B2's numbers exactly.

**Two things are left open for the owner and are not taken here** — §7.

## 4. Two departures from the architecture's sketch

**The row takes `ContaminationMixture<'_>`, not `&[ContaminationView]`.** The two halves of the
mixture sit in different tiers — the fraction frozen before the loop, the frequency moving with
it — and holding them in one type is what lets a single construction check that they describe the
same run, instead of the row rediscovering it per observation. It also gives the uncontaminated
case a name (`uncontaminated()`) rather than an empty slice. This is the same latitude arch §3
already took for `LogErrorSpreadTable`.

**The row reads `spreads_of`, which exponentiates.** The mixture needs `m`, not `log m`, so B2's
storage decision now costs an `exp` per error-side term. The argument that decided it — *spec
§3.3's closed form takes no logarithm inside its loop* — is precisely what C1 voids, since §8 puts
one there by specification. The change is one accessor wide and is raised rather than taken (§7).

## 5. What the reviews changed

Four review agents ran on the committed step, each in its own worktree, covering reliability,
errors, naming, idiomatic, smells, refactor safety and defaults. What they found, and what it
changed:

**Two silently-wrong-answer hazards, both fixed.**

- **`f64::max` swallows a `NaN`.** `x.max(MIN_BASE_ERROR)` returns the floor when `x` is `NaN` —
  so a `NaN` `q_sum` would have come back as **the most confident error probability the module
  admits**, and the row finite, plausible and wrong. Its sibling's `f64::clamp` propagates the
  `NaN` instead, which is why the difference was invisible: the doc paragraph explaining the
  hazard had been copied from the sibling and was true only there. `charged_error` now asserts
  that `q_sum` is finite and at or below zero, and that the scale is finite and positive —
  `scale` is a public field, so `from_fitted_rate`'s own guard is bypassable. A `q_sum` above
  zero was measured returning a *positive* log-likelihood of +48.90.
- **The row checked the spread table's stride and not its height.** The inner loop walks a
  column by striding and `zip` stops at the shorter of the two, so a table filled at another
  ploidy truncated the walk in silence: a tetraploid handed a diploid's table scored its last
  genotype `0.0` against its own `−0.863`, **and that genotype won**. This hole predates C1 —
  B2 added the stride check without it — and is fixed here, with a `genotype_count()` accessor
  to check against.

**One check that was eager on one half and lazy on the other.** The mixture's allele count was
checked against the locus before the walk; its read-group count was checked only when some
observation happened to name a group past the end. So a locus whose reads all came from the
first few groups passed, and a mismatch surfaced at whichever locus first reached further, or
never. Both halves are now checked up front, which makes both accessors' own panics unreachable
in a run — so they are exercised directly instead.

**One test was pinning a property on the function the row had stopped using.** Spec §12's tenth
test — the calibration scale reproduces the fitted rate — ran through the capped reading. §3.2's
property is *the average charged error equals the measured rate*, and the cap is what breaks it,
so the test was making a claim only the other function can keep. Its fixture sits well inside the
cap, so it passed either way; that is exactly why it had to be moved rather than left to fail one
day.

**Two pairs of look-alike methods, renamed and narrowed.** `charged_error` /
`calibrated_error` and `spreads_of` / `column_of` each had identical signatures and different
units, and a probe substituting one for the other compiled with no diagnostic. Reading `log m`
where `m` is meant charges 2.73 times too much where the spread applies and *divides by zero*
where it does not, since `NO_LOG_ERROR_SPREAD` is `0.0` and the linear form is `1.0` — giving a
log-likelihood of `+inf` rather than a panic. `column_of` is now `log_spreads_of` and
`pub(crate)`; `calibrated_error` is now `charged_error_capped_at_half` and `pub(crate)`.

**Six naming findings on the new code**, all applied: `contaminated` named a fraction and read as
a predicate; `from_somebody_else` was presented as the complement of `from_this_individual` and
is `c·q(o)` rather than `c`; `wrong_and_this_individuals` had no noun and held an intermediate;
`frequency_of` did not say *whose* frequency where two allele frequencies are in play at one
locus; `log_explained` was called `log_mixture` five lines from where it was read; and
`none()` borrowed `Option`'s word for a case the docs call *uncontaminated*.

**One route to a hidden default, closed.** `ContaminationMixture::new(&[], &[])` compiled and
produced a value indistinguishable from the named constructor, so a caller building both halves
from a fit that returned nothing could land on the clean case through a constructor whose name
says the opposite. It is now refused, with a message naming the constructor that means it.

**Two things the reviews checked and found sound, worth recording** so they are not re-checked:
every pair of the row's six parameters is a distinct type, so no transposition compiles — rustc
even offers `help: swap these arguments`; and `ContaminationMixture` has no `Default` impl, which
is the right call for a type whose empty value is a modelling claim.

## 6. What the tests pin

108 tests in the module, against 81 at B2 — 27 added. Beyond the sweep and the injected defects
of §2:

- **A hand-computed contaminated case**, every term written out: a diploid, two reference reads
  at a summed log error of −6, one alternative read at −7, 3% contaminated, the contaminant
  carrying the alternative at 1 in 1,000.
- **What the mixture does to a call, and what it depends on.** The heterozygote's lead over the
  reference homozygote falls from **6.019 nats** with no mixture to 5.981 at a contaminant
  frequency of 1 in 1,000, 5.376 at 1 in 100, and **2.131 at 1 in 2**. So a rare contaminant
  allele buys almost nothing and a common one buys 3.9 nats, 17 Phred. That is the same lever
  §3.6 measures for what the mixture *costs* the aggregation contract — 0.14 nats at 1 in 1,000
  against 1.89 at 1 in 2 — pointing the same way: the fixture where contamination changes a call
  most is the fixture where pooling reads costs most.
- **Each read group is charged its own fraction.** Three alternative reads from a 5%-contaminated
  library against the same three from a clean one: the reference homozygote finds them 8.57 nats
  — 37 Phred — less surprising in the contaminated library, and the heterozygote moves 0.12 nats
  the other way. A row that averaged the two fractions, or read the first group's for every
  observation, passes every other test in the file.
- **An explicit zero fraction and no mixture at all give bit-identical rows.**
- **The linear column is the bases the log column is the logarithm of** — `exp(0)` exactly `1.0`,
  and the other within a unit in the last place of `ERROR_SPREAD_BASES`.

## 7. Open for the owner — two questions, neither blocking C2

**Should `charged_error_capped_at_half` and `MAX_BASE_ERROR` survive at all?** After C1 neither
has a consumer outside this module's tests. The argument for keeping them is that they record
what production does and what we deliberately do not; the argument against is that a `pub(crate)`
method with a byte-identical signature to the one that *is* correct is a hazard whose only
defence is its name. **Recommendation: delete both**, and move the reasoning into
`MIN_BASE_ERROR`'s doc, which is where a reader asking "why is there no ceiling?" will look. Not
done here because it retires an A2 decision, and retiring one quietly inside a C1 commit is
exactly the move this plan's reviews keep catching.

**Should the error-spread table store `m` rather than `log m`?** B2 decided `log m` with a stated
argument — §3.3's closed form takes no logarithm inside its loop — and C1 is what makes that
argument false, since §8 puts one there by specification. Storing `m` would be exact and free and
would remove an `exp` per error-side term; the log form would then have no consumer at all.
**Recommendation: switch it**, as its own small commit so the change is visible, before D1 adds a
second reader. The whole change is `fill_log_error_spreads` writing `1.0`/`3.0`, the type losing
its `Log` prefix, and the accessor pair collapsing to one.

**And one thing the reviews found that belongs to the plan, not to the code.** Spec §3.6 requires
that *the run's output must still carry the fraction used, per sample*, because a genotype
computed at `c = 0.03` and one at `c = 0` are otherwise indistinguishable. **No step of this plan
owns that**, and `ContaminationMixture` currently offers nothing to read it from — `fraction_of`
returns the same `0.0` for *absent* as for *measured clean*, which is the distinction the rest of
this module is careful about (`ContaminationView::of_estimate` returns an `Option` for exactly
it). It needs a home before Milestone C closes.

## 8. Validation

All in the container, on the committed tree:

- `cargo fmt --all -- --check`: clean.
- `cargo clippy --lib --all-features --tests -- -D warnings`: clean. **The repo-wide
  `--all-targets --all-features` run is red on `main`**, in `examples/ng_duplicated_class_harness.rs`
  and `benches/freebayes_bookkeeping.rs` — pre-existing and out of this step's scope.
- `cargo test`: **4,298 passed, 0 failed, 14 ignored**; 108 of them in
  `ng::calling::likelihood`.

One build fix rode in ahead of this step and is on `main` as well: `cargo test` did not compile
at all, because `examples/ng_ssr_thin_stratum_gate.rs` names a module gated behind the
`bench-fixtures` feature and nothing in `Cargo.toml` said so.
