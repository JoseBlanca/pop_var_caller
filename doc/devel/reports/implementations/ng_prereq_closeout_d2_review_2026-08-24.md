# ng — the review the error-rate scale's denominator never had

**2026-08-24**, branch `ng-prereq-closeout`. Commit `021505f1` — the per-read-group numerator and
denominator of the read likelihood's error-rate scale — shipped with **no independent review at
all**: seven review agents died on API overload across that session, so its author's own three
hand-checks were the whole of it
([`ng_calling_prerequisites_d2_2026-08-24.md`](ng_calling_prerequisites_d2_2026-08-24.md) §4).

Four reviews ran here, each in its own worktree, against the merged commit: **the site-set claim**,
**the fixed-point arithmetic**, **whether the tests could fail**, and **whether the prose is true**.
This report says what they found and what was done about it.

---

## 1. The four findings that were defects, in the order they cost

### 1.1 The depth cap *does* divide the two site sets — 2.7% at 300× (owner's decision, open)

**The claim under review.** §3.2 requires the scale's denominator to average over exactly the reads
its numerator was fitted from. The commit argued the per-position depth cap does not break that,
because the histogram's thinning to 124 reads is hypergeometric on counts and never looks at a
read's quality.

**That argument is right per site and is not the whole of it.** The cap removes reads only from
*deep* positions, so it changes how much weight each position carries in the average: a 500-read
position casts 500 votes in the denominator and 124 in the population the numerator came from. Deep
positions are not a random sample of the genome — they are where reads pile up from elsewhere and
where mapping quality collapses.

**Measured on real reads** (`examples/ng_minted_error_means.rs`). On HG002's 100 benchmark regions
at 300×, where the fit sees 70,288,390 of 172,616,054 read-positions:

| | geometric mean of the minted error |
|---|---|
| as the accumulator counts it | 2.9055 × 10⁻⁴ |
| with each position thinned to the cap first | 2.9862 × 10⁻⁴ |
| ratio | **0.9730 — 2.7%, 0.12 Phred** |

On the tomato cohort — 63 accessions from 2.5× to 28.6×, median 10.7× — it is nothing: measured on
the deepest accession, 228,468,065 of 228,492,796 read-positions are under the cap and the mean
moves by a factor of 1.0000.

**Done:** §3.2 and the module doc now state the divergence with its size, in place of the sentence
that denied it. **And decided (owner, 2026-08-24): the accumulator does not thin, and the average
stays over every read** — because that is the population the scale is applied to at calling time.
The 2.7% is carried knowingly; §3.2 records both options and why this one. Revisiting costs one
multiply per site.

### 1.2 An `i64` running sum would have saturated on one ordinary deep sample, not on four hundred

The commit widened the sum from `i64` to `i128` and justified it at "four hundred human-scale
samples in one read group", reasoning from a billion **reads**.

**What the counter holds is a read at a *position*.** An observation contributes its `num_obs` at
every locus it appears at, so a sample's total is its covered length times its depth — measured,
172,616,054 over 571,984 bases on HG002 at 300×, which is 301.8 a base. That is about 150× more
than the paragraph assumed:

- a human genome at 30× is 9.3 × 10¹⁰ read-positions, 7.9 × 10¹⁷ scaled units — an `i64` holds
  **twelve** such samples;
- **the same genome at 300× is 7.9 × 10¹⁸ on its own — 86% of `i64::MAX`, for one sample.**

So the widening fixed a defect an ordinary run reaches, not a hypothetical one. **Done:** the type's
doc, and the previous report's paragraph, both corrected. **And it is now a test** —
`the_sum_is_wide_enough_for_a_run_that_an_i64_would_have_saturated` accumulates 1.288 × 10¹²
read-positions at 8 nats each and asserts the mean is exactly −8; an `i64` pins at `i64::MIN` and
returns −6.827. Swapping the field back to `i64` had left all 23 tests green.

### 1.3 Four tests could not have failed, and one of them guarded the step's own headline property

| test | why it could not fail | fixed by |
|---|---|---|
| `each_read_groups_mean_is_over_its_own_reads`'s ordering assertion | its fixture's read groups already arrived ascending, so the sort it names is a no-op on it — deleting the sort left all 23 tests green | the higher-numbered group is written first |
| `only_generic_loci_are_counted` | visits the empty case first, so a build that emptied the scratch buffer *after* the kind gate — leaving the previous locus's totals standing — passed | the generic locus goes first |
| `only_generic_loci_are_counted`, again | `LocusKind` has **three** arms and it names two; a gate rewritten as "everything but a bundle" poured every repeat-tract read into the denominator and stayed green | new test `a_repeat_tract_contributes_nothing` |
| nothing, for a supplied inbreeding coefficient | every accumulator fixture builds through a helper that hard-codes `Fitted`; moving the fold into the `Fitted` arm left **274** tests green while a supplied-`F` run shipped an empty denominator | two assertions on `a_supplied_inbreeding_coefficient_collapses_the_windowed_table` |

Each fix was run red against its mutation and green against the real code, both directions checked.

**And one guarding test earned its place.** The commit changed a shared fixture's `q_sum` from
`0.0` to distinguishable decimals so that the six-order merge test could see an order-dependent
fold. Two reviewers reconstructed the counterfactual independently: with a genuinely `f64` running
sum in place *and* the old zero fixture, `three_shards_merged_in_every_order_give_the_same_counters`
passes; restore the awkward fixture and it fails on the last bit. The account given for that fixture
change is confirmed by measurement.

### 1.4 A sentence in the walk still asserted the decision that had been overturned

`minted_ln_read_error`'s own doc — the function the new module links to, and the first place a
reader following that link lands — still said *"the scale wants the arithmetic one"* and *"the
accumulator needs its own sum of `exp` of this answer, minted where this is"*. Both were false the
moment the owner settled the question, and the commit did not touch that file. It also still said
the accumulator "does not exist yet". **Done:** rewritten, with the measured ratio and the
mate-overlap mechanism behind it.

## 2. Corrections to prose, with what each would have cost a reader

Eleven claims were checked against the code they cite and found wrong. The ones that would have
misled:

- **"an observation no read is behind carries a `q_sum` of exactly zero, which converts exactly"** —
  the rounding bound rests on this, and it is not a fact about the fold. The fold never builds a
  read-less observation at all (both column paths create one only when a read arrives), which is a
  *stronger* reason; but the module's sibling fixture does build one, with a non-zero `q_sum`. Now
  stated as it is, and pinned: a read-less observation carrying −1.5 moves the sum and not the
  count, and the test says so rather than a guard hiding it.
- **`fold_into`'s doc justified its sort by `f64` non-associativity** — in a function whose sum is
  an exact integer. Measured: deleting the sort changes no total. The sort is kept because the
  function *states* it returns ascending order and a caller reading the scratch vector directly
  should get it; the doc now says which of those two things it is for.
- **"the eight bytes are free"** — `{i128, u64}` is 32 bytes against `{i64, u64}`'s 16, because
  `i128` aligns the struct to 16. The field widens by 8 and the struct by 16.
- **"§3.3's identity is stated to a tolerance rather than bitwise"** — §3.3's *is* bitwise; it is
  §3.6's that is not. Wrong section cited.
- **"the six-order test above"** in `merge` — that test is below `merge`, not above. Now named
  instead of pointed at.
- **`posterior_engine.rs:2357`** for production filling its mixture table — 2357 is the
  genotype-shape-cache comment; the mixture branch starts at 2361.
- **"four fixed-point conversions' worth of rounding"** — eight, four in each of two shards. The
  miss is 1.226 × 10⁻⁹ against a 4.768 × 10⁻⁷ bar, 389 times inside it, both confirmed.
- **The plan's test table claimed a test that does not exist** for a property that does not hold —
  "the totals equal the histogram's own read count". They cannot: the histogram thins and the
  accumulator does not (§1.1). Replaced with what is actually checked, and on what.
- **§3.2 promised "§12 gains a test that reports both per read group on a real cohort"** and §12
  gained nothing. §12 now carries the measurement as a *change measurement*, done, with its result.
- **The previous report's `cargo doc` baseline of 24** — it is 23. A wrong baseline is exactly what
  lets a twenty-fourth in.
- **The previous report said "three of four review agents died"** in one paragraph and "seven" in
  the next, and "three of four relaunches" implies one survived, contradicting "no independent
  review" two lines below. Corrected to seven.

## 3. What the reviews checked and did not fault

**The site sets themselves are identical, and this was attacked rather than confirmed.** Both paths
run behind one `LocusKind::Generic` gate in `add_locus`, before its inbreeding-mode branch; both
iterate `complete_observations()`, the same iterator and not a lookalike; both count `num_obs`.
Measured across four run modes — one library and two, `F` fitted and `F` supplied — the totals are
identical in all four. The ploidy grain is right: the histogram keys `(read group, ploidy)`, the fit
gathers a group's ploidies into one scan and returns one rate per group, and the fit skips no cell —
a fixture with one haploid site against fifty diploid ones came back with all 51 counted and no
threshold to trip.

**And on real reads**: the reads the accumulator counts and the reads the walk emitted at those loci
agree exactly — 172,616,054 on HG002 at 300×, and on all 63 tomato accessions, with no locus left
unruled on in any of the 64 runs.

**Two numeric worries were chased and are not defects.** `i128 as f64` above 2⁵³ does lose bits —
584 scaled units on a 1.18 × 10¹⁹ sum — but the loss is *relative* to a sum the mean then divides
by the read count, so the worst mean miss across every magnitude from 2¹⁰ to 2⁶² scaled units is
3.55 × 10⁻¹⁵, which is 7.5 × 10⁻⁹ of the bound. And `x / 2²⁰ / reads` against
`x / (2²⁰ · reads)`: zero differences in seven million probed pairs, because 2²⁰ is a power of two.

**One correction to the bound's own wording.** It is *at most* 2⁻²¹ on the mean, not *under*: a
`q_sum` sitting exactly half a unit off the grid attains it at one observation and still attains it
at twenty million. The doc said "stays under"; a second assertion now pins the attained case.

## 4. One test gap left open on purpose

`the_fixed_point_rounding_does_not_accumulate_in_the_mean` wrote its tolerance as
`1.0 / PARTS_OF_ONE / 2.0` — which moves with the constant it is checking, so coarsening the grid
from 2⁻²⁰ to 2⁻¹⁹ left it green. **Fixed**: the tolerance is the literal `4.768…e-7`, and
coarsening now reddens it. That is the only mutation of this class that survived and it no longer
does.

## 5. Mutation totals, across the four reviews

**44 mutations run, 30 killed, 11 survived, 3 changed no behaviour at all** — 10 from the site-set
review, 10 from the arithmetic one, 21 from the one whose whole subject was the tests, and 3 from
the fact-check of the prose.

The eleven survivors are **eight distinct mutations**, three of which two reviews found
independently. **Six of the eight are closed by the changes made here** — the `Fitted`-arm fold, the
moved `out.clear()`, the "everything but a bundle" gate, the deleted sort, the `i64` swap, and the
coarsened fixed-point grid. Each was re-run on this branch against the test written for it and goes
red.

**Two are left open, and neither is a hidden defect.** Routing the running sum through `f64` is
exact below 2⁵³ scaled units, so no fixture can distinguish it; it becomes real one human sample
later, which is where the `i128` test above now stands. And swapping `saturating_add` for
`wrapping_add` changes nothing reachable: at `i128` saturation needs 1.8 × 10³⁰ read-positions. Both
are named in §8 rather than papered over.

Every survivor was proved to change behaviour by a program printing two different outputs on one
input, not by argument. The three behaviour-neutral ones: reordering the fold above the read-group
histogram call in `add_locus` (disjoint state), reversing a locus's read-group slice before
`fold_into` (integer addition commutes over distinct keys), and deleting the sort as far as the
*totals* are concerned — that one changes the returned vector's order, which is why it is a
survivor for the public function and neutral for the accumulator.

## 6. The reversion test, answered plainly

**Reverting only the integration** — the fold in `add_locus`, the merge arm, the accessor's
contents, leaving the type and its own tests — reddens **exactly one test of 23**,
`merge_carries_every_counter_and_both_tables`. Everything else stays green, including the six-order
merge test, whose new assertion then compares one empty map against another. That single by-value
assertion is the whole integration guard, which is what its own comment claims.

**Reverting the whole commit** — file, fold, merge arm, accessor and its tests together — leaves
the rest of the suite green. **That is the state of play and not a defect**: nothing consumes
`minted_errors()` anywhere in `src/`, `examples/`, `tests/` or `benches/`. The consumer is the read
likelihood's error-rate scale, in
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md) A2.

## 7. Validation

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo clippy --example ng_minted_error_means` | clean |
| `cargo test --lib ng::parameter_estimation::generic::calibration` | **9 passed** (was 7) |
| `cargo test --lib ng::parameter_estimation::generic::accumulators` | **16 passed** |
| `cargo test --lib ng::locus_generation::pileup::minted_error_census` | **4 passed** |
| `cargo test --lib` | **4,187 passed, 0 failed, 14 ignored**, 925 s — `main`'s 4,181 plus the six tests added across both commits |
| `cargo doc --no-deps` | 23 unresolved links, 12 redundant-explicit-link-target warnings — `main`'s baseline, unchanged |

## 8. What the reviews leave open

- ~~**The depth-cap decision** (§1.1)~~ — **settled 2026-08-24: the accumulator does not thin.**
- **A read group with a borrowed error rate.** Below 10,000 sites a group gets the mean of the
  other groups' rates rather than its own, while its denominator stays its own reads. §3.2's
  sentence about one site set does not describe that case; the module doc now names it, and whether
  it is right is a sentence the owner owes the spec. A capture panel, or a minor library in a
  multi-library sample, reaches it.
- **`saturating_add` on both fields is untested in both directions** — swapping it for
  `wrapping_add` left all 23 tests green. At `i128` neither is reachable (saturation needs
  1.8 × 10³⁰ read-positions), so this is a question about what a future narrowing should do, not
  about this build. Left as it is, and named.
- **`fold_into` has no direct test** — it is exercised only through `GenericAccumulators`.
- **No fixture in either module runs a locus the walk produced.** Every one is a hand-built
  `SequenceObservation`. What covers the walk is `ng_minted_error_means`, on real reads, and it is
  not a test.
