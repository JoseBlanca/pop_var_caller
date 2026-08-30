# ng calling loop — E2c: the repeat tract's scoring parameters

**Step:** E2c of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the assembly between the
run's fitted parameters and the STR row. **The step and its plan entry were both created by this
work**, on the owner's ruling of 2026-08-26.
**Design authority:** [`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §4.2, §4.3,
§4.4, §4.5; [`arch/read_likelihoods.md`](../../ng/arch/read_likelihoods.md) §4.1, §4.2.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

The repeat-tract row has been shipped and merged since `calling_read_likelihoods.md`'s H2, and
nothing outside its own tests could build what it takes. `SsrScoringContext::new` had no production
caller, `fill_reachable_lengths` had none, and the outlier weight had no source anywhere. A new
module, [`calling/inference/repeat_tract_parameters.rs`](../../../../src/ng/calling/inference/repeat_tract_parameters.rs),
now reads the run's fitted slippage and substitution numbers into one scoring context per
`(read group, candidate)`, builds the length support the junk term is spread over, and supplies the
outlier weight — so a tract's candidates plus a run's parameters produce an `SsrLocusParameters`
the row accepts. **What it does not do is wire that into the driver**, which is E3's; a repeat tract
is still refused at `call_locus`'s front door, unconditionally, and the refusal's message now names
E3 rather than the step that has already landed.

## 2. The three answers it owes, and what it answers

Two of the three fitted numbers a scoring context carries can be missing on perfectly good data,
and the outlier weight — which is not on a context at all, but beside them on `SsrLocusParameters`
— is not fitted anywhere. None of the three had a rule.

| what is missing | what happens | warrant |
|---|---|---|
| the candidate's stratum is not in the slippage fit | `StutterModel::hipstr_shipped()` | `Defaulted` |
| this read group has no fitted substitution rate at that stratum | `DEFAULT_SSR_SUBSTITUTION_RATE` (0.001) | `Defaulted` |
| the outlier weight, which nothing fits | `DEFAULT_OUTLIER_WEIGHT` (0.01) | *not in any warrant* — §2.3 |

### 2.1 A candidate the fit never reached

`NoSlippage::NoSuchStratum`'s own documentation says a caller has to have an answer, because *"a
candidate several repeats from its reference tract's length can land here on perfectly good data"*,
and `NoSlippage` names four different absences: no such stratum, this read group's slippage group
put no read here, an unknown read group, and a slippage group past the fit's own rows. **All four
take one answer**, because there is one rung below the fit and inventing a second would be a policy
nothing in the three calling documents states.

**But the four are not alike, and the first draft threw the reason away while claiming it was what
told them apart.** Two are ordinary; the other two — an unknown read group, and a slippage group
past the fit's own rows — are what `NoSlippage` itself calls *"the run is not what it claims"*: the
parameters and the reads came from different runs. A single count would let that arrive as routine,
so `cells_whose_read_group_the_fit_does_not_describe` counts it apart. A run should act on that
number rather than record it.

**The alternative considered and not taken** is borrowing a neighbouring stratum's numbers — the
same locus, the same chemistry, one repeat count away, which is biologically closer than HipSTR's
shipped row. It is not taken because *which* neighbour is a rule nobody has written, and
`stratum_fits` deliberately does not blend: a borrow invented here would make ng's slippage numbers
depend on an implementation step's choice. **This is banked for the owner rather than settled.**

**What the constant costs where it binds**, said plainly: `hipstr_shipped()` makes expansion and
contraction equal (0.05 each), where every real fit here is contraction-biased — tomato
dinucleotides sit at 0.83 shorter against 0.17 longer. So a defaulted cell scores a read that lost
a repeat and one that gained one alike, which no measured tract does. The cell says `Defaulted`, and
that is the whole of the protection.

### 2.2 A cell with no fitted substitution rate

The pre-pass emits this rate as `Provenance::FittedHere` or not at all: *"there is no rung below
it for this parameter"*. This module is that rung.

**What reaches it is not the case the first draft named.** The draft said *a stratum that compared
no bases has nowhere to fall* — but that is the case the emitter's own paragraph ends by ruling
out: *"Nothing observed can reach it: a read reaches a table only through a complete witness, which
compares its bases."* What actually reaches this code is the absence
`FrozenParameters::ssr_substitution_rate_at` documents — **the fit has no entry for this
`(read group, candidate stratum)` at all**, ordinarily because the candidate sits several repeats
from its tract's length and its stratum is not in the fit. That is the same ordinary case §2.1
describes, one parameter over, and the module's own test reaches it by deleting map entries.

**The number is the SNP/indel path's default, by definition rather than by coincidence:**
`DEFAULT_SSR_SUBSTITUTION_RATE` is written as `generic::DEFAULT_ERROR_RATE`, so editing that
constant moves this one. That is deliberate at the point where **neither** rate was measured — a
run should not default its two error parameters to two different guesses — and it ties nothing
that spec §4.3 forbids tying: wherever either rate is measured, its own measurement is used.

**The argument for the number is thinner than the first draft made it.** `parameter_prepass_ssr.md`
§4.5 requires the two rates to agree to within a quarter-Phred *where a stratum barely slips*,
which is a statement about low-slippage strata that **were** measured — not about strata that were
not measured at all, which is the condition this constant is reached under. And base quality inside
tracts is systematically worse than outside them (§4.1), so 0.001 is very likely optimistic at a
tract, by an amount nothing here establishes.

### 2.3 Why the outlier weight is in no warrant

Spec §4.5 settles the number: 0.01, *inherited from production and declared inherited*, with no
source in the parameters fit. What it does not settle is whether that inheritance enters the
locus's warrant, and **it must not**.

A warrant is per `(read group, candidate)`; the outlier weight is one run-wide constant that is
defaulted everywhere. Folding it in would make **every** repeat tract's call `Defaulted` in every
run, and the distinction §4.4 says the warrant exists to carry — *"a genotype resting on a direction
split borrowed from two repeat counts away is distinguishable in the run's output from one resting
on a fit"* — would be gone at every tract. The same line is already drawn one level down:
`PART_REPEAT_SHARE_OF_WHOLE` is a placeholder inside every stutter model this module builds **from
a fit** — a defaulted cell's model comes from `StutterModel::hipstr_shipped()`, whose part-repeat
shares are its own literals — and no provenance mentions it either.

**What the constant is owed instead is a line in the run's output**, which is where §3.6 already
puts the contamination fraction. That is E2b's neighbourhood and is banked there.

## 3. The shape, and the one borrow that decides it

**The contexts borrow the stutter models, so the two cannot live in one struct.** `TractScoringFits`
owns the models, the rates and the per-cell warrants; `scoring_contexts` builds the contexts from a
borrow of it, and they live as long as that borrow.

**It owns the motif too, and that is the review's Major.** The first draft had `scoring_contexts`
take a *second* `&Motif`, never compared against the one `gather_for_locus` looked the models up
under. Nothing stopped a caller gathering a mononucleotide tract and scoring it as a dinucleotide:
the context would then report an unreachable mass computed under one period beside a stutter model
looked up under another — measured at `1.18e-14` against the `2.00e-2` the model actually loses, a
factor of 1.7 × 10¹², with no panic. And since the context's motif also drives the emission, the
wrong motif then reached the row. **The fix removes the disagreement rather than checking for it**:
the motif is stored at gather time and `scoring_contexts` no longer asks. It is the same argument
`SsrScoringContext::new` already makes for taking the unreachable mass from the distribution rather
than from its caller — a fact taken twice is a fact that can disagree. Everything the type owns is a buffer that
`gather_for_locus` clears and refills, so one per worker allocates on the first few tracts and then
stops — **except the contexts**, which cannot outlive one locus and so cost one allocation of
`read groups × candidates` per tract. `#![forbid(unsafe_code)]` closes the usual escape.

**The table covers every read group of the run**, not the ones whose reads reached the tract, for
two reasons: `SsrScoringContextTable::of` indexes by `ReadGroupId` directly, and the row asserts
`contamination.fraction_of_each_read_group.len() == contexts.read_group_count()` (`likelihood/ssr.rs`),
so the two halves of the mixture must be on one axis. So a tract costs `read groups × candidates`
stratum lookups — 6 at one library and six candidates, 6,000 at a thousand libraries. **Nothing in
this repository reaches the second** (the tomato panel is 63 accessions of one library each) and no
measurement here says what it costs against the row's own work at that size. Banked.

## 4. What the tests can fail on, and the fixture shape they were built against

Twenty-nine tests, ten of them added after the reviews (§9). The fixture is **three read groups
against two candidates** — deliberately unequal,
because at an equal shape a table filled read-group-major and one filled candidate-major are the
same length and the same set of cells, so a transposition passes every shape check. The two
candidates are **6 and 11 repeats** of one dinucleotide, five apart, so they land in different
strata; every read group is in its own slippage group; and the slippage level and the substitution
rate differ on **both** axes, keyed by the candidate's repeat count rather than by its position.

**Two mutations were run rather than reasoned about**, both against
`every_cell_carries_its_own_read_groups_and_its_own_candidates_numbers`:

- reading the cells candidate-major while `gather_for_locus` filled them read-group-major leaves all
  six models present and fails at read group 0 / candidate 1, which then carries **read group 1's
  slippage level, 0.067, where its own is 0.064**;
- hoisting the lookup out of the candidate loop — one lookup per read group, at candidate 0's repeat
  count — fails at the same cell with **0.044 against 0.064**.

**The measurements the prose carries**, each from a test that asserts it:

- the reachable-length support for this fixture holds **41 lengths, every whole number of bases from
  2 to 42**, against the two the candidates themselves spell (12 bases and 22);
- the unreachable mass is the candidate's own: **2.02 in 10,000 at 6 repeats against 1.85 in a
  million at 11**, about 109-fold;
- twenty reads of the 11-repeat tract score the three diploid genotypes at **−161.35, −16.37 and
  −2.51 nats** in table order, so the homozygous 11-repeat call wins by 14 nats over the
  heterozygote;
- scoring both candidates under one stratum's numbers moves the heterozygote by **0.221 nats**
  (0.96 Phred) on two levels 0.044 and 0.064 apart — small, and not zero, which is what says the
  candidate axis reaches the row's arithmetic rather than only its parameter table.

**`gathering_a_second_tract_leaves_nothing_of_the_first` was built against a trap the first draft
walked into.** Both tracts reach 41 lengths, so a length buffer appended to rather than cleared
would still be a legal width; what separates them is *which* lengths — 2 to 42 against 20 to 60 — so
the test asserts the first tract's own candidate length, 12, is gone.

**Three of the prose's numbers were wrong in the first draft and the tests caught all three**: the
support was written as 39 and is 41, the unreachable-mass ratio as 1,100 and is 109, the
candidate-axis shift as 2.99 nats and is 0.221. Every figure above was re-read off a run.

**Three of the four figures above are asserted as thresholds rather than as values**, and the
report should not read as though they were pinned: the three genotype scores are asserted only as
*the homozygous 11-repeat call wins*, the mass ratio as `> 50`, the shift as `> 0.1`. Only the
support's 41 lengths and its two ends are asserted exactly. A change that halved either separation
would pass — which is the right bar for a test about *which axis reaches the row*, and the wrong
one for a reader who takes the numbers as pinned.

## 5. The release-held assertion battery

Eight checks outside a test module. All eight were downgraded to `debug_assert!` in one run and
`cargo test --release --lib ng::calling::inference::repeat_tract_parameters --all-features` was run:
**eight tests failed, one per check**, so each is reached by a test that fails without it.

**They do not call disjoint entry points — three of the eight enter `locus_parameters` and two
enter `tract_candidates` — and the first draft gave that as the reason no pair shadows another.**
What actually shows it is the run itself: with all eight downgraded, every one of the eight failed
with *"test did not panic as expected"* rather than with another check's message, so no check is
catching another's fixture. Two of the three messages in `locus_parameters` were reworded for the
same reason: *"reached the row"* appeared in both, so a `should_panic` string could have been
satisfied by the wrong one.

One further check — that the run has at least one read group — **is a `debug_assert!` with its
reason written beside it**: `FrozenParameters` refuses an empty calibration list at construction and
that list is the axis this counts, so no test can reach it. It guards that refusal being relaxed
later. And one is an unconditional `panic!` rather than an assertion — the message an ungathered
`TractScoringFits` gets when asked for a motif — with its own test.

`a_narrower_candidate_set_is_refused_by_the_rows_parameters` exists only because the first battery
run showed that check unreached.

## 6. The sentences this step retires

**Six sentences named step E2 as what a repeat tract was waiting for**, which landed and did not
gather these parameters. All six are in `summarise_condition.rs` — counted at `HEAD`, `grep "step
E2"` returns six lines and `calling/mod.rs` returns none. **Eight places were edited**: those six,
plus two that said the assembly was unwritten without naming a step (`generic_evidence_of`'s doc,
and `LocusInference::weakest_provenance` in `calling/mod.rs`).

Only one of the eight names E3 in so many words, and it is the one a reader will actually hit: the
front-door refusal's own panic message. The rest say what is missing — this driver's route from a
tract's evidence to the row — without a step label.

The plan carried the same class of error in **three** places, all corrected: the intro's gap count,
the E3 entry's "only on the candidates", and the Checkpoint E verification row.

## 7. What is deliberately not here

- **The contaminant seed** — spec §4.5.1's third term, and the one field of `SsrLocusParameters`
  left `None`. It is E2d, whose plan entry this work added. **`locus_parameters` refuses a run whose
  fit found contamination**, naming that step, rather than handing back the two-term form: the
  first draft made that a doc comment saying a caller *must not*, which is the same guarantee with
  no mechanism, and the failure it prevents is silent because the two-term row returns perfectly
  plausible numbers.
- **The driver's wiring**, which is E3.
- **Where a tract's candidates and their repeat counts come from.** `tract_candidates` takes the
  counts because a tract's repeat count is not its byte length divided by the period — an
  interrupted tract holds fewer whole repeats than its bases suggest. Today they come from a
  fixture; when the repeat-tract half of candidate selection lands they come from it, and the
  signature does not change.

## 8. Validation

Run with `./scripts/dev.sh` from the `ng-calling-loop` worktree.

| gate | before | after |
|---|---|---|
| `cargo test --lib` | 4,786 passed / 0 failed / 14 ignored | **4,815 passed / 0 failed / 14 ignored** |
| `cargo test --release --lib ng::calling --all-features` | 725 passed | **754 passed** |
| `cargo test --test ng_calling_loop_allocation --features dhat-heap` | 1 passed | 1 passed |
| `cargo fmt --all -- --check` | exit 0 | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 | exit 0 |
| `cargo doc --no-deps --lib` | 28 unresolved links, exit 101 | **28, exit 101** — none in this step's files |

**`cargo doc` exits 101 and always has**: the crate denies broken intra-doc links, so the 28
pre-existing ones are hard errors rather than warnings. It is in the table because the check that
matters is *this step added none*, and the first draft's row read as though the command passed.

## 9. What three reviews found, and what changed because of it

Three agents in worktrees detached at `a27c3bbd` with the diff applied: one on arithmetic and
control flow, one on tests and mutation, one on design conformance and claim-checking. Every
finding was applied.

**One Major on the code**, and it is §3's: `scoring_contexts` took a motif it never checked.

**Two blind spots in the tests, each hiding a defect that reaches a genotype.** The mutation
reviewer ran 16 mutations: 9 killed, **5 that changed behaviour and survived**, 2 that changed
nothing. The arithmetic reviewer ran 12: 9 killed, 3 survived. They came from one accident and one
ordering:

- **The period is 2 in every fixture, and 2 is also the candidate count.** Replacing
  `motif.ssr_period()` with a hard-coded 2 left all nineteen tests green — and on a trinucleotide
  that mutation defaults **six cells of six**, because the period keys *both* lookups. Closed by
  `a_tract_of_another_period_is_scored_under_its_own_period`, which is also the first time this
  module has been run at any period but 2.
- **The reuse test gathered the fitted tract first and the defaulting one second.** In that order,
  dropping `substitution_rate.clear()`, `warrant.clear()`, or either counter reset leaves all
  nineteen tests green: the stale entries sit past the ones the second tract reads, and a stale
  counter of zero adds nothing. With the order reversed, each deletion now fails — a stale rate
  reads 0.001 where the fit says 0.0016, a stale warrant reads `Defaulted` where the cell is
  fitted, a stale counter reads 3 where the second tract defaults nothing.
- **One test helper set both direction shares' provenance from one argument**, so `warrant_of`'s
  two conjuncts were indistinguishable and dropping either left the suite green. The helper now
  takes both, and a cell that fitted its contraction bias and read its fall-off off a curve is a
  test.

**Seven mutations were re-run after the fixes and all seven now fail a test** — the hard-coded
period, the two dropped `clear()`s, the dropped counter resets, `warrant_of` ignoring the fall-off
share, the rate's own provenance replaced by `FittedHere`, and the narrowed unknown-read-group
matcher. One of the fixes had
the same accident in it: the first version of
`a_library_the_fit_does_not_describe_is_counted_apart` produced **two** cells of each kind of
absence, so narrowing the matcher to the ordinary one still reported the right number. The fixture
now gives four against one.

**Eleven of seventy-one claims were wrong**, all of them explanations rather than figures — every
counted number in the first draft was re-derived and correct, including both mutation levels, the
109-fold mass ratio, the 41-length support and all three genotype scores. The wrong ones are
corrected in place above and each says what it replaced: the substitution rate's fall-through
mechanism (three places), what tells the four slippage absences apart (three places), two plan
claims, `PART_REPEAT_SHARE_OF_WHOLE` being in *every* model, "refused by name" for a refusal that
did not exist, "four numbers a context carries" mixing two containers, the empty fold being
unreachable, and the battery's "disjoint entry points".

**Two smaller things the reviews closed.** An ungathered `TractScoringFits` — reachable, since it
derives `Default` and is public — answered `Provenance::FittedHere` from folding zero cells, which
is the strongest warrant on the ladder from having read nothing; it now refuses, as does asking it
for contexts, which divided the cell count by a stride of zero. And `ReadGroupId(group as u32)` is
now a checked conversion.

## 10. Banked for the owner

1. **Should a candidate off the fit borrow a neighbouring stratum instead of taking HipSTR's shipped
   row?** The shipped row is symmetric where every real fit is contraction-biased. A within-locus or
   nearest-count borrow is closer to the data and is a rule nobody has written. Recommendation:
   leave it as it is until a tract is called on real data and the share of defaulted cells is
   counted — the answer turns on how often this binds, which nothing measures yet.
2. **The run's output should say the outlier weight was inherited**, beside what it already owes for
   the contamination fraction (§3.6). E2b is the step that opens that channel. The count of cells
   whose read group the fit does not describe belongs in the same line: it is the one absence here
   that means the parameters and the reads came from different runs.
3. **The per-tract cost is on the run's read-group axis** — `read groups × candidates` stratum
   lookups per tract, whoever covered it — and is unmeasured above 63 libraries.
