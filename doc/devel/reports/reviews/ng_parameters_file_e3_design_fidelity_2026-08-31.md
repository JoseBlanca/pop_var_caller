# E3 — design fidelity and truth-of-prose review

Reviewed in `/Users/jose/devel/pop_var_caller-e3-rev2`, a detached checkout of `6e434561`
(E2) with `tmp/e3_step.patch` applied and nothing else (`git apply --check --reverse`
passes, so the tree is exactly the patch).

**Verification run.** `cargo test --lib` on the patched tree: 5,563 passed, 0 failed.
`cargo clippy --lib`: clean. `cargo doc --no-deps` before the patch: **25** `error:
unresolved link`, **23** `warning: redundant explicit link target`. After the patch:
**25** and **23**, and the diagnostic *locations* are a byte-identical set, so nothing was
traded. Every intra-doc link the patch adds resolves: with `--document-private-items` the
100 unresolved links include none in `to_toml.rs` or `run_report.rs`, which covers the two
links on the private `no_stratum_was_fitted`. **No new diagnostic.**

Three mutations were run to size the test-strength findings; both survivors are recorded
below with the command that produced them.

---

## Summary

E3 delivers one of the three things `PROJECT_STATUS.md` assigns it, and delivers it well:
the two empty repeat-tract tables in a defaults run's file now say what they mean, in the
reader's language, beside the tables themselves. The two other items the same
`PROJECT_STATUS.md` bullet names as E3's are untouched.

The largest finding is not in the design but in the evidence: **the predicate the whole
step exists for is not pinned by any test.** Hard-coding `strata_with_slippage: 0` — which
makes `every_tract_falls_back()` return `true` unconditionally — passes all 5,563 library
tests.

Counts: **0 Blockers, 6 Majors, 7 Minors.**

---

## 1. Does E3 discharge its brief?

The brief has three parts. Taking them in order.

### Trace — partly, and the doc comment overstates what was traced

The traced answer is right and is stated: `StratumFits::at` looks the read group up
*before* the stratum
(`src/ng/parameter_estimation/joint/stratum_fits.rs:771-778`), so a defaults run's cell
gets `NoSlippage::NoSuchStratum` and
`src/ng/calling/inference/repeat_tract_parameters.rs:347` hands it
`(StutterModel::hipstr_shipped(), Provenance::Defaulted)`. That is the behaviour, and the
patch names it correctly in three places.

But see **Major 4**: the test that claims to have *run* the trace does not run it.

### Make the gap visible in the run's report — one of three items delivered

`PROJECT_STATUS.md:2111` is the bullet that assigns E3 its work, and it assigns three:

> Those three are E3's — its brief is *make the gap visible in the run's report rather than
> silent.*
>
> **(a) An empty slippage table reads as *no read group put a read in any stratum*** …
> Same for an empty substitution-rate table, silently 0.001. **(b) `[fitted_from]`, "What
> these numbers were fitted from", heads the second section of a file where nothing was
> fitted from anything.** **(c) No line says no fit ran**, where "0 of 9 parameter groups
> were fitted from your data" would serve a partial fit too.

(a) is done. (b) and (c) are not touched by the patch. I generated a defaults run's file
from the patched tree (`ParametersFile::of_run` over `a_runs_read_groups()`) and confirmed
both:

- `[fitted_from]` still carries the note *"What these numbers were fitted from. A run whose
  reference, samples or read groups do not match these is refused…"* followed by a
  `reference_digest`, a sample list, a read-group table and a `[fitted_from.census]` block
  of eleven digest terms — on a run that fitted nothing from anything.
- No line anywhere in the file says no fit ran. The reader still has to read the whole
  file and notice that every warrant says `defaulted`.

See **Major 1**.

### Do not invent a fallback — honoured

Nothing in the patch invents a number. The two notes render values read off
`StutterModel::hipstr_shipped()` and `DEFAULT_SSR_SUBSTITUTION_RATE`; no new constant is
introduced anywhere. This part of the brief is cleanly met, and the
`every_share_the_note_quotes_is_a_whole_percent` test
(`src/ng/calling/parameters_file/to_toml.rs:1900`) is a good guard on the one thing the
rendering could quietly get wrong.

### Stopping for a ruling — the citation does not check out

The brief: *"if the traced behaviour scores a tract rather than refusing it, say so and
stop for a ruling."* It does score. The patch's answer
(`src/ng/calling/parameters_file/defaults.rs:1040-1042`):

> Step E3's brief said to stop for a ruling if the trace scored the tract; the owner's
> ruling of 2026-08-31 is that reasonable numbers stand until the GIAB measurement exists,
> so what is owed is that the run says so — this is where it does.

That ruling is not in the tree. `PROJECT_STATUS.md` records exactly two rulings dated
2026-08-31 (`grep -n "RULED [^*]*"` gives them at lines 2112 and 2113): the inbreeding
coefficient defaults to 0, and *a read group the fit could not measure* keeps the stated
error rate. Neither is about slippage. The nearest sentence is inside the second
(`PROJECT_STATUS.md:2113`), and it is about `DEFAULT_ERROR_RATE`:

> 0.001 stands as a placeholder and **will be replaced by a value fitted from GIAB**,
> alongside §8's slippage measurement.

See **Major 2**.

---

## 2. Does Checkpoint E's bar hold? — the seven defaults, one at a time

Checkpoint E: *"a defaults run produces parameters, and says which of them were guessed."*
The plan's verification table restates it as *"a defaults run assembles and its report
names every guessed number."*

There are two candidate surfaces and only one is real today. `RunParameterReport` has no
consumer outside tests — its own doc says so (`src/ng/calling/run_report.rs:47-48`: *"The
output stage that will print it is step 10's and is not written, so today the only callers
are tests"*), and `grep` confirms every `.report(` call site in `src/`, `tests/` and
`examples/` is inside a `#[cfg(test)]` block. **So "a run says X" means today "the
parameters file the run writes says X", and nothing else.** Below, "✓" means the file says
it.

Walking `parameters_file/defaults.rs`'s module table (lines 20-27):

| # | default | can a run state it? | through what |
|---|---|---|---|
| 1 | base-quality multiplier, per read group | **✓** | `[base_quality_calibration].by_read_group` row, `warrant = "defaulted"`, with `origins::CALIBRATION_MULTIPLIER` above it. Caveat already recorded in E1 and in the file's own header: `defaulted` here does not pin the value to 1.0 |
| 2 | repeat-tract outlier weight | **✓** | `[stated_constants].repeat_tract_outlier_weight = { value = 0.01, warrant = "defaulted" }` + `origins::OUTLIER_WEIGHT` |
| 3 | tract ladder's fallback concentration | **✓** | `fallback_length_spectrum_concentration = { value = 1.0, warrant = "defaulted" }` + `origins::FLAT_CONCENTRATION` |
| 4 | contamination — absence | **✗** | **nothing.** See below |
| 5 | repeat-tract substitution rate | **✓ for a run that fitted none** (new in E3); ✗ for a partially-fitted run | `no_substitution_rate_was_fitted()` note above `substitution_rate_by_stratum = []` |
| 6 | slippage, per (stratum × slippage group) | **✓ for a run that fitted none** (new in E3); ✗ for a partially-fitted run | `no_stratum_was_fitted()` note above `slippage_by_stratum_and_group = []` |
| 7 | inbreeding coefficient, per sample | **✓** | `[inbreeding].by_sample` row, `warrant = "defaulted"`, + `origins::INBREEDING_COEFFICIENT` |

(The prior's seed, which the module header deliberately excludes from the seven, is also
stated: `[ordinary_site_prior].rung = "stated_heterozygosity"`.)

**So four of seven were already stated before E3, E3 adds two, and one — contamination —
still cannot be stated at all.** The gaps, named:

- **Contamination's absence is completely silent in a defaults run's file.** `to_toml`
  writes the `[contamination]` section *and its explanatory note* inside one
  `if let Some(contamination) = &self.contamination`
  (`src/ng/calling/parameters_file/to_toml.rs:181-198`), so when the section is absent the
  three-state explanation that would say what the absence means goes with it. In the file
  I generated, the word "contamination" appears nowhere except incidentally in the
  `[sequencing_batches]` note. The only cover is the file header's general line *"An
  absent key is not a zero. A missing section, a missing row and a missing key each mean
  the thing was not measured"* — which is the exact shape of argument E3's own comment
  rejects one section down (*"An empty table and a missing row are two different claims,
  and the paragraph above only covers the second"*,
  `src/ng/calling/parameters_file/to_toml.rs:314-315`). This is not in E3's brief, but it
  is squarely inside Checkpoint E's bar. **Minor 1.**
- **A partially-fitted run states neither #5 nor #6 at the run level.** Both notes fire
  only on an *empty* table. A run with 40 fitted strata out of 4,000 candidate ones writes
  no note, and `RepeatTractFitsUsed.fitted_substitution_rates` is a count with no
  denominator — nothing says how many cells wanted a rate — so a non-zero count cannot be
  read as coverage. What exists for that case is per locus
  (`TractScoringFits::cells_with_no_fitted_slippage`), which is the right grain for *how
  much of this tract fell back* and, as the patch itself argues, the wrong grain for the
  run-level question. `PROJECT_STATUS.md:2111` item (c) asked for the thing that would
  cover this — *"0 of 9 parameter groups were fitted from your data" would serve a partial
  fit too* — and it was not built. **Major 1** covers this.
- **The two surfaces test different predicates for the same claim.** The file's note fires
  on `slippage_by_stratum_and_group.is_empty()`, which is built from
  `each_stratum_and_group_with_numbers()` and therefore skips `None` cells; the report's
  `every_tract_falls_back()` fires on `strata() == 0`, which is `by_stratum.len()` and
  counts a stratum row whose every group is `None`. A run holding such a row writes the
  note *"no stratum was fitted at all"* into its file while its report says
  `every_tract_falls_back() == false`. I could not construct that state through
  `of_gathered_rows` (it only inserts strata that have at least one `Some`), so it is
  reachable only through `StratumFits::over` with an all-`None` outcome — theoretical, but
  it is two spellings of one rule. **Minor 2.**

**Verdict on the bar: it does not yet hold.** Six of seven can be stated by a defaults
run's file; contamination cannot, and no line of the file says a fit did not run.

---

## 3. Is `RepeatTractFitsUsed` in the right place and the right shape?

**Place: yes.** The argument in `src/ng/calling/run_parameters.rs:609-613` is correct and
is the reason this belongs on the report rather than on a locus — *"this answers did this
run fit any slippage at all, which no locus can"*. `RunParameterReport` is the only
run-level surface in the tree. Public fields match the module's existing convention
(`ReadGroupContamination`, `src/ng/calling/run_report.rs:189-228`, is all public fields
plus one derived predicate).

**Shape: two objections, one of which lands.**

*The module-doc contract.* `src/ng/calling/run_report.rs:15-17` still reads:

> Every number here is read off [`RunParameters`] unchanged. **Nothing is computed and
> nothing is summarised**: a report that averaged two of a sample's fractions would erase
> the one distinction the finer grain exists to express.

None of the three new fields is read off `RunParameters` unchanged. Two are cardinalities
(`by_stratum.len()`, `BTreeMap::len()`) and the third is a filtered walk over the
read-group axis joined against the fit's map. The patch updated the module doc's opening
sentence and left this one, and it also wrote *"Plain data with no arithmetic … four parts
read off four different runs"* at `src/ng/calling/run_report.rs:64-66`, which is now false
of the fourth. **Major 3.**

I do not think the fields should be removed. The sentence's actual argument is against
*aggregating values whose difference is the point* — averaging two read groups' fractions
— and a count of strata aggregates nothing a consumer wanted individually. The right fix
is to narrow the sentence to the claim it makes, and to say that this one part is derived
and from what.

*Zero as a state.* `strata_with_slippage: usize` uses `0` to carry a qualitatively
different claim, in a module whose sibling type expresses exactly that distinction with an
enum arm — `ContaminationUsed::NoneFitted`, documented at
`src/ng/calling/run_report.rs:151-156` as *"Absent, not a fitted zero."* And the zero is
doubly ambiguous: `StratumFits::strata`'s own doc
(`src/ng/parameter_estimation/joint/stratum_fits.rs:863-868`) says *"It cannot tell a
cohort with no repeat tracts from one where every stratum was refused"*, so a run over a
genome with no repeat tracts writes *"no stratum was fitted at all, so every repeat tract
in this run was scored under…"* — vacuously true, and alarming to read. **Minor 3.**
`every_tract_falls_back()` mitigates this for a caller that uses it; the public field does
not.

---

## 4. Is anything the prose asserts untrue?

Six Majors and four Minors. Taking the ones the brief named first.

### Checked and **true**

- **`StutterModel::hipstr_shipped`'s numbers.** `src/ng/alignment/stutter.rs:308-317`:
  `whole_repeat_longer_share: 0.05`, `whole_repeat_shorter_share: 0.05`,
  `part_repeat_longer_share: 0.01`, `part_repeat_shorter_share: 0.01`. The run-report
  doc's *"one read in twenty comes back a whole repeat short, one in twenty a whole repeat
  long, and one in a hundred each way for a part-repeat slip"*
  (`src/ng/calling/run_report.rs:270-273`) is right, and so are the four rendered values in
  the file (I read the generated file: *"5 reads in 100 report a whole repeat short and 5
  in 100 a whole repeat long, with 1 in 100 short and 1 in 100 long for a part repeat"*).
- **HipSTR's fitted values being contraction-biased.** `src/ng/alignment/stutter.rs:304-305`:
  *"Note the shipped row makes expansion and contraction **equal**. That symmetry is a
  starting point, not a claim — HipSTR's *fitted* values are contraction-biased."* The
  patch's paraphrase is faithful.
- **A defaults run's `read_groups_with_no_slippage_group` is empty, and why.**
  `RunParameters::of_defaults` builds
  `slippage_group_of_each_read_group` over `0..read_groups.len()`
  (`src/ng/calling/parameters_file/defaults.rs:384-394`), so
  `StratumFits::slippage_group_of` answers `Some(0)` for every read group and the filter
  keeps none. Correct, and the reason given is the right one.
- **What `NoSlippage`'s variants mean.** `NoSuchStratum` and `GroupPutNoReadHere` are
  marked *"Ordinary"*; `UnknownReadGroup` and `GroupNotInTheFit` are marked *"The run is
  not what it claims"* (`src/ng/parameter_estimation/joint/stratum_fits.rs:84-105`), and
  `repeat_tract_parameters.rs:340-345` does count the latter two apart. The patch's
  claims about them are accurate — with one incompleteness, **Minor 4** below.
- **`DEFAULT_SSR_SUBSTITUTION_RATE` is 0.001** (via
  `generic::DEFAULT_ERROR_RATE`, `src/ng/parameter_estimation/generic/mod.rs:309`), and
  *"Base quality inside tracts is usually worse than outside them, so on real reads that
  number is likely optimistic"* is the constant's own documented claim, sourced to
  `parameter_prepass_ssr.md` §4.1 (`repeat_tract_parameters.rs:126-129`).
- **The rounding-guard test's own arithmetic.** *"a shipped share of, say, 0.035 would be
  shown to a user as '4 in 100'"* — `(0.035 * 100.0).round()` is `4.0` in Rust
  (half away from zero). Correct.

### Majors

**Major 1 — E3 claims to have made the gap visible in the run's report; two of the three
things `PROJECT_STATUS.md` assigns to E3 are not built, and nothing records that.**
`PROJECT_STATUS.md:2111`, quoted in §1 above, names (a), (b) and (c) and says *"Those
three are E3's."* Only (a) shipped. The patch touches no documentation — no
`PROJECT_STATUS.md` entry, no impl report, no change to
`doc/devel/ng/impl_plan/parameters_file.md` — so a reader arriving at Checkpoint E has no
way to learn that (b) and (c) were consciously deferred rather than overlooked. Given that
the checkpoint bar is *"says which of them were guessed"* and (c) is precisely the line
that would say it, this is the finding I would act on first.

**Major 2 — `defaults.rs:1040-1042` cites an owner's ruling that is not recorded anywhere
the reader can find.** Quoted in §1. `PROJECT_STATUS.md` holds two rulings dated
2026-08-31 and neither is about slippage. The brief for E3 said to *"say so and stop for a
ruling"*, so this sentence is the load-bearing one for whether the step was allowed to
proceed at all — and it is the one sentence a reader cannot check. Either the ruling needs
its own `PROJECT_STATUS.md` bullet (in the shape the other two 2026-08-31 rulings take), or
the sentence should point at whatever does record it.

**Major 3 — the module doc still says nothing here is computed or summarised, and the
constructor doc still says the four parts are plain data.**
`src/ng/calling/run_report.rs:16` and `:64-66`. Argued in §3. Both sentences are false of
`RepeatTractFitsUsed`, and the patch edited the adjacent line in each case without
touching them.

**Major 4 — `a_defaults_runs_tract_is_scored_under_the_shipped_model_and_counted` does not
go through the assembly it says it goes through, and counts nothing.**
`src/ng/calling/parameters_file/defaults.rs:1063-1068` claims:

> This goes through the caller's own assembly, so it is the same lookup a locus makes:
> `TractScoringFits::gather_for_locus` asks `StratumFits::at` per `(read group,
> candidate)`, gets `NoSuchStratum`, and takes `StutterModel::hipstr_shipped` with
> `Provenance::Defaulted`.

The test body (lines 1069-1093) calls `run.ssr_slippage_fits().at(...)` directly and then
compares four literals against `StutterModel::hipstr_shipped()`. `TractScoringFits` does
not appear in it, `gather_for_locus` is never called, and nothing is counted despite the
test's name ending `_and_counted`. **This is the trace the brief asked for, asserted rather
than executed.** The underlying behaviour *is* covered — by
`repeat_tract_parameters.rs`'s own tests at lines ~1280-1375, which do drive
`gather_for_locus` and check `cells_with_no_fitted_slippage` — but for a *partially*
fitted run, not for a defaults run, and those tests predate E3.

**Major 5 — the run-report doc says the parameters file cannot distinguish the two states,
in the same commit that made the file distinguish them.**
`src/ng/calling/run_report.rs:268-270`:

> Nothing in the parameters file distinguishes *no read group put a read in that stratum*
> from *no stratum was ever fitted, so every tract falls back*, because both are the same
> empty table. This is what says which.

The other half of this patch is `no_stratum_was_fitted()`, whose entire purpose is to make
the file say which — and whose own comment
(`src/ng/calling/parameters_file/to_toml.rs:314-319`) argues the same distinction. So both
the premise and the conclusion (*"This is what says which"*, implying uniqueness) are
false as of this change.

**Major 6 — the note tells the user there is nothing in the file to edit, and the file
format says otherwise.** `src/ng/calling/parameters_file/to_toml.rs:478-480`:

> There is nothing in this file to edit: fit the run, or read the calls at repeat tracts
> as resting on somebody else's chemistry.

A hand-written `slippage_by_stratum_and_group` row *is* read back:
`to_run_parameters.rs:309-338` walks that table row by row and
`StratumFits::of_gathered_rows` (`stratum_fits.rs:548`) builds the fit from it. The
sentence also contradicts the paragraph printed 20 lines above it in the same file, which
spends five paragraphs explaining exactly how to write such a row
(`share_of_reads_that_slip`, `shorter_share`, `fall_off`,
`share_of_reads_that_slip_origin`, `shorter_share_and_fall_off_origin`,
`curve_fitted_on`). And it cuts against §7's stated purpose, which
`parameters_file/defaults.rs:10-12` gives as *"defaults a person can see and edit"*. The
honest third option — *write your own slippage rows here; the keys are described above* —
is the one the file exists to offer.

### Minors

**Minor 1 — a defaults run's file is entirely silent about contamination.** Argued in §2.
Not E3's brief, but it is the one remaining default that Checkpoint E's bar cannot be met
for.

**Minor 2 — the file's note and `every_tract_falls_back()` test different predicates.**
Argued in §2.

**Minor 3 — `strata_with_slippage == 0` carries a state, in a module whose sibling type
uses an enum arm for exactly that.** Argued in §3.

**Minor 4 — `read_groups_with_no_slippage_group` reports one of the two "the run is not
what it claims" absences, and the doc reads as though it reports both.**
`src/ng/calling/run_report.rs:287-291` says the state is *"the run is not what it claims,
which is why `NoSlippage` gives it a variant of its own and the locus counts it apart from
the ordinary absences."* The locus counts *two* variants apart —
`repeat_tract_parameters.rs:341-345` matches `UnknownReadGroup | GroupNotInTheFit` — while
the report's list is built from `slippage_group_of(group).is_none()` and therefore holds
only the first. `StratumFits::at`'s own comment
(`stratum_fits.rs:780-783`) calls `GroupNotInTheFit` *"the same class of fact as
`NoSlippage::UnknownReadGroup`"*. So a run can report loci with unknown-read-group cells
beside an empty `read_groups_with_no_slippage_group`, and nothing explains the difference.

**Minor 5 — `no_stratum_was_fitted`'s doc says the note "cannot come to disagree with the
model their tracts were actually scored under"** (`to_toml.rs:462-463`). The note reads its
shares from `StutterModel::hipstr_shipped()`; the model a tract is actually scored under is
chosen independently at `repeat_tract_parameters.rs:347`, which calls the same constructor.
The two agree because both name one function, not because they share a binding — so a
change at the fallback site alone would make the note wrong in silence. "cannot" is one
word too strong for name-coupling.

**Minor 6 — *"a whole repeat short"* is ambiguous, and the two readings differ.**
`whole_repeat_shorter_share = 0.05` is the total mass on contractions of *any* whole number
of repeats (`Regime::probability`, `stutter.rs:766-784`, is
`share × one_step × (1−one_step)^(size−1)`, summing to `share`). Reads exactly one repeat
short are `0.05 × 0.95 = 4.75` in 100. The note's *"5 reads in 100 report a whole repeat
short"* reads naturally as the second. *"a whole repeat or more short"* would settle it.

**Minor 7 — two claims in the file's user-facing note carry no number, in a project whose
own writing rule forbids that.** `to_toml.rs:476-478`: *"they are symmetric where real
slippage is usually biased towards the shorter tract. A PCR library slips more than this."*
CLAUDE.md: *"Never assert a property without its size, its subject, and its measure …
Words like real, large, significant … are placeholders for a number."* The claims are
inherited from `PROJECT_STATUS.md:2111` (*"wrong in both magnitude and shape for a PCR
library"*), which does not carry a number either. This tree does hold one usable figure —
`PROJECT_STATUS.md:1321` records HG002's known-homozygous loci at *"2.0% slippage at ≥6
repeats and a 3.4× direction split"* — and quoting it would turn "usually biased" into
something a geneticist can weigh against their own library. The `no_substitution_rate_was_fitted`
note next door is the model here: it names the number (0.001) and gives the direction of
the error.

---

## 5. Test strength — two surviving mutations

Both were run in the worktree and reverted; the tree is back to patch-only.

**Survivor A (Major-weight) — `strata_with_slippage: self.ssr_slippage_fits.strata()` →
`strata_with_slippage: 0`.** `cargo test --lib`: **5,563 passed, 0 failed.**

This is the field the type exists for. With it constant-zero, `every_tract_falls_back()`
returns `true` for every run in the project — including a fully fitted one — and nothing
notices. Reading the patch's tests explains why: three of them assert
`every_tract_falls_back()` is *true*, and none asserts it is ever false. The patch's own
`a_run_holding_substitution_rates_and_no_slippage_reports_both` was written after
discovering exactly this weakness on the neighbouring field (its doc says *"reporting a
constant zero there passed every other test in `ng::calling`"*) — the same check was not
then made for `strata_with_slippage`. **What is missing is one fixture with a fitted
stratum whose report says `strata_with_slippage > 0` and `!every_tract_falls_back()`.**

**Survivor B (Minor-weight) — the substitution-note guard reads the slippage table.**
Replacing `if tracts.substitution_rate_by_stratum.is_empty()`
(`to_toml.rs:363`) with `if tracts.slippage_by_stratum_and_group.is_empty()`:
`cargo test --lib parameters_file`: **192 passed, 0 failed.**

The two fixtures the file-level tests use have both tables empty (`of_defaults`) or both
non-empty (`a_file_using_every_shape`, whose `slippage_by_stratum_and_group` and
`substitution_rate_by_stratum` are both populated at
`parameters_file/mod.rs:1559` and `:1635`), so the two conditions are never separated. A
fixture with slippage fitted and no substitution rates — a state the patch's own
`a_run_holding_substitution_rates_and_no_slippage_reports_both` argues is reachable, in the
other direction — would kill it.

**Not a survivor:** `read_groups_with_no_slippage_group` is well pinned;
`a_read_group_the_slippage_fit_does_not_name_is_named_in_the_report` asserts the exact
vector and the walk direction it depends on.

---

## 6. What is good, and worth keeping as written

- The placement of both notes is right: each sits immediately above the key it is about, so
  a reader scanning for `slippage_by_stratum_and_group = []` meets the explanation on the
  way. I confirmed this on the generated file rather than from the source.
- Deriving the four shares from `hipstr_shipped()` instead of typing them, and then pinning
  the rounding assumption with `every_share_the_note_quotes_is_a_whole_percent`, is the
  right shape for a note that quotes a constant — and the test's failure message tells the
  next person what to do (*"the note rounds it and would show a reader a number the model
  does not hold"*).
- `unwrapped_comments` searching the sentence rather than the line, with the reason given,
  is a good answer to a real hazard in testing wrapped output.
- The `read_groups_with_no_slippage_group` walk direction, and the comment explaining that
  asking the map can never report what the map is missing, is correct and worth having in
  prose.
- The comment at `to_toml.rs:314-319` — *"An empty table and a missing row are two
  different claims, and the paragraph above only covers the second"* — is the clearest
  statement of the defect in the tree, and it names where the observation came from.
