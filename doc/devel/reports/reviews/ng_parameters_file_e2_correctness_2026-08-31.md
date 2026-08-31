# E2 correctness review — a run with no fit assembles from the defaults

Reviewed: the uncommitted step in `tmp/e2_step.patch`, applied to a detached worktree at
`/Users/jose/devel/pop_var_caller-e2-rev1` on base `accd45a1`.
Every number below was produced by running the code in that worktree; nothing was fixed and
nothing was committed.

**Verdict.** The step's arithmetic is sound and the two loosened assertions do not let a wrong
number through — that is proved rather than argued, below. There is **one Major finding**: the
parameters file a defaults run writes now carries, beside every sample the operator did not
name, a comment saying *"a run should not be able to write this line."* That sentence is the
step's own former justification, left behind in `to_toml.rs`, and it is printed into the artefact
the step exists to make writable. The rest are Minors: one guard branch the mutation suite cannot
pin, one vacuous assertion, three doc statements the code contradicts, and one gap that step F1
will walk into.

> **⚑ The patch I was handed is behind the author's own working tree, and four of these are
> already fixed there.** `tmp/e2_step.patch` touches three files;
> `/Users/jose/devel/pop_var_caller` has uncommitted edits to seven, including `to_toml.rs`. The
> Major (M-1) is fixed there, and so are three of the Minors (m-3, m-6, and the `# Panics` half of
> m-5) — each is marked below. **Still open in the author's live tree: m-1, m-2, m-4, m-7, the two
> remaining halves of m-5, and the F1 gap in §4.** Everything measured below was measured against
> the patch as given, in the isolated worktree.

---

## 1. Blocker

None.

---

## 2. Major

### M-1. The file a defaults run writes tells its reader the file cannot exist

> **Already fixed in the author's live working tree** (`to_toml.rs`, uncommitted): the constant is
> rewritten to *"nobody said how inbred this sample is, so it is scored as an outcrosser … if your
> cohort selfs, say so"*, with two `to_toml` assertions pinning it. The finding stands against the
> patch under review; no action is needed if that edit lands with the step.

**`src/ng/calling/parameters_file/to_toml.rs:537-542`** (the constant), reached from
**`to_toml.rs:546-548`** (`where_it_came_from` attaches it to every `Defaulted` warrant) and
**`to_toml.rs:242-250`** (the inbreeding rows).

```rust
/// An inbreeding coefficient nothing could be fitted for — **which the pre-pass has no
/// default for**, so a run should never write one.
pub const INBREEDING_COEFFICIENT: &str = concat!(
    "no coefficient was fitted for this sample, and inbreeding has no default: a run ",
    "should not be able to write this line"
);
```

Before this step that sentence was true, and `defaults.rs`'s own header quoted it as the proof
that a defaults run had no coefficient to write. The step deletes that header paragraph — the
owner's ruling of 2026-08-31 makes zero the run's declared default — but leaves the constant, and
`where_it_came_from` emits it beside **every** `defaulted` inbreeding row. A defaults run's rows
are defaulted by construction, so the sentence now appears in the ordinary case rather than in an
impossible one.

**Concrete input** — this is the step's own fixture, `defaults.rs:869-887`
(`a_defaults_run_writes_a_file_that_reads_back_as_the_same_run`), with `to_toml()` printed:

```
[inbreeding]
# How inbred each sample is, as a fraction in [0, 1) — one row a sample, counted
# over the reference positions that sample's reads covered.
by_sample = [
    { sample = "TS-1", inbreeding_coefficient = { value = 0.9, warrant = "supplied", observations = { covered_positions = 0 } } },
    # no coefficient was fitted for this sample, and inbreeding has no default:
    # a run should not be able to write this line
    { sample = "Ailsa Craig", inbreeding_coefficient = { value = 0.0, warrant = "defaulted" } },
]
```

Nothing machine-readable breaks: `validate()` accepts the file and it round-trips. What breaks is
spec §1.2 goal 3 — *"a person can read it and change one line"* — for the one file format a
defaults run produces. A user who runs without a fit, opens the file the run wrote, and reads down
the cohort-sized axis is told, once per unnamed sample, that the row in front of them is
impossible.

The three other origin notes in that module (`CALIBRATION_MULTIPLIER`, `FLAT_CONCENTRATION`,
`OUTLIER_WEIGHT`, `SUBSTITUTION_RATE`) all still read correctly for a defaults run; this is the
only one the ruling invalidated.

---

## 3. Minor

### m-1. The rate-count guard's *short* branch is no longer pinned by any test

`src/ng/calling/parameters_file/from_run_parameters.rs:186-189`.

The step splits one `assert_eq!` into a disjunction:

```rust
assert!(
    base_quality_rate_by_read_group.is_empty()
        || base_quality_rate_by_read_group.len() == read_group_count,
```

Mutating `==` to `<=` — which admits a rate set *shorter* than the run's read-group axis —
**survives the whole module: 184 passed, 0 failed** (mutation M15 below). The old
`assert_eq!(len, count)` was killed by
`a_rate_set_over_more_read_groups_than_the_run_has_is_refused` (5 rates, 3 read groups), which is
still the only test that reaches this guard, and it only exercises the *long* direction.

**It is not an equivalent mutant, and it is also not a correctness hole.** A short set whose
missing entries all belong to `Defaulted` read groups is newly admitted by the mutant, and it
writes exactly the file the honest set writes, because `warranted_value`
(`from_run_parameters.rs:588-592`) drops the count on a `Defaulted` warrant. A short set missing a
*fitted* read group is caught one frame later by the per-row guard at line 283. So the surviving
mutant proves the top guard is now redundant belt over the per-row brace rather than the thing
protecting a count — and the suite cannot tell the two apart.

The test that would pin it: a 3-read-group run handed a 2-entry rate set, expecting
`base-quality rates for 2 read groups and the run has 3`. I ran exactly that as a probe and the
unmutated code refuses it at `from_run_parameters.rs:186`, so the branch works; it is only
untested.

### m-2. A vacuous assertion in the step's flagship test

`src/ng/calling/parameters_file/defaults.rs:690`:

```rust
assert_eq!(coefficient.get(), 0.0);
assert_eq!(1.0 - coefficient.get(), 1.0);   // line 690
```

The second line cannot fail once the first has passed — `1.0 - 0.0` is `1.0` by IEEE 754, with no
rounding to lose. Its comment explains the genotype prior's `1 − F` weighting, which is worth
saying; the assertion under it asserts nothing about this code.

### m-3. `names_not_in` does not return names "in the order they were given" — *fixed in the live tree*

> The author's uncommitted `defaults.rs:311` now says *"Sorted, not in the order they were
> given"*, with a two-typo test at line 873. The finding below is against the patch.

`src/ng/calling/parameters_file/defaults.rs:289` (doc) against `defaults.rs:298` (code).
The doc promises the order of the statements; `by_sample` is a `BTreeMap<Box<str>, _>` and
`keys()` yields them sorted.

Measured: `nothing_said().and_this_sample("zeta", …).and_this_sample("alpha", …).names_not_in(…)`
returns `["alpha", "zeta"]`. The committed test
(`a_statement_naming_a_plant_the_run_does_not_have_is_reported`) names one bad sample, so it
cannot see the difference. Harmless today; it becomes a wrong error message the moment a caller
prints "the first name I could not find".

### m-4. A mistyped sample name reaches `of_defaults` and is silently ignored

`src/ng/calling/parameters_file/defaults.rs:337-375`. The step supplies `names_not_in` for a
caller to check with, and documents that checking is the caller's job — but `of_defaults`, the
door whose whole argument is that *"a caller cannot leave one out"*, does not call it.

Measured: `of_defaults(&three_lanes_two_plants, diploid(), &nothing_said().and_this_sample("Alisa
Craig", 0.9))` — one letter wrong — returns coefficients `[0.0, 0.0]` and no complaint. The
operator believes they set `F = 0.9` on a selfing landrace and every genotype at every locus is
scored under Hardy–Weinberg. Nothing between that typo and the wrong prior exists yet, because the
caller that would call `names_not_in` is step F1's or the CLI's.

### m-5. Three doc statements the code contradicts

- **`defaults.rs:212-215`** — the type doc: *"A sample this does not name takes
  `DEFAULT_INBREEDING_COEFFICIENT` and is marked `Defaulted`."* False whenever `everyone` is set:
  after `one_value_for_every_sample(0.9)`, a sample the map does not name takes 0.9 and is marked
  `Supplied`. That is precisely the state
  `a_per_sample_statement_lands_on_that_sample_and_overrides_the_run_wide_one` exercises. The
  private field's own doc at line 219 states the rule correctly.
- **`defaults.rs:334-335`** — *"On a run with no read groups or no samples, which
  `Self::of_gathered_values` refuses one frame later."* Arguments are evaluated left to right, so
  the third one, `SequencingBatches::all_together(read_groups)` (line 358), panics first —
  in `sequencing_batches.rs:197-206`, not in `of_gathered_values`. Same reason, different frame
  and a different message than the doc sends the reader to. *(This half is already gone from the
  author's live tree.)*
- **`defaults.rs:834-835`** — one sample and one library is *"the only shape where the read-group
  and sample axes are the same length."* Any run with one library a sample has equal axes; two
  plants sequenced on one lane each is the common case. **Still present** in the live tree, at
  line 882.

### m-6. The test comment overstates what one line is holding up — *fixed in the live tree*

> The sentence is gone from the author's uncommitted `defaults.rs`. Recorded because the
> measurement is what says the module is better protected than the comment claimed.

`defaults.rs:681-683`: *"a build that shipped the constant at 0.5 passed every test in this module
until this line was written."* As committed, setting `DEFAULT_INBREEDING_COEFFICIENT = 0.5` fails
**three** tests (mutation M1) — `a_run_with_no_fit_takes_every_default`,
`a_stated_coefficient_is_supplied_and_an_unstated_one_is_defaulted` and
`a_statement_naming_a_plant_the_run_does_not_have_is_reported` — because the latter two also
compare against the literal `0.0`. The claim may describe a moment during drafting; it does not
describe the module.

### m-7. `is_none_or` is a weakening that cannot be reached

`from_run_parameters.rs:295-296`. Replacing `rate.is_none_or(…)` with `rate.is_some_and(…)`
survives (184 passed, mutation M18), and is **provably equivalent**: `||` short-circuits, so the
second disjunct runs only when `calibration.provenance != Defaulted`, and the assertion at line
283 guarantees `provenance != Defaulted ⟹ rate.is_some()`. Where `rate` is always `Some`, the two
combinators agree. No defect — but the spelling reads as though a missing rate may disagree about
its warrant, which the guard above already forbids.

---

## 4. The two loosened assertions: what they now admit

This is the question the step is riskiest on, so it was settled by running code rather than by
reading. Fixture: `from_run_parameters.rs`'s own `a_fitted_run`, which is the mixed run —
read group 0 `FittedHere`, read group 1 `Defaulted` (its fitted rate came back zero and
`ReadGroupCalibration::from_fitted_rate` refused it), read group 2 `Borrowed`.

| input | before E2 | after E2 | what the file says |
|---|---|---|---|
| rate set covering 0, 1, 2 (the honest one) | accepted | accepted | rg0 812,344 reads; rg1 no count; rg2 640,918 reads |
| rate set covering 0, 2 and a **stranger at key 9** — same cardinality, read group 1's entry gone | panic at the missing lookup | **accepted** | byte-identical `to_toml()` to the honest one |
| rate set covering 0, 2 only (2 entries, 3 read groups) | panic | panic at line 186 | — |
| **empty** rate set, run has a fitted calibration | panic | panic at line 283, naming `FittedHere` | — |
| empty rate set, every calibration `Defaulted` (the defaults run) | panic | **accepted** — the point of the step | no counts at all |

**No count is lost and none arrives from the wrong read group.** The lookup is
`rate_by_read_group.get(&id)` on a `BTreeMap<ReadGroupId, _>`, so a row can only ever read the
rate filed under its own id; the stranger at key 9 is simply never read. And the only row that can
have no rate is one whose warrant is `Defaulted`, which `warranted_value` writes without an
`observations` table at all. The row is therefore the same row whether a rate was offered or not —
which is what the code comment claims, and it is true.

Two consequences of the loosening are worth recording rather than fixing:

- The *rationale* in both comments — *"a short list … writes some other read group's count beside
  a multiplier"* (line 182) — was never achievable through this argument, before or after the
  step: it is a keyed map, not a positional list. The check is a cardinality proxy for *"these
  rates came from a different fit"*, and the step makes it a weaker proxy (see m-1).
- The mutation that reverts the loosening (M14, strict equality) **is killed**, by the step's own
  round-trip test. The loosening is load-bearing and pinned.

### The gap step F1 will hit

`calibration_rows` admits a missing rate **only** under `Defaulted`. A run scoring from a supplied
parameters file has `Supplied` calibrations (spec §2.1 demotes a file fitted under another census)
or the file's original `FittedHere`/`Borrowed` ones — and has no rate map at all, because the file
carries counts and warrants but never the fitted error rates. Measured:

```
thread '…probe_a_supplied_calibration_with_no_rate' panicked at from_run_parameters.rs:282:
read group 0's calibration is Supplied and no rate was offered for it; only a `Defaulted`
calibration can have none …
```

Impl plan F1 is *"one writer, three sources — supplied file, defaults, or fit — after assembly the
run cannot tell them apart, and writes the file beside its VCF unconditionally."* Two of the three
sources work; the supplied-file one panics after the last locus. F1 will have to either widen this
guard again or reconstruct `Estimate`s from `RunParametersFromFile::reads_behind_each_calibration`
— the very "caller inventing an `Estimate` to satisfy a lookup" the step's comment argues against.
Flagging it here so the choice is made deliberately at F1 rather than discovered by a panic.

---

## 5. `of_defaults`, field by field, against spec §8 and §5

All nine hold. Traced through `of_gathered_values` and the consumers, not just asserted.

| field | what `of_defaults` puts there | checked against |
|---|---|---|
| `calibration_by_read_group` | `ReadGroupCalibration::defaulted()` × read groups — scale 1.0, `Defaulted` | §8 "has one": scale of one, marked `Defaulted`. `charged_error` at scale 1 is the geometric mean of the reads' own errors, pinned at `defaults.rs:427` |
| `contamination_by_read_group` | **empty** | §5 row 1 and §8 "absence is the default". `view()` takes `FrozenParameters::uncontaminated`, the plain read-likelihood formula |
| `sequencing_batches` | `all_together` — one batch, `defaulted: true` | §3.4's default; both axes are derived from the same `ReadGroups`, so `of_gathered_values`' two count checks cannot disagree |
| `inbreeding_coefficient_by_sample` | 0.0, or what the run declared | owner's ruling 2026-08-31. `1 − F = 1`, so `calling_priors.md` §7's heterozygote weighting is inert — the spec reference is accurate (§7, line 800) |
| `prior_seed` | `seed_from_moments(None, None)` → `SeedRegime::FallbackDiversity`, `alpha_alt_total = ExpectedHeterozygosity::SPECIES_FALLBACK = 1e-3` | `seed_generic.rs:254-259`; the constant is pinned to `1e-3` at `types.rs:1791` |
| `ssr_slippage_fits` | `StratumFits::over(&[], every read group → group 0)`; `stated_concentration = STATED_FLAT_CONCENTRATION = 1.0`, warrant `Defaulted`; no strata, no length spectra | §8's flat concentration; `length_spectrum_at` falls to `stated_flat(1.0)` (`stratum_fits.rs:695`) |
| `ssr_substitution_rate` | empty | §8: the default is taken at the tract. `TractScoringFits::gather` (`repeat_tract_parameters.rs:350-357`) counts the cell into `substitution_defaulted` and takes `DEFAULT_SSR_SUBSTITUTION_RATE = 0.001` |
| `ploidy` | the run's own | §3.2 |
| `repeat_tract_outlier_weight` | `defaulted()` — `DEFAULT_OUTLIER_WEIGHT = 0.01`, `Defaulted` | §3.8; the literal is pinned at `defaults.rs:458` |

### The slippage-group reasoning holds

`StratumFits::at` (`stratum_fits.rs:771-778`) looks the read group up **first** and the stratum
second:

```rust
let group = *self.slippage_group_of.get(&read_group).ok_or(NoSlippage::UnknownReadGroup)?;
let row   = self.by_stratum.get(&stratum).ok_or(NoSlippage::NoSuchStratum)?;
```

so declaring every read group into group 0 means the empty `by_stratum` answers `NoSuchStratum`,
and `TractScoringFits::gather` (`repeat_tract_parameters.rs:341-347`) adds the cell to
`slippage_defaulted` but **not** to `slippage_defaulted_by_an_unknown_read_group`, which is what
`cells_whose_read_group_the_fit_does_not_describe` reports. Declaring nothing would have put every
cell of every tract into the second counter. The step's claim is correct, and the three mutations
that break it (declare nothing / declare group 1 / declare only as many groups as there are
samples) are all killed.

Two notes on the choice, neither a defect:

- Group 0 for everybody is the same default the joint walk takes with `SLIPPAGE_PER_READ_GROUP`
  unset (`examples/ng_joint_records_walk.rs:775-784`), so the doc's claim about it is accurate.
- Which group is chosen is **inert today** — `by_stratum` is empty, so `0..n` distinct groups would
  answer `NoSuchStratum` too. It starts to matter at E3, when a shipped default row exists; one
  pooled group is then the shape that lets one shipped row cover the run, which looks like the
  right choice to have made.

### The one field the round-trip test does not follow — checked, and covered elsewhere

`a_defaults_run_writes_a_file_that_reads_back_as_the_same_run` asserts eight of the nine fields
through the TOML; `sequencing_batches` is checked only in memory, by the first test. That is not a
gap: mutation M23 forces `sequencing_batches_of`'s `batching_was_declared: !batches.is_default()`
to `true` — a defaults run's file would then claim somebody declared its batching — and it is
**killed** by two pre-existing tests, `a_run_that_declared_no_batching_writes_the_flag_false` and
`of_run_writes_a_single_sample_single_read_group_run`. The `defaulted` flag survives the trip and
is pinned.

---

## 6. Mutation results

25 mutations. Every run was verified to print
`CWD /Users/jose/devel/pop_var_caller-e2-rev1` and, where the text changed,
`Compiling pop_var_caller v0.1.0 (/Users/jose/devel/pop_var_caller-e2-rev1)`; the two touched files
were restored from a pristine copy after each and `diff`ed clean before the next. Suite:
`cargo test --lib ng::calling::parameters_file`, baseline **184 passed**.

Short names: `every_default` = `a_run_with_no_fit_takes_every_default`; `no_stratum` =
`a_defaults_runs_tracts_find_no_stratum_rather_than_an_unknown_read_group`; `supplied/defaulted` =
`a_stated_coefficient_is_supplied_and_an_unstated_one_is_defaulted`; `overrides` =
`a_per_sample_statement_lands_on_that_sample_and_overrides_the_run_wide_one`; `typo` =
`a_statement_naming_a_plant_the_run_does_not_have_is_reported`; `round_trip` =
`a_defaults_run_writes_a_file_that_reads_back_as_the_same_run`; `missing_rate` =
`a_rate_set_missing_one_of_the_runs_read_groups_is_refused`; `wider_fit` =
`a_rate_set_over_more_read_groups_than_the_run_has_is_refused`.

| # | mutation | outcome |
|---|---|---|
| M1 | `DEFAULT_INBREEDING_COEFFICIENT` 0.0 → 0.5 | killed — `every_default`, `supplied/defaulted`, `typo` |
| M2 | `of_each_sample`: per-sample lookup keyed to a name no sample carries | killed — `overrides`, `typo`, `supplied/defaulted`, `round_trip` |
| M3 | `of_each_sample`: run-wide value takes precedence over per-sample | killed — `overrides` |
| M4 | `of_each_sample`: a stated coefficient marked `Defaulted` | killed — `overrides`, `supplied/defaulted`, `round_trip` |
| M5 | `of_each_sample`: an unstated coefficient marked `Supplied` | killed — `supplied/defaulted`, `typo`, `round_trip` |
| M6 | `of_each_sample`: `Supplied` arm `observations` 0 → 7 | killed — `supplied/defaulted` |
| M7 | `of_each_sample`: `Defaulted` arm `observations` 0 → 7 | killed — `supplied/defaulted` only (the file drops the count, so the round-trip cannot see it) |
| M8 | `names_not_in`: predicate inverted | killed — `typo` |
| M9 | `of_defaults`: no read group declared into a slippage group | killed — `no_stratum` |
| M10 | `of_defaults`: slippage group 0 → 1 | killed — `no_stratum` |
| M11 | `of_defaults`: slippage map covers samples, not read groups | killed — `no_stratum` |
| M12 | `of_defaults`: one calibration a sample, not a read group | killed — `no_stratum`, `round_trip`, `every_default` |
| M13 | `of_run`: rate-count guard always passes | killed — `wider_fit` |
| M14 | `of_run`: strict equality restored (the E2 loosening reverted) | killed — `round_trip`. **The loosening is load-bearing and pinned.** |
| M15 | `of_run`: `is_empty() \|\| len == n` → `len <= n` (a short rate set admitted) | **SURVIVED** — see m-1. Not equivalent (it admits inputs the original refuses) but safe: every newly admitted input writes the correct file, because the per-row guard is what protects a count. |
| M16 | `calibration_rows`: missing-rate guard always passes | killed — `missing_rate` |
| M17 | `calibration_rows`: missing rate legal everywhere **but** `Defaulted` | killed — `missing_rate`, `round_trip` |
| M18 | `calibration_rows`: `is_none_or` → `is_some_and` | **SURVIVED** — provably equivalent, proof in m-7 |
| M19 | `calibration_rows`: every count written as zero | killed — 4 tests, including `a_file_read_into_a_run_and_written_back_is_the_file_that_was_read` |
| M20 | `calibration_rows`: the no-rate fallback count is 7 (the author's known equivalent mutant) | **SURVIVED** — confirmed equivalent, and the code comment's "passes all 184 tests of this module" reproduces exactly |
| M21 | `calibration_rows`: every row takes the **first** rate in the map rather than its own | killed — 6+ tests |
| M22 | `calibration_rows`: warrant-agreement guard always fires | killed — 6+ tests |
| M23 | `sequencing_batches_of`: every batching written as `declared` | killed — `a_run_that_declared_no_batching_writes_the_flag_false`, `of_run_writes_a_single_sample_single_read_group_run` |
| M24 | `of_defaults`: seed taken from a stated diversity rather than the fallback rung (`NeutralShape` at the same `alpha_alt_total`) | killed — `every_default` |
| M25 | `of_defaults`: the batching's two axes transposed (samples for read groups) | killed — `every_default`, `round_trip`, `no_stratum` |

**22 killed, 3 survived — one a real (safe) weakening (M15), two provably equivalent (M18, M20).**

---

## 7. Numbers in the patch, re-derived

| claim | where | verdict |
|---|---|---|
| the prior multiplies its heterozygote branch by `1 − F`, `calling_priors.md` §7 | `defaults.rs:191-193` | correct — §7, line 800 of the spec |
| `DEFAULT_SSR_SUBSTITUTION_RATE` is 0.001 | module table, `defaults.rs:27` | correct — pinned at `repeat_tract_parameters.rs:2474` |
| `ExpectedHeterozygosity::SPECIES_FALLBACK` is the seed's bottom rung | `defaults.rs:697-700` | correct — `seed_generic.rs:254`; constant is `1e-3`, pinned at `types.rs:1791` |
| `STATED_FLAT_CONCENTRATION` where a run fitted no stratum | `defaults.rs:703` | correct — `median_concentration` returns it on an empty map (`stratum_fits.rs:988`) |
| three read groups over two samples; `TS-1` is sample 0 | `defaults.rs:783-784` | correct — `group_by_sample` is first-seen order (`read_groups.rs:406-409`) |
| "putting 7 in its place passes all 184 tests of this module" | `from_run_parameters.rs:311` | correct — reproduced exactly (M20) |
| "a build that shipped the constant at 0.5 passed every test in this module until this line was written" | `defaults.rs:681-683` | **overstated** — three tests catch it (m-6) |
| one sample and one library is "the only shape where the read-group and sample axes are the same length" | `defaults.rs:834-835` | **wrong** — any run with one library a sample (m-5) |
| `names_not_in` returns names "in the order they were given" | `defaults.rs:289` | **wrong** — sorted; measured `["alpha", "zeta"]` for `zeta` then `alpha` (m-3) |
| `of_gathered_values` refuses the empty run "one frame later" | `defaults.rs:334-335` | **wrong frame** — `SequencingBatches::all_together` panics first (m-5) |

---

## 8. What was checked and found clean

- No count can be lost or cross-attributed by either loosened assertion (§4, proved by running).
- All nine `RunParameters` fields against spec §8 and §5 (§5 table).
- The `NoSuchStratum`-not-`UnknownReadGroup` reasoning, traced to the two counters in
  `TractScoringFits` that keep them apart (§5).
- `of_defaults` at one sample and one library, and at three read groups over two samples — the two
  ends `CLAUDE.md` §"cohort size" asks about. Both assemble; a defaults run at 3,000 samples is
  the same code with a longer inbreeding vector and no per-sample allocation that grows faster.
- The parameters file a defaults run writes passes `validate()` and round-trips through TOML back
  into equal parameters, for eight of the nine fields.
