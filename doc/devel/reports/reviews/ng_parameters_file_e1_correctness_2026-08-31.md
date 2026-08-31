# E1 — correctness review

Step **E1** of `doc/devel/ng/impl_plan/parameters_file.md` ("the defaults, compiled in"),
reviewed against `doc/devel/ng/spec/parameters_file.md` §§3.3, 3.8, 5, 8.

Worktree: `/Users/jose/devel/pop_var_caller-e1-rev1`, detached at `8877316f`, with
`tmp/e1_step.patch` applied (3 files modified, 60 insertions, 5 deletions; plus the new
282-line `src/ng/calling/parameters_file/defaults.rs`).

## How the runs were made

Every build and test ran through the **worktree's own** `scripts/dev.sh`, which computes
`PROJECT_DIR` from its own location and so bind-mounts the worktree:

    /Users/jose/devel/pop_var_caller-e1-rev1/scripts/dev.sh bash -c \
      'export PATH=/usr/local/cargo/bin:$PATH
       cd /Users/jose/devel/pop_var_caller-e1-rev1 || exit 9
       echo "CWD $(pwd)"; cargo test --lib -- ng::calling:: ng::parameter_estimation::'

Each run printed `CWD /Users/jose/devel/pop_var_caller-e1-rev1` and, after every mutation,
`Compiling pop_var_caller v0.1.0 (/Users/jose/devel/pop_var_caller-e1-rev1)` — both checked
automatically by the driver (`tmp/mutate.sh`), which refuses a run whose `sed` changed
nothing and prints the compile count. The scoped suite is `ng::calling::` plus
`ng::parameter_estimation::`: **1,852 tests, 0 failures** at baseline before and after the
whole mutation table. (`cargo test --lib` over the *whole* crate has 19 pre-existing failures
in `ssr::catalog` and `ng::tandem_repeat` that are unrelated to this step, so the scope was
narrowed to a set that is green at baseline.)

After each mutation the file was restored from its own byte-for-byte backup and
`git status --short` re-checked: `M likelihood/mod.rs`, `M parameters_file/mod.rs`,
`M parameters_file/validate.rs`, `?? parameters_file/defaults.rs` — the applied patch and
nothing else.

---

## Findings

### BLOCKER — the run's own writer produces a file the new rung refuses

**Where.** The new rung is `src/ng/calling/parameters_file/validate.rs:485-499`. What it
refuses is produced by `src/ng/calling/likelihood/mod.rs:358-370`
(`ReadGroupCalibration::from_fitted_rate`) via `src/ng/calling/run_parameters.rs:206-221`
(`RunParameters::assemble`) and written out by
`src/ng/calling/parameters_file/from_run_parameters.rs:252-294` (`calibration_rows`).

**What happens.** `from_fitted_rate` sets `scale = rate / mean_minted_error` and copies the
**rate's own** provenance. `src/ng/parameter_estimation/generic/fallback.rs:144-153`
(`resolve_error_rates`) mints an `Estimate<ErrorRate>` whose provenance is
`Provenance::Defaulted` and whose value is `DEFAULT_ERROR_RATE = 0.001`
(`src/ng/parameter_estimation/generic/mod.rs:309`) for any read group that (a) has fewer than
`MIN_SITES_TO_FIT = 10_000` clean sites of its own, (b) has no sibling read group above that
floor to borrow from, and (c) was handed no `fallback_error_rates`. So the calibration comes
out `Defaulted` at a scale that is `0.001 / (geometric mean of that library's minted error)`
— a number that is 1.0 only by coincidence. `calibration_rows` writes it with
`warrant = "defaulted"`, and `validate` then refuses the file the run just wrote.

**Proved by running code**, not by reading. A temporary probe test added to
`from_run_parameters::tests` (removed afterwards; the tree is back to the applied patch)
handed `RunParameters::assemble` three `Defaulted` rates of 0.001 over the module's own
`the_runs_minted_totals()` fixture, projected the file with `ParametersFile::of_run`, and
called `validate()`:

    PROBE row rg=0 value=0.12499999995166093 warrant=Defaulted observations=None
    PROBE row rg=1 value=0.12499999995166093 warrant=Defaulted observations=None
    PROBE row rg=2 value=0.24999999990264823 warrant=Defaulted observations=None
    PROBE: the file a run wrote is REFUSED: the parameters file cannot be used:
    base_quality_calibration.by_read_group[read_group = 0].error_probability_multiplier
    is 0.12499999995166093, and its warrant is `defaulted`, ... which is a multiplier of 1.0;
    ... change the warrant beside it to `supplied`

The **same probe with the new rung disabled** (`if false && …`, mutation M7) prints
`PROBE: the file a run wrote is ACCEPTED`. So this is a regression E1 introduces, not a
pre-existing break the step merely exposes to a test.

**Reach.** `MIN_SITES_TO_FIT` is 10,000 *sites*, so the ordinary way in is a run over a small
target — one gene, a panel, a small contig — or a shallow single-sample run where no read
group clears the floor. That is the corner `CLAUDE.md` names as the hardest and commits the
caller to handling. Note also that the failure is not a wrong number: it is the *whole run's*
parameters file becoming unreadable, so the run cannot be re-used, re-checked, or fed back in.

**Mitigating, and only that.** Nothing in production calls `RunParameters::assemble` today —
all 33 call sites are inside `#[cfg(test)]` blocks, which
`to_run_parameters.rs:1209` states in as many words. So the defect is latent until the fit is
wired to the caller, which is exactly when nobody will be looking for it.

**Why the step missed it.** The new rung's own comment
(`validate.rs:476-484`) enumerates what can carry a multiplier — *"a multiplier can be fitted
from this read group's own reads, borrowed from its sample's other read groups, or handed over
in a file (`ReadGroupCalibration::from_fitted_rate` copies the rate's own warrant)"* — and
concludes *"the other three warrants are all reachable here and none is checked"*. It counts
`FittedHere`, `Borrowed` and `Supplied` and does not notice that the same copying makes
`Defaulted` reachable too, at a scale that is not one.

**What spec §5 says.** Row three of §5's table reads *"a read group's error rate could not be
fitted | scale 1.0, warrant `Defaulted`"*. So the rung is right and the producer is wrong: a
read group that could not be fitted should take `ReadGroupCalibration::defaulted()` — scale
one — rather than a ratio built from a constant nobody measured against that library.
`from_fitted_rate` should return `None` for a `Defaulted` rate the same way it already does
for a zero one, and `assemble`'s `unwrap_or_else(ReadGroupCalibration::defaulted)` would then
do the right thing with no other change. (Choosing between that and relaxing the rung is the
owner's; the rung matches the spec and the producer does not, so I recommend fixing the
producer.)

---

### MINOR — the module's documentation is the step's deliverable and rustdoc renders none of it

**Where.** `src/ng/calling/parameters_file/defaults.rs:1-70`, reached by
`mod defaults;` at `src/ng/calling/parameters_file/mod.rs:167`.

The file's non-test content is entirely `//!` documentation — spec §8's four defaults, the
three kinds of thing they are, and where each constant lives. But `mod defaults;` is private
and the module has no items, so **`cargo doc` skips it entirely**: the prose renders nowhere
in the crate's documentation, and its intra-doc links are never checked. Under
`cargo doc --document-private-items` they turn out to be **eight broken links** — every one of
them:

    defaults.rs:18  [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`]   no item in scope
    defaults.rs:19  [`DEFAULT_OUTLIER_WEIGHT`]                 no item in scope
    defaults.rs:20  [`STATED_FLAT_CONCENTRATION`]              no item in scope
    defaults.rs:55  [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`]   no item in scope
    defaults.rs:55  [`ReadGroupCalibration`]                   no item in scope
    defaults.rs:57  [`DEFAULT_OUTLIER_WEIGHT`]                 no item in scope
    defaults.rs:59  [`STATED_FLAT_CONCENTRATION`]              no item in scope
    defaults.rs:59  [`StratumFits`]                            no item in scope

The names are imported only inside `#[cfg(test)] mod tests`, so the module scope the `//!`
docs are written against has none of them. The module's stated job — *"this module documents
the set"* and *"the numbers themselves stay where their users are"* — depends on those links
working, since the links are the only route from here to the constants. A `#[allow(unused)]
use` block at module scope, or `pub(crate) use` re-exports, fixes it. (The crate has 25 other
broken intra-doc links at baseline, so `cargo doc` was already not clean; these eight are new
but do not change the count under a plain `cargo doc`, which never looks.)

### MINOR — two of the three float-spelling changes are unpinned, and the third message was left behind

The step changed refusal messages from `{}` to `{:?}` so that "the number they quote is
spelled the way the writer spells it in the file (`1.0`, not `1`)", and
`defaults.rs:242-247` asserts exactly that. Two of those `{:?}` can be reverted with no test
noticing (mutations M11 and M13 below), because both refusals are exercised only at values
whose `Display` and `Debug` spellings are identical — `3.5` for the concentration and
`1.0 + 0.25 = 1.25` for the multiplier. The very edit the rung exists to catch — a person
changing a defaulted `1.0` to `2.0` — is the case where the two spellings differ (`is 2`
against `is 2.0`), and no test covers it. Using a whole number in the fixtures (say
`DEFAULT_ERROR_PROBABILITY_MULTIPLIER + 1.0` and `STATED_FLAT_CONCENTRATION + 2.0`) would pin
both.

Separately, the **outlier-weight rung was left on `{}`** for both its value and its constant
(`validate.rs:831-839`), while its two siblings now use `{:?}`. The step's own test passes for
that rung only because `0.01` spells the same either way. Three parallel rungs, two spellings.

### MINOR — one new test duplicates two existing ones and adds no coverage

`defaults.rs::contamination_defaults_to_a_missing_section_rather_than_to_zeros` asserts four
things, and every one of them is already asserted elsewhere:
`to_run_parameters.rs:946 an_absent_contamination_table_gives_an_uncontaminated_run` covers
the first three, and `to_toml.rs:1755-1759` covers `!to_toml().contains("[contamination]")`.
Both mutations aimed at this claim (M14, M15) killed the pre-existing tests as well. Not a
defect — worth knowing only because this is the one of the six new tests that pins nothing new.

---

## What I checked and found nothing wrong with

- **The three constants and the numbers quoted about them.** `DEFAULT_ERROR_PROBABILITY_MULTIPLIER
  = 1.0`, `DEFAULT_OUTLIER_WEIGHT = 0.01` (`likelihood/ssr.rs:83`), `STATED_FLAT_CONCENTRATION
  = 1.0` (`joint/stratum_fits.rs:352`) — all three match spec §8's first bullet, which names
  "1.0 — no calibration", "0.01, inherited, §3.8" and "1.0".
- **The T1 fixture's arithmetic.** Minted errors `[1e-4, 1e-2, 10^-1.3]` are Phred 40, 20 and
  13, so "spanning Phred 40 to Phred 13" is right; their geometric mean is
  `10^(-7.3/3) = 3.69e-3`, which is 3.7 billion times `MIN_BASE_ERROR = 1e-12`, so the guard
  assertion is doing real work and the equality is not asserting the floor. `log_scale()` is
  `ln 1.0 = 0.0` exactly, so `assert_eq!(…, 0.0)` is exact rather than lucky.
- **The claim about the fifth number.** `defaults.rs:41-47` says a repeat-tract cell with no
  slippage fit takes `StutterModel::hipstr_shipped` and `Provenance::Defaulted` and that the
  count is kept — both true, at `inference/repeat_tract_parameters.rs:347` and `:209`.
- **The other three routes to a `Defaulted` calibration.** `ReadGroupCalibration::defaulted()`
  always gives scale 1.0 (M2/M3 confirm nothing else does). The spec §2.1 demotion in
  `bindings.rs:546` is `weaker_of(warrant, Supplied)`, and `Defaulted` ranks *below* `Supplied`
  (`parameter_estimation/mod.rs:100-107`), so a demoted file's `Defaulted` 1.0 stays
  `Defaulted` at 1.0 and is accepted — covered by
  `validate::tests::the_warrant_a_demoted_file_carries_on_the_bottom_rung_is_accepted`, which
  passes. A `Supplied` rate does give a scale away from one, but no rule constrains `supplied`,
  so that is fine.
- **Rung ordering.** Setting `fallback_length_spectrum_concentration` to `defaulted` at 1.0 on
  the full fixture — which *has* fitted strata — is accepted, and correctly so: only
  `fitted_here` with no fitted stratum is refused (`validate.rs:684-690`).
- **`cargo fmt --check`**: clean. **`cargo clippy --lib --all-targets`**: no new warning; the
  ones present are pre-existing dead code in test and example targets.

---

## Mutation table

15 mutations. Each was applied on the host, compiled and run in the container against the
1,852-test scoped suite, then restored and the restore verified against `git status`.
"Killed by" names the tests among the **six new ones** where they fired, plus the total number
of failures the mutation caused.

| # | file:line | mutation | outcome |
|---|---|---|---|
| M1 | likelihood/mod.rs:229 | `DEFAULT_ERROR_PROBABILITY_MULTIPLIER = 1.0` → `2.0` | **killed**, 75 failures. New tests: `a_read_group_with_no_fitted_rate_…`, `a_defaulted_value_that_is_not_the_binarys_own_number_is_refused`, `each_of_the_three_constants_is_accepted_…`, `contamination_defaults_to_a_missing_section_…`. Also `validate::tests::the_file_a_run_writes_is_accepted` |
| M2 | likelihood/mod.rs:386 | `defaulted()` scale → `DEFAULT_ERROR_PROBABILITY_MULTIPLIER * 1.5` | **killed**, 28 failures. New test: `a_read_group_with_no_fitted_rate_is_charged_what_its_reads_were_minted_with` |
| M3 | likelihood/mod.rs:387 | `defaulted()` provenance `Defaulted` → `FittedHere` | **killed**, 9 failures. New test: `a_read_group_with_no_fitted_rate_…` |
| M4 | likelihood/ssr.rs:83 | `DEFAULT_OUTLIER_WEIGHT = 0.01` → `0.02` | **killed**, 37 failures. New tests: `the_outlier_weight_a_run_inherits_is_the_number_production_uses`, `contamination_defaults_to_a_missing_section_…` |
| M5 | joint/stratum_fits.rs:352 | `STATED_FLAT_CONCENTRATION = 1.0` → `2.0` | **killed**, 4 failures. New test: `a_run_that_fitted_no_stratum_seeds_every_tract_from_the_flat_rung` |
| M6 | joint/stratum_fits.rs:508 | no-stratum warrant `Defaulted` → `Supplied` | **killed**, 3 failures. New test: `a_run_that_fitted_no_stratum_…` |
| M7 | validate.rs:485 | new rung short-circuited (`if false && …`) | **killed**, 1 failure — `a_defaulted_value_that_is_not_the_binarys_own_number_is_refused`, and only it |
| M8 | validate.rs:486 | new rung's `!=` → `==` (refuses the constant, accepts the edit) | **killed**, 47 failures. New tests: `a_defaulted_value_that_is_not_…`, `each_of_the_three_constants_is_accepted_…`, `contamination_defaults_to_…` |
| M9 | validate.rs:494 | new rung quotes the constant `{…}` instead of `{…:?}` | **killed**, 1 failure — `a_defaulted_value_that_is_not_…` |
| M10 | validate.rs:495 | new rung's message says `` `fitted_here` `` where it said `` `supplied` `` | **killed**, 1 failure — `a_defaulted_value_that_is_not_…` |
| M11 | validate.rs:696 | concentration rung quotes the **value** `{value}` instead of `{value:?}` | **SURVIVED** — 1,852 passed |
| M12 | validate.rs:697 | concentration rung quotes the **constant** `{…}` instead of `{…:?}` | **killed**, 1 failure — `a_defaulted_value_that_is_not_…` |
| M13 | validate.rs:491 | new rung quotes the offending **value** `{}` instead of `{:?}` | **SURVIVED** — 1,852 passed |
| M14 | to_run_parameters.rs:199 | absent `[contamination]` → dense table of unmeasured zeros | **killed**, 5 failures. New test: `contamination_defaults_to_a_missing_section_rather_than_to_zeros` (also 4 pre-existing) |
| M15 | to_toml.rs:176 | `[contamination]` header written unconditionally | **killed**, 28 failures. New test: `contamination_defaults_to_…` (also 27 pre-existing) |

Two survivors, M11 and M13, and they are the same defect from two sides — the Minor above.
Every other mutation died, and M7 shows the new rung has exactly one guard: no pre-existing
test covers it, which is what a new rung should look like.
