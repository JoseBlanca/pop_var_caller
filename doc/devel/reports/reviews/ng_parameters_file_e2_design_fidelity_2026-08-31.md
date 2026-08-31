# E2 — design fidelity and truth of the prose

*Review of the uncommitted step `tmp/e2_step.patch` applied to `accd45a1` in the isolated
worktree `/Users/jose/devel/pop_var_caller-e2-rev2`. Design authority:
`doc/devel/ng/spec/parameters_file.md` §2, §3.3–§3.8, §5, §7, §8, §13, plus the two owner's
rulings of 2026-08-31 recorded in `PROJECT_STATUS.md`. There is no architecture document by
design.*

**Verdict: request changes.** 1 Blocker, 8 Majors, 7 Minors. The assembly itself is right —
all nine fields of `RunParameters` take what §8 and §5 require, and I could not fault the
choice of defaults. Every finding but two is a sentence that says something the code or the
spec refutes, which is the same failure mode E1 had.

## What was verified green

- `cargo doc --no-deps` in the worktree: **25 `error: unresolved link` and 23
  `warning: redundant explicit link target` before the patch, 25 and 23 after.** No new
  diagnostic; every intra-doc link the patch adds resolves.
- `cargo test --lib`: **5,555 passed, 0 failed, 13 ignored.** `cargo fmt --check` clean,
  `cargo clippy --lib --all-features` silent.
- The equivalent-mutant claim at `from_run_parameters.rs:310-314` is **true**: replacing
  `rate.map_or(0, …)` with `rate.map_or(7, …)` gives 184 passed / 0 failed on
  `ng::calling::parameters_file`.

---

## 1. Does `of_defaults` deliver spec §8 for every field?

`RunParameters` (`src/ng/calling/run_parameters.rs:97`) has nine fields. Walking them against
`of_defaults` (`src/ng/calling/parameters_file/defaults.rs:337-378`):

| field | what `of_defaults` puts there | §8 / §5 | verdict |
|---|---|---|---|
| `calibration_by_read_group` | `vec![ReadGroupCalibration::defaulted(); n]` — scale `DEFAULT_ERROR_PROBABILITY_MULTIPLIER` = 1.0, `Defaulted` | §8 "**has one**: the base-quality calibration scale (1.0 — no calibration) … marked `Defaulted`"; §5 row 3 | right per the spec — **but see M8**, it lands on the opposite side of the 2026-08-31 error-rate ruling |
| `contamination_by_read_group` | `Vec::<ContaminationView>::new()` | §8 "**absence is the default**: contamination"; §5 row 1 "the contamination table is absent … must not write zeros" | right |
| `sequencing_batches` | `SequencingBatches::all_together(read_groups)`, one batch, `is_default() == true` | §3.4 "written even where no contamination was fitted, because it is a fact about the run rather than about the fit" | right, with a gap — see m1 |
| `inbreeding_coefficient_by_sample` | `DeclaredInbreeding::of_each_sample(...).map(|e| e.value)` — 0.0 `Defaulted` where unstated, the stated value `Supplied` where stated | owner's ruling 2026-08-31; §3.5 "At least one is required" | right |
| `prior_seed` | `seed_from_moments(None, None)` | verified against `seed_generic.rs:254-260`: gives `SeedRegime::FallbackDiversity` at `ExpectedHeterozygosity::SPECIES_FALLBACK` (with `NEUTRAL_ALPHA_REF` on the reference side, which the inline comment does not mention but does not misstate) | right; §8 does not cover the seed and the module header says why |
| `ssr_slippage_fits` | `StratumFits::over(&[], {every read group → group 0})` | §8's third bullet's GIAB defaults do not exist (§12 q1); §3.7 "which slippage group each read group's reads are drawn under — the run's own declaration" | right, and the group-0 declaration is the right call — see below |
| `ssr_substitution_rate` | `BTreeMap::new()` | header table row 5: "a default taken at the tract, not written in the file" | right |
| `ploidy` | the argument | §3.2 "a property of the run rather than of the fit" | right |
| `repeat_tract_outlier_weight` | `RepeatTractOutlierWeight::defaulted()` — `DEFAULT_OUTLIER_WEIGHT` 0.01, `Defaulted` | §3.8, §8 first bullet | right |

**Nothing here is wrong or unjustified as a value.** Two things about the assembly are worth
recording as correct rather than as findings, because both were non-obvious:

- **Declaring every read group into slippage group 0 is right, and the doc's reason checks
  out.** `StratumFits::at` (`stratum_fits.rs:761-795`) looks up `slippage_group_of` *first*
  and returns `NoSlippage::UnknownReadGroup` when it misses, before it ever consults
  `by_stratum`. So an undeclared read group would answer `UnknownReadGroup` — "the run is not
  what it claims" (`stratum_fits.rs:94-97`) — on every cell of every tract. Declaring makes
  the answer `NoSuchStratum`, which that enum's own doc calls "**Ordinary**". The claim at
  `defaults.rs:345-354` is true.
- **"the same default the joint walk uses" is true.** `examples/ng_joint_records_walk.rs:775-785`
  puts every read group in group 0 unless `SLIPPAGE_PER_READ_GROUP=1` is set.

---

## 2. Is `DeclaredInbreeding` the right shape?

**For.** The ruling has three states, not two: nothing said, one value for the run, a value
for a named sample. A bare `(Option<InbreedingF>, BTreeMap<…>)` pair passed as two arguments
would let a caller hand over a per-sample map with no run-wide value and no name for that
combination; the type names it. The name-keyed join is a real invariant with a real scar —
D2's Blocker one level down was a positional join — and `of_each_sample` is the one place it
can be got wrong. And the split between the values (which `RunParameters` keeps) and the
warrants (which only the file keeps) needs *one* function serving both, or the two lists drift.
Three constructors plus one projection is not more than the ruling asked for.

**Against.** `names_not_in` (`defaults.rs:289-303`) has no caller and cannot have one yet: it
exists to let something refuse a mistyped name off a command line, and E2 has no command line
— the command surface is deferred (spec §11, "A command that writes the defaults without
running a caller … belongs with the rest of `pop_var_caller_exp`'s subcommands"). A public
method shaped for an argument nobody has designed is a guess, and it has already gone wrong:
its doc claims an ordering the code does not have (finding M4), and the one test that touches
it passes only because it has a single missing name.

**Verdict: the type is right and one of its four methods is not yet earned.** Keep
`nothing_said` / `one_value_for_every_sample` / `and_this_sample` / `of_each_sample`; either
drop `names_not_in` until the caller that refuses a typo exists, or fix M4. The concrete cost
of keeping it as written is a false comment on an unexercised path.

---

## 3. Are the two loosened assertions in scope for E2, or scope creep into F1?

**Verdict: in scope, narrowly — and having taken the writer on, the step then owed the rest of
what the writer prints, which is where the Blocker comes from.**

The reasoning. E2's claim is "a run with no fit and no supplied file assembles `RunParameters`
from the defaults" — assembly, which needs no writer. The loosening exists only so the step's
last test, `a_defaults_run_writes_a_file_that_reads_back_as_the_same_run`, can run, and that
test is F1's by the plan's own words ("**F1. One writer, three sources.** Supplied file,
defaults, or fit"). Two things pull the other way and they win:

- **The loosening is strictly a tightening for every warrant but one.** The count assertion at
  `from_run_parameters.rs:180-192` goes from "exactly `n` rates" to "`n` rates or none", and
  `calibration_rows` gains a per-read-group assertion (`from_run_parameters.rs:279-296`) that a
  missing rate is legal only under `Defaulted`. The retained refusal test still fires, and I
  checked why: the read group whose rate the test removes is `ReadGroupId(2)`, whose fixture
  rate is `Provenance::Borrowed` (`from_run_parameters.rs:863-866`), so the doc claim at
  `from_run_parameters.rs:2294` is true.
- **E1's Blocker was this caller refusing a file it had just written**, for the second time on
  this plan. Landing a defaults state the writer refuses would have re-armed that trap for F1
  to trip over. Proving the state is writable at the moment it is created is cheap insurance.

So the loosening stays. What does not follow is that the *file* it produces was checked, and
it was not — see B1.

---

## 4 & 5. Findings

### Blocker

**B1. `to_toml`'s `origins::INBREEDING_COEFFICIENT` is unchanged, so the file a defaults run
writes tells its reader that the line above it cannot exist.**
`src/ng/calling/parameters_file/to_toml.rs:537-542`, still reading:

> `/// An inbreeding coefficient nothing could be fitted for — **which the pre-pass has no`
> `/// default for**, so a run should never write one.`
> `pub const INBREEDING_COEFFICIENT: &str = concat!(`
> `    "no coefficient was fitted for this sample, and inbreeding has no default: a run ",`
> `    "should not be able to write this line"`
> `);`

`where_it_came_from` (`to_toml.rs:546-548`) prints that string above every `defaulted` row. I
ran the step's own fixture through `to_toml` and read the artefact. It contains, verbatim:

```
by_sample = [
    { sample = "TS-1", inbreeding_coefficient = { value = 0.9, warrant = "supplied", observations = { covered_positions = 0 } } },
    # no coefficient was fitted for this sample, and inbreeding has no default:
    # a run should not be able to write this line
    { sample = "Ailsa Craig", inbreeding_coefficient = { value = 0.0, warrant = "defaulted" } },
]
```

**How loudly it contradicts itself: twice, in the same file.** Besides the comment above, the
file's own header note (`to_toml.rs:94`) instructs the reader:

> "On most of them there is no built-in number for it to be, so writing one is a claim about
> this caller that no build makes — **do not mark an inbreeding coefficient or a substitution
> rate `defaulted`**."

The file then marks an inbreeding coefficient `defaulted`, three screens down, and
`validate()` accepts it — the step's own test asserts as much. So a user who follows spec
§1.2 goal 3 ("A person can read it and change one line") opens the first artefact a defaults
run produces and finds it instructing them that it should not exist.

`PROJECT_STATUS.md` names this text as E2's job in as many words: *"**Two texts say the
opposite today and are true of the fit rather than of the run**, so E2 makes each say which it
is about rather than deleting either: `parameter_estimation::generic::fallback`'s header …
and `to_toml`'s `origins::INBREEDING_COEFFICIENT`, which prints 'a run should not be able to
write this line' beside a defaulted coefficient in the file and now can."* Neither was
touched. Three sentences need rewriting: `to_toml.rs:537-538` (the doc), `to_toml.rs:539-542`
(the printed origin), and the clause in `to_toml.rs:94`.

### Major

**M2. `of_defaults`'s `# Panics` names a refuser that is never reached.**
`defaults.rs:332-335`:

> "On a run with no read groups or no samples, which [`Self::of_gathered_values`] refuses one
> frame later and for the same reason."

`of_gathered_values` is not reached. Its third argument is
`SequencingBatches::all_together(read_groups)`, and Rust evaluates arguments before the call.
I caught the panic on `ReadGroups::of_lanes(&[])` and the message is:

> "every read of a run belongs to a read group and a run has at least one, so a run whose
> read-group table is empty is one whose read groups went missing"

which is `checked_axes` at `src/ng/parameter_estimation/joint/sequencing_batches.rs:432-436`,
not `of_gathered_values`' "a set of parameters covering none is one whose read-group axis went
missing" (`run_parameters.rs:340-344`). Secondary: "or no samples" is unreachable through this
door — `ReadGroups` derives its sample list from its read groups, which the same file's
`checked_axes` records in a comment at `sequencing_batches.rs:438-443` ("no input reaches
this, so no test can").

**M3. The mutation comment at `defaults.rs:680-684` is false, and I ran the mutation.**

> "comparing a run's coefficient to `DEFAULT_INBREEDING_COEFFICIENT` moves both sides
> together, and a build that shipped the constant at 0.5 passed every test in this module
> until this line was written."

Shipping the constant at 0.5 and deleting only that one line fails **three** tests of this
module: `a_run_with_no_fit_takes_every_default`, `a_stated_coefficient_is_supplied_and_an_
unstated_one_is_defaulted`, and `a_statement_naming_a_plant_the_run_does_not_have_is_reported`
(181 passed / 3 failed). The first of those is the test the comment sits in — its very next
line, `assert_eq!(coefficient.get(), 0.0)` at `defaults.rs:686`, kills the mutant on its own.
The rationale the comment gives is sound in general; the measurement it offers as evidence is
not, and the guard adds nothing the following line does not already provide.

**M4. `names_not_in` does not return names "in the order they were given".**
`defaults.rs:289`:

> "**Every sample this names that the run does not have**, in the order they were given"

The source is `self.by_sample.keys()` (`defaults.rs:299`) over a `BTreeMap<Box<str>,
InbreedingF>` (`defaults.rs:221`), which iterates in lexicographic order of the name. Given
`and_this_sample("zeta", …).and_this_sample("alpha", …)` the result is `["alpha", "zeta"]`.
The one test that covers it (`defaults.rs:757-770`) has a single missing name, so it cannot
see the difference. Either drop the clause or store insertion order.

**M5. "the rule the projection out of the file already applies to every warranted number" is
applied to exactly one warrant.** `defaults.rs:254-256`:

> "**`observations` is zero on every row, whichever warrant it carries** … That is the rule the
> projection out of the file already applies to every warranted number."

`warranted_value` (`from_run_parameters.rs:578-593`) writes
`(warrant != Warrant::Defaulted).then_some(observations)` and its own comment eleven lines up
says "`Supplied` is **deliberately** *not* treated this way — see `of_run`'s three rules." The
reader side (`to_run_parameters.rs:457-462`) maps an absent count to zero and otherwise keeps
what the file says. So under either reading of "projection out", the rule covers `Defaulted`
and nothing else.

The consequence is visible in the artefact and is its own smaller defect: a coefficient an
operator declared writes `observations = { covered_positions = 0 }` (quoted in B1 above),
while the file's header (`to_toml.rs:88-91`) tells the reader that a `supplied` number
carrying `observations` "came that way from another run's file, and those counts are that
run's". Here they are this run's, and they are zero.

**M6. The module's own section on parameters with no default now contradicts the table thirty
lines above it — and the same patch rewrote that table.**
`defaults.rs:98-101`:

> `//! # The two parameters with no default, and the two reasons differ`
> `//!`
> `//! **The slippage numbers are owed a measurement; the inbreeding coefficient is forbidden a`
> `//! default.** They look alike in the file — both are simply absent — and they are not alike.`

The inbreeding coefficient is no longer forbidden a default and is no longer absent from the
file: `defaults.rs:29` now reads "| the inbreeding coefficient, per sample |
[`DEFAULT_INBREEDING_COEFFICIENT`], zero |". The heading, the lead sentence, and the module's
opening line at `defaults.rs:2` ("together with **the two parameters** that have no default")
all need to become one. The bullet under the heading was rewritten; its heading and lead were
not.

**M7. `DEFAULT_INBREEDING_COEFFICIENT`'s doc reintroduces the sentence shape E1's review
already found untrue.** `defaults.rs:187-189`:

> "It is the same kind of default as a base-quality multiplier of one: **not a guess at how
> inbred a plant is**, but the arithmetic that declines to correct for it."

Thirty lines earlier the same module argues the opposite about that very comparison
(`defaults.rs:38-42`): "**A multiplier of one declines to recalibrate; it does not abstain
from a claim.** It leaves every read's error probability at what the instrument minted, which
asserts the instrument was right." `PROJECT_STATUS.md` records this as an E1 Major: *"among
them 'a multiplier of one asserts nothing about the chemistry', which `read_likelihoods.md`
§3.2 and `validate.rs` eleven lines above the rung both refute."* F = 0 is likewise a claim —
that the cohort outcrosses — and on a selfing crop it is a large and wrong one, which the very
next paragraph of the same doc comment says plainly ("a landrace at `F = 0.9` scored at zero
is told every homozygous stretch of its genome is a surprise"). The two paragraphs contradict
each other; the second is the true one.

**M8 (design, for the owner rather than the author). A defaults run takes its reads at the
quality the instrument claimed, which is what the owner ruled against on 2026-08-31 for a read
group the fit could not measure.** The ruling's reason was general: *"a library's real error
rate is never its reported sequencing quality, because the quality scores describe base calling
while the reads also carry mismapping, chimeras and damage — so a library nothing could be
fitted for should be charged a stated rate rather than taken at its word."* A defaults run is
the case where **nothing** was fitted for **anybody**, and `of_defaults` gives every read group
a multiplier of exactly 1.0. Spec §8 does say 1.0, and `of_defaults` could not compute the
alternative — the `0.001 / mean minted error` form needs a `MintedReadErrors` accumulator that
only a pre-pass produces, and a defaults run never reads a read. So the behaviour may well be
forced. What is missing is that nobody says so: the module header spends a full ⚑ paragraph
(`defaults.rs:45-69`) on exactly this asymmetry and does not mention that the door this step
adds lands on the other side of it. Either a sentence in `of_defaults` saying why a defaults
run cannot be charged the stated rate, or an owner's ruling that it should be (which would
need a mean-minted-error estimate from somewhere).

**M1. `parameter_estimation::generic::fallback`'s header is unchanged.**
`src/ng/parameter_estimation/generic/fallback.rs:1-19` still opens "What each parameter falls
back to when its own data will not carry it — **and which of them are allowed to fall back at
all**" and states "**The inbreeding coefficient has one rung and it is not a default.** … it
is fitted or it fails." Read as a rule about the fit that is still true; read as it is written
— an unqualified statement about the caller — it is now false, and E2 was named as the step
that scopes it: *"E2 makes each say which it is about rather than deleting either."* Quieter
than B1 because it is not printed into an artefact a user reads, but it is the other half of
the same instruction. One clause — "the *fit* has one rung" — would do it.

### Minor

- **m1. `of_defaults` cannot be told the run's sequencing batching.** §3.4 makes the batching a
  fact "about the run rather than about the fit", the same class as ploidy, which *is* an
  argument. `of_defaults` fabricates `all_together`. Harmless today — contamination is absent,
  so `FrozenParameters::uncontaminated` never reads it, and `is_default()` marks it honestly in
  the file — but a direct-mode run that knows its batching has no way to state it through this
  door. `defaults.rs:359`.
- **m2. The round-trip test asserts eight of the nine fields.** `sequencing_batches` is checked
  on the assembled run (`defaults.rs:718-720`) but not after the trip through TOML, unlike the
  other eight. `defaults.rs:820-870`.
- **m3. "all 184 tests of this module" is the right number for the wrong module.**
  `from_run_parameters.rs:312`. 184 is `ng::calling::parameters_file` entire; the comment sits
  in `from_run_parameters`, whose own test module has 33. Say which.
- **m4. `hardy_weinberg` is the comparator, not the default prior.** `defaults.rs:687-688`
  attributes the `1 − F` weighting to `hardy_weinberg`, which
  `genotype_prior/mod.rs:72` describes as "the comparator". The claim is true of the production
  path too (`hardy_weinberg.rs:12`: "the same two-branch inbreeding mixture"), so this is a
  citation nit, not a wrong statement. The §7 citation at `defaults.rs:183` is exact —
  `calling_priors.md:800`, inside §7 (749–878), reads "The prior multiplies its heterozygote
  branch by `(1 − F)`".
- **m5. Two documents still say the coefficient is forbidden a default and neither was touched:**
  spec `parameters_file.md` §8 (whose "Three cases" is now four) and
  `doc/devel/ng/impl_plan/parameters_file.md:315` ("it is **forbidden** a default rather than
  owed a measurement, and it is what step E2 turns on"). E1's precedent was to leave the spec
  and record the ruling in `PROJECT_STATUS.md`, which is already done; the plan line is the one
  that now reads as a live instruction.
- **m6. The table row's third column describes the wrong state.** `defaults.rs:29` gives "what
  that is" as "a value the *run* declares" for the row whose whole subject is what a run takes
  when it declared nothing. The declared case is `Supplied` and never appears in a defaults
  table.
- **m7. "the same rule the parameters file's four bindings settled on" points at the wrong
  section.** `defaults.rs:208-210`. §6's four bindings are reference digest, sample list,
  read-group table and census terms, and the sample-list binding is "in order, **by name**".
  Where names-over-positions was actually settled for this quantity is §3.5: "the file writes
  the name beside the value, because the order is the run's and a file that carried only an
  order would be silently wrong against a re-ordered sample list."

## 6. Documentation diagnostics

Both baselines held. Measured in the worktree, before and after `git apply`:

| diagnostic | before | after |
|---|---|---|
| `error: unresolved link` | 25 | 25 |
| `warning: redundant explicit link target` | 23 | 23 |

Every `[link]` the patch adds resolves, including
`crate::ng::parameter_estimation::joint::stratum_fits::NoSlippage::NoSuchStratum`,
`…::NoSlippage::UnknownReadGroup`, `RunParameters::of_defaults`,
`crate::ng::calling::run_parameters::RunParameters::assemble` and `::of_gathered_values`, and
the same-module `[DEFAULT_INBREEDING_COEFFICIENT]` and `[DeclaredInbreeding]`.

## Claims checked and found true

Recorded so they are not re-checked: the `NoSuchStratum`-not-`UnknownReadGroup` argument
(`stratum_fits.rs:761-778`); "the same default the joint walk uses"
(`ng_joint_records_walk.rs:775-785`); "the read group whose rate went missing is `Borrowed`"
(`from_run_parameters.rs:863-866`); `seed_from_moments(None, None)` giving
`SeedRegime::FallbackDiversity` at `SPECIES_FALLBACK` (`seed_generic.rs:254-260`); "`assemble`
… refuses a run with no read groups" (`run_parameters.rs:198-203`, reached with two empty
maps); "putting 7 in its place passes all 184 tests" (ran it); and "it is a pure function of
those, so the two calls cannot disagree" (`of_each_sample` reads only `self` and
`read_groups`, both by shared reference).
