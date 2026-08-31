# E1 — design fidelity and truth of prose

**Step:** E1, "The defaults as named constants with their origin"
(`doc/devel/ng/impl_plan/parameters_file.md`, Milestone E).
**Reviewed in:** `/Users/jose/devel/pop_var_caller-e1-rev2`, detached at `8877316f` with
`tmp/e1_step.patch` applied.
**Spec:** `doc/devel/ng/spec/parameters_file.md` §§2, 2.1, 3.3, 3.5–3.8, 5, 8, 12 q1, 13.

## How the builds were run, and which findings rest on one

The command form in the brief tests the wrong tree. `/Users/jose/devel/pop_var_caller/scripts/dev.sh`
bind-mounts only its own directory, so `cd /Users/jose/devel/pop_var_caller-e1-rev2` fails inside
the container and cargo compiles the main repo. Every measurement below was taken with the
**worktree's own copy** of the script,
`/Users/jose/devel/pop_var_caller-e1-rev2/scripts/dev.sh`, and each run's log names the worktree:
`Documenting pop_var_caller v0.1.0 (/Users/jose/devel/pop_var_caller-e1-rev2)`. The `defaults`
tests exist only in the patched worktree and six of them ran, which is a second proof the right
tree was built.

Findings that rest on a build: **Major 4** (the broken links, from
`cargo doc --no-deps --document-private-items`) and **item 5** (the two baselines). Everything
else is read off the source.

---

## 1. Does the step deliver spec §8?

> **§8:** "the default for every parameter is a named `pub const` in the source with its origin
> recorded beside it, in this repo's existing convention".

**No — not for every parameter, and the module's own prose says otherwise.** Walking the nine
fields of `RunParameters` (`src/ng/calling/run_parameters.rs:97`), what a run with no fit would
put in each, and whether a compiled-in default with a recorded origin exists:

| field | what a no-fit run puts there | compiled-in default with origin? | named by E1? |
|---|---|---|---|
| `calibration_by_read_group` (:98) | `ReadGroupCalibration::defaulted()` — scale 1.0, `Defaulted` | **yes**, `DEFAULT_ERROR_PROBABILITY_MULTIPLIER` (`likelihood/mod.rs:229`), new here | yes |
| `contamination_by_read_group` (:99) | an empty axis; `view()` reads it as uncontaminated | **absence**, correctly no constant | yes |
| `sequencing_batches` (:105) | `SequencingBatches::all_together` (`sequencing_batches.rs:172`, doc: "**Every read group in one batch — the default**") | a named constructor, not a `pub const` — but it is a structure, not a number | no |
| `inbreeding_coefficient_by_sample` (:106) | **nothing — there is no answer** | **no.** `generic/fallback.rs:14` states the rule: "**The inbreeding coefficient has one rung and it is not a default.** It may be *supplied* … and otherwise it is fitted or it fails." No `InbreedingF` constant exists anywhere in `src/`. | **no** |
| `prior_seed` (:107) | `seed_from_population_moments(None, None)` → `SpectrumSeed::new(NEUTRAL_ALPHA_REF, ExpectedHeterozygosity::SPECIES_FALLBACK, SeedRegime::FallbackDiversity)` (`genotype_prior/seed_generic.rs:255`) | **partly.** `SPECIES_FALLBACK` (`types.rs:853`) is a `pub const` with a full origin. `NEUTRAL_ALPHA_REF` (`seed_generic.rs:36`) is a private `const`, so §8's letter is not met. And the default is marked by `SeedRegime::FallbackDiversity`, **not** by `Provenance::Defaulted` — the file has no warrant on this section at all, only a `rung` key | **no** |
| `ssr_slippage_fits` (:108) | an empty `StratumFits`; the stated concentration falls to `STATED_FLAT_CONCENTRATION`, `Defaulted` | **yes** (`stratum_fits.rs:352`); the slippage numbers themselves have none, correctly | yes |
| `ssr_substitution_rate` (:109) | an empty map; each tract cell then takes `DEFAULT_SSR_SUBSTITUTION_RATE` and `Provenance::Defaulted` (`inference/repeat_tract_parameters.rs:130`, `:356`) | **yes**, with an unusually candid origin ("the argument for the number is thinner than it looks") | **no** |
| `ploidy` (:110) | the run's own declaration — §3.2 calls it "a property of the run rather than of the fit", so no default is owed | n/a | no |
| `repeat_tract_outlier_weight` (:119) | `DEFAULT_OUTLIER_WEIGHT`, 0.01, `Defaulted` | **yes** (`likelihood/ssr.rs:82`) | yes |

**The one field with no answer is the per-sample inbreeding coefficient.** Everything else either
has a compiled-in default with its origin recorded, or legitimately needs none.

That gap is not E1's invention — the spec's §8 does not list the coefficient either, and its three
cases ("has one" / "absence is the default" / "has one and it has to be measured") have no slot for
"may not be defaulted at all". But E1 is the step whose job is the inventory, and the inventory it
wrote asserts completeness it does not have. See Major 2.

---

## Blocker

None.

---

## Major

### Major 1 — "A multiplier of one asserts nothing about the chemistry" is false, and the spec says the opposite

`src/ng/calling/parameters_file/defaults.rs:23`, repeated at
`src/ng/calling/likelihood/mod.rs:213-218` ("this one asserts nothing about the chemistry at all").

A multiplier of one makes the model charge every read exactly the error probability the instrument
minted. That is a claim about the chemistry: *this library's reported quality scores are right.*
`read_likelihoods.md` §3.2 exists precisely because that claim is doubtful, and its own comparison
table ends with the row

> | **ours, today** | yes … | nothing; the reported quality is used as reported |

against four competitor callers of which it says "**all four refuse to take it as reported**". So a
defaulted multiplier puts the caller back in the row the spec identifies as the deficiency §3.2's
calibration was written to remove — the opposite of asserting nothing.

The project's own code says it more briefly. `validate.rs:459`, eleven lines above the new rung:

> "Above one is legitimate and **common**: it says the instrument was optimistic."

If a multiplier above one is common, then one is a substantive claim about a library, not the
absence of a claim.

What *is* true is the weaker sentence already in the same paragraph — "A run that took it is not
guessing at a quantity; it is declining to change one" — and the spec's own framing, §8's last
paragraph: "A tomato PCR library taking a human PCR-free slip rate is a guess in a way that a
scale of one is not." The fix is to drop "asserts nothing about the chemistry at all" (and the
table's `asserts nothing` cell) for something like *declines to recalibrate*, which is what the
code does and what §8 licenses.

### Major 2 — the inventory claims to be complete and is not

`defaults.rs:1` — "**The numbers a run takes from the binary when nothing measured them**, in one
place"; `:11` — "# The four, and the three different kinds of thing they are"; `:38` — "# The
fifth number, which has no default and is not one of these".

The structure of the header is a closed set: four, plus a named fifth exception. A reader takes
that as the whole list. It is not:

- **the per-sample inbreeding coefficient is absent, and it is the one parameter that genuinely
  cannot be defaulted.** `generic/fallback.rs:14` — "The inbreeding coefficient has one rung and
  it is not a default … a cohort's diversity divides by `1 − F`, so a wrong constant would be
  amplified rather than absorbed." A defaults run needs one row a sample (§3.5: "At least one is
  required"), and nothing in the tree can supply one. This is the sixth number, it has no default,
  and unlike the slippage numbers it is not waiting on a measurement — it is *forbidden* a default.
  It belongs in the header beside the fifth, and its absence is exactly the state CLAUDE.md's
  range rule calls out: *emit it as absent* is a legitimate answer, a silent fitted zero is not.
- **the repeat-tract substitution rate has a compiled-in default that is not named.**
  `DEFAULT_SSR_SUBSTITUTION_RATE` (`inference/repeat_tract_parameters.rs:130`) is a `pub const`
  with its origin recorded, taken with `Provenance::Defaulted` at
  `repeat_tract_parameters.rs:356`, and counted (`cells_with_no_fitted_substitution_rate`). It is
  a field of `RunParameters` (:109) and a section of the file (§3.7's last bullet). It meets §8
  already — the header simply does not list it.
- **the ordinary-site prior's seed has one too, and it is marked differently.** With no fit the
  seed is `NEUTRAL_ALPHA_REF` and `ExpectedHeterozygosity::SPECIES_FALLBACK`
  (`seed_generic.rs:255`), and what says so is `SeedRegime::FallbackDiversity`, **not**
  `Provenance::Defaulted`. E1's brief is "Each marked `Defaulted` when used"; here the mark is a
  different enum in a different section of the file, with no warrant and so no `validate` rung
  possible. Whether that is acceptable is an owner's question, but the header should not be silent
  about the one default in the file that is marked by something other than a warrant. (`NEUTRAL_ALPHA_REF`
  is also a private `const`, which is §8's `pub const` requirement missed on a technicality.)

None of these need new code in E1. They need the header to say four-plus-two-plus-one rather than
four-plus-one, and to name the inbreeding coefficient as the parameter with no default at all.

### Major 3 — "whose `defaulted` constructor is the only thing that produces one" is false

`defaults.rs:55-56`:

> - [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`] beside [`ReadGroupCalibration`], **whose `defaulted`
>   constructor is the only thing that produces one**;

`to_run_parameters.rs:136`, in the sibling file, builds the struct literally:

```rust
calibration_by_read_group[at] = Some(ReadGroupCalibration {
    scale: row.error_probability_multiplier.value,
    provenance: row.error_probability_multiplier.warrant.into(),
});
```

So every run that reads a parameters file with `warrant = "defaulted"` gets a
`Provenance::Defaulted` calibration at a scale of 1.0 **without** the constructor being called.
`scale` and `provenance` are `pub` fields, which `charged_error`'s own `# Panics` section already
notes can be bypassed.

This matters beyond pedantry: the file's projection is the very path the new `validate` rung
guards, so the sentence is wrong about the one route the step exists to protect. (The charitable
reading — "the only thing that *reads the constant*" — is true, and is what the sentence should
say.)

### Major 4 — eight of the module's doc links do not resolve; it is the only file in `parameters_file/` with any

Measured with
`cargo doc --no-deps --document-private-items` in the worktree. `defaults.rs` is the single
worst file in the crate at **8** unresolved links; no other file under
`src/ng/calling/parameters_file/` produces one.

```
src/ng/calling/parameters_file/defaults.rs:18:55
18 | //! | the base-quality multiplier, per read group | [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`], one | …
   |    no item named `DEFAULT_ERROR_PROBABILITY_MULTIPLIER` in scope
```

and the same for `DEFAULT_OUTLIER_WEIGHT` (:19, :57), `STATED_FLAT_CONCENTRATION` (:20, :59),
`ReadGroupCalibration` (:55) and `StratumFits` (:59). The cause is that the module has no
module-level `use` at all — every import sits inside `#[cfg(test)] mod tests` — so nothing named in
the header is in scope for rustdoc. Only the one link written as a full path
(`crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight`, :58) resolves.

The plain `cargo doc --no-deps` gate does **not** catch this, because `mod defaults;` is private
and rustdoc does not document private items by default. Both baselines therefore hold (item 5,
below) while eight links in the new file are dead. This is the same evidence as finding 4 in
question 4: a doc-only private module cannot carry short-path intra-doc links.

---

## Minor

### Minor 1 — "Spec §8 sorts them by what a default *means*" attributes this module's taxonomy to the spec

`defaults.rs:13`. §8's own sort has three cases and they are not these:

> - **has one**: the base-quality calibration scale (1.0 — no calibration), the repeat-tract
>   outlier weight (0.01, inherited, §3.8), the flat concentration (1.0). All marked `Defaulted`.
> - **absence is the default**: contamination.
> - **has one, and it has to be measured before it exists** …

§8 puts all three constants in **one** case. The module's three kinds — *asserts nothing* /
*inherited guess* / *a model state* — split that first case in two, which is the module's own
reading (defensible, and §8's closing paragraph gestures at it) but is not "how §8 sorts them".
Same problem with the citation at `likelihood/mod.rs:218`, "(spec §8, **first bullet**)": the first
bullet is exactly the one that declines to distinguish the multiplier from the outlier weight.

### Minor 2 — "inherited" is false for the flat concentration, and §3.8 covers only the outlier weight

`defaults.rs:19-20, 27-31`: "**Two are inherited guesses** … the outlier weight is production's
0.01 … and the fallback concentration is one chromosome's worth of belief spread flat … which is
exactly why spec §3.8 writes **them** into a file a person can edit".

Two problems.

- `STATED_FLAT_CONCENTRATION`'s own origin (`stratum_fits.rs:330-351`) records no inheritance
  from anywhere: it is "one chromosome's worth of belief, spread flat … the same quantity and the
  same reading `ALPHA_REF` carries on the ordinary-site path". That is a stated uninformative
  prior, chosen here — structurally the same *kind* of thing as the multiplier of one, and
  emphatically not "production's 0.01". The table's `inherited guess` cell says something about
  its provenance that is not true.
- §3.8's subject is singular and explicit: "**The one that exists today is** the share of
  repeat-tract reads that came from nowhere the model can explain: `DEFAULT_OUTLIER_WEIGHT`". The
  flat concentration is written by §3.7 ("the run's stated concentration — the bottom rung"), not
  by §3.8.

### Minor 3 — the `to_bits` justification is not true of `f64`

`defaults.rs:91`: "Compared on `to_bits`, because a multiply by a scale that is *near* one would
pass an `==` on some inputs and not others."

For ordinary finite non-zero doubles, `a == b` and `a.to_bits() == b.to_bits()` give the same
verdict; they differ only on `NaN` and on `+0.0` against `-0.0`. A scale near-but-not-equal to one
that rounded back to the same double would pass **both**, and one that did not would fail **both**.
`to_bits` is the right habit here (it is what C3 established for the float round trip) but the
reason given for it is wrong.

### Minor 4 — "All three refusals say the same thing … one shape of message rather than three" is false

`defaults.rs:187-188`. The three messages the test exercises:

- multiplier (`validate.rs:493`): "is `1.25`, and its warrant is `defaulted`, which says no rate
  was fitted for this read group … which is a multiplier of `1.0`; a number you changed is one the
  run was handed, so change the warrant beside it to `supplied`"
- outlier weight (`validate.rs:834`): "is 0.02, and its warrant is `defaulted`, which says this run
  inherited the compiled-in 0.01; a number you changed …"
- concentration (`validate.rs:696`): "is `defaulted` at `3.5`, and the constant a run falls back to
  is `1.0`; a number somebody chose is `supplied`"

The third has a different clause order, a different closing sentence, and does not tell the reader
what to change. What the test actually pins is much weaker and is fine: each names its key, quotes
its constant, and contains the word `supplied`. The doc comment should claim that, not one shape of
message.

### Minor 5 — the step Debug-formats two of the three constants and leaves the third on `Display`

The patch changed the concentration refusal from `{value}` / `{STATED_FLAT_CONCENTRATION}` to
`{value:?}` / `{STATED_FLAT_CONCENTRATION:?}` and wrote the new multiplier refusal with `{:?}`,
under a stated principle (`defaults.rs:242`): "**Spelled as the file spells it**, which is `Debug`
for a float — the writer formats every value with it". That principle is right — the writer is
`to_toml.rs:627`, `format!("{value:?}")`.

The outlier weight's refusal (`validate.rs:834-838`) was left on `Display` for both the value and
the constant. The test passes only because `0.01` prints identically either way; a user who edited
the weight to `1e-7` would be told "is 0.0000001" while their file contains `1e-7`. Either bring it
into line or say why it is left out.

### Minor 6 — the contamination test restates two existing tests

`defaults.rs:164-180`. Its three assertions are already covered:

- the projection to an uncontaminated run — `to_run_parameters.rs:2117`,
  `an_absent_contamination_table_is_not_a_table_of_zeros` (C5's row 1), which additionally shows
  what the run *reports* and that the longhand form is refused;
- `!to_toml().contains("[contamination]")` — `to_toml.rs:1755-1759`, verbatim, and again at
  `mod.rs:2114-2118` through the parsed document.

The assertion is not vacuous (`section()` at `to_toml.rs:383` does write `[contamination]`
literally), but the test adds nothing the suite did not already hold, and its doc comment claims a
role — the §5 row-1 fixture — that C5's test already occupies.

### Minor 7 — "spec §1.2 goal 3 invites — copy the file your run wrote and change one line"

`defaults.rs:186`, and the same phrasing at `validate.rs:474`. Goal 3 is "A person can read it and
change one line … should not need a tool, a schema, or this document." The *copy the file your run
just wrote* half is §7's third bullet ("Editing starts from something"). Minor because the two
together do say what the comment says; worth pinning to §7 since that is where the sentence is.

---

## 2. Does the new `validate` rung belong to E1, or is it scope creep into C2?

**Verdict: it belongs to E1.** The step should keep it, and say in the commit message that it does.

**The case for C2.** C2 is the step that introduced `validate` and, by the owner's decision of
2026-08-28, "**also owns refusing a file that parses and means nothing**". The plan enumerates what
that covers — ranges, an empty sample list, a contamination table with no measured row, a
measurement with two zero counts — and a defaulted multiplier that is not 1.0 is not on the list.
C2 is marked ✅ and its checkpoint has been passed; adding a refusal to it reopens a closed step,
and a reader tracing "why does this refusal exist" from the plan will not find it under E1 either.

**The case for E1, which is stronger.** Three reasons.

1. **E1's own deliverable is unenforceable without it.** The step is "Each marked `Defaulted` when
   used", and C5 established that the mark is the *whole* of the difference — a defaulted 1.0 and a
   fitted 1.0 multiply every read by the same number. A mark nothing checks on the way in is a mark
   a hand-edit silently detaches from its number.
2. **The rung cannot precede the constant.** C2 had no `DEFAULT_ERROR_PROBABILITY_MULTIPLIER` to
   compare against; E1 is the step that names it. A refusal written at C2 would have had to inline
   `1.0`, which is the second spelling the header rightly objects to.
3. **The precedent is already established, and it is not C2.** The identical rung for the outlier
   weight landed in `dc1e1fd2` ("an edited outlier weight now reaches the score") and the one for
   the flat concentration in `21adc757` ("the tract ladder's bottom rung carries its own warrant") —
   neither under C2, both alongside the work that made the constant reachable. Verified with
   `git show <commit> -- src/ng/calling/parameters_file/validate.rs`. The multiplier's rung landing
   with E1 completes a set of three the same way the other two arrived.

The change is small (24 lines in `validate.rs`), sits inside a function C2 already wrote, and adds
no new refusal *category* — it is the third instance of one C2 shipped. That is completion, not
creep.

---

## 3. Is the new module's placement right?

**Verdict: no. The header should move into `parameters_file/mod.rs`, and each test to the module
that owns the behaviour it pins. `defaults.rs` as landed does not earn a module.**

The module contains no code — only `//!` documentation and `#[cfg(test)] mod tests`. Three
consequences, and the first two are measured:

- **Its documentation is invisible.** `parameters_file` is `pub mod`
  (`src/ng/calling/mod.rs:65`); `defaults` inside it is private, so `cargo doc --no-deps` never
  renders the header. The "one place" it claims to be is a place no reader of the generated docs
  can reach.
- **Its links cannot resolve** (Major 4). Eight of nine are broken, and they are broken *because*
  the module has no code: a module-level `use` to fix them would be an unused import in every
  non-test build. A doc-only module and short-path intra-doc links are incompatible.
- **Its tests are homeless.** Four of the six assert on other modules' types
  (`ReadGroupCalibration`, `RepeatTractOutlierWeight`, `StratumFits`) or call
  `ParametersFile::validate`; the two `validate` tests directly mirror
  `an_outlier_weight_whose_warrant_no_run_could_mean_is_refused`, which lives in `validate.rs`
  beside the code it checks. The contamination one duplicates two existing tests (Minor 6).

The header itself is worth keeping — the four-way table and the "where each constant lives" list
are the useful part of E1, and `mod.rs`'s header is already the module's map, is public, and has
the imports to make the links resolve. The tests belong beside the three checks in `validate.rs`
and beside `ReadGroupCalibration` in `likelihood/mod.rs`.

If the module is kept anyway, it needs (a) full paths in every link, and (b) an explanation of why
a private doc-only module is the right home for a summary the docs will not show.

---

## 4. `cargo doc --no-deps`: no new diagnostic

Counted in `/Users/jose/devel/pop_var_caller-e1-rev2` with the worktree's own `dev.sh`, patch
reverted then re-applied, both runs logging `Documenting pop_var_caller v0.1.0
(/Users/jose/devel/pop_var_caller-e1-rev2)`:

| | `error: unresolved link` | `warning: redundant explicit link target` | total `^error`/`^warning` lines |
|---|---|---|---|
| baseline (`8877316f`, patch reverted) | 25 | 23 | 51 |
| with the patch | 25 | 23 | 51 |

Identical, line for line (`diff` of the two diagnostic lists is empty). **The step adds nothing to
either baseline.** The eight broken links in Major 4 sit below this gate because `mod defaults;` is
private; they appear only under `--document-private-items`, where the crate carries a pre-existing
108 and `defaults.rs` contributes the largest single share of any file.

`cargo test --lib ng::calling::parameters_file`: 176 passed, 0 failed, 2 ignored — including the six
new `defaults::tests`. `cargo clippy --lib --all-targets` adds no warning naming any file in the
patch.

---

## What is right about the step

Worth saying, because the findings above are all about prose:

- The constant is placed correctly — beside `ReadGroupCalibration`, not re-declared in the new
  module — and `defaulted()` now reads it, so the two spellings the header warns about do not
  exist.
- The `validate` rung is the right shape and closes a real hole: the pre-existing guard
  (`a_warranted_value`, `validate.rs:1071`) only fires when a `defaulted` value carries
  `observations`, and a defaults run writes its multiplier with none — so before this patch the
  edit the step is worried about passed cleanly.
- The three claims most at risk of being wrong are right: `charged_error` really is
  `(scale · exp(q_sum/n)).max(MIN_BASE_ERROR)` (`likelihood/mod.rs:479`); the fixture's three reads
  really do span Phred 40 to Phred 13; `repeat_tract_parameters` really does give a cell with no
  slippage row `StutterModel::hipstr_shipped` and `Provenance::Defaulted` (`:347`) and count them
  (`:209`).
- The `{:?}` change to the concentration refusal is a genuine fix, matching the writer's own
  formatter at `to_toml.rs:627`.
