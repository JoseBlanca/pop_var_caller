# ng parameters file — the tract ladder's bottom rung carries its own warrant

**Date:** 2026-08-30
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone D — **its own
commit, before D3**, on the owner's ruling, the way the outlier weight's wiring landed before C2's
projection
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §2.1, §3.7, §13 test 5
**Code:** `StratumFits` in
[stratum_fits.rs](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs), and the three
places in `parameters_file/` that read or write it

---

## 1. The hole this closes

The tract ladder's bottom rung states one concentration for the whole run: the median of the
concentrations this run's strata fitted, or `STATED_FLAT_CONCENTRATION` where it fitted none. The
file writes it with a warrant beside it — and **the warrant was not carried anywhere**. The
projection *in* read only the number, and the writer worked the warrant out again from whether any
stratum was fitted.

That is fine for the two states a fit can be in and wrong for the one a *file* can be in. **Spec
§2.1 demotes every number of a file fitted under another census to `supplied`, wholesale**, and
§13's fifth test is *same calls, every warrant `Supplied`*. A demoted file written back out
re-emerged as `fitted_here`, because a run that took the median itself and a run handed one look
identical from the strata beside them — which is the whole of what a warrant is for.

**Owner's ruling of 2026-08-30: carry the warrant through `StratumFits`**, and do all of it — the
writer stops re-deriving, and `validate`'s rule keyed on the re-derivation moves or goes.

## 2. What moved

**`StratumFits` gained `stated_concentration_warrant: Provenance`**, beside the number it is about.
`over` sets it from whether the fit produced any stratum spectrum, which is what the writer used to
do afterwards; `of_gathered_rows` takes it as a sixth argument, so a file's own warrant travels in.

**The writer copies it** rather than deriving it. **The reader passes the file's through.**

## 3. `validate`'s rule: what survived and what could not

The old rule was `fitted_here` **exactly** where the file names a fitted stratum spectrum. Half of
that is still true and half of it refused the file spec §2.1 produces:

- **`fitted_here` still requires a fitted stratum spectrum** — a file naming none has nothing for
  the number to be the median of.
- **`defaulted` now requires the value to be `STATED_FLAT_CONCENTRATION`**, which is the rule the
  outlier weight already carries, for the same reason: *defaulted* is a claim about which constant
  was used.
- **`supplied` and `borrowed` are free.** That is the point. Nothing about the strata beside the
  number can say whether it was handed over, and the old rule refused every demoted file.

`the_warrant_a_demoted_file_carries_on_the_bottom_rung_is_accepted` is the case that could not pass
before.

## 4. Two fixtures were wrong about what a run produces

Both "smallest run" fixtures cleared the strata and set the warrant to `defaulted` while leaving
the value at **3.5**. A run that fitted no stratum states the compiled-in constant — `validate`'s
own copy of the fixture already did this. Both now do.

## 5. The file's own prose said the thing that stopped being true

The header told a reader that this key's warrant "is decided by whether this file holds any fitted
stratum spectrum". After this change it is not decided by anything the file holds: a demoted file
says `supplied`. Both notes rewritten — the header's rule, and the section note, which now says
what a `supplied` bottom rung means and tells an editor to write it.

**One sentence in the same paragraph was corrected while it was open**, from the geneticist's read
of D1: *"a larger number moves the prior less"* had no subject and read backwards — what moves less
is the reads' pull on the prior. It is *"a larger number makes the prior harder for the reads to
move"* now.

## 6. Tests

Two added:

| test | what it holds |
|---|---|
| `the_bottom_rungs_warrant_is_the_files_and_not_worked_out_from_its_strata` | three warrants through file → run → file, **with the fitted strata left in place** — which is what the existing round trip cannot see, since its fixture is `fitted_here` beside fitted strata and a re-derivation gives the right answer for the wrong reason |
| `the_warrant_a_demoted_file_carries_on_the_bottom_rung_is_accepted` | `supplied` and `borrowed`, with and without fitted strata beside them |

## 7. Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib ng::calling::parameters_file`: **160 passed, 0 failed** (158 before).
- `cargo test --lib`: **5,531 passed, 0 failed** (5,529 before).
- `cargo doc --no-deps`: **25 unresolved-link errors**, the baseline.
- The golden file regenerated and the diff read: the two rewritten notes, and nothing else.

## 8. Mutation testing

**Five mutants, every one fails a test** — and the first is the behaviour this commit replaces:

| mutant | tests it fails |
|---|---|
| the writer re-derives the warrant from the strata (the old behaviour) | 1 |
| the reader drops the file's warrant and says `fitted_here` | 2 |
| the fit marks the flat constant `fitted_here` too | 2 |
| `validate` stops refusing `fitted_here` with no fitted stratum | 1 |
| `validate` stops refusing `defaulted` away from the constant | 1 |
