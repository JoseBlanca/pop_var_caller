# ng parameters file — F1: reading the artefact as a geneticist

**Date:** 2026-08-31
**Step:** [parameters_file.md](../../ng/impl_plan/parameters_file.md) F1
**How it was run:** one agent, told to be a geneticist and not a Rust reviewer, in a worktree
detached at `ede29317` with the step applied. Nothing in the tree writes a parameters file to disk
in an ordinary run, so it **generated five** by adding a temporary `#[ignore]`d test and removing it
afterwards (verified gone).

**This seat has now found the most valuable defects on every step of Milestones E and F**, and F1
is no exception: three Blockers, none of which a compiler could see, and one of them a sentence
this step introduced.

---

## The files it read, and its own count of the seven groups

| file | how | its count | the file's claim |
|---|---|---|---|
| fitted everything | the golden `every_shape_as_written.toml` | 7 of 7 | ✓ |
| fitted nothing | real `RunParameters::of_defaults`, 3 lanes / 2 samples | 0 of 7 | ✓ |
| partly fitted | every-shape with contamination and substitution rates blanked | 5 of 7 | ✓ |
| one sample, one read group | real `of_defaults` | 0 of 7 | ✓ |
| a defaulted multiplier of 4.0 | every-shape, read group 1's multiplier set to 4.0 | 7 of 7 | ✓ |

**Every headline count is arithmetically right**, checked group by group by hand. The defects are
elsewhere.

---

## Blocker — the file told the reader their reads were left alone, above a number that made them four times worse

```
# not calibrated: this read group's reported qualities are taken at face
# value, because no usable error rate could be fitted for it
{ read_group = 1, error_probability_multiplier = { value = 4.0, warrant = "defaulted" } },
```

The origin comment is attached to the **warrant** and never to the value, and by the owner's ruling
of 2026-08-31 a read group nothing could be fitted for is *charged the pre-pass's stated rate*
rather than believed — so the multiplier is that rate over the library's own mean reported error
and 4.0 is an ordinary value (`a_run_whose_rates_were_defaulted_writes_a_file_its_own_reader_accepts`
pins 4.0 and 2.0 for libraries reporting 2.5 × 10⁻⁴ and 5 × 10⁻⁴). A reader concludes lib4 was
scored at the quality the instrument reported; it was scored **six Phred less confident on every
base**, which is the difference between a het call surviving and not.

The file's own header made the same claim in the other direction — *"the base-quality multiplier at
1.0"* — so a reader meeting `value = 4.0, warrant = "defaulted"` would take the file for corrupt.

**Fixed.** The section carries the explanation once (a run of thousands of read groups must not pay
for it a row at a time), the row note is one line pointing at it, and the header no longer lists the
multiplier among the numbers with a built-in value.

## Blocker — the four warrant words were never defined, and `borrowed` is the one the headline rests on

The vocabulary appears once, as a list — *"a `warrant`: fitted_here, borrowed, supplied or
defaulted"* — and `borrowed` then appears as a value twice and is explained nowhere. *Borrowed from
what?* A reader is left with a coefficient of 0.17 for a named plant, 9,411,027 covered positions
beside it, and no way to know whether 0.17 describes that plant or its lane-mates.

**This mattered more after F1 than before it**, because the new headline counts `borrowed` as
fitted. §1.2 goal 3 says a person should not need a tool or the spec; on the one word carrying the
headline, they did.

**Fixed.** All four are defined where the headline introduces them.

## Blocker — the new empty-census note promised a demotion that does not happen

```
# **This run had no census** … a run that reads this file will find a
# disagreement at the first line and report every number in it as supplied
```

`census_disagreement` zips the two term lists and then looks for a surplus on either side. **Two
empty lists agree.** And a run with no census passes `None` and compares nothing. So the sentence
was false for exactly the reader the line above it names — a defaults run, or **any direct-mode
run**, which is the whole reason this file format exists.

`bindings.rs`'s own doc scopes it correctly (*"so a **psp-mode** run reading this file…"*); the
file's prose dropped the scope. **Fixed**, and the note now gives both answers and says which run
gets which.

## Major — on the file making the strongest claim, the reader could not check it

*"All 7 groups … were fitted"* and then no list. Counting `[section]` headings gives **nine**, and
nothing said that two are outside the denominator or that `repeat_tracts` counts as three. **Fixed**
— the all-fitted arm names the seven too.

## Major — the file stated a rule about zeros that its own defaults run breaks

*"a zero means it was measured and found to be zero"*, in the position of greatest authority, a
hundred lines above `{ sample = "TS-1", inbreeding_coefficient = { value = 0.0, warrant =
"defaulted" } }`. **Fixed** — the rule now says *a zero under a `fitted_here` or `borrowed` warrant*,
and names the counter-example.

## Major — the contamination note lands inside `[base_quality_calibration]`

Content right, address wrong: a reader scanning section by section reads it as a remark about
calibration. TOML cannot head an absent section, so **a blank line was added** to detach it. Not a
full fix.

## Major, and NOT fixed — the fallback is disclosed only when a table is *completely* empty

`every_shape_as_written.toml` names 3 read groups and 3 strata, writes **one** substitution-rate
row, and its section never mentions a fallback: eight of nine cells silently take the stated 0.001.
Same for slippage. The reader concludes their repeat tracts were fitted.

**E3 decided the opposite deliberately** — a note saying *some tracts fell back* would be true of
almost every run and would tell a reader nothing about theirs — and recorded why. **Raised at
Checkpoint F rather than overturned here.** The reader's counter-argument is that the partly-covered
case is the common one and is the one with no note.

## Major, and NOT fixed — the sections' prose describes keys that are not there

A defaults run's `[repeat_tracts]` spends fifty lines defining `share_of_reads_that_slip_origin`,
`curve_fitted_on`, blends and curves, none of which appears in a file whose three tract tables are
all `[]`. The new headline is what makes this survivable. Conditional section prose is beyond F1.

## Major, and the same gap the other two reviewers found — "fitted from your data" cannot tell this run's fit from a previous run's

**Fixed**, by the same change: the headline says *fitted from reads* and a derived sentence names
the three groups that can say **whose** and the four that cannot.

---

## Minors, recorded and not fixed

`format_version` and `ploidy` carry no comment; the calibration rows are the only per-read-group
rows without the library name, so acting on a library means mapping by hand through
`[fitted_from].read_groups`; the prior's two concentrations have no unit where the repeat-tract one
does; the slippage curve's own vocabulary (`shape`, `intercept`, `slope`, `bend`, `centre_repeats`,
`held_out_error`, `strata`, `cells`, `rise_shape`, `curve_weight`, `reach`) is undefined;
`fitted_from_reads_of`'s two values are never explained; no legal range is stated for the
multiplier; nothing tells an editor that the stale origin comment above a row they changed becomes
false; and the one-sample file says nothing about being one sample — it still offers `borrowed`,
impossible with one read group and one sample.

## The goal-3 test: could it raise one library's error rate?

**Yes, from the file alone**, with one lookup and three things the file did not say: which read
group is lib3 (the row does not name the library); that deleting the stale origin comment is the
reader's job; and that forgetting to change the warrant is not caught — `{ value = 2.0, warrant =
"defaulted" }` is accepted and the run then reports a typed number as a compiled-in default.

## What it could not explain without reading the Rust

The seven groups; what `borrowed` borrows from; why a `defaulted` multiplier need not be 1.0; that
two empty censuses agree; and every field of the slippage curve's inline tables. **The first four
are fixed.**
