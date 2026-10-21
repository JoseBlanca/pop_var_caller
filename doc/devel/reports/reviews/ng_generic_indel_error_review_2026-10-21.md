# Code Review: generic-indel-error (commit f800df67)
**Date:** 2026-10-21
**Reviewer:** rust-code-review skill, single reviewer (no category fan-out; see Scope)
**Scope:** f800df67 — on the SNP/indel path, a wrong read that only an indel could explain is charged without the read group's calibration scale
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the PR-sized diff of commit f800df67 on branch `generic-indel-error`. I
  also read the net diff of `src/calling/likelihood/generic.rs` against `main`, because f800df67
  partly undoes the research commit 18dd5712. I checked the change's callers and callees across
  the whole crate.
- **Reviewed against:** f800df67, with HEAD at 9f51ac0e (a merge from main that does not touch
  these files).
- **In-scope files:** `src/calling/likelihood/generic.rs`, `src/calling/likelihood/mod.rs`,
  `src/calling/mod.rs`, `src/calling/inference/summarise_condition.rs`, `src/run/records.rs`,
  `doc/devel/ng/spec/read_likelihoods.md` (§3.3, §3.5, §5.3, §10, §12 item 20).
- **Also checked for consumers:** `src/calling/inference/`, `src/run/` (including `explain.rs`),
  `src/calling/quality/`, the `*_parity.rs` files, `src/parameter_estimation/`, `tests/`,
  `examples/`.
- **Deliberately out of scope:** the research scripts under `benchmarks/giab/src/`. The report
  `doc/devel/reports/ng_generic_indel_error_2026-10-21.md` is used only as the source of the
  numbers.
- **Deviation from the skill's procedure, stated plainly.** The dispatching agent asked for a
  focused correctness review, delivered to `tmp/review/`, with no source edits. So:
  - I did not fan out one sub-agent per category in its own worktree;
  - I did not save the report to `doc/devel/reports/reviews/` or update `PROJECT_STATUS.md`;
  - I did not run mutations, because each one needs a source edit. Each mutation result below is
    argued from the fixture's arithmetic, not measured, and says so.

  The categories covered here are reliability, errors, naming, refactor_safety, float_portability
  and extras ("does the diff match its stated intent").

## 2. Verdict

**Approve-with-changes.** I found no correctness defect. The code does what the amended spec says
in every case I checked:

- **the "every carried allele" rule** — the cell is `OneIndelAway` only when every carried allele
  is one indel away;
- **the mixed case** — one carried allele a substitution away and another an indel away is
  charged the scaled error;
- **contamination** — only this individual's share of the charge drops the scale; the
  contaminant's term `c·q(o)` is untouched;
- **partial reads** keep the scale;
- **the `MIN_BASE_ERROR` floor** is applied identically to the scaled and unscaled charge;
- **ploidy above 2 and loci with more than two alleles** — the fill loop does not depend on
  either.

What needs changing:

- two decisions the spec states explicitly have no test that would fail if they were reversed (M1);
- several doc comments, the spec's own §3.2, and the report's status note were not updated to
  match (Mi1 to Mi3).

## 3. Execution status

- `cargo test --release -j 2 --lib calling::likelihood -- --test-threads=2` — exit 0,
  `test result: ok. 274 passed; 0 failed; 0 ignored; 0 measured; 4969 filtered out`.
- **Not run:**
  - `cargo fmt --check`, `cargo clippy`, `cargo doc` — the author reports no clippy warnings in
    the touched files;
  - the full `calling::` suite and the byte-identity checks on GIAB and tomato — the author
    reports both;
  - mutation runs — they need source edits, which this review was told not to make.
- Findings labelled "Needs verification": 1 (Mi4, whether a scale below one occurs on real data).

## 4. Open questions and assumptions

1. **Can a fitted scale be below one, and is the amendment wanted there?**
   `ReadGroupCalibration::from_fitted_rate` (`src/calling/likelihood/mod.rs:397-402`) puts no
   lower bound on the scale, and the tests use scales of 0.37 and 0.5. Where the scale is below
   one, dropping it charges an indel-explained wrong read a *larger* error probability than the
   fit's scale would, which makes that read cheaper to explain as an error and so pushes the
   caller towards explaining the indel away. That is
   the opposite of what the amendment's "Why" paragraph argues for. It is still the `--defaults`
   behaviour, so it may well be acceptable, but the spec should say which. Affects Mi4.
2. **Assumption: a candidate allele may be the empty sequence** (a deletion projected without an
   anchor base). `differ_by_one_indel` handles that correctly — `""` against `"ACG"` is one
   indel — but no test pins it. Affects the Nits.

## 5. Top 3 priorities

1. **M1: pin the two untested decisions.** The partial-read rule (partials keep the scale) and the
   heterozygote of two indel alleles each survive a mutation that reverses them.
2. **Mi1: finish the doc trail.** Spec §3.2 states the scale applies to every read and has no
   pointer to the exception. Three doc comments in `generic.rs` still describe a wrong read's
   charge as always scaled.
3. **Mi2: correct the report's status note.** It says the owner ruled for the report's §7 text,
   but that text makes partial reads unscaled too, and what was built keeps their scale.

## 6. Findings

### Blocker

None.

### Major

#### M1: src/calling/likelihood/generic.rs:1865, :1909, :3662 — two decisions the spec states have no test that would fail if reversed
**Confidence:** High for the partial-read case; Medium for the heterozygote case (argued from the
fixtures, not run).

**(a) Partial reads keep the scale** (spec §5.3, new paragraph; the commit message says the
same). Every partial-read test builds its row through `row_with_partials`, which passes
`uncalibrated()`, a scale of 1 (`generic.rs:3662`, call at `:3672`). At a scale of 1 the
scaled and unscaled charges are bit-identical — that is what
`at_a_scale_of_one_the_minted_error_is_the_charged_error` asserts. So a change that charges an
unexplained partial its unscaled minted error would pass every test in the file. For example,
`fill_generic_emissions` (`:849-855`) could be changed to divide `minted_error` instead of
`charged_error`. Before this commit nobody would write that line. Now a sibling buffer of exactly
that type sits beside it, so it is a one-identifier slip, and the spec presents this case as a
deliberate decision.

**(b) A heterozygote whose two alleles are both one indel away.** The tests reach the
`OneIndelAway` cell only through homozygotes (`hom_ref` and `hom_insertion` at `:1865-1900` and
`:1909-1935`). Suppose the fill were written as "`OneIndelAway` only when the genotype carries a
single distinct allele". Then:
- the expected values in both tests stay the same;
- the existing mixed-case assertion still sees `OtherwiseDifferent`.

So that mutation survives, and so does any variant that treats ploidy above 2 differently. The
spec's rule ("when **every** allele a genotype carries…") covers this case explicitly, and it is
the case the word "every" exists for.

**Fix.** Add the two tests in §8:
- `a_partial_no_allele_fits_keeps_the_scale_at_a_scale_of_four`;
- `a_heterozygote_of_two_indels_charges_a_reference_read_unscaled`, with a tetraploid row.

### Minor

#### Mi1: doc/devel/ng/spec/read_likelihoods.md:410-416, :668 and src/calling/likelihood/generic.rs:652, :663, :782 — prose still says every wrong read is scaled
**Confidence:** High.

The amendment added notes to §3.3, §3.5, §5.3, §10 and §12. Five places that state the old rule
were not touched:

- **Spec §3.2** ("The decision: keep the read's own number, set its scale from the fit",
  `:410-416`) still gives `ε_read(after) = ε_read(as minted) × scale(read group)` for every read.
  §3.2 is the section a reader goes to for what the scale applies to.
- **The comparison table at spec `:668`** still charges a read the genotype cannot explain
  `… + log scale − log m`.
- **`genotype_log_likelihood_row`'s doc:**
  - `generic.rs:652` defines the error branch as `ε̄/m`;
  - `generic.rs:663` says §3.3 charges a wrong read `q_sum + n·(log scale − log m)`.
- **`fill_generic_emissions`' doc** (`generic.rs:782-786`) lists what the emission cache holds as
  "the charged error `ε̄` … scaled by the read group's calibration" and the partial verdicts. It
  omits the new unscaled minted error.

A reader who arrives at either function or at §3.2 learns the pre-amendment rule.

**Fix.**
- **§3.2:** after the formula, add one sentence: "Except for an error only an indel could have
  made, which §3.5's 2026-10-21 amendment charges unscaled."
- **Spec `:668`:** add "(no `log scale` for an indel error, §3.5)".
- **`generic.rs:652` and `:663`:** add "— without `scale` where every carried allele is one indel
  away ([`WrongReadCharge::OneIndelAway`])".
- **`fill_generic_emissions`:** add a bullet: "per complete observation, the same `ε̄` without
  the scale — what an indel error is charged".

#### Mi2: doc/devel/reports/ng_generic_indel_error_2026-10-21.md:7-12 vs :407 — the status note says the owner ruled for text that differs from what was built
**Confidence:** High.

The status note says "the owner ruled for §7's amendment. It is in `read_likelihoods.md` §3.5 and
built". §7's proposed text (`:407-408`) says "The same holds for a partial read that no allele
fits, when every carried allele is an indel away". The built rule and the adopted spec (§5.3)
keep the scale for partial reads. The report's §4 (`:297-299`) also recommends that partials be
treated the same way.

The "byte-identical to `nomult`" claim is still true, because the prototype did not cover partial
reads either. But the report now tells a reader that partials are unscaled, and they are not.

**Fix.** In the status note, add: "except that partial reads keep the scale (spec §5.3: deciding
that a fragment is one indel from an allele would need an alignment the partial rule avoids, and
complete reads alone restored every lost indel)."

#### Mi3: src/calling/likelihood/generic.rs:2384 — the production differential's reconciliation no longer accounts for every difference
**Confidence:** High that it is latent. It does not fire today.

The test's own doc says "every difference attributed, none unexplained" and lists three
differences. `scale_effect` takes `n·log scale` back out of every wrong read unconditionally. The
fixture holds only the alleles `A`, `C`, `G`, so no cell is `OneIndelAway` and the test passes. If
anyone extends the fixture with an indel allele, the reconciliation will be off by
`n·log 2.5` per indel-explained read. That person will be debugging the test, not the code. The
test now differs from production in a fourth way that it neither states nor reconciles.

**Fix.** Make it consult the cell, as `plain_log_likelihood_row` (`:2975-2982`) already does:

```rust
let scale_effect: f64 = supported
    .iter()
    .filter(|o| counts[usize::from(o.allele.get())] == 0)
    .filter(|o| spread_table.at(GenotypeIdx(genotype as u32), o.allele).is_scaled())
    .map(|o| f64::from(o.num_reads) * log_scale)
    .sum();
```

Also add a line to the doc: "and §3.5's 2026-10-21 amendment drops the scale for indel errors,
which this SNP-only fixture never reaches."

#### Mi4: doc/devel/ng/spec/read_likelihoods.md:782-817 — the amendment's argument assumes the scale is above one
**Confidence:** Low; Needs verification (Open question 1).

"Why" argues that the scale inflates indel errors, from HG002's fitted 6.2 to 14. The rule
itself, though, applies at any scale. Below one it raises an indel error's charge above what the
fit set for substitutions. That moves calls in the opposite direction from the motivation: a
fitted run there would explain indels away more readily than with the scale applied.

**Fix.** Two options:
- add one sentence to "Ends of the range" saying the rule is meant at scales below one too (the
  charge falls back to the `--defaults` one), if that is the intent;
- or restrict the exemption to `scale > 1`, i.e. charge `min(scaled, unscaled)`.

The first is the smaller change and keeps the rule free of anything estimated.

#### Mi5: src/calling/likelihood/generic.rs:155, :232-340, src/calling/mod.rs:1510 — "spread" names now hold a charge class
**Confidence:** Medium.

The cells are `WrongReadCharge` values, and they answer two questions: the divisor and whether the
scale applies. The surrounding names still say they hold a spread:
- `fill_error_spreads`, `ErrorSpreadTable`, `spreads_for`, `spreads_of`;
- the scratch field `error_spreads` and its accessor `error_spreads_mut`;
- the loop variable `spreads` (`:1026`).

The project's rule is that names say what the value is. The type doc already says the table
"stores [`WrongReadCharge`] cells", which is the giveaway.

**Fix.** Either rename (`fill_wrong_read_charges`, `WrongReadChargeTable`, `charges_for`,
`charges_of`), or, if the churn is unwanted now, leave the names and record that in the
`ErrorSpreadTable` doc. Renaming is mechanical, and the compiler catches every site.

### Nits

- **`generic.rs:48`:** "on ordinary ground" — the report and spec say "generic ground" or "invariant
  GIAB ground", and the 0.4-to-13 range includes repeat contexts, which "ordinary" suggests are
  excluded. "The parameters fit's error rate" (`:45`) reads awkwardly; "the parameter pre-pass's
  fitted error rate" matches the spec.
- **`src/calling/mod.rs:1760-1767`:** the scratch lost its unwritten-cell sentinel. Previously an
  unfilled cell was `NaN` and poisoned the row loudly; now it reads as `OtherwiseDifferent`, a
  plausible answer. The fill's exact-length assertion and the two-poison test
  (`every_cell_is_written_over_whatever_the_buffer_held`) cover the fill itself, so this is
  acceptable. The comment states it honestly.
- **`differ_by_one_indel` tests (`generic.rs:1841`)** cover the prefix case (`ACGT`/`ACG`), the
  suffix case (`ACGT`/`TACGT`), repeats and two runs. They do not cover an empty allele (`""`
  against `"A"` should be `true`; `""` against `""` should be `false`). Worth one line each;
  see Open question 2.
- **`ReadGroupCalibration::minted_error` (`likelihood/mod.rs:535`)** repeats `charged_error`'s two
  input assertions and its `exp`. The fill now takes two `exp`s per complete observation instead
  of one. The cost is once per sample per locus, not once per pass, so it is negligible.
  `charged_error` could be written as `(self.scale * exp(mean)).max(F)` beside
  `exp(mean).max(F)` from one shared helper, keeping the bit-identity at a scale of one.

## 7. Out of scope observations

- **`src/calling/inference/summarise_condition.rs:1497-1535` (`weakest_warrant_at_the_locus`)**
  folds every read group's calibration provenance into a generic locus's warrant. At a
  two-allele locus of reference and one indel, with no partial reads, the scale now enters no
  genotype's likelihood:
  - each homozygote charges the other allele's reads unscaled;
  - the heterozygote explains every read.

  The call's reported warrant can therefore be weaker than what it rests on. That is the same
  over-report the comment there rejects for repeat tracts. The impact is small (a provenance
  label). Follow-up: decide whether the fold should skip read groups whose only wrong reads
  were indel-explained, or note the approximation in that doc.
- **`plain_log_likelihood_row` (`generic.rs:2950-2990`)**, the independent §3.3 oracle, writes
  the error charge as `q_sum + n·log scale` with no `MIN_BASE_ERROR` floor. The row floors it.
  The floor is unreachable at admitted qualities, as the `charged_error` doc measures, so this is
  a documentation point, not a defect. It predates this change.

## 8. Missing tests to add now

**`score_partials` / `fill_generic_emissions` — partial reads keep the scale (M1a)**

```rust
/// Spec §5.3, as amended 2026-10-21: an unexplained partial keeps the scale even where a
/// complete read showing the same indel would not.
#[test]
fn a_partial_no_allele_fits_keeps_the_scale_at_a_scale_of_four() {
    // Reuse a_partial_no_carried_allele_explains_is_charged_its_own_quality's alleles and
    // partial (3 reads, q_sum = −21), but build the row with `calibrated(4.0)`.
    // Expect hom_ref == 3·(−7 + ln 4), and NOT 3·(−7).
}
```

Bug caught: an unexplained partial charged `minted_error` instead of `charged_error`, i.e. the
report's §7 wording built instead of the spec's.

**`fill_error_spreads` — every carried allele, heterozygote and ploidy 4 (M1b)**

```rust
#[test]
fn a_heterozygote_of_two_indels_charges_a_reference_read_unscaled() {
    // Reference ACGT, insertion ACGGT, deletion AT: each alternative is one indel from the
    // reference and from each other (ACGGT − CGG = AT).
    let alleles = locus(&[b"ACGT", b"ACGGT", b"AT"]);
    let table = diploid(3);
    let out = spreads(&alleles, &table);
    let ins_and_del = genotype_carrying(&table, &[0, 1, 1]);
    assert_eq!(table_over(&out, &table).at(ins_and_del, AlleleId(0)),
               WrongReadCharge::OneIndelAway);
    // And at a tetraploid carrying both, 2 + 2 copies.
    let tetraploid = GenotypeTable::build(Ploidy::try_new(4).unwrap(), 3);
    let mut out4 = vec![WrongReadCharge::default(); tetraploid.view().genotype_count() * 3];
    fill_error_spreads(&alleles, &tetraploid.view(), &mut out4);
    assert_eq!(table_over(&out4, &tetraploid)
                   .at(genotype_carrying(&tetraploid, &[0, 2, 2]), AlleleId(0)),
               WrongReadCharge::OneIndelAway);
}
```

Bug caught: an "indel only at a homozygote" or "only one distinct carried allele" implementation,
which every current test accepts.

**`differ_by_one_indel` — empty alleles (Nit)**

`differ_by_one_indel_treats_an_empty_allele_as_one_deletion`: `assert!(differ_by_one_indel(b"", b"A"))`,
`assert!(!differ_by_one_indel(b"", b""))`.

## 9. What's good

- **The charge is one enum, not two parallel tables** (`WrongReadCharge`, `generic.rs:42-71`). The
  divisor and the scale decision cannot be filled from different loci or read at different strides.
  The doc says this is why.
- **The mixed-case rule is enforced by "all equal to the first", not by an `any`/`all` pair**
  (`generic.rs:206-227`). The same three lines produce both main's old substitution rule and the
  new indel rule, so the `--defaults` byte-identity follows from the structure.
- **`differ_by_one_indel` caps the common suffix at what the common prefix leaves of the shorter
  sequence** (`generic.rs:407-421`). That is what makes `GAAAAT`/`GAAAAAAT` come out right inside
  a repeat, and the test pins exactly that case. Without the cap, prefix 5 plus suffix 5 would be
  10, which is not 6.
- **The two-poison test** (`every_cell_is_written_over_whatever_the_buffer_held`) is the right
  replacement for a `NaN` sentinel on an enum with no invalid value. Its fixture
  (`ACGT`/`ACCT`/`AT`) reaches all three cell values.
- **`minted_error` is checked bit-for-bit against `charged_error` at a scale of one**
  (`at_a_scale_of_one_the_minted_error_is_the_charged_error`). That test, not a tolerance, is what
  guarantees `--defaults` output is unchanged.

**Numbers in the diff, re-derived against the report:**
- "0.4 to 13 in 10,000": report table, 0.4 with no repeat to 12.8 at a 5-copy dinucleotide;
- "scale 6 to 14": 6.2 / 6.4 / 14.1;
- "30 to 95 in 10,000": report summary, 30 to 45 at 30× and 60 to 95 at 300×;
- missed indels 61/50, 26/18, 10/7: report `:258-259`;
- tomato, 46 records and 43 genotypes among 6,888: 36 added + 10 removed, `:270`;
- the §10 slippage figures 19/50/121/203/274 against 22/50/130/194/292;
- "8 fewer missed for 3 more false at 10×, 2 for 1 at 30×";
- "tenfold" (4.4 against 0.4) and "up to twofold" (8.3 against 5.2, and 12.8 against 6.8).

All CHECKED-CORRECT. The test comment's `q + ln 4 − ln 3` is confirmed by the test passing.

## 10. Commands to re-verify

- Ran: `cargo test --release -j 2 --lib calling::likelihood -- --test-threads=2` (274 passed).
- New: after adding §8's tests, the same command, plus a check that each new test fails when its
  target is reversed:
  - in `fill_generic_emissions`, divide `minted_error` instead of `charged_error` for partials;
  - in `fill_error_spreads`, return `first` only when `carried_copies` has a single non-zero
    entry.

---

## 11. Resolution (author, 2026-10-21)

Every finding was applied; nothing was declined.

- **M1** — three tests added: `a_partial_no_allele_fits_keeps_the_scale_where_a_complete_read_drops_it`,
  `a_heterozygote_of_two_indels_charges_a_reference_read_unscaled` (diploid and tetraploid) and
  `an_empty_allele_is_one_deletion_from_any_other`. **Mutations run, not argued**: charging
  unexplained partials `minted_error` fails the first; classing a cell `OneIndelAway` only at one
  distinct carried allele fails the second (and one existing test); dropping the rule in the
  assembly fails three tests. Each mutation was reverted from a saved copy and the file checked
  identical.
- **Mi1** — spec §3.2 gains the exception after its formula; the comparison table notes "no
  `log scale` for an indel error"; the docs of `genotype_log_likelihood_row` and
  `fill_generic_emissions` name the unscaled charge and the cache's new entry.
- **Mi2** — the report's status note says partial reads keep the scale, and why.
- **Mi3** — the production differential's scale correction filters on `is_scaled()`, with a note.
- **Mi4** — spec §3.5 states the rule holds at any scale, because the fitted rate carries no
  information about indels in either direction; every fit seen is above one (2.7 to 19.6).
- **Mi5** — renamed: `fill_wrong_read_charges`, `WrongReadChargeTable`, `charges_for`, `charges_of`,
  the scratch field `wrong_read_charges` and its accessor; prose "error-spread table" became
  "wrong-read charge table"; the architecture doc follows.
- **Nits** — wording fixed; `charged_error` and `minted_error` share one helper,
  `geometric_mean_error`, keeping them bit-identical at a scale of one.
- **Out of scope** — the warrant fold is recorded as a follow-up in the report's status note.

After the changes: `calling::` 1,301 passed; no clippy warning in the touched files; calls
byte-identical to the reviewed commit on three GIAB configurations, defaults and fitted.
