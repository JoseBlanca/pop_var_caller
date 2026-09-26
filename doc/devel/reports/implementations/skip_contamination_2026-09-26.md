# `estimate-parameters --skip-contamination`

*Implementation, review and fixes, 2026-09-26. Branch `skip-contamination`. Asked for by the owner
after the 2,169-sample run on kimura showed contamination to be the step that needs the most
memory: in many cohorts the estimate is not needed.*

## What it does

`--skip-contamination` (off by default) makes the SNP/indel fit return, for every (sample, read
group), `NotIdentified { reason: NotAsked }` without choosing markers. Per
`doc/devel/ng/spec/parameters_file.md` §3.4 a run where no read group has a fraction writes **no
contamination table**, which a calling run reads as uncontaminated — the plain read likelihood, not
a fraction fitted at zero. `RunParameters::assemble` already does this for any `NotIdentified`
reason, and a test pins it (`a_run_where_nothing_was_identified_is_uncontaminated`).

- `JointFitConfig::estimate_contamination` (default `true`);
- `NotIdentifiedReason::NotAsked`, with its message;
- `contamination::not_identified_anywhere`, shared by the refused panel and the skipped run;
- the command's help, and `ContaminationUsed::NoneFitted`'s doc, which now says a skipped run
  reports the same: **the file does not record that contamination was skipped rather than
  unidentifiable; only the progress log does.** Recording it would be a change to the
  parameters-file spec, left to the owner.

With the default nothing changes; the cross-platform checksum test passes unchanged — though its
fixture, like every command-level fixture, refuses contamination (too few markers), so that test
cannot see this switch.

## Tests

- `contamination::tests::a_fit_told_not_to_estimate_contamination_returns_no_fraction_anywhere` —
  a 30-sample panel with one sample contaminated: skipped, all 30 read groups are `NotAsked`;
  estimated, fractions come back.
- The command's parse test: off by default, on with the flag, and — because no fixture file can
  show it — the choice reaching the fit's config (`ordinary_position_config`). Flipping the
  negation fails this test (checked by hand, file restored).

## Review

One reviewer in its own worktree. Blocker: nothing tested that the flag reaches the fit, and no
command-level fixture could — fixed as above. Minors applied: the shared helper; "another
individual" rather than "another plant"; help text pointing at the printed `contamination:` line
instead of an unmeasured size claim; the `NoneFitted` doc. Nit not applied: moving the new fit
test into `fit.rs`.

## Validation

In the dev container: fmt and clippy `-D warnings` clean; targeted tests 62 passed; the full
suite before the review's fixes 4,880 passed, 3 failed — the pre-existing failures in
`examples/ng_generic_loci_dump.rs` and `examples/ng_ssr_loci_dump.rs`.
