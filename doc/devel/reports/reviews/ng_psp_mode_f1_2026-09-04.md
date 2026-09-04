# psp mode F1 — review, and what was done about each finding

**Date:** 2026-09-04
**Reviewed:** the working-tree diff of
[run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) step F1
**Implementation report:** [ng_psp_mode_f1_2026-09-04.md](../implementations/ng_psp_mode_f1_2026-09-04.md)
**Branch:** `ng-psp-mode`

Two reviewers, each in its own worktree detached at `91d10740` with the step's diff applied, one
on correctness and the spec's clauses, one on the report and on duplication. Between them they
ran fifteen deliberate defects against the scoped tests.

## The answer

**Three defects that a test should have caught survived one, and all three are closed.** The
shipped code was right in every case; what was missing was a test that could tell. The one
finding that changed behaviour rather than tests is the depth wording, and both reviewers reached
it independently.

## The findings, and what was done

### Applied

**A number was called depth and is not depth.** The report printed
`{n} loci read, {d:.1} reads a locus` and two doc comments called it the sample's mean depth. The
quantity is the record head's `reads-compared-with-reference`, which excludes reads a filter
turned away, reads the per-position cap discarded, reads that covered the locus and produced no
observation, and reads whose witness stopped inside it. At a repeat tract 40 reads cover but 22
anchor, it says 22. **The line now reads `reads a locus compared with the reference`**, and the
section's doc says what the number leaves out and why it is still worth printing.

**Three of the command's own refusals were untestable through the command.** Deleting
`refuse_an_output_that_cannot_be_written`, `refuse_an_output_whose_parameters_file_is_this_run_s_input`
and the `--ploidy` flag's journey into the parameters check left all nine of the new command's
tests green — the same gap direct mode's own review measured and wrote down, reintroduced.
**Six tests now drive `run_call_from_psps` itself**, and with those three calls deleted three of
them fail. One of the six had to be tightened first: a missing output directory is *also* named
by the writer, after the whole cohort has been called, so the test now asserts the command's own
sentence.

**Nothing named which sample a set of tallies belonged to.** Reversing the pairing between spent
sources and open files left all 25 `psp_caller` tests green, because both fixture samples read
exactly two loci and three reads. **The fixture now stores two loci for one sample and one for
the other**, the assertion names each sample beside its count, and the pairing carries the same
`assert_eq!` on names the source-building loop forty lines above already carried.

**The union-of-keys rule in the read-filter comparison was unpinned.** Collecting the settings to
compare from the first file only, instead of from every file, left all 24 report tests green: the
one test that could distinguish them put the sample *with* the key first. **The test now asserts
both orders**, and the second half fails under that defect.

**Counting at the decode rather than at the hand-over was untested.** Moving both counters above
the two refusals — the over-count the code's own comment says it avoids — left all 24 report
tests green. **Both refusal tests in `psp_source` now assert the counters**, and the mutation
fails both.

**A parameters file drops a read-filter key a later ng recorded.** The reader kept only keys its
own build's list names, so two files written months apart under different policies compared equal
on the difference that matters. **It now matches on the key prefix**, with a new
`READ_FILTER_PROVENANCE_PREFIX` beside the list and a test holding every listed key against it.

**Two messages changed in the lift.** `ParametersWouldBeOverwritten` was reworded and
`ParametersNotWritten` lost the clause direct mode's own field doc calls load-bearing —
*keep the VCF and re-run to recover its parameters*, without which a `set -e` pipeline throws away
a complete, correctly-headed file. **Both restored to direct mode's wording.**

**`header_for` and the report printer were copied between the two commands**, byte-identical
bodies. One of them encodes a decision that survived a mutation in direct mode until a
command-level test existed. **Both lifted into `calling_run`.**

**The fixture lift reached one of the two commands its own doc claimed.** `call-from-alignments`
still had a private cohort, and it gave both samples no reads — so its only end-to-end command
test drove a cohort with nothing in it. **It uses the shared one now**, and the shared module's
doc says what the two copies differed in.

**Two wordings that said the wrong thing.** The psp report's ground line said *bases asked for*,
where nobody asked and there is no `--regions`; it now says *as every file's header records
them*. The command's `# Errors` list named an order the code does not follow, because the round
width has to wait for the psp directories to be expanded; it now names the real one.

### Recorded, not changed

**A psp against a reference with a different contig count gets the catalog's refusal rather than
its own** — still a refusal before a block is decoded, and reordering would mean building the
segmentation after the cohort check, which needs the segmentation. In the implementation report's
*what is knowingly left*.

**The per-sample section has no cap** at the thousand-sample end. Direct mode's has the same
shape; capping one alone would make the two disagree.

**The disagreement line cannot be provoked from the shipped CLI**, because `generate-psps`
exposes no read-filter flag. It fires across ng builds whose compiled defaults differ, which is a
cohort walked over months — the cohort psp mode exists for. Said so in the section's doc.

**The walked branch's restructuring changes behaviour only for a cohort with no samples**, which
`--alignment`'s `required = true` makes unreachable and no test builds.

**Rendered-string comparison of parameter values** treats `Integer(20)` and `String("20")` as
equal. Reachable only through an overflow arm needing a read length above 9.2×10¹⁸.

## Mutations, and where they stand now

| mutation | before | after |
|---|---|---|
| `loci_read += 1` → `+= 2` | caught | caught |
| `mean_reads_a_locus` returns the sum | caught | caught |
| each sample's counts under another's name | **survived** | caught |
| the disagreement line never fires | caught | caught |
| the command's three refusals deleted | **survived** | caught |
| compared keys taken from the first file only | **survived** | caught |
| loci counted at decode rather than hand-over | **survived** | caught |
| "loci read" prints the compared-read sum | caught | caught |
| the two counts in the samples line swapped | caught | caught |
| the "no locus" list names the wrong samples | caught | caught |
| the psp ground line prints the walked partition | caught | caught |
| the "not recorded" arm dropped | caught | caught |
| the sample classification branches swapped | caught | caught |
| the disagreement line lists only the first value | caught | caught |

## After the fixes

`cargo test --lib` in the container: **6,214 passed, 0 failed, 15 ignored**. Both oracles still
green — direct mode byte-identical across the lifts, and psp mode's VCF byte-identical to direct
mode's apart from `##commandline`, 599 records on the six tomato accessions.
