# Review: the misread-after-a-run table, fitted per read group (2026-11-03)

**Scope.** Branch `str-context-error`, `main...321b744f`: the prototype term (8d0c745c), the
stutter-benchmark script fix (7d07d447), the fit spec (b3734f26) and the fitted table, runs of one
base only (321b744f). Run through the `code-review` skill at level high, with sub-agents; the
findings were then judged and acted on in the main session. Designs:
`doc/devel/ng/spec/str_context_error.md`, `doc/devel/ng/spec/str_context_error_fit.md`.

**Verdict.** One correctness defect, present since the prototype, fixed: a run before an STR tract
could be lengthened by flank bases the candidate deletes or already spells. Four smaller findings
fixed. Five left alone, each for a stated reason. The fix changes no call on HG002 chromosome 20 at
30×.

## 1. Fixed

| # | Where | What was wrong | Fix |
|---|---|---|---|
| 1 | `calling/likelihood/ssr.rs`, `complete_emission` | The run context was read over the whole reference flanks, even where the scored sequence covers flank bases: an allele that deletes flank bases, or a read and allele respelled over a wider flank extent. A left flank ending in `A`×8 with two bases deleted is a run of six, but was read as eight and charged. | `with_flanks_outside` cuts the flanks by the extent scored, in all three branches. Test `a_deleted_flank_base_does_not_lengthen_the_run_before_the_tract`, which fails without the cut. |
| 2 | `run/gatherer.rs`, `CensusPlan::writer_for` | The misread counter was attached by the caller after `writer_for`, which exists so that every census producer builds the same writer (`parameter_prepass_joint_records.md` §7.12). A second producer built through it would have written empty misread tables. | `writer_for` takes the reference and attaches the counter itself. |
| 3 | `parameter_estimation/misread_after_runs.rs`, `contexts_from` | Near a contig end, the end was found by binary search over fetches, about 17 per contig end; an assembly of many short scaffolds pays that per scaffold per sample. | The out-of-bounds error carries the contig's length; one more fetch cuts the stretch at the end. The fitted HG002 table is byte-identical. |
| 4 | `calling/parameters_file/to_run_parameters.rs` | A further-on row with a non-zero `as_the_run` was accepted, although scoring never reads it there. | Refused, naming the row. Test `a_further_on_misread_row_with_an_excess_as_the_run_is_refused`, which fails without the refusal. |
| 5 | `run_context.rs`, `ssr.rs`, `validate.rs`, `census_file.rs` | The run gate joined the flanks and the allele into a new vector per observation and candidate; a doc comment sat on the wrong function; the census comment gave 2.5 kB a read group for the misread counts. | The gate walks the three slices chained; the comment moved; the size is 1.1 kB. |

## 2. Left alone

| Finding | Why not now |
|---|---|
| `letters_after_runs` allocates three vectors and recomputes the run contexts per read, candidate and stretching. | No cost has been measured: chr20 at 30× calls in 12 s against main's 17 s, and with `--defaults` 1m23s against 2m11s (one run each, uncontrolled). Caching per candidate and direction needs scratch threaded through the emission. Worth doing if a profile ever shows it. |
| `bundle_emission` copies the flanks and clones the 19-entry table per emission. | As above. |
| The reverse-direction emission cache is filled even where no candidate is scored per direction. | As above; a per-row flag would avoid it. |
| The gate asks for a run within one base of 8, while a stretched candidate can lengthen a run by up to ten units. | Only a read slipping a run from 7 bases or fewer to 8 or more, in a candidate without a run of 7, goes uncharged, and the charge after runs of 8 to 10 is 2 to 8 in 1,000. Widening the gate would score most mononucleotide tracts twice. |

## 3. Measured after the fixes

HG002 chromosome 20 at 30× against GIAB, missed / false calls (psps of the runs-only build):

| route | before the review | after |
|---|---|---|
| `--defaults` | 747 / 453 | 747 / 453 |
| fitted | 838 / 356 | 838 / 356 |

The fixes change the genotype qualities of 24 records under `--defaults` and 22 fitted, of 102,580,
and no site, allele or genotype. All 5,081 library tests and every integration test pass; clippy
reports nothing in the changed code.
