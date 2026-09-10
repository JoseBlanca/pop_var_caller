# Code Review: the census lives inside the psp — Milestone C

**Date:** 2026-09-10
**Branch:** `census-vs-psp-perf`
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone C
**Implementation report:** [ng_psp_census_pair_milestone_c_2026-09-10.md](../implementations/ng_psp_census_pair_milestone_c_2026-09-10.md)

---

## How this review is run

**Milestone A's arrangement, unchanged through B.** The skill asks for one agent per category, each
in its own worktree; the build tree is 17 GB and the disk has 129 GB free, so that does not fit.
Instead one read-only agent covers grouped categories, forbidden to edit any file or run `cargo`,
and asked to **name** the mutations it wants rather than run them. The orchestrator runs them, one
at a time, restoring from a backup and proving the restore with `diff` before the next.

---

## C1 — the criteria from the header

**Reviewed against:** the working tree over `d3cfff1a`, three files. One read-only agent over four
grouped categories — reliability and errors, idiom and smells, naming and module structure, refactor
safety — plus the diff's own quantitative claims. **No blockers.** It confirmed the refusal order is
preserved for all eleven existing callers of `segments_over`, and re-derived the byte-identity
oracle from the artefacts on disk rather than taking it on trust: the before and after parameters
files share a SHA-256, so it is identity and not equal length.

### Findings

**M1 — the shared tail tells a header-driven caller to move a flag it ignores.** The refusal for
*this reader asks for tracts below what the catalog holds* names the flag that asked, because on a
walk that is what the person typed. On the new entry point nobody typed it. Measured, with a catalog
built at twenty copies: *"`--min-copies` asks for repeats the catalog … does not hold; raise it"* —
a sentence whose only actionable half is a knob this command does not read.
*Fixed:* `segments_cut_with` takes a `CriteriaSource`, and criteria that came from psp headers get
the general catalog refusal, which names the file the reader can change. A test asserts the message
names the catalog and does not name the flag; the mutation that reinstates the old rendering fails
it.

**M2 — four statements in new prose that the code contradicts.** *"The ground and what it is cut
with are the psps' own"* (they are the first psp's, and the check that makes one psp enough is
`CohortCensusEvidence::new`'s, elsewhere); *"the conversion below checks the same thing again"* (the
conversion checks nothing about the catalog); a doc listing *"finding the catalog"* four lines above
saying the path arrives resolved; and a test doc claiming the pre-C1 fit *"cut the ground
differently"* where it would in fact have been refused as built under another selection.
*All four fixed*, and the first now names which check is load-bearing, because the reader who moves
this code at C2 is the one who would otherwise drop it as redundant.

**M3 — a mechanism claim more confident than the code supports.** The comment leaned on the fit
refusing a selection built under other criteria. `fit_a_cohort` compares a digest of the kept
**generic** positions only, so it catches a difference only where one of those moves; at tomato's
1-in-400 keep rate a criterion that retypes one short tract can pass it. *Fixed*, and it is now
stated as the argument for the step rather than as a net.

**M4 — the catalog-path default written out four times.** `--catalog`'s promise —
*defaults to `<reference>.repeats.parquet`* — is part of the user-facing contract, and a fifth copy
arrives with `regenerate-census`. *Fixed:* one `run_ground::catalog_path_for`, which
`GroundRequest::catalog_path` now calls.

**M5 — two assertions satisfied by any successful fit.** The `--catalog` test asserted a sample
count and a non-empty term list, both true of every fit in the file. *Fixed:* it compares the file
fitted from the moved catalog against the same cohort fitted with its catalog where a run looks by
default — one catalog read from two paths is one answer.

**M6 — nothing pinned that the ground comes from the cohort**, and no fixture in the module could:
every walk covers whole contigs, where the psps' ground and the whole reference are the same
positions. *Fixed:* the fixture takes a BED, and a cohort walked over the first 400 bases of a
600-base contig is fitted with no ground given. The mutation that takes the ground from the
reference fails that test and nothing else.

**M7 — nothing pinned the refusal order.** A missing catalog is reported before a backwards period
range, and the test named for that order calls the conversion directly and never meets a catalog.
*Fixed* with a test in `call_from_alignments` that sets both.

**Minor, all fixed:** a bare participle for an opened catalog (`let open = …`), the `# Errors` list
in an order the code does not raise them in, and a module doc that did not mention the psp header
is read for its criteria.

### Recorded, not fixed

- **The catalog the fit reads is never compared with the catalog the psps name.** Spec §6 says
  `--catalog` is checked against the header's digest "as today"; no such check exists today. It
  belongs in C2's cohort opener, where `SegmentationInputs::first_difference` already names the
  field.
- **The `Segmentation` this command builds is never read** beyond the criteria handed into it. It
  is worth keeping for the refusals it raises, but it materialises every typed region of the
  analysed ground to do it — the one cost here that grows with the genome. C2 rewrites the block.
- `generate-census` still takes its criteria from its own flags, the shape `estimate-parameters`
  has just left. Milestone D replaces it.

### The mutations

| mutation | outcome |
|---|---|
| the criteria taken from the five flags again | 2 tests fail |
| the copy floor alone taken from `--min-copies` | 2 fail |
| `--catalog` ignored, the sibling path always read | 2 fail |
| the ground taken from the reference, not the cohort | 1 fails — the new ground test alone |
| the coarse-catalog refusal naming `--min-copies` again | 1 fails |
| the catalog checked after the flags are converted | 1 fails — the new refusal-order test |
| the ground cut with `StrRepeatCriteria::default()` | 0 here, 4 elsewhere in the lib suite |
| the criteria not recorded in `Segmentation::build` | 0 here, 5 elsewhere |

The last two are the reviewer's prediction confirmed: on this command the cut is unused and only
the record is read, so no test of `estimate-parameters` can see either, and `call-from-alignments`
and `generate-psps` are what catch them.

**One mutation reported a false pass and had to be run again.** The `--catalog` edit was written
against a line the review's own fix had rewritten, so it matched nothing and the tree was never
mutated — indistinguishable, in the log, from a mutation no test catches. The harness asserts its
match count; sending its output to `/dev/null` in the same command as the test run is what hid the
assertion.
