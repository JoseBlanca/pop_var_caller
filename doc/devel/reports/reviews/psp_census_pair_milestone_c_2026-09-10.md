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

---

## C2 — `--psp` in; `--census` and the five criteria flags out

**Reviewed against:** the working tree over `07fe0e91`, seven files. One read-only agent over the
same four grouped categories, plus the diff's own claims. **Two blockers**, eight should-fix, and
an accounting of every refusal the deleted opener made.

### Findings

**B1 — a rustdoc link to the function the step deleted**, in `census_freshness.rs`. Broken
intra-doc links are denied in `Cargo.toml`, so `cargo doc` fails on it, and the lib suite cannot
see that: rustdoc lints only fire under rustdoc. *Fixed*, along with a second, non-linking stale
mention in `RunError`'s own documentation. **And a finding beyond the step:** `cargo doc` is
already red on this tree — 40 unresolved links, none from this branch — so it is a gate the
milestone baselines have never included.

**B2 — the new laziness test could not fail.** It compared the census reader's byte counter
against the trailers' size; that counter counts *section* reads, and reading a cohort's censuses
reads no section, so the number was zero under every implementation — including one that took each
trailer whole. The assertion reduced to *the psps have a non-empty trailer*, and the test's own doc
claimed the opposite.
*Fixed* with three assertions over two instruments: the psp's own `trailer_bytes_read` must stay at
zero (no trailer taken whole), the census reader's must stay at zero (no section decoded), and then
asking for a section must move it off zero (the evidence is backed by the file, not resident). The
mutation that reads each trailer whole and decodes it resident fails it.

**M1 — `&mut OpenPspCohort` for a function that seeks nothing.** It read each psp's footer, which
opening already read, and then opened the file again by path. *Fixed:* a read-only sibling accessor,
and the census reader takes `&OpenPspCohort` — which is what lets a cohort be read while something
else holds it, the shape a thousand-sample run wants.

**M2 — the `--psp` listing rule written out twice**, with a third copy due at `regenerate-census`,
and the new doc asserting the equality in prose. *Fixed:* one `psp_inputs` module with five tests,
and each command dresses its refusal in its own words. The sibling's more helpful empty-directory
message is now both commands'.

**M3 — four prose statements the code contradicted:** "six tests went" (five did, three of them
about the pairing); "a few hundred bytes" for a head read of up to a megabyte a sample; "the flags
are still on the command line, and plan step C2 removes them", in C2; and `CohortRefusal`'s three
causes described as two. *All fixed.*

**M4 — two names that said the wrong thing:** `censuses` holding psp paths, `open` holding census
evidence — the second a leftover of the deleted `OpenCensusCohort`. *Fixed at all seven sites.*

**M5 — a refusal test asserting only the outer variant.** *Fixed:* it names the psp whose census is
missing, so C3 has something to tighten rather than replace.

### What the reviewer accounted for, and it is the useful half

Every one of the deleted opener's seven refusals was traced to where it lives now: four moved and
three of those got stronger, and the three census-versus-psp identity refusals are structurally
gone. **One input is newly unrefused and unreachable today** — a trailer holding another sample's
census — and the reviewer traced what would happen if it ever arose: not wrong numbers, but a
sample-name mismatch at call time, because the parameters file is bound to a run by name.

**Laziness was traced hop by hop** rather than assumed: `backed` → `renumbered`, which rewrites the
directory's keys and keeps the backed variant → `keys`/`holds`/`len` answering from the directory →
the first disk byte only in `fill`. That is what made B2 findable.

**And the cost that is not in any test:** the fit now holds every psp open for the run, each
keeping its block index — about 336 kB a sample at whole-genome scale, roughly 340 MB at a
thousand. Spec §5 asks for that shape; it is named in the implementation report because it is the
fit's first per-sample resident cost that grows with the genome.

### The mutations

| mutation | outcome |
|---|---|
| each trailer taken whole and decoded resident | 1 test fails — the rewritten laziness test |
| a census that will not read skipped instead of refused | 2 fail |
| each psp paired with another psp's path | 6 fail |
| the evidence's sample order reversed | 1 fails |
| the read-group table's file column taking the wrong path | 3 fail |

Two of the reviewer's five were predicted to survive. The first was B2 and is now caught; the
second — the file column — was predicted to survive because nothing read it, which is also why the
byte-identity oracle could not see it: the path never reaches the parameters file. One assertion
pins it.

---

## C3 — the refusal, before the reference

**Reviewed against:** the working tree over `a43ced4d`, four files. One read-only agent over the
usual grouped categories, and asked in addition to read the message itself as the person who meets
it when a 60-sample fit stops. **No blockers**, twelve should-fix and minor findings, and twelve
mutations with predictions.

**A process note the reviewer raised, and it is a real one.** The working tree changed under the
review: a mutation of mine was applied to `census_freshness.rs` while it was reading the file, so
its report is against the diff as first captured. It handled it by giving every mutation as an exact
`old → new` text pair rather than a line number. **Running mutations in the same worktree a review
is reading is a mistake to stop making** — the two must not overlap in time, or the review must be
given a copy.

### Findings

**S1 — the message inflected one of its three number-bearing phrases.** *"1 of this cohort's 1
samples cannot be fitted as they stand"* at a cohort of one, which is the low end of the range this
caller is built for and the case a person is likeliest to meet by hand. *Fixed everywhere, with a
test that reads the singular case whole.*

**S2 — the only arithmetic in the message was asserted by nothing.** The first line is what a
reader acts on before any other, and the mutation that changes *of this cohort's N samples* to
something else broke no test. *Fixed:* the line is asserted whole, and that mutation now fails
three tests.

**S3 — "Rebuild them" pointed at the psps that cannot be rebuilt.** In a cohort with both faults
the unreadable rows printed directly above it. *Fixed:* the command names the set it is for, a
separate line says the unreadable psps will not be mended by rebuilding, and the command is an
`Option` that exists only when something is stale — so the type cannot say *rebuild* without one.
The mixed case had no test and has one.

**S4 — the unreadable row printed its path twice and dropped the cause.** It stored the outermost
message; the fault hangs off `#[source]`. *Fixed* with `format_error_chain`, and the fixture's error
now carries the row's own path so the duplication would show.

**S5 — "before anything else is read" was false**: the cohort opener has already read every header,
footer and block index. *Fixed* to the claim that is true and load-bearing — before the reference.

**S6 — two doc comments still described a census file beside its psp**, one of them two lines above
the new code. *Fixed.*

**S7 — the constant naming a command that does not exist yet is the right carrier**, and printing
`regenerate-census` rather than today's `generate-census` is right. Two gaps the reviewer named:
nothing checks the printed name is a subcommand clap accepts, and plan step D1's task list does not
mention the constant. *A test now pins what can be pinned today* — that the name is not the command
which writes a census file — and D1 turns it into a parse. The plan gap is recorded for the
checkpoint rather than edited in.

**Minor, all fixed:** three-field tuples where row structs belong; an error variant whose name was
false in one of its two cases; a command line built eagerly on every run; `push_str(&format!(…))`
per argument; an unquoted path in a line whose whole value is that it is pasted; the report tests
splitting the cohort tests in two; and a doc sentence that did not parse.

### Recorded, not fixed

- **A psp that will not open never reaches this report**: the cohort opener refuses at the first
  one, so spec §4.1's *every sample is examined* holds for censuses and not for that fault.
- **The version-word surgery is now written in two test modules.** The reviewer asked for one
  fixture beside `a_census_this_build_wrote`; it is left where it is, and named here so the next
  step can move it rather than add a third.

### The mutations

| mutation | outcome |
|---|---|
| the run stops at the first stale psp | 4 tests fail |
| the reference read before the judgement | 1 fails |
| the header's arithmetic changed | 3 fail |
| the sample-name column dropped | 6 fail |
| the cause dropped from each line | 5 fail |
| an unreadable psp reported as one to regenerate | 2 fail |
| the stale lines left ungrouped | 2 fail |
| the cause chain not walked | 2 fail |
| `--psp` dropped from the command | 2 fail |
| `--catalog` dropped from the command | 2 fail |
| the old command name printed | 1 fails |

Two of the reviewer's predicted survivors — the header arithmetic and the sample column — are now
caught. The third, the constant's own value, is caught only against the old command's name, which
is the most that can be asserted before the new command exists.
