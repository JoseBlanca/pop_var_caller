# ng — the census lives inside the psp

*Design spec, 2026-09-09. **No code yet — this settles the design.** Companion plan:
[`impl_plan/psp_census_pair.md`](../impl_plan/psp_census_pair.md). **It reverses
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1's ruling of
2026-08-13** — one census file per sample, beside the psp — which shipped as code; §6.1 carries the
supersession note. An earlier revision of this document, the same day, kept the sidecar and made its
location mandatory; the owner's ruling below replaced it before any code was written.*

*Vocabulary, once. A **psp** is one sample's stored evidence: what its reads showed at every position
of the ground that was walked. A **census** is the much smaller object a parameters fit reads: what
that sample showed at a fixed set of positions and repeat tracts chosen for the whole run, so the fit
can ask the same question of every sample. The **selection** is that fixed set. The **trailer** is
the psp format's closing payload — a run of bytes the writer supplies when it seals the file, that
the container stores without interpreting, and that can be replaced without touching the records.
A census is **fresh** when the psp's trailer holds one in the format this build reads, recording the
selection this run rebuilds; otherwise the psp **needs its census regenerated**.*

---

## 1. What it is

**One file per sample.** The census is written into the psp's trailer when the walk seals the file,
and read from there by the fit. There is no census file. The commands a user runs, in order:

| step | command | in | out |
|---|---|---|---|
| 1 | `generate-psps` | alignment files | `<sample>.psp`, census inside |
| 2 | `estimate-parameters` | psps | the parameters file a calling run scores with |
| 3 | `call-from-psps` | psps and that file | the VCF |

And one repair command, run only when step 2 says to: **`regenerate-census`** — psps in, each one's
trailer replaced with a fresh census.

**What changes from today.** Step 1 writes `<sample>.census` beside the psp
([`generate_psps.rs:676-735`](../../../../src/pop_var_caller_exp/generate_psps.rs)) and the psp's
trailer is empty ([`psp_writer_line.rs:244`](../../../../src/ng/run/psp_writer_line.rs)); the census
moves into that trailer and the second file goes. Step 2 takes `--census` paths and six flags that
must equal what the psps were walked under
([`estimate_parameters.rs:69-141`](../../../../src/pop_var_caller_exp/estimate_parameters.rs)); it
will take `--psp` and read those settings out of the psp headers. `generate-census` and its
`--output-dir` ([`generate_census.rs:97`](../../../../src/pop_var_caller_exp/generate_census.rs))
become `regenerate-census`, which writes into the psp.

**Goals.**

- **The user names psps and nothing else, at every step.** Nothing else exists to name.
- **A psp whose census is missing or stale stops step 2 with a report that names every such
  sample, the reason, and the command to run.** Step 2 does not rebuild on its own.
- **Nothing a census depends on is typed twice.** The repeat criteria, the catalog and the analysed
  ground come from the psp header; the flags that duplicated them go.
- **A cohort whose psps disagree on the catalog or the criteria is refused by every command that
  opens one**, `call-from-psps` included.

**Non-goals.** A census file beside the psp, in any form — reversed, §3. Rebuilding inside step 2 —
rejected, §4. A flag to skip the census at step 1 — rejected, §3.3. Storing a larger census than the
fit uses — deferred, §7. Changing the psp format: the trailer is already there, and what goes in it
is the writer's business ([`psp_file_format.md`](psp_file_format.md) §3.4). What the census
*holds*, how it is encoded, and the fit itself — unchanged.

**It does not:** open an alignment file after step 1; read a psp's records at step 2; write any
file but the psp; accept a census path from anyone.

---

## 2. Why the census is kept at all — the cost of not having one

**The census is a cache, and rebuilding it means reading the whole psp.** It can be rebuilt with no
alignment file — that is what today's `generate-census` does — but the rebuild decodes every record
([`census_from_psp.rs:241-252`](../../../../src/ng/run/census_from_psp.rs)). It cannot skip records
by their head, because the read-error calibration totals it accumulates are summed over every
generic locus, and the per-read quality sums they need are in the record body
([`census.rs:2383-2396`](../../../../src/ng/parameter_estimation/joint/census.rs)).

**Measured, 2026-09-09**, with `examples/ng_census_read_vs_psp.rs` (branch `census-vs-psp-perf`; the
plan lands it), one thread, files in the page cache, the census at the density a whole-genome run
has — two million positions over the genome, 1 in 400 on tomato and 1 in 1,500 on human:

| sample | ground | reads a position | psp | read the census | rebuild it from the psp |
|---|---|---|---|---|---|
| tomato SRS3394606 | 32 Mb | 22.6 | 10.4 MB | 0.0002 s | 0.25 s |
| tomato DRR000741 | 10 Mb | 98.5 | 158 MB | 0.0002 s | 3.0 s |
| HG002, the 50× set | 4.13 Mb | 33.5 | 43.7 MB | 0.0002 s | 1.41 s |
| HG002, the 300× set | 5.08 Mb | 300.9 | 126 MB | 0.0001 s | 3.25 s |

The rebuild runs at about 45 MB of psp a second whichever corner it is in. **At 50× a human genome
is about 32 GB of psp a sample, and its rebuild about a quarter of an hour a sample**, on one
thread. Reading the census is a few milliseconds. Those numbers are why the census exists, why step
2 must not rebuild it unasked, and why the rebuild is a command that can be spread over machines.

---

## 3. Where the census lives — decided: in the psp's trailer

**Decision (owner, 2026-09-09): the psp is one file; the census is its trailer.** Three options were
live. **(a)** The trailer. **(b)** A file beside the psp, wherever the user says — today. **(c)** A
file beside the psp, same stem, nowhere else — this document's first revision.

**Why (a).** Simpler for the user: one file to make, copy, archive, name. And the reasons for a
separate file turned out to be about cases that do not arise. A census must be rebuilt *without*
re-walking the psp only when what the census holds changes — its format, or the selection's
compiled-in seed, budget and cap — and the owner's judgement is that neither is common once psps
exist at scale. Every other change that invalidates a census invalidates the psp too: the repeat
criteria partition the ground into tract loci and generic loci before a record is written
([`region_typing/mod.rs:279-299`](../../../../src/ng/region_typing/mod.rs)), so different criteria
mean different records, and the reference, the catalog and the analysed regions likewise.

**What (a) costs, accepted.** Regenerating a census rewrites the psp's tail: `replace_trailer`
truncates at the trailer offset and writes the new trailer and footer forward, leaving the header,
the blocks and the index untouched ([`psp/trailer.rs:60-103`](../../../../src/ng/psp/trailer.rs)).
Between the truncation and the new footer's last byte the file has no footer and every reader
refuses it — a partly written psp *should* be unusable, and the owner accepts that a process killed
in that window costs the sample a re-walk. And a psp on read-only storage cannot have its census
regenerated in place; it has to be copied somewhere writable first. The rebuild's read pass (§2)
costs the same in every option, so nothing here is a time argument.

**What is gained besides the file count.** No pairing. A census cannot be lost, copied apart from
its psp, or built from a different psp than the one beside it, so the identity a census carries
today — the digest of its psp's header and the record count
([`PileupIdentity`, `census_file.rs:83-92`](../../../../src/ng/parameter_estimation/joint/census_file.rs))
— has nothing left to check and is written absent. What remains to check is in §4.2.

### 3.1 The walk writes the census into the trailer as it seals the file

`PspWriter::finish` takes the trailer's bytes ([`writer.rs:672`](../../../../src/ng/psp/writer.rs)),
so the census is encoded in memory once the last locus has gone past and handed to `finish`. Today
the order is the other way round — the writer line is finished first, then the census
([`gatherer.rs:506-507`](../../../../src/ng/run/gatherer.rs)) — because the census needed the
psp's header digest, which it no longer does. Holding the census in memory until the end is what
happens today ([`CensusWriter`](../../../../src/ng/parameter_estimation/joint/census.rs)
accumulates and `finish` encodes); a whole-genome census is a few megabytes of positions plus tens of
megabytes of tracts — 1.03 bytes a position and about 35 a tract, measured on the tomato accession
at every position kept.

**Building the census during the walk is cheaper than not building it and rebuilding afterwards**
— 1.28 s against 1.40 s on six tomato accessions over 200 kb
([`ng_fit_stage_b_2026-09-05.md`](../../reports/implementations/ng_fit_stage_b_2026-09-05.md)) —
and it is always built. There is no walk without a census.

### 3.2 One `.partial` file, one rename — simpler than today

`generate-psps` writes the psp to a `.partial` name and renames it whole; today it does the same for
the census and reasons about which of the two to rename first
([`generate_psps.rs:713-735`](../../../../src/pop_var_caller_exp/generate_psps.rs)). With one file
there is one rename and nothing to order.

### 3.3 No `--skip-census` — decided

It would save a tenth of a percent of the psp's size and a small share of the walk, and manufacture
the one state the design exists to remove: a psp that looks complete and costs a quarter of an hour a
sample to make usable, discovered at step 2. `generate-psps` has no such flag today and keeps none.

### 3.4 Appending to a psp discards its census — recorded, not handled

`PspWriter::append` truncates at the index offset and discards the index, the trailer and the
footer ([`writer.rs:420-441`](../../../../src/ng/psp/writer.rs)). No shipped command calls it, and
the owner's ruling is that there is no user case for appending: regenerate the sample instead. If it
is ever used, the appended psp has an empty trailer, which §4.2 reports as *no census*.

---

## 4. When a census is missing or stale — decided: step 2 refuses and reports, it does not rebuild

Two options were live: rebuild inside step 2, silently — what §6.1 of the records spec wrote, and
what the owner first proposed; or refuse, naming every sample that needs regenerating and the
command that does it.

**Decision (owner, 2026-09-09): refuse.** "If we can't run the estimate, we stop, inform the user,
and the user has to decide to regenerate the census." The rebuild is a quarter of an hour a sample
and invisible: an `estimate-parameters` over thirty human samples that quietly rebuilt would run for
eight hours with nothing to show that anything unusual was happening. It is also per-sample and
distributable, and a rebuild inside step 2 would run serially on one machine. Refusing makes the
cost a decision. **The refusal does not estimate how long regeneration will take** — owner: that
depends too much on the hardware.

### 4.1 Every sample is examined before anything is refused

Today the cohort opener returns at the first census it cannot open or check
([`census_cohort.rs:184-215`](../../../../src/ng/run/census_cohort.rs)): a cohort with three stale
censuses reports one. **The whole cohort is examined and the report names all of them**, grouped by
cause, so one run tells the user the whole regeneration job.

### 4.2 What makes a census stale, and where each is caught

| cause | what is compared | cost |
|---|---|---|
| no census | the footer's trailer length is zero ([`reader.rs:228-247`](../../../../src/ng/psp/reader.rs)) | the footer, already read at open |
| a census of an older format | the version word at the trailer's front against `VERSION` ([`census_file.rs:75`](../../../../src/ng/parameter_estimation/joint/census_file.rs)) | one short read at the trailer offset |
| recorded under other settings | the twelve settings the census records against the ones this run records under, naming the first that differs ([`what_the_run_says_about_every_census_in_a_cohort`](../../../../src/ng/run/census_freshness.rs) against [`CensusPlan::recording_terms`](../../../../src/ng/run/gatherer.rs)): seven say which positions were chosen, five say in what units the evidence was written down. The digest of the kept positions is the backstop behind it ([`fit_a_cohort`](../../../../src/ng/run/census_fit.rs)) | the reference read and the selection rebuilt — 4 to 19 s in §2's measurements |

The first two are judged before the reference is read, so a refusal for a missing census is
immediate. The third needs the reference read and the selection rebuilt, so it is taken after
those and before anything is fitted: one verdict a sample, naming the setting that differs and the
command that regenerates it. **Its instruction needs no hedge because §6's two file checks come
first** — a run pointed at another reference or catalog is refused as having named the wrong file,
so a census that still differs from the run does not match its own psp, or was written by a build
that chooses positions differently, and regenerating is the fix for both.

**An older format is named as such, not as damage**: today a version word this build
does not know is reported as *malformed* ([`decode_census`, `census_file.rs:456-461`](../../../../src/ng/parameter_estimation/joint/census_file.rs)),
which sends the user looking for corruption when the truth is that this build changed what a census
holds. The version is read first and reported as *built by an older version; regenerate*.

### 4.3 What the refusal says

Per sample: the sample name, the psp path, the cause in the words of §4.2. Then one line with the
command, built from what was given — `regenerate-census --psp <the same arguments>` — so it can be
copied. Fresh samples are counted, not listed.

---

## 5. Reading the census out of the psp without reading the psp

**The fit's reader opens a census by path today and reads only its header and directory, seeking
to each section when it is asked for** ([`open_census`, `census_file.rs:494`](../../../../src/ng/parameter_estimation/joint/census_file.rs);
[`SampleCensusEvidence::backed`, `census.rs:1280`](../../../../src/ng/parameter_estimation/joint/census.rs),
whose section reads open the file and seek to the section's offset, l.1326-1330). That is what lets
a cohort too large to hold be opened lazily, and it must survive. **The reader takes the psp's path
and the trailer's extent — offset and length from the footer — and every seek adds the offset.** The
psp's own `trailer()` reads the payload whole ([`reader.rs:228-247`](../../../../src/ng/psp/reader.rs))
and is not what the fit uses; a thousand-sample cohort read whole would be thirty gigabytes of
census in memory.

**The psp's `read_header` does not read the footer.** Step 2 needs the header for the settings (§6)
and the footer for the trailer's extent; `PspReader::open` reads both
([`reader.rs:85`](../../../../src/ng/psp/reader.rs)) and holds one descriptor, which is the shape
`OpenPspCohort::open` already uses for a cohort ([`psp_caller.rs:111`](../../../../src/ng/run/psp_caller.rs)).

---

## 6. Where the settings come from — decided: the psp header, never a flag

**Every setting the selection depends on is in the psp header already.** `SegmentationInputs`
carries the catalog's own header, the repeat criteria the walk routed under, and the analysed
regions ([`segmentation_inputs.rs:25-40`](../../../../src/ng/segmentation_inputs.rs)). Step 2 and
`regenerate-census` read them from there.

**And they must check that the cohort agrees on them, which no cohort *opener* does today.** The
openers compare the analysed regions only ([`psp_caller.rs:677-690`](../../../../src/ng/run/psp_caller.rs);
[`census_cohort.rs:222-236`](../../../../src/ng/run/census_cohort.rs)).
`SegmentationInputs::first_difference` names the first of the three fields that differs, in the
order a person should fix them — catalog, criteria, regions
([`segmentation_inputs.rs:42-55`](../../../../src/ng/segmentation_inputs.rs)) — **and outside its own
tests it is called from one place: `PspVariantCaller::open`, which compares each psp against the
segmentation the run built** ([`psp_caller.rs:336-346`](../../../../src/ng/run/psp_caller.rs)). So a
calling run already refuses a cohort typed two ways; what it names is one sample and the run, rather
than the two samples that disagree, and the commands that build no run segmentation — the fit and
the repair — have no such refusal at all. A step 2 that took the criteria from the first psp's
header and never asked the rest would build a selection the other samples' censuses cannot match,
and learn it twenty seconds later from the fit's digest check, reported as *another selection*.

**Decision (owner, 2026-09-09): a cohort whose psps disagree on the catalog or the criteria is a
hard failure for every command that opens one** — `estimate-parameters`, `regenerate-census`, and
`call-from-psps` alike, since two psps typed under different criteria cannot be called together
either. The check goes where the analysed-regions check already is, `OpenPspCohort::open`, which
`call-from-psps` and today's `generate-census` already call
([`call_from_psps.rs:494`](../../../../src/pop_var_caller_exp/call_from_psps.rs),
[`generate_census.rs:413`](../../../../src/pop_var_caller_exp/generate_census.rs)); step 2 opens its
cohort the same way. Written once, refusing with the sample and the field named, before anything
else is read.

**Flags removed from `estimate-parameters`:** `--census`, `--min-copies`, `--min-period`,
`--max-period`, `--max-str-len`, `--min-purity`. **Kept:** `--reference` (the FASTA's path is not in
the header, and the selection needs the reference to know where it is sequence at all),
`--catalog` (defaulting to `<reference>.repeats.parquet`, and checked against the header's catalog
digest as today), `--output`, `--force`, `--ploidy`, `--inbreeding`. **Added:** `--psp`, once per
file or a directory, as `call-from-psps` takes it. `regenerate-census` takes `--psp`, `--reference`
and `--catalog`, and nothing else.

**Trap for the coder.** `estimate-parameters` builds its ground from the five flags as a
`run_ground::RepeatRouting` ([`estimate_parameters.rs:312-322`](../../../../src/pop_var_caller_exp/estimate_parameters.rs)),
and `segments_over` converts that to a `StrRepeatCriteria` itself, through `routing_criteria`
([`run_ground.rs:242-281`](../../../../src/pop_var_caller_exp/run_ground.rs)). The header holds the
finished `StrRepeatCriteria`, and there is no conversion back. Split `segments_over` after its call
to `routing_criteria`, so the flag-driven commands and the header-driven ones share everything past
that line.

---

## 7. The selection's two counts, and why storing more would be free — deferred

The budget of positions and the cap of tracts per stratum are compiled-in constants: about two
million and five thousand, with a fixed seed ([`CensusSelection::SHIPPED`, `gatherer.rs:173`](../../../../src/ng/run/gatherer.rs)).
No command takes them as flags. Nothing here changes them; a build that changes them makes every
census stale under §4.2's third row, which is the intended consequence.

**Recorded because it decides what a later knob would cost.** The selection nests. A position is
kept when `hash(contig, position, seed) < threshold`, the threshold rises with the budget, and the
hash depends on the position and the seed alone ([`loci.rs:316-376`](../../../../src/ng/parameter_estimation/joint/loci.rs)):
the set kept at two million is a strict superset of the set kept at half a million, and a fit
wanting the smaller filters the stored one and gets exactly the positions a fresh selection would
choose. The per-stratum cap is the lowest-hash prefix of the same kind
([`StratumSample`, `strata.rs:58-70`](../../../../src/ng/repeat_catalog/strata.rs)). A census stored
at ten times the budget serves any smaller one without a rebuild, for about one part in a hundred of
the psp. **Deferred** to [`parameter_prepass_runs.md`](../impl_plan/parameter_prepass_runs.md)
Checkpoint C's knobs question, where the fit-side knob it would need is owned.

---

## 8. `regenerate-census` — the repair command

**Psps in; each one's trailer replaced with a fresh census; nothing else written.** It replaces
today's `generate-census`.

- **`--psp`** once per file or a directory; `--reference`; `--catalog` defaulting as in §6. The
  criteria and the ground come from the headers; the cohort's agreement on them is checked first.
- **No `--output-dir`, no `--force`.** The census goes into the psp; replacing is the command's
  job. A psp that refuses the write — read-only storage — fails that sample with the path.
- **Fresh psps are skipped, and the run says so.** Freshness is §4.2 whole — all three rows,
  since this command rebuilds the selection anyway — so a run stopped part-way and started again
  does only the samples still owed, and a build whose selection constants changed regenerates every
  one without being told to.
- **Per sample: one read pass, then `replace_trailer`.** The pass is §2's quarter of an hour at
  50× human; the write is milliseconds. One psp open at a time, as today
  ([`generate_census.rs:29-33`](../../../../src/pop_var_caller_exp/generate_census.rs)).
  Progress to stderr as each sample finishes, in the words of the final report, as `generate-psps`
  does.
- **Spread it the way step 1 is spread**: one invocation per sample or per directory, on as many
  machines as there are. This is why it is a command and not a mode of step 2.

---

## 9. What will bite the coder

- **The writer line finishes before the census does today** ([`gatherer.rs:506-507`](../../../../src/ng/run/gatherer.rs)),
  and `PspWriterLine::finish` passes an empty trailer ([`psp_writer_line.rs:244`](../../../../src/ng/run/psp_writer_line.rs)).
  The census has to be finished and encoded first and handed through the line to
  `PspWriter::finish`. The line runs the writer on its own thread; the bytes cross that seam once.
- **No shipped command checks that a cohort's psps agree on the catalog or the criteria** — §6.
- **The cohort opener stops at the first failure** ([`census_cohort.rs:184-215`](../../../../src/ng/run/census_cohort.rs)).
  §4.1 needs a loop that judges every psp and refuses once, with the list. Most of that module —
  the stem rule `psp_beside` (l.263), the identity check, `CensusInCohort` — goes with the sidecar.
- **The lazy census reader is path-based** ([`backed`, `census.rs:1280`, `1326-1330`](../../../../src/ng/parameter_estimation/joint/census.rs)).
  §5: path plus a base offset, added to every seek. The resident reader, `read_census`, is a
  parity oracle for it and is unaffected.
- **An old-version census is `Malformed` today** ([`census_file.rs:456-461`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).
  §4.2 wants the version word read and reported before that check.
- **`PileupIdentity` is already `Option`** in the census header
  ([`CensusFile.pileup`, `census_file.rs:224`](../../../../src/ng/parameter_estimation/joint/census_file.rs)),
  so writing it absent needs no format change; the type and `freshness`/`freshness_by_header` are
  deleted once nothing reads them.
- **Tests that look for `<sample>.census`.** `generate_psps/tests.rs` l.479-487, 699, 975, 1029
  and the block at 1128; `census_from_psp`'s and `census_cohort`'s; `estimate_parameters/tests.rs`
  builds `--census` argument lists (`args_over`, l.58); `generate_census/tests.rs` — fifteen tests
  — moves whole to the renamed command, and its byte-for-byte test (l.111) is the oracle that the
  move changed nothing.
- **Everything that spells the old shape.** [`cli.rs:88-107`](../../../../src/pop_var_caller_exp/cli.rs)
  and its subcommand docs; `SUBCOMMAND` in both command modules;
  [`scripts/ng_fit_stage_end_to_end.sh:5-9`](../../../../scripts/ng_fit_stage_end_to_end.sh), which
  calls `generate-census`, passes `--census`, and diffs census files; check
  `scripts/ng_census_route_cost.sh` and `scripts/ng_census_agreement_mutations.sh`;
  `PROJECT_STATUS.md`'s pipeline line; and `examples/ng_census_route_cost.rs`, which writes census
  files through `write_psp(path, Some(census))`.

---

## 10. Cross-cutting

**Memory.** Step 2 reads one header and one footer per psp, then census sections lazily at the
trailer's offset — as today, plus one offset. The walk holds its census until the end, as today.
`regenerate-census` holds one psp open at a time.

**Errors.** Every refusal names the sample, the file and the cause; none is a panic. A refusal at
step 2 leaves nothing changed on disk. `regenerate-census` leaves the samples already done whole; a
sample interrupted inside `replace_trailer`'s window is a psp with no footer, which every reader
refuses and the report names — the accepted cost of §3.

**Concurrency.** Nothing new. Both commands are per-sample and are spread by invocation.

---

## 11. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the trailer as a replaceable payload | `PspWriter::finish(trailer)`, [`writer.rs:672`](../../../../src/ng/psp/writer.rs); `replace_trailer`, [`trailer.rs:103`](../../../../src/ng/psp/trailer.rs); the footer's trailer extent | called as they are; the census is the payload |
| the census's bytes | `write_census` / `decode_census`, [`census_file.rs:237`, `:456`](../../../../src/ng/parameter_estimation/joint/census_file.rs) | unchanged; written into memory instead of a file, read at an offset |
| the settings a census depends on | `SegmentationInputs` in the header, [`segmentation_inputs.rs:25`](../../../../src/ng/segmentation_inputs.rs) | read instead of the deleted flags |
| whether a cohort agrees on them | `SegmentationInputs::first_difference`, [`segmentation_inputs.rs:55`](../../../../src/ng/segmentation_inputs.rs) | called for the first time outside its tests, from `OpenPspCohort::open` |
| listing a cohort of psps and refusing one that is not one | `OpenPspCohort::open`, [`psp_caller.rs:111`](../../../../src/ng/run/psp_caller.rs) | every command opens its cohort here |
| rebuilding one census | `census_from_psp`, [`census_from_psp.rs:202`](../../../../src/ng/run/census_from_psp.rs) | called as is by `regenerate-census`, its identity argument dropped |
| the selection digest | `fit_a_cohort`'s check, [`census_fit.rs:116-139`](../../../../src/ng/run/census_fit.rs) | the third freshness row; shared with `regenerate-census`'s skip |
| the parity oracle | `each_census_it_writes_equals_the_one_the_walk_wrote`, [`generate_census/tests.rs:111`](../../../../src/pop_var_caller_exp/generate_census/tests.rs) | becomes: a psp copied, regenerated, and identical to the original **byte for byte, whole** — header, blocks, index, trailer and footer |

---

## 12. Deferred, with a recommended home

- **A footer rebuilt from an intact index**, so a psp torn inside `replace_trailer`'s window costs a
  tail rewrite rather than a re-walk. The index has no magic of its own — only the footer says where
  it starts ([`footer.rs:19-33`](../../../../src/ng/psp/footer.rs)) — so this is a scan, not a
  read. [`psp_file_format.md`](psp_file_format.md) §6.5.
- **A fit-side budget or cap**, served by filtering a generously stored census (§7) —
  [`parameter_prepass_runs.md`](../impl_plan/parameter_prepass_runs.md) Checkpoint C.
- **The census's per-position rule under a widened record.** A record wider than one base — one
  the sample's own deletion widened — contributes its depth at every kept position it covers and
  its alleles at none ([`add_generic`, `census.rs:2439-2510`](../../../../src/ng/parameter_estimation/joint/census.rs)).
  Measured 2026-09-09 with `examples/ng_census_locus_spans.rs`: 8,855 of 367,599 positions
  carrying non-reference evidence on HG002 at 50×, and 2,270 of 75,606 on the tomato accession —
  about 1 in 35 — read back as covered with nothing disagreeing.
  [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §2 says a spanning
  deletion lands in the fifth allele bucket; the code never reaches that bucket for one. Not this
  document's to settle.

---

## 13. Resolved decisions

- The census is the psp's trailer; there is no census file — §3. *Rejected:* a file beside the psp,
  anywhere; a file beside the psp, same stem only.
- Step 2 refuses a missing or stale census and reports every one, without a time estimate — §4.
  *Rejected:* rebuilding inside step 2.
- The selection's settings come from the psp header; six flags go — §6.
- A cohort disagreeing on catalog or criteria is refused by every command that opens one — §6.
- No `--skip-census` at step 1 — §3.3.
- `regenerate-census` replaces `generate-census`, writes into the psp, has no `--output-dir` and
  no `--force`, and skips fresh psps by the full freshness judgement — §8.
- Appending to a psp has no user case and stays off the command surface — §3.4.

**Open.** None.
