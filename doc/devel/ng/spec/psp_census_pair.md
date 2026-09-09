# ng — the psp and its census are a pair

*Design spec, 2026-09-09. **No code yet — this settles the design.** Companion plan:
[`impl_plan/psp_census_pair.md`](../impl_plan/psp_census_pair.md). It amends one sentence of
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1 — the one saying a
stale census is rebuilt silently at the fit — which was never built; the code refuses, and this
document makes that refusal the design.*

*Vocabulary, once. A **psp** is one sample's stored evidence: what its reads showed at every position
of the ground that was walked. A **census** is the much smaller file a parameters fit reads: what
that sample showed at a fixed set of positions and repeat tracts chosen for the whole run, so the fit
can ask the same question of every sample. The **selection** is that fixed set. A census is
**fresh** when it was built from the psp beside it, under the selection the run will rebuild, in the
format this build reads; otherwise it is **stale**. The **stem** of a file is its name without the
extension: `HG002.psp` and `HG002.census` share the stem `HG002`.*

---

## 1. What it is

**A rule about files, and three commands that follow it.** The rule: every psp has its census
beside it, same directory, same stem, and there is nowhere else a census can be. The commands a
user runs, in order:

| step | command | in | out |
|---|---|---|---|
| 1 | `generate-psps` | alignment files | `<sample>.psp` **and** `<sample>.census`, always both |
| 2 | `estimate-parameters` | psps | the parameters file a calling run scores with |
| 3 | `call-from-psps` | psps and that file | the VCF |

And one repair command, run only when step 2 says to: **`regenerate-census`** — psps in, a fresh
census beside each.

**What changes from today.** Step 1 already writes both files from one pass over the reads
([`generate_psps.rs:676-735`](../../../../src/pop_var_caller_exp/generate_psps.rs)). Step 2 today
takes `--census` paths and six flags that must equal what the psps were walked under
([`estimate_parameters.rs:69-141`](../../../../src/pop_var_caller_exp/estimate_parameters.rs)); it
will take `--psp` and read those settings out of the psp headers. Today's `generate-census` takes an
`--output-dir` that may be anywhere ([`generate_census.rs:97`](../../../../src/pop_var_caller_exp/generate_census.rs));
`regenerate-census` writes beside the psp and nowhere else.

**Goals.**

- **The user names psps and nothing else, at every step.** No command takes a census path.
- **A stale or missing census stops step 2 with a report that names every affected sample, the
  reason, and the command to run.** It does not rebuild on its own.
- **Nothing a census depends on is typed twice.** The repeat criteria, the catalog and the analysed
  ground come from the psp header; the flags that duplicated them go.

**Non-goals.** Storing the census inside the psp — considered and rejected, §3. Rebuilding a stale
census inside step 2 — considered and rejected, §4. A flag to skip the census at step 1 — rejected,
§3.3. Storing a larger census than the fit uses — deferred, §7. What the census *holds*, how it is
encoded, and the fit itself — unchanged, and settled in their own documents.

**It does not:** open an alignment file; read a psp beyond its header at step 2; write anywhere but
beside a psp; accept a census path from anyone.

---

## 2. Why the census is kept at all — the cost of not having one

**The census is a cache, and rebuilding it means reading the whole psp.** A census can be rebuilt
from its psp with no alignment file — that is what today's `generate-census` does — but the rebuild
decodes every record of the psp ([`census_from_psp.rs:241-252`](../../../../src/ng/run/census_from_psp.rs)).
It cannot skip records by their head, because the read-error calibration totals it accumulates are
summed over every generic locus, and the per-read quality sums they need are in the record body
([`census.rs:2383-2396`](../../../../src/ng/parameter_estimation/joint/census.rs)).

**Measured, 2026-09-09**, with `examples/ng_census_read_vs_psp.rs` (branch `census-vs-psp-perf`;
the plan lands it), one thread, files in the page cache, the census at the density a whole-genome
run has — two million positions over the genome, which is 1 in 400 on tomato and 1 in 1,500 on
human:

| sample | ground | reads a position | psp | read the census | rebuild it from the psp |
|---|---|---|---|---|---|
| tomato SRS3394606 | 32 Mb | 22.6 | 10.4 MB | 0.0002 s | 0.25 s |
| tomato DRR000741 | 10 Mb | 98.5 | 158 MB | 0.0002 s | 3.0 s |
| HG002, the 50× set | 4.13 Mb | 33.5 | 43.7 MB | 0.0002 s | 1.41 s |
| HG002, the 300× set | 5.08 Mb | 300.9 | 126 MB | 0.0001 s | 3.25 s |

The rebuild runs at about 45 MB of psp a second whichever corner it is in. **At 50× a human genome
is about 32 GB of psp a sample, and its rebuild about a quarter of an hour a sample**, on one
thread: 0.34 s a megabase over 3.0 Gb. Reading the census instead is a few milliseconds.

Nothing below changes those numbers. They are why the census exists, why step 2 must not rebuild it
without being asked, and why the rebuild is a command of its own that can be spread over machines.

---

## 3. Where the census lives — decided: beside the psp, same stem, nowhere else

Three options were live. **(a)** Inside the psp, in its trailer — the format has a payload slot at
the end of the file that ng writes zero bytes into today ([`psp_writer_line.rs:244`](../../../../src/ng/run/psp_writer_line.rs));
`replace_trailer` rewrites it by truncating at the trailer offset and writing forward, leaving the
blocks untouched ([`psp/trailer.rs:60-103`](../../../../src/ng/psp/trailer.rs)). **(b)** Beside the
psp, in any directory the user names — today's `generate-census --output-dir`. **(c)** Beside the
psp, in its directory, same stem — `generate-psps`'s own convention, made the only one.

**Decision (owner, 2026-09-09): (c).**

**Why not (a).** The rebuild costs the same either way — §2's full pass — so what separates them is
what happens to the 32 GB file when the 30 MB census is rewritten. With a sidecar, the psp is never
opened for writing. With the trailer, `replace_trailer` truncates the psp before it writes, so
until the new footer lands the file has no footer and every reader refuses it; a failed truncation
means the sample has to be walked from its reads again ([`psp/trailer.rs:66-79`](../../../../src/ng/psp/trailer.rs)).
And on storage that is read-only — an archive, an object store, which is where 32 GB files end up
— a trailer cannot be rewritten at all, where a sidecar can be written wherever the psp can be
copied. The owner's words: what the sidecar buys is *safety and reach*.

**Why not (b).** A fit handed census paths cannot tell *this sample's census is missing* from *this
sample is not in the cohort*: it fits the censuses it was given, one sample short, and nothing says
so. Handed psps, it knows the cohort, and "is the census beside this psp fresh?" is a question it
can ask per sample (§4). A census that may live anywhere also has to be *named* somewhere, and that
naming is the flag this design removes.

**What (c) costs.** A sample is two files to move instead of one. `sample.*` is the unit. A pair
copied apart is not silent: the census carries the digest of its psp's header, and the fit checks
it (§5).

### 3.1 The walk always writes both

`generate-psps` writes the census from the same pass that writes the psp, and building it there is
cheaper than walking without it and rebuilding afterwards — 1.28 s against 1.40 s on six tomato
accessions over 200 kb ([`ng_fit_stage_b_2026-09-05.md`](../../reports/implementations/ng_fit_stage_b_2026-09-05.md)).
A whole-genome census is a few megabytes of positions plus tens of megabytes of tracts — 1.03 bytes a
position and about 35 a tract, measured on the tomato accession's census at every position kept —
against a psp of gigabytes.

### 3.2 A pair is written census first, psp second — unchanged

Both files are written to a `.partial` name and renamed once whole; the census is renamed first,
so a run that dies between the two renames leaves a new census beside an old psp, which the
identity check refuses, rather than a finished psp beside a stale census, which would be trusted
([`generate_psps.rs:713-735`](../../../../src/pop_var_caller_exp/generate_psps.rs)). Nothing here
changes that; `regenerate-census` writes its one file the same way.

### 3.3 No `--skip-census` — decided

Rejected because it saves a tenth of a percent of the psp's size and a small share of the walk, and
manufactures the one state this design exists to remove: a psp that looks complete and costs a
quarter of an hour a sample to make usable, discovered at step 2. `generate-psps` has no such flag
today and keeps none. Someone who wants the space back deletes the file, and step 2 will say so.

---

## 4. When the pair is not whole — decided: step 2 refuses and reports, it does not rebuild

Two options were live. **Rebuild inside step 2** when a census is missing or stale, silently —
what [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1 wrote, and
what the owner first proposed for this design. **Refuse**, naming every sample that needs
regenerating and the command that does it.

**Decision (owner, 2026-09-09): refuse.** "If we can't run the estimate, we stop, inform the user,
and the user has to decide to regenerate the census."

**Why.** The rebuild is a quarter of an hour a sample and invisible: an `estimate-parameters` over
thirty human samples that quietly rebuilt would run for eight hours with nothing to show that
anything unusual was happening. It is also per-sample and distributable, and a rebuild inside step
2 would run serially on the one machine step 2 runs on. Refusing makes the cost a decision.

### 4.1 Every sample is examined before anything is refused

Today `open_census_cohort` returns at the first census it cannot open or check
([`census_cohort.rs:184-215`](../../../../src/ng/run/census_cohort.rs)): a cohort with three stale
censuses reports one, and the user finds the other two on the next run. **The whole cohort is
examined and the report names all of them**, grouped by cause, so one run tells the user the whole
regeneration job.

### 4.2 What makes a census stale, and where each is caught

| cause | what is compared | caught today at |
|---|---|---|
| no census beside the psp | the file exists | the cohort opener, one sample at a time |
| a census of an older format | its version word against `VERSION` ([`census_file.rs:75`](../../../../src/ng/parameter_estimation/joint/census_file.rs)) | `open_census`, as *malformed* — indistinguishable from a damaged file |
| built from another psp | the digest of the psp header it names against the psp beside it ([`freshness_by_header`, `census_file.rs:176`](../../../../src/ng/parameter_estimation/joint/census_file.rs)) | the cohort opener |
| built under another selection | the digest of the kept loci it carries against the loci the run rebuilds ([`fit_a_cohort`, `census_fit.rs:116-139`](../../../../src/ng/run/census_fit.rs)) | the fit, after the cohort has opened and the reference has been read |

The first three are judged from the two files' headers and cost nothing. **They are checked first,
before the reference is read** — reading the reference and rebuilding the selection is 4 to 19
seconds in the measurements of §2 — so a refusal for a missing census is immediate. The fourth
needs the rebuilt selection and stays where it is; its message changes to say what to run.

**An older format is named as such, not as damage.** A census that opens as *malformed* because its
version word is old tells the user to look for corruption; the truth is that this build changed
what a census holds. The version is read before anything else and reported as *built by an older
version; regenerate*.

### 4.3 What the refusal says

Per sample: the sample name, the census path, the cause in the words of §4.2. Then one line with
the command, built from what was given — `regenerate-census --psp <the same arguments>` — so it
can be copied. The samples that are fresh are counted, not listed.

---

## 5. How a psp and a census are known to be a pair

**The stem identifies; the digest verifies.** `<stem>.census` is the census of `<stem>.psp` and
of nothing else. Whether it really was built from that psp is what the census's stored identity
says: the MD5 of the psp's header bytes as they stand in the file, and the psp's record count
([`PileupIdentity`, `census_file.rs:83-92`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).

**Only the header half is checked, and that is a known gap.** A fit compares the header digest —
one short read of the psp — and not the record count, because the count sits inside the compressed
blocks and reading it means decompressing the psp at every fit
([`freshness_by_header`, `census_file.rs:176`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).
What that leaves unchecked is a psp whose header is unchanged and whose records are not, which
only `PspWriter::append` can produce, and no shipped command calls it. The owner, 2026-09-09: there
is no user case for appending; regenerate instead. **Deferred**: a record count in the psp footer
would close the gap for one cheap read; it is a psp format change and belongs to
[`psp_file_format.md`](psp_file_format.md).

---

## 6. Where the settings come from — decided: the psp header, never a flag

**Every setting the selection depends on is in the psp header already.** `SegmentationInputs`
carries the catalog's own header, the repeat criteria the walk routed under, and the analysed
regions ([`segmentation_inputs.rs:25-40`](../../../../src/ng/segmentation_inputs.rs)). Step 2 and
`regenerate-census` read them from there.

**And they must check that the cohort agrees on them, which nothing does today.** The cohort
openers compare the analysed regions only ([`census_cohort.rs:222-236`](../../../../src/ng/run/census_cohort.rs);
`OpenPspCohort`'s `AnalysedRegionsDiffer`). `SegmentationInputs::first_difference` names the first
of the three fields that differs, in the order a person should fix them — catalog, criteria,
regions ([`segmentation_inputs.rs:42-55`](../../../../src/ng/segmentation_inputs.rs)) — **and is
called only from its own tests.** A step 2 that took the criteria from the first psp's header and
never asked the rest would build a selection the other samples' censuses cannot match, and learn
it twenty seconds later from the fit's digest check, reported as *another selection*.

**Decision (owner, 2026-09-09): a cohort whose psps disagree on the catalog or the criteria is a
hard failure for every command that opens one** — `estimate-parameters`, `regenerate-census`, and
`call-from-psps` alike, since two psps typed under different criteria cannot be called together
either. The check goes where the analysed-regions check already is, `OpenPspCohort::open`
([`psp_caller.rs:677-690`](../../../../src/ng/run/psp_caller.rs)), which `call-from-psps` and
today's `generate-census` already call ([`call_from_psps.rs:494`](../../../../src/pop_var_caller_exp/call_from_psps.rs),
[`generate_census.rs:413`](../../../../src/pop_var_caller_exp/generate_census.rs)); step 2 opens
its cohort the same way. Written once, refusing with the sample and the field named, before
anything else is read.

**Flags removed from `estimate-parameters`:** `--census`, `--min-copies`, `--min-period`,
`--max-period`, `--max-str-len`, `--min-purity`. **Kept:** `--reference` (the FASTA's path is not in
the header, and the selection needs the reference to know where it is sequence at all),
`--catalog` (defaulting to `<reference>.repeats.parquet`, and checked against the header's catalog
digest as today), `--output`, `--force`, `--ploidy`, `--inbreeding`. **Added:** `--psp`, once per
file or a directory, as `call-from-psps` and today's `generate-census` take it. `regenerate-census`
takes the same set less `--output`, `--force`, `--ploidy` and `--inbreeding`.

**Why the criteria in particular cannot be a flag.** The repeat criteria decide how the ground is
partitioned: a stretch is a tract segment, a bundle, a satellite, or `Generic`, which is what is
left over ([`region_typing/mod.rs:279-299`](../../../../src/ng/region_typing/mod.rs)). A position
inside a tract is emitted once, as part of that tract's locus; there is no per-position generic
evidence for it anywhere in the psp. So the criteria are a property of the psp — change them and
the psp is wrong, not the census — and a flag that has to equal a header field is a way to build a
census that silently disagrees with the file it belongs to. The census selection is cut by the
same criteria as the psp ([`select_kept_loci`, `loci.rs:827-865`](../../../../src/ng/parameter_estimation/joint/loci.rs)
takes its generic domain from the catalog's `Generic` segments), and that is the only census a psp
can support.

**Trap for the coder.** `estimate-parameters` builds its ground from the five flags as a
`run_ground::RepeatRouting` and takes only the analysed regions from the cohort
([`estimate_parameters.rs:312-322`](../../../../src/pop_var_caller_exp/estimate_parameters.rs)).
`segments_over` takes that `RepeatRouting` and converts it to a `StrRepeatCriteria` itself, through
`routing_criteria` ([`run_ground.rs:242-259`, `270-281`](../../../../src/pop_var_caller_exp/run_ground.rs)).
The header holds the finished `StrRepeatCriteria`, and there is no conversion back. So the
segmentation needs an entry point that takes the criteria already made — split `segments_over`
after its call to `routing_criteria`, so the flag-driven commands and the header-driven ones share
everything past that line.

---

## 7. The selection's two counts, and why storing more would be free — deferred

The budget of positions and the cap of tracts per stratum are compiled-in constants: about two
million and five thousand, with a fixed seed ([`CensusSelection::SHIPPED`, `gatherer.rs:173`](../../../../src/ng/run/gatherer.rs)).
No command takes them as flags, so they cannot go stale by being retyped. Nothing here changes
them.

**Recorded because it decides what a later knob would cost.** The selection nests. A position is
kept when `hash(contig, position, seed) < threshold`, the threshold rises with the budget, and the
hash depends on the position and the seed alone ([`loci.rs:316-376`](../../../../src/ng/parameter_estimation/joint/loci.rs)):
the set kept at two million is a strict superset of the set kept at half a million, and a fit
wanting the smaller set filters the stored one and gets exactly the positions a fresh selection
would choose. The per-stratum cap is the lowest-hash prefix of the same kind
([`StratumSample`, `strata.rs:58-70`](../../../../src/ng/repeat_catalog/strata.rs)). So a census
stored at ten times the budget serves any smaller one without a rebuild, and costs about one part
in a hundred of the psp.

**Deferred, with a home.** A fit-side budget or cap needs a knob that does not exist, and the
question of which selection knobs become flags was reopened at
[`parameter_prepass_runs.md`](../impl_plan/parameter_prepass_runs.md) Checkpoint C and belongs
there. The nesting is what makes the answer cheap when it is wanted.

---

## 8. `regenerate-census` — the repair command

**Psps in; a fresh census beside each; nothing else written.** It replaces today's
`generate-census`, whose three stated cases — psps written before censuses existed, a census lost
or built under settings since changed, and wanting a larger census — are all *regenerate*.

- **`--psp`** once per file or a directory; `--reference`; `--catalog` defaulting as in §6. The
  criteria and the ground come from the headers.
- **No `--output-dir`.** The census goes beside the psp. A directory that refuses the write fails
  that sample with the path; there is no alternative location.
- **No `--force`.** Replacing is the command's job, and replacing a fresh census costs time and
  changes no bytes — the two producers agree byte for byte
  ([`each_census_it_writes_equals_the_one_the_walk_wrote`, `generate_census/tests.rs:111`](../../../../src/pop_var_caller_exp/generate_census/tests.rs)).
- **Fresh pairs are skipped, and the run says so**, unless **`--all`** is given. Freshness is the
  same judgement step 2 makes (§4.2, the first three rows) — one function with two callers — so a
  run stopped part-way and started again does only the samples still owed. `--all` exists for the
  case the judgement cannot see: the selection constants changed in a new build while the format
  version did not.
- **One psp open at a time**, as today ([`generate_census.rs:29-33`](../../../../src/pop_var_caller_exp/generate_census.rs)).
  Per-sample progress to stderr as each finishes, in the same words as the final report, as
  `generate-psps` does.
- **Spread it the way step 1 is spread**: one invocation per sample or per directory, on as many
  machines as there are. This is why it is a command and not a mode of step 2.

---

## 9. What will bite the coder

- **No shipped command checks that a cohort's psps agree on the catalog or the criteria.** Only the
  analysed regions are compared; `first_difference` has the check and no caller outside its tests
  (§6). Both commands need it.
- **The cohort opener stops at the first failure.** `open_census_cohort` uses `?` on `open_census`
  and `return Err` on a stale identity ([`census_cohort.rs:184-215`](../../../../src/ng/run/census_cohort.rs)).
  §4.1 needs a loop that judges every pair and refuses once, with the list.
- **The stem rule exists in both directions, in two files.** Census → psp is `psp_beside`
  ([`census_cohort.rs:263-276`](../../../../src/ng/run/census_cohort.rs)); psp → census is
  `census_path_for` ([`generate_psps.rs:847`](../../../../src/pop_var_caller_exp/generate_psps.rs)).
  Step 2 now goes psp → census. One helper, both callers.
- **An old-version census is `Malformed` today** ([`decode_census`, `census_file.rs:456-462`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).
  §4.2 wants the version read and reported before the magic-and-version check turns it into damage.
- **Order of work in step 2.** Today: open the cohort, then read the reference, then rebuild the
  selection. The pair judgement must come before the reference read, or a missing census is
  reported after twenty seconds of FASTA.
- **The fit's own refusal.** `fit_a_cohort` returns `AnotherSelection` when the rebuilt loci do
  not match the digest the censuses carry ([`census_fit.rs:139`](../../../../src/ng/run/census_fit.rs));
  the command renders it as *fitting the cohort*. It has to say *regenerate*.
- **The tests that build arguments.** `estimate_parameters/tests.rs` constructs `--census`
  argument lists (`args_over`, l.58) and has `a_census_without_its_psp_is_refused` (l.217); the
  new shape is *a psp without its census is refused, and the report names it*.
  `generate_census/tests.rs` — fifteen tests — moves whole to the renamed command; the byte-for-byte
  test (l.111) is the oracle that the rename changed nothing.
- **Everything that spells the old names.** [`cli.rs:88-107`](../../../../src/pop_var_caller_exp/cli.rs)
  and its subcommand docs; `SUBCOMMAND` in both command modules;
  [`scripts/ng_fit_stage_end_to_end.sh:5-9`](../../../../scripts/ng_fit_stage_end_to_end.sh), which
  calls `generate-census` and passes `--census`; check `scripts/ng_census_route_cost.sh` and
  `scripts/ng_census_agreement_mutations.sh`; `PROJECT_STATUS.md`'s pipeline line; and the §6.1
  sentence this document amends.

---

## 10. Cross-cutting

**Memory.** Step 2 reads one header per psp and no records — one short read a sample, as today.
`regenerate-census` holds one psp open at a time. Neither grows with the cohort beyond the list of
verdicts.

**Errors.** Every refusal names the sample, the file and the cause; none is a panic. A refusal at
step 2 leaves nothing on disk. A refusal inside `regenerate-census` leaves the samples already done
in place and whole, and the one that failed absent, never a stump (§3.2).

**Concurrency.** Nothing new. Both commands are per-sample and are spread by invocation.

---

## 11. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the pair's identity check | `freshness_by_header`, [`census_file.rs:176`](../../../../src/ng/parameter_estimation/joint/census_file.rs) | called as is; its verdict becomes one of §4.2's causes |
| the psp header, cheaply | `psp::header_and_its_digest`, [`psp/mod.rs:192`](../../../../src/ng/psp/mod.rs) | called as is |
| the settings a census depends on | `SegmentationInputs` in the header, [`segmentation_inputs.rs:25`](../../../../src/ng/segmentation_inputs.rs) | read instead of the deleted flags |
| whether a cohort agrees on them | `SegmentationInputs::first_difference`, [`segmentation_inputs.rs:55`](../../../../src/ng/segmentation_inputs.rs) | called for the first time outside its tests |
| listing a cohort of psps and refusing one that is not one | `OpenPspCohort::open`, [`psp_caller.rs:111`](../../../../src/ng/run/psp_caller.rs), and `generate-census`'s directory expansion | the same expansion; the cohort agreement checks it makes |
| rebuilding one census | `census_from_psp`, [`census_from_psp.rs:202`](../../../../src/ng/run/census_from_psp.rs) | called as is by `regenerate-census` |
| writing beside the psp, whole or not at all | the `.partial`-then-rename in [`generate_psps.rs:676-735`](../../../../src/pop_var_caller_exp/generate_psps.rs) | the same shape for one file |
| the parity oracle | `each_census_it_writes_equals_the_one_the_walk_wrote`, [`generate_census/tests.rs:111`](../../../../src/pop_var_caller_exp/generate_census/tests.rs); [`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh) | must keep holding after the rename |

---

## 12. Deferred, with a recommended home

- **A record count in the psp footer**, so the pair check covers both halves of the identity for
  one read — [`psp_file_format.md`](psp_file_format.md) §3.3, the footer.
- **A fit-side budget or cap**, served by filtering a generously stored census (§7) —
  [`parameter_prepass_runs.md`](../impl_plan/parameter_prepass_runs.md) Checkpoint C's knobs
  question.
- **The census's per-position rule under a widened record.** A record wider than one base — one
  the sample's own deletion widened — contributes its depth at every kept position it covers and
  its alleles at none ([`add_generic`, `census.rs:2439-2510`](../../../../src/ng/parameter_estimation/joint/census.rs)).
  Measured 2026-09-09 with `examples/ng_census_locus_spans.rs` (same branch as §2's probe): 8,855
  of 367,599 positions carrying non-reference evidence on HG002 at 50×, and 2,270 of 75,606 on the
  tomato accession — about 1 in 35 — read back as covered with nothing disagreeing.
  [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §2 says a spanning
  deletion lands in the fifth allele bucket; the code never reaches that bucket for one. Not this
  document's to settle; it belongs to that spec.

---

## 13. Resolved decisions & open questions

**Resolved.**

- The census lives beside the psp, same stem, nowhere else — §3. *Rejected:* the psp's trailer;
  any directory the user names.
- Step 2 refuses a stale or missing census and reports every one — §4. *Rejected:* rebuilding
  inside step 2.
- The selection's settings come from the psp header; six flags go — §6.
- No `--skip-census` at step 1 — §3.3.
- `regenerate-census` replaces `generate-census`, writes beside the psp, has no `--output-dir`
  and no `--force`, and skips fresh pairs unless `--all` — §8.
- Appending to a psp has no user case and stays out of the command surface — §5.

- A cohort whose psps disagree on the catalog or the repeat criteria is refused by every command
  that opens one, `call-from-psps` included — §6.
- The refusal at step 2 does **not** estimate how long regeneration will take — owner, 2026-09-09:
  it depends too much on the hardware. It names the samples and the command, nothing more.

**Open.** None.
