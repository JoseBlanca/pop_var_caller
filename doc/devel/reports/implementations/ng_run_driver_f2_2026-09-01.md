# F2 — every run writes the parameters it used, beside its VCF

**Date:** 2026-09-01. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone F step F2. **Design:** [`parameters_file.md`](../../ng/spec/parameters_file.md) §7,
with §3.2, §6 and §8 in the margins.

---

## What landed

`call-from-alignments` now writes two files: the VCF it was asked for, and beside it the
parameters that produced it. Nothing turns this off — spec §7 makes writing unconditional, *"the
run that most needs its parameters recorded is the one whose operator did not think to ask"* —
and a `--defaults` run writes one exactly as a run given a file does.

```
$ pop_var_caller_exp call-from-alignments --reference ref.fa \
      --alignment a.cram --alignment b.cram --output calls.vcf.gz --defaults
output	calls.vcf.gz
parameters	calls.parameters.toml
  groups_of_numbers_the_file_says_were_fitted	0 of 7
…
```

Almost all of the machinery existed: `ParametersFile::of_run`, `write_beside_the_vcf`,
`beside_the_vcf`, `ReadsBehindEachCalibration`'s three constructors,
`DeclaredInbreeding::of_each_sample`, `CensusIdentity::of_a_run_with_no_census`. **What F2 built
is the wiring and the two decisions those pieces deliberately left to a run driver.**

### The two decisions

**When the file is assembled, and when it is written — they are different moments.**

`ParametersFile::of_run` holds its checks in release and its own note says why the driver has to
choose: *"this runs after the last locus, so it discards a cohort's calling work where
`RunParameters::assemble`'s equivalent checks discard a startup… Still open, and it is the run
driver's."* Every one of those checks is a startup question — that this run's read-group table,
its parameters and its inbreeding estimates were all minted from the same inputs — and nothing
about the file changes while the run calls: it records what the run was **configured with**, not
what it found. So the file is **assembled before the first read is decoded**, where a panic costs
a startup, and **written after the VCF is renamed into place**, because spec §7's three purposes
(a run reproducible from its own output, a defaults run auditable, an edit that starts from
something) are all about a run that finished, and a parameters file standing beside a VCF that
does not exist answers none of them.

**A run may not write its parameters over the file it was handed.**

The second thing `write_beside_the_vcf` leaves to the driver: *"Spec §7 tells a user to copy the
file their run wrote and edit a line, and a re-run whose supplied file and whose VCF share a stem
writes over the file it was handed. Whether a driver should refuse that is the driver's."*
`--parameters calls.parameters.toml --output calls.vcf.gz` is the natural next command after a
first run, and it would destroy whatever the person changed by hand while the numbers came back
looking ordinary. **Refused before anything is read**, naming both files and what to do instead.
The comparison resolves each path's directory, so `<dir>/sub/../calls.parameters.toml` does not
slip past it.

### What else changed

**`##parametersFile` is filled**, which F1 deliberately left off because it had no file to name.
By **name and not by path**: the two are siblings by construction, so a path would say the same
thing at greater length and would be wrong the moment somebody moved the pair.

**`run_parameters` now hands back three things rather than one** (`TheRunsNumbers`).
`RunParameters` keeps what *calling* reads — one bare multiplier a read group, one bare
coefficient a sample — and a file has to say how much data stood behind each and under what
warrant. Neither is recoverable from the parameters afterwards, so both travel from wherever the
numbers came from. On the defaults path they are `ReadsBehindEachCalibration::nothing_was_fitted`
and `DeclaredInbreeding::of_each_sample`, built from the same two arguments
`RunParameters::of_defaults` takes so the two cannot disagree; on the supplied path they are what
the file recorded, because **a run that read a file writes back what it read** — a run that
dropped them would write a file claiming its supplied numbers rest on nothing.

**The summary names the second file and says what it rests on.** A file a fit wrote and a file a
defaults run wrote are the same shape by design, which is what makes a defaults run auditable and
also what makes the two look alike on the page; `what_the_run_fitted` is the file's own answer,
and the run prints it as a count of groups.

---

## What was measured

In the development container, at this step's working tree.

| check | before F2 | after F2 |
|---|---|---|
| `cargo test --lib` | 5,857 passed, 13 ignored | **5,868 passed**, 13 ignored |
| `cargo fmt --check` | exit 0 | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 | exit 0 |
| `cargo doc --no-deps` | 26 unresolved-link errors, 23 redundant-link warnings | unchanged |

**11 tests added**: 10 on the command surface — two of them driving `run_call_from_alignments`
itself — and one in `parameters_file` pinning that the file the run writes is as readable as the
VCF beside it.

**Spec §7's first purpose is an assertion, not prose.**
`what_a_defaults_run_writes_reads_back_into_the_numbers_it_scored_with` builds a real read-group
table from two fixture alignment files, assembles the file `of_run` would write, writes it beside
a VCF path, reads it back through `from_toml` and `to_run_parameters_for`, and checks the ploidy,
the library and plant counts, every coefficient and every multiplier survive the round trip. It
goes back in through the *binding* door, so what it shows is that a run's own file passes its own
refusals — a file naming another reference or another cohort's read groups is turned down there.

`a_defaults_runs_file_says_it_fitted_nothing` pins the other half: 0 of 7 groups, with the
denominator asserted non-zero so the count cannot say nothing.

---

## What the review found and what was done about it

Two agents reviewed the uncommitted step in isolated worktrees. **The second built a 4.6 kb
reference, its catalog and two 24× samples and drove the whole thing end to end**, which is what
produced everything below.

**The round trip is exact and a hand edit reaches the calls.** Feeding a run's own parameters file
back with `--parameters` over the same cohort gave a byte-identical VCF body *and* a
byte-identical parameters file. Setting one sample's inbreeding coefficient to 0.9 by hand and
marking it `supplied` moved exactly one field of one record — that sample's homozygote from
**GQ 63 to GQ 76**, the right sample and the right direction — and the run wrote the edit back
faithfully, dropping the "nobody said how inbred this sample is" comment that was no longer true.

**Four defects, all fixed here:**

- **The summary's count claimed a fit that never happened.** The line read
  `groups_of_numbers_fitted_from_this_cohort`, and `call-from-alignments` fits nothing — the
  number is repeated back from whatever a supplied file asserts. Marking one calibration row
  `fitted_here` by hand in a file a `--defaults` run had produced made the summary say `1 of 7`
  when no number in that run came from the cohort's reads. `what_the_run_fitted`'s own
  documentation says so in as many words, and the file's own banner phrases it correctly. The
  label is now `groups_of_numbers_the_file_says_were_fitted`.
- **The parameters file was written mode `0600` beside a VCF at `0644`.** Measured on every run:
  a colleague on a shared directory could read the calls and not the file saying what they rest
  on, which defeats §7 for everyone but the person who launched the run. It came from
  `NamedTempFile`, which creates `0600` and keeps it through the rename; the write now uses
  `File::create` and a rename, which is what the VCF's own sink does, so both files take their
  mode from the process's `umask` rather than from which crate created them. The test asserts the
  two **agree**, rather than asserting a fixed number.
- **A plain `--defaults` re-run destroyed a hand-edited parameters file in silence.** The new
  refusal guards the `--parameters` route; this is the other route to the same loss, where no flag
  names the file. Writing stays unconditional — §7 leaves no room for a switch — but a run that
  replaces an existing parameters file now says so on stderr.
- **A run whose VCF was complete and correct exited 1 saying only that the parameters could not
  be saved.** In a `set -e` pipeline that VCF is thrown away, which is the opposite of what the
  message meant. The run still fails — it did not produce all of its output — but the message now
  names the finished VCF and says to keep it.

**And one thing the file said about itself that was not true.** Of
`fallback_length_spectrum_concentration` it claimed *"a run handed this file by another run marks
it `supplied`"*; the round trip writes back `defaulted`, and the **code is right**: spec §2.1
settles that a supplied file keeps its warrants, because demoting on every read is demoting on
every direct-mode run under another name and would break the two-mode oracle. The sentence now
says what §2.1 says. The golden file `every_shape_as_written.toml` was regenerated and its diff is
that sentence and nothing else.

**Recorded and not acted on:**

- **The per-sample and per-read-group comment blocks repeat verbatim, once a row.** At two samples
  that is the file's main virtue; at the tomato cohort's 63 it is 378 identical lines, and at the
  committed ceiling of 3,000 samples about 18,000 lines of one repeated paragraph wrapped around
  3,000 numbers. Spec §9 already prices this axis as *"the largest artefact a run writes that is
  not a VCF"* at a thousand samples and says the shape is not settled; stating the paragraph once
  above the table is the obvious move and it is `to_toml`'s, not this step's.
- **Nothing tells a person their hand edit reached the calls** beyond the calls themselves. The
  new fitted-groups line covers the fit-versus-default question and not the edited-versus-shipped
  one.
- **266 comment lines to 50 data lines.** The reviewer's judgement was that this is right for a
  first reader and awkward for a fiftieth, and that they would not trade the prose away.

## The correctness pass, and the one thing it found that was wrong

19 mutations across two passes. **The first pass killed 5 of 14, and every one of the nine
survivors was the same gap: nothing in the tree called `run_call_from_alignments` or `report`.**
All the tests reached the private helpers directly, so deleting the refusal's *call site* left
17 of 17 green; so did a header naming `calls.vcf` instead of `calls.parameters.toml` on every
VCF a run writes; so did handing `ParametersFile::of_run` an axis of the wrong length, which
panics at startup on any real run. **Two end-to-end tests now drive the command itself** over a
reference, a catalog and two samples built on disk, and close that cluster.

**And one real defect: a `--parameters` run wrote a file saying it had no census.**
`ParametersFile::of_run`'s own contract says the opposite for that path — *"a run that read a
file fitted under other terms and writes its parameters out again has to write back the terms it
read, not its own"* — and the first draft passed `CensusIdentity::of_a_run_with_no_census()` on
both paths. The consequence is not a wrong number but a **silent loss of provenance one hop
through direct mode**: a psp fit writes a file with its census terms and `fitted_here` warrants;
somebody re-runs it in direct mode; the new file keeps the warrants and records no census; a
later psp run over the same cohort and the same census now finds a disagreement and demotes every
number to `supplied`, where reading the *original* file would have kept it. That is precisely the
two-mode divergence spec §2.1 exists to prevent, and the file was internally contradictory
besides — numbers marked as fitted from reads, beside a record saying no census produced them.
The census now travels with the other two things `TheRunsNumbers` carries.

**Three further findings, all fixed:**

- **A symlinked `--parameters` slipped past the refusal**, probed: `--parameters handy.toml`
  pointing at `calls.parameters.toml` was admitted, and the run then renamed its output over the
  link's target. Both sides are now resolved through the file system rather than by their
  directories alone, which also closes the case-insensitive-volume gap — on macOS's default
  `CALLS.vcf.gz`'s sibling *is* `calls.parameters.toml`, and a byte-wise name comparison said
  otherwise.
- **A `--parameters` file that does not exist was refused with the wrong message**, telling the
  person to copy a file that is not there. Only a file that exists can be overwritten, so only
  one that exists is compared; a mistyped name now gets the message about the file not existing.
- **The round-trip test was order-blind.** Its two samples both scored at the default
  coefficient, so permuting the per-sample axis passed — which a surviving mutation made visible.
  One sample now selfs and one outcrosses, with the fixture asserting the two differ before it
  compares them.

**One gap left open and recorded.** The summary's own text is still unpinned: a mutation that
swapped its two numbers, or printed the VCF's path on the `parameters` line, survives. Pinning it
means either capturing stdout or returning the lines rather than printing them, and the summary
is what every read-the-artefact pass looks at first — so it is recorded here rather than closed
by a refactor.

## Deviations from the plan, absorbed

**One, and it is the census above.** The plan's F2 is one sentence and this is otherwise it; the
two decisions in *What landed* are ones the code's own documentation asked the driver to make.
The census was not a decision at all — it was a field written wrongly, found by review, fixed.

## What F2 does not do

- **A finished VCF can name a parameters file that is not there.** The header is fixed before
  the walk and the parameters are written after `finish()`, so a failed parameters write ships a
  complete VCF whose `##parametersFile` dangles. `ParametersNotWritten` says so and names both
  files; the alternative — writing the parameters first — trades this for a parameters file
  beside a VCF that does not exist, which §7 answers none of its purposes with.
- **No `##fileDate`, and nothing in the VCF header says `--defaults` in as many words.** What the
  header does say is which file to open, and that file's own opening line says what was fitted.
  Whether the VCF should repeat it is the emission spec's, not this step's.
- **The file is built as one `String` before any of it is written.** `write_beside_the_vcf`'s own
  note prices that at up to 79 MB at 3,000 samples, taken after the last locus is called;
  nothing at the tomato cohort's 63. Recorded there, not changed here.
