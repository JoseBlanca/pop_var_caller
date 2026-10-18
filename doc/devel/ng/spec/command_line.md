# The `pvc` command line: what users do, and what the interface should be

*Design spec, 2026-10-15; decisions settled 2026-10-16. **No code changes yet.** It inventories
the command line as it stands on `main` (e4bd1e94), maps it onto what users of the caller
actually do, and settles what to keep, rename, move under `pvc dev` and remove, and what a run
prints and writes beside its output. The owner's rulings are in §5; §6 lists the changes that
needed none. This is the specification the implementation and the user documentation are written
from.*

---

## 1. What this is

`pvc` grew one subcommand at a time, each written for the experiment or milestone that needed
it. The result works, but it was never designed as one interface. There are nine subcommands
and about 114 options, a user has to read the source documents to know which subcommands are
for them, and the same idea is spelled differently from one subcommand to the next.

This document starts from the user's side: what a person running the caller wants to do (§2),
and the commands each of those tasks needs today (§3). It then lists every subcommand and
option with a proposed place for it (§4). The owner's rulings are in §5.

**Two readers matter for the interface, and they need different things:**

- **A user** — a geneticist calling a cohort. They need few subcommands, each named for a task,
  sensible defaults, and the less common settings kept out of the way.
- **A developer of the caller** — who needs the diagnostic tools and the experiment settings.
  These should stay callable, but not crowd the help a user reads.

## 2. What users do

The caller can run in two ways, and the tasks follow from them:

- **Direct mode** reads every sample's alignment files (BAM or CRAM) and writes a VCF in one
  run. Nothing is stored in between. Every sample's file is open for the whole run.
- **psp mode** splits the same work in two. First each sample's alignments are read once and
  what they show is stored in a *psp* file, one per sample. Then the cohort is called from those
  files. A sample added later costs one new psp; the other samples are not read again. Fitting
  the model's parameters (base-quality calibration, contamination, inbreeding, population
  diversity, repeat stutter) is only possible in psp mode, because the fit reads a small sample
  of positions every psp stores in its tail.

The tasks, the owner's four first:

| | Task | Mode |
|---|---|---|
| T1 | Call one or several samples from their alignments | direct |
| T2 | Store one or several samples' evidence as psp files | psp, first stage |
| T3 | Call one or several samples from their psp files | psp, second stage |
| T4 | Estimate the model's parameters from psp files | psp |
| T0 | Prepare a reference: build the repeat catalogue once | both |
| T5 | Find out why a locus was or was not called (developers) | both |

T0 is not a task a user sets out to do, but every user meets it as a required first step: every
run needs the repeat catalogue, a file listing where the reference's tandem repeats are, built
once per reference.

T5 is a developer's task, not a user's. It is how a developer traces why a locus was or was not
called, to check whether something in the caller went wrong. Its tools therefore belong with the
other developer tools (§5.3).

## 3. The tasks, as commands today

**T0 — prepare the reference.**

```
pvc repeat-catalog --reference genome.fa
```

Writes `genome.fa.repeats.parquet` beside the reference. Every other command finds it there.

**T1 — call from alignments.**

```
pvc call-from-alignments --reference genome.fa \
    --alignment s1.bam --alignment s2.bam \
    --output calls.vcf.gz --defaults
```

`--defaults` or `--parameters <file>` is required. **A direct-mode run cannot fit its own
parameters**, so without a parameters file it calls with the constants compiled into the binary:
no base-quality calibration, no contamination, no inbreeding. Today it also writes
`calls.parameters.toml` beside the VCF, holding those numbers; §7 replaces that file with the
run's log. (§5.1 is the decision this raises.)

**T2 — store evidence as psps.**

```
pvc generate-psps --reference genome.fa --alignment s1.bam --output-dir psps/
```

One `<sample>.psp` per sample, named from the alignment file's `@RG SM` tag. Samples are read
one after another; a cohort is spread over cores or machines by running one command per sample.

**T4 — estimate parameters.**

```
pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
```

**T3 — call from psps.**

```
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --output calls.vcf.gz
```

**T5 — explain a call (developers).** `--explain-loci <regions.bed> <out.tsv>` on either calling command
writes, for every locus in the BED, the candidates, the genotype likelihoods, the quality and
its penalties, and what became of the record. `pvc inspect-psp` dumps what one psp stored over a
BED, to tell whether a variant's reads ever reached the file.

## 4. Every subcommand and option, with a proposed place

**Placement** is one of:

- **keep** — in the main help, as it is;
- **advanced** — kept, but listed under an *Advanced* heading in the help;
- **dev** — a subcommand moved under `pvc dev` (§5.3), or an option listed under a *Developer* heading at the end of its subcommand's help;
- **change** — kept, with the change stated;
- **remove**.

### 4.1 Subcommands

| Subcommand | Task | What it does | Proposed |
|---|---|---|---|
| `repeat-catalog` | T0 | Scans the reference once for tandem repeats and writes the catalogue | keep |
| `call-from-alignments` | T1 | Direct mode: alignments → VCF | keep |
| `generate-psps` | T2 | Reads each sample's alignments and writes its psp | keep |
| `estimate-parameters` | T4 | Fits the model's parameters from psps → parameters file | keep |
| `call-from-psps` | T3 | psp mode: psps → VCF | keep |
| `inspect-psp` | T5 | Dumps what one psp stored over a BED, as a table | dev (§5.3) |
| `regenerate-census` | — | Rewrites the parameter-fitting sample stored in old psps' tails, for psps written by earlier builds | remove (§5.3) |
| `type-regions` | — | Writes how the reference is divided into repeat and non-repeat ground, as a table | dev (§5.3) |
| `estimate-contamination` | — | Contamination estimate on its own, from alignments, for comparing methods | dev (§5.3) |

The task names are already clear and consistent (`call-from-…`, `generate-…`, `estimate-…`),
so **I propose keeping all five main names.** Renaming would break every script and document
written so far for no gain a user would notice.

### 4.2 Options shared by several subcommands

Options that mean the same thing should be spelled the same, default the same, and be described
by the same sentence everywhere. These are the shared ones.

| Option | Used by | Default | What it controls | Proposed |
|---|---|---|---|---|
| `--reference` | all but `inspect-psp` | required | The FASTA every input was aligned to | keep |
| `--catalog` | all calling and fitting commands | `<reference>.repeats.parquet` | Where the repeat catalogue is | keep |
| `--alignment` (repeated) | direct mode, `generate-psps` | required | BAM/CRAM files; samples named by `@RG SM` | keep |
| `--psp` (repeated, file or directory) | `call-from-psps`, `estimate-parameters` | required | psp files | keep |
| `--output` / `--output-dir` | all | required | Where the result goes | change: `-` writes the VCF to stdout (§7.4) |
| `--log` | all user subcommands | `pvc_<output name>.log` beside the output | Where the log goes | new (§7.3) |
| `--tmp-dir` | subcommands that write scratch files | the output's directory | Where scratch files go | new (§7.5) |
| `--quiet` | all user subcommands | off | Terminal shows warnings and errors only | new (§7.1) |
| `--regions` | direct mode, `generate-psps` | whole genome | BED of the ground to call or store | keep |
| `--threads` | direct mode, `call-from-psps`, `repeat-catalog` | 0 = all cores; **1** in `repeat-catalog` | Parallelism; output never depends on it | change (§6) |
| `--force` | `generate-psps`, `estimate-parameters`, `repeat-catalog` | off | Overwrite an existing output | change (§6) |
| `--ploidy` | calling commands, `estimate-parameters` | 2 | Copies of each chromosome | keep |
| `--parameters` / `--defaults` | calling commands | one required | Where the model's numbers come from | §5.1 |
| `--build-index-if-missing` | direct mode, `generate-psps` | off | Writes a `.bai`/`.crai` if absent | advanced |

### 4.3 Calling options (`call-from-alignments`, `call-from-psps`)

| Option | Default | What it controls | Proposed |
|---|---|---|---|
| `--min-site-quality` | 1 | Records with a lower QUAL are not written; 0 writes everything | keep (move out of Advanced) |
| `--paralog-fdr` | 0.01 | False-discovery target of the hidden-dup filter, which removes calls better explained by a hidden duplication — two copies of a sequence collapsed onto one in the reference — than by a variant; 0 turns it off | change: rename `--hidden-dup-fdr`, advanced |
| `--paralog-filter-tag` | off | Keep those records, marked with a VCF filter, instead of dropping them | change: rename `--hidden-dup-tag`; the filter becomes `hiddenDup` (was `hiddenParalog`); advanced |
| `--max-candidate-alleles` | 6 | The most alleles a SNP/indel locus is called over. **Generic path only**: repeat loci use their own cap of 32, fixed in the code (`src/calling/allele_candidates/ssr.rs`) | change: rename `--max-generic-alleles`, and the help says repeat loci are not affected; advanced |
| `--max-cohort-locus-span` | 100 bp | The widest locus assembled, and so the longest deletion called | change: rename `--max-locus-span` |
| `--max-reads-per-position` | 1,000 | Reads of one read group kept at one position | dev; follows `--max-depth` unless set (§5.4) |
| `--max-read-group-depth` | 1,000 | A locus where any read group is deeper than this is not called | change: becomes `--max-depth`, which also sets the read cap (§5.4) |
| `--explain-loci BED TSV` | off | Writes why each locus in the BED was called as it was | dev (§5.3) |
| `--cohort-locus-builder-regions-len` | chosen from cohort size | How much genome one round of locus building covers; trades memory for speed | change: rename `--locus-builder-regions-len`; dev (§5.3) |
| `--psp-prefetch-bytes` (psps only) | 8 GiB | How much of the psps to read ahead into the OS file cache | dev (§5.3) |

### 4.4 What counts as a repeat

Six options decide which stretches of the reference go to the repeat (STR) caller rather than
the SNP/indel caller: `--min-copies` (default `8,6,6,6,5,4`, one value per motif length 1–6),
`--min-period` (1), `--max-period` (6), `--max-str-len` (100 bp), `--max-str-bundle-len`
(100 bp) and `--min-purity` (0.8). They appear on three subcommands today:

- `generate-psps` — **they are recorded in every psp**, so the cohort's files all carry them;
- `call-from-alignments` — direct mode has no psp, so it needs them here;
- `call-from-psps` — **where they must equal what the psps recorded, or the run is refused**.
  `estimate-parameters` already reads them from the psp headers instead.

**Agreed (owner, 2026-10-16):** keep them under their own help heading on `generate-psps` and
`call-from-alignments`, and **remove them from `call-from-psps`**, which reads them from the
psps as `estimate-parameters` does. A user can then not get them wrong at calling time.

#### 4.4.1 Sending more ground to the generic path: `--call-satellites` and `--generic-path-only`

**Asked for by a user; owner's rulings, 2026-10-16: two options, both off by default, both under
*Advanced*, and a run given both is refused.**

**Why a user wants them.** A run sorts the genome into four kinds of ground, and a default run
leaves one of them uncalled on purpose:

| Ground | What it is | Called by, in a default run |
|---|---|---|
| ordinary sequence | everything that is not a tandem repeat | generic path |
| STR tract | one tandem repeat, up to `--max-str-len` (100 bp) | STR path |
| STR bundle | repeats within 15 bp of each other, together up to `--max-str-bundle-len` (100 bp) | STR-bundle path |
| satellite | a tract or a bundle longer than those limits | nothing; the run report counts it as refused |

GATK calls on satellites, which is what the user noticed. The two options send more of this ground
to the generic path, the SNP/indel model:

| Setting | STR tracts | STR bundles | Satellites | Catalogue needed |
|---|---|---|---|---|
| default | STR path | STR-bundle path | not called | yes |
| `--call-satellites` | STR path | STR-bundle path | **generic path** | yes |
| `--generic-path-only` | **generic path** | **generic path** | **generic path** | no |

The owner expects the generic path's calls on repeats to be worse than the STR paths'. The options
exist for a user who prefers a call there to no call. `--call-satellites` changes only ground a
default run does not call at all; `--generic-path-only` gives up the STR model everywhere.

**Given both, a run is refused** before anything is read: `--generic-path-only` already sends
satellites to the generic path, so the pair has no meaning of its own. The message says so and
names the one to keep.

What the options mean for each subcommand:

- **`generate-psps` and `call-from-alignments`** take them. With `--generic-path-only`,
  `--catalog` and the six repeat settings of §4.4 are refused too, since no ground is sorted.
  `--call-satellites` keeps them: the catalogue and the length limits still decide what is a
  satellite.
- **The psp header records which of the three settings a psp was written with**, as it records
  the repeat settings, so `call-from-psps` and `estimate-parameters` follow the psps, and a cohort
  mixing psps written under different settings is refused.
- **`estimate-parameters`** with `--generic-path-only` psps fits no stutter, since there are no
  STR loci to fit it from, and writes a parameters file whose stutter section is empty (settled
  below; the ruling said without one); a calling run under that setting reads no stutter
  numbers. `--call-satellites` changes nothing in the fit: satellites never fed the stutter fit.
- **One limit stays under both.** A generic locus wider than `--max-locus-span` (100 bp) is not
  assembled. Inside a long satellite, reads that disagree about its length chain positions into
  loci wider than that, so **what these options add on satellites is mostly SNPs and short
  indels; the satellite's own length variation mostly stays uncalled.** The run report says how
  many loci were refused for width, and the user documentation says this plainly.

**Before `--call-satellites` could become the default**, its calls on satellites have to be
compared with GATK's on GIAB over that ground. Until then it stays off.

**Settled in implementation (2026-10-18,** plan [`call_satellites.md`](../impl_plan/call_satellites.md)**)**:
the three questions this section left open, a change to one of its rulings, and three things
the rulings did not say:

- **No catalogue at all, not an empty one.** Under `--generic-path-only` a run opens no catalogue:
  each analysed span is one generic region. An empty catalogue would have worked, but it is a file
  the user would have to build or be handed, and a psp would record a catalogue that played no
  part. A `--catalog` named for psps written with the option is refused by `call-from-psps` and
  `estimate-parameters`, which read none for them.
- **A satellite's edges are always ordinary sequence.** When the ground is typed, a satellite
  absorbs every tract and bundle member within 15 bp of it, so ordinary sequence, a contig end or
  a BED edge lies on each side, never an STR segment. Under `--call-satellites` each satellite is
  clipped to the analysed regions and merged with the ordinary sequence either side into one
  generic region; left as a region of its own, an indel across either old edge would be called by
  neither side. The STR edges further out are the junctions every run already has.
- **The parameters file is written with its stutter section empty, not without it** (owner,
  2026-10-18, changing the ruling above). The section is required by the file's format; an empty one
  is what a fit over zero repeat tracts already writes, every reader accepts it, and the file says in
  words that nothing was fitted. Leaving it out would have meant a new file format for a file one
  section shorter. The file records no `[repeat_routing]`, since no catalogue was asked anything.
- **`--call-satellites` keeps the fit unchanged by leaving satellite loci out of the census.** The
  census's positions never lay on satellites, but it also tallies every ordinary-sequence locus's
  read errors for the base-quality calibration, and with satellite loci in that tally each
  library's error multiplier moved in its seventh figure. Psps written with and without the option
  now fit to one file.
- **The psp format went back to 1.0** with the new header key, and every earlier psp is refused
  ([`psp_file_format.md`](psp_file_format.md), *Format 1.0 again*).
- **The run report** says how many satellite bases `--call-satellites` handed to the generic path,
  in both calling modes, beside the count of loci too wide to assemble.

### 4.5 `estimate-parameters`

| Option | Default | What it controls | Proposed |
|---|---|---|---|
| `--inbreeding` | fitted | States the inbreeding coefficient for every sample instead of fitting it | keep |
| `--skip-contamination` | off | Treats every read group as uncontaminated instead of estimating it | keep |
| `--str-param-estimates-at-once N` | 1 | How many groups of repeat tracts are fitted at once; memory against speed | change: rename `--str-strata-at-once`, advanced |
| `--threads` | **absent** | — | change: add it (agreed, owner, 2026-10-16). The fit already runs on every core through the global thread pool, so today it can only be narrowed with `RAYON_NUM_THREADS` |

### 4.6 `generate-psps`

`--force` (overwrite existing psps) and the shared options of §4.2 and §4.4. There is no
`--threads` option, and that was decided: samples are read one at a time, and a cohort is spread
by running one command per sample (`run_streaming.md` §5.2). **Proposed: keep it that way**, and
say in the help how to spread a cohort.

### 4.7 `repeat-catalog`

`--output`, `--threads` (default 1 — each thread holds a whole contig), `--force`, and six
settings for what the scan writes down: `--min-period`, `--max-period`, `--min-copies`
(`5,5,4,4,4,3`, deliberately looser than calling's so calling can choose within it),
`--min-flank-bp` (15), `--match-reward` (2), `--mismatch-penalty` (7).

**Agreed (owner, 2026-10-16):** keep `--output`, `--threads` and `--force`; put the six scan
settings under *Advanced*. A user should never need them, and a catalogue built with them changed is refused by
nothing, so they are a way to get quietly different results.

## 5. Decisions (all ruled, 2026-10-16)

### 5.1 Fitted parameters in direct mode

**Ruled (owner, 2026-10-16): left as it is.** A parameters file is made only by
`estimate-parameters`, from psp files. Direct mode calls with a parameters file or with
`--defaults`, and gets no route of its own from alignments to a parameters file.

**Why.** A fit needs many repeat tracts in every group of tracts it estimates stutter for, so a fit
over a few Mb of genome would not give good estimates, and the documentation will not suggest one.
The only other source of the numbers is a pass over the whole genome before calling, and a user
who is going to read the whole genome once to fit is better off writing psps in that pass, which
gives them the fit and psp-mode calling for the same cost.

**What follows for the documentation.** A parameters file is bound to its cohort: a calling run
refuses one that lists samples or read groups it does not have, or lacks ones it has
(`parameters_file.md` §6). A file fitted on a cohort therefore serves only that cohort, and a
cohort that has been fitted already has its psps. So **in practice direct mode calls with the
defaults**, and the user documentation says so plainly: direct mode is the quick route, with the
numbers compiled into `pvc`; psp mode is the route to calls made with fitted numbers. A user who
wants to change a default number can still start from the TOML block a `--defaults` run writes into
its log (§7.3) and pass the edited copy with `--parameters`.

**`--defaults` stays required** when no parameters file is given, so that calling with unfitted
numbers is a deliberate choice, and its refusal message names `estimate-parameters` as the way to
fitted ones.

### 5.2 One `call` subcommand, or two

`call-from-alignments` and `call-from-psps` could become one `pvc call` that takes either
`--alignment` or `--psp`. The calling options are the same in both. The ones that differ are
`--regions`, `--build-index-if-missing` and the repeat settings (alignments only) and
`--psp-prefetch-bytes` (psps only).

**Ruled (owner, 2026-10-16): two.** Each help page then lists only options that apply, and the
documentation can describe each mode on its own page. One subcommand would save a word of typing
and add a class of "this option does not apply to this input" errors.

### 5.3 What to do with the developer tools

`regenerate-census` repaired psps written by earlier builds. **Ruled (owner, 2026-10-16): it is
removed**, with its code. Before a first release no user has such psps, and a psp written by an
older build is rewritten with `generate-psps`; the refusal messages of `estimate-parameters` and
`call-from-psps` say so in place of naming `regenerate-census`.

`type-regions` and `estimate-contamination` are tools for developing and comparing the caller. `inspect-psp` and the calling commands' `--explain-loci` are the developer's tools for
tracing why a locus was or was not called (T5).

**Ruled (owner, 2026-10-15): the other three subcommands move under one `dev` subcommand.** A user
sees one line saying developer tools exist; a developer finds them all with `pvc dev --help`:

```
$ pvc --help
Commands:
  repeat-catalog         Build the reference's repeat catalogue (once per reference)
  call-from-alignments   Call samples from BAM/CRAM files and write a VCF
  generate-psps          Store each sample's evidence as a psp file
  estimate-parameters    Fit the model's parameters from psp files
  call-from-psps         Call samples from psp files and write a VCF
  dev                    Tools for developing and debugging the caller

$ pvc dev --help
Commands:
  inspect-psp             What one psp stored over a BED
  type-regions            How the reference is divided into repeat and non-repeat ground
  estimate-contamination  Contamination estimate on its own, for comparing methods
```

Each tool's name gets one word longer (`pvc dev inspect-psp`). The scripts under `scripts/` and
`benchmarks/` that call these tools change with it, as do the two scripts that call
`regenerate-census`.

**Options for developers cannot move under `dev`**, because they belong to a user subcommand.
They are listed under a *Developer* heading at the end of that subcommand's help:
`--explain-loci` on both calling commands, `--locus-builder-regions-len` on both, and
`--psp-prefetch-bytes` on `call-from-psps`. Hiding them (clap's `hide`) was rejected: an option
that `--help` never shows is findable only by someone who already knows it exists.

The user documentation does not describe the `dev` tools.

### 5.4 Two depth limits with the same default — one user option

Both limits exist for the same thing: a collapsed repeat, where many copies of a sequence are
mapped onto one copy in the reference and pile up thousands of reads (a median of 5,669 reads a
position in the 3× tomato samples). They act differently:

- `--max-reads-per-position` (the *read cap*) keeps at most 1,000 reads of one read group at a
  position and discards the rest. It is applied while alignments are read, so it bounds the
  memory and, in psp mode, what a psp stores. `generate-psps` has it.
- `--max-read-group-depth` (the *depth ceiling*) refuses to call a locus where any sample has a
  read group deeper than 1,000 reads, for every sample. It is applied when calling.
  `generate-psps` does not have it.

**At the defaults, the read cap never changes a call.** The ceiling counts the reads the cap
discarded (`src/run/depth_ceiling.rs`, `is_exceeded_by`). So a position the cap thinned holds
1,000 kept reads plus at least one discarded, which is over the ceiling, and the locus is dropped.
The cap's thinned reads are never called; its only effect is to keep the memory small on the way
to the drop. This is from the code; no run has checked it.

The two separate only when a user sets them apart. Someone with 2,000× data must raise both, or
the ceiling drops every locus the cap thinned. Someone who wants thinned calls instead of dropped
loci raises the ceiling alone.

**Options.**

- **(a) Keep both, under *Advanced*,** with help text saying how they interact.
- **(b) One user option, `--max-depth`** (default 1,000): a locus where a read group is deeper
  than this is not called, and reading keeps no more than this many reads of a read group at a
  position. Both limits follow the one number, so they cannot be set out of step. The separate
  read cap stays as a developer option for experiments.

**Ruled (owner, 2026-10-16): (b).** A user has one question — how deep is too deep to trust — and (b) asks it
once. Under (a), a user with 2,000× data who raises only the cap still has every deep locus
dropped, and nothing in the run says why. Thinned calls instead of dropped loci remain possible
under (b), through the developer option.

## 6. Changes that need no ruling

These make the interface consistent. I would make them in the implementation step unless the
owner objects.

- **`--threads` everywhere it can apply, meaning the same thing:** 0 = all cores. Add it to
  `estimate-parameters`. `repeat-catalog` keeps a default of 1, because each thread holds a whole
  contig in memory, and its help says so.
- **`--force` on every command that writes a file.** The calling commands overwrite an existing
  VCF without asking today, while `generate-psps` and `estimate-parameters` refuse. Proposed: all
  refuse unless `--force` is given.
- **One sentence per shared option.** Today `--catalog`'s help is three paragraphs on one
  command and one line on another. The longer explanations move to the user documentation; the
  help keeps one or two sentences and the default.
- **Fix the `estimate-parameters` description.** It says the inbreeding coefficient is
  *declared rather than fitted*; the `--inbreeding` option says it is fitted unless stated. The
  option is right.
- **Remove internal references from help text.** Help strings name spec files
  (`doc/devel/ng/spec/explain_loci.md`), the code's internal step numbers ("step 3's walk"), and
  dates of measurements. A user has none of these.
- **Define *psp* in the top-level help**, since every user meets the word.
- **A real top-level description.** `pvc --help` opens with "Top-level CLI for the `pvc`
  binary", which is a code comment leaking into the help. It should say what the caller does
  and show the two workflows of §3 in a few lines.
- **List subcommands in the order a user runs them**: `repeat-catalog`, `call-from-alignments`,
  `generate-psps`, `estimate-parameters`, `call-from-psps`. Today the list opens with
  `type-regions`, a developer tool, and puts `regenerate-census` between the two psp steps.

## 7. What a run tells the user: the terminal and the log

What a run prints is as much a part of the interface as its options. Today each subcommand prints
what its author needed while building it: a calling run prints progress lines on stderr and, at
the end, a report of per-sample counts on stdout; nothing is kept once the terminal closes.

**The owner's ruling (2026-10-15) sets the shape.** A run has two outputs for people:

- **the terminal** — short: which operation is running, on how many samples, how far it has got,
  and how it ended;
- **a log file** — complete: the input files by name, every setting and parameter needed to
  reproduce the run, and everything the run has to report.

**And a calling run no longer writes `calls.parameters.toml`.** Only `estimate-parameters`
writes a parameters file, because only it produces parameters. The parameters a calling run used
go into its log instead (§7.3). This reverses `parameters_file.md` §7, which made every run write
the file; that section is to be marked superseded when this spec is adopted.

### 7.1 Conventions both outputs follow

These are the common practice of command-line tools (the *Command Line Interface Guidelines*,
clig.dev, and the GNU coding standards), and what `samtools`, `bcftools` and GATK do:

- **Everything meant for a person goes to stderr. stdout is kept for data.** A tool whose stdout
  carries only its output can be put in a pipeline; one that mixes messages into it cannot.
  `inspect-psp` already writes its table to stdout, and the calling commands gain `--output -`
  for the VCF (§7.4). So the run report moves from stdout to stderr and the log.
- **Errors say what went wrong and what to do about it**, in one or two sentences, prefixed
  `error:`. Warnings are prefixed `warning:`. The exit status is 0 on success, 1 when the run
  failed, 2 when the command line was wrong (clap's own).
- **The terminal adapts to where it is going.** When stderr is a terminal, the progress line
  rewrites itself in place. When it is redirected to a file — a cluster job, `2> run.err` — each
  progress update is a new line, at a slower pace, with no control characters.
- **No colour.** It adds nothing a prefix does not, and breaks redirected output.
- **`--quiet`** on every subcommand: errors and warnings only. There is no `--verbose`; the
  detail lives in the log.

### 7.2 The terminal

A run says, in order:

1. **What it is doing, on how much.** One line naming the operation, the number of samples, the
   ground and the threads, and one line saying where the log is.
2. **Progress**, for anything longer than a few seconds: where it is in the genome, the
   percentage of the ground done, and the time elapsed.
3. **Warnings**, one line each. A warning that would repeat per sample or per locus is counted
   and printed once, with the detail in the log.
4. **How it ended**: what was written, how long it took, and the log's path again.

A direct-mode calling run, for example (the numbers are illustrative):

```
pvc call-from-alignments: calling 63 samples over the whole genome (782 Mb), 16 threads
log: pvc_calls.log
calling: SL4.0ch06:18,783,657  47%  1h12m
warning: 3 samples have read groups deeper than --max-read-group-depth at 1,204 loci; see the log
done: 4,068,400 records written to calls.vcf.gz in 2h31m
log: pvc_calls.log
```

The first line for each subcommand:

| Subcommand | First line |
|---|---|
| `repeat-catalog` | `pvc repeat-catalog: scanning genome.fa (12 contigs, 782 Mb), 1 thread` |
| `call-from-alignments` | `pvc call-from-alignments: calling 63 samples over the whole genome (782 Mb), 16 threads` |
| `generate-psps` | `pvc generate-psps: storing 2 samples' evidence over regions.bed (12.4 Mb), one after another` |
| `estimate-parameters` | `pvc estimate-parameters: fitting parameters from 63 psp files, 16 threads` |
| `call-from-psps` | `pvc call-from-psps: calling 63 samples from psp files over the whole genome (782 Mb), 16 threads` |

`generate-psps` also says which sample it is on (`sample 2 of 2: TS-104`), since it reads them
one after another.

### 7.3 The log file

**Written by every run of a user subcommand, always**, in the directory of the main output (the
working directory when the VCF goes to stdout). Its default name (owner, 2026-10-16):

| Subcommand | Log |
|---|---|
| `call-from-alignments`, `call-from-psps` | `pvc_calls.log` |
| `estimate-parameters` | `pvc_call_parameters.log` |
| `generate-psps` | `<sample>.psp.log`, one per sample, beside `<sample>.psp` |
| `repeat-catalog` | `<catalogue>.log`, beside the catalogue: `genome.fa.repeats.parquet.log` |

A file that has a log of its own takes the file's whole name plus `.log`, so the two are
plainly related and list together. `generate-psps` writes one log per sample because a cohort's
samples are usually written by separate commands running at once, and one shared log would
interleave them.

`--log <path>` puts the log elsewhere. Since the calling and fitting logs have fixed names, a
second run into the same directory finds the first one's log: it is refused or replaced together
with the output it belongs to (`--force`), and two runs that must share a directory each name
their own `--log`.

**It is written as the run goes**, a line at a time, so a run that crashes or is killed leaves a
log that shows how far it got.

What it holds, in this order:

1. **What ran.** `pvc` version and the full path of the binary (the version alone does not tell
   two development builds apart), the command line exactly as typed, the working directory, the
   start time with its time zone, the host name, the thread count. *(The git commit was first
   proposed here. Dropped in step 3: a build script that records it changes the crate's build
   inputs at every commit, and Cargo then recompiles the whole crate, minutes each time.)*
2. **The inputs, by name.** The reference, with its digest; the repeat catalogue, with its
   digest; each alignment or psp file with the sample it holds and its read groups; the regions
   BED, with the total ground it names.
3. **Every setting's value, including the ones left at their default.** A default can change
   between versions; a log that records only what was typed cannot reproduce a run made with an
   older one.
4. **The parameters the run scored with** (calling runs): where they came from — a file, named,
   with its digest, or the defaults compiled into this version — and **the full values, as a
   TOML block in the parameters-file format**. Copied out of the log, that block is a valid
   `--parameters` file, so the one thing the dropped `calls.parameters.toml` was good for, reusing
   or editing a run's numbers, is kept.
5. **Events, each with a timestamp**: stages started and finished, progress, every warning in
   full.
6. **The report**: the per-sample counts the terminal no longer prints — loci dropped as too deep,
   loci too long to assemble, what the hidden-dup filter removed.
7. **How it ended**: the finish time, elapsed time, peak memory, and the exit status; on failure,
   the full error.

**Plain text, for a person to read.** Header and settings as `name = value` lines, events as
`<timestamp> <level> <message>`. A machine-readable format such as one JSON object per line was
considered and rejected: the reader of a log is a person working out what a run did, and the
parts a program needs again — the parameters — are already in a format a program reads.

### 7.4 The VCF on stdout

`--output -` on both calling commands writes the VCF to stdout, so a run can feed `bgzip`,
`bcftools` or anything else directly:

```
pvc call-from-psps --reference genome.fa --psp psps/ --parameters cohort.parameters.toml \
    --output - | bcftools view -i 'QUAL>20' -Oz -o filtered.vcf.gz
```

- **What goes to stdout is uncompressed VCF**, as `bcftools` writes by default. Compression is
  the next tool's job; a user who wants `.vcf.gz` without a pipe names a file.
- **`--output` stays required.** `-` must be typed, so a run that forgets `--output` is refused
  instead of pouring a genome's VCF into the terminal.
- **The file-output guarantee does not carry over.** A named VCF appears whole or not at all,
  because it is written as `<name>.tmp` and renamed at the end. A stream cannot be taken back: a
  run that fails halfway has already sent half a VCF downstream. The exit status is then 1, and
  the log says where it stopped. This is the usual contract of a tool writing to a pipe.
- **With the hidden-dup filter on, nothing reaches stdout until the whole genome is
  called**, since no record's fate is known before then (the records wait in the temporary
  directory, §7.5). A downstream tool simply waits.

### 7.5 Temporary files

A run writes two kinds of files that are not its final output, and they must go to different
places:

- **Scratch files** — data a run parks and reads back, such as the records the hidden-dup
  filter holds until every record has been scored (several times a compressed VCF's size). These
  go to **the temporary directory**.
- **An output in the making** — the VCF, psp, catalogue or parameters file written as
  `<name>.tmp` and renamed into place when complete, so that a reader never sees half a file.
  **These stay beside the output.** A rename is only instant and all-or-nothing within one
  filesystem; moved to a temporary directory on another disk, the final step becomes a copy, and
  the guarantee is gone.

**`--tmp-dir <dir>` on every subcommand that writes scratch files.** By default it is the
directory the output goes to (the working directory when the VCF goes to stdout). Inside it the
run creates a directory of its own, `pvc_tmp_<random>/`, holds every scratch file there, and
removes it when the run ends, whether it succeeded or failed. The user's directory itself is never
removed. A directory of its own per run is what lets several runs share one `--tmp-dir`, a
cluster's scratch disk for example, without colliding. The log records where the temporary
directory was and how much the run put in it at most, so a user can size the scratch space for the
next run.

### 7.6 The VCF header

The VCF already records `##source` and `##commandline`. **It also names its log and the
parameters' source** (`##pvcLog=pvc_calls.log`, `##pvcParameters=defaults` or the file and its
digest), so a VCF separated from its log still says how its numbers were chosen. `##source` reads
`pvc <version>` (it read `ng <version>`). It carries no commit, for §7.3's reason, and because a
header that changed at every commit would move the cross-platform checksums and the identity
oracle's at every commit too.

### 7.7 Logging in the code

Messages are printed today by about sixty `println!` and `eprintln!` calls spread across the
subcommands and the library. **They are replaced by a logging framework: the `tracing` crate**,
the standard one in Rust, with `tracing-subscriber` to send its messages on.

**How it works.** Code anywhere in the crate states an event and its importance —
`info!("calling {n} samples …")`, `debug!(…)`, `warn!(…)` — and says nothing about where it goes.
The binary, at start-up, installs two destinations, each with its own threshold of importance:

| Level | Used for | Terminal | Log file |
|---|---|---|---|
| `error` | the run cannot continue | yes | yes |
| `warn` | something the user should check | yes | yes |
| `info` | what the terminal shows (§7.2): the operation, progress, the end | yes | yes |
| `debug` | inputs, settings, parameters, the report, stage timings | no | yes |
| `trace` | developer detail | no | only when asked for |

`--quiet` raises the terminal's threshold to `warn`. A developer can turn on `trace` with the
`RUST_LOG` environment variable, which `tracing-subscriber` reads; no option is added for it.

**Rules that keep this working:**

- **The library never prints.** Only the binary decides where messages go. A `println!` in
  `src/` outside `src/main.rs` and `src/cli/` is a defect, and a test can check for it.
- **One event, one message.** The terminal line and the log line come from the same event; there
  is never a separate "log version" of a message to keep in step.
- **A progress update is an ordinary event with its own target**, `pvc::progress`. The terminal
  destination draws it in place when stderr is a terminal, each update over the last, and lets
  the next ordinary message leave the last update standing and start below it; redirected, each
  update is a line. The log records each update as a line. *(Built this way in step 1 rather than
  with `indicatif`, as first proposed: drawing one line in place is a few lines of the terminal
  destination, and `indicatif` would add a dependency and the coordination with `tracing` it
  needs.)*

**Why `tracing` and not the older `log` crate.** Both would do the two-destination arrangement.
`tracing` also has *spans*: a stage such as "fit contamination" is opened and closed once, and
every message inside it carries the stage's name and its timing. That is what the log's stage
timings (§7.3, item 5) need, and with `log` they would be written by hand.

Building this is the first step of the implementation, before any message is reworded.

## 8. Out of scope

- The VCF's content and header beyond §7.6, which `vcf_output.md` owns.
- The parameters file's format, which `parameters_file.md` owns.
