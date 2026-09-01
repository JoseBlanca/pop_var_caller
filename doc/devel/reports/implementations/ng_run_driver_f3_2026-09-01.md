# F3 — the run report: what the run refused, and why

**Date:** 2026-09-01. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone F step F3. **Design:** [`cohort_merge.md`](../../ng/spec/cohort_merge.md) §3.3 (the
refusal counts reach the user), [`parameters_file.md`](../../ng/spec/parameters_file.md) §8 (*"the
file must not hide which"* numbers have an honest default), with §2's warrant ladder underneath
it, [`run_streaming.md`](../../ng/spec/run_streaming.md) §8 (the finish-time
read-filter tally).

---

## Why a run has to say anything at all

**A VCF cannot distinguish ground the caller examined and found nothing at from ground it never
spoke for.** That is `cohort_merge.md` §3.3's argument for the failed-locus count being *counted*
and not merely dropped — *"in the VCF alone, that absence is indistinguishable from analysed and
found nothing"* — and it holds for every other kind of nothing this caller produces: a repeat
tract it has no generator for, a locus where the allele cap left nobody callable, a sample whose
reads were all filtered away. Each looks, in the file, exactly like a quiet stretch of genome.

F1's implementation report records what that costs. A run over a 60-base `AT` tract at 24 reads a sample
**printed six zeros and exited successfully** — indistinguishable from a clean genome. What
follows is the run saying which.

## What a run now prints

This is the fixture run's own output, verbatim — two samples with no reads over 300 bases of the
shared fixture reference:

```
calls: /tmp/…/calls.vcf
parameters: /tmp/…/calls.parameters.toml
records written: 0
loci called: 0 — 0 written, 0 establishing no variant and so left out
analysed ground: 2 intervals, chr1:1-100 … chr2:1-200 — 300 bases asked for, in 2 typed regions
  called: 300 bases (100.0%)
  not called — repeat tracts this caller has not built yet: 0 bases (0.0%)
  not called — tandem arrays longer than this run types as callable: 0 bases (0.0%)
loci the merge declined to assemble for being too wide: 0
loci where the allele cap left no sample callable: 0
samples: 2 — 0 whose reads the caller used, 0 whose reads the filters took all of, 2 that contributed none
  no read reached the caller: zeta, alpha — either the sample has none over this ground, or none of its ground was walked; written ./.
numbers behind the calls: 0 of 7 groups the file says were fitted
  taken from constants or supplied, not measured here: the base-quality calibration, contamination, the inbreeding coefficients, the ordinary-site prior's seed, repeat-tract slippage, repeat-tract length spectra, repeat-tract substitution rates
```

A run that refused something adds, under its count, up to five spans with their lengths and the
bound they were measured against:

```
loci the merge declined to assemble for being too wide: 8
  chr1:100-160 (61 bases)
  …
  … and 3 more
  the bound is --max-cohort-locus-span 50; raise it and call again if the lengths above cluster just past it
```

Five rules shape it, and each was a choice:

- **The ground is named, not only counted.** Two runs in review printed *"analysed ground: 3000
  bases"* with nothing on either page saying which chromosome — and the report's whole argument is
  that a VCF cannot say what ground the caller spoke for.
- **Every count is a share of a stated whole, and the whole is what the walk covered.** A repeat
  tract is typed and walked whole even where a BED asks for part of it (spec §4.2: findings whole,
  generic clipped), so the bases handed to the walk can exceed the bases asked for. Dividing by
  the latter printed **200.0%** on a BED of 120 bases inside two tracts. The shares are of the
  three parts' own sum, and both totals are printed together whenever they differ.
- **A refusal that did not happen gets a count and no advice.** A line telling somebody how to fix
  a thing that did not occur is a line they read and discard. Where a refusal *did* happen it
  shows a handful of spans **with their lengths** and the bound in force **as the flag a person
  types** — which is what §3.3 says a non-zero count should lead a reader to, and none of the
  three could be done from the page before.
- **A reason that did not fire is not printed.** A read filter has nine reasons and a run trips
  two or three; the other six as zeros are six numbers to scan past.
- **The groups of numbers that were *not* fitted are named, not counted.** A run whose
  contamination and slippage are compiled-in constants is a different claim from one whose
  base-quality calibration is, and *five of seven* says neither.

## The three counters the report needed, and why none existed

**Three counters, and each was recorded as owed to this step before it began** — the read-filter
tallies in `arch/run_streaming.md` §3.4 and in `CohortWalkTallies`'s own documentation, the other
two in the milestone's handover. The report cannot state its arithmetic without them.

### The per-read-group read-filter tallies

Recorded as needing *"a change to the generator, not an accessor"*, and that was right — but the
change is one the generator already makes for the *aggregate* cursor counts, one axis finer.
Read-filter tallies belong to a cursor from the moment it is made
([`run_streaming.md`](../../ng/spec/run_streaming.md) §8), and a cursor is **rebuilt at every
chromosome change**, so a walk had already lost every contig but its last: a run over twelve
chromosomes would have reported the twelfth's drop rates as the whole run's. The generator now
takes a retiring cursor's per-read-group counts at the boundary — beside the line that already
took its aggregate ones — and sums the live cursor in when asked.

Reached through a **defaulted `LocusGenerator` trait method**, so a boxed generator can be asked;
that is the same gap `LocusGenerator::counts` closed for the ten locus counters, whose own note
records that *"a walk that emitted nothing for a covered region counted the truncations that
explained it into a struct nobody could see"*. Kept off `PileupGeneratorCounts`, because that
struct is printed verbatim by the dump tools against byte-identical committed baselines, and this
is a fact about the *reader* rather than about the loci.

Spec §8 names the failure if the sum is skipped, and it applies one contig at a time exactly as it
does one worker at a time: drop rates under-report **"silently, since every number stays
plausible"**.

### `LocusCounts::regions_handled_bp`

The ground partitioned in *regions* — `regions_handled` plus the two unhandled counters sum to
`regions_in` — and in bases only on the unhandled side. So a run could say it handled 9,000 of
10,000 typed regions and not how much genome that was, and **typed regions differ in length by
orders of magnitude**: half the regions can be a twentieth of the ground. One counter, incremented
where `regions_handled` already is, completes the partition.

### Contigs named rather than numbered

`GenomeRegion`'s `Display` writes `contig 0:15-15`, and its own documentation says why: a region
is a position with no reference beside it, so it cannot spend a contig table it does not have. A
**run** has one. `RunReport` therefore names every span it shows, and the `Display` is unchanged —
it still has nothing to name a contig with, and the fix belongs where the table is.

**The fourth is not built**, and is recorded rather than forgotten: `RunError::RecordNotWritten`
still renders its locus through that same `Display`, so the one message whose job is to say how
far a partial file got prints `contig 0:15-15`, whose chromosome a reader cannot match
against the VCF they are holding. It is an error type, reached from a path that has no contig
table in hand; giving it one is a change to `RunError`'s shape.

## What was measured

In the development container, at this step's working tree.

| check | before F3 | after F3 |
|---|---|---|
| `cargo test --lib` | 5,868 passed, 13 ignored | **5,885 passed**, 14 ignored |
| `cargo fmt --check` | exit 0 | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 | exit 0 |
| `cargo doc --no-deps` | 26 unresolved-link errors, 23 redundant-link warnings | unchanged |

**17 tests added**: 15 in `ng::run::report` and 2 in the pileup generator, plus one ignored test
that prints the report for a person to read (`the_report_a_person_sees`).

**The report is lines, not printed output, and that is why it is testable.** F2's implementation report records its correctness pass finding the summary the one part of this command a mutation could change with the whole
suite green — swapping its two numbers, or printing the VCF's path on the parameters line, both
survived. Every rule above is now an assertion.

**Two things reading the printed report changed.** Two different lines opened on the word
`parameters` — one the file's path, one the count of fitted groups — so a reader scanning for
either could land on the other; the second is now `numbers behind the calls`. And a sample with no
reads was given a filter line reading `0 reads kept, 0 dropped`, which the line naming it as
having no reads already said.

## What the reviews found

Two agents reviewed the uncommitted step in isolated worktrees. **The artefact pass drove nine
runs on genuinely different inputs** — ordinary ground, a contig that is almost all repeat tract,
a sample every read of which is a duplicate, a sample with no reads, a BED cutting across typed
regions, and a 120-base deletion that made the too-wide refusal fire.

**The one that mattered: the report told a false thing about the VCF.** A sample the caller could
not use was printed as *"each still carries a genotype, from the prior alone"* — which is what the
**loop** does and not what the **file** writes. `vcf_output.md` §7.1 no-calls a sample whose
likelihoods are flat, F1 implemented it, and every such sample is `./.` in the VCF. It was
reassurance that stopped a reader opening the file, and the file said the opposite.

**And "no reads here" collapsed three problems and was false for two of them.** A sample whose
720 reads were all duplicates was called *no reads* four lines above the line saying it had 720;
so was a sample whose reads exist over the ground but whose ground no generator walked. A
geneticist checks the duplicate marking in one case and the sample sheet in the other, so the
three are now told apart.

**The advice named a flag that did not exist.** `--max-cohort-locus-span` was not an option of
this command; typing what the report told you to type got a clap error. Both bounds are flags now
(`--max-cohort-locus-span`, `--max-candidate-alleles`), each refused when it is not a bound, and
the advice quotes the value in force beside each refused span's length.

Two smaller ones, both fixed: read groups were **numbered where contigs were named** — now
`library rg1 of zeta` — and *"satellite, which this caller will never call"* claimed a permanent
refusal from a threshold that is a tunable 100 bases, over 83.4% of one reviewed contig.

**The correctness pass ran 14 mutations and killed 4.** Ten survived, and five of them were the
same hole: **the step's hardest change had no test at all.** Deleting the harvest at retirement,
deleting the live sum, double-counting the live cursor, dropping the filled slot from the set's
sum, and replacing the whole call with an empty vector each passed all 5,880 tests. Deleting the
live sum alone means **every single-contig run reports zero drops**, because the last cursor never
retires.

Two tests close it, beside `cursor_tallies_are_taken_from_a_chromosome_before_it_is_retired`,
which is the same claim one axis coarser. One walks two contigs with a plain read and a
duplicate-flagged read on each and asserts `(1, 1)` after the first and `(2, 2)` after the second —
**measured against the mutation**: deleting the retirement harvest makes it read `(1, 1)` where
`(2, 2)` is expected. The other is the cross-check the review asked for, and it holds: the reads
the filters kept, summed per read group, equal the reads the cursor decoded.

**And one meaning defect the pass found by reading.** `other_sample` was counted as a dropped
read, against its own field's documentation: *"Not a drop, and deliberately outside every other
counter here … counting it as a drop would make a shared file look like a low-quality one."* On a
cohort sharing one multi-sample BAM it would have dominated every sample's drop count. It is
reported on its own terms now, and both the count and the reasons come from **one exhaustive
destructuring**, so a counter added to `ReadFilterCounts` later must be routed here or the build
stops — which is `ReadFilterCounts::add`'s own guarantee, and which naming the fields one by one
had given up.

**Recorded and not acted on:**

- **`regions_handled_bp` and `analysed_bases` can each be mutated with the suite green.** Every
  report test builds `LocusCounts` by hand, and the one end-to-end run that would show them is the
  ignored printing test. Closing it means an end-to-end assertion over a fixture whose ground is
  not all one kind.
- **The repeat-tract generator has a reader of its own and does not override
  `read_filter_counts`.** `SsrGenerator` keeps its own `retired_cursor_counts` and retires at a
  boundary exactly as the pileup generator does, so the trait default's justification — *"a
  generator with no reader of its own"* — does not describe it. No run is affected: both tract
  slots are unfilled and `generic_path_generators` refuses them. **The moment one is filled the
  report under-reports by whatever share of the ground is tracts, silently.**
- **A contig the cursor read nothing on contributes a fabricated `None` read group**, which
  `AlignmentCursor::read_group_counts`'s own comment calls "a wart worth knowing about". Invisible
  today, because an all-zero entry names no reason; it becomes visible on a shared BAM, where the
  fabricated entry carries a non-zero `other_sample`.
- **`reset_read_group_counts` would desynchronise the two harvests**, since it resets the
  per-group tally and not the aggregate. It has no production caller today; the cross-check test
  above would catch it.

## What F3 does not do

- **It is not a log and holds no spans beyond a handful.** Where the refused loci themselves
  surface — a sidecar, a BED — is `cohort_merge.md` §14's open question 5. This states the counts
  and the ground they cover, which is what §3.3 fixes as reaching the user.
- **The assembly check is not in it.** `CohortWalkTallies::assembly_check` says whether every
  sample's contig checksums were compared and agreed, or whether none could be; a run report
  should say which, and where it belongs is beside the identity lines a run does not yet print
  (no `##fileDate`, no caller version).
- **Nothing is printed while a run is going.** Single-threaded over a cohort of CRAMs that is a
  long silence with no way to tell a slow run from a hung one — F1's report records it, still
  open, and a progress line is not a report.
