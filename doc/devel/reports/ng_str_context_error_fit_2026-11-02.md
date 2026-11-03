# The misread-after-a-repeat table, fitted per read group — built and measured, 2026-11-02

*Revised 2026-11-03 for the owner's ruling of §7. Branch `str-context-error`, from main 92b00dd8. Designs:
[`../ng/spec/str_context_error.md`](../ng/spec/str_context_error.md) (the term) and
[`../ng/spec/str_context_error_fit.md`](../ng/spec/str_context_error_fit.md) (the fitted table),
both approved by the owner on 2026-11-02. The prototype and the spec are committed (8d0c745c,
b3734f26); the fitted table is not.*

## 1. The answer

**The fitted table works and calls as well as the prototype's hand-made one.** On HG002 at 30× this build
removes 18 to 21 false calls inside STR segments on chromosome 20 and 26 to 34 on the tandem-repeat
benchmark; at 5× all three builds are within a handful of calls of each other; the two simulations
do not move.

**Charging after repeats of 2 to 4 bases, which the owner asked to try, buys nothing.** The same fit
with every one of those entries set to zero calls within 6 calls of it in all four HG002 runs, the
misread bases there carry good base quality — what stutter written as mismatches looks like, not a
sequencing misread — and the charge makes most STR loci score each read twice. **On 2026-11-03 the
owner ruled the table covers runs of one base only, neither counting nor charging other repeats**
(§7); the branch now does that, and its fitted route is the *fitted, runs only* column below.

HG002 against GIAB, vcfeval, genotype match: missed / false calls. *Prototype* is the hard-coded
table of runs of one base; *fitted table* is this build; *fitted, runs only* is this build's own fit
with every entry after a repeat of 2 to 4 bases set to zero.

| run | main | prototype | fitted table | fitted, runs only |
|---|---|---|---|---|
| chr20, 30×, `--defaults` | 753 / 473 | 746 / 453 | 747 / 452 | — |
| chr20, 30×, fitted | 841 / 373 | 837 / 357 | 837 / 355 | 838 / 356 |
| chr20, 5×, `--defaults` | 24,519 / 4,619 | 24,519 / 4,619 | 24,520 / 4,620 | — |
| chr20, 5×, fitted | 24,913 / 4,346 | 24,908 / 4,341 | 24,908 / 4,341 | 24,908 / 4,342 |
| tandem-repeat benchmark, 30×, `--defaults` | 14,245 / 2,185 | 14,212 / 2,151 | 14,212 / 2,151 | — |
| tandem-repeat benchmark, 30×, fitted | 14,512 / 1,919 | 14,478 / 1,890 | 14,487 / 1,893 | 14,485 / 1,895 |
| tandem-repeat benchmark, 5×, `--defaults` | 24,028 / 2,537 | 24,030 / 2,540 | 24,033 / 2,542 | — |
| tandem-repeat benchmark, 5×, fitted | 24,300 / 2,389 | 24,306 / 2,392 | 24,309 / 2,392 | 24,303 / 2,391 |

Inside STR segments on chromosome 20 at 30× (missed / false): main 115 / 139 with `--defaults` and
108 / 118 fitted; this build 109 / 118 and 104 / 100.

## 2. What was built

**The count, in the parameter pre-pass** (`src/parameter_estimation/misread_after_runs.rs`). At
every generic position the walk hands the census, the reference around it decides which entry of the
table a forward read and a reverse read there belong to. Where one direction crossed a repeat and
the other did not, the first direction's reads are counted as *through* reads and the second's as
*arriving first*, with how many of each showed the repeat continuing or another wrong base, per read
group. The reference is read once per 64 kb stretch, so the count is linear in the ground walked. The
totals travel in the census, whose format goes from 4 to 5: **every existing psp has to be
regenerated** before it can be fitted from, and the freshness check says so.

**The rule deciding a base's entry** (`src/calling/likelihood/run_context.rs`, `contexts_along`) is
the same function for counting, on the reference, and for charging, on the candidate allele. As
built and measured on 2026-11-02 it covered repeats of period 1 to 4, 43 entries; since the ruling
of §7 it covers runs of one base only, 19 entries.

**The fit** (`estimate-parameters`): per read group and entry, the through rate minus the
arriving-first rate, zero where that is under two standard errors, and the default where either group
has fewer than 1,000 reads. It logs how many entries it fitted and checks each next-base entry
against its misread bases' base qualities.

**The parameters file** carries the table as `[misread_after_runs]`, one row per read group and
entry, each excess a warranted number; format 7 → 8, and a version-7 file reads with every read
group on the defaults. **Calling** hands each read group's table to the STR scoring.

**The defaults** are now this pre-pass's own fit on HG002 chromosome 20 at 30×, regenerated from its
parameters file, replacing the prototype's hand numbers.

Tests: all 5,081 library tests and every integration test pass (debug build), with new ones for the rule, the count, the census round trip, the
fit, the check, the file and the scoring under a read group's own table; each fails without its
piece. The three cross-platform checksum tests moved and were re-recorded **after checking that the
VCFs are byte-identical except the header line naming the parameters file's checksum** — only the
parameters file moved. They were re-recorded on arm64 Linux (the dev container); the macOS half of
"the same bytes on every platform" has not been re-checked. `cargo clippy --all-targets`: no warning
in the changed code.

## 3. The fitted table against the prototype's

Chromosome 20, 30×, the next base after a run of one base read as one more base of the run, per
10,000 reads (± standard error):

| run length | fitted, 30× | fitted, 5× | prototype |
|---|---|---|---|
| 8 | 21 ± 2 | 19 ± 5 | 10 |
| 12 | 205 ± 7 | 193 ± 17 | 193 |
| 16 | 639 ± 14 | 655 ± 34 | 902 |
| 20 | 1,051 ± 26 | 1,143 ± 67 | 1,702 |
| 23 or more | 1,174 ± 18 | 1,186 ± 44 | 2,036 |

**The 5× fit agrees with the 30× fit within its error at every entry**, so 64 Mb at 5× is already
enough; a whole genome at 3× has about ten times as many reads.

**Above 15 bases the fitted rates are 30 to 40 in 100 lower than the prototype's**, against the 11 to
20 the spec predicted (§3.4) from forward reads that started inside the run and so did not cross it.
The prototype counted only reads that crossed the whole run. In nats a read the difference is
`ln(902/639)` = 0.34 at 16 bases and 0.55 above 22, against the 3 to 6 the term moves; the calls in
§1 show it costs nothing measurable.

After repeats of 2 to 4 bases the fit agrees with the spec's measurement (§2.1): the only next-base
entries above 50 in 10,000 are 2-base repeats of 20 bases or more (82 ± 13 at 20 to 29 bases, 214 ±
14 at 30 or more, against 124 and 233 there). Every further-on entry of periods 2 to 4 is below 70 in
10,000, and most below 25 — but on millions of reads, so most clear two standard errors.

**The base-quality check passes for runs of one base and flags 2-base repeats.** After 2-base repeats
of 16 to 29 bases the bases read as the repeat continuing claim an error rate of 2 to 5 in 1,000
(about quality 25), while the reads through the repeat show them 13 and 21 times in 1,000. Bases the
sequencer is sure of, that differ by direction: a read one repeat unit longer, written by the mapper
as mismatches, which the stutter model already prices.

## 4. The simulations

Neither simulation plants a misread after a repeat, so nothing should move.

- **The flank simulation**: identical scores at 1 sample 3× and 30× and 20 samples 3×; at 20 samples
  30×, 2 of 4,810 records lose a rare candidate allele after a 3- or 4-base repeat, no genotype
  changes, and the score gains one right haplotype (876 → 877).
- **The stuttering simulation**: identical in every route at one sample, 30×; in the 50-sample 3×
  cohort, of about 150,000 sample-segment cells, at most 8 move in any route (false calls 7,388 →
  7,395 with `--defaults`, 387 → 387 fitted).
- **The fit on simulated reads**: at one sample 30×, 1 of the 27 entries with enough reads charges
  anything; in the cohort, 16 of 655, 2.4 in 100 — what a two-standard-error rule passes by chance
  where nothing is there (2.3 in 100). The cohort's 1,495 entries short of reads, 50 read groups on a
  4 Mb genome at 3×, fall back to the defaults.

## 5. The range

- **One sample, few reads**: the table is fitted per read group over the whole run, and §3 shows a
  5× chromosome fits it as well as 30×. A small region or genome falls back to the defaults entry by
  entry, marked `defaulted`.
- **Thousands of samples**: 19 entries of seven totals per read group in the census (1.1 kB) and in
  the parameters file (19 rows); the fit is per read group and independent.
- **At 5× the term costs a few heterozygous calls** inside long runs with one reference read
  (`str_context_error.md` §6.3); fitting does not change that.
- **Not measured**: a real cohort at 3× (tomato), and timing under controlled conditions.

## 6. Where the build departs from the spec as approved

All three are written into `str_context_error_fit.md`.

1. **A run of one base won over a longer repeat of another period** in the rule of §2.3. Without
   it, a 19-base `TTTC` repeat six bases back took the entry from an 8-base `T` run right before a
   base, and the base lost the run's charge; about 280 records on chromosome 20 moved against the
   prototype. With it the rule reproduced the prototype exactly on chromosome 20 at 30×, except one
   record beside a gap of `N`s, which the prototype counted as a run. *Moot since §7's ruling: with
   runs only, the rule is the prototype's plus the `N` guard.*
2. **An excess under two standard errors is charged as zero** (§4.1): 6 ± 27 in 10,000 is no
   evidence of anything, and charging it costs a second scoring.
3. **The base-quality check compares with the rate the through reads show**, not with the
   substitution rate plus the excess (§4.2): the first build warned at every short repeat, because
   the misread bases include the background misreads both directions share.

## 7. The decision: which periods to charge

**Decided by the owner on 2026-11-03: runs of one base only, and the other repeats neither charged
nor counted** — they showed no improvement, and biologically there is no reason they should: the
misread comes from the sequencer losing phase over many identical bases in a row. The branch was
changed to that the same day (19 entries; the census and parameters-file layouts lost the period;
the section is now `[misread_after_runs]`), with every test passing and the cross-platform checksums
re-recorded after checking the calls are byte-identical.

**The runs-only build on HG002 chromosome 20 at 30×** (missed / false; inside STR segments in
brackets), regenerating its psps and fitting from them:

| route | the build with periods 2 to 4 | fitted, runs only (§1) | runs-only build |
|---|---|---|---|
| `--defaults` | 747 / 452 (109 / 118) | — | 747 / 453 (109 / 119) |
| fitted | 837 / 355 (104 / 100) | 838 / 356 (105 / 101) | 838 / 356 (105 / 101) |

**Its counts are not exactly the earlier build's**: a position where a run stood on one side and a
2-to-4-base repeat on the other had an entry in both directions and was skipped (the spec's §3.2:
a position is counted only where one direction has an entry); it is now counted. Each entry rests on
slightly more reads — 37,030 through reads after 16-base runs, against 36,753 — and its excess
moves by much less than its standard error (0.0642 against 0.0639, ± 0.0013). Every other fitted
parameter is identical. The calls differ from the *fitted, runs only* column in genotype qualities
and allele frequencies at 310 of 102,580 records, in the alternative allele at one (chr20:401,686,
where the alternative allele replaces the `C` between two runs of `A` with an `A`, and the
runs-only build adds one more `A` to it), and in none of the
scores. Under `--defaults` the built-in table no longer charges after repeats of 2 to 4
bases, which moves one false call.

What I had recommended, and the evidence it rested on, kept for the record:

**Charging only after runs of one base, and keeping the count and the parameters-file rows for
periods 2 to 4 with their charge set to zero.**

- **The gain is nil.** On HG002 the fitted table with and without its period-2-to-4 entries differs
  by at most 6 calls in any run, in both directions (§1, last two columns).
- **The cost is a second scoring at most STR loci.** A candidate is scored once per read direction
  wherever a charged repeat sits in or beside it, and every common 2-base repeat of 12 bases or more
  has a charged further-on entry. The fitted chromosome-20 call at 30× took 18 s against 12 s with
  those entries zero; another session was compiling at the time, so that is a hint, not a
  measurement.
- **What the data say about those bases**: the check of §3 — high-quality bases that differ by
  direction, the signature of stutter, not of a misread.
- **Keeping the count** costs nothing at calling, and a library where a period-2 effect is real would
  show it in its parameters file.

The next steps are a timing measurement at six threads, then a tomato run to check a 3× cohort with
real misreads.
