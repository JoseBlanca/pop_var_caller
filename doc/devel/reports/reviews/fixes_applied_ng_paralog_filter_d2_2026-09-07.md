# Fixes applied — ng_paralog_filter_d2

**Date:** 2026-09-07
**Review:** [ng_paralog_filter_d2_2026-09-07.md](ng_paralog_filter_d2_2026-09-07.md) (0 Blocker /
11 Major / 13 Minor, Request-changes; two sub-agents in isolated worktrees)
**Branch:** `ng-paralog-filter`

## What was done

**All eleven Majors and all thirteen Minors applied. Nothing deferred.** The five surviving
mutations are re-run and closed — four killed by strengthened tests, one closed by the type system.

## The code

### The clock said one thing and measured another (M-code-1, M-code-2)

Two descriptions of the timing span, both wrong, and the one an operator reads was the wronger.

`SpillFile::parking_began` is renamed **`named_at`**, because that is when it is taken: `beside`
creates nothing, so no record has been parked yet, and `beside` is also used purely as a path helper
in two test modules where "parking began" means nothing. The field's doc now says the clock starts
when the spill is *named*, just before the calling loop, and not when its first record arrives.

`WhereTheTimeWent::calling`'s doc now says what is inside the span — the calling loop and the
spill's own flush — and lists what is already done before it starts: the reference read, the
segmentation built, every input opened and index-checked, the run parameters fitted or read, and the
header's metadata assembled. The printed line changed to match:

```text
time from the start of the calling pass, 988ms in all: calling 934ms (95%), … — the run's startup
before the calling pass, which reads the reference and the catalogue and opens every input, is not
counted here
```

**This matters more than a doc usually does**, because the step's conclusion is a *share*: a reader
who believed the old doc would put the file-opening cost inside `calling` and read the scoring share
as smaller than it is.

### The header assembly was charged to nobody (Mi3)

`began_writing` moved above the header work, so assembling the header — which copies one `String` per
sample, a thousand allocations on a thousand-sample cohort — is inside `writing`, where it belongs:
it is the header being written. The line says "in all", and now the four figures are.

### Units that could not carry a small run's numbers (Mi1, Mi2)

Two helpers, `how_long` and `how_big`. Before them the line read `0.01s in all: calling 0.00s (43%),
… 0.00s (31%), … 0.00s (5%), … 0.00s (21%)` — four numbers that sum to zero, said to sum to a
hundredth — and every spill under 51 kB printed `0.0 MB`, so one holding forty thousand bytes and one
holding nothing said the same words. **On a cohort run both were fine, and on the small runs an
operator uses to check their wiring both were nonsense.** Durations now scale through µs, ms and s;
sizes through bytes, kiB, MiB and GiB. `MB` was also the wrong name for a count divided by 1024².

### A `None` size said nothing at all (Mi7)

The size is there so an operator can size a filesystem, and an `if let` with no `else` left silence
where the metadata read failed — which reads as an empty spill. It now says the size could not be
read, and that the file itself was written and read back normally.

### Names, exports and destructuring (Mi4, Mi6, Mi8, Mi9)

`FilteredRun::spill_bytes` → **`spill_bytes_on_disk`**, restoring what the value is.
`WhereTheTimeWent` joins `mod.rs`'s re-exports beside `FilteredRun`, which carries it as a public
field. The duration binding is `scoring_wall`, so it no longer shares a scope with `scoring`, the
`ParalogScoringContext`. And `WhereTheTimeWent` **is now destructured where it is read**, which is
what closes the surviving mutation below.

## The tests, and the five mutations that survived

| mutation | before | now |
|---|---|---|
| the `calling` and `writing` durations swapped, shares left in place | survived | **killed** |
| the `fitting` and `scoring` rows swapped whole | survived | **killed** |
| `calling` clocked inside `fit_score_and_write_the_calls` (83 ns) | survived | **killed** |
| the size divisor 1000² instead of 1024² | survived | **killed** |
| the size unit named `MB` while dividing by 1024 | not tried | **killed** |
| a fifth pass added to `WhereTheTimeWent` and left out of the total | survived | **will not compile** |
| the size read before the flush | killed | **killed** (re-run properly — see below) |
| the shares divided by the calling pass rather than the total | killed | killed |
| the total forgetting `writing` | killed | killed |

**Three changes did that work.**

**The fixture gives the calling clock something to measure.** `a_spill_and_an_output` now holds the
spill open for a known 20 ms before handing it over, and the test asserts `calling` is at least that.
`Instant` is monotonic, so it is a hard floor a clock started inside the function cannot reach — under
that mutation `calling` was 83 nanoseconds, which passed `> ZERO` and passes nothing now.

**A first draft also asserted the calling pass was the largest of the four, and that was flaky.** It
passed module-scoped and failed under the full suite, where the writing pass contends for the same
filesystem as six thousand other tests: `calling: 22.7ms, writing: 76.7ms`. It was deleted rather
than loosened, because it pinned nothing the per-pass comparison below does not pin exactly — and a
test that passes on an idle machine and fails on a busy one is worse than no test. **Found by
running the full suite rather than the module**, which is the only run whose conditions are the
ones a reader will meet; the three runs after the deletion all give 6,643.

**Each pass's number is checked against that pass's own duration.** The old test checked that the
four names appeared and that the four shares summed to about a hundred; both survive swapping two
passes' figures, which is the likeliest way this line goes wrong. `what_the_line_says_about` now
pulls the duration and the share printed for one named pass and compares them against that pass's
measured value. **Writing that helper found a defect in itself first**: anchored on the bare pass
name it matched "the calling pass" in the line's own opening clause, so it is anchored on the
separator before each name.

**The size fixture is large enough for the printed digits to be wrong.** 1,200 records rather than
two, so the spill exceeds the 64 KiB write buffer — which is what makes the "read before the flush"
mutation fail — and renders as `647.0 kiB` rather than `0.0 MB`, where a divisor of 1000 or a unit
named `MB` both move the printed text.

**Two tests added beside them**: a spill whose size could not be read says so rather than printing
nothing, and finishing pass one twice keeps the size the first call measured (the early return sits
above the `stat`, and a refactor hoisting the `stat` out of the match would clobber it undetected).

**One mutation is now closed by the type system rather than by a test.** Adding a fifth field to
`WhereTheTimeWent` no longer compiles: the constructor is missing a field and the destructuring in
`what_to_tell_the_operator` does not mention it. That is the better outcome — the previous version
compiled clean, printed the same line, and left the total seven seconds short.

**And one of the review's mutations was malformed.** "The size read before the flush" was written as
an *added* `stat` rather than a moved one, so the real post-flush assignment overwrote it and the
mutation was a no-op that read as a survival. Re-run as an actual move, it kills two tests.

## The report

Rewritten from the run of the committed code. The corrections:

| was | is |
|---|---|
| kept heterozygosity "453 of 13,361, 3.4%", "nine times" | **474 of 13,361, 3.55%**, **8.7 times** — the count matched the string `0/1` and dropped 21 multiallelic heterozygotes; and **11.5 times** compared like for like |
| "this cohort of six spans a factor of ten" | **6.6×** measured (26.08 against 3.96 reads a window); 9.7× is a ratio of two fitted parameters |
| "this run cannot say whether it is the duplication or the shallow samples' fits" | it can: the fits sit 0.72× to 1.35× of each sample's own median depth at these positions, and the run's own window dump says so |
| "the over-covered samples are homozygous reference and the heterozygous sample is at ordinary depth — the one thing that does not fit" | **wrong.** Against each sample's own median, the cluster reads 2.43 / 1.22 / **1.93** / 2.75 / 0.83 / 0.99 — three samples over-covered and **the heterozygous one among them**. The model's story holds there |
| "not what a statistical fluke looks like" | the null, computed: 200,000 random draws of 36 from the 2,311 written records never reached 34 in multi-member clusters (mean 12.6, highest 28) |
| the filter-on runs peak 15 MB lower, with a mechanism for it | **the sign did not replicate** — eight further replicates put them 10 MB higher. Both are 2–3% of 420 MB; the direction and the mechanism are withdrawn, the "below the noise" conclusion stands |
| "31 deletions" | 22 deletions and 9 equal-length substitutions |
| "0.6% of the 200,000 the run analysed" | of the **199,672** the run called |
| the spill "about 7 times a compressed output" | 6.6× the filtered run's own gzipped output, 7.8× the pre-filter run's — two different comparisons |
| the time extrapolated, the disk not | **about 350 GB of scratch** at 3,000 samples and five million records, from a measured 134 bytes a record plus 25 a sample — against about 14 hours of scoring. **The disk binds first, and spec §5 does not price it** |
| scoring "0.05 s … 3.2 to 4.0 µs … 13 to 17 hours" | **47 ms, 3.4 µs, about 14 hours** — the line now prints milliseconds, so the figure has two significant digits rather than one |

**And one thing the report had missed entirely, now its own section**: `SRS3394712` and
`SRS3394712_SRR7279484` are two sequencing runs of **one biosample**, so the plan's standing
six-accession slice is five plants. The benchmark's own
[`rename_dup_samples.sh`](../../../../benchmarks/tomato1/scripts/rename_dup_samples.sh) exists to stop
a cohort VCF seeing the duplicate column, and six of the 63 accessions are such pairs. Where those two
columns agree it is a consistency check on the pipeline, not two samples corroborating each other.

## Validation after the fixes

In the container, via the absolute path to `scripts/dev.sh`:

- `cargo test --all-features --lib --bins --tests` — lib **6,643 passed, 0 failed, 15 ignored**, in each of three consecutive runs
  (6,639 at D1's commit; the four are this step's). One integration failure,
  `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, which is `main`'s.
- `cargo test --lib --all-features ng::run::paralog_filter` — **159 passed**, re-run as the last
  thing before `git add`.
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` — **9 errors in three kinds**,
  counted by kind, none in this step's files.
- `cargo fmt --check` — **9 unique files**, back to the baseline. It was 10: `rustfmt` was run on
  this step's test file alone, and the other nine were left as they were.
- **The standing oracle**: six accessions at `--paralog-fdr 0` give 2,311 records and sha256
  `84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d`, unmoved.
- **The three file relations hold on the six-accession run**, with 36 records flagged.
