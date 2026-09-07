# Fixes applied — ng_paralog_filter_d1

**Date:** 2026-09-07
**Review:** [ng_paralog_filter_d1_2026-09-07.md](ng_paralog_filter_d1_2026-09-07.md) (0 Blocker /
9 Major / 12 Minor, Request-changes; two sub-agents in isolated worktrees)
**Branch:** `ng-paralog-filter`

## What was done

**All nine Majors and all twelve Minors applied. Nothing deferred.**

Five of the Majors were wrong sentences in the implementation report and four were blind spots in
the two measuring tools. One change went further than the review asked and adds library code: the
run now prints what each sample's coverage fit came to, which is what makes three of the findings
impossible to repeat rather than merely corrected once.

**Everything below was re-measured after the changes, on a run of the code being committed.** The
report was rewritten from that run rather than edited in place, because the first draft's figures
came from a run of a different build.

## The library change, and why a measurement step made one

**M3, M4 and Mi3 are one defect wearing three faces: the run printed nothing about the fit it had
just performed.**

The report needed the depth a record's window was compared against. Nothing printed it, so the first
draft used the obvious stand-in — the median window depth over the 217 records the run wrote — and
called it "one copy's depth". Variant sites are exactly where excess coverage collects, so that
median is **8.49 reads** where the fit's own one-copy level is **5.22**. The flagged record's window
therefore read as 2.24 copies in the report where the model sees about 3.5, near its winsor cap of
four. The number a reader carries away was wrong by a factor of 1.6, in the direction that makes a
flagged record look *less* duplicated than the model thinks, and **nothing in the run's output could
have corrected it** — the reviewer recovered the truth only by re-implementing the scorer and
inverting all 217 ratios.

Plan step D1 asks for **"the fit's outcome"**. The run was reporting the fit's *acceptance count*
(`over 1 of 1 sample(s) with a coverage model`), which is a different thing. So:

- **[`ParalogScoringContext::what_each_fit_came_to`](../../../../src/ng/run/paralog_filter/scoring_context.rs)**
  returns one `WhatTheFitCameTo` per sample — the fitted one-copy depth and σ₀ — or `None` where the
  fit was refused. **A small owned summary rather than the models themselves**, deliberately: a
  caller holding a `SingleCopyCoverageModel` could ask it for a copy number at a GC of its choosing
  and index the answer against the wrong sample, which is the shape of spec §6's trap 3 and the
  reason C1's review made those accessors private. Two numbers per sample cannot be misindexed into
  a score.
- **[`what_to_tell_the_operator`](../../../../src/ng/run/paralog_filter/finish.rs)** names every
  fitted sample up to ten and gives the spread past that — the lowest, the median and the highest of
  both quantities. A run is up to several thousand samples (spec §4) and a report is read by a
  person; printing the first ten and stopping would leave a reader believing that is the cohort.

At one sample it prints:

```text
  sample SRS3394712 fitted one copy at 5.22 reads a window, scatter 0.300
```

**σ₀ = 0.300 is the number the reviewer's inversion had reconstructed**, to three decimals — so the
reconstruction is corroborated by the run, and the run no longer needs it.

**Two tests, and seven mutations against them, all seven killed:**

| mutation | verdict |
|---|---|
| the ends of the spread swapped | killed |
| the median taken as the lowest | killed |
| the cap raised past any cohort, so the spread branch never fires | killed |
| the scatter printed where the depth belongs | killed |
| every line naming the first sample | killed |
| a sample's depth and scatter swapped | killed |
| the fit reporting σ₀ where the one-copy depth belongs | killed |

**The fixtures give each sample a different fitted depth for exactly this reason.** A cohort fitted
at one depth cannot tell a line that prints each sample's own number from one that prints the first
sample's number twelve times — which is this plan's recurring fixture failure (C1's Blocker: every
coverage fixture on a single GC bin; C2's: every record on contig 0). `a_fittable_histogram_peaking_at`
shifts the peak bin, and both tests assert the spread is real before asserting what is printed.

**The mutation runner itself was fixed mid-sweep**, which is worth recording because this plan has
been bitten by its own tooling four times. Its first pass reported one mutation as NOT-APPLIED: the
pattern contained a `/` — `values.len() - 1) / 2` — which closed perl's `s///` early. **The runner
reported that as "not a result" rather than as a pass**, which is what it was built to do; both
sides now travel through the environment instead of into the script text. It also restores before it
judges, verifies the restore landed by comparing against the backup, and treats a run that produced
no verdict line as an error rather than a pass. After the sweep both source files were confirmed
byte-identical to their pre-sweep backups.

## The tools

### M6 — the relations could not see the filter losing its own marks

Both relations strip the two INFO keys and the four header lines before comparing, so a filter that
stopped writing them satisfied both. Now they are **asserted**: present on the dropping and tagging
runs, **absent** on the off run, and a record carrying one of the two keys without the other is a
failure. An on-run whose records carry no ratio at all is a failure too.

### M7 — a flagged record was compared against nothing

Added a third relation: **a flagged record is its off-run record with the verdict spliced in.**
Undoing the two INFO keys and the `FILTER` splice must give back the off run's line exactly. The
splice's inverse is read off
[`patch.rs`](../../../../src/ng/run/paralog_filter/patch.rs)'s `append_filter_column` — the id is
joined with `;` unless what was there was `PASS` or `.`, in which case it replaced the column — so
a flagged record left carrying only `hiddenParalog` must have an off-run value of `PASS` or `.`, and
anything else is reported as the splice having replaced what it should have joined.

### M8 — `(CHROM, POS, REF, ALT)` is not an identity

Replaced with a positional walk plus a length check. All three runs call the same cohort over the
same ground and write in the same order, and tagging removes nothing, so the off run and the tagging
run are the same list and can be walked together. **Exact under duplicate keys, and shorter.** The
length check is itself something the tool lacked: nothing previously noticed a short tagging file.

### M9 — an empty run report read as a quiet section

`sed` prints nothing and exits 0 when its marker is absent. The harness now fails, naming the
likely cause (the report's first line reworded).

### The Minors

| finding | what changed |
|---|---|
| Mi1 percentile convention | stated in the report: the value at sorted index `int(p·n)`, and that three percentiles differ under linear interpolation |
| Mi2 correlation over-read | now gives `r² = 0.208` and Spearman 0.364, and attributes the shortfall to the ratio's convexity in copy number and to the read count pushing the other way — **not** to `F`, which cannot move the ranking at one sample |
| Mi3 depth table not reproducible | the report has a section giving the `NG_WINDOW_COVERAGE_FILE` recipe, the row format, that the two fields are `f32` bit patterns written as decimal `u32`, and that it must be the tagging run's dump |
| Mi4 depth compared with a per-position figure | the report says the 8.49 is the median at written variant sites, is biased upward, and is not comparable with `CLAUDE.md`'s figure; the unbiased number is the fit's 5.22 |
| Mi5 peak resident from n = 1 | reported as a range over nine runs, 260.9 to 267.7 MB, with the observation that the dropping run was highest every time |
| Mi6 `mark before this run` | replaced by an assertion that the mark is zero, with a comment saying what a future edit would break |
| Mi7 `ru_maxrss` units | converts on macOS, exits naming the platform on anything else |
| Mi8 unread `.comparable` files | only the off run's is written, since it is the only one whose hash means anything |
| Mi9 unattributed binary | the harness prints the binary's path and build time |
| Mi10 substring match on `FILTER` | matches whole `;`-separated ids |
| Mi11 a zero-record run | prints a warning saying nothing below it is a measurement of anything |
| Mi12 five copies of one shift | a length mismatch is now reported and stops, rather than falling through to a positional diff |

**Two further guards added while fixing the above**, neither asked for: the off run carrying a
calibration line, or a filtering run lacking one, is now a failure rather than a printed note; and
the harness prints a marked line when nothing was flagged, so the empty case cannot be read as a
result in the harness's output as well as in the checker's.

## The report

Rewritten from the run of the committed code. The five Major corrections:

| was | is |
|---|---|
| "the two are separated by 0.14 on the ratio axis" | 13.4586 apart; **0.1384** is the kept record's distance to the cut |
| "what differs is the GC" | the read counts differ too — `AD 0,15` against `AD 0,19`, worth `4 × 0.6831 = 2.73` nats, a fifth of the gap |
| "the model's expected one-copy depth is not the same number at the two places" | with the size, and with the counterfactual: the kept record would have scored about 19.4 and been flagged at the other GC — carried as a reconstruction, explicitly, since the GC multiplier is still not printed |
| "19.01 reads against the run's median of 8.49 — 2.24 times one copy's depth" | 19.01 against the **fitted** one-copy level of 5.22 — about 3.5 copies, near the winsor cap of four |
| "the allele split and the inbreeding coefficient move it too" | the allele split does; `F` is locus-invariant by construction and at one sample cannot move any record relative to any other |

**The arithmetic in the second row was checked independently rather than taken from the review.** The
scorer enumerates seven carrier configurations, `(T,m) = (3,1),(4,1),(4,2),(6,1),(6,3),(8,1),(8,4)`
([`locus_score.rs`](../../../../src/ng/paralog/locus_score.rs), pinned by
`default_carrier_configs_are_the_seven_expected`), so `m/T` never exceeds ½ and `ln 0.5` really is
the duplication story's best branch for a hom-alt read; `ε = 0.01`
([`model_params.rs`](../../../../src/ng/paralog/model_params.rs) `DEFAULT_PSEUDOCOUNT_VAF`) gives the
other term.

## Validation after the fixes

In the container, on this worktree, via the absolute path to `scripts/dev.sh`:

- `cargo test --all-features --lib --bins --tests` — lib **6,639 passed, 0 failed, 15 ignored**
  (6,637 at `5e7f20da`; the two are this step's). One integration failure,
  `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, which is `main`'s and fails at
  the merge base too.
- `cargo test --lib --all-features ng::run::paralog_filter` — re-run as the last thing before
  `git add`.
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` — **9 errors in three kinds**,
  counted by kind, the same three the merge base has and none in this step's files.
- `cargo fmt --check` — dirty on the same **9 unique files**, none of them this step's.
- **The standing oracle**: six accessions at `--paralog-fdr 0` give 2,311 records and sha256
  `84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d`, unmoved by the report line.
- **The falsification suite**: ten broken inputs, all ten now caught, and an untouched triple still
  passes.
