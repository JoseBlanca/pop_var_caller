# Code Review: ng_paralog_filter_d4

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), one sub-agent in an isolated worktree
**Scope:** step D4 of the hidden-duplication filter plan — what the filter scores and flags by kind
of record, and the comparison against production's own filter
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the two by-kind tables and the production comparison in the step's
  implementation report, and
  [`scripts/ng_paralog_filter_by_kind.py`](../../../../scripts/ng_paralog_filter_by_kind.py). **No
  library code changed in this step** — confirmed by the reviewer against `git diff HEAD`.
- **One sub-agent**, in its own detached worktree at `4e061476`, which rebuilt both binaries and
  **reproduced both ng runs and both production runs from scratch** rather than reading the
  author's outputs.

### 2. Verdict

**Request-changes.** Zero Blockers, five Major, nine Minor.

**Every cell of both tables reproduced exactly** — written, scored, flagged and all five quantiles
across six kinds and two runs, re-derived with the reviewer's own classifier and index arithmetic —
as did the whole production comparison and the tract headline counts. **What failed was again the
explanation**, and one of the failures makes the step's conclusion *stronger* rather than weaker.

### 3. Execution status

No library code in this step, so the gate is `4e061476`'s: lib 6,643 tests pass, 0 failed, 15
ignored; one integration failure which is `main`'s; clippy 9 errors in three kinds; `fmt` dirty on
9 unique files.

### 4. Findings

#### Major

**M1: "far beyond anything in either run" is not quantified, and is false when quantified**

**Confidence:** High.

The report said a tract would need coverage "far beyond anything in either run" to reach the cut on
coverage alone. Measured: at σ₀ = 0.092 a zero-read record reaches 4.65 at **1.372 copies, 4.04 σ₀
above one copy**, and the highest-scoring tract in the run sits at **1.310, 3.37 σ₀**. **The gap is
0.68 σ₀ — about 17 reads in a 364-read window**, and the run came 83% of the way. The arm looks flat
in nats because it is flat and then abrupt: 1.3 at 1.30 copies, 3.5 at 1.35, 6.2 at 1.40, 56.7 at
two.

**Why it matters:** the report used the flatness to argue the coverage-only arm is safely inert.
It is not a safe zero, and that is a stronger argument for §8's allele term than the one made.

**M2: two different records were quoted as one**

**Confidence:** High.

"The deepest tract window measured is 1.42 times its sample's fitted one-copy scale, and its ratio
is 1.48." The deepest window is **chr6:150,466,707 at 1.4156× the scale, whose ratio is −0.6252** —
on the plateau. The ratio of 1.4843 belongs to **chr11:9,161,348 at 1.2033×**. Reproduced by the
orchestrator on its own run.

**The GC multiplier decides, and this is the fourth step in a row to trip on it**: the deep-but-quiet
tract sits at GC 0.295 where a window is expected to carry more, the loud one at GC 0.457 where it
is expected to carry less. **Raw depth over the one-copy scale is not the axis the score reads.**

**M3: "at six samples only one pair coincides" — no pair coincides**

**Confidence:** High.

All 21 tomato tract ratios are distinct. The claim was a misreading of the analysis script's own
last line, which printed "1 of them share one ratio exactly" to mean *the largest group has size
one* — that is, no repetition at all, with an arbitrary tie-break value attached. **The script
produced the error and the script is fixed**; it now says "no two tracts share a ratio" in words.

**M4: the mechanism given for the constant is the wrong one, and the right one is the finding**

**Confidence:** High.

The report said the ratio is constant "where a sample's coverage is at one copy". It is that same
constant for **every copy number from about 0.2 to about 1.1** — a tract at half a copy scores what
a tract at one copy scores. What sets the plateau's width is **σ₀, not the coverage**: at 0.092 the
carrier hypothesis at 1.5 copies is more than five standard deviations away and its branch
underflows, so the ratio answers the prior. Confirmed directly from the run's window dump: the 469
tracts within 0.0001 of the plateau span raw depths of **0.202 to 1.272 times the scale**.

**The report's reading suggested an arm that works and finds nothing. What is true is an arm blind
over a band a copy wide that then turns on abruptly** — which is what a reader deciding about §8
needs, and which is why M1's margin is thinner than the nat count suggests.

**M5: "agrees wherever it speaks" overstates the scope**

**Confidence:** High.

Production removed **6,223 records** over the whole of `regions.bed`; **17 of them — 0.27% — fall
inside the two intervals ng ran.** The body of the report was explicit about the restriction; the
summary sentence was not.

#### Minor

- **Mi1** — the biallelic SNP span on HG002 is **432.53** nats, not "over 440".
- **Mi2** — **p99 is the maximum by construction wherever n ≤ 100**, since `int(0.99·n)` is the last
  index; that is five of six tomato rows and two of six HG002 rows, and nothing said so.
- **Mi3** — "the ninetieth percentile equal to the median equal to the lowest" holds at the table's
  two decimals; at four, p90 is −0.6566 against −0.6572.
- **Mi4** — "share one ratio exactly" cannot be established from a VCF that prints four decimals.
  True in substance and not shown by the file.
- **Mi5** — "the indel arm is 22 of 284" contradicts the same paragraph's "19 indels": 22 counts the
  3 equal-length substitutions, which the report's own reasoning excludes.
- **Mi6** — 5.23 against a cut of 6.2500 is **1.02 nats short**, so "within a nat" is wrong.
- **Mi7** — **601 = 21 + 580 pools two incomparable runs**: different species, cohort size, cut and
  depth — and the two disagree about the plateau, which is the reason not to add them.
- **Mi8** — "equal-length multi-base substitutions" flatters the kind: **8 of the 9 on tomato change
  exactly one base**, so they are mostly SNPs written with a shared flanking base, and "their depth
  is not biased at all" is asserted rather than measured.
- **Mi9** — the script's modal-ratio line is meaningless when nothing repeats, which is what
  produced M3.

#### Why-claims flagged as asserted rather than measured

- "a deletion depresses depth across its own span in exactly the samples that carry it" and "their
  depth is not biased at all" — the stated reasons for the kind split. Neither is measured here; the
  split stands on the outcome, which is.
- "A deletion's ratios reach the far tail and a tract's do not, **because** a deletion carries the
  allele signal and a tract does not." Too strong: a tract's coverage arm reaches 56.7 at two copies,
  so this is a statement about the copy numbers these tracts had.

### 5. What is right, and was checked

Both by-kind tables cell by cell; the tract counts and the plateau value; the deletion-versus-tract
contrast; the entire production comparison including the 1,867 shared records, the 17/27/17/0 split,
the 9-never-called and 10-called-and-kept breakdown, both calibrations, and the two clusters holding
15 of the 17. **The BED conversion was checked and is right**: `$2 > lo && $2 <= hi` is the correct
half-open zero-based to one-based mapping, and no record sits on any of the four boundaries, so all
four candidate conventions give the same 1,879.

**The classifier's `STR` key is the filter's own key** — `pass_one.rs` chooses the row shape on
`record.is_repeat_tract()`, and `encode.rs` writes `STR` under the same condition, so the two
cannot disagree. That was worth checking, because the whole tract row depends on the analysis
agreeing with the code about what a tract is.
