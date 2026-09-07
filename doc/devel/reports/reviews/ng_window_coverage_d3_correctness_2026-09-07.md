# Review — window coverage D3: reliability, correctness, and every number the step claims

**Date:** 2026-09-07
**Reviewed:** commit `7fae9525` ("wip: D3 for review") against `bdf80ada`
**Report under review:**
[`ng_window_coverage_d3_2026-09-07.md`](../implementations/ng_window_coverage_d3_2026-09-07.md)
**Raw measurements read:** `/Users/jose/devel/pop_var_caller-window-coverage/tmp/d3/`
**Remit:** correctness and reliability, and a recomputation of every number the step quotes.
Naming, structure and prose went to a second reviewer and are not covered here.

---

## The short of it

**Every statistic in the report reproduces to the digit, and one of its two central inferences does
not follow from the data it cites.** The regression, the residual spreads, the standard errors, the
per-size means, the record counts, the "no psp byte changed" claim and the 10 kb readings all check
out exactly. What does not survive is:

- **the argument that the 26 MB term is not per sample** — the table it rests on cannot show that,
  though the slope in the table above it can, decisively (52 standard errors);
- **one sentence that its own numbers contradict** — "the difference is larger than either spread"
  at 10 kb, where the difference is 5.5 MB and the before arm's spread is 12.5 MB;
- **the claim that the two new tests pin the four footprint numbers** — before this review they
  pinned two of them; the other two were arithmetic on literals that would have gone on passing
  after the code they describe changed;
- **a per-sample term the price omits altogether**: the finalised-windows array, 24 bytes a
  retained covered position, all of it new in this plan.

Eight findings below, four Major and four Minor, no Blockers. All are fixed in this worktree —
tests strengthened, report and spec corrected — and the changes are listed at the end.

---

## 1. Findings

### Major 1 — CONFIRMED. The two new tests did not pin two of the four numbers they claim to

**What the report said:** "Every byte above is asserted by a test rather than derived in prose …
so a struct that grows a field fails a test rather than quietly making this report wrong."

**What the tests did.** `one_samples_measurement_costs_what_the_memory_report_says_it_does` never
called `observe`. Two of its four assertions were arithmetic on literals:

```rust
let sliding_buffer =
    (usize::try_from(shipped.window_bp).unwrap() + 1) * size_of::<CoveredPosition>();
assert_eq!(sliding_buffer, 8_016);          // 501 × 16, with 16 asserted three lines above

let held_back = 16_384 * size_of::<WindowMeans>();
assert_eq!(held_back, 262_144);             // 16_384 is a literal, not read from anything
```

Neither can fail for any change to the code it describes. Make `observe` buffer two entries per
position and `8_016` still passes. Replace `DepthBinWidth::AwaitingWindows(Vec::new())` with
`Vec::with_capacity(10_000)` — which would take the transient from 262 kB to 160 kB — and `262_144`
still passes. The guarantee the report claimed held only for the three `size_of` lines.

**Fixed.** The test now drives the accumulator over 10,249 covered positions and reads `len` and
`capacity` back off `positions`, `counts` and the held-back `Vec`. It fails if the buffer's
high-water mark moves, if the histogram is ever resized, or if the held-back list is reserved
rather than pushed.

### Major 2 — CONFIRMED. "So the term is not per sample" does not follow from the ground table

**What the report said:** 10 kb gives a 5.5 MB difference and 200 kb gives 33.2 MB — "Twenty times
the ground, six times the difference. **So the term is not per sample**".

**Why it does not follow.** Both rows are at **one sample**. At one sample a per-run term and a
per-sample term are the same number, so varying the ground while holding the cohort at one cannot
separate them. A term that were both per-sample *and* per-ground would produce this table
unchanged. What the table shows is that the term grows with the ground, and sub-linearly — which is
worth saying, and is not what was claimed.

**The correct argument is already in the report, one section earlier.** Were the 26.3 MB paid once
per sample, the fitted slope difference would be 26.3 MB a sample. It is 0.82 ± 0.49, which puts
26.3 at **52 standard errors above the point estimate**. That excludes a per-sample term of that
size outright — far more strongly than the ground table could.

**Fixed.** The section now leads with the slope argument and demotes the ground table to what it
actually shows.

### Major 3 — CONFIRMED. "The difference is larger than either spread" is false for one of the two arms

**Recomputed from `tmp/d3/small_before/time{1..5}` and `small_after/time{1..5}`:**

| arm | five readings, MB | mean | range |
|---|---|---|---|
| before | 17.8, 13.6, 5.3, 17.8, 12.1 | **13.32** | **12.5** |
| after | 18.6, 18.6, 19.7, 18.6, 18.6 | **18.82** | **1.1** |

The means, 13.3 and 18.8, and the 5.5 MB difference are all correct. But **5.5 is less than half the
before arm's 12.5 MB spread**, so the parenthesis contradicts the very numbers it quotes in the same
sentence. (Read as standard deviations rather than ranges the claim scrapes through — 5.15 against
5.5 — but the sentence quotes ranges.)

**Fixed.** The parenthesis now states both spreads and says the 10 kb row is the weaker of the two.

### Major 4 — NEW. A per-sample term is missing from the price altogether

`WindowCoverageInProgress::finalised` is a `Vec<(GenomePosition, WindowCoverage)>` — 24 bytes an
entry — holding every window the accumulator has produced that the merge has not evicted.
`SampleWindow::evict_before` drops those in the same call that drops the held records, over the same
stretch of ground. **The report's table does not mention it**, and the spec's §5 mentions it only as
"the ready deque", which is a different thing (the accumulator's own `ready` really is drained to
empty after every record observed, so it costs a record's worth, not a sample's).

Two consequences:

- **The look-ahead's own share was undercounted.** The extra half-window makes a sample hold 250
  more *positions*, and each retained position costs 48 bytes of held record **plus 24 bytes of
  finalised window** = 72. The look-ahead term is 250 × 72 = **18,000 bytes, not 12,000**.
- **The rest of the retained stretch is unpriced and is entirely new in this plan.** The merge held
  records over that stretch before this branch; it held no windows. At one record a base that is
  24 bytes for every base the merge retains — half again what psp mode pays for the records over
  the same ground. How long the stretch is depends on the organiser's release pace, which nothing
  in D3 measured.

**Fixed as far as it can be here.** The report's table now charges 18,000 for the look-ahead, adds
an explicit unpriced row for the rest, and says what would measure it. The
`a_held_record_costs_…` test pins the 24 bytes and the 72. **The unpriced row is a genuine open
item, not a wording fix** — it is the one thing in this step that could still change the per-sample
answer materially.

### Minor 1 — CONFIRMED. 8,016 is the live figure; the run pays 8,192

Measured by driving the accumulator: the buffer's high-water mark is **501 live entries in 512
allocated ones**. `VecDeque` grows by doubling and gives nothing back, so a sample pays 512 × 16 =
**8,192 bytes**, 176 more than the 8,016 the report charged. Small, but it is the allocation that a
memory budget sees. Both figures are now in the table and both are asserted by the test.

**Two related things the report should have said and now does:**

- **The 501 bound is a property of the stream, not of the type.** The buffer holds one entry per
  *observation*; `accumulator.rs`'s own module doc says a stream repeating one coordinate without
  advancing "would buffer without bound". What bounds it is that the merge feeds one record per
  covered position (`for_each_reported_depth`'s three rules make the region partition disjoint).
  The old test comment garbled this into a self-contradiction — asserting the bound while stating
  that repeats add entries.
- **12,000 (now 18,000) is an amortised marginal cost, not a high-water mark.** The held arrays
  double, so 250 more entries cost between nothing and about twice the figure.

### Minor 2 — CONFIRMED. 48 bytes is psp mode's, and the report did not say so

`held_bodies` is pushed only on the `Drawn::Kept` branch — psp mode. Direct mode pushes
`held_observations` instead, and `size_of::<SampleLocusObservations>()` is **120**, so a direct-mode
held record is 32 + 120 = **152 bytes of stack**, before its `Box<[u8]>` of reference bases and its
`Vec<SequenceObservation>` on the heap. The test's doc comment noted the split; the report's table
and headline did not, and the headline figure is the one that gets quoted.

The measurement runs are psp mode throughout, so the two halves of the report were at least
consistent — but "100.2 kB a sample" read as a property of the caller rather than of one mode.
**Fixed:** the headline says psp mode, the direct-mode figures are stated, and the test now asserts
both so neither can be read as the other.

### Minor 3 — CONFIRMED. "Nineteen times wider" needs its convention stated, and the 480 kB is misdescribed

**Two arithmetic slips in the prose:**

- **"0.82 ± 0.49 MB a sample, an interval nineteen times wider than the whole effect."** With
  ± 0.49 named in the same breath the reader computes 0.98 / 0.10 ≈ **ten**. Nineteen is the width
  of the *95%* interval, 2 × 1.96 × 0.4932 = 1.93, which the report never says. Now stated as
  "the 95% interval is 0.82 ± 0.97".
- **"the same 480 kB header that `run_streaming.md` §7.2 already names as three quarters of a human
  open sample."** §7.2 says the opposite shape: the human **open sample** is 480 kB, of which the
  **header** is 357 kB, and three quarters of *that header* is the contig list. 480 is not a header
  and is not three quarters of anything. The downstream arithmetic (480 → 123 → 223) was right;
  only the description was wrong.

**Also fixed in passing:** "the 262 kB transient … is a fifth of the budget" — 262 kB against
§7.2's 500 kB is **52%**, not 20%. The fifth belongs to the charged total, not the transient.

### Minor 4 — CONFIRMED. Three claims asserted as measured that were inferred

1. **"it is the reference the cache now reads once per cover … and the ground each cover spans"** —
   a mechanism claim about a 26 MB term. No arm was built with the reference read removed and no
   heap profile was taken. What the data establish is a cohort-independent term that grows with the
   ground. The reference read is the candidate whose size has that shape; that is an inference, and
   the report now says so and names the `dhat` run that would settle it.
2. **"the 39.5 MB … is dominated by the merge's held records and by the psp source's arena"** —
   contradicted three paragraphs later by "nothing here measures which term it actually is". Now
   stated once, as two candidates and no measurement.
3. **"301 MB, measured 2026-09-07"** (spec §4) and **"100.2 kB a sample in all, measured"** (spec
   §5) — the largest cohort measured was 63, and the per-sample figure is priced by construction.
   Now "priced by construction and pinned by tests", with the extrapolation labelled.

---

## 2. What the measurement can and cannot carry — the statistics, recomputed

**Every figure the report quotes reproduces exactly.** Ordinary least squares of peak against
sample count, 27 rows an arm from `tmp/d3/fine_sweep.txt`:

| | before | after | difference |
|---|---|---|---|
| slope, MB a sample | 39.458 (se 0.369) | 40.276 (se 0.327) | **0.818 ± 0.493**, t = 1.66 |
| intercept, MB | −11.827 (se 9.633) | 14.520 (se 8.544) | **26.347 ± 12.876**, t = 2.05 |
| residual standard error, MB | 38.707 | 34.328 | |

The report's 39.5 / 40.3, −11.8 / +14.5, 38.7 / 34.3, 0.82 ± 0.49 and 26.3 ± 12.9 are all correct.
**The ± are standard errors**, not standard deviations or intervals: 0.493 is
√(0.369² + 0.327²) and 12.876 is √(9.633² + 8.544²). The residual figures are residual standard
errors (√(RSS/25)), not sample standard deviations of the residuals (which are 37.95 and 33.66) —
worth naming, since a reader who recomputes with an n−1 denominator will get a different number.
The report now names all three conventions under the table.

**Is "the one thing the measurement can see" justified at 2.0 standard errors? Only barely, and it
is not robust.** Refitting with the smallest cohorts dropped:

| rows used | intercept difference | t |
|---|---|---|
| all nine sizes | 26.3 ± 12.9 | **2.05** |
| without N = 1 | 25.0 ± 15.1 | 1.66 |
| without N = 1, 2 | 28.3 ± 17.8 | 1.59 |
| without N = 1, 2, 3 | 17.3 ± 21.6 | 0.80 |

The sign is stable; the significance is not, and it leans on exactly the rows most exposed to the
instrument's blind spot (below). **26 MB should be read as an order of magnitude with the sign
known, not as a measured constant** — which is what the report now says.

**"Two runs of the same binary at 63 samples differ by as much as 134 MB": CONFIRMED**, and it is
the *after* arm — 2,589.0 against 2,454.8, a range of 134.2. The before arm's range there is 131.5.

**"The slope did not move" overstates what the fit says.** The point estimate moved by 0.82 MB a
sample, eight times the priced 0.11, and the priced value sits 1.4 standard errors below it. The
data are as consistent with the priced cost as with none. The honest statement is the report's own
next clause — "could not have been seen to" — and the heading now matches it.

### The instrument has a one-sided bias the report did not mention

`scripts/peak_rss.sh` falls through to the `/proc/<pid>/status` path in this container (no `time(1)`
is shipped) and samples `VmHWM` every 20 ms **while `kill -0` still succeeds**. A peak reached after
the last sample that lands inside the process's life is never seen, so **the reading can only err
low**. That is the likeliest explanation for the single-sample outliers — the before arm returning
5.3 MB once at 10 kb against 17.8 twice, and 23.0 / 23.1 / 31.7 at one sample on the 200 kb slice.
Those are precisely the rows the intercept difference leans on, and a low-biased before arm at small
N would manufacture exactly the positive intercept difference observed. The report now carries this
as a caveat on the instrument.

**One more design point about the run:** the arms were executed sequentially — all 27 before runs
(08:02–08:03), then all 27 after runs (08:03–08:04). Anything that drifted over that half hour sits
entirely inside the arm difference. Interleaving the two binaries within each repeat removes it, and
is recorded as what a rerun should do.

---

## 3. What I checked and found sound

**The per-construction arithmetic, all of it.**

| claim | recomputed | verdict |
|---|---|---|
| histogram 80,200 | 50 GC bins × (400 + 1 overflow) columns × 4 bytes = 80,200; `counts` is `vec![0; cell_count()]`, so `capacity == len` and it is never resized | **sound** |
| 8,016 | 501 × 16 — correct as the live figure; see Minor 1 for the allocated one | sound with a caveat |
| 12,000 | 250 × 48 — correct for the record half; see Major 4 for the window half | undercounted |
| 100,216 | 80,200 + 8,016 + 12,000 — adds up | sound |
| 262,144 | **verified by driving the accumulator**: `AwaitingWindows(Vec::new())` plus pushes, capacity 16,384 at 9,999 held windows, × 16 = 262,144. No `reserve` anywhere on that path | **sound and now observed** |
| "a fifth of 500 kB" | 100,216 / 500,000 = 20.0% | sound |
| "108 → 208 kB on tomato" | §7.2's tomato open sample is 108 kB; 108 + 100.2 = 208.2 | sound |
| "480 → 580 kB on human" | §7.2's human open sample is 480 kB; 480 + 100.2 = 580.2, over the 500 kB budget | sound |
| "480 to 123 so the total to 223" | §7.2: sharing the contig list takes 480 → 123; 123 + 100.2 = 223.2 | sound |
| "100 MB across a thousand, 301 MB across three thousand" | 100,216 × 1,000 = 100.2 MB; × 3,000 = 300.6 MB | sound |
| "spec §3.4's own estimate of 80 MB and 240 MB" | 80,200 × 1,000 and × 3,000 — the histogram alone, as stated | sound |
| "39.5 MB a sample is 79 times the 500 kB budget" | 39.5 / 0.5 = 79.0 | sound |
| "halving the bins halves 80.2 kB to 40.2" | 50 × 201 × 4 = 40,200 | sound |

*(The bytes-a-sample figures in this table are the report's originals; Minor 1 and Major 4 revise
8,016 → 8,192 and 12,000 → 18,000, which carry the total to 106,392 and the downstream figures to
214 / 586 / 229 kB and 106 / 319 MB. All of those were recomputed too.)*

**The three-size table is the nine-size sweep's means, and no runs are mixed.** Recomputed from
`fine_sweep.txt`: 25.93 / 59.13 at one sample, 216.87 / 240.27 at six, 2,447.87 / 2,539.00 at 63 —
the report's 25.9 / 59.1, 216.9 / 240.3, 2,447.9 / 2,539.0 exactly. The earlier four-repeat
three-size run under `tmp/d3/before/`, `after/`, `before_rep1..4/` and `after_rep1..4/` gives quite
different numbers (17.9 / 60.7 at one sample, 2,675.0 / 2,530.9 at 63) and **none of them appear
anywhere in the report**. No silent mixing.

**Both binaries produce identical VCF record counts — and at all nine sizes, not just three.**
From the `vcf_records` column of the six `memory_vs_samples.tsv` files: 21, 62, 108, 196, 245,
2,120, 2,587, 4,208, 5,492, identical in every one of the six runs. The report's 21, 245 and 5,492
are the three it quotes; the full agreement is stronger than what it claimed.

**"This plan changes no psp byte": CONFIRMED, two ways.** `git diff --name-only a33ada0f 7fae9525
-- src/ng/psp/` is empty — none of `block.rs`, `chain_ids.rs`, `footer.rs`, `header.rs`, `index.rs`,
`reader.rs`, `record.rs`, `segmentation_section.rs`, `trailer.rs` changed. And operationally the
question does not arise: `memory_vs_samples.sh` points both binaries at the same
`tmp/d3/psps63` directory, so the same files were read either way.

**The harness is sound in the two ways the report claims.** It names its binary (`bin=$1`, with an
`[ -x ]` guard) rather than taking the newer of two target directories, which for a two-build
comparison would have measured one build twice. And it takes the cohort at each size as a prefix of
one `sort`ed list, so a step from *N* to *M* adds samples and changes nothing else — `all=$(ls
"$psps"/*.psp | sort)` then `head -"$n"`.

**The accumulator code behind the numbers.** `counts` is allocated whole at `new` and only ever
written through `fold`'s `cell.saturating_add(1)` — never pushed to, never resized. The held-back
list is `std::mem::take`n at the fit and consumed by the fold loop, so it is genuinely freed at the
ten-thousandth window and the "then never again" holds. `Unfittable` latches, which is what keeps
the held-back list bounded for a sample with no positive median. The `ready` deque is drained to
empty after every record the cache observes, so it is correctly not a per-sample residency.

**One thing worth recording as sound because it looked wrong at first.** In `finalise_centre`'s
shrink-left step, entries sharing a coordinate are popped one at a time while
`distinct_positions_summed` is decremented only on the last of them. That is correct — the guard
`self.positions.len() == 1 || self.positions[1].position != front.position` — and it is what makes
the floor count distinct coordinates rather than observations, which is the whole point of that
counter. Repeated coordinates cannot arise from today's merge, but the code handles them, and the
501-entry bound is the one place where "cannot arise today" is load-bearing (Minor 1).

---

## 4. What I changed

**`src/ng/window_coverage/accumulator.rs`** — rewrote
`one_samples_measurement_costs_what_the_memory_report_says_it_does` so that it measures instead of
re-multiplying. It now builds the shipped configuration, drives `observe` over 10,249 covered
positions at one record a base, and asserts:

- `counts.len() == counts.capacity() == 50 × 401` at construction **and after ten thousand
  windows** — 80,200 bytes, never resized;
- the buffer's observed high-water mark: `positions.len()` peaks at 501 (8,016 bytes live) and
  `positions.capacity()` at 512 (**8,192 bytes allocated**);
- the held-back list observed one window short of the fit: `len == 9,999`, `capacity == 16,384`,
  262,144 bytes;
- the two whole-pass terms at 88,216 bytes.

The doc comment now states that the 501 bound is a property of the stream rather than of the type,
and why the ready deque is not priced.

**`src/ng/run/cohort_merge/observation_cache.rs`** — extended
`a_held_record_costs_what_the_memory_report_says_it_does` to pin what the report was reading as one
number and is three:

- psp mode's held record at 48 bytes (32 + 16), unchanged;
- **direct mode's at 152** (32 + `size_of::<SampleLocusObservations>()` = 120), asserted so the two
  cannot be confused;
- **a finalised window at 24 bytes**, and a retained position in psp mode at **72**;
- the look-ahead's own share at 250 × 72 = **18,000 bytes**, alongside the old 250 × 48 = 12,000.

Doc comment rewritten to say which mode each figure is, that 0.99 records a base comes from C3's
198,710 covered positions over 200 kb, and that the figure is an amortised marginal cost rather
than a high-water mark.

**`doc/devel/reports/implementations/ng_window_coverage_d3_2026-09-07.md`** — corrected:

- headline and table: **106.4 kB a psp-mode sample** (80,200 + 8,192 + 18,000 = 106,392), with an
  explicit unpriced row for the finalised windows over the rest of the retained stretch, and the
  three qualifications (mode, stream-dependence, amortisation) spelled out;
- downstream figures carried through: 21% of budget, 108 → 214 kB, 480 → 586 kB, → 229 kB,
  106 MB / 319 MB;
- §7.2 described correctly (480 kB open sample, 357 kB header, contig list three quarters of the
  header);
- ± labelled as standard errors, with each component given, and the residual figure named as a
  residual standard error;
- "nineteen times" restated as the 95% interval, 0.82 ± 0.97;
- "the slope did not move" → "the fit cannot resolve a change", with the 1.4-standard-error
  distance to the priced value stated;
- the `VmHWM` sampler's one-sided bias added as a caveat, and the sequential-arms design added as
  assumption 5;
- the "not per sample" section rebuilt on the slope (52 standard errors) with the ground table
  demoted to what it shows, and the 10 kb spreads stated correctly;
- the 26 MB mechanism marked as an inference with the `dhat` run named;
- the 39.5 MB attribution stated once, as two unmeasured candidates;
- the record-count claim widened to all nine sizes and the psp-unchanged claim given its evidence;
- the three-size table annotated as the sweep's means, with the earlier four-repeat run explicitly
  not quoted;
- the 262 kB transient corrected from "a fifth of the budget" to half of it again;
- a follow-up added for direct mode, which was priced but never run.

**`doc/devel/ng/spec/window_coverage.md`** — §4 and §5 brought in line: 319 MB at three thousand
(labelled as an extrapolation from a per-sample figure, largest cohort measured 63), 106.4 kB a
psp-mode sample, the unpriced finalised-windows term named, direct mode's 152 bytes recorded, and
the 26 MB attributed to the slope argument rather than asserted as a measured mechanism.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` was touched.

---

## 5. Verification

```
$ scripts/dev.sh cargo test --lib --all-features
test result: ok. 6380 passed; 0 failed; 15 ignored; 0 measured; 0 filtered out; finished in 52.44s

$ scripts/dev.sh cargo check --lib --tests --all-features
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 12.70s

$ scripts/dev.sh rustfmt --check --edition 2024 src/ng/window_coverage/accumulator.rs
(no output — clean)

$ scripts/dev.sh rustfmt --check --edition 2024 src/ng/run/cohort_merge/observation_cache.rs
Diff in .../observation_cache.rs:260:
Diff in .../observation_cache.rs:971:
Diff in .../observation_cache.rs:1061:
Diff in .../observation_cache.rs:1761:
```

**Four hunks in `observation_cache.rs`, the same four as before this branch** — lines 260, 971,
1061 and 1761, all far above the test at 1,902 that this step and this review touched. Not a
regression, and unchanged in count by either.

---

## 6. What is still open

1. **The finalised-windows term over the retained stretch is unpriced** (Major 4). It is 24 bytes a
   covered position and all of it is new in this plan. Measuring it needs the organiser's release
   distance, which nothing in D3 recorded — a `dhat` run at one sample, or an instrumented
   high-water mark on `finalised.len()`, would settle it.
2. **Where the 26 MB actually lives** (Minor 4.1). One `dhat` run at each of the two grounds, one
   sample, names the allocation instead of nominating a candidate.
3. **The intercept difference deserves a rerun with the arms interleaved** and, if it is to carry
   any weight, an instrument that cannot read low — the 20 ms `VmHWM` poll is the wrong tool for a
   three-second run whose peak may fall at the end.
4. **Direct mode has never been run under this measurement.** Its per-sample price is three times
   psp mode's on the look-ahead term alone, before the records' heap.
