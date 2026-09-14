# Review: B7, what converting the caller's maths to libm cost and what it moved

*Reviews [`portable_float_B7_measurements_2026-09-15.md`](../implementations/portable_float_B7_measurements_2026-09-15.md)
at its draft of 2026-09-15, before commit. Every figure below was recomputed from the files under
`tmp/` that the report names; the scripts used are in `tmp/review_B7/`. No cargo build or caller run
was made.*

## Verdict

**Commit after fixes.** The headline speed figures are right to the decimal: the fit costs 27% on
macOS and 29% on Linux, `call-from-psps` 1.7% and 2.0%. Every checksum, line count and VCF
comparison reproduces, and default-parameter calling is identical record for record. Three
statements go beyond their evidence, though:

- the direction of the noisy macOS bench passes is stated backwards;
- the `estimate-contamination` output check compares two files that contain no estimate;
- peak memory is summarised as inside a spread it mostly falls outside.

None of the three changes the recommendation at Checkpoint B. Each needs rewording before the
report is committed.

## Findings

### Major

**M1. §2.3, "In each, one pass sits near A3's figure for the same case and the other above it, so
the change there is if anything overstated."** The direction is reversed. A pass that ran slow
*raises* the unchanged mean, which *shrinks* the percentage. Measured against the unchanged pass
that sits near A3's figure (`tmp/ab2/benches_macos2/*.{1,4}.std.log`,
`tmp/baseline_A3/macos_benches/`):

| case | reported | against the unchanged pass near A3 | against A3's mean |
|---|---:|---:|---:|
| site quality, 63 samples | +10.8% | +16.9% | +15.1% |
| SNP/indel joint fit, 8 samples | +13.7% | +18.1% | +20.2% |
| SNP/indel joint fit, 32 samples | +8.3% | +11.8% | +11.7% |

So the three figures are *understated*, by 3 to 6 points. The claim that 32 samples is the cheapest
case of each half still holds on those figures. The lower end of the macOS joint-fit range in §1
("+8 to +21%") is probably nearer +10%. **Fix:** write "so the change there is if anything
understated. Against the pass near A3's figure it is +17%, +18% and +12%." Consider re-running
those two benches on macOS if the Checkpoint B decision leans on the low end.

**M2. §3.1, "`estimate-contamination` writes identical JSON from the two builds over the 20
regions."** The JSON files are identical, but they hold no estimate.
`tmp/measure_B7/outputs/contamination_{std,float}.json` report `"estimated": 0`,
`"not_identified": 4` and `"markers": 0`. Every sample has `"contamination": null`, with the reason
"too few positions where the cohort varies". The logs say "0 of 4 libraries … estimated over 0
varying positions". The only numbers that could differ are counts of positions. So this check
cannot show whether the converted `ln`/`exp` in `contamination.rs` moved anything, and the gap the
B5 review raised is still open. **Fix:** say the JSON is identical because four low-depth samples
give no estimate, and that the contamination estimate's parity is therefore unchecked; add it to §5.
The same logs time the two runs, 624 s unchanged against 942 s with the change, under the load §2.1
discards and on an unrecorded platform. That is +51%, above the fit's 27–29%. Either quote it in §5
as a single loaded run, or time the command interleaved.

**M3. §1, "peak memory, every command | within the unchanged build's own spread", and §2.4,
"falls inside or a few per cent from A3's range … in both directions."** Compared with A3's ranges
(`tmp/measure_B7/*/runs.tsv`, `fit.tsv` against `tmp/baseline_A3/{linux,macos}_{runs,runs160}`,
`macos_fit20`, `redo/linux_fit20_macos_psps`), 13 of 14 command-and-platform cells fall outside
A3's range. Only macOS `call-from-psps` over 20 regions is inside. The largest departures:

- Linux `generate-psps`, 160 regions: 409.7–491.4 MiB against 406.1–438.8, so 12.0% above;
- Linux `generate-psps`, 20 regions: 6.8% above;
- Linux `call-from-psps`, 20 regions: 5.7% above;
- macOS `call-from-psps`, 160 regions: 5.1% below;
- macOS `call-from-alignments`, 20 regions: 5.0% above.

The departures go both ways, and one build's own repeats spread by up to 20% (Linux
`generate-psps`, 410–491 MiB). That looks like noise, not a cost, but "within the spread" is not
what the files show. **Fix:** in §1 write "no consistent direction; within 12% of A3's range". In
§2.4 give the largest departure and its command, and say the runs came from a different session
than A3's.

### Minor

**Mi1. §2.2, "Nor, in effect, does `generate-psps`: on Linux the change was the faster build, and
on macOS its fastest run is 0.3% above the unchanged build's slowest."** The report calls
`call-from-psps` "slower by a small, consistent amount" because its change-build runs all sit above
the unchanged ones (a 0.9% gap on macOS). The macOS `generate-psps` runs show the same pattern: all
three change-build runs are above all three unchanged runs, and the change is slower in every round
(`tmp/ab2/gen_macos/gen.tsv`). On Linux the pattern reverses: all three change-build runs are
faster, by 1.0% in the mean. Also, `call-from-alignments` on Linux is slower in 4 of 5 rounds.
Complete separation of three runs against three happens by chance 1 time in 20. Here it happened on
both platforms, in opposite directions, for a command that does little floating-point maths. That
suggests a ~1% effect of the binary itself, and it also bounds how much the 1.7–2.0% for
`call-from-psps` can be read as libm's cost. **Fix:** give `generate-psps` as "+0.8% on macOS and
−1.0% on Linux, each with the builds' ranges apart, so a difference of about 1% in either direction
arises without a maths cost". Replace "does not move" for `call-from-alignments` with the figures.

**Mi2. §2.1, the load description.** Checked against `tmp/ab2/hostload.tsv`:

- "The two change-build fit rounds that overlapped the longest burst took 578.6 s and 579.1 s."
  The longest burst ran from 22:31 to 22:44, about 12 minutes. In it, the media-analysis or
  PDF-preview service was at 40% or more of a core in most samples. It overlapped change round 1
  (579.1 s) and unchanged round 1 (451.9 s), not change round 0. Change round 0 (578.6 s) saw only
  single-sample spikes. The conclusion survives in a better form: the change's round inside the
  burst took 579.1 s and its round outside 578.6 s, and the unchanged round inside it, 451.9 s,
  sits inside the unchanged range.
- "with nothing else running but macOS's own background services." One sample shows `Visual`
  (the editor, from its path) at 113.7% during Linux `call-from-psps`, and `claude` reaches 35%.
- "in bursts." During the `generate-psps` runs the PDF-preview service held about one core
  continuously: 19 of the 23 samples show it or the media-analysis service at 40% or more. It was
  also near 98% through most of the Linux passes of the site-quality, delimiter, pileup and psp
  benches (`benches_linux/load.txt`).
- "host load 1.7 to 11." Within the runs the minimum is 2.36. The 1.72 sample was taken 5 s before
  the first run.

**Fix:** correct the burst sentence as above, and name the editor and the Claude session. Say the
PDF-preview service ran through the `generate-psps` and Linux bench passes.

**Mi3. §2.1, "at a load of 26 on 18 cores the unchanged fit ranged from 280 s to 496 s between
runs."** 280 s is A3's minimum, from a quiet session two days earlier (`tmp/baseline_A3/macos_fit20`).
The loaded session has a single unchanged macOS fit, 496.24 s (`tmp/ab/macos/ab.tsv`). No file
records a load of 26. **Fix:** "the unchanged macOS fit took 496 s under load, against A3's 280–288
s; on Linux 493–567 s against 448–459 s (`tmp/ab/*/ab.tsv`)". Drop the load figure or cite its
source.

**Mi4. §2.3, "All 11 draws hash identically in the two trees, and every SNP/indel fit runs exactly
eight passes in both."** The check was run only on macOS (`probe_{std,float}_macos.txt`). The
unchanged Linux tree draws through glibc, which A2 found differs from libm on about one `exp`
argument in ten, so the Linux draws were not checked. The eight passes are the bench's cap
(`PASSES = 8`), and every fit reports `converged false`, so the pass count says nothing about the
work done inside each pass. **Fix:** say "on macOS". Say the eight passes are the cap. Either run
the check on Linux or say the Linux draws are assumed the same, because the draws round to integer
counts.

**Mi5. §2.3, macOS, "no case moved beyond A3's own run-to-run spread."** 7 of the 27 cases did,
though none by more than 1.4%. Examples: pileup coverage 10 at −1.4% against A3's 0.7%, psp
`walk_full/shallow_10_reads` at +0.9% against 0.1%, and `walk_heads/shallow` at +0.9% against 0.3%
(`benches_macos2/table.md`). **Fix:** "all 27 within 2.6% of A3's mean; only the three psp cases
whose A3 passes were 3.8–5.5% apart moved by more than 1.4%". Replace the "within noise" cells in
§1 with those bounds.

**Mi6. §2.2, "the fit's cost is the same (A3: +31% macOS, +30% Linux)", and "sits a little
above".** On macOS the change measures 27.1% against A3's 31%. The unchanged mean here (298.4 s)
is 5% above A3's (284.3 s), while the change's mean (379.1 s) is only 1.6% above the libm build's
(373.2 s). Plan §8 keeps B7 because "the real conversion can inline libm and the substitution
cannot", and the report never tests that. The bench tables can: on macOS the change's cost is at or
below the libm build's in all 18 joint-fit and site-quality cases. Examples: SNP/indel with 8
samples +13.7% against +23%, site quality with 63 samples +10.8% against +21%. On Linux the two
agree within about 3 points, except six alleles (+22.6% against A3's unreliable +56%). A3's figures
came from a different session, so this is suggestive only. **Fix:** replace "the same" with the
numbers. Replace "a little above" with "5% above". Add one sentence comparing the benches with A3's
libm build and naming the inlining question as unresolved.

**Mi7. §2.3, "Every Linux figure A3 flagged as unreliable is now settled."** A3's +13% and +56%
measured the libm build. B7 measures a different build, so it replaces those figures rather than
settling them. **Fix:** "A3 flagged two Linux figures, for four and six alleles, as unreliable.
For the change they are +8.6% and +22.6%, with its two passes 0.7% and 0.0% apart."

**Mi8. §6 omits two departures from the plan.**

- A3 timed all three calling commands over 20 regions as well as 160. The interleaved session
  timed only 160 regions, so the 20-region timings exist only in `tmp/measure_B7/*_runs/`, under
  the discarded load.
- The plan asks for the parity oracle to run at B7. §3.3 cites B4 and B5 results instead. Those
  still hold at `cd056a3f`, because B6 changed `src/` only by a comment rewrap and a test-module
  `allow`, but the report should say so.

### Nits

- §2.3 Linux delimiter: "2.1% to 3.5% apart" should read 3.4%; the 3.5% belongs to `depth/412`,
  which moved only +0.5%. "had been 0 to 3.2% apart" should read 0.3 to 3.2%
  (`benches_linux/table.md`; `delim_linux/`, recomputed).
- §2.3, the lto note: "`crate::float` is inlined across the crate boundary only where link-time
  optimisation runs." `crate::float`'s functions are `#[inline]` and live in the caller's own crate.
  What needs LTO is libm's bodies, which are not `#[inline]` (`libm-0.2.16/src/math/exp.rs:85`,
  `log.rs:75`). "Overstates" is expected, not measured; write "can overstate".
- §2.3, the `log_add_exp` note: the hidden-duplication filter also runs inside
  `call-from-alignments` (`src/cli/call_from_alignments.rs`). Say whether the threshold question
  stays open for the planned performance follow-up.
- §2.3, site quality: the allele sweep runs at 1,000 samples (its 2-allele median equals the
  1,000-sample one). The +23 to +26% at six alleles is a figure at 1,000 samples.
- §3.1, "the change's records are also the same on the two platforms." So are the unchanged
  build's. All four record sets per command match: `call-from-alignments` MD5 `be3e091e`,
  `call-from-psps` `796521a6`, over 120,538 records.
- §3.2: on macOS, the change's fit differs from the unchanged build's in 57 of 574 lines
  (`tmp/ab2/macos/fit_{std,float}.toml`). Worth one clause beside Linux's 56.
- §3.2, "one step lower", and §3.4, "steps". Define once: one step is the gap to the adjacent
  `f64`.
- §3.1, "the aligner's recorded-answer test covers one of its two algorithms." Name them
  (the test covers algorithm 3; the caller ships 4u, per the B2 review) and say what a
  recorded-answer test is.
- §2.3, "They use few transcendental calls." Give A1's count or drop it.
- §4, the failing example tests: cite the log. "That step fails on `main` today" is inferred from
  arm64 runs; CI runs x86_64. Write "would fail" or cite a CI run.
- §6, "eight to ten minutes": the Linux fits took 446–593 s, 7.4 to 9.9 minutes.
- §2.1, "Linux ran first and macOS second." That holds for the commands; for the benches macOS ran
  first.
- §1's test row (4,795 and 4,796) must be updated with §4's counts when they arrive.

## Verified and correct

- **§2.2 command table and §1 speed rows.** Every range, mean and percentage matches
  `tmp/ab2/{linux,macos}/ab.tsv` and `gen_{linux,macos}/gen.tsv`. The files hold no duplicate
  round+command+build rows; the repeats are only in `run.log`'s echo. Round orders alternate as
  stated. Fastest change against slowest unchanged: fit +23.6% on macOS and +26.3% on Linux (the
  report rounds to 24% and 26%); `call-from-psps` +0.9% on macOS, 0.0% on Linux.
- **Session times and loads.** 22:04:38 to 00:57:27 for the command and `generate-psps` runs, and
  01:07 for the delimiter. Loads while macOS ran the commands and fit: 6.88–22.0, median 18.72.
  Linux median 8.36. Bench loads 1.16–5.36 on macOS and 1.56–6.58 on Linux. `generate-psps` start
  and end loads 2.29–5.07. The 1,671% CPU sample exists (23:10:42).
- **§2.3 bench tables.** Every percentage and run-to-run figure matches `table.md` on both
  platforms. The four readings in "What to read from these" hold, apart from M1: 32 samples is the
  cheapest case of each half, SNP/indel costs more on Linux and repeat tracts less in every case,
  the site-quality trends, and change passes within 1.7%.
- **Draw probe.** `b7_draw_probe.rs` uses the bench's constants, seeds, spans, depth and
  starting points. The draw hashes are identical in the two trees on macOS; the fit hashes differ.
- **Linux delimiter re-run.** 13 cases from −2.3% to +3.5%, and the depth cases −2.3, −1.9,
  +0.4, −0.5 and +0.5%. Against A3, the depth cases read 1.5–3.7% slower, and pileup and psp
  cases were within 2.0%.
- **Oracle.** The five checksums appear in `oracle_{linux,macos}.log`.
- **160-region VCFs.** 120,538 records, identical across the two builds and the two platforms for
  both commands (`--defaults`, per `ab2.sh`).
- **Parameter files.** `0fb3e50d` on all six `measure_B7` fits and on both `ab2` last-round fits.
  Unchanged `dfdefb5d` on Linux, `9d8c5539` on macOS, 7 lines apart. Unchanged against change on
  Linux, 56 of 574 lines; change against the libm build, 8 lines, all in read-group multipliers,
  inbreeding coefficients and the two concentrations. The `ln 3` probe file is `3f5e5e7f`.
- **§3.3.** 6,735 records, 30 differ: AF 9, PARALOG_POST 19, PARALOG_LR 2. No QUAL, FILTER, site
  or per-sample field differs. Every differing value is exactly one unit apart in its last printed
  place (two cases carry into the digit before). The parity records at B4 and B5 are `1e40bd1e`
  on both platforms.
- **§3.4.** The figures match the B2 and B4 commit messages.
- **Appendix A.** SHA-256 prefixes of `tmp/ab/bin/{std,float}_{macos,linux}` match, and so do
  `binaries.txt` and the `measure_B7` `build.txt` files. `c6a4394b` differs from `75722336` only
  in documentation and one example. Release profile `lto = "fat"`, `codegen-units = 1`; `profiling`
  and `soak` `lto = false`.
- **Scope, §5.** It covers x86_64, GIAB, `estimate-contamination` speed and larger or deeper
  cohorts. The B2, B3, B4 and B5 notes are each addressed (with the corrections in Mi4, M2 and the
  nits).
- **Not verified.** The claim that the unchanged macOS binary was "rebuilt during this step to the
  same checksum", the example-test failures, and the Linux-only test.

## Fixes applied

| finding | resolution |
|---|---|
| M1 noisy macOS passes, direction | "understated"; the costs against the pass near A3's figure (+16.9%, +18.1%, +11.8%, and +10.1% for the repeat-tract fit at 32 samples) are in §2.3, and §1's low end is marked as about +10% |
| M2 `estimate-contamination` | §3.1 now says both files hold no estimate (all four samples refused, 0 varying positions); §5 records output and speed as unchecked |
| M3 peak memory | §1 and §2.4 say no consistent direction, within 12% of A3's range from another session, with the 20% spread of one command's own repeats |
| Mi1 calling commands at ±1% | §2.2 gives the per-round pattern: `generate-psps` slower on macOS, faster in every Linux round; `call-from-alignments` slower in 4 of 5 Linux rounds |
| Mi2 load description | burst attributed to round 1; editor and session named; the PDF-preview service's near-continuous core stated; load from 2.4 |
| Mi3 load of 26 | replaced by the one loaded fit (496 s) against A3's quiet 280–288 s, citing `tmp/ab/*/ab.tsv` |
| Mi4 draw probe | "on macOS"; eight passes is the cap and no fit converged |
| Mi5 macOS other benches | "all 27 within 2.6%; seven beyond A3's own spread, none by more than 1.4%" |
| Mi6 against the libm build | numbers given (27.1% against 31%, 29.1% against 30%; unchanged 5% above A3's); bench comparison added, inlining question stated as unsettled |
| Mi7 "settled" | reworded as replaced, with the numbers |
| Mi8 §6 | 20-region commands not re-timed, parity oracle not re-run at `cd056a3f`, both added |
| nits | delimiter spread 2.1–3.4% and A3's 0.3–3.2%; LTO note names libm's bodies and says "expected, not measured"; filter runs in both calling commands; allele sweep at 1,000 samples; unchanged build's cross-platform identity stated; macOS 57 lines; *step* defined; the aligner's two algorithms and *recorded-answer test* explained; "few transcendental calls" dropped; 7.4–9.9 minutes; overlap order stated per phase; §1 test row matches §4. The delimiter's change stays +3.5% (the 3.4% nit applied to the passes' spread) |
| not verified: example tests | re-run on macOS on both trees; the same three fail (`tmp/ab2/example_tests_{std,float}_macos.log`); CI on x86_64 stated as inferred |
| not verified: macOS rebuild | the author's rebuild of the unchanged macOS binary gave `5252cff30207d56c`, the timed copy's checksum; recorded in this session's shell output only |
| not verified: Linux-only test | `src/psp/writer.rs:1924`, `#[cfg(target_os = "linux")]` |
