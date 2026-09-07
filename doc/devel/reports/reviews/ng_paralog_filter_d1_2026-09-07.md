# Code Review: ng_paralog_filter_d1

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), two sub-agents in isolated worktrees
**Scope:** step D1 of the hidden-duplication filter plan — the filter run at one sample, and the two
tools that measure it
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the uncommitted working tree on branch `ng-paralog-filter` at `5e7f20da` —
  [`scripts/ng_paralog_filter_runs.sh`](../../../../scripts/ng_paralog_filter_runs.sh),
  [`scripts/ng_paralog_filter_relations.py`](../../../../scripts/ng_paralog_filter_relations.py),
  and the step's implementation report.
- **Out of scope:** everything under `src/` (D1 was a measurement step when the review was
  dispatched; the report line it later gained is covered under *Findings applied after the review*
  below), production, and steps A1–C5.
- **Categories dispatched:** two, and neither is one of the skill's usual per-category checklists,
  because this step's deliverable is **numbers plus two measuring tools** rather than library code:
  1. **adversarial re-measurement** — re-run the step independently and check every figure in the
     report against the reviewer's own run;
  2. **the tools as code** — falsify them with deliberately broken inputs.
- **Both sub-agents worked in their own detached worktrees at `5e7f20da`**, using the release binary
  built from that commit rather than rebuilding.

### 2. Verdict

**Request-changes.** Zero Blockers, nine Major, twelve Minor.

**The counted figures were all right and the explanations were nearly all wrong.** The
re-measurement reviewer confirmed twenty-one separate quantities exactly — record counts, the whole
calibration line, the flagged record field by field, six ratio percentiles, the window depths, the
Pearson correlation, the record-kind and genotype counts, the six-accession oracle's hash — and then
found **five Major defects, every one of them a claim about *why*.** This is `CLAUDE.md`'s and the
reporting skill's named failure mode arriving in its purest form: not a wrong number, but a wrong
mechanism attached to a right number.

**And the two measuring tools could not see seven of the sixteen failures they were shown.** The
tools reviewer built sixteen broken inputs; nine were caught. Four of the seven misses were the kind
that would let a real defect through silently.

### 3. Execution status

| command | with this step | at `5e7f20da` |
|---|---|---|
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6639 passed; 0 failed; 15 ignored` | `6637 passed; 0 failed; 15 ignored` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, three kinds | same 9, same kinds |
| `cargo fmt --check` | dirty on 9 files, none this step's | the same 9 |

One integration test fails at both — `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`,
which is `main`'s. `--all-targets` is unavailable (`main` is red); `cargo doc` and `cargo audit` not
run.

**Clippy was counted by error kind, not by total**, per the failure this plan logged on 2026-09-07:
`grep -E "^error: " | grep -v "could not compile" | sort | uniq -c`, which can see a kind the
baseline does not have. The three kinds are `can be more succinctly written as a byte str` (3),
`the following explicit lifetimes could be elided` (3) and `this operation has no effect` (3). `fmt`
was counted by **unique file** — nine — and not by the twenty-eight `Diff in …:line:` lines, which is
the same error in a different measure.

**Findings labelled "Needs verification": 0.** Every finding came from a re-run or a broken input in
a reviewer's own worktree.

**Falsification totals: 16 broken inputs built, 9 caught, 7 missed.**

### 4. Open questions and assumptions

1. **Should the run print what the fit came to?** Raised by findings M3, M4 and Mi8 together: the
   report's central mechanism could only be recovered by inverting the run's own ratios, because
   nothing printed the fitted one-copy depth or σ₀. **Answered during synthesis and applied** — see
   *Findings applied after the review*. Plan step D1's own words are "the fit's outcome", and the run
   was reporting only the fit's acceptance count.
2. **`--paralog-fdr 1e-12` writes `target_fdr=0.0000` in the header**, which a reader would take for
   the off run's setting. Not this step's code
   ([`vcf/header.rs`](../../../../src/ng/vcf/header.rs), C4's), and not reachable at any sane target.
   Filed as an out-of-scope observation.

### 5. Top 3 priorities

1. **M1** — the report says two records are 0.14 apart on the ratio axis. They are 13.46 apart.
2. **M4** — "2.24 times one copy's depth" is not the model's one copy; the model reads that window
   as about 3.5 copies.
3. **M6–M9** — four ways the measuring tools pass on a broken file.

### 6. Findings

#### Major — the report's mechanism claims

**M1: `ng_paralog_filter_d1_2026-09-07.md` — "0.14 apart on the ratio axis" is the wrong quantity**

**Confidence:** High. **Categories:** re-measurement.

The report said the cut separated two records "0.14 apart in ratio". The two ratios are 21.5702 and
8.1116 — **13.4586 apart**, a factor of about e¹³·⁵. 0.1384 is the distance from the *kept* record to
the *cut*, which is a different statement. As written it told the reader the filter split two records
whose scores nearly agree.

**Why it matters:** it is the report's headline illustration of soft calibration, and it inverts what
the run shows.

**Suggested fix:** state both quantities and name each.

---

**M2: the read counts were left out of the explanation, and they are a fifth of the gap**

**Confidence:** High. **Categories:** re-measurement.

The report explained the two records' difference as GC alone. Its own table two paragraphs above
shows the other difference: `AD 0,15` against `AD 0,19`. Both are fully hom-alt, so under the
real-variant story each read costs `ln(1 − ε) = −0.01005` at `ε = 0.01`
([`model_params.rs`](../../../../src/ng/paralog/model_params.rs) `DEFAULT_PSEUDOCOUNT_VAF`), and
under the duplication story's best branch `ln 0.5 = −0.6931`. Four extra reads cost the kept record
`4 × 0.6831 = 2.73` nats.

**Checked independently during synthesis**: the seven enumerated configurations are
`(T,m) = (3,1),(4,1),(4,2),(6,1),(6,3),(8,1),(8,4)`
([`locus_score.rs`](../../../../src/ng/paralog/locus_score.rs) `enumerate_carrier_configs`, pinned by
`default_carrier_configs_are_the_seven_expected`), so `m/T` never exceeds ½ and no better branch
exists. The sub-agent's arithmetic holds.

---

**M3: the GC mechanism was asserted with no size**

**Confidence:** High. **Categories:** re-measurement.

The report said "the model's expected one-copy depth is not the same number at the two places" and
stopped, so nothing told the reader whether a GC step of 0.086 can move a score by 13 nats. The
reviewer re-implemented the scorer and inverted all 217 ratios: expected one-copy depth 5.51 reads at
GC 0.453 against 6.33 at GC 0.367, copy numbers 3.45 and 3.02, worth 10.7 nats — and the
counterfactual the report gestured at, that **the kept record would have scored about 19.4 and been
flagged had it sat at the flagged record's GC**.

**Corroborated by the run**, after the change below: the reconstruction's σ₀ of 0.300 is what the run
now prints, to three decimals.

---

**M4: "2.24 times one copy's depth" is a model quantity the number is not**

**Confidence:** High. **Categories:** re-measurement.

The report divided the flagged record's window depth by the **median depth of the 217 records the run
wrote**, and called the result a multiple of one copy. Variant sites are where excess coverage
collects, so that median (8.49) is a biased subset and not the fit's one-copy level. The run now
prints the level: **5.22 reads**. The model reads that window as about 3.5 copies, near its winsor
cap of four — not 2.24.

**Why it matters:** it is the sentence a reader carries away, and it was wrong by a factor of 1.6 in
the direction that makes the flagged record look less duplicated than the model thinks it is.

---

**M5: the inbreeding coefficient cannot move the ranking at one sample**

**Confidence:** High. **Categories:** re-measurement.

The report attributed part of the ranking to "the allele split and the inbreeding coefficient".
`ParalogScorePrecompute` is built once per pass from a cohort-length slice of coefficients, so `F` is
locus-invariant by construction; at one sample it is one number applied identically to all 217
records and cannot move any record relative to any other. The allele-split half is right.

#### Major — the measuring tools

**M6: `ng_paralog_filter_relations.py` — the relations cannot see the filter losing its own INFO keys and header lines**

**Confidence:** High. **Categories:** tools.

Both relations strip `PARALOG_LR`, `PARALOG_POST` and the four `##` lines before comparing, so a
filter that stopped writing them satisfied both. Demonstrated: with both keys stripped from every
record of the dropping and tagging runs and three of the four header lines deleted, the tool reported
both relations HOLD and exited 0. The one count that would have shown it was printed and never
compared against anything.

---

**M7: a flagged record could be corrupted in any column and both relations still held**

**Confidence:** High. **Categories:** tools.

Relation one removes the flagged records from the off side; relation two removes them from the tag
side. Between them nothing compared a flagged record's QUAL, its remaining INFO, its FORMAT or any
sample column. Demonstrated: replacing the flagged record's sample column with `0/0:1:1:1,0` and its
QUAL with `0.1` left both relations holding, exit 0.

**Why it matters:** the tagging run's flagged records are the only place a ratio survives, and every
per-record number in D1's report is read out of them.

---

**M8: records were matched by `(CHROM, POS, REF, ALT)`, which is not an identity**

**Confidence:** High. **Categories:** tools.

Spec §3.2 states that a repeat tract may share a position with the generic locus owning its anchor
base. Two records sharing all four columns make relation one wrong in **both** directions — a false
failure when one twin is flagged and the filter correctly keeps the other, and a false pass when the
filter wrongly removes both. **Latent on this data**: over 8,510 records from six accessions there is
no repeated `(CHROM, POS)` at all, which is why this is Major and not a Blocker.

---

**M9: `ng_paralog_filter_runs.sh` — the run report section can be silently empty**

**Confidence:** High. **Categories:** tools.

`sed -n '/^hidden-duplication filter:/,$p'` prints nothing and exits 0 when the marker is absent, so
rewording the run report's first line would print an empty section and exit 0, taking π, the cut, the
convergence flag and the drop count out of every report written from that output while the record
counts above still looked healthy.

#### Minor

- **Mi1** — the percentile convention is unstated, and three of the eleven quoted percentiles change
  at the printed precision under linear interpolation.
- **Mi2** — "0.456 linear correlation" is over-read: `r² = 0.208`, Spearman 0.364, and the two real
  causes of the shortfall are the ratio's convexity in copy number and the read count pushing the
  other way — not `F` (M5). GC is *not* the missing piece: restricted to the flat part of the GC
  curve the correlation is 0.434, essentially unchanged.
- **Mi3** — the window-depth table cannot be reproduced from the report: it needs
  `NG_WINDOW_COVERAGE_FILE` set on the run, the rows are `f32` bit patterns written as decimal `u32`,
  and it must be the tagging run's dump. None of that was stated.
- **Mi4** — "deeper than three reads a position" compares a variant-site median with a per-position
  figure.
- **Mi5** — peak resident quoted to 0.1 MB from a single run.
- **Mi6** — `mark before this run` is structurally always zero, and its comment described a guard
  that was not doing anything.
- **Mi7** — `ru_maxrss` is kilobytes on Linux and bytes on macOS; a host run would print a figure a
  thousand times too large.
- **Mi8** — `drop.comparable` and `tag.comparable` were written and never read.
- **Mi9** — the harness never said which binary produced its numbers, though it chooses between two
  target directories by modification time.
- **Mi10** — the `hiddenParalog` count matched a substring of column seven rather than a whole id.
- **Mi11** — a run that wrote zero records reported three clean zeros and exited 0.
- **Mi12** — on a length mismatch the difference report printed five copies of the same one-record
  shift.

#### Out of scope observations

- **`--paralog-fdr 1e-12` prints `target_fdr=0.0000`** in the header
  ([`vcf/header.rs`](../../../../src/ng/vcf/header.rs)), which reads as the off run's setting. C4's
  code, not reachable at any sane target.
- **Nine `cargo fmt` files and nine clippy errors are `main`'s**, unchanged by this step.

### 7. What the tools survived

Recorded because a tool that correctly fails on a broken input is a result worth keeping. Of the
sixteen broken inputs, these nine were caught by the version under review:

| broken input | verdict |
|---|---|
| one genotype byte changed in a kept record of the dropping run | both relations fail, exit 1 |
| one QUAL changed in a non-flagged record of the tagging run | relation two fails, exit 1 |
| a non-flagged record deleted from the dropping run | both fail, the length mismatch stated |
| a spurious `PARALOG_LR` added to the off file | relation one fails, exit 1 |
| a record whose INFO is exactly the two filter keys, against `.` in the off run | holds — the round trip gives back `.`, which is what the writer emits |
| an unrelated key whose name starts with `PARALOG_` | holds — whole names are compared |
| a flagged record carrying `EMNoConv;hiddenParalog` | recognised as flagged; the id test splits on `;` |
| the dropping and tagging files swapped on the command line | both fail, exit 1 |
| the harness pointed at a nonexistent CRAM | prints the run's error and exits 1 |

Three further harness behaviours were checked and are sound: `##parametersFile=` really carries only
the basename and really is the whole difference a rename makes; a missing parameters file aborts the
harness through `set -e`; and the container has `bash`, `python3`, `shasum` and `sha256sum` but no
`/usr/bin/time`, with `/bin/sh` as dash, exactly as the script's header says.

### 8. Findings applied after the review

Recorded here rather than only in the fix-application report, because one of them adds library code
that this review did not see:

**The run now prints what each sample's coverage fit came to** —
`ParalogScoringContext::what_each_fit_came_to` and a report line naming every fitted sample up to
ten, giving the spread past that. It closes M3, M4 and Mi3 at source rather than in prose: the
quantity the report got wrong is now a number the run states. **Two tests, seven mutations, all
seven killed.** The change adds no behaviour and moves no VCF byte — the standing oracle still gives
2,311 records at sha256 `84ad19c2…`.
