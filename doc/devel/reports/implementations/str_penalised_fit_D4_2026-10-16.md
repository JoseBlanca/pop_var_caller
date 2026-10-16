# D4: the oracle and the cross-platform digests re-recorded (2026-10-16)

Plan step D4 of [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md) (spec §4 item 6:
"the oracle moves, and the move is explained before its checksums are re-recorded"). Branch
`str-penalised-fit`, binary built from `932fc9b0` (the next commit, `2b3416bd`, is documents only).
Measured on macOS, arm64.

Two sets of checksums pin this caller's output, and both had gone stale on this branch:

- **the identity oracle** (`scripts/promote_ng_oracle.sh`): four tomato2 accessions over 20 regions,
  called with default parameters, fitted, and called again with the fitted file. Its checksums sit
  in `scripts/promote_ng_oracle.baseline`, last re-recorded at plan step B4;
- **the cross-platform digests** (`src/cli/cross_platform_digests.rs`): three small synthetic
  cohorts (2 samples, 2 samples with sequencing errors, 21 samples), each walked, fitted and called
  in a unit test. Last recorded on 2026-10-12 (`e7c9106e`).

## 1. The answer

**Every move is explained by a commit, and both sets are re-recorded.** The new fit moves the fitted
parameters file and the calls made with it; nothing called with default parameters moves, beyond the
notes in the defaults run's parameters file. One move was not described in the commit that caused it:
**the stutter prior (C8) left six thin three-base repeat classes without an F of their own** (§2.3).
It is the existing rule for an F the data do not place, not a defect, but the owner should know of it.

The full suite passes except the three tests that check a panic only debug builds raise, as expected
in release: 5,362 passed, 3 failed, 12 ignored (`cargo test --release --features bench-fixtures
--no-fail-fast`).

## 2. The oracle

### 2.1 What moved

| checksum | B4 baseline | now | moved by |
|---|---|---|---|
| default-parameter VCFs (both routes), windows, histograms | — | unchanged | — |
| `from_alignments.parameters.toml` (the defaults run's file) | `3d1719c2…` | `4f272d7e…` | notes only: C8, C9, `7265efd0` |
| `fitted.parameters.toml` | `419869a0…` | `5ed794a0…` | C7, `0b12941d`, C8, C9, C10, `7265efd0` |
| `from_fit.comparable` (calls made with the fitted file) | `e4af91bd…` | `61b09b60…` | C7, C8, C9, C10 |

Each change was bracketed by its own oracle run when it was committed (`tmp/oracle_C7`, `oracle_Cfix`,
`oracle_C8`, `oracle_C9`, `oracle_C10`), so each difference below is pinned to one commit by
comparing neighbouring runs. Today's run (`tmp/oracle_D4/`) gives the same calls as C10's, byte for
byte, and a fitted file that differs from C10's only in one note (§2.2). The performance commits after
C10 move nothing, as their end-to-end check found (`tmp/e2e/`, the same `5ed794a0…` before and after).

### 2.2 The two parameters files

**The defaults run's file** moves only in comment lines; every line that is not a comment is
identical to B4's. C8 added a paragraph on the stutter prior, C9 a sentence saying the stated F is a
median over strata of 8 or more loci, and `7265efd0` reworded the prior's paragraph (the centre is
printed 0.1, not 0.10, and the sentence on rows below one expected slipped read now says those rows
carry no error on any of their three numbers).

**The fitted file**, B4 against now (345 → 391 lines that are not comments):

- **Repeat-tract rows** (C7, then `0b12941d`, C8, C9). Slippage rows 37 → 49 and allele-length rows
  15 → 49, 40 of them with an F. On the 37 classes both fits wrote, the slip rate now is a median 0.88
  of the climb's (0.24 to 1.60). The climb's flat tail of borrowed values is gone: one-base classes of
  16 to 25 repeats held 0.184 each and now follow a curve, from 0.106 at 16 repeats down to 0.045 at
  20 and 21 and up to 0.123 at 25. The smallest slip rate in the file is 0.0017; before C8's prior it
  was 4.2 × 10⁻¹⁸. The fallback F moved 0.526 → 0.723 (C7) → 0.710 (C8) → 0.6135 (C9) → 0.6131 (C10).
- **The order of the cohort** (C10). Samples and read groups are listed in name order, so read groups
  are renumbered, and each library's numbers move with its renumbered row.
- **Every SNP/indel number in its last digits** (C10): the same psps are now summed in a different
  order, and that fit stops at its pass limit, so the order shows. Each sample's F moves by at most
  0.0003 (SRS1839219 0.96645 → 0.96617); each library's error multiplier by at most 4 parts in 100,000
  (2.67999 → 2.67990); the ordinary-site prior's reference concentration 1.8699 → 1.8730; the three
  fit starts' pass counts by 1 or 2.
- **Notes** (C8, C9 and `7265efd0`, as in the defaults run's file).

The substitution rates, census digests and every other table did not move.

**`7265efd0` and the fall-off's error.** That commit removes the fall-off's standard error from any
row with fewer than one expected slipped read but at least one repeat slipped. No oracle row is such:
between C10's file and today's, only the one note differs.

### 2.3 Six three-base classes lost their own F at C8

At C7 46 of the 49 allele-length rows carried an F; from C8 on, 40 do. The six that lost it are all
three-base classes, with 1 to 5 loci each (13 loci together), at 7, 8, 9, 10, 11 and 17 repeats; at
C7 their F ran from 0.54 to 0.83.

What happened, from the two runs' logs (`tmp/oracle_Cfix/fit.log`, `tmp/oracle_C8/fit.log`): the
stutter prior raised those classes' slip rates from almost nothing (5 × 10⁻⁶ to 0.0017, and 8 × 10⁻¹⁰
at 17 repeats) to 0.013–0.026, because a class with no slipped read now leans on one read's worth of
calling's default stutter. With stutter now available to explain small length differences, the
three-base F curve went from stiff to loose (its bend penalty 1.2 × 10³ → 3.4 × 10⁻³), and in each of
the six classes the F estimate ran to within 1% of an end of its search. The rule built in milestone A
then writes no F (the run log says "F absent: its estimate stopped within 1% of an end of its
search"), and calling uses the fallback F (0.613) at those tracts.

This is the designed behaviour for an F the data cannot place, and it affects classes of at most 5
loci. C8's commit message described its oracle move as "repeat-tract rows only", which is true, but
did not mention these six. I record it here rather than stopping, because its cause is certain (C8 is
the only code change between the two runs) and the mechanism is in the log; **if the owner wants a
class with a few loci to keep an F at an end of the search rather than fall back, that is a design
question for spec §2.10, not something D4 changes.**

### 2.4 The calls made with the fitted file

B4 against now, records keyed by chromosome, position and reference allele:

- **6,683 → 6,680 records.** Gone: a one-base deletion in an A homopolymer at SL4.0ch01:11857915
  (QUAL 22.1, a 0/1 in one sample at GQ 3) at C7; two SNVs 5 bases apart at SL4.0ch01:85900664 and
  85900669 (QUAL 9.3, hidden-duplication posterior 0.51–0.52, at that filter's edge) at C10. One
  record gained an allele at C7: the repeat tract at SL4.0ch01:80960949 became multi-allelic.
- **Three records change a sample's genotype**, all repeat tracts, all at C7:
  SL4.0ch01:79435466 (one sample ./. → 1/1), 80960949 (one sample 0/0 → 2/2, the new allele) and
  81048133 (one sample 1/1 → 0/0, QUAL 53.8 → 34.9).
- **QUAL moves in 747 records**, from −19.8 to +21.8, median −0.1. 426 of them are repeat tracts,
  moved by the new stutter, shares and F (C7: 427 records, −20.9 to +29.6; C8: 326, −7.8 to +4.3; C9:
  2). The other 321, and most of the small moves elsewhere, are C10's: with every SNP/indel number
  moved in its last digits, AF moves in 5,852 records by a median 5 × 10⁻⁶ (at most 0.005), the
  hidden-duplication score in 5,840 by a median 0.0003, sample GQs by 1 or 2, and QUAL by −0.2 to
  +0.3.
- **The hidden-duplication filter's prior** (the header line `##paralogFilter`, its estimate of how
  common duplicated sites are) moved 0.050935 → 0.050928 at C7 and → 0.050664 at C10. The filter fits
  that prior over all records, so a change at the repeat tracts moves every record's posterior
  slightly: at C7, 675 SNP/indel records changed their posterior alone, by at most 3.4 × 10⁻⁵.

### 2.5 The new baseline lines

```
4f272d7ed368a743eb2295c44cc5749c  from_alignments.parameters.toml
5ed794a02799828688b98edbd7661be7  fitted.parameters.toml
61b09b604ee5fb2d7196150824c2c996  from_fit.comparable
```

The other five lines are unchanged. The baseline file carries a paragraph naming these moves; checked
against `tmp/oracle_D4/digests.txt` with the script's own comparison, every line matches.

## 3. The cross-platform digests

The three tests' fitted files and VCFs were read out of the failure messages the suite printed at
each stage (`tmp/test_merge.log`, `suite_milestone_B.log`, `suite_milestone_C.log` for C7,
`suite_C_final.log` for C10, and today's run) and checked against the checksums printed beside them.
The state they were recorded at was rebuilt from `e7c9106e` and reproduces all six recorded checksums
on macOS, so the recorded files could be compared too (`tmp/d4/digests/`).

**What moved, in order:**

1. **The rename** "length spectrum" → "allele length variation" (`aa42d629`): two keys and four note
   paragraphs. Nothing else. The allele window at ±9 (`38b1fe56`) arrived at the same time and moved
   nothing here, because the climb wrote no repeat-tract row for these fixtures; spec §3 names both as
   the cause, and only the rename is.
2. **Format 3** (B2, B4): `format_version = 3`, the fallback written as F 0.5 `defaulted` in place of
   the concentration 1.0, and notes.
3. **The penalised fit** (C7) writes what the climb did not: one slippage row and one allele-length
   row for the fixtures' two-base tract at 10 repeats, and the fallback F becomes that class's F,
   `fitted_here` (0.254 for the two-sample cohorts, 0.897 for 21 samples). The calls follow, at the
   repeat tract only.
4. **The stutter prior** (C8) moves that row's slip rate from 5.8 × 10⁻⁷ to 0.0042 (two samples) and
   from 5.4 × 10⁻⁸ to 0.00032 (21 samples), and its fall-off from 0.29 to 0.049 (two samples) and from 0.47 to 0.050 (21
   samples).
5. **Notes** (C9, `7265efd0`). C10's ordering moves nothing: the samples were already in name order.

No number outside the repeat-tract rows moved, and in the VCFs only the repeat tract at chrV:201:

| fixture | QUAL | AF | GQ, the heterozygote | GQ, the others |
|---|---|---|---|---|
| two samples | 36.7 → 90.1 | 0.249902 → 0.250006 | 34 → 75 | 45 → 47 |
| with sequencing errors | 37.4 → 90.4 | 0.249924 → 0.250010 | 34 → 74 | 42 → 44 |
| 21 samples | 37.3 → 78.6 | 0.023718 → 0.023811 | 24 → 65 | each up 2 or 3 |

No genotype moved, and the SNPs did not move.

| constant | was | now |
|---|---|---|
| `FITTED_PARAMETERS_MD5` | `a5b4997b…` | `b673dc9611c7126b5c35e9d7629493b3` |
| `CALLS_MD5` | `b1a1a907…` | `dae7e77c280d7f240349cc8433207a81` |
| `WITH_ERRORS_FITTED_PARAMETERS_MD5` | `4859eacd…` | `3a016fb7a40bc91eb5bb9ae2adc6e286` |
| `WITH_ERRORS_CALLS_MD5` | `f65f586f…` | `d315b4e4d0de2f21fd34cd673bbdfa8b` |
| `OF_21_SAMPLES_FITTED_PARAMETERS_MD5` | `0573095b…` | `5abcd04153d5acc32463fcd751970024` |
| `OF_21_SAMPLES_CALLS_MD5` | `5c2414e4…` | `b23487bd69f9d39f1943fbe9350a7a3f` |

Recorded on macOS (arm64) only; the Linux container was not run. All four tests in the module pass,
the pool-width test included (1, 4 and 7 threads give the same bytes).

## 4. Where the evidence is

- Oracle run: `tmp/oracle_D4/`, log `tmp/oracle_D4.log`.
- Fixture files at each stage and the scripts that compared them: `tmp/d4/` (`extract_digests.py`,
  `section_diff.py`, `vcf_moves.py`, `digests/`); the rebuilt `e7c9106e` tree in
  `tmp/d4/src_e7c9106e/` and its build in `tmp/d4/target_e7c9106e/`.
- Suite log: `tmp/d4/suite_D4.log`.
