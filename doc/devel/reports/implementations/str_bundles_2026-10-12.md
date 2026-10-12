# STR bundles called as one STR locus: the four-way comparison

*Implementation report, 2026-10-12. Branch `str-bundles`, merged into main after this report. The
design is [`../../ng/spec/str_bundle_loci.md`](../../ng/spec/str_bundle_loci.md); the measurement
that chose it is [`../../ng/impl_plan/str_bundles_measurement.md`](../../ng/impl_plan/str_bundles_measurement.md)
§4; each step's review is in [`../reviews/`](../reviews/) (`str_bundles_*_2026-10-1*.md`).*

---

## 1. What changed

**An STR bundle — two or more repeats within 15 bases of each other — is now called, as one STR
locus over its hull.** Before, nothing called a bundle: on GIAB HG002's tandem-repeat benchmark
regions that left 2,422 truth variant records in 2,473 bundles uncalled.

- **Typing.** Each bundle member carries its motif. A bundle longer than `--max-str-bundle-len`
  (default 100 bases; 0 turns bundle calling off) is typed as a satellite and not called.
- **Runs.** Stutter is placed in *runs*: stretches of pure, whole copies of one motif of one to six
  bases, at least the typing's copy floor, found once from the hull's reference bases, chosen
  longest first. The typing's members overlap and nest too often (734 of HG002's 2,410 bundles) to
  carry stutter themselves.
- **Reads.** A read that crosses the hull is spelled over it, as over an STR segment. A read that
  ends inside is kept as a partial, with the flank it held recorded. psp format 2.1.
- **Calling.** Each candidate allele is parsed into the runs and fixed sequence. A read is scored
  as the allele with any set of per-run slips, each run at its own period's stutter rate and its
  own count; a partial read says nothing past where it ends. The prior multiplies the runs' seeds.
  The periodicity check is off at bundles.
- **Record.** One multi-allelic record over the hull; `RU` and `PERIOD` per run (header
  `Number=.`); the run counts per called allele in a new per-sample field, `REPCNB`.
- **STR segments are unchanged** in every benchmark run against main, apart from two changes they
  gained on the way: main's fix for partial reads held on the right flank (d02430db), and, through
  the merge, the held flank recorded for whole-tract partial reads (14 more right cells at 50
  samples, PCR, 30×).

## 2. The comparison

Four ways to call bundles, on the same reads:

- **new bundle caller** — this work;
- **our STR caller** — today's STR path with each bundle typed as one STR segment, motif of its
  longest member (an experiment build from main, not merged);
- **generic path** — bundles typed as generic sequence, coalesced with the regions either side (an
  experiment build, not merged);
- **main** — bundles uncalled.

**The stuttering benchmark**, extended to lay 500 bundles drawn to resemble HG002's (median hull 36
bases against 34; one bundle in eight holds a mono- or dinucleotide run of 16 copies or more), true
stutter, no fit. A cell is one sample at one bundle; "right" counts cells carrying a variant (about
10,435 at 50 samples, 208 at one), false calls those carrying none (about 14,565 and 292).

| bundles: right / false calls | new bundle caller | our STR caller | generic path | main |
|---|---|---|---|---|
| 50 samples, PCR, 30× | 10,407 / 1 | 10,265 / 154 | 10,071 / 756 | 398 / 0 |
| 50 samples, PCR, 5× | 7,929 / 82 | 7,479 / 241 | 7,859 / 759 | 347 / 6 |
| 50 samples, PCR-free, 30× | 10,422 / 0 | 10,341 / 11 | 10,369 / 49 | 398 / 0 |
| 50 samples, PCR-free, 5× | 8,159 / 76 | 7,770 / 122 | 8,149 / 195 | 349 / 15 |
| 1 sample, PCR, 30× | 207 / 0 | 158 / 0 | 206 / 0 | 5 / 0 |
| 1 sample, PCR, 5× | 99 / 0 | 83 / 0 | 115 / 0 | 5 / 0 |
| 1 sample, PCR-free, 30× | 206 / 0 | 160 / 0 | 208 / 0 | 5 / 0 |
| 1 sample, PCR-free, 5× | 91 / 0 | 72 / 0 | 100 / 0 | 3 / 0 |

In bundles holding a mono- or dinucleotide run of 16 copies or more (50 samples, PCR, 30×; 1,580
cells with no variant), the generic path makes 690 false calls, our STR caller 100, the new caller 1.

**GIAB HG002** (one sample, PCR-free), per benchmark region: 1,131 bundles with a truth variant, 1,342
without.

| HG002 bundles: right / wrong / missed; false calls | 30× | 50× |
|---|---|---|
| new bundle caller | 862 / 151 / 118; 20 | 895 / 136 / 100; 20 |
| our STR caller | 650 / 142 / 339; 27 | 687 / 122 / 322; 24 |
| generic path | 901 / 184 / 46; 16 | 934 / 170 / 27; 18 |
| main | 55 / 100 / 976; 1 | 56 / 104 / 971; 1 |

## 3. What it means

- **Where stutter matters — PCR libraries, cohorts — the new caller is the only one of the four that
  does not take stutter for alleles.** At 50 samples, PCR, 30×, it makes 1 false call where the
  generic path makes 756 and our STR caller 154, and it gets the most cells right.
- **Our STR caller would not have been enough.** Modelling a bundle as one repeat of one motif misses
  about a quarter of HG002's variant bundles (339 of 1,131 at 30×).
- **Where it is behind the generic path:**
  - **HG002**: 39 fewer right at 30× (862 against 901), 72 more missed, 33 fewer wrong. 42 of its 118
    misses are bundles over 100 bases, which the size cap leaves uncalled. The owner's judgement
    (2026-10-12): calling those with the generic path would bring its PCR false calls back.
  - **One sample at 5×**: 99 right against 115. An allele seen by fewer than two reads that cross
    the whole bundle is never offered as a candidate; partial reads score candidates but do not
    nominate them, as at segments.
- **Cost.** Calling the 500 bundles of a 50-sample, 30× run takes about 6.8 s of CPU on one thread
  (on a loaded host): about 14 ms a bundle, against about 4 ms before partial reads were scored and
  at most 2 ms a segment. Most of it is comparing partial reads letter by letter against every
  candidate; sharing that work between candidates is the next lever if cohort size makes it matter.

## 4. What it does not cover

- **Simulated bundles are the easy shape**: pure members separated by gaps. Real bundles whose
  repeats abut, overlap or carry impurities (about 30% of HG002's) are not simulated. Next
  (owner, 2026-10-12): extend the simulator to them (spec §10 step 5).
- **Partial reads nominating candidates** is the low-depth lever, to be measured on bundles and
  segments together at 3–5× (owner, 2026-10-12).
- **The parameter fit** was not run; every benchmark run used the simulator's true stutter.
