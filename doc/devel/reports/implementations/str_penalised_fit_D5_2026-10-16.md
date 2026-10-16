# D5: calling with the new fit's parameters against the climb's, HG002 and the GIAB trio (2026-10-16)

Plan step D5 of [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md) (spec §4 item 7).
For the owner's decision at checkpoint D: whether calling with parameters from the new repeat-tract
fit is at least as good as calling with parameters from the fit it replaces (the climb), before the
climb is deleted.

## 1. The answer

**The two fits' parameters call the same to within a handful of tracts and one indel.** On HG002's
tandem-repeat truth, the new fit's parameters get 1 and 3 more homopolymer tracts right than the
climb's (30× and 50×, of about 3,850 called), and in the period-2-and-longer class 10 fewer at 30×
and 1 more at 50× (of about 2,850). On the GIAB trio, SNP calls are identical and the new fit finds one true indel
fewer at each depth (303 against 304 of 330 at 30×, 312 against 313 at 50×), with the same false
positives.

**The new fit is also about 60 times faster on one sample:** `estimate-parameters` on HG002's
tandem-repeat ground at 30× took 46 s with the new fit and 2,927 s (49 minutes) with the climb, the
rest of the command being the same code.

### HG002, tandem-repeat tier: right genotypes of the tracts called

Scored by `benchmarks/lib/tract_qual_experiment.py` against GIAB's tandem-repeat truth, as in
`str_partial_read_phase_2026-10-12.md` §4.

| depth | class | new fit | climb | compiled defaults (context) |
|---|---|---|---|---|
| 30× | homopolymer | 3,462 of 3,835 (0.9032) | 3,461 of 3,839 (0.9020) | 3,515 of 3,869 (0.9092) |
| 30× | period 2 and longer | 2,556 of 2,854 (0.8965) | 2,566 of 2,856 (0.8994) | 2,573 of 2,871 (0.8971) |
| 50× | homopolymer | 3,568 of 3,884 (0.9194) | 3,565 of 3,885 (0.9183) | 3,598 of 3,905 (0.9223) |
| 50× | period 2 and longer | 2,609 of 2,877 (0.9072) | 2,608 of 2,880 (0.9059) | 2,619 of 2,889 (0.9075) |

Of 6,975 tracts, 46 change verdict from the climb's parameters to the new fit's at 30× and 27 at
50×. At 30× the largest group is 14 homopolymers that go from a spurious heterozygote to right; the
losses are spread (8 period-2-and-longer tracts whose true allele is no longer offered, 5 and 4
homopolymers with no record or a collapsed heterozygote). Full list: `tmp/d5/out/tables.txt`.

### GIAB trio: recall and precision, pooled over HG002, HG003 and HG004 (QUAL ≥ 30)

| depth | type | new fit | climb | compiled defaults (context) |
|---|---|---|---|---|
| 30× | indels | 303 found, 2 false, 27 missed | 304, 2, 26 | 306, 2, 24 |
| 30× | SNPs | 2,004, 3, 57 | 2,004, 3, 57 | 2,020, 6, 41 |
| 50× | indels | 312, 1, 18 | 313, 1, 17 | 313, 1, 17 |
| 50× | SNPs | 2,023, 3, 38 | 2,023, 3, 38 | 2,031, 6, 30 |

The one indel lost at 30× is HG003's; at 50×, HG004's. Per-sample rows: `tmp/d5/out/tables.txt`.

## 2. What the comparison held fixed

- **The same psps** for both fits (written once by the new binary), and **one calling binary**
  (HEAD `8c6ce6f5`) calling with `call-from-alignments --parameters` in both arms.
- **One fit per sample and depth**: a calling run refuses a parameters file whose samples are not
  its own, and both benchmarks call one sample at a time, so every fit is a one-sample fit, the
  hardest case for the fit. **Six fits per arm, not the eight intended:** the trio's HG002 fit and
  the tandem-repeat tier's HG002 fit were written to `HG002_<depth>` and `hg002_<depth>`, which
  macOS's case-insensitive filesystem treats as one file, so the trio's HG002 was called with the
  fit made on the tier's ground (same sample, wider ground). Both arms did the same, so the
  comparison between them is unaffected; the trio's HG002 rows are not a fit of its own 100
  regions.
- **Nothing but the repeat-tract numbers differs** between the arms' files: with comments and every
  `[repeat_tracts…]` table removed, all six pairs of files are identical.
- The climb's parameters came from `pvc` built at `9c57b942` (the last commit before the switch),
  which already writes parameters-file format 3; the calling binary reads it unchanged (checked on the
  C6 oracle: the same VCF record for record).

## 3. Context the table shows but this step does not judge

The compiled defaults (no fit) call more SNPs right than either fitted arm on the trio: 2,020
against 2,004 at 30×, with 6 false positives against 3. That difference comes from the SNP/indel
fit at one sample, which both arms share and this branch does not change; it is not a property of
either repeat-tract fit.

## 4. How it was run

`tmp/d5/run_d5.sh all` (the steps and their commands are in `tmp/d5/README.md`), every step under the
memory guard. Wall: about 2 hours, almost all of it the climb's two HG002 fits. Largest peak of any
step: 720 MiB.
