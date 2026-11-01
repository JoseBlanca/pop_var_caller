# The neighbour-evidence pass, phase two: solving one window

*Implementation report, 2026-11-01, branch `neighbour-pass`. Spec:
[`doc/devel/ng/spec/neighbour_evidence_pass.md`](../../ng/spec/neighbour_evidence_pass.md).*

## What exists now

The pass can re-call every sample in a window of records, given each record's `PL`, `GP` and `AF`
and each sample's inbreeding coefficient. Nothing reads a VCF or writes one yet; that is phase
three. The code is in `src/neighbour_evidence/`:

| file | what it does | spec |
|---|---|---|
| `records.rs` | one record as the pass holds it: each sample's likelihoods from `PL` and allele shares from `GP`, as `f32`; the background template's shares from `AF` | §3.2, §3.3 |
| `templates.rs` | each sample's 16 templates in a window, by agreement over at most 500 of its informative records spread over the window | §3.3 |
| `fit.rs` | each sample's `S`, `R` and `F`, fitted on the run's first windows and then frozen | §3.5 |
| `windows.rs` | where each window's core and margins end, under the memory budget | §2.2, §6 |
| `refine.rs` | one window solved: every sample's walk, the posteriors at the core's records, and the rule that withholds a call from a sample with no reads and no evidence nearby | §4.2 |
| `model.rs` | the walk (phase one); its posterior now sums the pair states as `Q̃ᵀ · W · Q̃` | §3.6 |

**The background template is checked against the Python prototype.** `hmm.py` gained
`background_shares`, and the reference file gained a fourth window that walks two templates plus
the background; the Rust walk matches it to 1e-9 like the other three. The first three windows are
byte-identical to before.

## Where the code had to be more precise than the spec

Each is now written into the spec.

- **A record where a sample has no reads is walked, not skipped.** §7 said skipping it was exact.
  It is not: the transition across two gaps differs from the transition across their sum, because a
  state that changes class and changes back redraws its templates, while one long step counts them
  as kept. On the prototype, removing one such record from a nine-record window moved the other
  records' posteriors by up to 5 in 10,000. So every record with at least two samples' rows is a
  site of every sample's walk.
- **How several fitting windows combine.** Each window adds its twelve grid walks' log-likelihoods
  to a running total per grid point; the largest total is walked with; `F` is re-estimated from that
  point's walk over every site the fitting windows walked, and freezes with `S` and `R`.
- **Which 500 records templates are chosen over**: spread evenly over the window, not the first
  500, which would sit in the left margin.
- **How the budget shrinks a window**: the core first, down to the margins' length, then both
  together to the 500-kb floor. At 3,000 samples and tomato's density that gives the 1.5-Mb window
  §6 quotes.
- **The no-read rule at its edges**: a record at the same position does not count as evidence on
  either side, and at a contig's right end the one-sided allowance needs the header's contig length.

The forward vectors stay `f64`: at 16 templates and the background a walk keeps 306 numbers a
record, about 98 MB a thread over a 40,000-record window. Halving that is a change for when memory
is measured.

## A first look on the founder simulator

Twenty samples from 4 founders, inbred (95% of each genome in autozygous runs), 4 Mb in two contigs,
PCR-free libraries; mapped with minimap2 and called by pvc with fitted parameters; then the windows
solved by `examples/ng_neighbour_pass_probe.rs` with 600-kb cores and 400-kb margins, so each contig
holds four windows. Scored by `doc/devel/ng/research/neighbour_evidence/score_founder_probe.py` at
SNP records whose position and alleles match the truth's. One unit is one sample at one record.

| depth | sample-sites | where the sample has reads: wrong, calling → pass | where it has none: calling writes `./.`; the pass calls, and is wrong at |
|---|---|---|---|
| 3× | 29,020 | true non-reference: 282 → 15 of 13,938; true reference: 5 → 6 of 13,534 | 1,548, wrong at 23 |
| 1× | 10,960 | true non-reference: 73 → 80 of 4,825; true reference: 1 → 13 of 2,025 | 4,110, wrong at 463; 289 of them true reference called non-reference |

- **The fit finds the simulation's settings.** At 3× every sample chose `S` = 300 kb, the simulated
  value; `R` came out 1 or 3 Mb against 2 Mb simulated, and `F` 0.86 to 0.99 against 0.95.
- **At 3× the pass removes 19 in 20 errors at true non-reference genotypes and calls every site
  without reads**, wrong at about 1 in 70 of them.
- **At 1× it adds errors where samples have reads**, 19 more over 6,850, and its calls at sites
  without reads are wrong at 1 in 9; a fifth of the true-reference ones become non-reference. pvc
  itself wrote 1,063 records there against 3,585 at 3×, so each sample had a few hundred informative
  records over the whole 4 Mb to choose templates by, and no sample's fit reached the 5,000 records
  that freeze it. This is the whole-cohort-at-1× case the spec already leaves to measurement (§11 Q8);
  phase four's 1× cells decide what to do about it.
- **Not scored**: pvc SNP records with no truth record at that position — 206 records at 3×, where
  calling and the pass are both wrong about 15 times in 100. Some are false records, some real SNPs
  the truth writes inside a repeat's allele; an exact-match scorer cannot tell which.

These are one seed and one cohort shape. They show the code runs end to end on mapped reads and is
not wrong in an obvious way; they are not the measurement §8 asks for.

## Tests

25 in `src/neighbour_evidence/`, run with `cargo test --lib neighbour_evidence`. Besides the 1e-9
reference, they check template choice (agreement, ties by sample order, absent samples as equal
shares, the 500 spread over the window), the fit (clamping, freezing, `F` high for an inbred sample
and low for an outbred one), the windows (every record in exactly one core and in order at 1 to
10,000 samples, the 1.5-Mb window at 3,000, the budget exceeded and reported past the floor, no plan
before the right margin is read), and the solve (errors fewer and no-read sites filled on a drawn
inbred cohort at 1×, a sample that stops having reads withheld, only the core returned, the same
bits at one thread and four).
