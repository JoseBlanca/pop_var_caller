# psp record head — H3: what the two new head fields cost, measured

**Date:** 2026-09-04
**Plan step:** [psp_head_compared_reads.md](../../ng/impl_plan/psp_head_compared_reads.md) Milestone H, step H3
**Spec:** [psp_head_compared_reads.md](../../ng/spec/psp_head_compared_reads.md) §4
**Branch:** `ng-psp-mode`

## The answer

**The change costs 3.9 % of the compressed file at 10.25 reads a position and 8.5 % at 280.32.**
In raw bytes it is exactly what §4 predicted — one byte a record at low depth, two at high — and
the kind's move contributes none of it. What §4 did not predict is how little compression removes
at depth: **81 % of the added bytes disappear at ten reads a position and only 30 % at 280**,
because the compared-read count tracks depth and varies record to record, where the non-reference
count is almost always zero and a column of zeros costs nothing.

**So the field is cheapest where the flat floor already answered the keep rule, and dearest exactly
where it is needed.** Whether 8.5 % of a deep file is worth the body decodes it saves is a question
for the reader that consumes it, which is the successor plan's and does not exist yet.

Measured on two production psps rebuilt with the current `pileup`, at the 100 kb default grid and a
32 kB window, compressed bytes a record. Every row is the same records and the same bodies, with the
head written a different way:

| what the file carries | tomato SRR7279481, 10.25 reads a position | HG002 chr21, 280.32 |
|---|---:|---:|
| bodies alone, no head scalars | 4.752 | 14.870 |
| + the head as it stood before this milestone | 4.944 (+4.05 %) | 15.831 (+6.46 %) |
| **+ the head as it stands now** | **5.138 (+8.12 %)** | **17.178 (+15.52 %)** |
| the change itself | **+3.91 %** | **+8.51 %** |
| raw bytes a record added | +1.00 | +1.93 |
| the head compressed on its own | 0.253 → **0.509** | 2.425 → **3.945** |

7,687,686 records and 56,627 records. The command is
`./target-container/release/examples/ng_psp_head_encoding tmp/h3/<psp> --label <name>`, and the two
runs as they printed:

```
sample	tomato				mean-reads-a-position	10.25	records	7687686
encoding                        blocks  head-B/rec  uncomp-B/rec  comp-B/rec  heads-alone  vs-varint  vs-no-scalars  vs-old-head
no head scalars                    160       2.195        24.589       4.752        0.161     -7.51%         +0.00%       -3.89%
varint, head before 2026-09-04     160       6.195        29.589       4.944        0.253     -3.76%         +4.05%       +0.00%
varint                             160       8.195        30.589       5.138        0.509     +0.00%         +8.12%       +3.91%
fixed 4/4/1/4/4/4                  160      23.195        45.589       5.433        0.630     +5.75%        +14.34%       +9.89%
fixed 4/2/1/2/2/4                  160      17.195        39.589       5.324        0.590     +3.62%        +12.04%       +7.68%
fixed 2/1/1/2/2/2                  160      12.195        34.589       5.233        0.552     +1.86%        +10.14%       +5.85%

sample	hg002-chr21			mean-reads-a-position	280.32	records	56627
no head scalars                     15       4.660        38.214      14.870        1.800    -13.43%         +0.00%       -6.07%
varint, head before 2026-09-04      15       8.661        43.215      15.831        2.425     -7.84%         +6.46%       +0.00%
varint                              15      11.588        45.143      17.178        3.945     +0.00%        +15.52%       +8.51%
fixed 4/4/1/4/4/4                   15      25.660        59.214      18.016        4.476     +4.88%        +21.15%      +13.80%
fixed 4/2/1/2/2/4                   15      19.660        53.214      17.783        4.324     +3.52%        +19.59%      +12.33%
fixed 2/1/1/2/2/2                   15      14.660        48.214      17.540        4.167     +2.10%        +17.95%      +10.79%
```

**The fixed-width rows say the same thing they said at Milestone D2 and by a wider margin:**
variable-length wins on both samples and at every width tried, now that the head has two more
fields to pay for.

**And the skipping walk is faster, not slower, than the figure on record.** On a store ng wrote
itself — 8,105,483 loci over tomato SRR7279481's whole genome at 9.7 reads a record, keeping one
record in a hundred over seven timed rounds — **2.930×**, 0.355 s against 1.041 s, stepping over
99.0 % of the body bytes.

## Two things this step found that the plan did not anticipate

**1. The probe had not been updated since the chain ids joined the head, so every head-cost figure
in the specs is older than the plan thought.** `ng_psp_head_encoding.rs` writes the head by hand and
checks itself against the shipped `BlockBuilder`'s bytes. That check failed on the first run of this
step — the harness wrote 2,489,991 bytes where the builder wrote 2,678,719 over one block — because
the harness never wrote the **chain ids' live-set changes**, which are part of the head and arrived
at Milestone E4. It now writes them through the real `LiveSetWriter`, and the check passes on both
corpora. So the figures H3 was asked to re-take were stale for two reasons, not one, and the head
they described had four scalars and no changes.

**2. The 9.2 % and the 5.8 % cannot be re-taken, and nothing replaces them exactly.** They compare
this format against spec §4.3's *no-head* row, whose bodies code coverage and chain ids as
differences from the previous record — a format nothing has ever implemented, so there is nothing to
measure it on. What replaced them in the spec is a narrower quantity that is exact: **the same
bodies with and without the head's own bytes**, which is the table above. It is a different
denominator and the spec now says so, in place, rather than letting one number be read as the
other's successor.

## The corpora had to be rebuilt, and one of them is not the same size

**The psps the figures were first taken on are refused by this tree's reader.** Their header
predates production's `kind` field, so `PspReader::open` rejects them naming it. Both were rebuilt
with the current `pileup`:

- **tomato SRR7279481** — the benchmark CRAM, whole genome. **7,687,686 records**, which is exactly
  what Milestone D2's corpus had, so this arm is the same corpus re-made.
- **HG002 chr21** — the same alignment the old psp came from, restricted to chr21.
  **56,627 records** against D2's **74,623**. Not the same corpus: the pileup that wrote the old one
  was a different build. Its depth, which is what the arm is for, is 280.32 reads a position against
  D2's 279.99.

## Changes made

**[`examples/ng_psp_head_encoding.rs`](../../../../examples/ng_psp_head_encoding.rs)** — the head is
six fields, and three arms are new:

- The two new scalars in both the varint and the fixed-width writers, with the fixed arm's widths
  extended (one byte for the kind, four for the compared count at the width the format allows).
- **The chain-id live-set changes, written for every arm** through `LiveSetWriter`, because they are
  the head's last field and are not what the arms differ in.
- **`VarintBeforeTheKindAndTheDenominator`** — the head as it stood before this milestone,
  reconstructed exactly: four scalars, and the locus-kind tag written at the end of the *body*,
  which is where it lived. The bodies are otherwise identical, so the difference between that row
  and the varint row is the whole cost of the change, measured on one corpus in one run rather than
  by building two versions of the code.
- **`NoHeadScalars`** — the head's scalars removed and its changes kept, as the denominator. Its doc
  says plainly that it is **not** spec §4.3's "no head" row and must not be quoted as one.
- The rows are collected before any is printed, so each can be shown against a reference measured
  after it, and a row whose widths cannot hold the corpus prints its refusal rather than `inf`.
- `put_head`'s six values became a struct: they are all `u64`, and a pair transposed at a call site
  would be invisible.

**Documents.** [`psp_file_format.md`](../../ng/spec/psp_file_format.md) §4.3: the 2.06× replaced by
the 2.930×, with what the old store was and why it read high; the 9.2 %/5.8 % replaced by the table
above with its different denominator named; the 0.077 bytes a record replaced by the head compressed
on its own at both depths, saying that the old figure left the changes out. The fixed-versus-varint
paragraph now states that the question was settled by measurement at D2 rather than left open.
[`psp_record_encoding.md`](../../ng/spec/psp_record_encoding.md) §6: the paragraph saying *how much
of the 2.06× survives at depth is unmeasured, and nothing can measure it until Milestone F opens a
file* — both halves false — replaced by the measurements, with the deep end still named as a bound.
[`arch/psp_file_format.md`](../../ng/arch/psp_file_format.md): the narrowed open question records
that the shallow end is no longer an upper bound and the deep end still is.
`src/ng/psp/record.rs` and `src/ng/psp/walk.rs`: the same figures where the code quotes them.

## Validation results

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — **6,157 passed**, 0 failed, 14 ignored; `cargo test --example
  ng_psp_head_encoding` — 4 passed.
- **The probe's own oracle passes on both corpora**: its varint arm reproduces `BlockBuilder`'s
  block payloads byte for byte, which is what makes the fixed-width rows a comparison against ours
  rather than against something nobody wrote.
- Every number above was re-taken after the last refactor of the probe and reproduces exactly.

## Left standing

- **No ng-written store at 280 reads a record exists**, so the deep end of the skip's value is still
  an upper bound taken on a converted store. Building one needs a deep human alignment walked by
  `generate-psps`.
- **`psp_record_encoding.md` §2.3's "the length costs 1.4 % of the file at three reads a position
  and 3.3 % at 279"** is from the same old regime and is not re-taken here: isolating one head field
  needs an arm per field, which is more than this step's question.
