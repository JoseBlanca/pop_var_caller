# An insertion at a region's last base belongs to whatever follows, not to the SNP/indel path

**Date:** 2026-09-02. **Follows** [C4](ng_ssr_loop_c4_2026-09-02.md), which measured the defect.
**Owner's ruling**, on the example C4 printed: *"the A belongs to the str, not to the snp/indel.
It would be different if it were an indel that crosses the boundary, but this is not the case.
This is an indel inside the str that happens to be coded inside the generic segment, but it does
not belong in there."*
**Module:** [`src/ng/locus_generation/pileup/open_record.rs`](../../../../src/ng/locus_generation/pileup/open_record.rs).

---

## The rule, and why it is the clip's rule rather than a second one

An insertion's anchor base is **unchanged**. It is in the record because a VCF record needs a
base, not because the variant touches it; the inserted bases sit *between* that base and the next.
So an insertion anchored on a region's last base has all of its content past the region's end —
and the bound that already truncates a deletion's tail there leaves nothing of it.

**What is untouched**: an indel that genuinely starts inside the generic ground. A deletion that
removes the last base and continues into the tract still opens its record — its anchor is earlier,
the deleted bases include one of the region's own, and the existing clip gives it the part inside.
That is the owner's second case, and it stays a SNP/indel call.

**Two places, one predicate.** `inserted_past_the_regions_last_base` is asked where an event opens
a record — so no record is opened for such an insertion alone — and again where a read's allele
bases are composed, so its bases stay out of a record something else opened. The second is the one
that matters: the read almost always has a `Match` at the same anchor, so a record exists there
either way, and what changes is that it now holds the reference base alone.

## What it costs and what it buys, measured

Three GIAB samples, each on its own 100-interval confident BED against its own truth VCF,
`--defaults`. The middle column is C4's run; the last is the same binary with this rule.

| depth | class | | no tract calling | tract calling | + this rule |
|---|---|---|---:|---:|---:|
| 30× | indels | recall | 0.673 | 0.915 | **0.900** |
| | | precision | 0.987 | 0.816 | **0.983** |
| | SNPs | recall | 0.974 | 0.980 | 0.980 |
| 50× | indels | recall | 0.676 | 0.939 | **0.933** |
| | | precision | 0.987 | 0.809 | **0.981** |
| | SNPs | recall | 0.979 | 0.984 | 0.984 |

**Duplicated indel records: 62 at 30× and 66 at 50×, now 0 at both.** False indel calls fall from
68 to 5 at 30× and from 73 to 6 at 50× — below the 3 the run had before tracts were called at all,
plus a handful.

**It costs five true indels at 30× and two at 50×**, which is the honest half. Of the 27
insertions anchored on a generic region's last base on HG002's ground, 25 were also written by the
tract path and are pure duplicates; the other two are an insertion at a **bundle's** edge, which
nothing can call either way, and one at a tract the repeat path called as reference. Those are now
silence rather than a call.

## Against the bar

The production caller on the same samples and regions is 0.987 SNP and 0.930 indel at 30×. ng is
**0.980 and 0.900**, from 0.974 and 0.673. And against ng's own best — 0.946 indel before the
region clip landed — it is 4.6 points short, so the clip is still not fully paid back.

## Tests — 1 new

`an_insertion_anchored_on_the_regions_last_base_contributes_no_allele`, written against
`process_position` directly. **Its first draft was written against the generator and was
vacuous**: that fixture's reference is a run of `A`s, so an inserted base left-aligns to the
contig's start and no record opens at the anchor whatever the rule says. The draft passed before
the rule existed. It is recorded in the surviving test's own comment.

Both halves are asserted, and the second is what stops the rule being *drop insertions*: one base
further inside the region the same insertion still reaches the record, as ordinary SNP/indel
ground.

## What is left open

**A run whose analysed ground simply ends at a region boundary loses an insertion there**, because
the walk is given one region at a time and cannot see whether another follows. On HG002's
benchmark ground that costs nothing — all 27 such insertions are followed by a typed region — but
it is a real hole for a BED that ends mid-repeat, and closing it means the generator knowing its
region's successor.

## Validation

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` 6,014 passed / 0 failed / 14 ignored, all in the container. The benchmark
figures are the shipped release binary, scored with `benchmarks/giab/src/score_ng_recall.sh`.
