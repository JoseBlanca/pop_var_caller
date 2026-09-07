# A2 — the decode reads a window of the reference, not a chromosome

**Date:** 2026-09-07
**Branch:** `ng-cram-window`
**Plan:** [cram_reference_window.md](../../ng/impl_plan/cram_reference_window.md), Milestone A step 2
**Design:** [alignment_cursor.md](../../ng/spec/alignment_cursor.md) §10 points 2–4;
[`FORK.md`](../../../../vendor/noodles-cram/FORK.md) §5

## Verdict

**The CRAM decode now takes the bases each slice declares it needs, and nothing else.** Every
mapped slice is decoded through `records_over_window`; the repository handed to noodles is empty
from this commit on, so a slice that fell back to it would fail loudly rather than quietly find a
chromosome. 6,280 library tests pass, up 2 — and, for the first time in this tree, a test can see
a reference read at the wrong offset.

## What landed

`decode_container_at` asks each slice for `reference_span()` — the contig and the first and last
position its records touch, readable from the slice header before any block is decoded — fetches
exactly that span through the reader A1 gave it, and decodes with `records_over_window`. A slice
with no single reference takes the existing call, which needs no external bases at all.

**One deviation from the plan, forced by the borrow checker.** The plan put the window buffer in
`RecordScratch`. It cannot go there: the decoded records borrow the window for as long as they
live, while `scratch` is borrowed *mutably* to build each record from them, so the two cannot be
fields of one struct. The window is its own local, reused across the container's slices exactly as
the plan intended.

## The fixture this needed, and why every existing one was blind

**Every CRAM fixture in this tree is written against an all-`A` reference**, and the in-memory
accessor the cursor tests hand out is all-`A` too. That is what makes them agree by construction,
and it is also what makes them **unable to fail for a wrong window**: a CRAM stores each read as
its *differences from the reference*, so decoding against the wrong offset of an all-`A`
chromosome rebuilds exactly the same read — and the slice's stored MD5 matches as well, since
every window of that reference has the same digest.

So A2 adds `indexed_cram_over_a_varied_reference`: a FASTA of pseudo-random `{A,C,G,T}` from a
fixed-seed xorshift, so the fixture is identical on every machine, with a CRAM written against it.
`a_cram_decoded_against_per_slice_windows_rebuilds_the_reads_its_bam_twin_holds` compares every
read's name, position, bases and qualities against the BAM twin of the same records — a BAM being
the right oracle because it stores its bases literally and reads them back without consulting a
reference at all. Filter #8 is off in that test on purpose: the fixture's reads are all-`A`
against a random reference, so the mismatch filter would drop every one and the comparison would
be two empty lists agreeing.

**Mutation-verified.** Fetching the span one base late fails it — measured twice, at 304 passed /
1 failed in `ng::read`, and the mutation reverted and the revert confirmed by grep each time.

**What actually caught the shift is worth recording, because it is not what the test asserts.**
The message is `reference sequence checksum mismatch` — the MD5 the slice header stores over its
own span, which noodles validates against whatever bases it is handed. So on a file carrying that
digest the window is guarded twice and the digest fires first. **The sequence comparison is what
remains when it does not**: CRAM v3 §8.5 says an all-zero checksum is not to be validated, and a
writer may store one, in which case wrong bases would reach the caller with nothing else
objecting.

**A window that falls short is refused, and by the bounds check rather than the digest.**
`a_reference_window_that_falls_short_of_its_slice_is_refused` needs an accessor that lies, because
nothing in the honest path can produce a short window — the decode asks for exactly the declared
span, and a real reference either serves it or errors. The double serves one base fewer than
asked. The assertion requires the refusal to name the *window*, narrowly: accepting the digest
instead would let the only guard silently become the one some real files switch off.

## Validation

Run in the container on this tree.

- `cargo test --lib` — **6,280 pass**, 0 failed, 15 ignored, in 50.97 s.
- `cargo test --lib --tests` — as above; every integration binary passes except
  `ng_calling_loop_calls_genotypes`, `main`'s pre-existing contamination failure.
- `cargo clippy --lib --tests --all-features -- -D warnings` — nothing in the changed files; the
  three `run/cohort_merge/` lifetime errors are `main`'s.
- `cargo fmt --all -- --check` — the 29 drifting hunks are `main`'s, in files this step does not
  touch.

**The existing CRAM oracles pass unchanged** — `a_run_of_regions_through_one_cram_cursor_matches_a_linear_scan`,
`the_cram_cursor_oracle_is_not_vacuous`, `a_read_past_the_first_container_carries_its_own_read_group`
and `t8`. They cannot see a window offset, for the reason above, but they do prove the windowed
path serves the same reads in the same order across regions and container boundaries.

## What is not proven here

**Nothing has been decoded from a real file yet.** Every number above is from fixtures whose
contigs are 100 and 200 bases. The whole-genome tomato CRAM — 112,140 slices, spans averaging
7,075 bases — is B1's, together with the resident-memory figure this whole plan exists for and a
VCF compared line for line. Until then the saving is arithmetic.

**`OpenReference` still holds its repository**, unread from this commit. A3 deletes it.
