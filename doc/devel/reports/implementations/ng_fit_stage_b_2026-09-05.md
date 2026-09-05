# The fit stage — Milestone B: `generate-census`, and what each route to a census costs

**Date:** 2026-09-05
**Plan steps:** [parameter_prepass_runs.md](../../ng/impl_plan/parameter_prepass_runs.md) Milestone B, steps B1, B2 and B3
**Spec:** `parameter_prepass_joint_records.md` §6.1; `parameter_prepass_census_sites.md` §2, §3
**Branch:** `ng-psp-mode`

## The answer

**`generate-census` builds each stored psp's census without opening a single alignment file, and
what it writes is what the walk wrote.** On the six tomato accessions over the two 100 kb
intervals, all six censuses are byte-identical between the two routes.

**Building the census during the walk is the cheaper of the two, by about a tenth of the walk's
time**, and costs a little more memory: **1.28 s against 1.40 s, and 192 MB peak resident against
188 MB**, over three repetitions that moved by 0.02 s and 1 MB.

## What was measured, and on what

Six tomato accessions (`SRR7279481`, `488`, `501`, `533`, `536`, `537`) over the two 100 kb
intervals of `tmp/c1_two_regions.bed`, at about three reads a position. The ground is the one a
real run walks — `run_ground::segments_over`, which cuts it into **318 segments**; asking the
catalog directly at the storage floors gives several times as many over the same ground, and a
cost measured over that segmentation is a cost no run pays.

Each route runs in a process of its own, because the peak resident memory of two routes in one
process is the larger of them and says nothing about either
([`scripts/ng_census_route_cost.sh`](../../../scripts/ng_census_route_cost.sh),
[`examples/ng_census_route_cost.rs`](../../../examples/ng_census_route_cost.rs)).

| route | wall time over the work | peak resident memory |
|---|---|---|
| the walk writes both files | 1.28, 1.28, 1.28 s | 192.3, 192.2, 191.3 MB |
| the walk writes the psp; a second pass builds the census | 1.40, 1.39, 1.41 s | 188.6, 188.4, 188.4 MB |

The clock starts after the setup — reading the reference, building the segmentation, choosing the
selection — because that work is the same whichever route runs and including it would dilute what
the measurement is about. It took 3.14 s here.

**The memory difference is 3.5 MB in 190, under 2%, and this report does not say what it is.**
Both routes hold the same 198,182-position selection, so it is not the selection; beyond that
would be a guess.

## One figure that must not be read as a general one

**The census keeps 198,182 of the 200,000 analysed bases here — 99 in 100.** The budget is about
two million positions and this ground is 200 kb, so the threshold keeps very nearly everything.
On a whole tomato genome the same budget keeps roughly **1 position in 400**.

So this measurement charges the census something like four hundred times the share of the walk it
would carry on a whole-genome run, and **the same is true of the file sizes**: 1,305,915 bytes of
census against 3,586,053 of psp, about a third, is a fact about a 200 kb BED and not about the
format.

What that does to the comparison is worth being explicit about. The two routes build the same
census, so what separates them is the second route's extra pass over the stored psp. That pass
scales with the psp rather than with the census, so **the gap measured here is not one that
shrinks away on bigger ground** — but this run cannot say how it grows, and nothing here should be
quoted as if it could.

## A difference the check caught before any timing was believed

The first run of the comparison reported **all six censuses different**, while the fixture tests
said the two producers agree. The difference was exactly sixteen bytes at one offset in every
file: the md5 of the psp header, which a census carries to name the psp it was built from.

The cause was the harness, not the routes. It recorded its own command line — including the route
word — into each psp's provenance, so the two routes wrote psps whose headers differed by one
character, and the censuses correctly named two different files. The harness now records a
constant command line, and the comparison is what it claims to be.

**This is why the script compares the files rather than only timing them.** A timing comparison
between two different outputs measures nothing, and the failure looked exactly like a real defect
in the second producer.

## What the command says about each file it read

Per sample, in one line shared by the progress note and the final report: the psp read and how
many stored loci it holds, the census written and its size, and how many kept positions and kept
tracts this sample has a read at — **each beside how many the selection holds**, because 900
positions is good or bad depending on whether the selection kept a thousand or two million.

**A sample whose census holds no read at any kept locus is named rather than omitted.** That is a
legitimate outcome — the walk covered ground the selection kept nothing in, or no read reached
what it did keep — and a run that said nothing about it would leave somebody hunting for a file
written exactly as asked. The fixture cohort's second sample has no reads at all and exercises it.

## Two things the command does not do

**There is no `--regions`**, for the reason `call-from-psps` has none, and here it is sharper than
convenience: the digest of the analysed regions travels in every census as one of its recording
terms, so a selection made over other ground produces censuses the cohort cannot be fitted from —
and the disagreement would surface hours later, at the fit.

**It does not fit anything.** Nothing yet reads a census back at command level; that is Milestone C.

## Also fixed here

`generate-psps`' own help said it writes no census beside the psp and does not refuse to overwrite
one already there. Both shipped before this plan started — the census in Milestone G of the
psp-mode plan, the refusal with `--force` earlier still.

## Validation

`cargo test --lib` in the container: **6,246 passed, 0 failed, 15 ignored** — the 6,236 Milestone A
left, plus this milestone's ten. `cargo fmt` and
`cargo clippy --all-targets --all-features -D warnings` are clean.
