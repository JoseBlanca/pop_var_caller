# psp mode — Milestone G: the census beside the psp, and the walk is once

**Date:** 2026-09-04
**Plan steps:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone G, steps G1 and G2
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §2, §5.2; `parameter_prepass_joint_records.md` §6.1; `parameter_prepass_census_sites.md` §3, §5.1
**Branch:** `ng-psp-mode`

## The answer

**`generate-psps` writes two files a sample from one pass over the alignment files.** Measured on
the six tomato accessions over the two 100 kb intervals: **6 samples, 3,592,149 bytes of psp and
1,305,915 bytes of census, in 5 seconds**, and the run names both files and both sizes for every
sample.

The census is fed **at the walk's yield point**, where a locus is handed over — not by the psp
writer. Feeding it from the writer would make the census a function of what was *stored* rather
than of what the walk *saw*, and would leave a consumer that iterates a gatherer without writing
a psp building no census at all.

## What the numbers behind the selection are, and where they come from

Which positions and which tracts a run keeps is a function of a seed and two counts
([`CensusSelection::SHIPPED`](../../../../src/ng/run/gatherer.rs)). Two are the design's own
figures rather than choices made here:

- **about two million generic positions** — *"the knob is a number of positions, and about two
  million of them is the default"* (`parameter_prepass_census_sites.md` §5.1), sized to yield
  about ten thousand segregating sites;
- **five thousand tracts a stratum** — where `parameter_prepass_joint_loci.md` §6's first
  question closed, measured on a tomato archive.

**The seed is a compiled-in constant, and that is what makes this command's own advice work.**
`generate-psps` tells a person to spread a cohort by running it once per sample; two invocations
that seeded differently would keep **disjoint** sets of positions and their samples could not be
pooled at all. A constant makes them agree across processes by construction rather than by
somebody typing the same number twice —
`two_invocations_keep_the_same_census_positions` is the test.

**Whether these three become flags is not settled**, and it is the same open question Milestone C
recorded about the read filters and the locus-generator knobs.

## The reference is now read on the walk's own thread, observed

Choosing census positions needs to know where the reference is sequence at all: a position inside
a run of `N` has no base to compare a read against. That is a pass over the FASTA, and the reader
`generate-psps` used hands its pass to a background thread and takes the result from a **shared
cache** — which cannot hand back one caller's *observations*. So a second entry point was added
(`read_reference_observing_or_creating_fai`) that reads on the calling thread and tells an
observer about every base. **It is the same single pass, not an extra one**; what it gives up is
overlapping that pass with opening the alignment files.

## A defect this milestone's own test found before it shipped

**The census would have named a psp that does not exist.** `PspWriter::create` records the
compression level into the header *before* encoding it, so the header the gatherer holds is not
the header in the file — one line of TOML, and every byte of the digest differs. A census built
from the gatherer's copy names a header no psp carries, and
[`Freshness`](../../../../src/ng/parameter_estimation/joint/census_file.rs) would have answered
*rebuild* for ever, on every census, silently. That is precisely the failure naming the pileup
exists to prevent.

**The fix is that the writer hands back what it wrote**: `WriteStats` now carries
`header_digest`, the md5 of the header exactly as it went into the file — sixteen bytes rather
than the header itself, which runs to 16 MB at the format's ceiling. The test that caught it
rebuilds the expected identity from the psp on disk rather than from anything the gatherer holds.

## What the pair does when a run stops

Both files go to a scratch name and are renamed only once whole, so a stopped walk leaves neither
at the sample's own path and a stopped **re-walk** leaves the pair it was replacing intact.

**Two renames cannot be one**, so a run that dies between them leaves this sample's new census
beside its old psp — and the run says so, because the second rename's failure ends it. What makes
that state safe rather than silent is the identity: the census names the header and record count
of the psp it was built from, so a fit reaching that pair refuses it. The order is census first
for that reason — renaming the psp first would leave a **finished-looking psp** beside a stale
census, and the psp is what a calling run reads and would trust.

## One figure that must not be read as a general one

On this slice the censuses come to 1.3 MB against 3.6 MB of psp — **about a third**. That ratio is
an artefact of the corner: the analysed ground is 199,672 bases and the budget is two million
positions, so the selection keeps essentially every position. On a whole genome the budget is far
below the ground and the census is a small fraction of it. The ratio measured here is a fact about
a 200 kb BED, not about the format.

## What is knowingly left

- **A census already in the output directory does not trigger the `--force`-less refusal**; only
  a psp does. A census without its psp is useless, so what this can replace is a file nothing can
  read.
- **The fit does not read these files yet.** Reading a census back, building one *from* a psp, and
  the §7.12 byte-for-byte census-equality oracle are the next plan's
  (`parameter_prepass_runs.md`, unwritten), which is what this milestone's scope says.

## How it is verified

Eight tests: four in the gatherer (both files written and the census naming its psp; no plan means
no census; a census that cannot be written fails the walk; what the census holds is what the
sample showed) and four at the command (a census a sample naming its own psp; one run's samples
select alike; **two invocations** select alike; a stopped walk leaves no scratch file; the report
names both).

`cargo test --lib` in the container: **6,229 passed, 0 failed, 15 ignored**. `cargo fmt --check`
and `cargo clippy --all-targets --all-features -- -D warnings` clean. Both standing oracles still
green — direct mode byte-identical, and psp mode's VCF byte-identical to direct mode's apart from
`##commandline`, 599 records.
