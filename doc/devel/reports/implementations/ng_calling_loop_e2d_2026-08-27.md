# ng calling loop — E2d: the contaminant seed at a repeat tract

**Step:** E2d of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the third term of a
repeat tract's read-likelihood mixture.
**Design authority:** [`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §4.5.1, and
§3.6 for the ordinary-site term it is deliberately unlike;
[`spec/population_diversity.md`](../../ng/spec/population_diversity.md) §4 for the length spectrum.
**Date:** 2026-08-27. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

**A run whose parameter fit found contamination now has its repeat tracts called instead of
refused.** A contaminating read at a tract shows a length that is a real allele in some
individual, and until this step the model had to explain it as slippage — inflating the apparent
slip rate — or drop it to the outlier floor. Now it has a term of its own: how common each length
this tract can show is in the contaminating population.

## 2. Where the numbers come from, and the one that had to be turned inside out

The mixture is three terms: the sample's own copies, a flat outlier share, and `c · seed(o)` —
the read group's fitted contamination fraction times how common the length observation `o` showed
is among contaminants.

**The fraction is the pre-pass's and needed no work.** The seed did. The genotype prior builds its
belief about a tract **per candidate**, placing each candidate by how many whole repeats it holds
as an offset from the reference tract's count; `seed(o)` asks for a probability **per observed
length**, in bases. Converting one into the other is the calling loop's job because it is the only
place holding both the candidate table and the tract's reachable-length support.

**Each candidate's share is added at the support entry its bases land on**, which settles the
three cases the sibling module's documentation asks about: two candidates spelling one length sum
into one entry rather than each taking the full share; a length no candidate reaches gets nothing
and its reads fall to the outlier floor, which is where they went before; and a read that ran out
inside the tract takes the seed's mass at or above what it witnessed.

## 3. A tract's mixture is frozen where an ordinary site's moves, and that changes a cost claim

**At an ordinary site the contaminant's half of the mixture is the cohort's own frequency for the
allele an observation shows**, so the loop rewrites it at every pass and the genotype-likelihood
table has to be assembled again each time.

**At a tract it is not.** The spec weighed the cohort's own per-locus frequencies for the job and
refused them for exactly that reason — contamination must not move from one pass to the next — and
took the joint repeat fit's **length spectrum** for the tract's stratum instead. It is specific to
the locus, because it is indexed from that tract's own reference length, and it is fitted before
calling starts.

**So a contaminated tract's table is built once, like every other tract's.** The driver's per-locus
question is therefore *does this locus's table move as the loop iterates*, which is contamination
**and** the SNP/indel path — not contamination alone. A contaminated tract needs no per-batch
contaminant tables and no reassembly, and is given neither.

## 4. What it buys, measured on the fixture

Three samples at a dinucleotide tract called over 6 whole repeats and 7. The middle sample carries
two copies of the 6-repeat tract and shows twenty reads of its own plus **four** at the 7-repeat
length that came from another individual.

| the run | that sample's call |
|---|---|
| no fraction fitted | **`0/1`** — a second allele |
| fitted fraction 5 in 100 | `0/1` — not enough mass |
| fitted fraction 8 in 100 | **`0/0`** — somebody else's DNA |

**The middle row is the one that matters**, and it is asserted: the fixture cannot be satisfied by
a model that reads the fraction as a flag rather than as a number. Without any fraction, the other
two explanations cannot carry four reads — slippage to exactly one repeat longer runs about one
read in two hundred at this stratum's fitted numbers, and the outlier term, spread flat over every
length the tract can reach, is smaller still.

**Where the window is.** At one, two and three reads every run calls `0/0` — the slippage term
alone covers them. At five, no fraction below about 20 in 100 recovers the homozygote. Four is
where 5 in 100 and 8 in 100 differ.

## 5. What the reviews found, and what changed

Three review agents in isolated worktrees: arithmetic and control flow; tests and mutation; design
conformance and claim-checking. **2 Blockers, 8 Majors, and 7 of 42 checked claims wrong.** Every
finding was applied.

### No arithmetic defect, and no test that could have found one

Both the arithmetic and the design reviewers confirmed the seed sums to one on every path, that
every candidate's own length is in the support by construction, and that skipping the per-batch
tables at a tract is safe. **The mutation reviewer then showed the whole suite stayed green with
the computation removed.** Two mutations, both run:

- replace the fit's length spectrum with a uniform shape — **green**;
- re-centre the spectrum on the last candidate instead of the reference — **green**.

Two accidents caused it. The unit fixture that asserted the seed's values ran on a fit carrying
**no length spectrum at all**, so the prior answered from the ladder's flat bottom rung and "half
each" was both the right answer and the mutant's. And the fixtures that *did* have a fitted
spectrum asserted only genotypes — which did not move, because the spectrum was
`[0.10, 0.25, 0.45, 0.15, 0.05]`, whose upper tail falls by a factor of three at each step, so a
shape read one repeat off centre gave the two candidates the **same** pair of shares.

**Both are closed by one fixture and one repair.** The spectra are now
`[0.10, 0.30, 0.44, 0.11, 0.05]` and its siblings, with no two adjacent pairs in the same ratio.
And a new test asserts the seed's own numbers over three candidates chosen so that no shortcut
reproduces them: a clean 6-repeat tract of twelve bases, an **interrupted** one of twelve bases
holding four whole repeats, and a clean 7-repeat tract of fourteen. Twelve bases carries the first
two summed — `0.54/0.65` — and fourteen carries `0.11/0.65`.

Five mutations that ran green now fail:

| mutation | now caught by |
|---|---|
| a uniform shape instead of the fit's spectrum | the seed's own numbers |
| the spectrum re-centred on the last candidate | the seed's numbers, and a new refusal |
| each share written at its length rather than added | the two twelve-base candidates |
| a candidate placed by its repeat count rather than its bases | the same, which cannot even complete |
| the driver treating a contaminated tract as a moving table | the contaminant tables it must not size |

### Two lookups became one

The gather and the driver each looked the length spectrum up, keyed identically — so the run could
report one rung of the tract ladder while scoring against another if either key changed. The first
fix was a comparison; the better one, and what shipped, is a single lookup handed to both, so the
two cannot disagree. The comparison is gone with the possibility.

### Nine sentences the change made false

Two `# Panics` sections still promised a refusal that had been deleted, and the field documentation
on the emission-cost counter still said a contaminated locus assembles its table once per pass.
Also corrected: a claim that a tract's interrupted candidates cannot be told apart by the reads,
where the spec says the read likelihood separates them by about 28 Phred per distinguishing base
and it is only *this* term that cannot; and an assertion message claiming the locus's warrant
covers the contamination fraction, which it does not.

### One number the reviews asked to be stated with its size

The prior's shares are normalised over **this locus's candidates**, so the spectrum's mass at
lengths no candidate carries moves onto the ones that do. On the shipped fixture the fit puts 0.44
of a stratum's chromosomes at the reference length and the seed says 0.8 — so at a fitted fraction
of 5 in 100 the mixture credits 0.040 of a read at that length to the contaminant where the
spectrum alone would credit 0.022, **1.8 times the fitted weight**. The row's own sum-to-one
contract forces some normalisation; what is now written down is what it costs.

## 6. Validation

All under `./scripts/dev.sh`:

| gate | result |
|---|---|
| `cargo test --lib` | **4,896 passed / 0 failed / 14 ignored** |
| `cargo test --test ng_calling_loop_calls_genotypes` | **15 passed** |
| `cargo test --test ng_calling_loop_allocation --features dhat-heap` | 1 passed |
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo doc --no-deps --lib` | 27 unresolved links, all pre-existing |

**The release-held assertion battery ran on the three checks this step adds**: each downgraded to
`debug_assert!` in one run, `cargo test --release --lib ng::calling --all-features` run, **four
failing tests over three checks** — the two-direction one is reached from both sides. Source
restored and the restore verified.

## 7. Banked for the owner

- **No contaminated fixture has a read that ran out inside the tract.** The seed's *ordering* is
  load-bearing only there — such a read takes a suffix of it — and the row's own tests cover that
  arithmetic, so what is untested is the seam rather than the sum.
- **Every contaminated fixture gives every library the same fraction**, and only one library sends
  a read. Nothing here would notice the fraction list arriving reversed; the row's own tests cover
  its indexing.
- **No candidate sits outside the fitted spectrum's span**, so the floor that a candidate several
  repeats from its tract's length takes is not reached through this path.
