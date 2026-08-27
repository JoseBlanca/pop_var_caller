# ng calling loop — E3b: ng calls genotypes at a repeat tract

**Step:** E3b of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the driver's
repeat-tract branch, and the integration fixture in front of it.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §1, §5.0.1,
§8, §9; [`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §4 for the row;
[`spec/population_diversity.md`](../../ng/spec/population_diversity.md) §4, §5 for the prior.
**Date:** 2026-08-27. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

**ng calls genotypes at a repeat tract.** Until this step the driver turned every tract away at
its front door: its emission build, its row assembly and its evidence accessor all took the
SNP/indel path's per-sample evidence, and a tract's evidence is a different shape. Five places
now branch, the refusal is gone, and a tract's prior is seeded from its stratum's fitted **length
spectrum** — how that stratum's chromosomes are spread over tract lengths — with the rung of the
three-rung ladder it came from travelling onto the locus's record. The integration fixture runs
the same path at a tract that it already ran at a SNP:

```text
per-sample observations + supplied candidates → shape_ssr_locus → call_locus → LocusInference
```

## 2. The five branches, and what each one does

| place | SNP/indel | repeat tract |
|---|---|---|
| evidence accessor | `generic_evidence_of` | new `tract_evidence_of`, returning the samples, the tract's motif and flanks, and the candidates' repeat counts |
| `fill_what_no_pass_recomputes` | fills the emissions only | fills the **whole likelihood table**, through `fill_tract_likelihood_table` |
| `assemble_genotype_likelihood_table` | folds the emissions into a per-genotype row | returns early — there is nothing left to fold |
| `weakest_warrant_at_the_locus` | the calibrations of the read groups whose reads reached the locus | the tract's own `(read group, candidate)` parameter table's warrants |
| the per-locus prior seed | `fill_locus_concentration` off the run's frequency spectrum | `fill_ssr_seed` off the stratum's length spectrum |

**Why the table is filled in the emission build rather than in the assembly, which is the one
structural choice here.** The tract's row builder computes its emissions and assembles them into
a per-genotype row in one call, and there is no seam between the two — that is
`likelihood::ssr`'s shape, not a decision taken here. So the fill had to live in one of the two
functions. Putting it in the assembly would have meant threading the repeat-tract emission model
through `TableReassembly` and `run_frequency_loop`, which are generic over the model's *scratch*
rather than over the model; putting it in the emission build needed no signature change, because
that function already took the model and had never used it.

**And it changes no invariant.** Where nothing is contaminated a tract's row reads no allele
frequency either, so the table written once before the loop is the table every pass reads —
which is exactly what Milestone D's cost invariant says. Where something *is* contaminated a
tract is refused inside that same call.

**The function was renamed, because on the tract path its old name was half true.** It is now
`fill_what_no_pass_recomputes` — what the two paths have in common — where `build_locus_emissions`
named only what the SNP/indel path does.

## 3. What a tract still cannot be

**A contaminated one.** The third term of the tract's read-likelihood mixture — how common the
length an observation showed is in the contaminating population — is not built, and
`TractScoringFits::locus_parameters` refuses a run whose fit found one rather than handing back
the two-term form. That refusal now fires from inside a scored tract rather than from a driver
that refused every tract, and it names the step that supplies the missing term. **It is a
mechanism rather than a doc comment for a measured reason**: the two-term row returns perfectly
plausible numbers with the fitted fraction silently dropped.

**A bundle of tracts.** The calling seam sends a bundle down the repeat path deliberately, and
nothing on that path scores one. It is refused by name, and the message says which of the two
failures it is — a gap in what the repeat path covers, not a locus routed to the wrong read
model. Nothing constructs a bundle into calling today.

**Selected rather than supplied.** The repeat-tract half of candidate selection is unwritten, so
a tract's candidates and each candidate's repeat count are stated by whoever calls the loop.
Every fixture says so in its own doc comment, and so does the test binary's module documentation,
because a supplied candidate set read as a selected one is a claim about a step that does not
exist.

**Why the repeat counts travel at all**: how many whole repeats a candidate holds is not its byte
length divided by the motif's. An interrupted tract — one whose repeat is broken by a
substitution — holds fewer, and the count is what picks the stratum whose slippage numbers score
the candidate and the offset the prior's length spectrum reads at.

## 4. The one point where this step and a spec section disagreed — ⚖ ruled 2026-08-27

**`population_diversity.md` §5 wanted a repeat tract in a run carrying no repeat-tract parameters
refused by name; §4.4 wants the tract ladder to always answer.** The two meet at exactly one run:
one whose fit produced no length spectrum anywhere. The driver took §4.4's side — such a tract is
called, and the rung on its record says `StatedFlat` — and put the question to the owner.

**Ruled: keep it.** *"Refusing turns a whole class of runs into a hard failure for a condition the
output already states."* The rung is what that spec's third goal asks for — a call resting on a
stated constant has to be distinguishable from one resting on a measurement without re-running
anything — and it carries that distinction at every tract of such a run, where a refusal would
carry it by producing nothing.

**No code changed on the ruling**; it is what this step built. §5 and §6 of the spec now record the
ruling rather than the refusal, and the plan's E3b entry carries it too. The refusal §5 used to ask
for would have been one predicate in `call_locus` — no length spectrum at any stratum and none at
any motif period, which `StratumFits::strata_with_a_length_spectrum` and
`periods_with_a_pooled_length_spectrum` answer between them.

## 5. What the reviews found, and what changed

Three review agents in isolated worktrees cut from `72cc100c` with the working tree applied as a
patch: arithmetic and control flow; tests and mutation; design conformance and claim-checking.
**No Blockers. 5 Majors, 11 Minors, 5 Nits, and 13 of 44 checked claims wrong** — every counted
figure was right and every failure was an explanation of why or where, which is the fifth step
running with that shape. Every finding was applied.

### The arithmetic was sound and the fixtures could not have shown it

**The finding worth the whole review**: nothing in the step could fail if the tract's prior seed
were wrong, and nothing could fail if the driver ignored the repeat counts it goes to such
lengths to thread through. Both were demonstrated rather than argued — one reviewer replaced the
supplied repeat counts with `bases.len() / period` and ran the suite green; another took the
reference repeat count from the *last* candidate instead of the first, which is the mistake the
code's own comment warns about, and ran green.

Both survived for the same reason: **every fixture gave each sample twelve to twenty reads all of
one length**, where the likelihood separates the genotypes by tens of nats and the largest shift a
wrong prior can produce is about four. And **every candidate of *n* repeats spelled exactly `2n`
bases**, so counting the bases gave the right answer by accident.

The fixtures were rebuilt to remove six coincidences, and four mutations that had run green now
fail, each caught by one test:

| mutation | now caught by |
|---|---|
| repeat counts derived from the candidates' bases | the interrupted-tract fixture: two candidates of twelve bases holding six and four whole repeats |
| the reference repeat count taken from the last candidate | three tests — the seed's own numbers, the interrupted tract, and the tract whose reads decide nothing |
| the row count charged where the candidate count belongs | the cost fixture, now three samples against two candidates |
| only the reads that spanned the whole tract counted | the same, now holding five observations of which one ran out |

**The other coincidences removed**: one library became three, against two candidates, in two
slippage groups with different numbers; the length spectra stopped being palindromes and no two of
their classes share a weight; the two strata stopped sharing their direction split, their fall-off
and their substitution rate.

**Two new fixtures exist only to make the prior testable**: one asserts the seed as numbers —
`20 × [0.45, 0.15]` = `[9.0, 3.0]`, the fitted weights at offsets 0 and +1 from the reference tract
length — and one gives a sample no reads at all, so its posterior is its prior and the call is
whatever the length spectrum says. Against the right seed the three genotypes come out 0.577,
0.346 and 0.077 and the call is `0/0`; against the seed a wrong reference count produces it is
`1/1`.

### Four claims that were wrong about the code rather than about a number

- *"Only the read groups whose reads are actually here"* — the fold's first rule, and no longer
  true at a tract. It is now qualified, and the departure from `read_likelihoods.md` §4.4 is
  stated with what narrowing it would cost. **It also gained a test**: a run of three libraries
  where only one sends a read and the fit describes only that one comes back `Defaulted`.
- *"The two paths read two different sets of parameters, and neither reads the other's"* — false,
  because a tract's **site quality** is still folded against the run's ordinary-site prior seed.
  That is pre-existing (`calling_quality.md` §8 leaves a tract's quality to a document that is not
  written), but this step is what makes it reachable, so it is now named at the call site.
- *"Two things are checked"* on the evidence's ordering contract — three, since this step.
- *"Its six vectors are cleared and refilled at each tract"* — four.

### Three checks the reviews asked for, and each is reached by a test

`assemble_genotype_likelihood_table`'s tract arm now asserts that the run is uncontaminated rather
than relying on another module to have refused first — because the day the tract's third mixture
term lands, a contaminated tract falling through that early return would be scored against no
contaminant frequency at all, silently. A bundle is refused with its own message. Both are
unreachable through the driver and both are reached by a test that calls the function directly;
downgrading either to `debug_assert!` fails exactly one test in `--release`.

## 6. Validation

All under `./scripts/dev.sh`:

| gate | result |
|---|---|
| `cargo test --lib` | **4,890 passed / 0 failed / 14 ignored** (4,874 before this step) |
| `cargo test --test ng_calling_loop_calls_genotypes` | **14 passed** (10 before) |
| `cargo test --test ng_calling_loop_allocation --features dhat-heap` | 1 passed |
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo doc --no-deps --lib` | 27 unresolved links, all pre-existing — none added |

**The release-held assertion battery ran twice**, once for the checks this step's first half added
and once for the two the reviews asked for. Every new `assert!` outside a test module was
downgraded to `debug_assert!` in one run and `cargo test --release --lib ng::calling
--all-features` was run: **five checks, five failing tests, one per check**. The source was
restored and the restore verified before moving on.

## 7. Banked for the owner

- **The tract's site quality is folded against the ordinary-site prior seed.** Pre-existing and
  named in `calling_quality.md` §8 as a gap; this step is what makes it reachable. Everything else
  on the tract path now reads a tract quantity.
- **A tract's warrant fold runs over every library of the run**, not the ones whose reads reached
  it — so a tract in a run of many libraries can be reported `Defaulted` on account of a library
  that sent it nothing. It is the conservative direction, it is now tested, and narrowing it means
  narrowing what the parameter table covers rather than changing the fold.
- **Nothing exercises a tract whose scratch rows and run samples differ**, because a tract rules no
  sample out, so no such tract exists. The map is read rather than assumed; the code comment says
  so, so a green suite is not read as evidence.
