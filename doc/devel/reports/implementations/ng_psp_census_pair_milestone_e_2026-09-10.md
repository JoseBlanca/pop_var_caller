# The census lives inside the psp — Milestone E: the sidecar's last traces

**Date:** 2026-09-10
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone E
**Spec:** [psp_census_pair.md](../../ng/spec/psp_census_pair.md) §3, §8, §11
**Branch:** `census-vs-psp-perf`, on top of `9fdb6d6a`
**Milestone D's report:** [ng_psp_census_pair_milestone_d_2026-09-10.md](ng_psp_census_pair_milestone_d_2026-09-10.md)
**Review:** [psp_census_pair_milestone_e_2026-09-10.md](../reviews/psp_census_pair_milestone_e_2026-09-10.md)

---

## What this milestone is for

Milestones A to D moved the census into the psp's trailer, taught the fit to read it there, and
built the repair for a psp whose census is stale. What is left is everything outside that path
that still spells the old shape: the scripts that ran the four commands, the machinery of the
sidecar that nothing calls any more, the two measuring probes, and the words — subcommand help,
`PROJECT_STATUS.md`, and the report for the plan as a whole.

**The baseline every step is judged against is the one Milestone C's report records** — the tree
is red in four ways, none of them this plan's — re-measured in each step's gate table.

---

## E1 — the scripts

**Committed:** see `git log` for `test(ng): E1`.

### What it does

**`scripts/ng_fit_stage_end_to_end.sh` was broken, not merely out of date.** It invoked
`generate-census`, a subcommand that no longer exists, told it `--output-dir`, and then checked
that a `.census` file sat beside every psp. Since A2 the walk writes no such file and since D1 the
command is `regenerate-census`, so the script failed on every run at its second step.

It now runs psp mode as it stands — three commands and a repair:

    1. generate-psps        the six CRAMs      ->  one .psp a sample, census sealed inside
    2. estimate-parameters  those psps         ->  cohort.parameters.toml
    3. call-from-psps       psps + that file   ->  two VCFs, one --defaults and one fitted
    4. regenerate-census    a copy of each psp ->  its trailer rebuilt, and the file whole again

**Step 4 is the plan's D4 oracle reaching real reads**, which is what its own line asks for. Each
walked psp is copied into a directory of its own, the copy's census is dropped, the repair rebuilds
it from the copy's own records, and the copy has to be the walked file again — header, blocks,
index, trailer, footer, every byte.

### The deviation: one new file the plan did not name

**`examples/ng_psp_drop_census.rs`**, thirty lines of program around one library call.

The oracle needs a psp that is *owed* a rebuild, and a copy of a walked psp is not one: a psp whose
census is the one this run would write is skipped and its records are never read (spec §8), so a
copy handed straight to the command comes back untouched and `cmp` passes without a rebuild having
happened. The plan's own note on D4 records this. Emptying a trailer is what makes the copy owed
one, and **nothing that ships does it**: the file has to be cut at the trailer's offset and a
footer re-encoded behind the cut with its trailer length zeroed, which is `replace_trailer`'s job
and not `dd`'s. `--force` on the repair would have been the other route and spec §8 rules it out.

So the example is `replace_trailer(psp, b"")` with a command line around it. It reads every
footer before it writes anything, so a run naming one unreadable file leaves every psp as it found
it; it refuses a psp whose trailer is already empty, because emptying it would arm nothing; and a
failure part-way names the psps it had already emptied, since `replace_trailer` can leave the one
it failed on with no footer at all.

### The oracle, and it reproduces the last recorded run

**On the same six accessions and the same ground as Milestone A** — the first six CRAMs of
`benchmarks/tomato1/crams` by name (`SRR7279481` to `486`, which the reads name `SRS3394606`,
`SRS3394711`, `SRS3394712`, `SRS3394712_SRR7279484`, `SRS3394713`, `SRS3394714`), over the first
two 100 kb intervals of `benchmarks/tomato1/regions.bed`, at about three reads a position.

| | Milestone A, 2026-09-09 | this run |
|---|---|---|
| census inside the psps | 1,545,479 bytes | **1,545,479** |
| parameters file | 38,124 bytes | **38,124** |
| records, `--defaults` | 2,275 | **2,275** |
| records, fitted | 2,082 | **2,082** |
| called only by the defaults | 196 | **196** |
| called only by the fit | 3 | **3** |
| of the records both called, differing in a genotype | 87 of 2,079 | **87 of 2,079** |
| genotypes differing | 113 of 12,474 | **113 of 12,474** |
| psps | 8,465,826 bytes | 8,465,862 |
| every copy identical to its walked psp, whole | — | **6 of 6** |

**That is stronger than "nothing moved".** Milestone A's figures were measured when this script
still called `generate-census` and the fit read census *files*; this run's fit reads each psp's
trailer instead. The parameters file coming back the same size and all five calling counters the
same across that change is C1 and C2's byte-identity re-confirmed end to end on real reads.

**What was compared is a size and five counts, not bytes.** The script does not `cmp` the
parameters file or the VCFs against Milestone A's, and could not: those artefacts were not kept.

**The one figure that moved is 36 bytes of psp, and it is the run's own command line.** A psp's
header records the argv verbatim (`current_command_line()` in `generate_psps.rs`), and the header
is text, so one extra character in the command line is one extra byte in the file. Measured, not
inferred: walking one accession into `tmp/e1_len_a` and again into `tmp/e1_len_aaaaaa` — paths five
characters apart — gave psps of 1,787,534 and 1,787,539 bytes. The gap here is 36 bytes over six
psps that share one command line, so six characters. **Milestone A's report does not record the
output path it used**, so which token was six characters shorter cannot be confirmed; what can be
said is that the census inside those psps is identical to the byte, so nothing about their content
moved.

### Two other scripts, and what was found in them

**`scripts/ng_census_route_cost.sh` and `examples/ng_census_route_cost.rs`** name the two ways to
end up with a census. Their prose said "both ship" and named `generate-census` as the second. Both
halves are now wrong, and in different ways:

- the second route is `regenerate-census`'s pass, which writes into the trailer rather than into a
  file of its own;
- and **no command produces the psp that route starts from** — one with an empty trailer.
  `generate-psps` always builds a census and spec §3.3 rules out a flag to skip it. The harness
  makes that state deliberately, in order to price what repairing such a psp costs.

Two claims went stale with them and are corrected: the census the harness writes is no longer "24
bytes a sample short of what `generate-census` writes", because nothing that ships writes a census
with a pileup identity — the two `write_census` calls outside tests, the walk's and the repair's,
both pass `None`. And the reason the after-the-walk route's census-file write stays inside the
clock has changed: it is no longer "work its real counterpart does" in the sense of writing a file,
but a write of the same census of the same size that `replace_trailer` makes into the psp.

**`scripts/ng_census_agreement_mutations.sh`** names neither the old command nor `--census`, so the
plan's E1 line left it alone. Reading it turned up something else: **its fourth case's stated
expectation was false.** It says the fourth defect — one read's minted error arriving one step
off — cannot be caught, because a census holds depth codes and allele counts and no per-read
quality. That stopped being true when the fit stage's Milestone C put the per-read-group minted
read-error totals into the census: `CensusWriter::add_locus` folds every complete observation's
`q_sum` into a per-read-group sum and `write_census` encodes it.

Measured rather than reasoned — the case was applied on its own and the tree restored afterwards:
**it fails two of the module's three tests**, the two that compare the two censuses. The third does
not fail, and not for the reason a first draft of this paragraph gave: it is built on the same
sample with reads as the first, and it never compares the two censuses at all — it checks the
fixture's walk reached its kept loci.

### What the review changed

One blocker and eleven smaller findings; the review report has them in full. The blocker:
**the script printed the run's headline answer and asserted nothing about it.** A fitted calling
pass that fell back to the compiled-in defaults, or two passes emitting header-only VCFs, both
printed four lines of zeros and exited 0 — and on this cohort every one of those counters has a
non-zero true value.

Two checks close it, and the split matters. The comparison now stops the run if either pass emitted
no record at all; it deliberately does **not** demand that the two VCFs differ, because at one
sample, or on ground where the fit lands near the defaults, calling twice and getting the same file
is a correct answer, and a script that required a difference would fail a good run. What catches a
fitted pass that used the defaults is a check on the invocation instead: each pass's report line
says *numbers behind the calls: N of 7 groups the file says were fitted*, and the defaults pass must
say 0 while the fitted pass must say more than 0.

### The mutations

Eleven cases, each applied on its own to a restored tree and run end to end over two accessions and
one 100 kb interval, with the diff printed at the apply step. The baseline over that small cohort
passes.

| the defect | what happened |
|---|---|
| the copies are never emptied | **caught** — the arming comparison stops the run |
| …and the arming check is disabled too | **survives**, and it is inherent — see below |
| the repair never runs | **caught** — every copy still carries an empty trailer |
| the repair is pointed at the originals | **caught** — the fresh originals are skipped, the copies untouched |
| both comparisons address the copy, so each is a file against itself | **caught** — the arming comparison then says "already identical" |
| the fitted pass quietly scores with the defaults | **caught** — by the fitted-groups check, which is why it exists |
| neither pass runs and both VCFs are empty files | **caught** |
| the fit never runs | **caught** — the parameters file is checked for content |
| the census total is read off the per-sample lines | **caught** — those lines are printed twice a sample, and the parse insists on one match |
| one psp is left uncopied | **caught** — the copies are counted against the psps |
| the arming tool reports success and empties nothing | **caught** — the arming comparison |

**One of these was a faulty mutation first time and the driver showed it.** The case meant to read
the census total off the per-sample lines dropped the pattern's leading anchor but kept the words
`bytes are census`, which those lines do not carry — they say `bytes are its census` — so it
matched exactly what it had matched before and the run passed. Re-aimed at the real words, it is
caught. Printing the apply step's diff is what made a mutation that changed nothing distinguishable
from one nothing catches.

**The surviving case is not a defect to fix.** If the copies are never emptied *and* the check that
would notice is disabled, the repair skips every psp, the copies stay byte-identical to the
originals, and the final comparison passes. The before-comparison *is* the guard, and a guard
cannot guard itself.

### The gate

**The four gates come back as the baseline set, item for item.**

| gate | at `9fdb6d6a` (baseline) | after E1 |
|---|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | 6,738 lib tests pass, 0 failed, 15 ignored; one target red, `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` | **the same** |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 11 errors of 5 kinds in 6 files | **the same 11, in the same 6 files** |
| `cargo check --all-targets --keep-going` | 4 examples do not compile | **the same 4** — the new one compiles |
| `cargo fmt --check` | 4 files | **the same 4** |

**One thing had to be put right to get there.** The first run of the gate listed a fifth file for
`cargo fmt`: the new example. It was formatted with `rustfmt` on its own — running `cargo fmt`
over the tree would have rewritten the four files that were already unformatted before this branch
— and the check is back to those four.
