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

---

## E2 — the sidecar's machinery deleted

**Committed:** see `git log` for `refactor(ng): E2`.

### What it does

The sidecar was the `<sample>.census` file that used to sit beside each psp. Since Milestone A the
census is the psp's own trailer, so the two cannot come apart, and everything that existed to check
that pairing has nothing left to check. Gone: `PileupIdentity`, the `Freshness` verdict and the two
functions that produced it, `CENSUS_FILE_EXTENSION` and `census_path_for`, and three tests asserting
that a trailer's census names no psp — a property the type system now holds, there being nothing
left to name one with.

**It is a refactor and not a deletion, because `PileupIdentity` was not dead code.**
`regenerate-census` reports each sample by the record count that type carried, so `census_from_psp`
returns a plain `records: u64` instead; the header digest beside it had no reader outside tests.
`write_census` and both census readers lose the argument every shipped caller passed `None` for.

### The census format did not change, and no psp was invalidated

**The identity occupied one flag byte** saying present or absent, and every shipped writer writes
absent. That byte is still written, always zero, so every census this build has written is byte for
byte what it was and `VERSION` stays at 4. Removing it would have moved every field of the
directory, cost a format version, and made every psp already on disk unreadable — to save one byte
a sample.

**What a census with that byte set now gets is a refusal**, and that is a change from the first
draft of this step, which stepped over the 24 bytes behind it. The review's argument is the one
that decided it: such a census names a psp by a digest and a record count that this build no longer
compares against anything, so reading past them is accepting a file it can say nothing true about.
Nothing in the shipped commands can produce one — the fit and the repair both read a psp's trailer,
and the walk has written the byte zero since Milestone A — so the refusal is reachable only from a
census file kept from an older build.

### Three things went newly unused and were not deleted

`psp::header_digest`, `psp::header_and_its_digest` and `WriteStats::header_digest`. Their only
reader was the identity. They are left in place with ⚠ notes because they live in the psp writer
rather than in the sidecar, and dropping the last of them changes what the walk reports about
itself — **a decision for the owner at Checkpoint E**, where it is listed.

### The decision the plan asked for, either way

`every_census_in_the_cohorts_psps` is **kept**. Production calls the two halves so that a command
can judge each census between them; this is the whole-cohort read that the tests of that assembly
are written against, and deleting it would put the composition by hand into five test call sites in
three modules. The reason is now in its own doc comment, so the next reader does not hunt for a
production caller that was never there.

### What the review changed

One reading of the format, six sentences that survived the deletion and had become false, and one
gap in the tests. No finding was a wrong deletion: the reviewer checked the pre-change tree and
every reader of every deleted item was a test.

- **The unreachable branch.** The step-over arm could not be exercised by any writer in the tree,
  and the test that had covered the format's two flag values went with the type. Refusing instead
  makes the behaviour reachable from a hand-built census, and there is now a test that builds one:
  it finds the flag byte from `encode_header`'s own output rather than guessing an offset, so it
  cannot drift from the layout it pokes at.
- **The version had no test behind it.** Every other test compares a file's version word against
  `VERSION`, so all of them pass whatever it holds. A second constant beside it, and a test that
  the two agree, makes a bump a deliberate edit in two places.
- **Six false sentences**, in `gatherer.rs` (two), `generate_psps/tests.rs`,
  `regenerate_census/tests.rs`, `regenerate_census.rs` and `census_from_psp.rs`: each promised a
  check that this step deleted, in the present tense. Three more outside those files — the psp
  header's read-filter keys were justified by the census naming its psp, the writer line still said
  the census was encoded twice, and the route-cost harness's constant command line rested on the
  same mechanism. All rewritten around the reason that still holds, which in the first case is that
  a cohort's psps are refused unless their headers agree.

**Two documents outside this plan are now false and were left alone**: `doc/devel/ng/arch/run_streaming.md`
cites two deleted items by line, and `doc/devel/ng/spec/run_streaming.md` §6.1 justifies the psp
header's contents by a consumer that no longer exists. Editing another plan's spec is not this
plan's to do; both are listed at Checkpoint E.

### The mutations

Nine cases, each applied on its own to a restored tree, with the diff printed at the apply step.
The baseline over the census tests passes with 204.

| the defect | what happened |
|---|---|
| the header's flag byte is written present, with nothing behind it | **caught** |
| the flag byte is not written at all, so every later field lands one byte early | **caught** |
| a census that names a psp is read past instead of refused | **caught** — by the test added for it |
| the format version moves | **caught** — by the test added for it |
| the record count starts at one | **caught** |
| the count is of bodies decoded rather than of records read | **survives** |
| the repair reports no records for any sample | **caught** |
| the walk seals every psp with an empty trailer | **caught** |
| a census is read from one byte past its trailer | **caught** |

**The survivor is narrow and named rather than fixed.** The walk this producer uses always hands
over a body, so counting bodies and counting records give the same number on every fixture; the two
would part company only under a selective walk, which nothing asks for here. **And one of the nine
did not compile first time** — the driver printed the compiler's own output beside the result, which
is what kept it from being recorded as a defect nothing catches.

### The gate

| gate | at `9fdb6d6a` (baseline) | after E2 |
|---|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | 6,738 lib tests pass; one target red | **6,737 pass**, the same target red |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 11 errors of 5 kinds in 6 files | **the same 11, in the same 6 files** |
| `cargo check --all-targets --keep-going` | 4 examples do not compile | **the same 4** |
| `cargo fmt --check` | 4 files | **the same 4** |

**The one lib test fewer is arithmetic, not a loss**: three deleted with the type they exercised,
two added by the review's fixes.

---

## E3 — the two measuring programs checked

**Committed:** see `git log` for `docs(ng): E3+E4`, which carries E4 as well.

### What this step is, and is not

Both programs were already in `examples/` — they landed on this branch before the plan's first
commit — so E3 is **verification**, not landing. What it had to establish: that they still build
after E2's deletions, that they still run, and that the one which writes a census writes it
somewhere of its own and leaves the psp it was given alone.

### What was measured

**Both build**, and `cargo check --all-targets` lists the same four examples that did not compile
before this branch and no fifth.

**`ng_census_read_vs_psp` answers its question on a real psp.** On `SRS3394712.psp` — one tomato
accession over the two 100 kb intervals, 1,190,325 bytes holding 193,603 records — at a census
budget of one position in one:

| what was timed | one round |
|---|---|
| decode the census whole | **0.0003 s** |
| rebuild it from the psp, every record's body | 0.0440 s |
| rebuild it, declining the bodies a census does not need | 0.0494 s |
| every head, no body — the floor a psp route cannot go below | 0.0073 s |

**One round, so read the ratio loosely**: a second run of the same command gave 0.0440 against
0.0436 for the full-body walk, and the ratio it prints moved from 152 to 127. What the numbers
support is *two orders of magnitude*, not a figure.

**Declining bodies did not save anything here, and the program says why**: at this budget the
selection keeps 198,182 of 200,000 positions, so every record is at a kept position and the
selective walk built all 193,603 bodies — 1 in 1. The saving that arm exists to measure needs a
budget that keeps a small share of the ground, which is what a whole-genome run has and a 200 kb
run does not.

**It wrote its census to `tmp/ng_census_read_vs_psp/SRS3394712.k1.census`, 240,827 bytes** — the
same size `regenerate-census` wrote into that psp's trailer in E1's run — and the psp it was given
was not touched.

**`ng_census_locus_spans` answers its own.** Of the same psp's 193,603 records, **198 are wider
than one base** (1 in 1,000, the widest 42 bases), and of the 5,978 kept positions carrying
evidence that disagrees with the reference, **146 have that evidence dropped** because the record
covering them spans more than the base it is recorded at. That is the price of the census's
per-position rule, on this ground: about 1 kept position in 40 of those that carry non-reference
evidence.

### What changed, and what deliberately did not

**`ng_census_read_vs_psp`'s module doc opened on a false sentence** — *"The census beside each psp
is a cache"* — and its first repair was false in a subtler way: it said the harness writes the
psp's trailer bytes out, which it does not. It rebuilds the census with the shipped producer and
writes that. The two are the same bytes only at a budget of one in one, which is what the
walk-versus-rebuild oracle guarantees; at any other budget the census this program must time is
one no file on disk holds. That, and not tidiness, is why it writes a file of its own, and the doc
now says so.

**`ng_census_locus_spans` was not touched.** Nothing in it names a census file, a sidecar, an
identity or a freshness check, so E2 falsified nothing in it — and it is one of the four files
`cargo fmt --check` listed before this branch, so leaving it alone keeps that set at four.

---

## E4 — the words

**Committed:** with E3, above.

### What changed

**`generate-psps`'s help said a census goes beside each psp.** It now says a sample is one file and
the census is sealed into that psp's tail. The other three subcommand docs were already current
from D1; `regenerate-census`'s list of what it repairs was one cause short of the four the verdict
type has, and now names all four.

**`PROJECT_STATUS.md` has a new entry at the head of its current-focus block** — the pipeline as it
stands, what the plan makes impossible, and the figures E1 measured.

**The report for the plan as a whole** is
[`ng_psp_census_pair_2026-09-10.md`](ng_psp_census_pair_2026-09-10.md). It is written for someone
who did not follow the plan: what changed for a person running the caller, what it cost, what was
decided along the way, and what is left open.

### The judgement the plan asked for and this did not do

**`PROJECT_STATUS.md`'s stale pipeline line was not rewritten.** The plan's E4 names it, but it
sits inside a dated entry of 2026-09-05 describing a run that really did go
`generate-psps → generate-census → estimate-parameters → call-from-psps`. Rewriting it would
falsify a record of what happened. What went in instead is a parenthetical in the file's own house
style — it has two other entries marked as superseded the same way — pointing at the entry at the
head of the block. The review agreed the entry should stand and asked for exactly that marker,
because an arrow diagram is the most copyable thing on the page and carries no date inside itself.

### What the review changed

Three blockers, all of them claims in the whole-plan report that the code contradicts, and all
three were things I had reasoned rather than read:

- **"the reference is never opened" was inverted.** Two refusals happen at different moments: a psp
  with no census, or one of a format this build does not read, is caught on ten bytes and a seek
  **before** the reference is opened; a census recorded under other settings cannot be caught that
  cheaply, because the settings are a digest over a selection that has to be rebuilt from the
  reference first. The report had merged them into one refusal that costs nothing.
- **The read filters are not compared by any cohort opener.** `SegmentationInputs::first_difference`
  compares the catalog, the repeat-tract criteria and the analysed regions, and nothing else. The
  filters are recorded in every psp's header and the calling run's report names the ones that
  disagree — no command refuses over them. **This one had reached committed code**: E2's own review
  fix replaced a stale reason in `ReadFilterConfig`'s doc with a false one. It is corrected in this
  commit rather than by amending E2's.
- **The four-commands-to-three story was wrong about which work moved.** The walk already built the
  census in the pass it was making anyway; `generate-census` built a *second* copy from the stored
  psp, and an end-to-end run compared the two. What this plan removed is that second copy, not the
  work — which is also why the two-producer agreement survives as this plan's parity oracle rather
  than being something it had to invent.

And six smaller ones: a measurement quoted as six characters when five were measured, headline
figures with no source, a spec section named as stale in one place when it is stale in four, the
twelve settings split seven-and-five when the code says nine-and-three, "no reader left" where one
test still reads, and milestone letters doing work in a document written for someone who has never
seen them.

**`PROJECT_STATUS`'s new entry had three of its own**: it said a pair coming apart was invisible
*and* that there was a check for it, it pointed at "the last recorded run" when the last run
recorded in that file is a different one, and it said nothing a run produces moved while omitting
the 36 bytes of psp that did.

### The gate

E3 and E4 change no behaviour. The gate is the baseline's, item for item, with the same lib-test
count as E2's — 6,737 — and `cargo fmt --check` still listing the same four files, which is what
leaving `ng_census_locus_spans.rs` alone preserves.
