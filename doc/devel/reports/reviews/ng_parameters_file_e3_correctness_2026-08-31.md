# E3 correctness review — the slippage slot, and what a run does without it

**What was reviewed.** The uncommitted step in `tmp/e3_step.patch`, applied to base commit
`6e434561` in the detached worktree `/Users/jose/devel/pop_var_caller-e3-rev1`. Four files:
`src/ng/calling/run_report.rs`, `src/ng/calling/run_parameters.rs`,
`src/ng/calling/parameters_file/to_toml.rs`, `src/ng/calling/parameters_file/defaults.rs`.

**Baseline.** `cargo test --lib ng::calling` in the worktree's own container: 1,026 passed, 0
failed, 2 ignored.

**⚑ First, a caveat that changes how to read the rest — see §0.** The patch is a stale snapshot:
the main repository's working tree carries a substantially rewritten E3, and two of the three Major
findings below are already closed there. Everything else still stands.

**Headline.** The shipped behaviour is correct in every state I could construct. The two
conditional notes fire in exactly the four states they should, the read-group walk is over the
right axis, and every number the code and the comments quote re-derives. What the step is short of
is **test coverage on the half of each new value that its fixtures never reach**: no fixture
anywhere builds a run that *fitted* a stratum and then asks the report about it, and no fixture
has one repeat-tract table empty while the other is full. Thirty-two mutations were run; twenty-two
died, ten survived, and six of the ten are real holes rather than equivalent mutants.

There are **no Blockers**. Nothing here can produce a wrong genotype: the report field has no
consumer outside tests yet (`run_report.rs:48` — *"The output stage that will print it is step 10's
and is not written, so today the only callers are tests"*), and the file notes are comments the
parser discards.

---

## 0. ⚑ Read this first — the patch is behind the working tree

**`tmp/e3_step.patch` is a stale snapshot of E3.** The main repository's working tree (branch
`ng-parameters-file`, HEAD `6e434561`, which *is* the base I was given) carries a substantially
rewritten E3, and I reviewed the patch, as instructed, without touching it. The differences are not
cosmetic:

| | in the patch I reviewed | in the working tree today |
|---|---|---|
| `no_stratum_was_fitted()` | returns one paragraph quoting the model's **four direction shares** | returns **two** paragraphs quoting the section's own three numbers — `share_of_reads_that_slip`, `shorter_share`, `fall_off` — inverted out of the shipped model, plus the part-repeat mass |
| an uncontaminated run | writes no `[contamination]` section and no explanation | writes a note in the section's place, because *"the one reader who most needed it … was the only reader who never saw it"* |
| files touched | four `.rs` | those four plus `testdata/every_shape_as_written.toml` (a **byte-compared golden output** file, `to_toml.rs:1139`) and `PROJECT_STATUS.md`, neither of which is in the patch |
| tests | six new | eight, including two that close two of my three Major findings |

**Two of my three Majors are already fixed there, independently.** The working tree adds
`a_run_that_fitted_a_stratum_does_not_report_every_tract_falling_back` (§3.1) and
`each_empty_table_note_answers_for_its_own_table` (§3.2), and the second one's doc records the same
measurement I made from the other side: *"pointing the substitution-rate note's guard at the
slippage table passed all 192 tests of this module, because every other fixture has both tables
empty or both full."*

**What still applies to the working tree**, checked by reading it (not by running it — I did not
build or mutate the main repo):

- **§3.3, the read-group walk's lower bound.** `run_parameters.rs` is byte-identical between the
  patch and the working tree, and the only fixture that exercises the walk still drops read group
  **2** alone. `(1..len)` would still survive.
- **§1.2 and §3.6**, `GroupNotInTheFit` and what `every_tract_falls_back()` promises. The word
  `GroupNotInTheFit` appears nowhere in `run_report.rs` or `run_parameters.rs`, and
  `every_tract_falls_back` is unchanged at `run_report.rs:316-318`.
- **§3.5**, the unchecked slippage-map keys, and **§3.7**, `strata()` counting strata present
  rather than strata with numbers — both in code the rewrite did not touch.

**What no longer applies**, because the rewrite already says it better: three of §3.8's four prose
points. *"There is nothing in this file to edit"* has become *"There is nothing here to change, only
rows to add … one written by hand is read back like any other"*; *"a whole repeat short"* is now
explicitly *"one repeat **or more**"*; and *"A PCR library slips more than this"* is now *"A PCR
preparation generally slips more than a PCR-free one, by an amount that depends on how many cycles
it ran."* The mutation results in §2 are results about the patch and would need re-running against
the working tree — in particular the rewritten note is built from four derived quantities rather
than four direct reads, which is new arithmetic and new mutation surface.

**So: everything below is a true report on `tmp/e3_step.patch`, and §3.1 and §3.2 are already
closed in the code that will actually be committed.** If a review of the working tree is wanted,
that is a second run.

---

## 1. Answers to the three questions the brief asked

### 1.1 Is `read_groups_with_no_slippage_group` walking the right axis? — Yes.

`report()` walks `0..calibration_by_read_group.len()` and maps each index to `ReadGroupId(index)`.
That is the run's own read-group axis, and it is the right one:

- `report()` opens by asserting `read_groups.len() == calibration_by_read_group.len()`
  (`run_parameters.rs:553-561`), so the walk's length is the read-group table's length;
- `ReadGroupId` **is** an index into that vector — `report()` itself indexes
  `self.contamination_by_read_group[read_group.get() as usize]` two dozen lines above
  (`run_parameters.rs:594`), and `gather_for_locus` mints the same identity
  (`inference/repeat_tract_parameters.rs:328-330`);
- the direction is the only one that can find an absence, exactly as the code comment says: asking
  the map for its keys can never name a key it lacks. Mutation **M10** (replace the walk with
  `Vec::new()`) and **M05** (`0..0`) both die.

**The reverse direction — a slippage map naming a read group the run does not have — is neither
reported nor refused.** I built it: `StratumFits::over(&[], {0→0, 1→0, 2→0, 7→3})` on a
three-read-group run assembles without complaint and reports
`read_groups_with_no_slippage_group=[]`. `RunParameters::assemble` performs precisely this check
for the *contamination* map (`run_parameters.rs:226-236`, *"an estimate past the axis is dropped
in silence, which leaves a contaminated library uncorrected"*) and performs none for the slippage
map. This is pre-existing rather than introduced by E3, and it is harmless today because the extra
key is never looked up — the lookups are keyed by the run's own read groups. See §3.5.

### 1.2 Is `NoSlippage::GroupNotInTheFit` reachable, and is it reported? — Reachable, not reported.

Run rather than reasoned. Building `StratumFits::of_gathered_rows` with the declaration
`{0→0, 1→0, 2→5}` and **one** stratum fitted over **one** slippage group:

```
at(ReadGroupId(2), period 2, 11 repeats)
    → Err(GroupNotInTheFit { slippage_group: 5, groups_fitted: 1 })
at(ReadGroupId(0), period 2, 11 repeats)  → Ok(...)
```

Every candidate of read group 2 at every stratum therefore falls back to
`StutterModel::hipstr_shipped`, and it falls back through the *"the run is not what it claims"*
door. The run-level report on that same run says:

```
strata_with_slippage = 1
fitted_substitution_rates = 0
read_groups_with_no_slippage_group = []
every_tract_falls_back() = false
```

— that is, nothing. `TractScoringFits` counts `UnknownReadGroup` and `GroupNotInTheFit` **together**
in `slippage_defaulted_by_an_unknown_read_group` (`repeat_tract_parameters.rs:340-346`), and the
new field's doc borrows that framing (*"which is why `NoSlippage` gives it a variant of its own and
the locus counts it apart from the ordinary absences"*) while covering only one of the two variants
the locus counts.

**Why this is Minor and not Major: neither production assembly path can produce it.**

- The parameters-file path densifies it away. `to_run_parameters` sizes every stratum's group
  vector to `max(declared slippage group, row slippage group) + 1`
  (`to_run_parameters.rs:293-305`), so a hand-edited file naming group 5 with rows only for group 0
  comes back as `GroupPutNoReadHere { slippage_group: 5 }` — an *ordinary* absence — not
  `GroupNotInTheFit`.
- The fit's own path is consistent by construction. `gather_strata` derives its group count from
  the same map, `max(group) + 1` (`ssr_fit.rs:2213-2217`), and `derive_thin_strata` only promotes
  an outcome once `furnished_any` is true (`ssr_fit.rs:1271-1274`).

So `GroupNotInTheFit` is defensive-only today, reachable through the public API by a caller that
assembles its own outcomes. Recording it in the report is cheap and would close the gap between
what the field's doc claims and what it covers.

### 1.3 Do the two notes fire exactly when they should? — Yes, in all four states.

I crossed `(slippage table empty?) × (substitution table empty?)` over
`a_file_using_every_shape()` with each table cleared in turn, and read the produced text:

| slippage table | substitution table | slippage note | substitution note |
|---|---|---|---|
| full | full | absent | absent |
| **empty** | full | **present** | absent |
| full | **empty** | absent | **present** |
| **empty** | **empty** | **present** | **present** |

That is the correct table. A partially-fitted run gets neither note; a run with slippage and no
substitution rates gets the second and not the first.

The produced text of a defaults run, read where it lands, is right: both notes sit immediately
above the tables they explain, wrapped to 80 columns, and the file still parses (the pre-existing
`a_defaults_run_writes_a_file_that_reads_back_as_the_same_run` covers that).

**But nothing in the suite separates the two conditions**, which is finding M-2 below.

---

## 2. Mutation testing

Thirty-two mutations across the new report field, the `report()` builder, the two conditional
notes and the derived numbers. Scored against **`cargo test --lib ng::calling`** — the author's own
scope, 1,026 tests. Test names are abbreviated; all live in
`src/ng/calling/parameters_file/defaults.rs` unless marked `to_toml`.

| # | mutation | verdict |
|---|---|---|
| M01 | `strata_with_slippage: 0` (constant) | **SURVIVED** — real hole, §3.1 |
| M02 | `strata_with_slippage: …strata_with_a_length_spectrum()` | **SURVIVED** — real hole, §3.1 |
| M03 | `fitted_substitution_rates: 0` (constant) | KILLED — `a_run_holding_substitution_rates_and_no_slippage_reports_both` |
| M04 | `fitted_substitution_rates: …strata()` | KILLED — same |
| M05 | read-group walk emptied, `(0..0)` | KILLED — `a_read_group_the_slippage_fit_does_not_name_is_named_in_the_report` |
| M06 | filter `is_none()` → `is_some()` | KILLED — `a_defaults_runs_report_says_every_tract_falls_back`, `a_read_group_…_is_named_in_the_report` |
| M07 | walk starts at 1, dropping read group 0 | **SURVIVED** — real hole, §3.3 |
| M08 | walk drops the last read group | KILLED — `a_read_group_…_is_named_in_the_report` |
| M09 | walk runs one past the axis | KILLED — both report tests |
| M10 | walk replaced by `Vec::new()` | KILLED — `a_read_group_…_is_named_in_the_report` |
| M11 | `every_tract_falls_back`: `== 0` → `> 0` | KILLED — all three report tests |
| M12 | `every_tract_falls_back` → `true` (constant) | **SURVIVED** — real hole, §3.1 |
| M13 | `every_tract_falls_back` reads `fitted_substitution_rates` | KILLED — `a_run_holding_substitution_rates_and_no_slippage_reports_both` |
| M14 | `every_tract_falls_back` → `read_groups_with_no_slippage_group.is_empty()` | KILLED — `a_read_group_…_is_named_in_the_report` |
| M15 | slippage note never fires | KILLED — `a_defaults_runs_file_says_what_its_empty_repeat_tract_tables_mean` |
| M16 | slippage note always fires | KILLED — `a_fitted_runs_file_carries_no_such_note`, `to_toml::the_whole_shape_writes_the_documented_toml` |
| M17 | slippage note keyed to the **substitution** table's emptiness | **SURVIVED** — real hole, §3.2 |
| M18 | substitution note never fires | KILLED — `a_defaults_runs_file_says_what_its_empty_repeat_tract_tables_mean` |
| M19 | substitution note always fires | KILLED — `a_fitted_runs_file_carries_no_such_note`, `to_toml::the_whole_shape_writes_the_documented_toml` |
| M20 | substitution note keyed to the **slippage** table's emptiness | **SURVIVED** — real hole, §3.2 |
| M21 | `in_a_hundred` drops the `× 100` | KILLED — `a_defaults_runs_file_says_what_its_empty_repeat_tract_tables_mean` |
| M22 | `in_a_hundred` scales by 10 | KILLED — same |
| M23 | whole-repeat share replaced by the part-repeat share | KILLED — same |
| M24 | whole-repeat shorter ↔ longer exchanged | **SURVIVED — provably equivalent**, §2.1 |
| M25 | `hipstr_shipped()` → `hipstr_em_start()` | KILLED — same |
| M26 | substitution note quotes a literal `0.01` instead of the constant | KILLED — same |
| M27 | the two note *bodies* exchanged | KILLED — same |
| M28 | slippage note wrapped at `COMMENT_WIDTH` not `ROOM_AT_THE_MARGIN` | **SURVIVED** — real hole, §3.4 |
| M29 | slippage note emitted *below* its table instead of above | **SURVIVED** — real hole, §3.4 |
| M30 | both counts `+ 1` | KILLED — all three report tests |
| M31 | unnamed read groups returned in descending order | **SURVIVED** — §3.6 |
| M32 | note's gloss changed to "one base in a **hundred**" beside `0.001` | **SURVIVED** — real hole, §3.4 |

**22 killed, 10 survived.** One survivor is a provably equivalent mutant; the other nine are
coverage gaps, of which six matter.

**The six survivors that matter were re-run against the whole library suite**, not just
`ng::calling` — `cargo test --lib`, 5,576 tests — in case some fixture elsewhere in the tree
happened to cover them. It does not: M01, M07, M12, M17, M20 and M28 survive the full suite too.

### 2.1 The one equivalent mutant, with its proof

**M24** exchanges `shipped.whole_repeat_shorter_share()` and `shipped.whole_repeat_longer_share()`
in the format arguments. `StutterModel::hipstr_shipped` sets **both** to `0.05`
(`src/ng/alignment/stutter.rs:310-311`), and `in_a_hundred` is a pure function of the share, so the
two arguments evaluate to the same `5u32` and the formatted string is byte-identical. The mutant
cannot be killed by any test that reads the produced text; killing it would need a test asserting
the *call sites*, which is not worth a test. The same argument applies to exchanging the two
part-repeat shares (both `0.01`).

That the file cannot be wrong *about which direction is which* is not an accident here — it is the
consequence of the shipped row being symmetric, which the note itself is at pains to point out.

---

## 3. Findings

### 3.1 MAJOR — nothing reads `strata_with_slippage` or `every_tract_falls_back()` on a run that fitted a stratum

**Where:** `src/ng/calling/run_parameters.rs:615`, `src/ng/calling/run_report.rs:307`.

`strata_with_slippage: 0` (M01) and `every_tract_falls_back() { true }` (M12) both survive the whole
of `ng::calling`. Every fixture that reads either builds its `StratumFits` with
`StratumFits::over(&[], …)` — no strata — so the *false* branch of the type's headline question is
never exercised.

This is **the same hole the author already found and fixed once, one field over**. The doc comment
on `a_run_holding_substitution_rates_and_no_slippage_reports_both` (`defaults.rs:1217-1219`) says:
*"reporting a constant zero there passed every other test in `ng::calling`, because every other
fixture that reads this field is a run with nothing fitted at all."* That is now exactly true of
`strata_with_slippage` — and that test's own fixture still passes `StratumFits::over(&[], …)`, so it
did not close the sibling case.

**Concrete input that produces the wrong outcome.** A run whose `StratumFits` holds one stratum:

```rust
let fits = StratumFits::of_gathered_rows(
    BTreeMap::from([(ReadGroupId(0), 0), (ReadGroupId(1), 0), (ReadGroupId(2), 0)]),
    BTreeMap::from([(Stratum { period: 2, reference_repeats: 11 }, vec![Some(fitted)])]),
    BTreeMap::new(), BTreeMap::new(),
    STATED_FLAT_CONCENTRATION, Provenance::Defaulted,
);
```

Correct: `strata_with_slippage == 1`, `every_tract_falls_back() == false`. Under M01 or M12 the
run's one statement about repeat tracts reads *"every repeat tract of this run is scored under
another caller's shipped constants"* on a run that fitted its own, and no test objects. Since the
whole purpose of the field is to distinguish those two runs, the untested half is the half that
matters.

**Closing it costs one test**: build a run with one stratum, assert `strata_with_slippage == 1` and
`!every_tract_falls_back()`. That single test kills M01, M02 and M12.

### 3.2 MAJOR — nothing separates the two notes' conditions

**Where:** `src/ng/calling/parameters_file/to_toml.rs:320` and `:363`.

M17 (the slippage note keyed to `substitution_rate_by_stratum.is_empty()`) and M20 (the mirror)
both survive. The reason is that the suite has only two repeat-tract fixtures and both tables move
together in each: `a_file_using_every_shape()` has both full, and a defaults run has both empty. No
test ever puts a file in a state where the two conditions differ.

**Concrete input that produces the wrong outcome under M17.** A run that fitted slippage and no
substitution rates — a state the patch's own report test
(`a_run_holding_substitution_rates_and_no_slippage_reports_both`) establishes is reachable, and
which the section's own prose says is normal because *"the two are fitted separately — slippage per
`(stratum × slippage group)`, the rate per `(read group × stratum × ploidy)`"*. Take
`a_file_using_every_shape()`, clear only `substitution_rate_by_stratum`, and write it. The mutated
writer prints, immediately above three fitted slippage rows:

> This table is empty, which is not the same as a missing row: no stratum was fitted at all, so
> every repeat tract in this run was scored under the stutter model this caller ships …

That is a false statement about the reads, in the user-facing file, in the exact register the note
exists to correct. The shipped code does **not** do this — I verified all four states — but the
suite would not notice if a future edit made it.

**Closing it costs one test**: the four-state table in §1.3, asserting each note's presence equals
its own table's emptiness. That single test kills M17 and M20 (and, incidentally, M15, M16, M18 and
M19 as well).

### 3.3 MAJOR — the read-group walk's lower bound is unguarded

**Where:** `src/ng/calling/run_parameters.rs:620`.

M07 changes `(0..calibration_by_read_group.len())` to `(1..…)` and survives. The only fixture that
exercises the walk's boundaries — `a_read_group_the_slippage_fit_does_not_name_is_named_in_the_report`
— drops read group **2**, the last one, so the top of the range is pinned (M08 and M09 both die)
and the bottom is not.

**Concrete input that produces the wrong outcome.** The patch's own fixture with the declaration
`{ReadGroupId(1): 0, ReadGroupId(2): 0}` instead of `{ReadGroupId(0): 0, ReadGroupId(1): 0}` — that
is, the *first* library is the one the pre-pass never saw. Correct:
`read_groups_with_no_slippage_group == [ReadGroupId(0)]`. Under M07: `[]`, and the run reports
nothing about a library whose every tract is scored under HipSTR's constants through the
*"the run is not what it claims"* door.

Adding `ReadGroupId(0)` to the existing fixture's dropped set — expecting
`[ReadGroupId(0), ReadGroupId(2)]` — kills M07 and M31 at once.

### 3.4 MINOR — three properties of the produced notes that no test pins

All three survive because the only test that reads a defaults run's file uses `contains` on seven
phrases, and the only test that measures comment geometry reads a *different* file.

- **M28 — wrap width.** `no_comment_line_is_longer_than_the_prose_it_carries` (`to_toml.rs:1759`)
  runs over `a_file_using_every_shape().to_toml()`, which by construction contains neither new note.
  So `wrapped(&no_stratum_was_fitted(), COMMENT_WIDTH)` — 80 characters of text plus the `# `,
  giving 82-column lines — passes. Fix: run that same loop over a defaults run's file too.
- **M29 — placement.** Moving the note *below* the table it explains survives, because
  `unwrapped_comments` flattens the whole file into one string and the assertions are `contains`.
  Placement is the entire point of the note (the section's own paragraph is above the table and is
  what the reader misread). Fix: assert the note's first sentence appears before the line
  `slippage_by_stratum_and_group = [`.
- **M32 — the gloss.** The test pins `"the caller's stated 0.001"` but not the phrase beside it, so
  the note can read `0.001 — about one base in a hundred read wrong inside a tract` and pass. The
  gloss is the half a geneticist actually reads. Fix: add `"about one base in a thousand"` to the
  owed-phrases list.

### 3.5 MINOR — a slippage map naming read groups the run does not have is neither reported nor refused

**Where:** `src/ng/calling/run_parameters.rs:329` (`of_gathered_values`) and `:183` (`assemble`).

Measured: `StratumFits::over(&[], {0→0, 1→0, 2→0, 7→3})` on a three-read-group run assembles, calls
and reports `read_groups_with_no_slippage_group=[]`. Nothing anywhere objects. `assemble` runs the
identical check for the contamination map and refuses (`run_parameters.rs:226-236`).

Harmless today — the extra key is never looked up, because lookups are keyed by the run's own read
groups — and **pre-existing rather than introduced by E3**. Listed because the brief asked, and
because it is the other half of the evidence that "the parameters and the reads came from different
runs", which is the thing the new field exists to say.

### 3.6 MINOR — `every_tract_falls_back()` is false on a run where every tract does fall back

**Where:** `src/ng/calling/run_report.rs:298-307`.

The method's first sentence promises *"Whether every repeat tract of this run is scored under the
shipped stutter model"*; its second gives the mechanism, *"true exactly where no stratum carries
numbers"*. The two disagree whenever the fall-back is caused by an unnamed read group rather than
by an absent stratum. Measured, on a run with one fitted stratum and an empty declaration map:

```
strata_with_slippage = 1
read_groups_with_no_slippage_group = [ReadGroupId(0), ReadGroupId(1), ReadGroupId(2)]
every_tract_falls_back() = false
```

while `at()` answers `UnknownReadGroup` for every read group at every stratum — every tract does in
fact fall back. The same holds for the `GroupNotInTheFit` state of §1.2, where the sibling field is
empty too and *nothing* in the report shows the fall-back.

The information is not lost — a consumer reading both fields can work it out — but the method that
carries the name of the question does not answer it. Either narrow the first sentence to match the
mechanism ("whether the run fitted any slippage at all"), or widen the predicate.

**No consumer exists yet**, which is why this is Minor: `run_report.rs:48` records that the output
stage which will print this is step 10's and unwritten.

### 3.7 MINOR — `strata_with_slippage` counts strata *present*, not strata *with numbers*

**Where:** `src/ng/calling/run_report.rs:280` (field doc) and `run_parameters.rs:615`.

The field is documented as *"How many strata carry slippage numbers"* and is filled from
`StratumFits::strata()`, which is `by_stratum.len()` — how many strata are in the fit, whether or
not any slippage group in them has numbers. The file's table, by contrast, is written from
`each_stratum_and_group_with_numbers()`, which skips the empty cells
(`from_run_parameters.rs:485-486`). So the note's condition (an empty *table*) and the report's
condition (`strata() == 0`) are two different predicates for one claim, computed in two places, with
no test crossing them.

**In practice they cannot differ**, and the proof is worth recording: `derive_thin_strata` only
promotes a `Refused` outcome once `furnished_any` is true (`ssr_fit.rs:1271-1274`), and a `Fitted`
outcome exists because reads were fitted, so every stratum that reaches `by_stratum` has at least
one group with numbers. Worth one clause in the field's doc rather than a code change.

### 3.8 MINOR — four claims in the user-facing note

The note is prose a geneticist reads and argues with, so `CLAUDE.md`'s writing rules apply to it.

- **"There is nothing in this file to edit"** is not true of the table it sits above.
  `slippage_by_stratum_and_group` is read back by `to_run_parameters` and `validate` accepts
  hand-written rows (`validate.rs:714-745` checks the shares and the fall-off and nothing about
  where they came from). Supplying slippage rows by hand is exactly what a user with no fit would
  reach for, and the note tells them there is no such thing. Whether the project *wants* that is a
  design question; the sentence asserts a property of the file that the file does not have.
- **"A PCR library slips more than this."** No size, no subject, no measure — the pattern
  `CLAUDE.md` names first. There is a number in the tree that would carry the *shape* half of the
  claim: `Slippage`'s own doc records tomato dinucleotides at a shorter-share near 0.83, *"2,438
  shorter against 501 longer"* (`ssr_fit.rs:86-88`). Nothing in the tree measures the magnitude
  half.
- **"5 reads in 100 report a whole repeat short"** is the total mass of the whole-repeat-contraction
  branch across *all* step sizes, not the one-repeat slip the phrase describes. `StutterModel::new`
  derives `same_length_share` as one minus the four direction shares (`stutter.rs:282-287`) and
  `probability` distributes each direction's share over step sizes through a geometric
  (`stutter.rs:378-388`), so exactly one repeat short is `0.05 × 0.95 = 4.75` in 100. The same
  wording is already in `defaults.rs`'s module header, so the note is consistent with the tree —
  correcting it means correcting both.
- The type doc's *"a `(read group, stratum)` the substitution-rate fit never accumulated"*
  (`run_report.rs:257`) names two keys where `StratumKey` has three; the field doc four lines
  below says `(read group × stratum × ploidy)` correctly.

### 3.9 MINOR — a test's name and doc claim more than the test does

**Where:** `src/ng/calling/parameters_file/defaults.rs:1070`.

`a_defaults_runs_tract_is_scored_under_the_shipped_model_and_counted` neither scores a tract nor
counts anything. Its doc says *"This goes through the caller's own assembly, so it is the same
lookup a locus makes: `TractScoringFits::gather_for_locus` asks `StratumFits::at` …"* — it does not;
it calls `StratumFits::at` directly and then asserts four constants of
`StutterModel::hipstr_shipped`. The body is a superset of the pre-existing
`a_defaults_runs_tracts_find_no_stratum_rather_than_an_unknown_read_group`
(`defaults.rs:798`) over three candidate lengths instead of one, plus the constant pinning — which
is genuinely useful, because it is the anchor the derived note is checked against.

**The traced behaviour it claims is nonetheless real and is covered elsewhere**, so this is a naming
defect and not a coverage gap: `repeat_tract_parameters.rs:1304-1321` drives `gather_for_locus` and
asserts `*table.of(read_group, 1).stutter == StutterModel::hipstr_shipped()`,
`weakest_provenance == Provenance::Defaulted`, and
`cells_with_no_fitted_slippage() == READ_GROUPS`. The E3 brief's "say so and stop for a ruling"
condition — *does the traced behaviour score the tract?* — is answered yes, correctly, and the
owner's ruling stands.

### 3.10 Not a finding — the `as u32` cast, checked and dismissed

`run_parameters.rs:621` mints `ReadGroupId(group as u32)` with a plain cast, where
`gather_for_locus` uses the checked form and says why: *"Checked rather than cast: `group as u32`
would silently make read group 2³² into read group 0 and score a library against another's
polymerase"* (`repeat_tract_parameters.rs:325-330`). I was going to raise it, and then found that
the same file already casts this way twice, in `assemble` (`run_parameters.rs:209` and `:239`), both
predating this step. The new line matches its own file's convention, the value is unreachable at
4.3 billion read groups, and clippy is silent because `cast_possible_truncation` is not in
`Cargo.toml`'s lint table. Recorded so the next reviewer does not spend the same twenty minutes.

---

## 4. Every number the patch asserts, re-derived

| claim | where | check |
|---|---|---|
| `whole_repeat_shorter_share() == 0.05` | `defaults.rs:1091` | `stutter.rs:311` ✓ |
| `whole_repeat_longer_share() == 0.05` | `defaults.rs:1092` | `stutter.rs:310` ✓ |
| `part_repeat_shorter_share() == 0.01` | `defaults.rs:1093` | `stutter.rs:314` ✓ |
| `part_repeat_longer_share() == 0.01` | `defaults.rs:1094` | `stutter.rs:313` ✓ |
| the four survive `StutterModel::new` | — | `sanitize_direction_share` clamps to `[0,1]`; 0.05 and 0.01 pass unchanged (`stutter.rs:853`) ✓ |
| "5 … 100 … 5 … 100 … 1 … 1" in the note | `to_toml.rs:481-484` | `(0.05 × 100).round() = 5`, `(0.01 × 100).round() = 1` ✓ (see §3.8 for what "a whole repeat short" elides) |
| "one read in twenty … one in a hundred each way" | `run_report.rs:272-277` | 0.05 = 1/20, 0.01 = 1/100 ✓ |
| `DEFAULT_SSR_SUBSTITUTION_RATE` prints as `0.001` | `to_toml.rs:495`, `defaults.rs:1130` | `= DEFAULT_ERROR_RATE`; pinned at `repeat_tract_parameters.rs:2474` ✓ |
| "about one base in a thousand" | `to_toml.rs:495` | 0.001 ✓ |
| `fitted_substitution_rates == 2` | `defaults.rs:1263` | two distinct `StratumKey`s — `(rg 0, period 2 × 6, 2n)` and `(rg 2, period 3 × 9, 2n)` — so the `BTreeMap` has length 2 ✓ |
| `strata_with_slippage == 0` for a defaults run | `defaults.rs:1054` | `StratumFits::over(&[], …)` leaves `by_stratum` empty ✓ |
| `read_groups_with_no_slippage_group` empty for a defaults run | `defaults.rs:1060` | `of_defaults` declares all `read_groups.len()` groups into slippage group 0 (`defaults.rs:384-395`) ✓ |
| `at(rg, 2, r) == Err(NoSuchStratum)` for 9 (group, length) pairs | `defaults.rs:1078-1086` | `at` resolves the read group first (all declared), then misses `by_stratum` ✓ |
| `at(ReadGroupId(2), 2, 11) == Err(UnknownReadGroup)` | `defaults.rs:1206` | read group 2 absent from `{0→0, 1→0}`; the read-group lookup precedes the stratum lookup (`stratum_fits.rs:770-773`) ✓ |
| every shipped share is a whole percent | `to_toml.rs:1902-1917` | 0.05 × 100 = 5.000000000000001, 0.01 × 100 = 1.0000000000000002 — both within 1e-9 of an integer; the guard is not vacuous, it fails at e.g. 0.035 ✓ |
| `ROOM_AT_THE_MARGIN` gives ≤ 80-column comments | `to_toml.rs:457` | 78 + `"# "` = 80 ✓ (but see M28 — unguarded for these notes) |

---

## 5. What I would ask for before this lands

Of the three tests the patch needed, **the working tree already has two** (§0): a run with one
fitted stratum asked for its report, and the four-state note table. The third is still open:

1. ~~A run with one fitted stratum, asked for its report~~ — **already added** as
   `a_run_that_fitted_a_stratum_does_not_report_every_tract_falling_back`. Would have killed M01,
   M02, M12.
2. ~~The four-state note table of §1.3~~ — **already added** as
   `each_empty_table_note_answers_for_its_own_table`. Would have killed M17 and M20.
3. **Still open: the unnamed-read-group fixture extended to drop read group 0 as well as 2** —
   asserting `[ReadGroupId(0), ReadGroupId(2)]`. Kills M07 and M31, and the code it guards is
   byte-identical between the patch and the working tree.

And three one-line additions to the existing tests: run the comment-width loop over a defaults run's
file (M28); assert the note precedes `slippage_by_stratum_and_group = [` (M29); add
`"about one base in a thousand"` to the owed phrases (M32).

The design questions — whether `GroupNotInTheFit` should reach the run report (§1.2), whether
`every_tract_falls_back()` should mean what its name says (§3.6), and whether the note should keep
telling the reader there is nothing here to edit (§3.8) — are the owner's, not mine.

---

## 6. How this was run

- Worktree `/Users/jose/devel/pop_var_caller-e3-rev1`, patch applied with `git apply`, base
  `6e434561`.
- Every cargo invocation through the **worktree's own** `scripts/dev.sh`, each run confirmed to
  print `CWD /Users/jose/devel/pop_var_caller-e3-rev1` and, after each mutation,
  `Compiling pop_var_caller v0.1.0 (/Users/jose/devel/pop_var_caller-e3-rev1)`.
- Mutations applied by exact-string replacement with a uniqueness check on the anchor (a mutation
  whose anchor matched zero or two places was reported rather than run); the four touched files
  restored from a snapshot between every mutation, and verified byte-identical to the snapshot at
  the end.
- Probe tests and the mutation driver live in the worktree's `tmp/e3rev/`; the probes were removed
  from the tree before mutation scoring, so every verdict above is against the author's tests alone.
