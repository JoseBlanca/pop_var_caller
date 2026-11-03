# ng — fitting the misread-after-a-run table in the parameter pre-pass

**Approved by the owner on 2026-11-02 and built on branch `str-context-error`; narrowed to runs of one
base on 2026-11-03 (§1).** The design this implements is
[`str_context_error.md`](str_context_error.md): a read that has just been sequenced through a long
run of one base misreads the bases after it, most often reading the next base as one more base of
the run, and the STR path's substitution term has to price that per read direction. That document's
prototype hard-codes one table measured on HG002; this one says how each run fits the table itself.

Reads on: [`str_context_error.md`](str_context_error.md) §2–§3 (the effect and the decision),
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.1 (the STR substitution rate this adds to),
[`read_likelihoods.md`](read_likelihoods.md) §4.3 (where the rate is used), and
[`../../specs/design_principles.md`](../../specs/design_principles.md) principle 0 (the range).

## 1. What this does

**For each read group, the pre-pass counts how often reads misread a base just after a long run of
one base, separately for reads that were sequenced through the run and for reads that reached the
base from the other side, and `estimate-parameters` turns the difference into a table the STR path
charges.** The difference is the point. A true variant at the base shows in both groups of reads; a
misread caused by the run shows only in the reads that crossed it. So the table needs no truth set
and no called genotype, and a sample's own variants cancel out of it.

The owner's rulings this follows:

- fit the table per read group in the pre-pass, rather than hard-coding it (2026-11-01);
- **charge after runs of one base only, and do not count or report other repeats** (2026-11-03).
  The first build also counted and charged tandem repeats of 2 to 4 bases per unit, letting each read
  group's data decide (2026-11-02). Measured on HG002, that bought nothing: the fit with every such
  entry set to zero called within 6 calls of the fit with them, in all four runs at 30× and 5×
  (`../../reports/ng_str_context_error_fit_2026-11-02.md` §7), and the bases it charged there carried
  high base quality, which is what stutter written by the mapper as mismatches looks like, not a
  sequencing misread. Nor is there a reason to expect one: the misread comes from the sequencer
  losing phase over many identical incorporations in a row, which a run of one base causes and a
  repeat alternating between bases does not. So the table is 19 entries, all for runs of one base.

## 2. The table

### 2.1 What one entry of the table is

**An entry says, for reads that have just crossed a run of one base of a given length, how much more
often they misread a base at a given distance after it than reads arriving from the other side.**
The distance is whether the base is the next one after the run or one of the 19 after that.

At the next base there are two numbers, because the misread has a favourite outcome:

- **reading it as the run** — one more base of the run: after `A`×16, an `A`;
- **reading it as any other wrong base.**

Further on there is one number, any wrong base, because the measurement found no favourite there
(`str_context_error.md` §2).

Each number is an **excess**: the rate in reads through the run minus the rate in reads arriving
first, never below zero. It is added to the read group's fitted STR substitution rate, exactly as the
prototype adds its hard-coded excess.

### 2.2 The rows

| | run lengths, in bases |
|---|---|
| next base | each length from 8 to 22, then 23 or more |
| the 19 bases after | 8–11, 12–15, 16 or more |

That is 16 + 3 = 19 rows. **Runs shorter than 8 bases are not charged**: the measurements found the
two groups of reads within 1 in 1,000 of each other there (`str_context_error.md` §2), and §5 of this
document says what an uncharged base costs.

### 2.3 Which entry a base belongs to

**One rule decides a base's entry, for counting and for scoring alike**, so the pre-pass counts
exactly what calling later charges. It is a function of a sequence and a read direction:

1. Walk the sequence in the direction the read was sequenced: left to right for a forward read,
   right to left for a reverse one.
2. A **run** is a maximal stretch of one base, `A`, `C`, `G` or `T`. A stretch of `N` is no run: it
   is no sequence at all.
3. A base's candidate entries are every run of 8 bases or more that ended within the last 20 bases
   before it.
4. Of those, **the longer run wins, then the nearer one.** A base takes one entry or none — the
   measurement has no case of two effects adding, and this keeps the charge bounded. Longest first
   because the effect grows with length: twenty bases after a 16-base run the measured excess is 6
   in 100, at the base after an 8-base run 1 in 1,000.

This is the prototype's rule, except for the `N` guard. On HG002 chromosome 20 at 30× the build
reproduces the prototype's records exactly, except one beside a gap of `N`s that the prototype
counted as a run.

**Calling applies the rule to the candidate allele with the locus's reference flanks beside it**, as
the prototype does (`str_context_error.md` §3 point 1). **The pre-pass applies it to the reference**,
because what it counts are reads at reference positions.

## 3. Counting in the pre-pass

### 3.1 Where

**At every generic locus the pre-pass walks, at every position of it whose base has an entry in
either direction.** The census writer already visits every generic locus for the per-read-group
base-quality totals (`CensusWriter::add_locus`, `parameter_estimation/joint/census.rs`), not only the
two million positions it keeps, and this count goes beside those totals.

**Why generic positions and not the STR loci.** The bases after a run are mostly generic ground:
default typing makes a run of one base of 6 or more an STR segment, so the base just past its end
belongs to the generic path. An STR locus's observation holds the tract only, so the base after a
whole run is not in it. Fitting on generic positions and applying at STR loci also keeps the loci a
table is applied to apart from the ones it was fitted on.

**What a generic locus gives the count.** Each observation is a group of reads from one read group
showing the same bases, with how many of them were forward reads (`num_fwd`) and their summed log
error (`q_sum`). So at each position the count knows, per read group, how many forward and reverse
reads showed each base. It does not know where each read started — §3.4 says what that costs.

Only observations whose bases are as long as the locus are used, compared position by position with
the reference: a read carrying an insertion or deletion at the locus says nothing about a substitution
there.

### 3.2 Which reads are through and which arrived first

At a position `x`, take the entry the rule of §2.3 gives a **forward** read there, computed on the
reference. If it has one, and a **reverse** read at `x` has none:

- the position's forward reads are counted as **through** reads for that entry;
- its reverse reads are counted as **arriving first** for the same entry.

And the mirror image for an entry a reverse read has and a forward read does not. **A position where
both directions have an entry is not counted**: neither group is a clean comparison for the other
(`A`×16 `G` `A`×10 is one: the forward reads crossed sixteen `A`s, the reverse reads ten).

### 3.3 What is kept

Per read group and per entry, seven totals:

| total | meaning |
|---|---|
| `through_reads` | reads through the run at counted positions |
| `through_as_the_run` | of those, reads showing one more base of the run (next-base entries only) |
| `through_other_wrong` | of those, reads showing any other non-reference base |
| `arriving_first_reads` | reads arriving first at the same positions |
| `arriving_first_as_the_run` | of those, reads showing one more base of the run (next-base entries only) |
| `arriving_first_other_wrong` | of those, reads showing any other non-reference base |
| `through_as_the_run_log_error` | summed `q_sum` of the reads in `through_as_the_run` — the base qualities' own claim, for the check of §4.2 |

For a further-on entry the "as the run" totals stay zero and "other wrong" means any wrong base.

### 3.4 What this count does not see: reads that started inside the run

**Some forward reads at a position after a run started inside the run, or between it and the
position, and so did not cross all of it.** The count cannot tell them apart, because an observation
does not keep where its reads start. The scripts behind `str_context_error.md` §2 could, and left
them out.

They dilute the through group. With reads of 150 bases starting uniformly, the share of forward reads
at a position that started within the last `L + d` bases — a run of `L` bases ending `d` bases
before — is about `(L + d) / 150`: **about 1 read in 9 at the base after a 16-base run, about 1 in 4
twenty bases further on.** The fitted excess is lower than the true one by that share.

**What that costs at calling is small.** At the base after a 16-base run the excess becomes about
0.08 instead of 0.09, and the charge a misread pays changes by `ln(0.09 / 0.08)`, about 0.1 nats,
against the 3 to 6 nats a read the term exists to move (`str_context_error.md` §4). A correction would
need each read's start in the census; I recommend not making it. *(Measured 2026-11-02: above 15
bases the fitted rates are 30 to 40 in 100 below the prototype's, more than this predicts, and the
calls show no cost — the report's §3.)*

### 3.5 True variants at counted positions

**Not removed, because the comparison cancels them.** A heterozygous variant at a counted position
raises both groups' rates by about a half; a homozygous one by about one. Summed over the positions
of an entry, they move both rates by the same amount, and the difference is unchanged. They do add
noise: an entry fitted from few positions can be pushed either way by one variant, which is what the
floor of §4.1 is for.

Excluding positions that look variable would be worse: any rule that drops a position for showing
the run's base drops it more often when the misread is strong, and so biases the excess down.

### 3.6 Stored in the census

The totals are a new per-read-group table in the census's header, written after the base-quality
totals, as version 4 added those. **The census format goes from version 4 to 5**, so every existing
psp file must be regenerated before `estimate-parameters` reads it; the freshness check already
refuses an older census and says so. The merge across samples needs nothing new: a read group
belongs to one sample, so its table moves with it.

## 4. The fit

### 4.1 From totals to the table

For each read group and entry:

```
excess as the run  = max(0, through_as_the_run / through_reads
                            − arriving_first_as_the_run / arriving_first_reads)
excess other wrong = the same with the "other wrong" totals
```

with a binomial standard error for each, written to the parameters file beside it. **An excess less
than two standard errors above zero is charged as zero** (added 2026-11-02, after the first fit on
HG002): such an entry cannot be told from none, and charging it would cost every read near such a
run a second scoring (§6) for nothing the data showed.

**An entry with fewer than 1,000 reads in either group takes the default table's value** (§5), marked
`Defaulted` like any other parameter that fell back. At 1,000 reads an excess of 1 in 100 has a
standard error of about 3 in 1,000, so the floor is where an entry can tell a real excess from none.
Every entry is fitted independently: nothing forces longer runs to charge more, because at the floor
the data have already earned the right to say otherwise, and below it the default does.

### 4.2 The check against base qualities

**`str_context_error.md` §3 says the qualities check the table**, and the census carries what they
claim. For each fitted next-base entry, the fit compares the rate the base qualities imply for the
misread bases — `exp(through_as_the_run_log_error / through_as_the_run)` — with the rate the reads
through the run show for reading the base as the run, `through_as_the_run / through_reads`. It
writes both to the fit's log. *(Corrected 2026-11-02: the first build compared with the substitution
rate plus the excess, which leaves out the background misreads both directions share; at short runs,
where the excess is near zero, every entry then warned.)* Where an entry has 100 or more misreads and
the two differ by more than a factor of 3, the log carries a line naming the read group and the
entry. **The check never changes a value**: it exists to make a broken count visible. *(Changed
2026-11-03 from a warning to an information line: on the tomato2 cohort it fired for 43 of 59 read
groups, after runs of 8 to 12 bases in both directions — qualities better than the misread rate
after runs of 8 and 9, worse after 10 and 11 — on a dataset with nothing broken.)* On HG002 the
misread bases carry quality about 10, a rate of 1 in 10, against a measured rate of about 1 in 10
after a 16-base run.

## 5. The parameters file and the defaults

**The table travels as a new table of the parameters file**, `[misread_after_runs]`, one row per read
group and entry, holding both excesses, their standard errors, the read counts and the provenance.
The file's format goes from version 7 to 8; a version-7 file is read with the table absent, and every
entry then takes the default, as version 7 did for the table it added.

**The defaults are the table this pre-pass fits on HG002 chromosome 20 at 30×**, written into the
code once and regenerated by the same counting the fit uses. They replace the prototype's hand
numbers, which came from separate scripts. So `--defaults` and the fitted route charge from the same
definition.

**What the defaults assume.** They are one Illumina library's numbers. A read group from another
platform that falls back to them for a sparse entry is charged HG002's misread, which is a guess; the
provenance says which entries are guesses.

## 6. Calling

**The STR scoring reads each read group's table instead of the constants.** Concretely:

- `RunParameters` and its frozen view hold the table per read group, filled from the file or from the
  defaults;
- the per-locus gathering that already looks up each `(read group, candidate)` substitution rate
  hands the scoring context its read group's table;
- `run_context.rs` takes the table as an argument where it held constants, and its scoring applies
  the rule of §2.3;
- a candidate is scored once per read direction where a run of 8 bases or more, or 7 — one slipped
  base short of charged — sits in the candidate or beside it in its flanks, and the read group's
  table charges anything. A table of zeros costs no second scoring.

Unchanged from the prototype: reads that saw only part of a tract are scored at the flat rate, since
they carry no forward-read count; and the slippage re-fit's read attribution reads the forward
direction's score for every read (a small approximation the prototype already carries).

## 7. One sample and thousands, three reads and three hundred

- **One sample at 3×.** Entries are pooled over the whole genome per read group. The scripts behind
  `str_context_error.md` §2 saw 2,935 reads through 16-base runs in 10 Mb at 30×; a whole genome at
  3× gives about thirty times that, so the common entries clear the floor of §4.1 by far. The rare
  ones — runs of 22 or more — may not, and take the default.
- **A small region or a small simulated genome.** Most entries fall below the floor and take the
  default. A simulation that plants no misread fits an excess near zero where it does clear it,
  which is right for that data.
- **Thousands of samples.** The count is seven totals per entry per read group, about 1.1 kB a read
  group; the scoring cost is the prototype's, at most two scorings per distinct read sequence.
- **Three hundred reads.** Nothing changes in kind. The count's cost grows with reads, as every
  census count does.
- **The cost already measured.** At 5×, a heterozygous variant inside a long run with a single
  reference read can be called homozygous when that read crossed the run
  (`str_context_error.md` §6.3). Fitting the table does not remove this; it can only make the
  charge the right size.

## 8. Tests — each fails without the piece it covers

1. **The rule (§2.3).** One function gives a base's entry from a sequence and a direction; tests for
   a run of one base, the precedence order, a run of `N` and a two-base repeat taken as no run, and
   the prototype's cases unchanged.
2. **The count (§3).** A generic locus after `A`×12 with forward reads showing `A` at the next base
   and reverse reads showing the reference: the forward reads land in `through_as_the_run` of the
   12-base entry, the reverse reads in `arriving_first_reads`; a position with entries in both
   directions is not counted; an observation of the wrong length is skipped.
3. **The census file.** The table round-trips, and the version is 5.
4. **The fit (§4.1).** Totals give the excess; a negative difference gives zero; an excess within two
   standard errors of zero gives zero; a sparse entry takes the default and is marked `Defaulted`.
5. **The check (§4.2).** A table whose misreads claim quality 30 at a fitted rate of 1 in 10 logs a
   warning.
6. **The parameters file (§5).** The table round-trips; a version-7 file reads with every entry
   defaulted.
7. **Calling (§6).** A read group whose fitted table is all zero scores every read bit for bit as
   main does; the chr20:44,076,996 case of `str_context_error.md` §1 still flips with the default
   table; each read group scores under its own table.

## 9. Measuring it

Done on 2026-11-02 for the first build, which also covered repeats of 2 to 4 bases; the numbers are
in `../../reports/ng_str_context_error_fit_2026-11-02.md`.

1. **The table against the prototype's.** Fitted on HG002 chromosome 20 at 30×, entry by entry
   against `str_context_error.md` §2 — the size of §3.4's dilution.
2. **Calls.** HG002 chromosome 20 and the tandem-repeat benchmark at 30× and 5×, `--defaults` and
   the fitted route, against main and against the prototype (`str_context_error.md` §6.3).
3. **The simulations.** The flank and stuttering simulations, where nothing should change beyond
   noise.
4. **A cohort at 3×.** The 50-sample stuttering simulation is the one cohort already set up; its
   reads have no misread planted, so it checks that the fitted table stays near zero there.
