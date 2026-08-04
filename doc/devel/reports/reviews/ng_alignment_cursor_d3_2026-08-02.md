# ng — the alignment cursor, D3: review

*Two agents in separate worktrees over the same staged change (base `5ec6a60`): one on
correctness and design fit, one on whether the tests can fail.
Implementation: [ng_alignment_cursor_d3_2026-08-02.md](../implementations/ng_alignment_cursor_d3_2026-08-02.md).*

**Both agents found the same defect, independently, and both proved it with a test.** It is the
one this step stopped on, and the owner chose the fix.

---

## Blocker — a cohort would have reported one sample N times

The generator keeps its reader for a whole chromosome but is handed a `&SampleReads` on every
call, and the reader was keyed on the chromosome alone. So the sample argument was read once,
when the reader was opened, and ignored after that.

`ng_ssr_cohort_stutter` asks every sample about one repeat before moving on, so the first sample
would have answered for all of them: one plant's reads, N times, under N names, with no error
and a table of the right shape.

Proven both ways round — sample A has a read on the repeat, sample B's nearest is 45 bases away
and cannot reach it, and B reported A's read.

**Fixed** on the owner's decision: one generator per sample, and the generator refuses a sample
it was not opened for (`LocusGenerationError::ForeignSample`). `SampleReads::identity` is the new
interface that makes the comparison possible; it compares by files, not by sample name. The
cohort tool now builds a generator per sample. The same guard went onto the generic generator,
where no caller loops samples yet but the same trait makes the same mistake available.

Two tests each side: one that a foreign sample is refused, one that a generator per sample gives
each sample its own reads — because "refuses everything" would pass the first on its own.

---

## Major, fixed

- **The retirement tally could lose four of its five numbers, or overwrite instead of adding, with
  the whole suite green.** Both mutations survived, because only one chromosome was ever retired
  in a test and one retirement makes adding and overwriting identical. This matters beyond
  bookkeeping: the reuse-versus-jump ratio is the *only* evidence that the reader is being kept
  at all, so a fold that drops it would leave the feature detectable on the first chromosome and
  undetectable on every other. Fixed twice over — `CursorCounts` gained `AddAssign` so the five
  numbers cannot be folded partially at any of the four sites, and a chromosome-1 → 2 → 1 test
  distinguishes adding from overwriting.
- **A public method carried the wrong doc comment.** Two doc blocks had merged, so the summary
  rustdoc shows for `cursor_counts` described a different method, and the method it belonged to
  had none.
- **"Still strictly better than a query per locus" was asserted with no measurement**, on the one
  generator where the saving might have been negative — repeats are far apart, and the previous
  checkpoint's audit showed the ratio collapsing with distance. Replaced with the measured
  number: −28 % of CPU time, output byte-identical.

## Minor, fixed

The `+ Send` bound the previous review recommended was reverted on the generic generator, where
the agent showed it costs nothing and nothing needed it removed — only the STR generator has the
caller that cannot satisfy it. Two assertions in a chromosome test re-stated the test's own input
and could not fail; replaced with the one that carries weight. Stale statements about the STR
generator still using the old read path, in the read module and a test fixture's doc.

## Known and accepted

The reviewers listed things that are correct but unwitnessed, or reachable only through paths
this step does not add: the depth cap has never been exercised on a reader that is reusing (I
checked; it holds), the cursor error paths through the STR generator have no test, and
`delimit_segment_reads` — a bake-off diagnostic — now moves the shared reader, so calling it
"read-only" is true of the tallies but no longer of the reader.

## Deferred to the owner

The design documents for both generators now describe the old arrangement. The previous step's
review already raised `locus_generation_pileup.md`; `locus_generation_ssr.md` joins it, and so
does one line of `alignment_cursor.md` that names the STR generator as a caller of the old read
path. Raised at Checkpoint D — this skill does not edit design documents, and the plan does not
include a fold-in.
