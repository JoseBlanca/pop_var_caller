# Code-quality map: `src/locus_generation/` (2026-11-26)

Read-only review. Nothing was built or run. Every `path:line` below was checked against the
worktree at `f84129cf`. The production part of `pileup/` was reviewed by a helper agent; I
re-read its main citations (B1, A1, A2) and corrected one claim of its (see D4).

## Summary

1. **One read can stop a whole run.** A read with a deletion longer than 5,000 bp fails the
   generic walk with `RecordTooWide`. That error ends the sample's locus stream, so the run
   stops. No read filter removes such reads, although the error message says one should. Long
   reads produce such deletions; short reads rarely do. This is the only finding that can
   change calls (B1).
2. **`ssr.rs` already has four seams, and they should become files.** Its 4,494 production
   lines are a header of types plus three inline modules: read coordinates (350 lines), the
   per-read decision (2,350 lines, half of them tests) and the tally (740 lines). Moving each
   to `ssr/` and moving `ssr_bundle.rs` under it as well is a pure move. It also lets six
   `pub(super)` items and eight `pub` items become private (A1, A2).
3. **The two generators each keep the same per-chromosome cursor bookkeeping.** About 60
   lines each: the foreign-sample check, the counters of retired cursors and the per-read-group
   drop totals. A small shared struct removes one copy (D1).
4. **The STR path has eleven hand-written CIGAR walks.** Each one re-states which operations
   move along the reference. Two of them compute the same predicate (D2).
5. **Two counters cost speed and are never read.** `column_census` adds one write to a
   process-wide atomic at nearly every position; `minted_error_census` lost its only reader in
   `0f9281c5`. Deleting both removes about 220 lines (A3).
6. **A log warning says "the run will fail" when it will not.** The chain-id allocator gets a
   real cap in every run, contrary to the comment beside it, and two error messages name a
   `--max-active-reads` flag that does not exist (B4).
7. **The history in the comments outweighs the code.** In `generator.rs`, `open_record.rs`
   and `genome_walk.rs`, comment lines outnumber code lines. Much of it narrates a deleted
   "production" caller and plan steps (D3, C2, "Milestone C review") (R1).

`locus_generation` depends on `crate::parameter_estimation` only through three doc-comment
links (`pileup/minted_error_census.rs:5`, `:19`; `pileup/open_record.rs:3594`). It is not a
code dependency, so there is nothing to invert. Two of the three go away with A3.

---

## Bugs

### B1. A deletion longer than `max_record_span` stops the run
- **Where:** `pileup/open_record.rs:1788-1795` (`open_new`) and `:1706-1713` (`widen`); the
  message is at `pileup/errors.rs:71-75`.
- **What:** A read with CIGAR `100M6000D100M` in a generic region gives a record span of
  6,001, above the default cap of 5,000, so `open_new` returns `WalkerError::RecordTooWide`.
  A 4,900 bp deletion anchored 200 bp into an open record fails the same way in `widen`.
  The error becomes `LocusGenerationError::Walker`. That variant is documented as fatal
  (`mod.rs:862-866`), and `SampleLocusObservationsIterator` fuses on it (`mod.rs:1313-1316`).
  The helper agent found no deletion-length filter in `read/filtering.rs` or
  `bam::alignment_input::cigar_is_bad`. The repeat slide is capped so that it cannot trigger
  this (`open_record.rs:3173-3192`); the deletion's own width is not.
- **Proposal:** Drop such reads before the walk, counted under a new `ReadFilterCounts`
  reason. The alternative is to drop only the event in the walk and count it.
- **Downside:** Reads carrying a real large deletion are lost from the generic path. Today
  they are lost too, together with the run. A deletion-length filter is a new read-filter
  policy, so the owner should decide it.
- **Effect:** output changes only on inputs that fail today. **Cost:** M.
  **Confidence:** high that the error fires (likely); medium that it occurs at all on
  short-read data.

### B2. STR loci are far shorter than 65,535 bases only because of a default
- **Where:** `witness.rs:324-330` (`LocusLen::from_positions` saturates at `u16::MAX`);
  `ssr.rs:1931-1932` (a partial read's reach is clamped the same way).
- **What:** If a user raises `--max-str-len` or `--max-str-bundle-len` above 65,535, a tract
  that long gets `LocusLen` 65,535. Its partial witnesses then cover fewer positions than the
  tract has, and depth along it (`num_obs_along_locus`) is short at the far end. I found no
  check on either option's upper bound. With the defaults (100 bp each) this cannot happen.
- **Proposal:** Reject `--max-str-len` and `--max-str-bundle-len` above `u16::MAX` where the
  options are parsed.
- **Downside:** It rules out a configuration nobody runs today.
- **Effect:** none on supported inputs. **Cost:** S. **Confidence:** medium (suspected; I
  did not trace the option parsing end to end).

### B3. Run counters that miscount (calls unaffected)
- **Where:** `pileup/genome_walk.rs:1186-1197` returns before the update at `:1311-1313`;
  `genome_walk.rs:1309-1310` and `:1370-1372` compare a per-read-group cap
  (`sample_to_cap`, `:2092-2108`) with a whole-column depth.
- **What:**
  - `column_depth_high_water` never sees positions handled by the fast column path, which
    are most positions at 30×.
  - `positions_short_of_cap` compares a column's total depth with a cap that applies per read
    group. With two or more read groups it counts the wrong positions.
  - Both counters are diagnostics only.
- **Proposal:** Pass the depth out in `FastColumn::Emitted` and update the counter there.
  Compute the shortfall per read group inside `sample_to_cap`.
- **Effect:** counters only. **Cost:** S. **Confidence:** high.

### B4. A warning that says the run will fail, and two messages naming a flag that does not exist
- **Where:** `pileup/genome_walk.rs:842-857` (the comment says the allocator gets a cap it can
  never reach); `pileup/generator.rs:959-962` (the generator passes
  `config.max_active_reads`); `pileup/chain_id_allocator.rs:458-474` (the warning);
  `pileup/errors.rs:46-49`, `:63-65` (the two messages).
- **What:**
  - In every real run the allocator is capped at `max_active_reads` (32,768). Its one-shot
    warning fires at three-quarters of that and tells the user the run will fail. The run
    will not fail, because `admit_read` evicts reads first.
  - The two messages tell the user to change `--max-active-reads`, which no command offers.
  - The generator's doc does name this warning (`generator.rs:276-279`), but describes it as
    "inherited verbatim from production" rather than saying it is wrong.
- **Proposal:** Pass `u32::MAX` in the generator too, as the comment intends, or fix the
  comment and the warning. Either way, rewrite the two messages to name what the user can
  actually change.
- **Effect:** removes one wrong log line. **Cost:** S. **Confidence:** high.

### B5. Slow paths that only bite at extreme depth (speed, not output)
- **Where:**
  - `pileup/active_read_set.rs:526-536` (`worst_sampling_key`, called from
    `genome_walk.rs:1077-1079`)
  - `pileup/genome_walk.rs:2076-2084`
  - `pileup/open_record.rs:684-688` and `:1099-1106`
- **What:**
  - Once the 32,768-read hold limit is reached, every arriving read rehashes the name of
    every read held: 32,768 hashes per read.
  - Above 1,000 reads in one read group at one position, the list of dropped reads is
    searched linearly on every insert. At 12,792 reads, the deepest position the helper agent
    found measured, that is about 70 million comparisons per position.
  - At 300× neither path runs.
- **Proposal:** Compute each read's sampling key once, in `ActiveReads::admit`, and store it
  on the read (8 bytes, 256 KB at the limit). Push dropped ids without the check, then sort
  and deduplicate once when the record is finalised.
- **Effect:** speed only; the keys and counts are identical. **Cost:** S–M.
  **Confidence:** high on the cost, medium on how often real data gets there.

### B6. Every STR locus re-reads about 90 reference bases it has just released
- **Where:** `ssr.rs:4218-4222` (releases the reference before `start − flank_bp`);
  `ssr.rs:4406-4417` (the slide fetch starts `MAX_SLIDE + MAX_INDEL_UNIT + 1 = 107` bases
  before the tract); `raw_chrom_reader.rs:393-396` (`prepend_backward`).
- **What:** With the default flank of 15 bp (`DEFAULT_BUNDLE_THRESHOLD`), `cursor_for` drops
  the reference buffer before `start − 15`. Then `slide_reach` asks for `start − 107`. So each
  STR locus makes a small backward file read of about 92 bases.
  - Correctness is unaffected, because eviction is only a hint.
  - The cost is one extra read call per STR locus per sample.
- **Proposal:** Release from `start − max(flank_bp, 107)`.
- **Downside:** about 92 more bases held, which is negligible.
- **Effect:** speed only, likely small; measure before claiming a gain. **Cost:** S.
  **Confidence:** medium (suspected; not timed).

### B7. At very high depth, reads that only touch a flank use up the STR read cap
- **Where:** `ssr.rs:454-494` (`fetch_capped_reads` samples every read overlapping the tract
  plus its margin); `ssr.rs:985-992` (the CIGAR-only check that later classifies flank-only
  reads as `OutsideTract`).
- **What:** The cap of 1,000 reads per locus is drawn from every read overlapping the tract
  plus 15 bp each side. A read that ends before the tract or starts after it is then
  discarded as `OutsideTract`. With 150 bp reads and a 30 bp tract, about 1 overlapping read
  in 7 is such a read, so the cap starts to bind at about 700× depth. When it binds, the
  locus keeps fewer informative reads than the cap allows. At 300× the cap does not bind.
- **Proposal:** Apply the CIGAR-only `OutsideTract` test before offering a read to the
  reservoir, and count those reads separately.
- **Downside:** **This changes output** wherever the cap binds, because the kept set
  changes. It also changes what `reads_discarded_by_cap` means.
- **Effect:** output change above about 700×. **Cost:** S. **Confidence:** medium.

---

## Architecture

### A1. Split `ssr.rs` along the seams it already has
- **Where:** `ssr.rs:502-847` (`mod read_region`), `:857-3203` (`mod classify`, production
  part `:857-1971`), `:3209-3940` (`mod tally`), `:313-494` (the read cap: `SplitMix64`, the
  seeds, `Reservoir`, `fetch_capped_reads`). The generator is `:3948-4492`; the tests are
  `:4494-5821`.
- **What:** The file is already four modules. Three are written inline and use `pub(super)`
  for their boundaries, and `ssr/slide_reach.rs` already sits in the `ssr/` directory. So the
  seams exist; only the files do not.
- **Proposal:**
  - Move `read_region`, `classify` and `tally` to `ssr/read_region.rs`, `ssr/classify.rs` and
    `ssr/tally.rs`, each with its own tests.
  - Move the read cap to `ssr/read_cap.rs`.
  - Leave `MeasuredRepeat`, `SsrLocus`, the config, the counts and the generator in `ssr.rs`.
  - Optionally, move the flank-deletion helpers (`classify` `:1226-1723`, about 500 lines)
    into `ssr/flank_deletion.rs`. They are one topic: how a deletion at a tract's edge is
    measured.
  - Do **not** split the generator from its types. It is 550 lines about one idea.
- **Downside:** `git blame` needs `-C` to follow moved lines, and a reader of `classify_read`
  now opens a second file for its helpers. A move can also hide an edit, so make each move
  its own commit with nothing else in it.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### A2. Put the bundle generator under `ssr/` and shrink the public surface
- **Where:**
  - Items made `pub(super)` only for `ssr_bundle.rs`: `ssr.rs:4328` (`MeasuredReads`),
    `:4336` (`split_flanks`), `:4362` (`claim_the_one_locus`), `:4374` (`contig_name`),
    `:4383` (`fetch_locus`), `:4432` (`measure_reads`).
  - `pub` items nothing outside `ssr.rs`/`ssr_bundle.rs` uses: `Reservoir` (`:385`),
    `CappedReads` (`:437`), `fetch_capped_reads` (`:454`), `seed_for_segment` (`:367`),
    `SsrLocus` and its public fields (`:151`), `MeasuredRepeat` (`:45`), `MeasuredGround`
    (`:64`), `RepeatDelimiter` (`:3948`).
- **What:** `pub(super)` in `ssr.rs` means visible to all of `locus_generation`, including
  `pileup/`, although only the bundle generator uses those items. The `pub` items are part of
  the library's public API (`pub mod ssr`) and have no outside user. I grepped `src/`,
  `examples/`, `benches/` and `tests/`.
- **Proposal:**
  - Move `ssr_bundle.rs` to `ssr/bundle.rs`. A child module sees its parent's private items,
    so the six `pub(super)` items become private.
  - Keep the one outside path, `locus_generation::ssr_bundle::SsrBundleGenerator`
    (`run/walker.rs`), working through a `pub use`, or update that call site.
  - Narrow the eight `pub` items to `pub(super)` or private.
- **Downside:** A bundle is its own region kind, and a reader looking for its generator may
  not look under `ssr/`. The `mod.rs` doc can say where it is.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### A3. Delete the two measurement counters, one of which costs speed
- **Where:**
  - `pileup/mod.rs:461-506` (`pub mod column_census`)
  - `pileup/genome_walk.rs:1378-1463` (86 lines inside `process_position`)
  - `pileup/fast_column.rs:550` (an unconditional `FAST_COLUMNS.fetch_add`)
  - `pileup/minted_error_census.rs` (production lines 1-194, made `pub mod` at
    `pileup/mod.rs:112`) and its call sites at `generator.rs:1173-1179`,
    `fast_column.rs:423-435` and `open_record.rs:931`, `:975-981`
- **What:**
  - No code outside `pileup/` reads `column_census`; I grepped `src/`, `examples/` and
    `tools/`.
  - `FAST_COLUMNS` is a process-wide atomic incremented at every fast-path position by every
    worker thread, so all workers write the same memory location at nearly every position.
  - `minted_error_census` lost its only reader in commit `0f9281c5`. If switched on, it
    would also be wrong under parallel samples: it is one global table keyed by locus start
    alone.
- **Proposal:** Delete both modules and their call sites; git history keeps them. Keep the
  `PVC_FAST_COLUMN=0` switch (`genome_walk.rs:869-870`), which tests that the fast path and
  the general path agree.
- **Downside:** Re-measuring the minted error means would need the code restored from
  history.
- **Effect:** output none; slightly faster. **Cost:** S. **Confidence:** high.

### A4. The read query's extra margin may be unnecessary now (check before acting)
- **Where:** `pileup/generator.rs:681-707` (`query`, and its rationale).
- **What:** Each record anchored in a region is clamped at the region's last base
  (`open_record.rs:3007-3016`, `:3164-3192`), and the walk stops about one base past it
  (`genome_walk.rs:326-343`). So every read that can fold into a kept record overlaps the
  region itself, and the extra margin seems to fetch reads that are never used. The rationale
  in the comment predates the clamp.
- **Proposal:** Query the region alone, and keep the margin only for releasing the reference
  behind the walk. Adopt it only if the output-identity check shows byte-identical output.
- **Downside:** If some record does reach past the region's end, its support would be lost
  without any error.
- **Effect:** none expected; less cursor memory. **Cost:** S plus the check.
  **Confidence:** medium.

---

## Duplication

### D1. The two generators keep the same per-chromosome cursor bookkeeping
- **Where:**
  - `ssr.rs`: fields `:4006-4026`, `read_filter_counts` `:4159-4167`, `cursor_counts`
    `:4136-4142`, the sample check and retirement `:4208-4248`
  - `pileup/generator.rs`: fields `:885-898`, `cursor_counts` `:1011-1017`,
    `read_filter_counts` `:1029-1036`, the sample check `:1239-1251`, retirement
    `:1379-1392`
- **What:** The same five jobs are written twice, about 60 lines each:
  1. remember the sample and refuse another (`ForeignSample`);
  2. when the chromosome changes, add the retiring cursor's `CursorCounts` to a running
     total;
  3. add its per-read-group drop counts to a `BTreeMap`;
  4. report the total plus the live cursor's counts;
  5. adopt the sample only after a successful open.

  The `pileup` comment at `generator.rs:1231-1233` says the check was added after the
  mistake was made on the STR side.
- **Proposal:** One small struct in `locus_generation` (say `CursorLedger`) holding the
  sample and the two retired totals. It needs three methods: `check_sample`,
  `retire(&SampleCursor)` and `totals(Option<&SampleCursor>)`. Each generator keeps its own
  cursor, because the pileup cursor lives inside its walker.
- **Downside:** It couples the two generators' cursor lifecycles a little. They already
  share `SampleCursor`, so that is not new coupling.
- **Effect:** none. **Cost:** S–M. **Confidence:** high.

### D2. Eleven hand-written CIGAR walks in the STR path
- **Where:** `ssr.rs:530` (`read_footprint`), `:554` (`ref_to_read`), `:649`
  (`tract_span_from_input`), `:1355` (`reference_bases_deleted_in_the_flanks`), `:1402`
  (`deletion_takes_flank_bases_at_a_tract_edge`), `:1458` (`through_the_flank_indels`, two
  loops at `:1474` and `:1515`), `:1555` (`indel_in_the_guarded_flank`), `:1598`
  (`flank_bases_deleted`), `:1683` (`aligned_bases_beyond_the_tract`), `:1710`
  (`unaligned_bases_at_the_ends`).
- **What:**
  - Each walk re-states which CIGAR operations move along the reference and which move along
    the read. They do not all agree on `Skip` (`N`): some treat it as a deletion, others as a
    match. That is harmless on DNA-seq, where `N` does not occur.
  - **Two of the walks compute the same predicate.** `deletion_takes_flank_bases_at_a_tract_edge`
    (`:1417-1418`) tests exactly the two conditions `flank_bases_deleted` tests
    (`:1613`, `:1616`). The only difference is that the second returns `None` above 255
    deleted bases.
- **Proposal:**
  - Merge the pair: make `flank_bases_deleted` return the raw `(u32, u32)`, and let both
    callers derive what they need from it. Cost S.
  - Optionally, add one iterator over a CIGAR yielding each operation with its reference
    range and read range. The nine walks become short folds over it. Cost M.
  - Do not merge these with `pileup/cigar_cursor.rs`. That code answers a different question
    (the events at one position), is tuned for the hot loop, and should evolve apart.
- **Downside:** The iterator refactor touches the code that decides how a flank deletion is
  measured, which took several measured rounds to get right. Do it only with the classify
  tests and the output-identity check green.
- **Effect:** none if done exactly. **Cost:** S (the merge), M (the iterator).
  **Confidence:** high on the merge, medium on whether the iterator earns its keep.

### D3. The right-junction insertion rule is written once on each path
- **Where:** `ssr.rs:696-701` (an insertion after the tract's last base is the tract's when its
  last base equals the tract's last base); `pileup/open_record.rs:2701-2710`
  (`claimed_junction_insertion`, the exact complement: the generic path claims it otherwise).
- **What:** One ownership rule is written in two places as two comparisons. Its doc
  (`open_record.rs:2672-2678`) says every insertion is owned by exactly one path only because
  the two agree. If one copy changes alone, insertions at a repeat's right edge are counted
  twice or by nobody.
- **Proposal:** One `pub(crate)` predicate in `locus_generation`, for example
  `insertion_at_right_junction_extends_repeat(last_inserted, repeat_last_base) -> bool`,
  called by both paths.
- **Downside:** It is a one-line comparison, and the shared function adds a hop for the
  reader. The value is that changing it on one side alone is no longer possible.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### D4. Small duplicates to leave alone
- **The base-quality conversion is written three times:** `ssr.rs:1868-1871`
  (`ln_p_err_sum`), `pileup/open_record.rs:3645-3656` (a table) and
  `read/prepared_read.rs:150-155`. All three use the same formula, `−q·ln10/10`. Merge them
  only if the `read` module's owner agrees; the gain is small.
- **`ssr.rs`'s `SplitMix64` is not `crate::splitmix`.** Its increment is `0x9e3779b97f4a7615`
  (`ssr.rs:331`). The standard constant, used in `src/splitmix.rs:9` and five other places,
  is `0x9e37_79b9_7f4a_7c15`; the two differ in one hex digit. The STR value is still odd,
  so the sequence is sound, but swapping in the shared function would change which reads the
  cap keeps. Leave the value. Fix the comment at `ssr.rs:318-320`, which calls it SplitMix64
  "ported verbatim", so nobody "fixes" it into an output change. (This corrects the helper
  agent's note that the STR copy is the same generator as the test fixture's.)
- **`pileup/read_sampling.rs:91-97`** is a third variant, the finaliser alone. It is
  different by design, and changing it changes the kept reads.
- **`SsrGenerator::next_locus` and `SsrBundleGenerator::next_locus`** each assemble the output
  (`ssr.rs:4306-4323`, `ssr_bundle.rs:166-179`). They differ in the region and the kind
  payload; a shared helper would take as many arguments as it saves lines.

### D5. Test-only methods that duplicate a production one
- **Where:** `mod.rs:207-238` (`non_reference_reads`) and `mod.rs:240-264`
  (`reads_compared_with_reference`), both `#[cfg(test)]`, beside the production
  `non_reference_and_compared_reads` (`mod.rs:276-289`).
- **What:** Two methods, 55 lines with their docs, exist only for tests. They repeat the
  production method's loop, which already returns both numbers.
- **Proposal:** Delete them and have the tests call the combined method. Keep the doc
  paragraphs that explain why partial reads are excluded, on the combined method.
- **Effect:** none. **Cost:** S. **Confidence:** high.

---

## Readability

### R1. Comments narrate a deleted caller and plan steps
- **Where:**
  - In `ssr.rs` (production part), "production" appears 27 times and "parity" 11 times. The
    module header `:1-16` says it is "a port of production `src/ssr/pileup/`" that "lands
    across the STR generator plan".
  - Section banners `:313-316`, `:431-433`, `:496-501`, `:849-854`, `:3204-3208` name plan
    steps (D1, D2a, D2b, D3).
  - `seed_for_segment` (`:363-366`) and `with_default_aligner` (`:4065-4067`) justify
    choices by "parity" with a caller that no longer exists.
  - In `pileup/`, `mod.rs:1-97` describes "eight files transcribed … unchanged".
    `generator.rs:252-266` says "Ten fields" on a 22-field struct. Eleven comments cite a
    nonexistent `ia/specs/…`, for example `pileup/decompose.rs:3`.
  - In `generator.rs`, `open_record.rs` and `genome_walk.rs`, comment lines outnumber code
    lines: 847 against 599, 1,977 against 1,621 and 1,033 against 1,004.
- **Proposal:** One pass per file. Delete provenance and plan-step narrative, and repair
  dead links. Keep every comment that explains *why* the current code does what it does,
  including the measurements. History stays in git and in `doc/devel/`.
- **Downside:** The owner values recorded reasoning, and length is not the target. Cut only
  text that describes code or reviews that no longer exist.
- **Effect:** none. **Cost:** M. **Confidence:** high.

### R2. `tract_end` means two different things in one file
- **Where:** `ssr.rs:108` (`MeasuredRepeat::end`, 1-based, inclusive) against `ssr.rs:986`,
  `:1067`, `:1405-1406` and the flank helpers (0-based, exclusive).
- **What:** Within one file, `tract_end` is sometimes the tract's last base and sometimes one
  past it. Every helper's doc has to say "coordinates as in …". The bundle helper already
  takes `tract: Range<u32>` (`:1260`), which carries the convention in its type.
- **Proposal:** Pass the 0-based tract as `Range<u32>` to the flank helpers, as
  `flank_deletion_the_delimiter_implies` already does, or rename the local
  `tract_end_exclusive`.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### R3. `classify_read` is a 280-line function with seven exits
- **Where:** `ssr.rs:944-1224`.
- **What:** Lines `:1015-1110` decide whether a read is spelled from its own alignment: the
  input span, the flank deletion, extending it through flank indels, the smallest form and
  the crossing test. The decision returns from the middle of the function, and its inputs
  are used again by the delimiter's branches below (`:1165`, `:1196-1217`).
- **Proposal:** Extract that block as `spelled_from_own_alignment(...) -> OwnSpelling`, a
  struct holding the input span, the flank deletion and whether the read crossed the tract.
  `classify_read` then reads as: outside the tract, no alignment needed, own spelling,
  delimiter.
- **Downside:** About eight values have to pass in or out. If the struct grows past that,
  leave the function as it is.
- **Effect:** none. **Cost:** M. **Confidence:** medium.

### R4. The tally's docs describe a three-field key; the key has five fields
- **Where:** `ssr.rs:3240` (`Support`), `ssr.rs:3254-3260` (`tally`), and the test helper doc
  at `ssr.rs:3605-3607`.
- **What:** They say observations are grouped by `(bases, read_witness, read_group)`. The real
  key (`ssr.rs:3223-3229`) also has `flank_deletion` and `held_flank`.
- **Proposal:** Refer to `ObservationIdentity` rather than listing the fields.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### R5. Two copies of the per-op event code in `cigar_cursor.rs`
- **Where:** `pileup/cigar_cursor.rs:319-417` and `:465-551`.
- **What:** The two are identical apart from `continue` against `return`. The linear loop can
  call the `#[inline(always)]` helper instead.
- **Effect:** none. **Cost:** S. **Confidence:** high.

---

## Looked at, left alone

- **`GeneratorSet::begin_region` decides `adjoins_repeat_ground` from the previous region
  only** (`mod.rs:1162-1166`). That is correct, because typed regions tile each contig with
  no gaps and in order.
- **A satellite's right junction:** the generic path does not claim an insertion whose last
  base equals the satellite's last base, and satellites are not called. By design, since
  satellites are out of scope.
- **The STR generator emits a locus even with zero reads.** Absent evidence is recorded as an
  empty locus. The generic path emits nothing without reads; the two shapes are intended.
- **`tally` sums `q_sum` in floating point in read order** (`ssr.rs:3310`). Read order is
  fixed by the reservoir, and the sum is rounded by `SummedLogError::from_nats`, so it is
  deterministic. On the STR path it is unread downstream
  (`calling/allele_candidates/ssr.rs:1001`).
- **`HashMap` in `tally`** is followed by a total sort (`ssr.rs:3363-3370`), so iteration
  order cannot reach output.
- **The STR reservoir's seed is the same for every sample at a locus** (contig name plus
  start, `ssr.rs:347-376`). It is deterministic and does not depend on thread count.
- **The `u32` casts of positions in the STR path** (`ssr.rs:375`, `:971`, `:3316`) truncate
  only beyond 4 Gb, which the BAM format cannot address.
- **At a contig start the left flank is empty** (`SsrLocus::fetch`, `ssr.rs:182`). Every
  index of `reference[left - 1]` is guarded by `left > 0` (`ssr.rs:1025`) or by
  `left < period` (`:1809`).
- **The two `while` loops of `in_smallest_form`** (`ssr.rs:1661-1676`) index the deleted
  flank bases correctly on both sides; I worked through the left and right cases.
- **The helpers in `slide_reach.rs` index the window safely**: `start - 1 - unit >= first`,
  because `start >= first + unit + 1` (`slide_reach.rs:85-98`).
- **The generic and STR read caps are different by design.** The generic cap is per read
  group per position and hashes read names; the STR cap is a per-locus reservoir. Do not
  merge them.
- **Do not merge `decompose()` (test-only) with `cigar_cursor`.** It is an independent
  reference the cursor is tested against.
- **`SampleLocusObservationsIterator`'s explicit `Drop`** (`mod.rs:1340-1369`) duplicates the
  field order on purpose; the comment says why.
- **`LocusGenerationError` variants with no `#[from]`** (`mod.rs:791-795`) are deliberate, so
  every error names its region.
- **u32 positions in the walk are checked against overflow** (`genome_walk.rs:1556-1564`).
- **`AHashMap` in the chain-id allocator is never iterated**, so its order cannot reach
  output.
- **Test setup:** `classify`'s tests build their aligner in three near-identical helpers
  (`ssr.rs:2016`, `:2132`, `:2182`). It is small, and each helper names its fixture, so
  leave it.
