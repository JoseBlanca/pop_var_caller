# ng — the alignment cursor, Milestone A: review fixes applied

*Companion to [the review](../reviews/ng_alignment_cursor_a_2026-08-02.md) and
[the implementation report](ng_alignment_cursor_a_2026-08-02.md). A section per step.*

---

## A1

**All findings applied**, with **two disputed on measurement** and one recorded as a stated
gap rather than a test. Nothing was deferred.

The suite went from 9 tests to 19, and four of the new ones failed when first written —
each failure a real property of the code or the fixture that had been assumed rather than
checked.

### The Blocker

**B1 — the typed-region walk could drop three quarters of its regions silently.**

Added `the_typed_walk_visits_every_region_the_stream_produced`, which counts the same typed
stream a second time without the walk in the way and compares three things: the region
count, the base pairs those regions cover, and the loci (one per base, since the fixture's
reads are pure Match and tile the contig).

**Verified by re-running the review's own mutation.** With
`walks.into_iter().flatten().take(1)`:

```
the walk saw 1 Generic regions; the stream produced 4
test result: FAILED. 18 passed; 1 failed
```

### The Majors

| finding | what was done |
|---|---|
| **M1** the discarded `stream.counts()` | `regions_in`, `regions_handled`, `loci_emitted` and **both unhandled counters** now live in `ProbeReport` and print. `the_walk_accounts_for_every_region_and_every_locus` asserts the spec §13.2 partition across three run shapes. |
| **M2** 3 of 16 printed keys pinned | `the_report_prints_exactly_the_documented_key_set` pins all twenty-one, in order. `each_printed_key_carries_its_own_counter` gives every field a distinct value so a transposition inside `render` cannot pass. The copy out of `PileupGeneratorCounts` is now one destructure, so a counter added upstream stops the build here. |
| **M3** untested contig filter | Hoisted to `contig_is_selected`, one body instead of two, with a unit test and an end-to-end test through **both** region sources. |
| **M4** the record-span knob | **Disputed on measurement — see below.** |
| **M5** untested fixture, false tiling claim | Three `assert!`s in `SyntheticSample::build`, a `stride()` method the tiling rule is stated in, and five tests including three `should_panic`. |
| **M6** the wrong `Cargo.toml` claim | Rewritten to what was measured. |
| **M7** the allocator guard's justification | Rewritten to say what it catches (a contradiction the operator typed) **and what it does not** (the silent case, because nothing in the output names the allocator). |

### Two findings disputed on measurement

**M4 — the record-span knob is unobservable on this fixture, and the proposed test asserted
an effect that does not occur.** The review suggested asserting that a 200 bp halo and the
default 5,000 bp one admit different reads. Measured at 500 bp regions, the two produce
**identical reports** — every counter, including `reads_admitted=425` and
`records_outside_region=1288`. That is not the knob failing to reach the generator: the walk
stops at the region's end and the read stream is lazy, so halo records are never pulled
unless something reaches back into the region, and nothing here does — these reads are 150 bp
of pure Match with no long deletion. Exercising the halo needs a fixture containing the thing
the halo exists for. **Recorded as a comment stating the gap**, because an assertion that
happens to hold would claim coverage this file does not have.

**M2's walk-side transposition cannot be caught here, and the code now says so.** Swapping
two assignments out of `PileupGeneratorCounts` stays green because seven of those eight
counters are zero on this fixture. The destructure stops a counter being *forgotten*; it
cannot stop two being *crossed*. Written down rather than left as implied coverage.

### Four assertions that failed when first written

Each one taught something the review and the implementation had both assumed:

1. **`regions_handled == regions_in` is false by design.** On the typed path, `4 of 7`
   regions reach a generator — the rest are repeat regions with no generator filled. The
   correct invariant is the three-way partition, which is why both unhandled counters were
   added rather than just the two the review named.
2. **An observation is a folded row, not a read.** At 10× the fixture yields **5,332**
   observations over 5,000 loci, not ~50,000: ten reads agreeing with the reference fold into
   one row. The assertion is now a band — surplus between 1 and the number of reads written —
   with the fold explained. The review's suggested `loci * 5` would have been wrong in the
   same way the original `> loci` was weak.
3. **The fixture does not tile at coverage 1.** 33 reads at a stride of 152 over 150 bp reads
   covers **4,950 of 5,000** positions. `build` now refuses any geometry whose stride exceeds
   its read length.
4. **The exhaustive `ProbeReport` literal earned its keep immediately** — adding the two
   unhandled counters failed to compile at exactly the test that pins the printed surface.

### One finding the orchestrator's own mutation testing added

After applying the fixes, the four mutations were re-run. Two still passed. One was M2's
transposition (above); the other was **loosening `contig_is_selected` to a prefix match** —
which my test cases missed, because none of them had a contig name *extending* the filter.
`chr21_KI270872v1_alt` is exactly that, and exactly what the anchor's reference is full of.
Added, and the mutation now fails.

### Validation

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test` (plain) — reaches the probe's **19 passed**.
- `cargo test --lib --tests --examples --all-features` — `2770 passed; 0 failed; 5 ignored`.
- `cargo bench --bench ng_generic_pileup_perf -- --test` — six cases, all `Success`.
- **The chr21 anchor, after every fix:** `loci=236081 observations=251786
  reads_admitted=54709` — unchanged. The new lines are additive:
  `regions_in=205875 regions_handled=102938 unhandled_not_implemented=102659
  unhandled_out_of_scope=278 loci_emitted=236081`, which partitions exactly and records for
  the first time that **half of chromosome 21's typed regions have no generator filled**.

### Not applied

Nothing. Two items are recorded in the review as out of scope for A1:

- `benches/psp_writer_perf.rs:386` panics under `cargo test --all-targets` (pre-existing,
  untouched by this branch, and the cause of a wrong measurement in this review).
- `examples/ng_generic_loci_dump.rs` holds a third copy of `write_fasta` / `write_bam`;
  converting it belongs to Milestone F1.
