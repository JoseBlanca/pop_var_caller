# Code Review: str_bundles_typing
**Date:** 2026-10-11
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step 2 of the STR bundle caller (typing only), uncommitted diff on `str-bundles` at 90e6f69e
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted diff (49 files) on branch `str-bundles`, HEAD `90e6f69e`.
  It gives each STR bundle member a cut span and a motif (`BundleMember`), on the scan route
  (`segment_criteria::bundle_member_from_bases`) and the catalog route
  (`repeat_catalog/segments.rs`, `bundle_member_from_row`); adds the bundle cap
  `--max-str-bundle-len` (default 100) and plumbs it through three subcommands, the psp header and
  the parameters file; drops a bundle on a contig's first or last base; and re-records three
  parameters-file checksums.
- **In-scope files:** `src/region_typing/{mod.rs,segment_criteria.rs}`,
  `src/repeat_catalog/{segments.rs,criteria.rs,builder.rs,parquet_file.rs}`,
  `src/psp/segmentation_section.rs`, `src/calling/parameters_file/*`, `src/cli/*` (the flag and
  dump), `src/cli/cross_platform_digests.rs`, `tests/ng_repeat_catalog_differential.rs`,
  `examples/ng_typed_region_dump.rs`, `scripts/ng_str_stutter_check_bundles.py`, and the five spec
  documents. Direct callees read: `repeat_catalog/reader.rs` (`rows_for`),
  `repeat_catalog/mod.rs` (`reach_into_window`), `repeat_catalog/row.rs`, `classify`,
  `split_bundles`, `locus_generation/mod.rs` (`begin_region`), the census digests.
- **Out of scope:** the simulator (reviewed in `str_bundles_simulator_2026-10-11.md`); the
  later steps (reads, psp bundle records, calling).
- **Categories dispatched.** By the owner's rule, no worktrees: two read-only reviewers in this
  worktree, neither building nor editing. *Reliability + refactor safety* (tests, mutants by
  reading), and *errors + defaults + naming + idiomatic + smells + module structure + doc quality*.
  The orchestrator ran every build, the mutation and the probes itself. `float_portability` and
  `unsafe_concurrency` do not apply (no floats reach output here, no concurrency). Their files are
  in `tmp/review_2026-10-11_str-bundles-typing/`.

## 2. Verdict

Approve-with-changes. The typing does what §5 asks on both routes and leaves every HG002 segment
and generic region untouched. One defect is new in kind: a windowed catalog read could type a chain
of short repeats as a bundle where the whole contig types it a satellite (M1). One design question
the typing surfaces must be settled before step 3: members overlap and nest (M2).

## 3. Execution status

All in the dev container, release build:

- `cargo test --release --no-fail-fast --lib` before any fix: 4,995 passed, 3 failed — the three
  debug-only `should_panic` tests known to fail in release.
- `cargo test --release --test ng_repeat_catalog_differential` before any fix: 12 passed.
- A probe (a temporary integration test, deleted after) for M1: see the finding.
- **The implementer's HG002 claims, re-derived from `tmp/bundles/typed.tsv` (main) and
  `tmp/bundles_step2/typed.tsv` (this change): CHECKED-CORRECT.** 2,473 → 2,410 bundles, 265 → 328
  satellites; the 20,183 `ssr_locus` and 70,608 `generic` rows identical in their first eight
  columns; the bundle and satellite spans the same set. 385 bundles end past every member's cut end
  and 7 start before every member's cut start, as claimed — but the docs say *last member's* cut
  end, and measured against the last-listed member it is 520 (M5).
- Not run: `cargo doc`, `cargo audit` (no dependency changed).

Needs verification: 0.

## 4. Open questions and assumptions

1. **(M2)** Members overlap and nest. Should a member wholly inside another be dropped at typing,
   or should step 3's parse take overlapping members? Recommendation below.
2. **(Mi8)** A bundle dropped at a contig end is not charged as a `FlankClamped` rejection in the
   tally, although §5.3 says it is dropped "as an STR segment is". Charge it, or say in §5.4 that
   the tally does not?

## 5. Top 3 priorities

1. **M1** — a windowed read types a long chain as a bundle; the whole contig types it a satellite.
2. **M2** — overlapping and nested members, before step 3 builds the parse on §4's picture.
3. **M5** — "first repeat to last" in the help, the docs and the note written into every
   parameters file, which is wrong for 520 of 2,410 HG002 bundles.

## 6. Findings

### Major

**M1: [src/repeat_catalog/mod.rs](../../../../src/repeat_catalog/mod.rs#L466) — a windowed catalog read types a chain of short repeats as a bundle where the whole contig types a satellite**
- **Confidence:** High (probed). **Categories:** reliability, orchestrator.
- **Problem:** `rows_for` pads each requested span by `reach_into_window` = the contig's longest
  tract plus the bundle radius. A bundle is a chain of repeats each within the radius of the next,
  so its hull is not bounded by any one tract. Probe: a contig with six repeats of 20 to 24 bases,
  10 apart, 180 bases first to last, longest tract 24. Whole contig, scan and file alike:
  `satellite [200, 380]`. A read over bases 205–215: scan `satellite [200, 380]`, **file
  `bundle [200, 254]`**. Before the cap the same read gave a 55-base bundle against a 181-base
  bundle — wrong extent, which predates this change; the cap makes it the wrong kind too, and a
  bundle is what step 4 will call.
- **Why it matters:** production types through `genome_segments(ReadScope::Regions)`, so typing
  depended on the BED. At the default cap it needs a contig whose longest tract is under about 85
  bases; HG002's typing does not move (checked after the fix), but a raised cap (§5.2, "for long
  reads") or a small contig reaches it.
- **Suggested fix:** pad by `max(longest tract, bundle cap) + radius`. A bundle at most the cap
  then has every member inside the padded window, and a longer chain is still seen longer than the
  cap. The satellite's extent can still be cut at the padded window, as before.

**M2: [src/region_typing/mod.rs](../../../../src/region_typing/mod.rs#L285) — members overlap and nest; §4, §6.1 and §7.1 assume ordered, disjoint runs**
- **Confidence:** High (measured on `tmp/bundles_step2/typed.tsv`). **Categories:** orchestrator,
  both reviewers.
- **Problem:** of HG002's 2,410 bundles, **551 hold one member wholly inside another** (e.g.
  `chr1:8007443`, a `(CT)n` member `[8007442, 8007506)` and a `T` member `[8007497, 8007506)`
  inside it), 734 have consecutive members that overlap, 6 list members whose cut starts are out of
  order, and in 520 the last-listed member does not reach the hull's end. §4 defines a gap of 0 to
  14 bases; §6.1 takes "the last member's motif at the right edge"; §7.1 parses an allele as each
  member's copies in order with gaps between. None has a reading for a member inside another.
- **Why it matters:** the typing is what hands step 3 its members. Building the parse on §4's
  picture would mis-parse about one bundle in four.
- **Suggested fix:** an owner decision (question 1). Recommendation: drop a member wholly
  contained in another at typing — the container already spells its bases — after measuring what
  it does to HG002's member counts; leave partial overlaps to the parse. Until then, document it
  (done in this review's fixes).

**M3: [src/region_typing/mod.rs](../../../../src/region_typing/mod.rs#L926-L933) — nothing pins the hull to the detected spans any more**
- **Confidence:** High. **Categories:** reliability.
- **Problem:** the diff replaced `hull == detected span` with "every member lies inside the hull",
  which a hull drawn from the members' cut spans also satisfies. No fixture puts a ragged end on
  the hull's edge, and the differential runs `resolve_features` on both sides.
- **Why it matters:** §5.1 chose the detected hull on purpose; the natural refactor — hull from
  `BundleMember`s — would turn the ragged bases of 385 HG002 bundles generic with every test green.
- **Suggested fix:** a fixture whose last member is eight `CAG` then `CA`; assert the hull ends on
  the `A` and the member two bases earlier.

**M4: [src/repeat_catalog/segments.rs](../../../../src/repeat_catalog/segments.rs#L300-L303) — the catalog route's no-clean-cut fallback never runs in a test**
- **Confidence:** Medium. **Categories:** reliability.
- **Problem:** `row.trimmed.unwrap_or(row.detected)` is pinned on the scan route by a unit test,
  but no differential contig holds a bundled row with no cut; replacing the fallback with an
  `expect` passes every test.
- **Suggested fix:** a differential contig with an impure bundled repeat whose `minimal_trim`
  fails, asserting the fixture holds such a row and that both routes agree.

**M5: [src/cli/generate_psps.rs](../../../../src/cli/generate_psps.rs#L217) and four more — "first repeat to last" misdefines the hull**
- **Confidence:** High. **Categories:** errors/naming/docs.
- **Problem:** the hull runs from the earliest detected start to the **furthest** detected end
  (`mod.rs:930-933`). The three subcommands' help, `type-regions`' help ("first member's start to
  last member's end"), `RepeatRouting`'s doc, the note `to_toml.rs` writes into every parameters
  file, the `SsrBundle` doc and the glossary say "to the last". Wrong for 520 of 2,410 HG002
  bundles. The specs' "385 of 2,410 end past their last member's cut end" is the count past
  *every* member's.
- **Suggested fix:** "from its first repeat's start to the furthest end of any of them"
  everywhere; correct the figure's wording; re-record the three parameters checksums in the same
  change, since the note is in the file.

**M6: [src/region_typing/mod.rs](../../../../src/region_typing/mod.rs#L802-L803) — "no hull built from the catalog's rows touches an end" is not guaranteed**
- **Confidence:** Medium. **Categories:** errors/docs.
- **Problem:** the builder's flank floor is measured on the **cut** span
  (`repeat_catalog/row.rs:104-112`), the hull on **detected** spans, and `--min-flank-bp` can be
  1. The same claim is in both specs. The code applies the rule on both routes, so only the docs
  are wrong.
- **Suggested fix:** say both routes apply the rule and that it rarely fires on the catalog's.

### Minor

- **Mi1** `src/region_typing/mod.rs:344` — a corrupt catalog row panics in `BundleMember::new`
  (the decoder checks no span), where a locus from the same row degrades to `NoCleanTrim`.
  Pre-existing in kind: `detected_interval` already underflows on a stored start of 0. Fix: validate
  spans in the decoder.
- **Mi2** `src/cli/typed_regions.rs:392` — `expect("motif bytes are ASCII")`; `Motif::new` checks
  only the length. Low confidence: depends on whether the FASTA reader admits non-ASCII bytes.
- **Mi3** `src/repeat_catalog/segments.rs:252` — the lookup's `expect` states its invariant only in
  the message. Fix: a `PANIC-FREE` comment naming why every key is present.
- **Mi4** `src/cli/run_ground.rs:318` — `CriteriaSource::TheFiveFlags` documents six flags.
- **Mi5** `src/region_typing/mod.rs:972` — the right-hand contig-end rule has no one-base-in case;
  `hull.end + 1 >= contig_len` would pass.
- **Mi6** `src/psp/segmentation_section.rs:548` — the `i64::MAX` refusal of the new key is never
  run; `--max-str-bundle-len 18446744073709551615` is a legal flag.
- **Mi7** `tests/ng_repeat_catalog_differential.rs:90-95, 820-825` — the scan side still ends in
  `..TypedRegionConfig::default()`, so a future config field defaults on one side only.
- **Mi8** `src/region_typing/mod.rs:972` — a dropped contig-end bundle is not charged to
  `rejected_by_reason` (question 2).
- **Mi9** `doc/devel/ng/spec/{str_bundle_loci,parameters_file}.md` — the census digests "move", but
  not what follows: a parameters file fitted before the cap is **demoted** to `supplied` by a run
  whose census was built after, although bundles never enter the census.
- **Mi10** `src/cli/cross_platform_digests.rs:130-136` (×3) — "its note two sentences" is garbled;
  the note gained a clause and a count.
- **Mi11** three `Args` structs — the flag's doc and `#[arg]` are copied three times; a shared
  `#[command(flatten)]` struct would hold the six routing flags once.
- **Mi12** `builder.rs`, `parquet_file.rs`, `cli/repeat_catalog.rs`, `criteria.rs` — "a catalog
  records the default cap" restated four times; the builder silently overwrites a caller's cap.
- **Mi13** `src/region_typing/mod.rs:331` — `BundleMember` holds bare `u64`s where the hull holds
  `Position`s.
- **Mi14** `src/cli/typed_regions.rs:384` — `members_json` is `pub` in a CLI module so an example
  can call it; it belongs beside `BundleMember`.
- **Mi15** `src/region_typing/mod.rs:845` — `resolve_features` takes seven positional parameters.
- **Mi16** no test pins member order where one member contains another (tied to M2).
- **Mi17** `src/region_typing/mod.rs:396` — the field doc does not name
  `DEFAULT_MAX_STR_BUNDLE_LEN`.

### Nits

The psp refusal message ("its typing is not one this build can name") is not quite true of the
catalog section it also guards; the `the_hull_and_its_first_member_agree_on_the_first_base`
message states the old invariant; backticked `` `Satellite` `` in user help; a doc reflow in
`segmentation_section.rs:630`; `crate::region_typing::BundleMember` spelled in full in
`segment_criteria.rs`; `typed_regions_cli.md`'s T4 still reads as current. Old builds ignore the
new psp key (no `deny_unknown_fields` on `WireRepeatCriteria`); the format-version bump of
`str_bundle_loci.md` §6.2 is where that closes.

## 7. Out of scope observations

- `repeat_catalog/reader.rs` `rows_for`: after M1's fix, a satellite built from a chain or a run of
  abutting short tracts reaching further than the padded window is still cut at it, so its extent
  depends on the BED (the probe: `[200, 346]` from the window against `[200, 380]` whole). An exact
  fix widens the read until no fetched row lies within the radius of the window's edge.
- `RegionKind::SsrBundle` can be built with 0 or 1 members (`run/callers.rs:4673`,
  `locus_generation/mod.rs` test); nothing enforces `>= 2`.
- Clippy `-D warnings` stops on two `type_complexity` errors in untouched files
  (`locus_generation/ssr.rs:2117`, `run/edge_deletion_gathering.rs:201`).

## 8. Missing tests to add now

- `a_window_over_one_end_of_a_bundle_types_the_whole_chain` (M1) — added.
- `the_hull_ends_past_a_ragged_last_member` (M3) — added.
- `genome_segments_bundle_member_without_clean_cut_keeps_detected_span_on_both_routes` (M4) —
  a bundled catalog row with `trimmed: None`; catches a dropped or panicking fallback.
- right-end one-base-in case in `a_bundle_at_a_contig_end_is_dropped_unless_capped` (Mi5) — added.
- `max-str-bundle-len-bp = u64::MAX` refused (Mi6) — added to
  `the_rules_refuse_each_broken_shape_naming_its_field`.
- `bundle_members_keep_detector_order_when_one_member_contains_another` (Mi16).

## 9. What's good

- Both routes describe a member through one shared pair of functions (`upper`, `minimal_trim`) and
  the differential compares the result, with a ragged member in the fixture
  (`tests/ng_repeat_catalog_differential.rs`).
- A psp with no cap is refused naming the key rather than read as 100
  (`psp/segmentation_section.rs`), while a parameters file with no cap is compared on the other
  axes only — two different answers for two different kinds of record, each stated.
- The cap test measures the hull from an uncapped run instead of pinning the detector's phase
  (`a_bundle_longer_than_the_cap_is_a_satellite_and_is_counted_as_one`).

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --no-fail-fast --lib`
- `scripts/dev.sh cargo test --release --test ng_repeat_catalog_differential`
- `tmp/review_2026-10-11_str-bundles-typing/verify3.sh` under `DEV_EXTRA_MOUNT` (HG002 typing and
  the stuttering benchmark's typing check).
