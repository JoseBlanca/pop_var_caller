# A2 — a sample lends its sections and does not hand them over

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step A2 — implementation plan 2.
**Design authority:** [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.1, §2.2; [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§6.2.
**Date:** 2026-08-14.

---

## 1. What this step is for

**A caller can no longer hold a piece of a sample's census.** `SampleCensusEvidence`'s two public
maps are gone; what a caller gets is `with_generic` and `with_strata`, which decode or borrow the
sections asked for, hand a closure borrows of them, and take them back when it returns. Spec §6.2
is explicit that this has to be a property rather than a convention: a call that *returned* a
section would let a file-backed run reassemble the whole file in memory, which is the outcome the
by-section design exists to prevent.

**And the repeat-tract evidence is split per stratum**, which arch §2.2 calls "the one storage-type
change this decision forces". A section is *one read group's tracts for one stratum*, so the
per-tract vectors are indexed within the stratum and the two counts that used to be keyed by
stratum are one number each.

## 2. What changed

| | before | after |
|---|---|---|
| `SsrEvidence` | one read group, every stratum | one read group, **one stratum** |
| its `covering_not_crossing`, `bases_compared` | `BTreeMap<Stratum, u64>` | one `u64` each |
| `GuardObservation::locus`, `TractDifference::locus` | index into the whole STR selection | index **within the stratum** |
| `SampleCensusEvidence.generic`, `.ssr` | two public maps | one private `Sections` |
| access | `sample.generic[&group]` | `sample.with_generic(&groups, \|sections\| …)` |

**The writer builds the stratum-local index beside the stratum it already carried.**
`CensusWriter::new` walks the kept tracts in genome order and numbers each within its own stratum,
so a section's tracts run in the same order the selection does — which matters because the fit sums
over a stratum's tracts and a sum of floats is not associative.

**Every declared read group gets a section for every stratum, whether or not a read reached one.**
That is what makes the enumeration the same in every sample, and §4 measures what it cost.

### Deviations from the architecture's signatures, and why

Both are recorded here rather than raised, because neither changes what the design promises — a
section is lent for the length of a call and nothing may retain one.

1. **`with_generic` takes a band of read groups, where arch §2.2 takes one.** A position's depth is
   the sum of a sample's read groups' depths ([`fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs)
   reads them together at every position), so a call lending one group at a time would leave the fit
   holding a partial sum for two million positions while it asked for the next. It is the same
   argument arch §2.2 already makes for lending a *band* of strata rather than one.
2. **The scoped calls return `R`, not `Result<R, CensusError>`.** Nothing can fail while every
   section is resident: `CensusError`'s two variants (arch §2.3) are `Malformed` and `Io`, both of
   which need the file. Writing `Result` now would put an `.expect` at every call site for an error
   that cannot occur. **The `Result` arrives with `Sections::Backed` in step B2**, which is the step
   that can fail.

A third thing is a bridge rather than a deviation: `generic_sections` and `ssr_sections` are
`pub(super)`, so `fit.rs`, `contamination.rs` and `ssr_fit.rs` can still read across every sample at
once while they are moved onto the cohort's own scoped calls in A3 and A4. **Nothing outside the
module has that door**, which is what the step was for; the four examples all go through
`with_generic` / `with_strata`.

## 3. What moved on real reads

Both oracles, before and after, at four container CPUs
(`tmp/run_oracle.sh`, `tmp/a2_before_*.txt` against `tmp/a2_*.txt`).

**No fitted number moved, on either cohort.** The whole diff is the size accounting, on three lines
a sample:

| | before | after |
|---|---:|---:|
| tomato, 8 accessions | 20.9 – 21.3 bytes a tract | **19.4**, every sample |
| the GIAB trio | 25.2 bytes a tract | **20.6**, every sample |

**The saving is the pair of map entries a stratum used to cost.** A `BTreeMap<Stratum, u64>` entry
holds a 16-byte key and an 8-byte value, and there were two such maps: 48 bytes a charged stratum,
now 16 — the two `u64` fields. Tomato charges 20 to 23 of its 26 strata and holds 349 tracts, so
about 1.7 bytes a tract; the trio charges all 32 of its own over 216 tracts, so 4.7. Both are what
the diff shows, to within the walked bits' own growth: a bit array a stratum rounds up to a byte
26 or 32 times instead of once, which is 13 bytes on tomato and 16 on the trio — 0.04 and 0.07
bytes a tract.

**The spread across tomato's accessions is gone, and that is the property rather than a saving.**
Before, a sample paid for the strata it happened to charge (20.9 to 21.3 across the eight); now
every sample holds every stratum's section and the number is 19.4 for all of them. Two samples
enumerating different sections would be two samples answering different questions.

*(At the scale the census actually runs at these terms barely register: tomato SL4.00 holds 462,701
kept tracts in 141 strata, so 141 sections cost 2,256 bytes of counts against 8.4 MB of offsets —
0.005 bytes a tract. The measurable difference above is a property of a 349-tract census.)*

## 4. Tests

| test | what it pins |
|---|---|
| `a_tract_is_numbered_within_its_stratum_and_not_across_the_selection` | **new.** Two twelve-base tracts of different motifs, one tract in each stratum: the second is locus 1 of the selection and locus 0 of its own stratum, so a writer keeping the genome-wide index writes a mismatch and a guard entry against a tract the section does not have |
| `two_strata_are_charged_separately_for_their_censored_reads_and_compared_bases` | rewritten: the two counts are now read off two *sections* rather than two keys of one map |
| `a_read_that_reached_a_tract_without_crossing_it_is_counted_against_its_stratum` | unchanged in what it asserts; the "no other stratum was charged" half is now structural — a section has nowhere to put another stratum's count |
| `read_groups_fold_by_addition` | built through `SampleCensusEvidence::resident` |

**The mutation check, run by hand.** Putting the genome-wide index back in the two `locus` fields of
`add_ssr` — the defect the new test exists for — fails
`a_tract_is_numbered_within_its_stratum_and_not_across_the_selection` and nothing else in the census
filter: §5 records the counts.

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo test --lib ng::parameter_estimation::joint::census` | `41 passed; 0 failed` (40 before) |
| the same, with the stratum-local index mutated back to the genome-wide one | `40 passed; 1 failed` |
| `cargo test --lib` | `3,590 passed; 0 failed; 11 ignored` (3,589 before) |
| the 88-second tomato oracle | §3 — no fitted number differs |
| the 74-second trio oracle | §3 — no fitted number differs |

**The two red gates are the two that were red before this branch's first commit**, neither in code
this plan touches: `cargo clippy --all-targets --all-features -- -D warnings`, and
`cargo test --all-targets`, which panics in `benches/psp_writer_perf.rs:386`.

**Nothing this step wrote is in the clippy gate's red set**, checked per target rather than by
reading the aggregate — the run aborts on the first target that fails, so which errors it prints
varies between runs. `cargo clippy --lib --all-features -- -D warnings` is **clean**, and so are
three of the four examples touched; `ng_joint_duplicated_in_fit` fails on a `needless_range_loop`
that `git stash` shows is there without this step's changes.
