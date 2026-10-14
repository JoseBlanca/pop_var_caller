# Calling STR bundles — the measurement that decides how

*Plan, 2026-10-11. Branch `str-bundles`, worktree `/Users/jose/devel/pop_var_caller-str-bundles`,
from main at a87ce7b6 (main then held the STR edge-deletion work,
[`../../reports/implementations/str_edge_deletion_2026-10-11.md`](../../reports/implementations/str_edge_deletion_2026-10-11.md)).
Nothing in `src/` is changed on this branch yet.*

---

## 1. The gap

**An STR bundle is two or more repeats within 15 bases of each other**, so none of them has the
clean flanks an STR segment needs (`typed_regions.md` §2.4; the radius is
`DEFAULT_BUNDLE_THRESHOLD`, 15). Its region is the hull of its repeats, first start to last end, so
the bundle as a whole does have clean flanks outside its outermost repeats; it lacks them only
between members.

**Nothing calls a bundle today.** The typing emits `RegionKind::SsrBundle`; the walk's generator
slot for bundles is empty in production (`src/run/walker.rs`, `GeneratorSet::new` near line 1700:
`GeneratorSlot::Unfilled(UnhandledReason::NotImplemented)`), so their reads are never measured
(dispatch in `src/locus_generation/mod.rs`, `begin_region` and `next_locus`); the generic walk stops
at them; and calling would set a bundle locus aside if one arrived (`src/run/callers.rs`,
`LocusOutcome::BundleSetAside`).

**Its size, on GIAB HG002's tandem-repeat benchmark regions** (the 50,000-region tier set; typed
with main's `ng_typed_region_dump` at default calling floors; `scripts/ng_str_bundles/`):

| | STR segments | bundles |
|---|---|---|
| regions | 20,183 | 2,473 |
| truth variant records starting inside them | 8,534 | 2,422 |

So about one truth variant per bundle against about 0.4 per segment: bundles are the more variable
ground. 2,207 of the 2,473 bundles have two members, 229 three, 37 more; the median bundle is 35
bases long, 90% are 71 or shorter, the longest 329. Another 23,118 truth records start in what ng
types as ordinary sequence — mostly repeats with too few copies to be STR segments — which the
generic path already calls.

## 2. The decision (owner, 2026-10-11)

**Measure the cheap alternative first — let the generic path call bundles — and build an STR-path
version only if the measurement says we must.** The three options considered:

1. **One compound STR locus per bundle**, analogous to segments: each read crossing the hull is
   spelled over it, anchored on the outer flanks. Needs a motif per member in the read model (a slip
   can be a unit of either repeat; the model's per-run slip placement for interrupted repeats is the
   nearest existing piece), one repeat count per member in the prior's ladder and in `REPCN`.
   Buys stutter modelled in every member, phase between members, one record. About a week, a guess.
2. **The generic path over the bundle** — the option to measure. Small change; no stutter model (PCR
   stutter appears as low-fraction indels the candidate bar must reject), no phase between members.
3. **Each member as its own segment with the gap as its flank** — ruled out: the gap is under 15
   bases by definition, the STR path's anchoring fails there, and a deletion spanning the gap is the
   edge-deletion problem on both sides.

## 3. The measurement

**Build.** Two ways to let the generic path walk bundles; I recommend the first:

- **Type a bundle's hull as ordinary sequence**, coalesced with the ordinary regions either side, in
  the typing partition (`src/region_typing/mod.rs`, where bundles are clustered and emitted). The
  walk's invariant that two generic regions are never adjacent (`genome_walk.rs`, near lines
  708–713) then holds by construction, and the pileup's edge handling sees nothing new.
- Fill the bundle generator slot with the generic generator. Keeps the typing, but puts three
  generic regions side by side, which the walk has never seen.

Either way the psp's recorded typing changes, so a run needs fresh psps and the
`cross_platform_digests` VCF checksums may move (re-record only after showing the change alone moves
them).

**Score.** `scripts/ng_str_accuracy/hg002_two_binaries.sh` runs main's binary and the experiment's
on HG002 at 30× and 50× and scores both (`benchmarks/lib/tract_qual_experiment.py`, which reads `*`
as no change since 9a6709d4). Two grounds:

- **`TIER_KIND=ssr_bundle`** — the bundles. Main writes nothing there, so main's row is the floor.
  **The ground must come from main's typing**: the script types with its own checkout's dump, which in
  the experiment no longer has bundles — build `tier.bed` with main's dump, copy it into the output
  directory, and the script reuses it. A bundle has no single period; `PERIOD_FOR_BUNDLES` (default
  2) sets the class it is scored in.
- **`TIER_KIND=ssr_locus`** — the segments, to show nothing there moved.

**What to compare it with.** On segments main gets 3,503 of 3,995 homopolymer tracts and 2,566 of
2,980 longer-period tracts right at 30× (3,595 and 2,613 at 50×); the scorer's accuracy column
divides by the tracts both the truth and the caller call (about 90%). Useful calibration, not built
yet: the generic path's accuracy on the benchmark tracts ng types as ordinary sequence (the 23,118
records above) — the ground the generic path already calls.

**Decision rule.** If bundles called this way come out close to the STR path on segments, option 2
is the answer. If well below, that number is what option 1 is worth, and option 1 gets a design.

**What it cannot show.** HG002 is one sample, high depth, PCR-free. The stuttering benchmark places
no bundles (its simulator picks repeats with clean flanks), so PCR stutter in bundles and cohort
sizes are unmeasured; adding bundles to the simulator is a separate piece of work, worth proposing
only if option 2 looks good on HG002 (`CLAUDE.md`: one sample to thousands, 3× to hundreds).

## 4. Outcome (2026-10-11): the generic path fails under PCR stutter; bundles get an STR bundle path

**Decision (owner, 2026-10-11): build an STR bundle path — one STR locus per bundle, a motif per
member — with the same 100-base cap as STR segments; longer bundles are typed as satellites.** The
generic-path change was not merged. It is kept as a patch, not committed.

**On HG002 the generic path looked good** (bundles typed as generic, coalesced with the generic
regions either side, in `resolve_features`). Per region, at 30×, against main's typing:

| HG002, 30× | truth variant: right | wrong genotype | missed | no truth variant: false call | left as reference |
|---|---|---|---|---|---|
| bundles, main (nothing calls them) | 55 | 94 | 976 | 1 | 1,341 |
| bundles, generic path | 901 | 178 | 46 | 16 | 1,326 |
| STR segments, period 2+, STR path | 2,566 | 299 | 115 | 24 | 4,793 |
| STR segments, homopolymer, STR path | 3,503 | 361 | 131 | 73 | 8,217 |

At 50×: bundles on the generic path 934 / 164 / 27 / 18. The STR segments' counts were identical with
and without the change at both depths. In 139 of the 178 wrong bundles the right allele was never
offered.

**HG002 is PCR-free, so the owner asked whether that was why.** It was. The same build was run on the
stuttering benchmark's reads (2026-10-10 first pass, true stutter, no fit) twice: as is, and with
every STR segment typed as generic (an experiment-only switch, `NG_EXPERIMENT_SEGMENTS_AS_GENERIC`).
A cell is one sample at one segment; false calls are among the cells where the sample carries no
variant (about 87,800 at 50 samples):

| 50 samples | variant cells right, STR path | generic path | false calls, STR path | generic path |
|---|---|---|---|---|
| PCR 30× | 61,370 of 62,127 | 58,477 of 62,290 | 79 | 9,097 |
| PCR 5× | 43,569 of 62,620 | 40,097 of 64,122 | 1,557 | 8,563 |
| PCR-free 30× | 61,803 of 62,127 | 61,055 of 62,154 | 11 | 560 |

GATK on the same PCR 30× reads: 9,448 false calls. 8,279 of the generic path's 9,097 sit on mono-
and dinucleotides of 16 copies or more (17,831 cells). At one sample, PCR 30×: 10 false calls
against 3 of 1,771. On HG002, 293 of the 2,473 bundles hold a mono- or dinucleotide run of 16 copies
or more, and 795 one of 12 or more — the ground where the generic path takes stutter for alleles.

**Why not each member as its own STR locus** (owner's question, after the size worry): the STR path
aligns a read against the repeat and 15 repeat-free bases either side, and inside a bundle that flank
is partly the neighbouring repeat (`SsrGeneratorConfig::check_flank_within`). A length change in the
neighbour moves the anchor; a slip in `(CA)8 GT (CA)10` can be counted by both loci; a deletion
across the gap is the edge-deletion problem between two STR loci. **The size worry is met by the
cap:** 2,410 of the 2,473 bundles are 100 bases or shorter; 10 are longer than 150.

## 5. Practical notes for the session

- **Rules that still hold** (owner): don't run the parameter fit; where STR calling is measured on
  the stuttering benchmark, call with the true stutter (`true_stutter.parameters.toml`). Don't
  reintroduce an 8-base overhang rule for crossing reads.
- **The dev container** mounts only its own checkout; HG002 needs
  `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks`. `scripts/dev.sh` does not pass the
  host environment: set variables inside it (`scripts/dev.sh env VAR=… cmd`). Long `bash -c` strings
  can be blocked: write a script under `tmp/` and run `scripts/dev.sh bash tmp/that.sh`. A script
  under `tmp/` must not find the repository from `dirname "$0"/../..` (that is the parent of the
  worktree); scripts under `scripts/` can.
- **Host:** zsh (an unquoted `$VAR` of several words is not split); no `python3` (use
  `uv run --no-project python`); `rm -rf` is declined — delete your own scratch from inside the
  container.
- **Builds:** the binary is `pvc`; copy each measured build to a stable place before the next build;
  never edit source while it compiles. `cargo test --release` fails three debug-only `should_panic`
  tests; use `--no-fail-fast` (a debug `cargo test --no-fail-fast` passes everything at a87ce7b6).
- **Disk** was 98% full during the last session. `--explain-loci` outputs are large (8.6 GB for 7
  builds over 16 runs); delete scratch you no longer compare against. Other sessions may be running
  in other worktrees (one was compiling in `stutter-multiplier` on 2026-10-11): check before deleting
  anything outside this worktree, and ask before stopping a container that is not yours.
- **The stuttering benchmark's source runs** are in
  `/Users/jose/devel/pop_var_caller-tract-ceiling/tmp/stutter_bench/first_pass/` (reads, truth,
  GATK and HipSTR calls; `scripts/ng_str_accuracy/copy_run.sh` copies from there). Not needed for this
  measurement, but needed by any follow-up on cohorts and PCR.
