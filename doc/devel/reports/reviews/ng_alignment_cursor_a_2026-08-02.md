# Code Review: ng alignment cursor — Milestone A

**Date:** 2026-08-02
**Reviewer:** rust-code-review skill (orchestrator + five isolated sub-agents)
**Scope:** the working-tree diff of each Milestone A step, reviewed before its commit
**Status:** Approve-with-changes — all findings applied, see
[the fix report](../implementations/ng_alignment_cursor_a_fixes_2026-08-02.md)

*One file per milestone, a section per step. Per-category audit trail in the gitignored
`tmp/review_2026-08-02_ng-alignment-cursor-a1/`.*

---

## A1 — the probe, its tests, and the shared fixture builders

### 1. Scope

- **What:** the uncommitted working-tree change for step A1.
- **Against:** branch `ng-generic-perf` at `924e680` plus the step's patch. Each sub-agent
  worked in its own git worktree, detached at that commit with the patch applied — the
  change is not in any commit, so a worktree created from `main` would otherwise have
  reviewed nothing.
- **In scope:** `examples/ng_generic_walk_probe.rs` (new), `examples/shared/synthetic_alignment.rs`
  (new), `benches/ng_generic_pileup_perf.rs`, `Cargo.toml`.
- **Out of scope:** everything under `src/` — A1 changes no library code.
- **Categories dispatched:** `reliability` (the step's deliverable is a test module),
  `refactor_safety` (170 lines moved between targets, and the bench's numbers are only
  comparable across commits if the fixture did not move), `naming` (a new shared module's
  whole public surface, plus prose that makes checkable claims), `module_structure` (the
  placement decision), `tooling` (a manifest change and two feature gates against CI),
  `smells`. Not dispatched: `errors` (no error type added; `parse_count`'s `String` is
  pre-existing), `defaults`, `unsafe_concurrency`, `extras` — no trigger.

### 2. Verdict

**Approve-with-changes.** The move is proven safe; the test module as first written was not
adequate to its own stated job, and the Blocker below is the reason this loop exists.

### 3. Execution status

Run by the orchestrator on the host and passed verbatim into every sub-agent prompt:

| command | result |
|---|---|
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib --tests --examples --all-features` (CI's) | `2770 passed; 0 failed; 5 ignored`, plus `9 passed` on the example |
| `cargo test --lib ng::` (debug) | `1471 passed; 0 failed; 2 ignored` |
| `cargo bench --bench ng_generic_pileup_perf -- --test` | exit 0, six cases |
| the chr21 anchor | `loci=236081 observations=251786 reads_admitted=54709` |

Not run: `cargo audit`, `cargo doc --no-deps`. **`cargo test --release` is red on a clean
tree** for reasons predating this branch — four tests assert on `debug_assert!` messages
that release compiles out — so the debug profile is the bar, and every agent was told so.

### 4. Open questions and assumptions

1. **Resolved during the review, against the orchestrator.** The `Cargo.toml` comment
   claimed, and said it had checked, that `cargo test --all-targets` builds examples without
   running their tests. Three sub-agents independently measured the opposite. The
   orchestrator's own measurement had been invalidated by a **pre-existing panic in an
   unrelated criterion bench** (see *Out of scope*) that aborted the run before it reached
   any example. The comment was rewritten to what is true: `--all-targets` and `--examples`
   both run them; `test = true` earns its place only for a plain `cargo test`. *Affects M6.*
2. **The `compile_error!` swap** (the probe's allocator guard) was flagged for scrutiny
   rather than assumed correct. Both `smells` and `tooling` confirmed the replacement works
   and matches the repo's only comparable target, and both found the *justification* dishonest
   in the same way. *Affects M7.*

### 5. Top 3 priorities

1. **B1** — the typed-region walk could drop three quarters of its regions with every test
   green. The defect class the probe exists to detect, reproduced inside the detector.
2. **M1** — the dispatcher's region accounting was fetched and thrown away one line later:
   the three numbers that make "every region was walked" checkable.
3. **M6** — a comment asserting a cargo behaviour, claiming it was checked, and being wrong.

### 6. Findings

#### Blocker

**B1: `examples/ng_generic_walk_probe.rs` — the typed-region walk can drop three quarters of
its regions with all nine tests green.**
*Category: reliability. Confidence: High — mutated and measured.*
Truncating the typed-region stream to its first region takes the fixture from
`regions=4 loci=4979` to `regions=1 loci=1601` — 68 % of the loci gone — and the suite stays
green. The recorded `loci=236081` comes from this path, not from the whole-contig bound. A
walk that skips regions prints a smaller `loci` **and a shorter `seconds`**, which is
indistinguishable from the speed-up this branch exists to find. The only guard,
`the_typed_region_walk_reaches_the_generator`, asserted `> 0` and `<= span`, all satisfied by
one region.

#### Major

**M1: `walk` discarded `stream.counts()`.** `let counts = stream.counts(); let _ = counts;`
— `regions_in`, `regions_handled`, `loci_emitted`, computed by the library and dropped.
*Category: reliability, convergent with smells.*

**M2: only 3 of 16 printed keys were pinned.** Renaming three others at once was green, as
was crossing two counter assignments, as was making `generic_region_bp` count regions rather
than base pairs. The test was named *"the printed keys are the interface"* and covered a
fifth of it. *Category: reliability.*

**M3: `contig_filter` — the argument the anchor run passes — had no test at all.** Inverting
it and loosening it to a prefix were both green. The predicate was also written out twice,
so the two region sources could disagree about which contigs they walk. *Category:
reliability.*

**M4: `PVC_PROBE_MAX_RECORD_SPAN` could be dropped on the floor with the suite green.**
*Category: reliability.* **Partly disputed on measurement** — see the fix report: the knob is
genuinely unobservable on this fixture, and the proposed test asserted an effect that does
not occur.

**M5: the shared fixture module had no tests, and its central claim was false at low depth.**
Span 5,000 / read length 150 / coverage 1 covers **4,950 of 5,000** positions. `read_len >
span` panicked with a bare `attempt to subtract with overflow` inside a private helper — in
a release profile, a wrapped length reaching a slice instead. *Category: reliability,
convergent with refactor_safety.*

**M6: a `Cargo.toml` comment asserted a cargo behaviour, said it had been checked, and was
wrong.** *Categories: naming, smells, tooling — three independent measurements.*

**M7: the allocator guard's justification claimed more than the guard delivers.** It cites
the silent both-on mismatch, but fires only on a contradiction the operator typed; and the
error it names — a measurement lying about its allocator — remains available in the two
combinations it permits, because nothing in the output names the allocator. *Category:
smells.*

#### Minor

- **Mi1** `synthetic_reads` destructured the geometry with `..` while also reading the
  elided field indirectly. *refactor_safety.*
- **Mi2** `..ProbeReport::default()` in the test whose job is to pin the report's surface —
  the struct most likely to gain a counter during this work. *refactor_safety.*
- **Mi3** Four verified-false doc claims: "the same nine the dump prints" (it is eight),
  "three of which are `Option<u64>`" (two), `coverage` documented off by a factor of
  `read_len`, and a usage string naming two of four knobs. *naming.*
- **Mi4** `refuse_an_ambiguous_allocator` — a `-> bool` predicate named as a command.
  *naming.*
- **Mi5** `PVC_PROBE_WHOLE_CONTIG=0` turned the knob **on**, while the three count knobs go
  to lengths to reject a zero. *smells.*
- **Mi6** `split_generic`'s duplication was justified by "an example cannot import another
  example" — refuted by the shared module the same change introduces. *smells, convergent
  with module_structure.*
- **Mi7** `dhat::Profiler` was created before argument validation, so a mistyped command
  wrote a `dhat-heap.json` over a real profile. Reproduced: 5,352 bytes. *tooling.*
- **Mi8** The `alloc-mimalloc` feature doc promised a compile error no target has ever had.
  *tooling.*
- **Mi9** Test helpers took positional bools selecting between two different pipelines —
  the shape `ProbeRun`'s own doc says it exists to prevent. *smells.*
- **Mi10** `parse_count` had no boundary or malformed-input coverage. *reliability.*

#### Nits

Literal indexing after a length guard; `..generic.clone()` struct-update in a test;
`pieces[0]` after `assert_eq!(len, 1)`; the bench holding its span twice; a stray double
blank line in the manifest.

### 7. Out of scope observations

- **`benches/psp_writer_perf.rs:386` panics under `cargo test --all-targets`** —
  `index out of bounds: the len is 3300000 but the index is 3300000`, in
  `psp_writer_phases/flush_block_one`. Pre-existing: the file is untouched by this branch and
  last changed in `21dac5b`. It aborts `--all-targets` before any example runs, which is what
  produced the orchestrator's wrong measurement in open question 1. CI does not run
  `--all-targets`, so it is invisible there. **Worth its own fix.**
- **A third copy of the fixture primitives exists.** `examples/ng_generic_loci_dump.rs`
  carries its own `write_fasta` (line-for-line identical to the shared file's) and
  `write_bam` in its test module. The seam was drawn at the composed `SyntheticSample`, which
  the dump cannot use. Converting the dump is Milestone F1's business, not A1's.
- **`ProbeReport::render` has no compile-time link to the field set** — a counter added to
  the struct is printed only if someone remembers a `writeln!`. Mitigated by the new key-set
  test, not removed.

### 8. Missing tests added now

Every one is in the fix report's table. The suite went from 9 tests to 19.

### 9. What's good

- **The fixture parity harness.** `refactor_safety` did not assert the move was safe — it
  built the pre-move body beside the post-move one and compared FASTA bytes, whole BAM byte
  streams, and every record field at all three bench geometries (6,666 / 20,000 / 66,666
  records, zero differences). That is the difference between "the bench should still be
  comparable" and knowing it is.
- **Verifying prose as if it were code.** `naming` checked every claim these files make
  against the thing claimed, and found five wrong — including one of the orchestrator's.
  The Milestone D review found the same category was the highest-yield one; it still is.
- **Mutation over reading.** `reliability` applied 16 mutations and reported the eight the
  suite caught alongside the eight it did not, which is what made B1 a measurement rather
  than a worry.

### 10. Commands to re-verify

```
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test                                    # now reaches the probe's 19 tests
cargo test --lib --tests --examples --all-features
cargo bench --bench ng_generic_pileup_perf -- --test
./target/release/examples/ng_generic_walk_probe <ref.fna> <HG002 30x bam> chr21
```

---

## A2 — `CursorError` and `AlignmentFile::contigs()`

### 1. Scope

The uncommitted working-tree change for step A2, against `a400f73` plus its patch, each
sub-agent in its own worktree. In scope: `src/ng/read/input/cursor.rs` (new),
`src/ng/read/input/open_bam.rs` (new field, accessor, `Debug`, five tests),
`src/ng/read/input/mod.rs` (one line). Out of scope: the rest of `open_bam.rs`.

**Three agents over five checklists**, paired because the diff is small: `reliability`;
`errors` + `naming`; `module_structure` + `smells`. Not dispatched: `defaults`,
`unsafe_concurrency`, `extras`, `refactor_safety`, `idiomatic` — no trigger on a types-only
addition.

### 2. Verdict

**Approve-with-changes.** All findings applied. The step's central Major was found
independently by all three agents.

### 3. Execution status

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib ng::` `1476 passed` at review time. `cargo doc --no-deps --lib` exits 101
— **13 errors pre-date the branch**, 2 were new here.

### 5. Top 3 priorities

1. **M1** — the docs say `contigs` is the reference's list; it is the file's, and the test
   named after the property cannot fail.
2. **M2** — `WrongChromosome` names two indices and no file.
3. **M3** — two mutations survive: `sq_md5s` built in reverse order, and the `Debug` contig
   count off by 999.

### 6. Findings

#### Major

**M1: the "reference's list" claim, in three places, none of them true.**
*Categories: naming, reliability, module_structure — three independent confirmations.*
The value is `contig_list(&header)` — a re-reading of the file's header — reconciled against
the reference under a rule that treats an absent `M5` as a wildcard. Proved twice by
experiment: a probe showed the file's digests surviving into `contigs()` against a `.fai`
reference carrying none (`PartialEq says equal: true`, structurally identical: `false`), and
a mutation populating the field from the reference left the test named after the claim
**passing**. Arch §2.1 and spec §8 carry the same sentence.

**M2: `WrongChromosome` carries no file identity.** `alignment_file.md` §4 states the house
pattern as "each variant carrying its path"; every `AlignmentFileError` variant has one, and
the sibling `Io` proves a cursor holds it. With 32–320 live cursors, two bare integers name
neither the file, the sample, nor the contig table they index. *Category: errors.*

**M3: two surviving mutations.** `sq_md5s` built in reverse order passed all three new tests
(the fixture's digests are all `None`, so every pair was `None == None`); the `Debug` contig
count `+ 999` passed everything. *Category: reliability.*

#### Minor

- **Mi1** `Io` names the mechanism, not the failed operation — against the rule the crate's
  other three error enums follow. *errors.*
- **Mi2** `sq_md5s_by_file`, cited as the justification for a duplicated field, **does not
  exist**. The real consumer is `SampleReads::assembly_inputs`. *naming, reliability,
  module_structure.*
- **Mi3** `(spec §4)` cites "Error model" / "What the cursor promises"; the claim's home is
  `alignment_file.md` §3.1 check 2. *naming.*
- **Mi4** 35,228 quoted without the mode qualifier spec §11.5 explicitly demands — the
  typed-region walk counts 34,633, and B3's counter assertions have to choose. *naming.*
- **Mi5** The module's rustdoc **summary line** describes a type the module does not contain;
  the "only the errors" disclosure sits 18 lines below the fold, where the index does not
  show it. *naming.*
- **Mi6** `cursor: ContigId` / `requested: ContigId` — two same-typed fields whose names do
  not disambiguate them. *naming.*
- **Mi7** Two `unresolved link` rustdoc errors from intra-doc-link brackets on filesystem
  paths. *tooling, via cross-category.*
- **Mi8** `contigs` and `sq_md5s` are duplicated state. **Deleting the field and having
  `check_assembly` take a `&ContigList` was tried and works** — suite green, ~18 lines
  shorter — but it changes a `pub` signature and eight call sites. *module_structure.*

### 7. Out of scope observations

- **`cargo doc` is red on this branch's parent** — 12–13 `unresolved link` errors in
  `src/ng/locus_generation/ssr.rs`, `src/ng/region_typing/`, `src/ng/tandem_repeat.rs`,
  `src/ng/types.rs` and five `src/ssr/` files. CI runs `cargo doc --no-deps --lib
  --all-features` with `RUSTDOCFLAGS: -D warnings`, so **the doc step fails independently of
  this work**. Worth its own fix.
- **A spec/arch contradiction that will bite at B1.** `spec/alignment_cursor.md` §10: *"The
  cursor is consumed when it rejects a region, so a half-valid cursor cannot exist."*
  `arch/alignment_cursor.md` §2.2: *"`move_to_region` leaves the cursor usable on any
  outcome."* Arch's own signature `move_to_region(&mut self)` cannot consume. A2 has no
  cursor, so nothing is wrong yet — **raised at Checkpoint A**.
- **`SampleReads::assembly_inputs` has no production caller** — only a test that `.count()`s
  it. Its whole reason for existing is `check_assembly`, which no shipping path invokes.

### 9. What's good

- **The same Major, reached three different ways** — by reading the gate's wildcard rule, by
  a probe printing both lists, and by mutating the field to match the docs and watching the
  wrong test pass. Convergence on a claim nobody had checked.
- **Recording the divergences that are *correct*.** One agent noted that arch §1.4's
  `#[error("reading {path}")]` would not compile and spec §8's `Io(std::io::Error)` carries
  no path, so the code is right to differ — flagging them so a later reader does not "fix"
  them back.

---

## A3 — `RecordReader` and its in-memory arm

### 1. Scope

The uncommitted change for A3 against `cd9fbd9` plus its patch. In scope:
`src/ng/read/input/record_reader/{mod,in_memory}.rs` (both new) and one line in
`src/ng/read/input/mod.rs`. **Two agents over four checklists**: `reliability`; and
`module_structure` + `naming` + `smells` paired.

### 2. Verdict

**Request-changes, then approve.** The step's own deliverable — an oracle — was the thing the
tests did not pin. **Thirteen mutations run, six survived.**

### 5. Top 3 priorities

1. **B1** — `begin_region`'s rewind is never observed from mid-script, so it can be deleted
   at either layer with the suite green. That is spec §7's abandoned-region case.
2. **B2** — the clone's fidelity is never checked beyond the read name; dropping
   `alignment_start`, the one field the forget rule compares, survives.
3. **M1** — the contract claims `read_next` reuses the buffer's allocations. The only arm
   does not, and the doc also claims the real arms have "the same cost shape".

### 6. Findings

#### Blocker

**B1: the rewind can be deleted and nothing fails.** Both tests exercising `begin_region`
called it on a reader already at position 0 or already drained. `if self.next_index >=
self.records.len() { self.next_index = 0 }` passed the whole suite — and that is exactly
spec §7's first named case: a caller that abandons a region leaves the reader mid-script, and
the next region would get a truncated one. *Category: reliability. Mutation-verified.*

**B2: the clone is only proved to carry the name.** `drain` compared name lists, so a clone
losing `alignment_start` — **the field the forget rule compares** — survived, as did one
keeping the name and garbage everywhere else. *Category: reliability. Mutation-verified.*

#### Major

**M1: "reusing the buffer's allocations" is false for the only arm.** `RecordBuf` derives
`Clone`, and a derived `Clone` gets the default `clone_from` (`*self = source.clone()`),
which drops the destination's buffers. Measured with a standalone program: `derive(Clone):
ptr_same=false cap_before=256 cap_after=64` against a manual impl's `ptr_same=true
cap_after=256`. The second claim — that a real arm has "the same cost shape" — is wrong the
other way: noodles' `read_record_buf` is genuine reuse. *Category: smells.*

**M2: the forwarding test cannot detect a skipped delegation**, though its doc says it can.
Replacing the enum's `begin_region` body with `Ok(())` left all eight tests green, because the
test repositioned a freshly built reader. *Category: smells. Mutation-verified.*

**M3: `other_sample_records`'s justification is untrue.** `RecordReader` does not implement
`RecordSource`, so the trait's no-default rule does not reach it; and "happened twice in this
module's history" relocated another module's incidents into a module created in this commit.
*Category: smells.*

**M4: the read-group guard was one record deep**, and `header()` was pinned only by contig
*count* — a header with the right number of contigs under different names would put every
read on the wrong chromosome. *Category: reliability. Both mutation-verified.*

#### Minor

- **Mi1** "stated once so two arms cannot drift" describes a hope; prose does not fail a
  build. *module_structure.*
- **Mi2** "A reader holds only its position" drops arch §1.3's **one-record pushback**, which
  plan C2 implements — so the contract was already scheduled to be false. *smells.*
- **Mi3** "a record replayed from memory" transposes spec §5's *read* into *record*; a
  replayed read sits above the filter and goes through **fewer** lines, not the same ones.
  *smells.*
- **Mi4** `next: usize` is a bare participle where the crate's two nearest analogues use
  `next_index`. *naming.*
- **Mi5** the four-line allow reason was duplicated verbatim in both files. *smells.*

### 7. Out of scope observations

- `src/ng/read/input/region_query.rs` says "The trait's default is `0`" of a trait method
  that **has no default** — the stale sentence the new doc echoed.
- Arch §1.3's own bullets contradict its prose about who keeps what.

### 9. What's good

- **Three of the brief's four judgement calls came out in the code's favour, with the reason
  measured rather than asserted** — including why `#[expect]` would *not* work here (the
  lint fires for the lib target and not the test target, so `--all-targets` reports the
  expectation unfulfilled and fails `-D warnings`).
- **The `clone_from` finding was settled with a standalone program**, not by reading noodles'
  source: `grep -rn "fn clone_from"` over noodles-sam returning nothing is suggestive; the
  pointer and capacity numbers are proof.
