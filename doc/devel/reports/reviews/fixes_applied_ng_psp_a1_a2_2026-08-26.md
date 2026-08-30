# ng psp store A1+A2 — fixes applied

*2026-08-26. Applied against [the review](ng_psp_a1_a2_2026-08-26.md) of `c4c163ee`, branch
`ng-psp-encoding`.*

---

## Accounting — every finding

| # | finding | severity | status |
|---|---|---|---|
| B1 | the offset `decode` returns was untested | Blocker | **Applied** |
| M1 | `deny_unknown_fields` refused a same-major file that added a key | Major ×3 categories | **Applied** |
| M2 | `encode` wrote a header its own reader refused (`> i64::MAX`) | Major | **Applied** |
| M3 | the writer's and reader's encoding tables were two lists | Major ×3 categories | **Applied** |
| M4 | `FOOTER_BYTES` did not pin the wire width | Major ×2 | **Applied** |
| M5 | six rules never exercised; three fixtures far from their boundary | Major ×2 | **Applied** |
| M6 | `check_basename` accepted `/` and `..` | Major ×3 | **Applied** |
| M7 | read errors had an empty `source()` chain | Major | **Applied** |
| M8 | a foreign file reported as a damaged one | Major | **Applied** |
| M9 | `Io` named neither the file nor the operation | Major | **Applied** |
| M10 | `OutOfOrder` did not name the file | Major | **Applied** |
| M11 | `RecordHead` states a shape E3 must break | Major | **Applied** (documentation, per the reviewer's own recommendation) |
| M12 | two spec defaults had no named constant; three limits unpinned | Major ×3 | **Applied** |
| M13 | no property or fuzz test for the parser | Major ×2 | **Applied** |
| M14 | TOML keys read as two different things four lines apart | Major | **Applied** |
| M15 | a test decoded the binary length prefix as UTF-8 | Major | **Applied** |
| — | double public path per type; constants outside the re-export list | Minor | **Applied** |
| — | `footer.rs` tied ng's width to production's `pub(crate)` constant | Minor | **Applied** |
| D1 | the 1 MiB header cap is an ~11,300-contig ceiling | Major | **Deferred — Checkpoint A** |
| D2 | `BlockIndexEntry`'s dropped field: spec §3.3/§4.1 against §6.2 | Minor | **Deferred — a spec contradiction, not a code fix** |
| D3 | the manifest carries no cardinality | Minor | **Deferred — already open from A1** |
| D4 | the compression level has no home in the manifest | Minor | **Deferred — Checkpoint A** |
| D5 | `toml::value::Datetime` in a public struct | cross-category | **Deferred** |
| — | a version overflowing `u16` classed as damage, not as unsupported | Minor | **Deferred** — the honest fix needs `UnsupportedVersion.found` to stop being `(u16, u16)`, which is an arch-level change. The message now says the parts must be numbers below 65535. |

**Nothing was disputed.** Every finding reproduced against the current code.

## Two adaptations worth naming

**M2's reader half cannot reach the rule, and the table says so.** A number above `i64::MAX` has no
TOML spelling a parser will take back, so the reader stops at the syntax before any rule runs. That
is not a gap — it is why the *writer* needs the rule. The test's two over-large rows assert the
writer's message and the reader's separately, with the reason written at the assertion.

**One reviewer suggestion was tried and reverted, and the reason is now a comment.** The errors and
reliability agents both noted that `decode` parses the body twice — once as a bare `toml::Table` for
the version, once into the wire types. Feeding the already-parsed table to the wire types instead
**breaks every header**: `writer.created` comes back out of a `toml::Table` as a string rather than
as a TOML datetime, and deserialisation fails. A header is about a kilobyte and is read once per
file open, so the second pass costs nothing that matters. The comment at the call site records that
the alternative was measured by trying it.

## Validation

Run in the container on the tree that was committed:

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp` | **37 passed, 0 failed** (26 before) |
| `cargo test --lib --bins --tests --examples` | **4,569 lib tests passed**, 0 failed, 14 ignored; every other target green |

`--all-targets` remains outside the gate: it is red on `main` for a pre-existing panic at
`benches/psp_writer_perf.rs:386`.

## The file format changed, and this was the moment for it

Four wire keys were renamed on the naming finding: `fixed` → `fixed-width-integer`, `ieee` →
`ieee-float`, `bytes` → `width-bytes`, `scale` → `steps-per-unit`, and `window-log` →
`look-back-window-log`. **No psp has ever been written**, so this costs nothing now and would cost a
format version later. The Rust names follow the wire ones, and `ContigEntry` became `ContigIdentity`,
which also removes its shadowing of `crate::fasta::ContigEntry`.

**This is a departure from [arch §3.3](../../ng/arch/psp_file_format.md), which names the variants
`Fixed`, `Ieee` and `FixedPoint { scale }`.** Raised at Checkpoint A with the rest.
