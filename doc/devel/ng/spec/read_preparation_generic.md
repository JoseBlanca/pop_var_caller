# ng — read preparation (generic path) — RETIRED, see [`read_preparation.md`](read_preparation.md)

*Retired 2026-07-25. Read preparation turned out to be a **generic-path-only** step (the STR path has
no read preparation — see below), so the "generic vs STR" split that made this a separate file
dissolved. Its content — the three modes, `PreparedRead`, the `process_read` port, the reuse map and
parity oracle — merged into the single [`read_preparation.md`](read_preparation.md). Nothing was lost;
this file is a redirect so older links resolve.*

**Why there is only one read-preparation spec.** Read preparation is a per-read, *locus-independent*
transform: it canonicalises a read against the reference around its own span and produces a
`PreparedRead` that serves every locus the read overlaps. The STR per-read operation is not that — it
aligns a read *against a specific tract* to read out an observation *about that locus*, so it is
observation generation, not read preparation. It lives in
[`locus_generation_ssr.md`](locus_generation_ssr.md) (the generator) and
[`alignment.md`](alignment.md) (the repeat-aware aligner it calls).

**Where each part went:** all of it → [`read_preparation.md`](read_preparation.md).
