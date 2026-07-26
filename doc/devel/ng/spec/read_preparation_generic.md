# ng — read preparation (generic path) — RETIRED, see [`read_preparation.md`](read_preparation.md)

*Retired 2026-07-25. Read preparation turned out to be a **generic-path-only** step (the STR path has
no read preparation — see below), so the "generic vs STR" split that made this a separate file
dissolved. Its content — the three modes, `PreparedRead`, the `process_read` port, the reuse map and
parity oracle — merged into the single [`read_preparation.md`](read_preparation.md). Nothing was lost;
this file is a redirect so older links resolve.*

**Why there is only one read-preparation spec.** Read preparation canonicalises the line-up the mapper
gave a read. The STR path throws that line-up away — it re-aligns every spanning read against the
tract — so canonicalizing the mapper's CIGAR there would be work nothing reads. Its per-read operation
is also a different kind of thing: it produces an observation *about one locus*, not a read. It lives
in [`locus_generation_ssr.md`](locus_generation_ssr.md) (the generator) and
[`alignment.md`](alignment.md) (the repeat-aware aligner it calls). See
[`read_preparation.md`](read_preparation.md) §1.

**Where each part went:** all of it → [`read_preparation.md`](read_preparation.md).
