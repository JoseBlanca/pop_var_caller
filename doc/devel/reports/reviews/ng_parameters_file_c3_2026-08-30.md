# Review — C3, the float round-trip oracle

*2026-08-30. One agent in an isolated worktree at `4360ccb5`, handed the step's diff as a patch.
**0 Blockers, 2 Majors and 5 Minors applied.** Module tests 126 → 130.*

The step found no defect: every float already round-tripped. So the whole of its value is whether
the oracle would have caught one — and that is what the review was asked to judge.

## What it established, and what it turned out not to be

**Both formatters recover every double.** The pass traced them apart in the dependency: the
artefact's writer is `format!("{value:?}")`, Rust's shortest representation that round-trips;
`serde`'s is `Display`, which never uses exponent form, so `f64::MAX` is 22 characters through one
and 310 through the other. They are genuinely different text for the same number, so exercising
both is real coverage rather than one path twice. And there is one float decoder underneath —
`toml`'s `DeFloat::to_f64` — which both the shape's own reader and `toml::Value` reach, so the two
smaller tests measure the decoding a run uses.

**And the spec's proposed fix fits neither writer.** §4 says that if the digits were lost, "the fix
is a serializer that formats floats for round-trip". The writer that produces the artefact *is*
that serializer already, by the choice made at B2 — so that half was settled by construction and
what these tests add there is a guard. The one that had genuinely not been checked is `serde`'s,
and it writes a golden test file rather than the artefact, so a digit lost there would have been a
defect in a fixture. **Recorded as owed to the spec rather than edited into it**; the spec is the
owner's.

## The two Majors were both the header describing itself wrongly

**"Every assertion below compares `to_bits()`"** — the last test's principal assertion is `==` on
whole files through the derived `PartialEq`, which is exactly the comparison the paragraph beside
it says cannot be trusted. The test is sound: it plants a negative zero and checks that one on
bits. The header was wrong about it.

**"the module's round-trip tests read files written by each"** — no test in the module reads a
file. It calls `toml::to_string` in memory beside `to_toml`. The conclusion held; the reason named
an artefact that plays no part, which would have sent the next reader to the wrong place to check
it.

## Five Minors, of which two changed what is covered

**A negative subnormal was missing**, and the sweep is a poor net for that class: the subnormal
exponent field is 1 of 2,048, so about five of the ten thousand swept values are subnormal at all.
One line closes it, and a hand-rolled replacement for the formatter is exactly where a sign lost
off a subnormal would come from.

**The whole-file test's comment claimed a curve's coefficients** and set none — the curves kept
their fixture values throughout. One now carries a hard value.

The rest were prose: `1.0000000000000002` is one unit in the last place above one, not two; the
sweep's own coverage sentence read as though it were the subnormal net when the table is; and the
module header now says plainly which of the brief's three adversarial categories it does not cover.

## What the pass confirmed rather than found

The sweep cannot run out of draws — a `u64` pattern is non-finite exactly when its eleven exponent
bits are all ones, 1 in 2,048, so about 19 of the ~10,019 draws needed are skipped against a budget
of 40,000. The generator is Marsaglia's xorshift64 triple, full period from a non-zero seed, and
cannot degenerate. A failure is reproducible from its own message: the bits are printed in hex.
And no float in the file escapes the one formatter — every other `to_string` in the writer is on an
integer field.
