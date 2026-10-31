# Fuzzing pvc's input readers

**What this checks: a file a user hands pvc, however malformed, gives an error message and never a
panic.** Each target writes the fuzzer's bytes over one input of a small real run and runs the
`pvc` subcommands that read that input, in-process. An `Err` is the right answer; a panic, an
abort or an allocation over the limit is a finding, and libFuzzer saves the input that caused it.

| target | the input it replaces | subcommands run |
|---|---|---|
| `bam` | the first sample's BAM, given as the *decompressed* stream; the target compresses it, so mutations land on records and not on deflate bytes | `call-from-alignments --build-index-if-missing` |
| `cram` | the first sample's CRAM, raw | `call-from-alignments --build-index-if-missing` |
| `fasta` | the reference FASTA and its `.fai`, joined by a line holding only `@@FAI` | `repeat-catalog`, then `call-from-alignments` with that catalog |
| `bed` | `--regions` | `call-from-alignments` |
| `parameters` | `--parameters` | `call-from-alignments` |
| `psp` | one sample's psp | `call-from-psps`, `dev inspect-psp` |
| `psp_fit` | one sample's psp | `estimate-parameters` — seconds a run, so give it `-timeout=600` |

## The fixture

`fixture/` is a run small enough to repeat a hundred times a second: two 10 kb stretches of GRCh38
as contigs `ctgA` and `ctgB`, the GIAB 5x reads of HG002 and HG003 that fall wholly inside them, and
the catalog, psps and fitted parameters pvc built from those. `corpus/` holds each target's seeds.
Both are rebuilt by `make_fixture.sh`, inside the dev container, from the GIAB inputs:

```
DEV_EXTRA_MOUNT=<giab_trio dir> ./scripts/dev.sh bash fuzz/make_fixture.sh \
    <GRCh38.fa> <giab_trio dir> target-container/soak/pvc
```

## Running

cargo-fuzz needs nightly, which the dev container does not ship, and the container is discarded
after each command. So nightly and `cargo-fuzz` are installed under the repository's `tmp/fuzz/`,
once:

```
./scripts/dev.sh bash -c 'export RUSTUP_HOME=$PWD/tmp/fuzz/rustup CARGO_HOME=$PWD/tmp/fuzz/cargo
  rustup toolchain install nightly --profile minimal
  cargo-binstall --no-confirm --install-path $CARGO_HOME/bin cargo-fuzz'
```

and every run sets the same two variables and puts `$CARGO_HOME/bin` on `PATH`. The target must be
named: without `--target`, cargo-fuzz picks musl, where the address sanitizer does not link.

```
cargo +nightly fuzz run --fuzz-dir fuzz -a --target aarch64-unknown-linux-gnu psp \
    tmp/fuzz-corpus/psp fuzz/corpus/psp -- -max_total_time=1200 -rss_limit_mb=2048 \
    -artifact_prefix=tmp/fuzz-artifacts/psp/
```

`-a` keeps debug assertions on. Keep the working corpus and the artifacts under `tmp/`.

**The `bam` target's memory grows by about 35 MB a run**, held across runs in one process, and
reaches 2 GB in about 55 runs; give it `-rss_limit_mb=9000` and keep `-malloc_limit_mb=2048`, so a
single oversized allocation is still caught.
