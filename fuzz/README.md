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

**`regressions/` holds every input that once crashed pvc**, one directory a target, named for
what it broke. Pass that directory as a further corpus (`fuzz/regressions/psp` beside
`fuzz/corpus/psp`) so each run checks them first; `make_fixture.sh` rebuilds `corpus/` and never
touches it. Those found on 2026-10-31, all now refused with an error:

| file | what it did |
|---|---|
| `psp/crashed_overlapping_blocks` | blocks copied over each other, each passing its checksum; two records on one base aborted `call-from-psps` in the merge |
| `psp/crashed_census_depth_bin` | seven bytes of the census changed; a depth code past the ladder aborted `estimate-parameters` |
| `psp/crashed_renamed_header` | the header's sample renamed by hand; the census's name disagreed and aborted `estimate-parameters` |
| `parameters/crashed_slippage_group` | `slippage_group = 4000000000` sized a table and the run was killed for memory |
| `bam/crashed_reference_number` | a record naming reference sequence 2.6 billion, in a header listing two; building the index asked for 316 GB |
| `bam/aborted_under_a_cap_reference_count` | a header declaring 2^31 reference sequences; noodles-bam reserved 38 GB, and pvc aborted under an 8 GB cap |
| `bam/aborted_under_a_cap_record_size` | a record declaring 4.26 GB; noodles-bam reserved it, and pvc aborted under a 4 GB cap |

**The `bam` target looked as if it leaked, and it did not.** Its memory rose by about 35 MB a run
and reached 2 GB in about 55 runs, while one valid BAM replayed 150 times in one process stayed at
437 MB. What rose was the high-water mark: noodles-bam reserved what a BAM header's
reference-sequence count, a reference name's length or a record's size declared before reading
the bytes, so a mutated input asked for gigabytes — and under a memory cap (`ulimit -v`, a
cluster's `h_vmem`) that aborted pvc itself. The fork in `vendor/noodles-bam` reads all three as
the bytes arrive (its `FORK.md`). Run the target in fork mode anyway, so one input's memory stays
out of the next run's:

```
target-container/aarch64-unknown-linux-gnu/release/bam tmp/fuzz-corpus/bam fuzz/corpus/bam \
    fuzz/regressions/bam -fork=4 -ignore_ooms=1 -ignore_timeouts=1 -ignore_crashes=1 \
    -max_total_time=2400 -rss_limit_mb=3000 -malloc_limit_mb=3000 -timeout=60 \
    -artifact_prefix=tmp/fuzz-artifacts/bam/
```
