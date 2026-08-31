//! **Does the shipped reader keep the 1.8× the format was sold on?**
//!
//! ```text
//! cargo run --release --example ng_psp_against_production -- \
//!     <a production .psp> <the ng store built from it> [--rounds R] [--only ng|production]
//! ```
//!
//! # The claim under test
//!
//! Spec [`psp_file_format.md`](../doc/devel/ng/spec/psp_file_format.md) §5.4 records the design's
//! headline speed figure: over 62 samples and 471,520,156 records, production's `.psp` took
//! **42.4 s** and this format **23.1 s** — 1.8× faster on the same records.
//!
//! **That number was taken on the measuring prototype** (`examples/psp_row_stream_roundtrip.rs`),
//! not on `src/ng/psp/`, which is the code that shipped. Nothing had ever timed the shipped reader
//! against production's. This does.
//!
//! # What it measures, and what it deliberately does not
//!
//! One full walk of each reader over one sample, every record built. Both arms count the records
//! they saw and both touch every field of every record they build — a walk whose result is thrown
//! away can be optimised into a walk that reads nothing, and the two readers would not be
//! optimised away by the same amount.
//!
//! `--only ng` and `--only production` walk one reader instead of two, for putting a sampling
//! profiler on one of them without the other's frames in the same profile. A run under either
//! reports a time and no ratio.
//!
//! **The two arms are interleaved round by round**, for the reason
//! `examples/ng_psp_skip_value.rs` gives: run as two phases the arms see different machine
//! conditions and the ratio picks up whatever changed between them. The reported figure is the
//! **cheapest round of each arm**, which is the round least disturbed by everything else on the
//! machine.
//!
//! **It is one sample, not 62.** §5.4's walk holds 62 readers open and advances them together; the
//! per-sample memory that shape exists to measure is `examples/ng_psp_open_cost.rs`'s question and
//! was settled by milestone H4. What is compared here is decode throughput per record, which is
//! what the 1.8× is a ratio of.
//!
//! ⚠ **The ng store is built from the production file, so it names the reads production named.**
//! Production names about 3.4 % of the reads ng will name (the owner's ruling of 2026-08-17), and
//! ng's chain-id changes ride in the record head. So this walk under-weights the part of ng's
//! record that will grow most, and the ratio it reports is an **upper bound** on ng's advantage at
//! this depth.

// **The allocator the shipped run uses, and it has to be asked for here.** `alloc-mimalloc` is a
// default feature, but a `#[global_allocator]` is per *binary*: an example that does not declare
// one links the system allocator however the feature is set. Measured on this harness, over the
// tomato accession's 7.69 M records, `--only ng --rounds 30` both times: **1.024 s under the
// system allocator against 0.807 s under mimalloc**, so the system allocator makes this reader
// look 27 % slower than it is. Every psp harness under `examples/` that omits this declaration has
// been timing the wrong allocator, `ng_psp_skip_value.rs` (milestone H5) among them.
#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::Instant;

use pop_var_caller::ng::psp::PspReader as NgPspReader;
use pop_var_caller::psp::PspReader as ProductionPspReader;

/// How many rounds each arm is timed over. **Seven**, for the reason
/// `examples/ng_psp_skip_value.rs` measured: at three rounds five repeats of one command gave an
/// 8 % spread of the ratio, at seven a 3 % one.
const DEFAULT_ROUNDS: usize = 7;

/// What one walk cost, and enough of what it saw to prove it happened.
struct WhatTheWalkCost {
    seconds: f64,
    records: u64,
    /// Summed over every record built, so the compiler cannot delete the walk — and so a run that
    /// silently read half a file is visible as half a total.
    observations: u64,
    /// Summed reads behind those observations: the depth axis, and which corner of the range this
    /// reading belongs to.
    reads: u64,
}

fn main() {
    let mut args = std::env::args().skip(1);
    let production = PathBuf::from(args.next().unwrap_or_else(|| {
        panic!(
            "usage: ng_psp_against_production <a production .psp> <the ng store built from it> \
             [--rounds R] [--only ng|production]"
        )
    }));
    let ng = PathBuf::from(
        args.next()
            .expect("the ng store built from that production file"),
    );
    let mut rounds = DEFAULT_ROUNDS;
    let mut only = Arm::Both;
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--rounds" => {
                let value = args.next().expect("--rounds needs a count");
                rounds = value
                    .parse()
                    .unwrap_or_else(|why| panic!("--rounds {value:?} is not a count: {why}"));
            }
            "--only" => {
                let value = args.next().expect("--only needs ng or production");
                only = match value.as_str() {
                    "ng" => Arm::Ng,
                    "production" => Arm::Production,
                    other => panic!("--only takes ng or production, not {other:?}"),
                };
            }
            other => panic!("unknown flag {other}"),
        }
    }
    assert!(rounds > 0, "--rounds needs at least 1");

    let (best_production, best_ng) = best_of_interleaved(rounds, only, &production, &ng);

    if only == Arm::Both {
        assert_eq!(
            best_production.records, best_ng.records,
            "the two stores hold different record counts, so the walks are not over the same \
             records — rebuild the ng store from this production file with \
             examples/ng_psp_parity.rs"
        );
    }

    println!("phase\tng-psp-against-production");
    println!("production\t{}", production.display());
    println!("ng\t{}", ng.display());
    println!("rounds\t{rounds}");
    println!("records\t{}", best_ng.records);
    println!(
        "reads-a-record\t{:.1}",
        best_ng.reads as f64 / best_ng.records.max(1) as f64
    );
    println!("production-seconds\t{:.4}", best_production.seconds);
    println!("ng-seconds\t{:.4}", best_ng.seconds);
    println!(
        "ng-is-faster-by\t{:.3}",
        best_production.seconds / best_ng.seconds
    );
    println!(
        "production-ns-a-record\t{:.1}",
        best_production.seconds * 1e9 / best_production.records.max(1) as f64
    );
    println!(
        "ng-ns-a-record\t{:.1}",
        best_ng.seconds * 1e9 / best_ng.records.max(1) as f64
    );
    println!("production-observations\t{}", best_production.observations);
    println!("ng-observations\t{}", best_ng.observations);
    println!("production-bytes\t{}", bytes_of(&production));
    println!("ng-bytes\t{}", bytes_of(&ng));
}

/// Which readers a run walks.
///
/// **`Both` is the only setting that yields a ratio**, and it is the default: two arms timed in
/// one process, alternating. The single-arm settings exist for one job — putting a sampling
/// profiler on one reader without the other's frames in the same profile — and a run under one of
/// them reports a time and no comparison.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Arm {
    Both,
    Ng,
    Production,
}

/// Time each requested arm once per round, alternating, and keep the cheapest of each.
///
/// **The cheapest round, not the mean**, and the arms alternate rather than running in two
/// phases: `examples/ng_psp_skip_value.rs` measured a phased version of the same shape giving
/// ratios of 2.93, 4.41 and 3.07 on one corpus under load.
fn best_of_interleaved(
    rounds: usize,
    only: Arm,
    production: &Path,
    ng: &Path,
) -> (WhatTheWalkCost, WhatTheWalkCost) {
    let mut best_production: Option<WhatTheWalkCost> = None;
    let mut best_ng: Option<WhatTheWalkCost> = None;
    for _ in 0..rounds {
        if only != Arm::Ng {
            let one = walk_production(production);
            if best_production
                .as_ref()
                .is_none_or(|best| one.seconds < best.seconds)
            {
                best_production = Some(one);
            }
        }
        if only != Arm::Production {
            let one = walk_ng(ng);
            if best_ng
                .as_ref()
                .is_none_or(|best| one.seconds < best.seconds)
            {
                best_ng = Some(one);
            }
        }
    }
    (
        best_production.unwrap_or_else(nothing_walked),
        best_ng.unwrap_or_else(nothing_walked),
    )
}

/// The arm a run was told not to walk: every total zero, so a reading taken from it is visibly
/// absent rather than plausibly small.
fn nothing_walked() -> WhatTheWalkCost {
    WhatTheWalkCost {
        seconds: 0.0,
        records: 0,
        observations: 0,
        reads: 0,
    }
}

/// Production's reader, over its own `.psp`, with the 1 MB `BufReader` its cohort driver uses.
fn walk_production(path: &Path) -> WhatTheWalkCost {
    let file = File::open(path).unwrap_or_else(|why| panic!("opening {}: {why}", path.display()));
    let mut reader = ProductionPspReader::new(BufReader::with_capacity(1 << 20, file))
        .unwrap_or_else(|why| {
            panic!(
                "{} is not a production .psp this build reads: {why}",
                path.display()
            )
        });
    let started = Instant::now();
    let (mut records, mut observations, mut reads) = (0u64, 0u64, 0u64);
    for record in reader.records() {
        let record = record.expect("a production record");
        records += 1;
        observations += record.alleles.len() as u64;
        for allele in &record.alleles {
            reads += u64::from(allele.support.num_obs);
        }
    }
    WhatTheWalkCost {
        seconds: started.elapsed().as_secs_f64(),
        records,
        observations,
        reads,
    }
}

/// ng's reader, over the ng store, every body built.
fn walk_ng(path: &Path) -> WhatTheWalkCost {
    let mut reader =
        NgPspReader::open(path).unwrap_or_else(|why| panic!("opening {}: {why}", path.display()));
    let started = Instant::now();
    let (mut records, mut observations, mut reads) = (0u64, 0u64, 0u64);
    for record in reader.records().expect("start the ng walk") {
        let record = record.expect("an ng record");
        records += 1;
        let built = record.record.expect("a full walk builds every body");
        observations += built.observations.len() as u64;
        for observation in &built.observations {
            reads += u64::from(observation.num_obs);
        }
    }
    WhatTheWalkCost {
        seconds: started.elapsed().as_secs_f64(),
        records,
        observations,
        reads,
    }
}

fn bytes_of(path: &Path) -> u64 {
    std::fs::metadata(path).map(|it| it.len()).unwrap_or(0)
}
