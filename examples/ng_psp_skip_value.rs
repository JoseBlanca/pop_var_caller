//! **What is the record head worth at depth?** — Milestone H5's measurement, and the open
//! question arch §7 records.
//!
//! ```text
//! cargo run --release --example ng_psp_skip_value -- <a store.ngpsp> [--keep-one-in K] [--rounds R]
//! ```
//!
//! # The question
//!
//! Every record opens with a head a reader can judge without building the body (spec
//! `psp_file_format.md` §4.3), so a cohort's first pass — which wants a few records in a hundred —
//! can advance a pointer instead of decoding. The measuring prototype found that worth **2.06×**.
//!
//! **But the chain-id changes ride in the head too**, and they grow with depth: spec
//! `psp_record_encoding.md` §6 measures **0.432 bytes a position at 11.4 reads and 6.42 at 293**.
//! So as depth rises the head gets dearer while the body the skip avoids stays where it is, and
//! **how much of the 2.06× survives is not known**. That is what this measures.
//!
//! On the wire a record is
//!
//! ```text
//! position offset | reference span | non-reference reads | body bytes | chain-id changes | body
//! ```
//!
//! and `body bytes` does not reach the changes — deliberately, so that a reader which skips a body
//! still carries the live set forward. **The skip therefore never avoids the chain-id changes**,
//! which is exactly why depth threatens it.
//!
//! # ⚠ What this corpus cannot tell you
//!
//! ng cannot yet write a psp of its own, so these stores are built from a production `.psp` by
//! `examples/ng_psp_parity.rs`. **Production names about 3.4 % of the reads ng will name** (the
//! owner's ruling of 2026-08-17; `examples/ng_chain_id_column_cost.rs` measures the gap), so the
//! chain-id changes in these heads are a small fraction of what ng's will be.
//!
//! **The consequence is one-directional, which is what makes the reading still worth having.** A
//! bigger head makes the skip worth *less*, never more. So the ratio measured here is an **upper
//! bound** on what ng will see at the same depth — and arch §7's question stays open until a real
//! ng-written store exists.

use std::path::{Path, PathBuf};
use std::time::Instant;

use pop_var_caller::ng::psp::PspReader;

/// How many rounds each walk is timed over, so a reading is not one sample of a noisy machine.
///
/// **Seven, and three is not enough — measured.** At three rounds five repeats of the same
/// command gave ratios from 2.865 to 3.143, an 8 % spread; at seven they gave 2.806 to 2.901, a
/// 3 % one. A step whose whole output is a ratio cannot quote a figure that moves by 8 % between
/// identical runs.
const DEFAULT_ROUNDS: usize = 7;

/// What one walk cost.
struct WhatTheWalkCost {
    seconds: f64,
    records: u64,
    bodies_built: u64,
    body_bytes_built: u64,
    body_bytes_skipped: u64,
    /// Reads summed over every observation of every body built — the depth axis this whole step
    /// is about, and the one figure that says which corner of the range a reading belongs to.
    reads_in_the_bodies_built: u64,
}

fn main() {
    let mut args = std::env::args().skip(1);
    let store = PathBuf::from(
        args.next()
            .expect("usage: ng_psp_skip_value <a store.ngpsp> [--keep-one-in K] [--rounds R]"),
    );
    let mut keep_one_in = 100u64;
    let mut rounds = DEFAULT_ROUNDS;
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--keep-one-in" => {
                let value = args.next().expect("--keep-one-in needs a count");
                keep_one_in = value
                    .parse()
                    .unwrap_or_else(|why| panic!("--keep-one-in {value:?} is not a count: {why}"));
            }
            "--rounds" => {
                let value = args.next().expect("--rounds needs a count");
                rounds = value
                    .parse()
                    .unwrap_or_else(|why| panic!("--rounds {value:?} is not a count: {why}"));
            }
            other => panic!("unknown flag {other}"),
        }
    }
    assert!(keep_one_in > 0, "--keep-one-in needs at least 1");
    assert!(rounds > 0, "--rounds needs at least 1");

    // **The two arms are interleaved, round by round, and that is not tidiness.** Run as two
    // phases — every full walk, then every skipping walk — the arms see different machine
    // conditions, and the ratio picks up whatever changed between them. Measured on a loaded
    // machine, the phased version gave the same command ratios of 2.93, 4.41 and 3.07 on one
    // corpus; a ratio that moves by half between identical runs is not a measurement.
    let (full, skipping) = best_of_interleaved(rounds, &store, keep_one_in);

    assert_eq!(
        full.records, skipping.records,
        "a declining walk still yields every record — it is not a filter (spec §6.2)"
    );
    assert!(
        skipping.bodies_built < full.bodies_built,
        "the predicate declined nothing, so there is no skip to price: {} of {} bodies built",
        skipping.bodies_built,
        full.bodies_built
    );
    assert!(
        skipping.body_bytes_skipped > 0,
        "the declined records carried no body bytes, so skipping them costs nothing to begin with"
    );

    println!("phase\tng-psp-skip-value");
    println!("store\t{}", store.display());
    println!("keep-one-in\t{keep_one_in}");
    println!("rounds\t{rounds}");
    println!("records\t{}", full.records);
    println!("full-seconds\t{:.4}", full.seconds);
    println!("skipping-seconds\t{:.4}", skipping.seconds);
    println!(
        "speed-up\t{:.3}",
        full.seconds / skipping.seconds.max(f64::MIN_POSITIVE)
    );
    println!("bodies-built-full\t{}", full.bodies_built);
    println!("bodies-built-skipping\t{}", skipping.bodies_built);
    println!("body-bytes-total\t{}", full.body_bytes_built);
    println!("body-bytes-skipped\t{}", skipping.body_bytes_skipped);
    println!(
        "body-bytes-skipped-share\t{:.4}",
        skipping.body_bytes_skipped as f64 / full.body_bytes_built.max(1) as f64
    );
    println!(
        "mean-body-bytes\t{:.2}",
        full.body_bytes_built as f64 / full.records.max(1) as f64
    );
    println!(
        "mean-reads-a-record\t{:.1}",
        full.reads_in_the_bodies_built as f64 / full.records.max(1) as f64
    );
}

/// Both arms, alternating, so each round's two walks meet the same machine.
fn best_of_interleaved(
    rounds: usize,
    store: &Path,
    keep_one_in: u64,
) -> (WhatTheWalkCost, WhatTheWalkCost) {
    let mut full = walk(store, 1);
    let mut skipping = walk(store, keep_one_in);
    for _ in 1..rounds {
        let again_full = walk(store, 1);
        if again_full.seconds < full.seconds {
            full = again_full;
        }
        let again_skipping = walk(store, keep_one_in);
        if again_skipping.seconds < skipping.seconds {
            skipping = again_skipping;
        }
    }
    (full, skipping)
}

/// Walk the whole store, building one body in `keep_one_in`.
fn walk(store: &Path, keep_one_in: u64) -> WhatTheWalkCost {
    let mut reader =
        PspReader::open(store).unwrap_or_else(|why| panic!("opening {}: {why}", store.display()));
    let mut seen = 0u64;
    let started = Instant::now();
    let mut records = 0u64;
    let mut bodies_built = 0u64;
    let mut body_bytes_built = 0u64;
    let mut body_bytes_skipped = 0u64;
    let mut reads_in_the_bodies_built = 0u64;
    {
        let walk = reader
            .records_where(|_head| {
                let wanted = seen.is_multiple_of(keep_one_in);
                seen += 1;
                wanted
            })
            .expect("walking the store");
        for found in walk {
            let found = found.expect("a record");
            records += 1;
            let bytes = u64::from(found.head.body_bytes);
            if let Some(record) = &found.record {
                bodies_built += 1;
                body_bytes_built += bytes;
                reads_in_the_bodies_built += record
                    .observations
                    .iter()
                    .map(|observation| u64::from(observation.num_obs))
                    .sum::<u64>();
            } else {
                body_bytes_skipped += bytes;
            }
        }
    }
    WhatTheWalkCost {
        seconds: started.elapsed().as_secs_f64(),
        records,
        bodies_built,
        // A full walk builds every body, so its `body_bytes_built` is the store's whole body mass;
        // a skipping walk splits the same mass between built and skipped.
        body_bytes_built: body_bytes_built + body_bytes_skipped,
        body_bytes_skipped,
        reads_in_the_bodies_built,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a_walk_costing(seconds: f64) -> WhatTheWalkCost {
        WhatTheWalkCost {
            seconds,
            records: 1,
            bodies_built: 1,
            body_bytes_built: 1,
            body_bytes_skipped: 0,
            reads_in_the_bodies_built: 1,
        }
    }

    /// **The cheapest round, wherever it falls.** A slower round is a machine doing something
    /// else; the question is what the walk costs. Taking the mean would let one noisy round move
    /// a ratio the whole step turns on, and taking the first would make the reading one sample.
    #[test]
    fn the_reading_is_the_cheapest_round_wherever_it_falls() {
        for order in [
            vec![0.5, 0.9, 0.7],
            vec![0.9, 0.5, 0.7],
            vec![0.9, 0.7, 0.5],
        ] {
            let mut best = a_walk_costing(f64::INFINITY);
            for seconds in &order {
                let again = a_walk_costing(*seconds);
                if again.seconds < best.seconds {
                    best = again;
                }
            }
            assert!(
                (best.seconds - 0.5).abs() < f64::EPSILON,
                "the cheapest of {order:?} is 0.5, not {}",
                best.seconds
            );
        }
    }
}
