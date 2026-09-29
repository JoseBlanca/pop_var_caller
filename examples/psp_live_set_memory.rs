//! **How many read identifiers does a sample's psp source hold, and what does decoding cost?**
//!
//! ```text
//! cargo run --release --example psp_live_set_memory -- [--window-bp N] [--stop-after-bp N] <psp>...
//! ```
//!
//! A cohort run holds each sample's records over about one building round of the genome, and
//! for every record it holds it keeps what that record's body needs to be decoded: the reads
//! live at its position (`doc/devel/reports/psp_pileup_memory_2026-09-28/report.md`). This walks
//! each file the same way — drawing every record, releasing whatever starts more than
//! `--window-bp` bases (default 16,000, the tomato run's round) behind the newest — and prints
//! the most identifiers the source held at once, with the record count it held at that moment.
//!
//! **Then it decodes every record once**, in a second walk that holds the same window and builds
//! each body just before releasing it, and prints the time. That is more building than a run
//! does — a run builds about one record in eight — so it bounds the decoding cost from above.
//!
//! It depends only on `PspSummarySource`'s public surface, so it runs unchanged on either way of
//! storing the live sets, which is the comparison it exists for.
//!
//! `--stop-after-bp N` ends each walk at the first record past base `N` of the file's first
//! contig, or on any later contig — enough to reach the ch00 pile-up at 1.5 Mb in a tomato file
//! without walking the rest of the genome.

use std::collections::VecDeque;
use std::path::{Path, PathBuf};
use std::time::Instant;

use pop_var_caller::psp::PspReader;
use pop_var_caller::run::cohort_merge::observation_cache::{Drawn, ObservationSource};
use pop_var_caller::run::psp_source::PspSummarySource;
use pop_var_caller::types::ReadGroupId;

fn main() {
    let mut window_bp: u64 = 16_000;
    let mut stop_after_bp = u64::MAX;
    let mut paths = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--window-bp" {
            window_bp = args.next().expect("a width").parse().expect("a number");
        } else if arg == "--stop-after-bp" {
            stop_after_bp = args.next().expect("a base").parse().expect("a number");
        } else {
            paths.push(PathBuf::from(arg));
        }
    }
    println!(
        "sample\trecords\tpeak_ids\trecords_at_peak\tids_per_record_at_peak\twalk_s\tbuild_all_s"
    );
    for path in &paths {
        let (records, peak_ids, records_at_peak, walk_s) =
            walk(path, window_bp, stop_after_bp, false);
        let (_, _, _, build_s) = walk(path, window_bp, stop_after_bp, true);
        let sample = path.file_stem().unwrap().to_string_lossy();
        println!(
            "{sample}\t{records}\t{peak_ids}\t{records_at_peak}\t{:.1}\t{walk_s:.2}\t{build_s:.2}",
            peak_ids as f64 / records_at_peak.max(1) as f64
        );
    }
}

/// Walk one file holding a sliding window; build each record before it is released if `build`.
fn walk(path: &Path, window_bp: u64, stop_after_bp: u64, build: bool) -> (u64, usize, usize, f64) {
    let mut reader = PspReader::open(path).expect("a psp that opens");
    let groups: Vec<ReadGroupId> = (0..reader.header().read_groups.len() as u32)
        .map(ReadGroupId)
        .collect();
    let mut source = PspSummarySource::over(&mut reader, &groups).expect("a walk");
    let mut held: VecDeque<(u32, u64, core::ops::Range<usize>)> = VecDeque::new();
    let (mut records, mut peak_ids, mut records_at_peak) = (0u64, 0usize, 0usize);
    let mut first_contig: Option<u32> = None;
    let started = Instant::now();
    loop {
        let drawn = match source.next_drawn(None) {
            Some(drawn) => drawn.expect("a record that reads"),
            None => break,
        };
        let Drawn::Kept { summary, body } = drawn else {
            unreachable!("a summary source keeps every body")
        };
        records += 1;
        let contig = summary.region().contig.get();
        let start = summary.region().start.get();
        let on_a_later_contig = contig != *first_contig.get_or_insert(contig);
        if stop_after_bp != u64::MAX && (on_a_later_contig || start > stop_after_bp) {
            break;
        }
        held.push_back((contig, start, body));
        // Release what the window has passed: another contig, or too far behind.
        let mut release_to = None;
        while let Some((c, s, _)) = held.front() {
            if *c == contig && s + window_bp > start {
                break;
            }
            let (_, _, passed) = held.pop_front().unwrap();
            if build {
                std::hint::black_box(source.build(passed.clone()).expect("a body that decodes"));
            }
            release_to = Some(passed.end);
        }
        if let Some(to) = release_to {
            source.release_before(to);
        }
        if source.held_live_ids() > peak_ids {
            peak_ids = source.held_live_ids();
            records_at_peak = source.held_records();
        }
    }
    if build {
        for (_, _, body) in held.drain(..) {
            std::hint::black_box(source.build(body).expect("a body that decodes"));
        }
    }
    (
        records,
        peak_ids,
        records_at_peak,
        started.elapsed().as_secs_f64(),
    )
}
