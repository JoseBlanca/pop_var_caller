//! **Gain/loss of the recommended delimiter vs the previous best, on human data.** Compares
//! unit_robust (algo 4u, recommended) against unit_slip (algo 4, the best we had before) over real
//! HG002 reads, and adjudicates every disagreement with an INDEPENDENT oracle: anchor each of the
//! locus's unique reference flanks in the read by sliding the whole flank (no alignment model), and
//! read off what the read actually shows — did it span both borders, or run off inside the tract?
//!
//! That answers "what did we gain or lose" concretely:
//!  - a read unit_slip called COMPLETE but unit_robust demoted to PARTIAL is a GAIN if the oracle
//!    confirms the far flank is not in the read (unit_slip fabricated an exact length), a LOSS if the
//!    oracle finds both flanks (a real complete downgraded to a bound).
//!  - on reads both call complete, whose measured length matches the oracle's flank-to-flank tract?
//!
//! ```text
//! ng_ssr_gain_loss <reference.fa> <sample.bam|cram> [contig ...]
//! ```

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::alignment::PerQualityEmission;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner;
use pop_var_caller::ng::locus_generation::LocusGenerator;
use pop_var_caller::ng::locus_generation::ReadWitness;
use pop_var_caller::ng::locus_generation::ssr::{
    RepeatDelimiter, SegmentDelimitations, SsrGenerator, SsrGeneratorConfig,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{Bp, ContigId};

const MIN_OVERLAP: usize = 12; // bp of a flank that must lie in the read to call it anchored
const MAX_RATE: f64 = 0.2; // tolerate SNP / sequencing error in the anchored flank

/// Slide the whole left flank so its 3' end lands at read position `t` (allowing the 5' end to hang
/// off the read start); the smallest-mismatch `t` with enough overlap, or None. Using the full flank
/// (not a short seed) makes the unique flank pin the boundary; the repetitive tract cannot match it.
fn place_left(read: &[u8], flank: &[u8]) -> Option<usize> {
    let l = flank.len();
    let mut best: Option<(f64, usize, usize)> = None;
    for t in MIN_OVERLAP..=read.len() {
        let overlap = t.min(l);
        if overlap < MIN_OVERLAP {
            continue;
        }
        let fs = &flank[l - overlap..];
        let rs = &read[t - overlap..t];
        let mism = fs.iter().zip(rs).filter(|(a, b)| a != b).count();
        let rate = mism as f64 / overlap as f64;
        if rate <= MAX_RATE {
            let key = (rate, usize::MAX - overlap, t);
            if best.is_none_or(|b| key < b) {
                best = Some(key);
            }
        }
    }
    best.map(|b| b.2)
}

/// Slide the whole right flank so its 5' end lands at read position `e >= start`.
fn place_right(read: &[u8], flank: &[u8], start: usize) -> Option<usize> {
    let l = flank.len();
    let mut best: Option<(f64, usize, usize)> = None;
    let mut e = start;
    while e + MIN_OVERLAP <= read.len() {
        let overlap = (read.len() - e).min(l);
        if overlap >= MIN_OVERLAP {
            let fs = &flank[..overlap];
            let rs = &read[e..e + overlap];
            let mism = fs.iter().zip(rs).filter(|(a, b)| a != b).count();
            let rate = mism as f64 / overlap as f64;
            if rate <= MAX_RATE {
                let key = (rate, usize::MAX - overlap, e);
                if best.is_none_or(|b| key < b) {
                    best = Some(key);
                }
            }
        }
        e += 1;
    }
    best.map(|b| b.2)
}

/// The oracle's verdict for one read: what its own sequence proves about spanning the tract.
#[derive(PartialEq, Clone, Copy)]
enum Truth {
    Spanned(usize), // both flanks anchored; the tract between them is this many bp (the true length)
    RanOffRight,    // only the left flank is in the read
    RanOffLeft,     // only the right flank is in the read
    Unknown,        // neither flank could be placed — the oracle cannot judge
}

fn oracle(read: &[u8], left_flank: &[u8], right_flank: &[u8]) -> Truth {
    if left_flank.len() < MIN_OVERLAP && right_flank.len() < MIN_OVERLAP {
        return Truth::Unknown;
    }
    let l = place_left(read, left_flank);
    let r = l.and_then(|t| place_right(read, right_flank, t));
    match (l, r) {
        (Some(t), Some(e)) if e >= t => Truth::Spanned(e - t),
        (Some(_), None) => {
            // left placed, right not after it — but maybe the right flank is simply absent; confirm
            // the read has no right flank anywhere.
            if place_right(read, right_flank, 0).is_none() {
                Truth::RanOffRight
            } else {
                Truth::Unknown
            }
        }
        _ => {
            if place_right(read, right_flank, 0).is_some() {
                Truth::RanOffLeft
            } else {
                Truth::Unknown
            }
        }
    }
}

fn is_complete(o: &Option<(ReadWitness, Vec<u8>)>) -> bool {
    matches!(o, Some((ReadWitness::Complete, _)))
}
fn is_partial(o: &Option<(ReadWitness, Vec<u8>)>) -> bool {
    matches!(o, Some((ReadWitness::Observed { .. }, _)))
}
fn measured_len(o: &Option<(ReadWitness, Vec<u8>)>) -> Option<usize> {
    match o {
        Some((ReadWitness::Complete, b)) => Some(b.len()),
        _ => None,
    }
}

fn make_gen<A: RepeatDelimiter>(
    fasta: &Path,
    contigs: &ContigList,
    aligner: A,
    bundle: Bp,
) -> Result<
    SsrGenerator<WindowedRefSeq, impl FnMut() -> WindowedRefSeq, A>,
    Box<dyn std::error::Error>,
> {
    Ok(SsrGenerator::new(
        WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()),
        {
            let fasta = fasta.to_path_buf();
            let contigs = contigs.clone();
            move || WindowedRefSeq::new(fasta.clone(), contigs.clone())
        },
        aligner,
        SsrGeneratorConfig::default(),
        bundle,
    )?)
}

#[derive(Default)]
struct Counts {
    reads: u64,
    // headline class totals
    old_complete: u64,
    new_complete: u64,
    old_partial: u64,
    new_partial: u64,
    // the demotions (old=complete, new=partial), adjudicated by the oracle
    demote_total: u64,
    demote_gain: u64, // oracle: the read ran off — old fabricated a complete, new is right
    demote_loss: u64, // oracle: the read spanned — a real complete was downgraded to a bound
    demote_unknown: u64, // oracle could not judge
    // shared completes: length vs the oracle's flank-to-flank truth
    shared_complete: u64,
    old_len_right: u64,
    new_len_right: u64,
    both_len_right: u64,
    len_adjudicable: u64,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: ng_ssr_gain_loss <reference.fa> <sample.bam|cram> [contig ...]");
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter: Vec<String> = args[3..].to_vec();
    match run(&fasta, &bam, &contig_filter) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn fold(c: &mut Counts, old: &SegmentDelimitations, new: &SegmentDelimitations) {
    for (o, n) in old.reads.iter().zip(new.reads.iter()) {
        c.reads += 1;
        c.old_complete += is_complete(&o.observation) as u64;
        c.new_complete += is_complete(&n.observation) as u64;
        c.old_partial += is_partial(&o.observation) as u64;
        c.new_partial += is_partial(&n.observation) as u64;

        // Demotion adjudication.
        if is_complete(&o.observation) && is_partial(&n.observation) {
            c.demote_total += 1;
            match oracle(&o.read_seq, &old.left_flank, &old.right_flank) {
                Truth::Spanned(_) => c.demote_loss += 1,
                Truth::RanOffLeft | Truth::RanOffRight => c.demote_gain += 1,
                Truth::Unknown => c.demote_unknown += 1,
            }
        }

        // Shared completes: whose exact length matches the oracle's flank-to-flank tract?
        if is_complete(&o.observation) && is_complete(&n.observation) {
            c.shared_complete += 1;
            if let Truth::Spanned(truth) = oracle(&o.read_seq, &old.left_flank, &old.right_flank) {
                c.len_adjudicable += 1;
                let ol = measured_len(&o.observation) == Some(truth);
                let nl = measured_len(&n.observation) == Some(truth);
                c.old_len_right += ol as u64;
                c.new_len_right += nl as u64;
                c.both_len_right += (ol && nl) as u64;
            }
        }
    }
}

fn run(
    fasta: &Path,
    bam: &Path,
    contig_filter: &[String],
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = info.contig_list();
    let sample = SampleReads::open_only_sample(
        std::slice::from_ref(&bam.to_path_buf()),
        &info,
        ReadFilterConfig::default(),
        true,
    )?;
    let walk_config = TypedRegionConfig::default();
    let bundle = Bp(walk_config.criteria.bundle_threshold);
    let e = PerQualityEmission::new();
    let mut old = make_gen(fasta, &contigs, SsrUnitSlipAligner::new(e), bundle)?;
    let mut new = make_gen(fasta, &contigs, SsrUnitRobustAligner::new(e), bundle)?;

    let mut c = Counts::default();
    for (index, entry) in contigs.entries.iter().enumerate() {
        if !contig_filter.is_empty() && !contig_filter.iter().any(|n| n == &entry.name) {
            continue;
        }
        let walk_ref = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
        let mut walk = TypedRegionIterator::over_contig(
            walk_ref,
            ContigId(index as u32),
            walk_config.clone(),
        )?;
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            old.begin_segment(region.region);
            new.begin_segment(region.region);
            let o = old.delimit_segment_reads(segment, &sample)?;
            let n = new.delimit_segment_reads(segment, &sample)?;
            fold(&mut c, &o, &n);
        }
    }
    if let Some(h) = verify {
        h.join()?;
    }

    let pct = |a: u64, b: u64| {
        if b == 0 {
            0.0
        } else {
            100.0 * a as f64 / b as f64
        }
    };
    println!("== unit_robust (recommended) vs unit_slip (previous best), HG002 real reads ==");
    println!("reads compared: {}", c.reads);
    println!("\nobservation classes (per read):");
    println!(
        "  complete:  unit_slip {:>7}   unit_robust {:>7}   ({:+})",
        c.old_complete,
        c.new_complete,
        c.new_complete as i64 - c.old_complete as i64
    );
    println!(
        "  partial :  unit_slip {:>7}   unit_robust {:>7}   ({:+})",
        c.old_partial,
        c.new_partial,
        c.new_partial as i64 - c.old_partial as i64
    );

    println!(
        "\nGAIN/LOSS — the {} completes unit_robust demoted to partial, judged by the flank oracle:",
        c.demote_total
    );
    println!(
        "  GAIN (oracle: read ran off — unit_slip fabricated a complete): {:>6}  ({:.1}%)",
        c.demote_gain,
        pct(c.demote_gain, c.demote_total)
    );
    println!(
        "  LOSS (oracle: read spanned — a real complete downgraded)     : {:>6}  ({:.1}%)",
        c.demote_loss,
        pct(c.demote_loss, c.demote_total)
    );
    println!(
        "  oracle could not judge (far flank not resolvable)            : {:>6}  ({:.1}%)",
        c.demote_unknown,
        pct(c.demote_unknown, c.demote_total)
    );

    println!(
        "\nLENGTH ACCURACY on the {} reads BOTH call complete ({} oracle-adjudicable):",
        c.shared_complete, c.len_adjudicable
    );
    println!(
        "  unit_slip   length == oracle tract: {:>6}  ({:.1}%)",
        c.old_len_right,
        pct(c.old_len_right, c.len_adjudicable)
    );
    println!(
        "  unit_robust length == oracle tract: {:>6}  ({:.1}%)",
        c.new_len_right,
        pct(c.new_len_right, c.len_adjudicable)
    );
    Ok(())
}
