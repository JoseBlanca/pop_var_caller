//! Does the joint fit's locus selection do what its spec says, on a real genome?
//!
//! Four questions, one run over one reference and its repeat catalog:
//!
//! 1. **How much of the reference is not sequence at all** — runs of `N` and other
//!    ambiguity codes. A position inside one is kept by a rule that cannot see it, never
//!    covered by a read, and then enters every per-locus rate as a denominator with no
//!    numerator.
//! 2. **Does the generic rule keep what it was asked for**, and does it scatter — the
//!    realised count against the target, the gap distribution, and the per-contig counts
//!    against each contig's share of the genome.
//! 3. **Is it order-independent on real coordinates** — the whole-reference sweep against
//!    a sweep over the same reference cut into shards that respect nothing.
//! 4. **How many STR loci does each stratum hold at the calling floors** — the count that
//!    decides what the per-stratum cap should be, whether any sampling is needed at all,
//!    and how much of the STR path the joint route would hold
//!    (`parameter_prepass_joint_loci.md` §6 question 1).
//!
//! ```text
//! ng_joint_loci_probe <reference.fa> <catalog.parquet> [generic-target] [str-cap]
//! ```

use std::path::PathBuf;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::joint::loci::{
    SelectableRegions, UnambiguousRuns, select_generic_positions, threshold_for,
};
use pop_var_caller::ng::reference_info::{
    ContigInfo, ReferenceSource, read_reference_info_observing,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegmentCriteria;
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::types::{Bp, ContigId, GenomeRegion, Position};

fn whole_reference_regions(contigs: &[ContigInfo]) -> SelectableRegions {
    SelectableRegions::new(
        contigs
            .iter()
            .enumerate()
            .map(|(i, c)| GenomeRegion {
                contig: ContigId(i as u32),
                start: Position(1),
                end: Position(c.length),
            })
            .collect(),
    )
    .expect("whole contigs never overlap")
}

/// Cut every span into pieces at coordinates the rule knows nothing about, so that a
/// selection that secretly depends on traversal order comes out different.
fn shard(regions: &SelectableRegions) -> SelectableRegions {
    let mut spans = Vec::new();
    for (i, span) in regions.spans().iter().enumerate() {
        let length = span.end.get() - span.start.get() + 1;
        // Three uneven pieces, offset per span so the cuts never line up.
        let end_exclusive = span.end.get() + 1;
        let clamp = |c: u64| c.clamp(span.start.get(), end_exclusive);
        let first = clamp(span.start.get() + (length / 3) + (i as u64 % 7));
        let second = clamp(span.start.get() + (2 * length / 3) + (i as u64 % 11));
        let mut cuts = vec![span.start.get(), first, second, end_exclusive];
        cuts.sort_unstable();
        cuts.dedup();
        for pair in cuts.windows(2) {
            if pair[0] < pair[1] && pair[1] <= span.end.get() + 1 {
                spans.push(GenomeRegion {
                    contig: span.contig,
                    start: Position(pair[0]),
                    end: Position(pair[1] - 1),
                });
            }
        }
    }
    SelectableRegions::new(spans).expect("cuts of disjoint spans stay disjoint")
}

fn main() {
    let mut args = std::env::args().skip(1);
    let reference = PathBuf::from(
        args.next()
            .expect("usage: <reference.fa> <catalog.parquet> [generic-target] [str-cap]"),
    );
    let catalog_path = PathBuf::from(args.next().expect("a catalog"));
    let generic_target: u64 = args
        .next()
        .map_or(2_000_000, |a| a.parse().expect("a position target"));
    let str_cap: usize = args
        .next()
        .map_or(1_000_000, |a| a.parse().expect("a per-stratum cap"));
    let seed = 42_u64;

    // ---- the reference, and where it is sequence -------------------------------------
    let started = Instant::now();
    let mut callable = UnambiguousRuns::default();
    let info = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: reference.clone(),
            fai: None,
        },
        &mut callable,
    )
    .expect("the reference reads");
    let names: Vec<String> = info.contigs.iter().map(|c| c.name.clone()).collect();
    let whole = whole_reference_regions(&info.contigs);
    let ambiguous_bases = callable.ambiguous_bases();
    let callable_spans = callable
        .into_selectable()
        .expect("maximal runs are disjoint");

    println!("reference          {}", reference.display());
    println!("contigs            {}", info.contigs.len());
    println!(
        "reference bases    {} ({:.1} Mb), read in {:.1} s",
        whole.total_length().get(),
        whole.total_length().get() as f64 / 1e6,
        started.elapsed().as_secs_f64()
    );
    println!(
        "  A/C/G/T          {} ({:.2}%)",
        callable_spans.total_length().get(),
        100.0 * callable_spans.total_length().get() as f64 / whole.total_length().get() as f64
    );
    println!(
        "  N or ambiguous   {} ({:.2}%), leaving {} unambiguous runs",
        ambiguous_bases,
        100.0 * ambiguous_bases as f64 / whole.total_length().get() as f64,
        callable_spans.spans().len()
    );

    // ---- the generic rule ------------------------------------------------------------
    println!("\n--- generic selection, target {generic_target}, seed {seed} ---");
    let over_whole = {
        let threshold = threshold_for(generic_target, whole.total_length());
        let started = Instant::now();
        let kept = select_generic_positions(&whole, &names, seed, threshold).expect("selects");
        println!(
            "over every reference base : kept {} ({:+.2}% of target) in {:.1} s",
            kept.len(),
            100.0 * (kept.len() as f64 - generic_target as f64) / generic_target as f64,
            started.elapsed().as_secs_f64()
        );
        kept
    };
    let over_callable = {
        let threshold = threshold_for(generic_target, callable_spans.total_length());
        let kept =
            select_generic_positions(&callable_spans, &names, seed, threshold).expect("selects");
        println!(
            "over A/C/G/T only         : kept {} ({:+.2}% of target)",
            kept.len(),
            100.0 * (kept.len() as f64 - generic_target as f64) / generic_target as f64,
        );
        kept
    };

    // How many of the naive selection's positions sit where no read can ever align.
    let mut in_ambiguous = 0_u64;
    {
        let mut spans = callable_spans.spans().iter().peekable();
        let mut current: Option<GenomeRegion> = spans.next().copied();
        for position in &over_whole {
            while let Some(span) = current {
                if span.contig.get() < position.contig.get()
                    || (span.contig == position.contig && span.end.get() < position.position.get())
                {
                    current = spans.next().copied();
                } else {
                    break;
                }
            }
            let covered = current.is_some_and(|span| {
                span.contig == position.contig
                    && span.start.get() <= position.position.get()
                    && position.position.get() <= span.end.get()
            });
            if !covered {
                in_ambiguous += 1;
            }
        }
    }
    println!(
        "  of the first selection, {} positions ({:.2}%) sit in N — kept, never covered",
        in_ambiguous,
        100.0 * in_ambiguous as f64 / over_whole.len() as f64
    );

    // Scatter: how close the kept positions sit to one another.
    let mut gaps: Vec<u64> = Vec::with_capacity(over_callable.len());
    for pair in over_callable.windows(2) {
        if pair[0].contig == pair[1].contig {
            gaps.push(pair[1].position.get() - pair[0].position.get());
        }
    }
    gaps.sort_unstable();
    let quantile = |q: f64| gaps[((gaps.len() as f64 * q) as usize).min(gaps.len() - 1)];
    let within = |d: u64| gaps.iter().filter(|&&g| g <= d).count();
    println!(
        "  gaps between kept positions: p10 {}, median {}, p90 {} bp",
        quantile(0.10),
        quantile(0.50),
        quantile(0.90)
    );
    println!(
        "  neighbours within 100 bp {} ({:.1}%), within 1 kb {} ({:.1}%)",
        within(100),
        100.0 * within(100) as f64 / gaps.len() as f64,
        within(1_000),
        100.0 * within(1_000) as f64 / gaps.len() as f64
    );

    // Uniformity across contigs. The expectations are scaled to the **realised** total,
    // not to the target: the realised count is itself binomial around the target, and
    // charging that shared departure to every contig would report an unbiased hash as
    // skewed.
    let mut per_contig = vec![0_u64; info.contigs.len()];
    for position in &over_callable {
        per_contig[position.contig.get() as usize] += 1;
    }
    let mut selectable_per_contig = vec![0_u64; info.contigs.len()];
    for span in callable_spans.spans() {
        selectable_per_contig[span.contig.get() as usize] += span.end.get() - span.start.get() + 1;
    }
    let mut worst = (0.0_f64, String::new(), 0_u64, 0.0_f64);
    let mut chi_square = 0.0_f64;
    let mut degrees = 0_usize;
    for (i, contig) in info.contigs.iter().enumerate() {
        if selectable_per_contig[i] < 10_000 {
            continue;
        }
        let expected = over_callable.len() as f64 * selectable_per_contig[i] as f64
            / callable_spans.total_length().get() as f64;
        let z = (per_contig[i] as f64 - expected) / expected.sqrt();
        chi_square += z * z;
        degrees += 1;
        if z.abs() > worst.0.abs() {
            worst = (z, contig.name.clone(), per_contig[i], expected);
        }
    }
    let degrees = degrees.saturating_sub(1);
    println!(
        "  per-contig uniformity: chi-square {:.1} on {} degrees of freedom ({:+.1} sd of its own \
         null), worst contig {} at {:+.1} sd",
        chi_square,
        degrees,
        (chi_square - degrees as f64) / (2.0 * degrees as f64).sqrt(),
        worst.1,
        worst.0
    );

    // Order independence, on real coordinates rather than a fixture's.
    let sharded = shard(&callable_spans);
    let threshold = threshold_for(generic_target, callable_spans.total_length());
    let from_shards = select_generic_positions(&sharded, &names, seed, threshold).expect("selects");
    println!(
        "  sharded sweep ({} spans against {}): {}",
        sharded.spans().len(),
        callable_spans.spans().len(),
        if from_shards == over_callable {
            "identical"
        } else {
            "DIFFERENT — the rule is not order-independent"
        }
    );

    // ---- the STR strata --------------------------------------------------------------
    // A stale catalog is not a reason to lose the generic half of the run, which needs no
    // catalog at all.
    let catalog = match RepeatCatalog::open_checking_against_reference(&catalog_path, &info) {
        Ok(catalog) => catalog,
        Err(e) => {
            println!("\n--- STR strata: skipped, the catalog does not open ---\n{e}");
            return;
        }
    };
    let calling = StrRepeatCriteria {
        classification: SsrSegmentCriteria::default(),
        min_flank_bp: Bp(30),
        max_str_len_bp: Bp(100),
    };
    println!(
        "\n--- STR strata at the calling floors {:?}, flank {}, satellite cap {} ---",
        (1..=6)
            .map(|p| calling.classification.min_copies.for_period(p))
            .collect::<Vec<_>>(),
        calling.min_flank_bp.get(),
        calling.max_str_len_bp.get()
    );
    let started = Instant::now();
    let (counts, sample) = catalog
        .sample_loci_per_stratum(&calling, ReadScope::WholeReference, str_cap, seed)
        .expect("the catalog serves the calling floors");
    println!(
        "read in {:.1} s; {} strata hold {} loci",
        started.elapsed().as_secs_f64(),
        counts.iter_sorted().len(),
        counts.total()
    );

    println!("\nperiod  repeats      loci     kept");
    for ((period, repeats), held) in counts.iter_sorted() {
        println!(
            "{period:>6}  {repeats:>7}  {held:>8}  {:>7}",
            sample.get(period, repeats).len()
        );
    }

    // What a cap would cost, and whether one is needed at all.
    println!("\ncap    loci kept   strata capped");
    for cap in [100_u64, 500, 1_000, 5_000, 20_000, 100_000] {
        let kept: u64 = counts
            .iter_sorted()
            .iter()
            .map(|(_, n)| (*n).min(cap))
            .sum();
        let capped = counts
            .iter_sorted()
            .iter()
            .filter(|(_, n)| *n > cap)
            .count();
        println!("{cap:>5}  {kept:>10}  {capped:>13}");
    }
    let largest = counts
        .iter_sorted()
        .into_iter()
        .max_by_key(|(_, n)| *n)
        .expect("at least one stratum");
    println!(
        "largest stratum: period {} at {} repeats, {} loci",
        largest.0.0, largest.0.1, largest.1
    );
}
