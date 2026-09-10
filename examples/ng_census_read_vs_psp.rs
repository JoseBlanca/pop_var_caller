//! **What does a fit save by reading a psp's census instead of rebuilding it from that psp's
//! records?**
//!
//! A census is a cache: everything in it can be recomputed from the psp it lives in
//! ([`parameter_prepass_joint_records.md`](../doc/devel/ng/spec/parameter_prepass_joint_records.md)
//! §6.1). That document argues the cache is worth keeping — *"rebuilding is a full pass, not a
//! seek"* — but never priced it. This does.
//!
//! **The cache lives in the psp's trailer** since `psp_census_pair.md` §3, and this harness does
//! not read it from there. It rebuilds the census with the producer that ships (`census_from_psp`)
//! and writes it to a scratch file of its own, **because `--keep-one-in K` means the census this
//! run has to time is usually not the one in the trailer**: at any K but 1 the budget is this
//! run's, and no file on disk holds that census. The two are the same bytes at K = 1, which is
//! what the walk-versus-rebuild oracle guarantees.
//!
//! ```text
//! cargo run --release --example ng_census_read_vs_psp -- \
//!     <reference.fa> <catalog.parquet> <regions.bed> <a.psp> [--keep-one-in K] [--rounds R]
//! ```
//!
//! # The four walks it times
//!
//! All four end with the same object a fit consumes, or with a stated part of it.
//!
//! - **the census** — decode it **whole**, which is what a fit that holds a cohort in memory
//!   pays per sample. A fit over a cohort too large to hold reads it lazily instead, a section at
//!   a time, and pays less; this arm is the upper end of the census route, not the shipped one.
//! - **psp, every body** — `census_from_psp`, the producer that ships. It decodes every record's
//!   body, because the read-error calibration totals are summed over *every* generic locus and
//!   the per-read-group base-quality sums they need live in the body.
//! - **psp, only the bodies a census needs** — the same walk, declining every body whose head
//!   says the record is at no kept census position, overlaps no kept repeat tract, and carries no
//!   read disagreeing with the reference. This
//!   is the skip the format's record head exists for, and the ceiling on what a psp route could
//!   cost if the calibration totals came from somewhere else.
//! - **psp, no bodies at all** — the floor: decompress every block, parse every head, build
//!   nothing. What no psp route can go below.
//!
//! # The knob that decides the answer
//!
//! `--keep-one-in K` sets the census budget to one position in `K` of the analysed ground. **The
//! whole-genome run keeps about 1 in 400** — two million positions over 800 Mb of tomato — while a
//! run over a few megabases keeps nearly everything, because the budget is a count and not a
//! rate. The census grows with `1/K`; the psp walk does not shrink with it. So a ratio
//! quoted without its `K` says nothing.

use std::error::Error;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::joint::census::{
    CensusWriter, NamedReadGroup, SampleCensusEvidence,
};
use pop_var_caller::ng::parameter_estimation::joint::census_file::{read_census, write_census};
use pop_var_caller::ng::parameter_estimation::joint::loci::UnambiguousRuns;
use pop_var_caller::ng::psp::{PspReader, RecordHead};
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, read_reference_observing_or_creating_fai,
};
use pop_var_caller::ng::region_typing::DEFAULT_MAX_STR_LEN;
use pop_var_caller::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use pop_var_caller::ng::repeat_catalog::RepeatCatalog;
use pop_var_caller::ng::run::{CensusPlan, CensusSelection, Segmentation, census_from_psp};
use pop_var_caller::ng::types::{GenomePosition, GenomeRegion, ReadGroupId};
use pop_var_caller::pop_var_caller_exp::run_ground::{self, GroundRequest, RepeatRouting};

/// How many rounds each walk is timed over.
///
/// **Three, and each walk is run in its own round rather than all rounds of one walk together**,
/// so a machine that gets slower part-way through hurts all four equally.
const DEFAULT_ROUNDS: usize = 3;

/// What one timed walk cost.
struct WhatItCost {
    seconds: Vec<f64>,
    bodies_built: u64,
    records: u64,
    bytes: u64,
}

impl WhatItCost {
    fn new() -> Self {
        Self {
            seconds: Vec::new(),
            bodies_built: 0,
            records: 0,
            bytes: 0,
        }
    }

    /// The fastest round, which is the one least disturbed by whatever else the machine was doing.
    fn best(&self) -> f64 {
        self.seconds.iter().copied().fold(f64::INFINITY, f64::min)
    }
}

fn main() -> ExitCode {
    let mut positional: Vec<String> = Vec::new();
    let mut keep_one_in: u64 = 1;
    let mut rounds = DEFAULT_ROUNDS;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--keep-one-in" => match args.next().and_then(|v| v.parse().ok()) {
                Some(k) if k >= 1 => keep_one_in = k,
                _ => {
                    eprintln!("--keep-one-in wants a count of at least 1");
                    return ExitCode::from(2);
                }
            },
            "--rounds" => match args.next().and_then(|v| v.parse().ok()) {
                Some(r) if r >= 1 => rounds = r,
                _ => {
                    eprintln!("--rounds wants a count of at least 1");
                    return ExitCode::from(2);
                }
            },
            other => positional.push(other.to_string()),
        }
    }
    let [fasta, catalog, bed, psp] = positional.as_slice() else {
        eprintln!(
            "usage: ng_census_read_vs_psp <reference.fa> <catalog.parquet> <regions.bed> \
             <a.psp> [--keep-one-in K] [--rounds R]"
        );
        return ExitCode::from(2);
    };
    match run(
        Path::new(fasta),
        Path::new(catalog),
        Path::new(bed),
        Path::new(psp),
        keep_one_in,
        rounds,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            let mut cause = error.source();
            while let Some(next) = cause {
                eprintln!("  caused by: {next}");
                cause = next.source();
            }
            ExitCode::FAILURE
        }
    }
}

fn run(
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    psp_path: &Path,
    keep_one_in: u64,
    rounds: usize,
) -> Result<(), Box<dyn Error>> {
    let setting_up = Instant::now();

    let mut callable = UnambiguousRuns::default();
    let with_checksums = Arc::new(read_reference_observing_or_creating_fai(
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
        &mut callable,
    )?);
    let unambiguous = callable.into_selectable()?;
    let contigs = with_checksums.contig_list();
    let _reference = OpenReference::new(Arc::clone(&with_checksums));

    let ground = GroundRequest {
        reference: fasta,
        catalog: Some(catalog_path),
        regions: Some(bed),
        routing: RepeatRouting {
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        },
    };
    let analysed = run_ground::analysed_regions(&ground, &contigs)?;
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();
    let segmentation = Arc::new(run_ground::segments_over(
        &ground,
        &analysed,
        &with_checksums,
    )?);

    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &with_checksums)?;
    let selection = CensusSelection {
        seed: CensusSelection::SHIPPED.seed,
        generic_target: (analysed_bases / keep_one_in).max(1),
        ssr_cap: CensusSelection::SHIPPED.ssr_cap,
    };
    let plan = CensusPlan::of_run(
        selection,
        &catalog,
        &analysed,
        &unambiguous,
        &with_checksums,
        &segmentation.inputs().repeat_tract_criteria,
    )?;
    drop(catalog);

    // The census this comparison reads, built from the psp by the producer that ships, so the
    // two routes are two ways to the same bytes and not two different objects. **It goes to a
    // scratch file of its own** because at any `--keep-one-in` but 1 it is not the census in the
    // psp's trailer — the budget is this run's — and because a harness must leave every psp it is
    // given exactly as it found it.
    let work = PathBuf::from(
        std::env::var("NG_WORK").unwrap_or_else(|_| "tmp/ng_census_read_vs_psp".to_string()),
    );
    std::fs::create_dir_all(&work)?;
    let built = census_from_psp(psp_path, &plan, &segmentation)?;
    let census_path = work.join(format!("{}.k{keep_one_in}.census", built.sample));
    let truth = {
        let mut bytes = Vec::new();
        write_census(&built.evidence, &mut bytes)?;
        std::fs::write(&census_path, &bytes)?;
        bytes
    };
    let psp_bytes = std::fs::metadata(psp_path)?.len();
    let census_bytes = truth.len() as u64;
    drop(built);

    println!("# psp: {}", psp_path.display());
    println!("# analysed bases: {analysed_bases}");
    println!("# keep one position in: {keep_one_in}");
    println!(
        "# census positions kept: {} of {analysed_bases}",
        plan.loci.generic().len()
    );
    println!("# census tracts kept: {}", plan.loci.ssr().total());
    println!("# psp bytes: {psp_bytes}");
    println!("# census bytes: {census_bytes}");
    println!("# setup: {:.2} s", setting_up.elapsed().as_secs_f64());

    let mut census = WhatItCost::new();
    let mut every_body = WhatItCost::new();
    let mut needed_bodies = WhatItCost::new();
    let mut no_bodies = WhatItCost::new();
    let mut reads_a_position = 0.0_f64;

    // The kept tracts as regions, so a head can be judged against them the way the census
    // writer judges a locus — replicated from `CensusWriter::new`, which does not expose it.
    let kept_tracts = kept_tract_regions(&plan, &with_checksums);

    for _ in 0..rounds {
        let at = Instant::now();
        let file = std::fs::File::open(&census_path)?;
        let mut file = std::io::BufReader::new(file);
        let mut read = read_census(&mut file)?;
        // **Every decoded value is touched**, so this times a census a fit can use rather than a
        // decode whose arrays are never read. The sum is printed for the same reason.
        let touched = touch_every_value(&mut read.census)?;
        census.seconds.push(at.elapsed().as_secs_f64());
        census.records = touched;
        drop(read);

        let at = Instant::now();
        let produced = census_from_psp(psp_path, &plan, &segmentation)?;
        every_body.seconds.push(at.elapsed().as_secs_f64());
        every_body.records = produced.records;
        every_body.bodies_built = produced.records;
        drop(produced);

        let at = Instant::now();
        let selective = selective_census(psp_path, &plan, &segmentation, |head| {
            head.non_reference_reads > 0
                || holds_a_kept_position(plan.loci.generic(), head.region)
                || kept_tracts
                    .binary_search(&(
                        head.region.contig.get(),
                        head.region.start.get(),
                        head.region.end.get(),
                    ))
                    .is_ok()
        })?;
        needed_bodies.seconds.push(at.elapsed().as_secs_f64());
        needed_bodies.records = selective.records;
        needed_bodies.bodies_built = selective.bodies_built;
        needed_bodies.bytes = {
            let mut bytes = Vec::new();
            write_census(&selective.evidence, &mut bytes)?;
            bytes.len() as u64
        };

        let at = Instant::now();
        let (counted, depth) = heads_only(psp_path)?;
        no_bodies.seconds.push(at.elapsed().as_secs_f64());
        no_bodies.records = counted;
        no_bodies.bytes = depth.round() as u64;
        reads_a_position = depth;
    }

    let report = |name: &str, cost: &WhatItCost| {
        let each: Vec<String> = cost
            .seconds
            .iter()
            .map(|second| format!("{second:.4}"))
            .collect();
        println!(
            "{name}: best {:.4} s  (rounds {})  records={} bodies_built={}",
            cost.best(),
            each.join(" "),
            cost.records,
            cost.bodies_built,
        );
    };
    report("census             ", &census);
    println!("# values read out of the census: {}", census.records);
    report("psp-every-body     ", &every_body);
    report("psp-needed-bodies  ", &needed_bodies);
    report("psp-no-bodies      ", &no_bodies);
    println!("# reads compared with the reference, a position: {reads_a_position:.1}");
    println!(
        "# psp bytes a base of analysed ground: {:.2}",
        psp_bytes as f64 / analysed_bases as f64
    );
    println!(
        "# seconds a megabase of analysed ground: every-body {:.4}, needed-bodies {:.4}, no-bodies {:.4}",
        every_body.best() * 1e6 / analysed_bases as f64,
        needed_bodies.best() * 1e6 / analysed_bases as f64,
        no_bodies.best() * 1e6 / analysed_bases as f64,
    );

    println!(
        "ratio psp-every-body / census = {:.0}",
        every_body.best() / census.best()
    );
    println!(
        "ratio psp-needed-bodies / census = {:.0}",
        needed_bodies.best() / census.best()
    );
    println!(
        "ratio psp-no-bodies / census = {:.0}",
        no_bodies.best() / census.best()
    );
    println!(
        "psp-needed-bodies builds {} of {} bodies, 1 in {:.0}",
        needed_bodies.bodies_built,
        needed_bodies.records,
        needed_bodies.records as f64 / needed_bodies.bodies_built.max(1) as f64,
    );
    println!(
        "the selective walk's census is {} bytes against the full walk's {}: {}",
        needed_bodies.bytes,
        census_bytes,
        if needed_bodies.bytes == census_bytes {
            "the same length"
        } else {
            "a different length"
        },
    );
    Ok(())
}

/// Read every value the census holds, and say how many there were.
///
/// **A decode whose arrays are never read is not the work a fit does**, and on packed arrays the
/// difference between decoding and reading can be the whole cost.
fn touch_every_value(
    census: &mut SampleCensusEvidence,
) -> Result<u64, pop_var_caller::ng::parameter_estimation::joint::census::CensusError> {
    let groups = census.read_groups();
    let strata = census.strata();
    let mut seen = 0_u64;
    let positions = census.with_generic(&groups, |lent| {
        let mut seen = 0_u64;
        for evidence in lent {
            let depth = evidence.depth();
            for index in 0..depth.len() {
                if !matches!(
                    depth.get(index),
                    pop_var_caller::ng::parameter_estimation::joint::census::DepthCode::NeverWalked
                ) {
                    seen += 1;
                }
            }
        }
        seen
    })?;
    seen += positions;
    for group in &groups {
        for stratum in &strata {
            seen += census.with_strata(*group, std::slice::from_ref(stratum), |lent| {
                let mut seen = 0_u64;
                for section in lent {
                    for locus in 0..section.len() {
                        seen += u64::from(section.offsets(locus).total());
                    }
                }
                seen
            })?;
        }
    }
    Ok(seen)
}

/// Whether any kept census position falls inside `region`.
fn holds_a_kept_position(kept: &[GenomePosition], region: GenomeRegion) -> bool {
    let first = kept.partition_point(|position| {
        (position.contig.get(), position.position.get()) < (region.contig.get(), region.start.get())
    });
    kept.get(first).is_some_and(|position| {
        position.contig == region.contig && position.position.get() <= region.end.get()
    })
}

/// The kept tracts as (contig, start, end) triples, ascending — what a head is compared against.
///
/// **A triple rather than a `GenomeRegion` because that type is deliberately not `Ord`**, and a
/// probe is not the place to give it an ordering the library has declined to.
fn kept_tract_regions(
    plan: &CensusPlan,
    reference: &pop_var_caller::ng::reference_info::ReferenceInfo,
) -> Vec<(u32, u64, u64)> {
    let contigs = reference.contig_list();
    let id_of = |name: &str| {
        contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .map(|index| index as u32)
    };
    let mut regions: Vec<(u32, u64, u64)> = plan
        .loci
        .ssr()
        .iter_sorted()
        .into_iter()
        .flat_map(|(_, segments)| segments.iter())
        .filter_map(|segment| {
            id_of(segment.chrom()).map(|contig| (contig, segment.start(), segment.end()))
        })
        .collect();
    regions.sort_unstable();
    regions
}

/// What a selective walk produced.
struct SelectiveCensus {
    evidence: SampleCensusEvidence,
    records: u64,
    bodies_built: u64,
}

/// A census built from a psp, declining every body the predicate does not want.
fn selective_census(
    path: &Path,
    plan: &CensusPlan,
    segmentation: &Segmentation,
    mut want: impl FnMut(&RecordHead) -> bool,
) -> Result<SelectiveCensus, Box<dyn Error>> {
    let mut reader = PspReader::open(path)?;
    let sample = reader.header().sample.clone();
    let declared: std::collections::BTreeMap<ReadGroupId, NamedReadGroup> = reader
        .header()
        .read_groups
        .iter()
        .map(|group| {
            (
                group.walk_local_id,
                NamedReadGroup {
                    declared_id: group.id.clone(),
                    library: group.library.clone(),
                },
            )
        })
        .collect();
    let mut writer: CensusWriter = plan.writer_for(sample, declared, segmentation);
    let mut records = 0_u64;
    let mut bodies_built = 0_u64;
    for streamed in reader.records_where(&mut want)? {
        let streamed = streamed?;
        records += 1;
        if let Some(record) = streamed.record.as_ref() {
            bodies_built += 1;
            writer.add_locus(record);
        }
    }
    Ok(SelectiveCensus {
        evidence: writer.finish(),
        records,
        bodies_built,
    })
}

/// Every head, no body — the floor a psp route cannot go below.
///
/// **It also measures the depth this psp was written at**, from the heads it is reading anyway:
/// reads compared with the reference, summed over records and divided by the bases they cover.
/// That is the number every figure here has to be quoted beside, and taking it from the file
/// removes any chance of labelling a measurement with a depth somebody remembered.
fn heads_only(path: &Path) -> Result<(u64, f64), Box<dyn Error>> {
    let mut reader = PspReader::open(path)?;
    let mut records = 0_u64;
    let mut reads = 0_u64;
    let mut bases = 0_u64;
    for streamed in reader.records_where(|_| false)? {
        let streamed = streamed?;
        records += 1;
        let span = streamed.head.region.len();
        reads += u64::from(streamed.head.reads_compared_with_reference) * span;
        bases += span;
    }
    Ok((records, reads as f64 / bases.max(1) as f64))
}
