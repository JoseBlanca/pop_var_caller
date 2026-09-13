//! **The committed golden repeat catalog, read by ng's own tests.**
//!
//! `tests/data/tandem_repeat/golden.ssr_catalog.bed.gz` is the catalog production's `ssr-catalog`
//! built — `trf-mod`, then its post-filter — on the synthetic STR-diversity reference beside it. Four
//! of ng's tests treat it as an oracle: the scanner's recall, the resident partition, the repeat
//! catalog's shipping stack, and the reference digest. They read it through production's
//! `ssr::catalog::io::CatalogReader`; production is being deleted, so this is the part of that
//! reader the tests use, written for the one file they read. The scanner's test moved to it at
//! promotion step C12, and the other three move at C19.
//!
//! The file is BGZF-compressed text: a `## key: value` metadata block, a `#`-prefixed column header,
//! and one tab-separated row per locus — `chrom start end motif purity_fraction ref_seq_start
//! ref_seq`, with 0-based half-open coordinates. Only what the tests read is kept.

use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use noodles_bgzf as bgzf;

/// The post-filter settings the golden catalog was built at, from its metadata block.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct GoldenCatalogSettings {
    /// The purity floor applied after the post-filter recomputed purity.
    pub(crate) min_purity: f32,
    /// The early floor on the detector's score.
    pub(crate) min_score: i32,
    /// The flank, in bases, embedded each side of a tract.
    pub(crate) flank_bp: u32,
    /// The radius, in bases, within which neighbouring tracts were dropped as a bundle.
    pub(crate) bundle_threshold: u32,
}

/// One locus of the golden catalog: where its tract is and what repeats in it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GoldenLocus {
    /// The contig's name.
    pub(crate) chrom: String,
    /// The tract's first base, 0-based.
    pub(crate) start: u32,
    /// One past the tract's last base.
    pub(crate) end: u32,
    /// The repeat unit, as the post-filter wrote it.
    pub(crate) motif: Vec<u8>,
}

/// The whole golden catalog, as the tests need it.
pub(crate) struct GoldenCatalog {
    /// The MD5 of the reference the catalog was built on, as its metadata block records it.
    pub(crate) reference_md5: String,
    /// The settings it was built at.
    pub(crate) settings: GoldenCatalogSettings,
    /// Its loci, in file order.
    pub(crate) loci: Vec<GoldenLocus>,
}

/// The path of a file in the tandem-repeat test data directory.
pub(crate) fn tandem_repeat_fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("tandem_repeat")
        .join(name)
}

/// Read the golden catalog. Any line that is not what the format promises is a failure, not a skip:
/// a silently shortened oracle would let a recall test pass on fewer loci.
pub(crate) fn golden_catalog() -> GoldenCatalog {
    let file = std::fs::File::open(tandem_repeat_fixture("golden.ssr_catalog.bed.gz"))
        .expect("the golden catalog is committed");
    let reader = BufReader::new(bgzf::io::Reader::new(file));

    let mut metadata: Vec<(String, String)> = Vec::new();
    let mut loci = Vec::new();
    for line in reader.lines() {
        let line = line.expect("the golden catalog decompresses");
        if let Some(entry) = line.strip_prefix("## ") {
            let (key, value) = entry
                .split_once(": ")
                .unwrap_or_else(|| panic!("a metadata line is `## key: value`: {line}"));
            metadata.push((key.to_string(), value.to_string()));
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let [chrom, start, end, motif, _purity, _ref_seq_start, _ref_seq] = fields[..] else {
            panic!("a golden locus row has seven columns: {line}");
        };
        loci.push(GoldenLocus {
            chrom: chrom.to_string(),
            start: start
                .parse()
                .unwrap_or_else(|_| panic!("a tract start: {line}")),
            end: end
                .parse()
                .unwrap_or_else(|_| panic!("a tract end: {line}")),
            motif: motif.as_bytes().to_vec(),
        });
    }

    let value = |key: &str| {
        metadata
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
            .unwrap_or_else(|| panic!("the golden catalog records `{key}`"))
    };
    GoldenCatalog {
        reference_md5: value("reference_md5").to_string(),
        settings: GoldenCatalogSettings {
            min_purity: value("min_purity").parse().expect("a purity floor"),
            min_score: value("min_score").parse().expect("a score floor"),
            flank_bp: value("flank_bp").parse().expect("a flank"),
            bundle_threshold: value("bundle_threshold").parse().expect("a bundle radius"),
        },
        loci,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The reader returns what the file holds**, counted by hand from it: sixteen loci on two
    /// contigs, and the four settings production's catalog header wrote.
    #[test]
    fn the_golden_catalog_reads_as_written() {
        let catalog = golden_catalog();
        assert_eq!(catalog.loci.len(), 16);
        assert_eq!(catalog.reference_md5, "8a23ad24367daeacae74b8cefdd7cdf4");
        assert_eq!(
            catalog.settings,
            GoldenCatalogSettings {
                min_purity: 0.8,
                min_score: 0,
                flank_bp: 50,
                bundle_threshold: 50,
            }
        );
        assert_eq!(
            catalog.loci[0],
            GoldenLocus {
                chrom: "ctg1".to_string(),
                start: 119,
                end: 155,
                motif: b"TA".to_vec(),
            }
        );
        assert_eq!(
            catalog.loci[15],
            GoldenLocus {
                chrom: "ctg2".to_string(),
                start: 69,
                end: 109,
                motif: b"CTAT".to_vec(),
            }
        );
    }
}
