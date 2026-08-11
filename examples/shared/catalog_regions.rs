//! How every dev tool gets its typed regions now: **from the reference's repeat catalog, not
//! from a scan.**
//!
//! Each of these tools used to build a walk over the FASTA and let it cut the reference into
//! typed regions as it went. The catalog holds those repeats already — found once per
//! reference by `pop_var_caller_exp repeat-catalog` — so a tool opens the file beside the
//! reference and reads them back. It costs a run nothing to scan, and every tool sees the
//! same repeats rather than each rediscovering them.
//!
//! **The precondition that comes with it:** a reference with no catalog beside it stops the
//! run, and the error names the command that writes one. That is deliberate — falling back to
//! a scan would answer the question a second way without saying so.
//!
//! Two differences from the walk are worth knowing before reading a tool's output:
//!
//! - a tandem repeat within 15 bases of a contig's end is not in the file, so a scan would
//!   have found a handful of loci per genome that a tool no longer reports (0 in tomato at
//!   the calling floors, 4 at the catalog's own);
//! - the end-of-run tally has no counter for repeats turned down for touching a contig's very
//!   first or last base, because the file holds none of them and a `0` would read as "this
//!   genome has none".

use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

/// One contig, end to end, as the region list the catalog is asked with.
///
/// `length` is the contig's own, from the reference's table — never a window's and never a
/// slice's, for the reason the flank test needs it: whether a tract has room beside it is a
/// question about the contig.
#[allow(dead_code)]
pub fn whole_contig(contig: ContigId, length: u64) -> GenomeRegion {
    GenomeRegion {
        contig,
        start: Position(1),
        end: Position(length),
    }
}

/// Build the reference's catalog, for a **synthetic fixture that has none**.
///
/// A real reference gets one from `pop_var_caller_exp repeat-catalog`, once, and every tool
/// reads it. A fixture a test writes has no such run behind it, so the test does what that
/// command does: streams the reference once and writes the catalog beside it.
#[allow(dead_code)]
pub fn build_catalog_beside(fasta: &std::path::Path) {
    use pop_var_caller::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use pop_var_caller::ng::repeat_catalog::{
        RepeatCatalogBuilder, StrRepeatCriteria, sibling_catalog_path,
    };
    use pop_var_caller::ng::tandem_repeat::ScanParams;

    let mut builder = RepeatCatalogBuilder::create(
        &sibling_catalog_path(fasta),
        StrRepeatCriteria::default(),
        ScanParams::default(),
    )
    .expect("a builder");
    let info = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: fasta.to_path_buf(),
            fai: None,
        },
        &mut builder,
    )
    .expect("the reference pass runs");
    builder.finish(&info).expect("the catalog is written");
}
