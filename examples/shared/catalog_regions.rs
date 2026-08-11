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
