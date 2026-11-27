//! The reference's contig table, and what is left of the reference-FASTA fetchers.
//!
//! - [`ContigEntry`] / [`ContigList`] — the reference's contig table (name, length, optional
//!   MD5). Shared across the FASTA, CRAM and psp formats because every input asserts which
//!   reference it is against.
//! - [`ChromRefFetchError`] — the error a reference read raises, and [`fetcher`]'s
//!   canonicalisation rule. The fetchers that module was named for served the deleted production
//!   caller and are gone; the readers in use are [`crate::raw_chrom_reader`] and
//!   [`crate::ref_seq`].

pub mod fetcher;

use std::collections::HashMap;

pub use fetcher::ChromRefFetchError;

use crate::types::ContigId;

// ---------------------------------------------------------------------
// ContigList
// ---------------------------------------------------------------------

/// One reference sequence (`@SQ`) entry: name, length, and optional MD5.
#[derive(Debug, Clone)]
pub struct ContigEntry {
    pub name: String,
    pub length: u64,
    pub md5: Option<[u8; 16]>,
}

impl PartialEq for ContigEntry {
    fn eq(&self, other: &Self) -> bool {
        if self.name != other.name || self.length != other.length {
            return false;
        }
        match (self.md5, other.md5) {
            (Some(a), Some(b)) => a == b,
            // Absent MD5 acts as a wildcard: a CRAM that omits M5 is
            // not contradicting a CRAM that carries one.
            (None, _) | (_, None) => true,
        }
    }
}

impl Eq for ContigEntry {}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ContigList {
    pub entries: Vec<ContigEntry>,
}

impl ContigList {
    /// First-difference report between two lists, intended for error
    /// messages. Returns `Ok(())` when the lists agree; otherwise an
    /// `Err` carrying a short string identifying which field
    /// disagreed and where.
    pub(crate) fn first_disagreement(&self, other: &Self) -> Result<(), String> {
        if self.entries.len() != other.entries.len() {
            return Err(format!(
                "@SQ list length differs ({} vs {})",
                self.entries.len(),
                other.entries.len()
            ));
        }
        for (i, (a, b)) in self.entries.iter().zip(other.entries.iter()).enumerate() {
            if a.name != b.name {
                return Err(format!(
                    "name disagreement at index {} ('{}' vs '{}')",
                    i, a.name, b.name
                ));
            }
            if a.length != b.length {
                return Err(format!(
                    "length disagreement at index {} (contig '{}': {} vs {})",
                    i, a.name, a.length, b.length
                ));
            }
            if let (Some(ma), Some(mb)) = (a.md5, b.md5)
                && ma != mb
            {
                return Err(format!(
                    "md5 disagreement at index {} (contig '{}')",
                    i, a.name
                ));
            }
        }
        Ok(())
    }

    /// A lookup from each contig's name to its [`ContigId`] — its position in this list —
    /// built once, for a caller that resolves a name per locus.
    pub fn ids_by_name(&self) -> ContigIdsByName {
        ContigIdsByName::of(self)
    }
}

/// **Contig name → [`ContigId`], built once from a [`ContigList`].** The id is the contig's
/// position in the list, the number every record carries.
///
/// A map rather than a scan because the census asks once per STR locus it keeps, in every
/// sample: on a scaffold-level assembly of 30,000 contigs a scan averages 15,000 string
/// comparisons a lookup, about 1.5 billion a sample for 100,000 loci. The map costs about a
/// megabyte at that contig count, and nothing on a reference of a few chromosomes.
///
/// A name the list holds twice resolves to its first position, as a scan from the front would.
#[derive(Debug, Clone, Default)]
pub struct ContigIdsByName {
    by_name: HashMap<String, ContigId>,
}

impl ContigIdsByName {
    fn of(contigs: &ContigList) -> Self {
        let mut by_name = HashMap::with_capacity(contigs.entries.len());
        for (index, entry) in contigs.entries.iter().enumerate() {
            by_name
                .entry(entry.name.clone())
                .or_insert(ContigId(index as u32));
        }
        Self { by_name }
    }

    /// The id of the contig named `name`, or `None` when the list has no such contig.
    pub fn id_of(&self, name: &str) -> Option<ContigId> {
        self.by_name.get(name).copied()
    }
}

// ---------------------------------------------------------------------
// MultiChromRefFetcher trait
// ---------------------------------------------------------------------

#[cfg(test)]
/// The deleted production pileup walker's multi-chromosome reference interface, kept for tests
/// only: the generic locus generator's mock reference
/// (`locus_generation::pileup::mock_reference`) implements it so code copied from that walker
/// can be tested unchanged. Errors are typed ([`ChromRefFetchError`]) so callers can route I/O
/// failures, range failures and contract violations apart.
pub trait MultiChromRefFetcher {
    /// Fetch `length` reference bases starting at the 1-based
    /// position `start` on chromosome `chrom_id`. Bytes are
    /// uppercase ASCII over `{A,C,G,T,N}` (canonicalised by the
    /// fetcher implementation).
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::OutOfBounds`] if the requested
    ///   window exceeds the chromosome length.
    /// - [`ChromRefFetchError::InvalidStart`] if `start_1based == 0`.
    /// - [`ChromRefFetchError::Io`] on any underlying FASTA I/O
    ///   failure or unknown `chrom_id`.
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError>;
}

#[cfg(test)]
/// Forwarding impl so callers may pass either an owned fetcher or a
/// shared reference. The fetcher's methods all take `&self`, so the borrow is
/// sufficient. (Production's pileup walker, deleted with production, was the
/// caller that needed it.)
impl<T: MultiChromRefFetcher + ?Sized> MultiChromRefFetcher for &T {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        (**self).fetch(chrom_id, start_1based, length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(name: &str) -> ContigEntry {
        ContigEntry {
            name: name.to_owned(),
            length: 100,
            md5: None,
        }
    }

    /// The id is the position in the list; an absent name has none; a name listed twice
    /// resolves to its first position, as a scan from the front of the list would.
    #[test]
    fn a_name_resolves_to_its_first_position_in_the_list() {
        let contigs = ContigList {
            entries: vec![entry("chr1"), entry("chr2"), entry("chr1"), entry("chr3")],
        };
        let ids = contigs.ids_by_name();
        assert_eq!(ids.id_of("chr1"), Some(ContigId(0)));
        assert_eq!(ids.id_of("chr2"), Some(ContigId(1)));
        assert_eq!(ids.id_of("chr3"), Some(ContigId(3)));
        assert_eq!(ids.id_of("chrUn"), None);
        assert_eq!(ids.id_of(""), None);
    }
}
