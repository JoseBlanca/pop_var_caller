//! **Where a calling run has got to**, printed to stderr every
//! [`LOCI_BETWEEN_PROGRESS_LINES`] loci.
//!
//! A calling run writes nothing a person can watch until it ends: with the hidden-duplication
//! filter on, every record is parked on a temporary file and the VCF is only written once the
//! whole genome has been called. On a cohort that takes hours, a silent console is
//! indistinguishable from a hung one. So the calling loop counts each locus it is handed and,
//! every five million, says where it is:
//!
//! ```text
//! calling: SL4.0ch06:18,783,657, 47% of the analysed bases, 10,000,000 loci, 4,068,400 records, 1h12m elapsed
//! ```
//!
//! **It sits on the calling thread, where the loci arrive one at a time in genome order**, so the
//! lines cannot interleave and the position only moves forward. The percentage is the analysed
//! bases before the locus over all the analysed bases — the whole genome, or the regions a BED
//! file named — and "records" counts the records handed on so far, before any filter. Counting a
//! locus is one increment and one comparison, and nothing here touches what the run writes.

use std::time::{Duration, Instant};

use crate::fasta::ContigList;
use crate::types::GenomeRegion;

/// How many loci pass between two progress lines.
///
/// **A locus is far sparser than a position.** On eight tomato accessions at about three reads a
/// position, 8 Mb of analysed regions gave 208,956 loci — one per 38 bases — so five million
/// loci is some 190 Mb of genome there. How dense loci are grows with the cohort and its depth,
/// since more reads show more differences from the reference, so the same count is a different
/// stretch of genome on another run.
pub(crate) const LOCI_BETWEEN_PROGRESS_LINES: u64 = 5_000_000;

/// The calling loop's counter of loci, and what it needs to say where a locus is.
pub(crate) struct CallingProgress<'a> {
    contigs: &'a ContigList,
    /// The analysed regions in genome order — the ground the merge hands loci over from.
    analysed: &'a [GenomeRegion],
    analysed_bases: u64,
    /// The analysed region the last locus fell in, and the bases of every region before it —
    /// a cursor, because loci arrive in genome order and never move back.
    region_index: usize,
    bases_before_region: u64,
    loci: u64,
    started: Instant,
}

impl<'a> CallingProgress<'a> {
    /// Start counting a run over `analysed`, whose contigs `contigs` names.
    pub(crate) fn new(contigs: &'a ContigList, analysed: &'a [GenomeRegion]) -> Self {
        Self {
            contigs,
            analysed,
            analysed_bases: analysed.iter().map(|region| region.len()).sum(),
            region_index: 0,
            bases_before_region: 0,
            loci: 0,
            started: Instant::now(),
        }
    }

    /// Count one locus, starting at `locus`; every [`LOCI_BETWEEN_PROGRESS_LINES`]th prints a progress line.
    pub(crate) fn locus_passed(&mut self, locus: GenomeRegion, records_written: u64) {
        self.loci += 1;
        if self.loci.is_multiple_of(LOCI_BETWEEN_PROGRESS_LINES) {
            let line = self.line(locus, records_written, self.started.elapsed());
            eprintln!("{line}");
        }
    }

    /// The last line of a run: how many loci it called and how long that took.
    pub(crate) fn finished(&self, records_written: u64) {
        eprintln!(
            "calling: done, {} loci, {} records, {} elapsed",
            with_thousands(self.loci),
            with_thousands(records_written),
            elapsed(self.started.elapsed()),
        );
    }

    fn line(&mut self, locus: GenomeRegion, records_written: u64, took: Duration) -> String {
        let contig = self
            .contigs
            .entries
            .get(locus.contig.0 as usize)
            .map_or("?", |entry| entry.name.as_str());
        let percent = match self.analysed_bases {
            0 => 0,
            all => self.bases_before(locus) * 100 / all,
        };
        format!(
            "calling: {contig}:{}, {percent}% of the analysed bases, {} loci, {} records, {} elapsed",
            with_thousands(locus.start.get()),
            with_thousands(self.loci),
            with_thousands(records_written),
            elapsed(took),
        )
    }

    /// The analysed bases before `locus` starts. Advances the cursor past every region that
    /// ends before it.
    fn bases_before(&mut self, locus: GenomeRegion) -> u64 {
        while let Some(region) = self.analysed.get(self.region_index)
            && (region.contig, region.end) < (locus.contig, locus.start)
        {
            self.bases_before_region += region.len();
            self.region_index += 1;
        }
        let within = match self.analysed.get(self.region_index) {
            Some(region) if region.contig == locus.contig && region.start <= locus.start => {
                locus.start.get() - region.start.get()
            }
            _ => 0,
        };
        self.bases_before_region + within
    }
}

/// `1234567` as `1,234,567`.
fn with_thousands(number: u64) -> String {
    let digits = number.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

/// A duration as `1h12m`, `12m05s` or `42s` — to the second, since a line comes every few
/// minutes at the least.
fn elapsed(took: Duration) -> String {
    let seconds = took.as_secs();
    let (hours, minutes, seconds) = (seconds / 3600, seconds / 60 % 60, seconds % 60);
    if hours > 0 {
        format!("{hours}h{minutes:02}m")
    } else if minutes > 0 {
        format!("{minutes}m{seconds:02}s")
    } else {
        format!("{seconds}s")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::ContigEntry;
    use crate::types::{ContigId, Position};

    fn contigs() -> ContigList {
        ContigList {
            entries: ["chr1", "chr2"]
                .into_iter()
                .map(|name| ContigEntry {
                    name: name.to_owned(),
                    length: 1_000,
                    md5: None,
                })
                .collect(),
        }
    }

    fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    #[test]
    fn numbers_are_grouped_in_thousands() {
        assert_eq!(with_thousands(0), "0");
        assert_eq!(with_thousands(999), "999");
        assert_eq!(with_thousands(1_000), "1,000");
        assert_eq!(with_thousands(41_250_000), "41,250,000");
    }

    #[test]
    fn durations_are_read_at_the_largest_unit() {
        assert_eq!(elapsed(Duration::from_secs(42)), "42s");
        assert_eq!(elapsed(Duration::from_secs(12 * 60 + 5)), "12m05s");
        assert_eq!(elapsed(Duration::from_secs(3600 + 12 * 60 + 59)), "1h12m");
    }

    /// Two analysed regions of 100 bases each on chr1 and one of 200 on chr2: a locus at the
    /// start of the chr2 region has half the analysed bases behind it, whatever lies between
    /// the regions.
    #[test]
    fn the_percentage_counts_only_analysed_bases() {
        let contigs = contigs();
        let analysed = [region(0, 1, 100), region(0, 501, 600), region(1, 101, 300)];
        let mut progress = CallingProgress::new(&contigs, &analysed);

        let line = progress.line(region(0, 551, 551), 3, Duration::from_secs(7));
        // 100 bases of the first region and 50 of the second, of 400.
        assert_eq!(
            line,
            "calling: chr1:551, 37% of the analysed bases, 0 loci, 3 records, 7s elapsed"
        );
        let line = progress.line(region(1, 101, 101), 3, Duration::from_secs(7));
        assert!(
            line.starts_with("calling: chr2:101, 50% of the analysed bases"),
            "{line}"
        );
        let line = progress.line(region(1, 201, 201), 3, Duration::from_secs(7));
        assert!(line.contains("75% of the analysed bases"), "{line}");
    }
}
