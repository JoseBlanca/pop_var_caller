//! End-to-end tests for the pileup walker plus shared helpers
//! (`MockFasta`, fixture builders) used by inner-module tests.
//!
//! Tests exercise scenarios from `ia/specs/pileup_walker.md`: SNP
//! folding, deletion anchoring, REF widening on overlap, mate
//! overlap, phase-chain lifecycle markers, eager closure, etc.
//!
//! **Copied verbatim from `src/pileup/walker/tests.rs`, and that is the
//! milestone's gate** (`locus_generation_pileup.md` §12): production's own
//! suite, unmodified, green against ng's copy of the walker. Anything here
//! that needed a *behavioural* touch would be a transcription error, not a
//! design change. The only edit is the one the type forces — the 24
//! `PreparedRead` literals name the new `read_group` field, which is exactly
//! what not making that type `#[non_exhaustive]` was for. Strip those 24 lines
//! and this file is byte-identical to production's.
//!
//! It dies with plan 3, when ng's walker starts emitting
//! `SampleLocusObservations` and its behaviour deliberately diverges; §12
//! classifies each test three ways for that moment.
//!
//! `copy_fidelity.rs`'s `the_eight_copies_are_still_productions` is what keeps
//! "verbatim" a checked property rather than a claim in a doc comment — and it checks
//! this file too, which is why it is not in it.

use std::sync::Arc;

use super::CigarOp;
use super::MateRole;
use super::PreparedRead;
use super::WalkerConfig;
use super::run;
use crate::fasta::{ChromRefFetchError, MultiChromRefFetcher};
use crate::ng::locus_generation::{ReadWitness, SequenceObservation, WitnessedLocusPositions};
use crate::ng::types::ReadGroupId;

// ---------------------------------------------------------------------
// MockFasta
// ---------------------------------------------------------------------

/// In-memory `MultiChromRefFetcher` implementation backed by a single
/// chromosome string. Tests inject their reference here directly
/// instead of building a real FASTA file.
#[derive(Debug, Clone)]
pub struct MockFasta {
    /// Reference bases, indexed 0-based for storage. The walker
    /// uses 1-based coordinates externally — `fetch` translates.
    chromosomes: Vec<Vec<u8>>,
}

impl MockFasta {
    pub fn new(chr0: &str) -> Self {
        Self {
            chromosomes: vec![chr0.as_bytes().to_vec()],
        }
    }

    /// Multi-chromosome variant: each entry is the literal bases of
    /// chromosome `i` (`chrom_id == i`). Used by tests that exercise
    /// chromosome-boundary behaviour.
    pub fn with_chromosomes(chroms: &[&str]) -> Self {
        Self {
            chromosomes: chroms.iter().map(|s| s.as_bytes().to_vec()).collect(),
        }
    }
}

impl MultiChromRefFetcher for MockFasta {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        let chrom_name = format!("chrom_id={chrom_id}");
        let chrom =
            self.chromosomes
                .get(chrom_id as usize)
                .ok_or_else(|| ChromRefFetchError::Io {
                    chrom_name: chrom_name.clone(),
                    source: std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("unknown chrom_id={chrom_id}"),
                    ),
                })?;
        if start_1based == 0 {
            return Err(ChromRefFetchError::InvalidStart);
        }
        let end_exclusive = start_1based + length;
        if (end_exclusive - 1) as usize > chrom.len() {
            return Err(ChromRefFetchError::OutOfBounds {
                chrom_name,
                chrom_length: chrom.len() as u32,
                start: start_1based,
                end: end_exclusive,
            });
        }
        let start_idx = (start_1based - 1) as usize;
        let end_idx = start_idx + length as usize;
        Ok(chrom[start_idx..end_idx].to_vec())
    }
}

// ---------------------------------------------------------------------
// Fixture helpers
// ---------------------------------------------------------------------

pub fn snp_read(qname: &str, alignment_start: u32, seq: &[u8], qual: &[u8]) -> PreparedRead {
    let len = seq.len() as u32;
    PreparedRead {
        chrom_id: 0,
        alignment_start,
        alignment_end: alignment_start + len - 1,
        cigar: vec![CigarOp::Match(len)],
        seq: seq.to_vec(),
        bq_baq: qual.to_vec(),
        mq_log_err: -3.0,
        mapq: 60,
        is_reverse_strand: false,
        qname: Arc::from(qname),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),
    }
}

pub fn paired_snp_reads(
    qname: &str,
    alignment_start_a: u32,
    alignment_start_b: u32,
    seq: &[u8],
    qual: &[u8],
) -> (PreparedRead, PreparedRead) {
    let mut a = snp_read(qname, alignment_start_a, seq, qual);
    a.mate_role = MateRole::FirstOfPair;
    let mut b = snp_read(qname, alignment_start_b, seq, qual);
    b.mate_role = MateRole::SecondOfPair;
    (a, b)
}

/// **ng's locus type, which these inherited tests now assert on directly.**
///
/// Until this commit they ran through `to_pileup_record`, a back-projection onto production's
/// `PileupRecord` — the adaptation spec §12 sanctions, taken at B2 as one reviewed function
/// rather than 67 hand-edited assertions. But that projection **merged the observations ng splits** (by
/// witness and by read group) and dropped three fields production has no counterpart for, and at
/// Milestone B its losses hid three live surfaces from the review. D1 removed the *differential's*
/// need for it; this removes the suite's, so nothing in the module sees the emitted type through a
/// lossy view any more.
pub type Locus = crate::ng::locus_generation::SampleLocusObservations;

/// Production's positional allele idiom, said in ng's terms — the two accessors that carry the
/// inherited assertions across.
///
/// **They panic rather than return an `Option`, deliberately.** These tests were written against
/// a type where `alleles[0]` and `alleles[1]` always existed, and the two ways ng differs are
/// exactly what a silent `None` would hide: a locus no read matched the reference at has **no**
/// reference observation (production creates the bucket regardless of support), and a locus whose
/// reference-matching reads split by witness or read group has **several**. A test that lands on
/// either is asking a question ng's type does not answer, and should be rewritten to assert on the
/// observations it means — so it fails loudly and names which case it hit.
impl Locus {
    /// The observation for reads matching the reference across this locus — production's `alleles[0]`.
    pub fn reference_observation(&self) -> &SequenceObservation {
        let matching: Vec<&SequenceObservation> = self
            .observations
            .iter()
            .filter(|observation| observation.matches_reference(&self.reference_bases))
            .collect();
        assert_eq!(
            matching.len(),
            1,
            "the locus at {:?} has {} observations of reference bases, not one: with none, no read \
             matched the reference here and production's zero-support REF bucket has no ng \
             counterpart; with several, they differ by witness or read group and there is no \
             single REF observation to ask for. Rows: {:?}",
            self.region,
            matching.len(),
            self.observations,
        );
        matching[0]
    }

    /// The 1-based anchor — production's `pos`, narrowed back to the `u32` these assertions
    /// compare against.
    pub fn anchor(&self) -> u32 {
        u32::try_from(self.region.start.get()).expect("a fixture anchor fits u32")
    }

    /// How many reference positions the locus covers — production's `ref_span()`, which read it
    /// off `alleles[0].seq.len()`. ng carries the region, so this is the region's length and does
    /// not depend on an observation existing.
    pub fn footprint_len(&self) -> u32 {
        u32::try_from(self.region.len()).expect("a fixture footprint fits u32")
    }

    /// The first observation whose bases are *not* the reference's, in emission order — production's
    /// `alleles[1]`. Note ng sorts observations by bases, where production's order was bucket creation,
    /// so "the first alt" is the alphabetically first, not the first seen.
    pub fn first_alt_observation(&self) -> &SequenceObservation {
        self.observations
            .iter()
            .find(|observation| !observation.matches_reference(&self.reference_bases))
            .unwrap_or_else(|| {
                panic!(
                    "the locus at {:?} carries no non-reference observation: {:?}",
                    self.region, self.observations,
                )
            })
    }
}

/// Drive `run` on a fixed input list, collecting the emitted loci.
pub fn drive_walker(reads: Vec<PreparedRead>, ref_fetcher: MockFasta) -> Vec<Locus> {
    drive_walker_with_summary(reads, ref_fetcher).0
}

/// Same as `drive_walker` but also returns the run's
/// `RunSummary`. Useful for tests that assert on counters.
pub fn drive_walker_with_summary(
    reads: Vec<PreparedRead>,
    ref_fetcher: MockFasta,
) -> (Vec<Locus>, super::RunSummary) {
    drive_walker_with_config(reads, ref_fetcher, &WalkerConfig::default())
}

/// Drive `run` with an explicit `WalkerConfig`. Used by tests that
/// need to override defaults (e.g., the per-column depth caps to
/// trip the truncation path with a small synthetic input).
pub fn drive_walker_with_config(
    reads: Vec<PreparedRead>,
    ref_fetcher: MockFasta,
    config: &WalkerConfig,
) -> (Vec<Locus>, super::RunSummary) {
    let mut walker = run(reads, &ref_fetcher, config);
    let records: Vec<Locus> = (&mut walker)
        .map(|r| r.expect("walker yielded error"))
        .collect();
    let summary = walker.summary();
    (records, summary)
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[test]
fn pure_ref_pileup_emits_one_record_per_position_with_only_ref_allele() {
    // Reference: ACGTA at positions 1..5
    // Two reads, each spanning all 5 positions, all REF.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 1, b"ACGTA", &[30; 5]);
    let records = drive_walker(vec![r1, r2], fa);
    assert_eq!(records.len(), 5);
    for (i, rec) in records.iter().enumerate() {
        assert_eq!(rec.anchor(), (i + 1) as u32);
        assert_eq!(rec.observations.len(), 1, "REF only at clean position");
        assert_eq!(rec.reference_observation().num_obs, 2);
    }
}

#[test]
fn snp_at_one_position_emits_record_with_two_alleles() {
    // Reference ACGTA. r1 ref-everywhere, r2 has SNP G→T at pos 3.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 1, b"ACTTA", &[30; 5]);
    let records = drive_walker(vec![r1, r2], fa);
    assert_eq!(records.len(), 5);
    let rec_pos3 = &records[2];
    assert_eq!(rec_pos3.anchor(), 3);
    assert_eq!(rec_pos3.observations.len(), 2, "REF + SNP");
    // First is REF (G), supported by 1 read.
    assert_eq!(&*rec_pos3.reference_bases, b"G");
    assert_eq!(rec_pos3.reference_observation().num_obs, 1);
    // Second is SNP (T), supported by 1 read.
    assert_eq!(&*rec_pos3.first_alt_observation().bases, b"T");
    assert_eq!(rec_pos3.first_alt_observation().num_obs, 1);
}

#[test]
fn deletion_record_has_extended_ref_span() {
    // Reference AAAATTTTGG. One read with CIGAR 4M3D3M starting at 1:
    //   AAAA at 1..4, deletion at 5..7, TGG at 8..10.
    // Expected: anchor record at 4 (one-before deletion start),
    // REF "ATTT" (4 bases), DEL allele "A" (anchor only).
    let fa = MockFasta::new("AAAATTTGG");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 9,
        cigar: vec![CigarOp::Match(4), CigarOp::Deletion(3), CigarOp::Match(2)],
        seq: b"AAAAGG".to_vec(),
        bq_baq: vec![30; 6],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let records = drive_walker(vec![r], fa);
    let anchor = records
        .iter()
        .find(|r| r.anchor() == 4)
        .expect("must emit anchor at deletion's preceding base");
    assert_eq!(anchor.footprint_len(), 4, "anchor + 3 deleted = 4");
    assert_eq!(
        &*anchor.reference_bases, b"ATTT",
        "REF over the deletion span"
    );
    let del = anchor
        .observations
        .iter()
        .find(|a| a.bases.as_ref() == b"A")
        .expect("DEL allele = anchor only");
    assert_eq!(del.num_obs, 1);
}

#[test]
fn deletion_record_does_not_double_count_ref_reads() {
    // Reference: ACGTAC (positions 1..6).
    // r1: pure-Match across 1..5 (5M). All REF.
    // r2: 1M3D1M starting at 2 → Match at 2, Deletion of 3 anchored
    //     at 2, Match at 6.
    //
    // The deletion record at pos=2 widens to span 4 (footprint
    // [2, 6)). r1 spans the whole record's footprint with REF
    // bases, so REF.num_obs should be 1 (one ref-spanning read);
    // before the fix it was 4 because r1 was re-folded once per
    // walker step inside the open record's footprint, multiplying
    // every five-scalar value by ref_span.
    let fa = MockFasta::new("ACGTAC");
    let r1 = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACGTA".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let r2 = PreparedRead {
        chrom_id: 0,
        alignment_start: 2,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(1), CigarOp::Deletion(3), CigarOp::Match(1)],
        seq: b"CC".to_vec(),
        bq_baq: vec![30; 2],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r2"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let records = drive_walker(vec![r1, r2], fa);
    let anchor = records
        .iter()
        .find(|r| r.anchor() == 2)
        .expect("anchor at deletion's preceding base");
    assert_eq!(anchor.footprint_len(), 4, "anchor + 3 deleted = 4");
    // **`reference_observation()`, not `observations[0]`** — the conversion at `dad9baf` used the
    // subscript here and it does not mean what production's `alleles[0]` meant. `finalise` sorts
    // observations by `bases`, this locus holds `"C"` (r2's deletion) and `"CGTA"` (r1's reference
    // match), and `"C"` sorts first — so index 0 was r2's observation, the same observation the `del` assertion
    // below finds by name and checks again. The reference observation was never inspected. Production's
    // index 0 was positional and *was* the REF bucket; ng's is lexicographic — same subscript,
    // different meaning, no compiler error, which is why `reference_observation()` exists and why
    // nothing in this module should reach an observation by index again.
    //
    // **What this assertion catches, measured, and what it does not.** It is live: adding each
    // read's contribution to its observation twice makes it read 2 against 1. But it does **not** catch
    // the bucket-accounting bug the test is named for, and neither did the version before this
    // fix — because that bug stopped being representable at B1. Rows are derived from
    // `folded_reads`, one entry per read, each carrying `num_obs: 1`, so an observation's `num_obs` is
    // the number of distinct reads folded into it and no error in the *buckets* can reach it.
    // Deleting the re-fold's `subtract_contribution` outright is caught by
    // `observations_are_the_buckets_when_one_read_group_witnesses_completely`, not here.
    //
    // So the name now over-promises: what is actually pinned is "one read contributes
    // `num_obs: 1` to the observation it lands in", plus `reference_observation()`'s own
    // demand that exactly one observation carry the footprint's reference bases. The
    // historical guard is B1's per-read derivation, structurally.
    let ref_allele = anchor.reference_observation();
    assert_eq!(
        ref_allele.num_obs, 1,
        "REF: 1 obs from r1 only; got {}",
        ref_allele.num_obs
    );
    assert_eq!(ref_allele.num_fwd, 1, "REF: forward strand count = 1");
    let del = anchor
        .observations
        .iter()
        .find(|a| a.bases.as_ref() == b"C")
        .expect("DEL allele = anchor base only");
    assert_eq!(del.num_obs, 1, "DEL: 1 obs from r2");
}

#[test]
fn refold_after_widen_clears_chain_id_from_old_bucket() {
    // Reference ACGTAC (positions 1..6).
    //
    // R0 (1M2D1M @ pos 1): Match(A @ 1), Del(CG @ 2..3),
    //     Match(T @ 4). seq "AT". Opens record at pos 1
    //     directly at span 3 (footprint [1, 4)).
    // R1 (5M @ pos 1) with a T→C SNP at pos 4: seq "ACGCA".
    //     At walker_pos 1, R1 folds into the record at pos 1
    //     with events overlapping [1, 4): three Matches at
    //     pos 1..3, all REF — R1 lands in the REF bucket "ACG".
    // R3 (1M1D1M @ pos 3): Match(G @ 3), Del(T @ 4),
    //     Match(A @ 5). seq "GA". At walker_pos 3, R3's
    //     deletion has footprint [3, 5), which overlaps the
    //     existing record at pos 1 (footprint [1, 4)) and
    //     extends past it — widens the record to span 4
    //     (footprint [1, 5)).
    //
    // At walker_pos 3, R1 is a contributor (Match @ 3) and
    // re-folds under the now-wider REF "ACGT": its event
    // window now includes Match(4) = C (SNP), so the new
    // allele seq is "ACGC" — a *different* bucket from REF
    // "ACGT" (the REF bucket's seq was extended in-place by
    // `widen()` from "ACG" to "ACGT").
    //
    // Walker invariant under test: after R1 re-folds out of
    // the REF bucket into "ACGC", no bucket may retain a chain
    // id without a backing observation. Chain ids are derived
    // from `folded_reads` at finalise (each read appears once,
    // under its current bucket), so a re-fold can never leave a
    // stale id behind. REF (`alleles[0]`) chain ids are dropped
    // unconditionally, so the meaningful check here is on the
    // ALT "ACGC" bucket R1 moved into: it must carry exactly
    // R1's one chain id (`len <= num_obs`, asserted below over
    // every allele). A stale, observation-less id would break
    // that bound and, downstream, hand the merger a
    // chain-anchored constituent with no observations.
    let fa = MockFasta::new("ACGTAC");
    let r0 = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 4,
        cigar: vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(1)],
        seq: b"AT".to_vec(),
        bq_baq: vec![30; 2],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r0"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let r1 = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACGCA".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let r3 = PreparedRead {
        chrom_id: 0,
        alignment_start: 3,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(1)],
        seq: b"GA".to_vec(),
        bq_baq: vec![30; 2],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r3"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let records = drive_walker(vec![r0, r1, r3], fa);
    let anchor = records
        .iter()
        .find(|r| r.anchor() == 1)
        .expect("anchor record at pos 1");

    // Universal invariant: in every emitted record, every
    // allele bucket's `chain_ids.len()` must be `<= num_obs`.
    // Each chain id represents at least one observation that
    // landed in this bucket; a chain id with no backing
    // observation is a leftover from a re-fold that the
    // walker forgot to clean up.
    for allele in &anchor.observations {
        assert!(
            allele.chain_ids.len() <= allele.num_obs as usize,
            "allele {:?} has chain_ids={:?} but num_obs={} — \
             stale chain ids exceed observations",
            std::str::from_utf8(&allele.bases).unwrap_or("<non-utf8>"),
            allele.chain_ids,
            allele.num_obs,
        );
    }
}

#[test]
fn insertion_record_has_alt_longer_than_ref() {
    // Reference AAAACGT. Read 1M2I5M = 1 M at pos 1 ("A"), 2-base
    // insertion ("XX"), 5 M from pos 2.
    let fa = MockFasta::new("AAAACGT");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(5)],
        seq: b"AXXAAACG".to_vec(),
        bq_baq: vec![30; 8],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let records = drive_walker(vec![r], fa);
    let anchor = records
        .iter()
        .find(|r| r.anchor() == 1)
        .expect("anchor at 1");
    let ins = anchor
        .observations
        .iter()
        .find(|a| a.bases.len() > anchor.footprint_len() as usize);
    assert!(ins.is_some(), "INS allele should be longer than REF");
    let ins = ins.unwrap();
    assert_eq!(&*ins.bases, b"AXX", "anchor + 2 inserted bases");
    assert_eq!(ins.num_obs, 1);
}

#[test]
fn forward_strand_count_recorded_correctly() {
    // Reference ACG. Two reads: forward, reverse. Both REF.
    let fa = MockFasta::new("ACG");
    let mut r1 = snp_read("r1", 1, b"ACG", &[30; 3]);
    r1.is_reverse_strand = false;
    let mut r2 = snp_read("r2", 1, b"ACG", &[30; 3]);
    r2.is_reverse_strand = true;
    let records = drive_walker(vec![r1, r2], fa);
    let rec = &records[0];
    assert_eq!(rec.reference_observation().num_obs, 2);
    assert_eq!(rec.reference_observation().num_fwd, 1);
}

#[test]
fn placed_left_is_per_record() {
    // Reference ACGTA. Two reads:
    //   r1 starts at pos 1, covers 1..5
    //   r2 starts at pos 3, covers 3..5
    // At record pos 3, r1 was placed_left (start=1 < 3) and r2 was not (start=3 == 3).
    //
    // **Adapted at B2, and the half that went is the point** (spec §12): this was
    // `placed_left_and_placed_start_are_per_record`, and ng no longer computes
    // `placed_start` at all — nothing consumes it, and it is a pure function of the read's
    // start against the anchor, so a later consumer re-derives it without touching the fold
    // (spec §6). `placed_left` **is** consumed — `vcf::qual_refine` turns it into the
    // read-position-bias term subtracted from QUAL — so that half stands unchanged, and
    // the count of reads that did *not* start left is asserted through `num_obs` rather
    // than through the field that is gone.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 3, b"GTA", &[30; 3]);
    let records = drive_walker(vec![r1, r2], fa);
    let rec3 = records.iter().find(|r| r.anchor() == 3).unwrap();
    assert_eq!(rec3.reference_observation().num_obs, 2);
    assert_eq!(rec3.reference_observation().placed_left, 1);
    // **A real check, not arithmetic on the two lines above it.** The first version of
    // this substitution asserted `num_obs - placed_left == 1` right after asserting
    // `num_obs == 2` and `placed_left == 1`, which is true by construction whatever the
    // walk did — the "test that cannot fail" pattern, introduced while removing a genuine
    // assertion. What the deleted `placed_start` half actually pinned is that the counter
    // is **per record**: r1 starts left of *this* record's anchor and not of the one at 1.
    let rec1 = records.iter().find(|r| r.anchor() == 1).unwrap();
    assert_eq!(
        rec1.reference_observation().placed_left,
        0,
        "at the record anchored on r1's own start, nothing is placed left of it"
    );
}

#[test]
fn uncovered_positions_produce_no_records() {
    // Reference ACGTACGTAC (10 bp). Reads at pos 1..3 and pos 7..9.
    // Positions 4..6 should produce no records.
    let fa = MockFasta::new("ACGTACGTAC");
    let r1 = snp_read("r1", 1, b"ACG", &[30; 3]);
    let r2 = snp_read("r2", 7, b"GTA", &[30; 3]);
    let records = drive_walker(vec![r1, r2], fa);
    let positions: Vec<u32> = records.iter().map(|r| r.anchor()).collect();
    assert_eq!(positions, vec![1, 2, 3, 7, 8, 9]);
}

#[test]
fn paired_mates_with_overlapping_positions_share_chain_id() {
    // Both mates start at the same position; the second admits
    // while the first is still in the active set. They collapse
    // onto the same chain id via the `pending_mates` hash lookup.
    //
    // The reads carry a non-reference base (`C` over a `A` reference)
    // so the shared chain id is observable on the ALT allele
    // (`alleles[1]`): the walker drops REF (`alleles[0]`) chain ids.
    let fa = MockFasta::new("AAAAA");
    let (m1, m2) = paired_snp_reads("pair", 1, 1, b"CCC", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec1 = records.iter().find(|r| r.anchor() == 1).unwrap();
    // **This assertion changed with the type, and the old one was checking nothing.** It read
    // `alleles[0].chain_ids.is_empty()` — "REF chain ids are dropped" — but both mates here carry
    // `C` over an `A` reference, so **no read matched the reference at all** and production's
    // `alleles[0]` was an empty bucket. An empty bucket's chain-id list is empty whatever the
    // rule does, so the assertion held for free.
    //
    // ng emits no reference observation here, which is the honest statement of the same input, and the
    // rule it was reaching for is checked where it can fail: `only_the_reads_that_departed_from_
    // the_reference_carry_a_chain_id` (in `open_record.rs`) and the dump's fixture of the same
    // name both put a genuinely reference-matching read beside a departing one.
    assert!(
        rec1.observations
            .iter()
            .all(|observation| observation.bases.as_ref() != b"A"),
        "neither mate matched the reference, so there is no reference observation: {:?}",
        rec1.observations,
    );
    assert_eq!(&*rec1.first_alt_observation().bases, b"C");
    assert_eq!(
        rec1.first_alt_observation().chain_ids,
        vec![0u64],
        "the two mates collapsed onto one chain id, which is what this test is about"
    );
}

/// **The pin on the per-column mate-overlap skip.**
///
/// `process_position` asks the active set, in O(1), whether any pair it holds still has
/// both alignments on the reference at this column, and skips the whole mate-overlap step
/// when the answer is no. A skip that fired on a column where a pair *is* present would
/// lose that reconciliation silently: the emitted record would keep both mates' BQ instead
/// of zeroing the lower one, and the only counter that would notice is this one.
///
/// So this walks a pair whose mates overlap on three bases and one solo read that overlaps
/// neither, and asserts the exact number of reconciled columns. Three, not "more than
/// zero": a skip that fired on one column of the overlap would still leave two.
#[test]
fn every_column_a_mate_pair_spans_is_reconciled_not_skipped() {
    let fa = MockFasta::new("AAAAAAAAAA");
    // Mates at 1..3 and 3..5 — one shared column, position 3.
    let (m1, mut m2) = paired_snp_reads("pair_edge", 1, 3, b"AAA", &[30; 3]);
    m2.bq_baq = vec![10; 3];
    // Mates at 6..8 and 6..8 — three shared columns, 6, 7 and 8.
    let (m3, mut m4) = paired_snp_reads("pair_full", 6, 6, b"AAA", &[30; 3]);
    m4.bq_baq = vec![10; 3];
    // A solo read over the same span as the second pair: depth without a partner, so the
    // skip has columns where it must *not* fire and columns where nothing shares an id.
    let solo = snp_read("solo", 6, b"AAA", &[30; 3]);
    let (_records, summary) = drive_walker_with_summary(vec![m1, m2, m3, solo, m4], fa);
    assert_eq!(
        summary.mate_overlap_positions, 4,
        "one shared column for the staggered pair plus three for the stacked one",
    );
}

/// The other half of the same claim: a column with no pair present must not be *charged*
/// for one. A pair whose mates never overlap — the ordinary paired-end case, and the
/// reason the skip is worth having — reconciles nothing anywhere.
#[test]
fn mates_that_never_overlap_reconcile_nothing() {
    let fa = MockFasta::new("AAAAAAAAAA");
    let (m1, m2) = paired_snp_reads("pair", 1, 5, b"AAA", &[30; 3]);
    let (_records, summary) = drive_walker_with_summary(vec![m1, m2], fa);
    assert_eq!(
        summary.mate_overlap_positions, 0,
        "the mates span 1..3 and 5..7, so no column holds both",
    );
}

#[test]
fn paired_mates_within_lookup_window_share_chain_id_across_active_set_exit() {
    // The first mate exits the active set well before the second
    // mate admits. Mate-pair tracking is governed by
    // `mate_lookup_window`, not by active-set residence — the
    // `pending_mates` entry stays alive across the first mate's
    // exit, and the second mate's later arrival (still within the
    // window) reuses the first mate's chain id.
    //
    // Fixture: m1 covers pos 1-3, m2 admits at pos 10. m1's exit
    // at walker_pos=4 must *not* drop the pending entry; the
    // 10 → 1 = 9 bp separation is well inside the default
    // `mate_lookup_window` of 10 000 bp.
    // Non-reference reads (`C` over `A`) so the shared chain id shows on
    // the ALT allele; the walker drops REF (`alleles[0]`) chain ids.
    let fa = MockFasta::new("AAAAAAAAAAAAAAAAAAAA");
    let (m1, m2) = paired_snp_reads("pair", 1, 10, b"CCC", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec1 = records.iter().find(|r| r.anchor() == 1).unwrap();
    let rec10 = records.iter().find(|r| r.anchor() == 10).unwrap();
    assert_eq!(rec1.first_alt_observation().chain_ids, vec![0u64]);
    assert_eq!(
        rec10.first_alt_observation().chain_ids,
        vec![0u64],
        "the two mates of a single pair must share one chain id"
    );
}

#[test]
fn paired_mates_separated_beyond_lookup_window_get_distinct_chain_ids() {
    // When the second mate arrives more than `mate_lookup_window`
    // bp past the first mate's `alignment_start`, the chain-id
    // allocator's `evict_stale_pending` walk has dropped the
    // pending entry by then. The second mate cannot match, mints
    // a fresh chain id, and is therefore treated as a separate
    // molecule. This is correct: when the pair is that far apart,
    // we have no read-level evidence that they're the same
    // physical fragment within our trustworthy-pairing window.
    //
    // Fixture: m1 at pos 1, m2 at pos 12_001. With default
    // `mate_lookup_window = 10_000`, m2 is past the eviction
    // threshold (1 + 10_000 + 1 = 10 002 first hits it).
    let n = 12_010_usize;
    // Non-reference reads (`C` over `A`) so the chain ids show on the ALT
    // allele; the walker drops REF (`alleles[0]`) chain ids.
    let fa = MockFasta::new(&"A".repeat(n));
    let (m1, m2) = paired_snp_reads("pair", 1, 12_001, b"CCC", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec_a = records.iter().find(|r| r.anchor() == 1).unwrap();
    let rec_b = records.iter().find(|r| r.anchor() == 12_001).unwrap();
    assert_eq!(rec_a.first_alt_observation().chain_ids, vec![0u64]);
    assert_eq!(
        rec_b.first_alt_observation().chain_ids,
        vec![1u64],
        "beyond the lookup window the pair-tracking entry has been evicted; \
         the second mate gets a fresh id"
    );
}

#[test]
fn mate_overlap_bq_tie_prefers_first_mate_not_earlier_position() {
    // Two paired mates overlapping at the same anchor with the
    // SAME alignment_start and the SAME BAQ-capped BQ. Distinguish
    // them by mq_log_err so the kept-mate's contribution to q_sum
    // is identifiable. Per spec §"Tie-breaking and disagreement":
    // BQ tie → prefer mate 1 (the read whose `mate_role` is
    // `FirstOfPair`).
    //
    // The previous tie-break used `alignment_start` and dropped
    // contributor `b` arbitrarily on equal starts, so the
    // first-of-pair distinction was ignored.
    let fa = MockFasta::new("ACG");
    let m_first = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACG".to_vec(),
        bq_baq: vec![30; 3],
        mq_log_err: -2.0, // distinct mq_log_err for the kept mate
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: MateRole::FirstOfPair,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let m_second = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACG".to_vec(),
        bq_baq: vec![30; 3],
        mq_log_err: -10.0, // distinct mq_log_err for the loser
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: MateRole::SecondOfPair,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    // First-mate appears AFTER the second mate in the input stream,
    // even though both have alignment_start = 1, to make sure the
    // tie-break is decided by `MateRole::FirstOfPair` rather than
    // by stream order or alignment_start.
    let records = drive_walker(vec![m_second, m_first], fa);
    let rec = &records[0];
    assert_eq!(rec.reference_observation().num_obs, 2);
    // Kept mate's contribution = max(ln_BQ(Q=30), -2.0) ≈ -2.0.
    // Zeroed mate contributes max(ln(1)=0, -10.0) = 0.
    // Sum ≈ -2.0. If the tie-break wrongly kept mate 2, sum would
    // be max(0, -10.0) + max(ln_BQ, -2.0) ≈ -2.0 too — but the
    // distinguishable case here is the chain id count: kept
    // mate's bq is non-zero, ln_q ≈ -2.0; loser's ln_q = 0 from
    // its zeroed BQ AND -10 mq, so its contribution to q_sum = 0.
    // Net q_sum ≈ -2.0, NOT ≈ -10.0 (which would be the case if
    // the tie-break wrongly kept mate 2).
    assert!(
        rec.reference_observation().q_sum > -3.0 && rec.reference_observation().q_sum < -1.0,
        "q_sum ≈ -2.0 (first mate kept); got {}",
        rec.reference_observation().q_sum
    );
}

#[test]
fn mate_overlap_zeroes_lower_bq_contribution() {
    // Two paired mates overlapping at the same position with
    // different BAQ-capped BQs. Match-only agree case under S7:
    // keeper carries the *summed* BQ (capped at 200), other is
    // zeroed. Both mates still count as observations.
    //
    // Mate 1 at pos 1, BQ=30, length 3.
    // Mate 2 at pos 1, BQ=10, length 3 (overlapping).
    let fa = MockFasta::new("ACG");
    let (m1, mut m2) = paired_snp_reads("p", 1, 1, b"ACG", &[30; 3]);
    m2.bq_baq = vec![10; 3];
    let records = drive_walker(vec![m1, m2], fa);
    let rec = &records[0];
    assert_eq!(rec.reference_observation().num_obs, 2, "both mates count");
    // q_sum at default mq_log_err = -3.0:
    //   keeper: max(ln_perr(40), -3.0) = -3.0  (MQ dominates)
    //   other:  max(ln(1)=0, -3.0)     = 0
    //   sum ≈ -3.0
    // (the BQ-summing change from S7 is invisible here because MQ
    // dominates; tests at low MQ_log_err pin the BQ math directly).
    assert!(
        rec.reference_observation().q_sum > -4.0 && rec.reference_observation().q_sum < -2.0,
        "expected q_sum ≈ -3 (MQ-dominated), got {}",
        rec.reference_observation().q_sum
    );
}

#[test]
fn mate_overlap_agree_keeper_carries_summed_bq() {
    // S7 agree case: when both mates' bases agree at walker_pos,
    // the surviving observation carries the *sum* of BQs (not the
    // higher mate's BQ as the original walker did). To make BQ
    // dominate q_sum so the change is observable, set a strongly-
    // negative mq_log_err so MQ never wins the max.
    let fa = MockFasta::new("A");
    let make = |is_first: bool, bq: u8| PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 1,
        cigar: vec![CigarOp::Match(1)],
        seq: b"A".to_vec(),
        bq_baq: vec![bq],
        mq_log_err: -100.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: if is_first {
            MateRole::FirstOfPair
        } else {
            MateRole::SecondOfPair
        },
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let m1 = make(true, 20);
    let m2 = make(false, 20);
    let records = drive_walker(vec![m1, m2], fa);
    assert_eq!(records.len(), 1);
    let rec = &records[0];
    assert_eq!(rec.reference_observation().num_obs, 2);
    // Combined BQ = 40. ln_perr(40) = -40 * ln(10) / 10 ≈ -9.21.
    // Keeper contribution: max(-9.21, -100) = -9.21.
    // Other contribution: max(ln_perr(0)=0, -100) = 0.
    // Total q_sum ≈ -9.21. Pre-S7 (Q=20 unsummed): ≈ -4.61.
    let q = rec.reference_observation().q_sum;
    assert!(
        q < -8.5 && q > -10.0,
        "q_sum should reflect summed BQ (≈ ln_perr(40) ≈ -9.21), got {q}",
    );
}

#[test]
fn mate_overlap_agree_combined_bq_caps_at_200() {
    // S7 agree case: the combined BQ is clamped to samtools'
    // MPLP cap of 200 (htslib/sam.c:5919-5921). 150 + 100 = 250
    // → capped at 200.
    let fa = MockFasta::new("A");
    let make = |is_first: bool, bq: u8| PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 1,
        cigar: vec![CigarOp::Match(1)],
        seq: b"A".to_vec(),
        bq_baq: vec![bq],
        mq_log_err: -100.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: if is_first {
            MateRole::FirstOfPair
        } else {
            MateRole::SecondOfPair
        },
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let m1 = make(true, 150);
    let m2 = make(false, 100);
    let records = drive_walker(vec![m1, m2], fa);
    let q = records[0].reference_observation().q_sum;
    // ln_perr(200) ≈ -46.05. Without the cap it would be
    // ln_perr(250) ≈ -57.56.
    assert!(
        q < -45.0 && q > -47.0,
        "q_sum should reflect cap-200 (≈ -46.05), got {q}",
    );
}

#[test]
fn mate_overlap_disagree_winner_bq_scaled_by_0_8() {
    // S7 disagree case: when mate bases disagree, the higher-BQ
    // mate keeps its BQ scaled by 0.8 (samtools' "we trust this
    // less" haircut at htslib/sam.c:5927); the loser is zeroed.
    let fa = MockFasta::new("A");
    let make = |is_first: bool, base: u8, bq: u8| PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 1,
        cigar: vec![CigarOp::Match(1)],
        seq: vec![base],
        bq_baq: vec![bq],
        mq_log_err: -100.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: if is_first {
            MateRole::FirstOfPair
        } else {
            MateRole::SecondOfPair
        },
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    // Mate 1 has REF "A" with higher BQ → winner. Mate 2 has SNP
    // "G" with lower BQ → loser, zeroed.
    let m1 = make(true, b'A', 30);
    let m2 = make(false, b'G', 20);
    let records = drive_walker(vec![m1, m2], fa);
    let rec = &records[0];
    let ref_allele = rec
        .observations
        .iter()
        .find(|a| a.bases.as_ref() == b"A")
        .expect("REF allele present");
    let snp_allele = rec
        .observations
        .iter()
        .find(|a| a.bases.as_ref() == b"G")
        .expect("SNP allele present");
    assert_eq!(ref_allele.num_obs, 1);
    assert_eq!(snp_allele.num_obs, 1);
    // Winner BQ = (30 * 0.8) as u8 = 24. ln_perr(24) ≈ -5.53.
    // Pre-S7 (Q=30 unscaled): ≈ -6.91.
    let q_ref = ref_allele.q_sum;
    assert!(
        q_ref < -5.0 && q_ref > -6.0,
        "REF allele q_sum should reflect scaled BQ=24 (≈ -5.53), got {q_ref}",
    );
    // Loser BQ zeroed → ln_perr(0) = 0 → max(0, -100) = 0.
    assert_eq!(
        snp_allele.q_sum, 0.0,
        "SNP allele's BQ was zeroed; q_sum should be 0",
    );
}

#[test]
fn record_widen_events_counter_only_increments_on_real_widens() {
    // Three pure-Match reads at consecutive positions on a 10-bp
    // reference. Every record opens fresh at span 1 and never
    // widens; the run's `record_widen_events` counter should
    // therefore stay at 0. The previous implementation summed
    // ref_span across all open records before/after each
    // process_position step and incremented when the after-sum
    // grew — but a freshly-opened record also grows the sum, so
    // the counter conflated opens with widens.
    let fa = MockFasta::new("ACGTACGTAC");
    let r1 = snp_read("r1", 1, b"ACG", &[30; 3]);
    let r2 = snp_read("r2", 4, b"TAC", &[30; 3]);
    let r3 = snp_read("r3", 7, b"GTA", &[30; 3]);
    let (_records, summary) = drive_walker_with_summary(vec![r1, r2, r3], fa);
    assert_eq!(
        summary.record_widen_events, 0,
        "no widening occurred; counter should be 0, got {}",
        summary.record_widen_events
    );
}

#[test]
fn paired_mate_indel_overlap_yields_single_observation() {
    // Both mates of a pair report the same insertion at the same
    // anchor. Per spec §"Mate overlap on indels": treat as one
    // observation, not two — assign the higher-BQ-proxy event to
    // the bucket and drop the other. The previous walker zeroed
    // the loser's BQ but still folded it as a separate
    // observation.
    //
    // Reference AAAACGT (length 7). Both mates' CIGAR is 1M2I5M
    // → anchor Match at 1, Insertion of 2 bp ("XX") at anchor 1,
    // five Matches over positions 2..6.
    let fa = MockFasta::new("AAAACGT");
    let cigar = vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(5)];
    let mate_a = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: cigar.clone(),
        seq: b"AXXAAACG".to_vec(),
        bq_baq: vec![30; 8],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: MateRole::FirstOfPair,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let mut mate_b = mate_a.clone();
    mate_b.mate_role = MateRole::SecondOfPair;
    mate_b.bq_baq = vec![20; 8]; // lower BQ → loser

    let records = drive_walker(vec![mate_a, mate_b], fa);
    let anchor = records
        .iter()
        .find(|r| r.anchor() == 1)
        .expect("anchor record at pos 1");
    let ins = anchor
        .observations
        .iter()
        .find(|a| a.bases.len() > anchor.footprint_len() as usize)
        .expect("INS allele present");
    assert_eq!(
        ins.num_obs, 1,
        "indel-overlap collapses to one observation; got {}",
        ins.num_obs
    );
    // Forward-strand count should also reflect a single observation.
    assert_eq!(ins.num_fwd, 1);
}

#[test]
fn record_emits_in_coordinate_order_across_reads() {
    // 100 reads at increasing starts. We just want to verify the
    // emitted record stream is monotonically ordered by pos.
    let mut chrom = String::with_capacity(200);
    for _ in 0..40 {
        chrom.push_str("ACGTA");
    }
    let fa = MockFasta::new(&chrom);
    let mut reads = Vec::new();
    for i in 0..50u32 {
        let start = i * 2 + 1;
        reads.push(snp_read(&format!("r{i}"), start, b"AC", &[30; 2]));
    }
    let records = drive_walker(reads, fa);
    // Emitted records' positions must be sorted ascending.
    for w in records.windows(2) {
        assert!(w[0].anchor() <= w[1].anchor(), "out-of-order emission");
    }
}

#[test]
fn out_of_order_input_is_a_hard_error() {
    let fa = MockFasta::new("ACGTACGTAC");
    let r1 = snp_read("r1", 5, b"ACG", &[30; 3]);
    let r2 = snp_read("r2", 1, b"ACG", &[30; 3]); // before r1 — invalid
    let err = first_walker_error(vec![r1, r2], fa);
    let msg = err.to_string();
    assert!(msg.contains("out-of-order"), "got: {msg}");
}

#[test]
fn chromosome_id_regression_is_a_hard_error() {
    // Spec invariant: `(chrom_id, alignment_start)` must be
    // monotonically non-decreasing across the input stream. The
    // within-chromosome regression case is already pinned by
    // `out_of_order_input_is_a_hard_error`; this test pins the
    // *cross-chromosome* case, which used to be silently accepted
    // because `flush_chromosome` reset `last_admitted_locus = None`
    // and the regressing read then admitted as a fresh start on
    // the smaller chrom_id. See finding M2 in
    // `ia/reviews/pileup_2026-05-09.md`.
    let fa = MockFasta::with_chromosomes(&["ACGTA", "ACGTA"]);
    let mut r1 = snp_read("a", 1, b"AC", &[30; 2]);
    r1.chrom_id = 1;
    let mut r2 = snp_read("b", 1, b"AC", &[30; 2]);
    r2.chrom_id = 0;
    let err = first_walker_error(vec![r1, r2], fa);
    let msg = err.to_string();
    assert!(
        msg.contains("out-of-order"),
        "error message should reference out-of-order semantics; got: {msg}",
    );
}

#[test]
fn forward_chromosome_change_is_accepted() {
    // Pin the legitimate forward case so the M2 fix doesn't
    // accidentally reject monotonic chrom_id transitions.
    let fa = MockFasta::with_chromosomes(&["ACG", "TTT"]);
    let mut r1 = snp_read("a", 1, b"AC", &[30; 2]);
    r1.chrom_id = 0;
    let mut r2 = snp_read("b", 1, b"TT", &[30; 2]);
    r2.chrom_id = 1;
    let records = drive_walker(vec![r1, r2], fa);
    assert!(
        records.iter().any(|r| r.region.contig.0 == 0),
        "chrom 0 records must be emitted",
    );
    assert!(
        records.iter().any(|r| r.region.contig.0 == 1),
        "chrom 1 records must be emitted",
    );
    // And they must come in chrom order.
    let chrom_ids: Vec<u32> = records.iter().map(|r| r.region.contig.0).collect();
    let mut sorted = chrom_ids.clone();
    sorted.sort();
    assert_eq!(chrom_ids, sorted, "records must emit in chrom_id order");
}

#[test]
fn chain_ids_are_unique_and_monotonically_allocated() {
    // Three solo reads admit cleanly; every chain id observed across
    // emitted records must be unique (no recycling) and the ids
    // appear in non-decreasing order of first reference position.
    // Reads carry a non-reference base (`C` over an `A` reference) so
    // the ids land on ALT alleles — the walker drops REF chain ids.
    let fa = MockFasta::new("AAAAAAAAAA");
    let reads = vec![
        snp_read("a", 1, b"CC", &[30; 2]),
        snp_read("b", 4, b"CC", &[30; 2]),
        snp_read("c", 7, b"CC", &[30; 2]),
    ];
    let records = drive_walker(reads, fa);

    // Collect every chain id from every allele observation.
    let mut all_ids: Vec<u64> = records
        .iter()
        .flat_map(|r| {
            r.observations
                .iter()
                .flat_map(|a| a.chain_ids.iter().copied())
        })
        .collect();
    let n_total_observations = all_ids.len();
    all_ids.sort_unstable();
    all_ids.dedup();
    assert_eq!(
        all_ids.len(),
        3,
        "three reads must produce exactly three distinct chain ids \
         (observed across {n_total_observations} allele observations)",
    );
    // Monotonic: starts at 0 and increases by 1.
    assert_eq!(all_ids, vec![0, 1, 2]);
}

#[test]
fn paired_mates_share_a_single_chain_id() {
    // First and second mates of a pair must collapse to one chain id.
    // Non-reference reads (`C` over `A`) so the id lands on the ALT
    // allele — the walker drops REF chain ids.
    let fa = MockFasta::new("AAA");
    let (m1, m2) = paired_snp_reads("p", 1, 1, b"CCC", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let mut all_ids: Vec<u64> = records
        .iter()
        .flat_map(|r| {
            r.observations
                .iter()
                .flat_map(|a| a.chain_ids.iter().copied())
        })
        .collect();
    all_ids.sort_unstable();
    all_ids.dedup();
    assert_eq!(
        all_ids.len(),
        1,
        "a paired read produces exactly one chain id shared across both mates",
    );
}

#[test]
fn chain_ids_persist_across_chromosome_boundaries() {
    // Chain ids are per-`.psp`-file unique, not per-chromosome. The
    // walker's chain-id allocator must NOT reset `next_id` on
    // chromosome change.
    // Non-reference reads (`C` over `A`) so the ids land on ALT alleles —
    // the walker drops REF chain ids.
    let fa = MockFasta::with_chromosomes(&["AA", "AA"]);
    let mut r0 = snp_read("a", 1, b"CC", &[30; 2]);
    r0.chrom_id = 0;
    let mut r1 = snp_read("b", 1, b"CC", &[30; 2]);
    r1.chrom_id = 1;
    let records = drive_walker(vec![r0, r1], fa);
    let mut all_ids: Vec<u64> = records
        .iter()
        .flat_map(|r| {
            r.observations
                .iter()
                .flat_map(|a| a.chain_ids.iter().copied())
        })
        .collect();
    all_ids.sort_unstable();
    all_ids.dedup();
    assert_eq!(
        all_ids,
        vec![0, 1],
        "chain ids must remain unique across chromosomes",
    );
}

#[test]
fn column_depth_cap_truncates_snp_only_column_when_over_cap() {
    // Five SNP-only reads anchored at pos 1, each spanning 5
    // bases. Every column has 5 contributors and only Match
    // events, so the SNP cap applies. With max_snp_column_depth=3
    // we expect every column to truncate to 3 contributors.
    let fa = MockFasta::new("ACGTA");
    let reads: Vec<_> = (0..5)
        .map(|i| snp_read(&format!("r{i}"), 1, b"ACGTA", &[30; 5]))
        .collect();
    let cfg = WalkerConfig {
        max_snp_column_depth: 3,
        max_indel_column_depth: 99,
        ..WalkerConfig::default()
    };
    let (records, summary) = drive_walker_with_config(reads, fa, &cfg);

    // 5 columns, all over-cap → 5 truncations.
    assert_eq!(summary.column_depth_truncations, 5);
    // No column should report num_obs > cap.
    for rec in &records {
        for allele in &rec.observations {
            assert!(
                allele.num_obs <= 3,
                "pos {}: num_obs {} should be capped at 3",
                rec.anchor(),
                allele.num_obs,
            );
        }
    }
}

/// **What survives a capped column, pinned — and the rule it pins changed on
/// 2026-08-05.**
///
/// This test was `column_depth_cap_keeps_first_n_of_admission_order`, and it asserted
/// that the walk keeps the first `cap` contributors in the order the active set holds
/// them: production's rule, and ng's until now. That rule made the emitted bases depend
/// on the container's storage order, and where arrival order survived it meant "keep the
/// leftmost-starting reads" — a subsample tilted towards early alignment starts.
///
/// ng now keeps the `cap` reads with the smallest
/// [`sampling_key`](super::read_sampling::sampling_key), a hash of the query name. So
/// the expected answer is no longer a fixed pair of names: it is worked out here from the
/// reads themselves, which is exactly the property being asserted — **the survivors are a
/// function of the reads and of nothing else.**
///
/// Five SNP reads at one position, each with a distinct read base so the survivors are
/// identifiable from the emitted observations. Cap = 2.
#[test]
fn column_depth_cap_keeps_the_smallest_sampling_keys() {
    let fa = MockFasta::new("A");
    let reads = vec![
        snp_read("r0", 1, b"C", &[30]),
        snp_read("r1", 1, b"G", &[30]),
        snp_read("r2", 1, b"T", &[30]),
        snp_read("r3", 1, b"A", &[30]), // matches REF
        snp_read("r4", 1, b"C", &[30]),
    ];
    let cfg = WalkerConfig {
        max_snp_column_depth: 2,
        max_indel_column_depth: 99,
        ..WalkerConfig::default()
    };

    // The two smallest keys, and therefore the two bases that must be folded — worked
    // out from the reads before the walk sees them.
    let mut ranked: Vec<_> = reads
        .iter()
        .map(|read| {
            (
                super::read_sampling::sampling_key(read),
                read.qname.to_string(),
                read.seq.clone(),
            )
        })
        .collect();
    ranked.sort();
    let survivors: Vec<&str> = ranked[..2]
        .iter()
        .map(|(_, name, _)| name.as_str())
        .collect();
    let mut expected_bases: Vec<Vec<u8>> =
        ranked[..2].iter().map(|(_, _, seq)| seq.clone()).collect();
    expected_bases.sort();

    let (records, summary) = drive_walker_with_config(reads, fa, &cfg);

    assert_eq!(records.len(), 1, "single column emitted");
    assert_eq!(summary.column_depth_truncations, 1);
    let rec = &records[0];
    assert_eq!(&*rec.reference_bases, b"A", "the locus's reference base");

    // Every folded read, as bases, one entry per observation counted. ng emits an
    // observation per allele that a read actually matched — including the reference
    // allele when a read matched it, and *no* reference observation when none did,
    // which is the spec §12 divergence from production this file already carried.
    let mut folded: Vec<Vec<u8>> = rec
        .observations
        .iter()
        .flat_map(|observation| {
            std::iter::repeat_n(observation.bases.to_vec(), observation.num_obs as usize)
        })
        .collect();
    folded.sort();

    assert_eq!(
        folded, expected_bases,
        "the cap kept the reads whose bases are {folded:?}; the two smallest sampling \
         keys belong to {survivors:?}, whose bases are {expected_bases:?}"
    );
    let total: u32 = rec
        .observations
        .iter()
        .map(|observation| observation.num_obs)
        .sum();
    assert_eq!(total, 2, "exactly the cap's worth folded, and no more");
}

#[test]
fn column_depth_cap_uses_indel_cap_when_any_indel_event_present() {
    // Four SNP-only reads + one indel-bearing read, all anchored
    // at pos 1. At column pos 1 the indel-bearing read contributes
    // an Insertion event — the column flips to "indel column" and
    // the tighter indel cap applies. At pos 2 onward only Match
    // events remain, so the (much higher) SNP cap applies and
    // does not fire.
    let fa = MockFasta::new("AAAACGT");
    let mut reads: Vec<PreparedRead> = (0..4)
        .map(|i| snp_read(&format!("snp{i}"), 1, b"AAAACGT", &[30; 7]))
        .collect();
    // Indel read: 1M 2I 5M starting at pos 1; anchor of the
    // insertion is pos 1.
    let indel = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(5)],
        seq: b"AXXAAACG".to_vec(),
        bq_baq: vec![30; 8],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("indel"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    reads.push(indel);

    let cfg = WalkerConfig {
        max_snp_column_depth: 99,  // far above 5; SNP-only cols don't fire
        max_indel_column_depth: 2, // below 5; indel col at pos 1 fires
        ..WalkerConfig::default()
    };
    let (_records, summary) = drive_walker_with_config(reads, fa, &cfg);

    // Exactly one column carried an indel event (pos 1), and that
    // column had 5 contributors > indel cap of 2 → one truncation.
    assert_eq!(
        summary.column_depth_truncations, 1,
        "indel cap should fire exactly once at the indel-anchor column",
    );
}

#[test]
fn column_depth_cap_does_not_fire_below_threshold() {
    // Two SNP-only reads, default config (caps 8000 and 250).
    // Far below either cap → no truncation, every contributor
    // folds.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 1, b"ACGTA", &[30; 5]);
    let (records, summary) = drive_walker_with_summary(vec![r1, r2], fa);

    assert_eq!(summary.column_depth_truncations, 0);
    for rec in &records {
        assert_eq!(
            rec.reference_observation().num_obs,
            2,
            "pos {}: both reads should fold (no truncation under default cap)",
            rec.anchor(),
        );
    }
}

// --- G1: adaptor-region per-base filter, walker integration ---------

#[test]
fn g1_walker_drops_match_observations_past_adaptor_boundary() {
    // Reference: ACGTACGT (positions 1..9 1-based, 8 bases).
    // Two reads on the same molecule, ancient-DNA shape:
    //   - r_fwd: forward-strand, M(8) at pos 1, REF on every base.
    //     Adaptor boundary at 5 (insert size = 4 < seq_len = 8) →
    //     positions 1..4 emit, positions 5..8 are dropped as adaptor.
    //   - r_rev: reverse-strand, M(8) at pos 1, REF on every base.
    //     Adaptor boundary at 4 (mate.start - 1) → positions 1..4
    //     are dropped, positions 5..8 emit.
    //
    // Net result at every position 1..8: exactly one read contributes
    // (the one whose strand has that position outside its adaptor).
    // Without G1, both reads would fold at every position and num_obs
    // would be 2 — this is the regression the filter prevents.
    let fa = MockFasta::new("ACGTACGT");
    let mut r_fwd = snp_read("pair", 1, b"ACGTACGT", &[30; 8]);
    r_fwd.is_reverse_strand = false;
    r_fwd.mate_role = MateRole::FirstOfPair;
    r_fwd.adaptor_boundary = Some(5);
    let mut r_rev = snp_read("pair", 1, b"ACGTACGT", &[30; 8]);
    r_rev.is_reverse_strand = true;
    r_rev.mate_role = MateRole::SecondOfPair;
    r_rev.adaptor_boundary = Some(4);

    let records = drive_walker(vec![r_fwd, r_rev], fa);
    assert_eq!(records.len(), 8, "every covered position emits one record");
    for rec in &records {
        // Each position is covered by exactly one of the two mates
        // after the adaptor filter applies. Without G1, mate-overlap
        // resolution would still cap to 1 (sum-and-cap on agreement
        // applies), but that path mutates BQ; here the filter cleanly
        // removes the adaptor base from the contributor list before
        // overlap resolution sees it.
        assert_eq!(
            rec.reference_observation().num_obs,
            1,
            "pos {}: exactly one mate is outside adaptor at this position",
            rec.anchor(),
        );
    }
}

// --- Error-variant coverage --------------------------------------
//
// M14–M17 of `ia/reviews/pileup_2026-05-11.md`: every `WalkerError`
// variant must have a regression test pinning the exact variant
// returned. Without these, a swallowed or mis-mapped error could
// silently degrade to `None` / `Ok(())` (data loss) or to the wrong
// variant (operator triage misled).
//
// (The pre-iterator `run_returns_channel_closed_when_receiver_dropped_mid_stream`
// regression has been removed: the pull-shape walker no longer
// owns or sends through a channel, and the `ChannelClosed` variant
// is gone.)

/// Iterate the walker until it yields its first error and return
/// it. Panics if the walker exhausts without erroring.
fn first_walker_error(reads: Vec<PreparedRead>, fa: MockFasta) -> super::WalkerError {
    run(reads, &fa, &WalkerConfig::default())
        .find_map(|r| r.err())
        .expect("walker did not surface any error")
}

#[test]
fn zero_ref_span_input_is_a_hard_error() {
    // M15. A malformed PreparedRead with alignment_end <
    // alignment_start must hard-error on admission.
    let fa = MockFasta::new("ACGT");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 3,
        alignment_end: 2,
        cigar: vec![CigarOp::Match(0)],
        seq: vec![],
        bq_baq: vec![],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("zero"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    assert!(
        matches!(err, super::WalkerError::ZeroRefSpan { .. }),
        "got {err:?}",
    );
}

#[test]
fn open_record_widening_past_max_record_span_errors() {
    // M16 widen path. A single read with a deletion that pushes the
    // open record's footprint past MAX_RECORD_SPAN must error.
    let ref_len = (super::DEFAULT_MAX_RECORD_SPAN as usize) + 50;
    let fa = MockFasta::new(&"A".repeat(ref_len));
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: ref_len as u32,
        cigar: vec![
            CigarOp::Match(1),
            CigarOp::Deletion(super::DEFAULT_MAX_RECORD_SPAN + 1),
            CigarOp::Match(1),
        ],
        seq: b"AA".to_vec(),
        bq_baq: vec![30, 30],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("wide"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::RecordTooWide { cap, .. } => {
            assert_eq!(cap, super::DEFAULT_MAX_RECORD_SPAN);
        }
        other => panic!("expected RecordTooWide, got {other:?}"),
    }
}

// M21: malformed-PreparedRead validation. The walker must reject
// the upstream contract violation as a typed error, not panic on
// out-of-bounds indexing.

#[test]
fn admit_rejects_seq_shorter_than_cigar_consumes() {
    // CIGAR consumes 5 read bases (5M), but seq has only 3.
    let fa = MockFasta::new("ACGTA");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACG".to_vec(),
        bq_baq: vec![30; 3],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("short"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::MalformedRead { reason, .. } => {
            assert!(
                reason.contains("CIGAR consumes 5 read bases but seq.len = 3"),
                "got reason: {reason}",
            );
        }
        other => panic!("expected MalformedRead, got {other:?}"),
    }
}

#[test]
fn admit_rejects_seq_bq_length_mismatch() {
    let fa = MockFasta::new("ACGTA");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACGTA".to_vec(),
        bq_baq: vec![30, 30, 30], // 3 instead of 5
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("bq_short"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::MalformedRead { reason, .. } => {
            assert!(
                reason.contains("seq.len (5) != bq_baq.len (3)"),
                "got reason: {reason}",
            );
        }
        other => panic!("expected MalformedRead, got {other:?}"),
    }
}

#[test]
fn admit_rejects_cigar_consuming_more_read_bases_than_seq_provides() {
    // CIGAR = 2M + 4I + 2M. Consumes 8 read bases; seq has 5.
    // The cursor would index seq[2..6] for the Insertion, which
    // is OOB. Admit-time check must catch it first.
    let fa = MockFasta::new("ACGTA");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 4,
        cigar: vec![CigarOp::Match(2), CigarOp::Insertion(4), CigarOp::Match(2)],
        seq: b"ACGTA".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("cigar_long"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::MalformedRead { reason, .. } => {
            assert!(
                reason.contains("CIGAR consumes 8 read bases but seq.len = 5"),
                "got reason: {reason}",
            );
        }
        other => panic!("expected MalformedRead, got {other:?}"),
    }
}

#[test]
fn fasta_fetch_failure_propagates_as_walker_error_fasta() {
    // M17. Reference shorter than the read — MockFasta::fetch
    // returns UnexpectedEof when the open_new fetch runs off the
    // chromosome's end. Must surface as WalkerError::Fasta with
    // the correct locus.
    let fa = MockFasta::new("AC");
    let r = snp_read("r", 1, b"ACGT", &[30; 4]);
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::Fasta {
            chrom_id,
            start,
            start_plus_len,
            ..
        } => {
            assert_eq!(chrom_id, 0);
            assert!(start_plus_len > start);
        }
        other => panic!("expected Fasta, got {other:?}"),
    }
}

#[test]
fn walker_iterator_returns_none_after_yielding_error() {
    // The iterator's terminal-on-first-error contract: once `next()`
    // returns `Err`, every subsequent call returns `None`. Mirrors
    // the previous push-shape behaviour, where `run` returned `Err`
    // and stopped emitting at the same moment.
    let fa = MockFasta::new("AC"); // too short to satisfy a 4-base read.
    let r = snp_read("r", 1, b"ACGT", &[30; 4]);
    let mut walker = run(vec![r], &fa, &WalkerConfig::default());
    // The error may be preceded by zero or more successful records
    // (depending on the order of fasta fetch vs. walker_pos). Drive
    // until we see the first Err, then assert exhaustion.
    let mut saw_error = false;
    for item in &mut walker {
        if item.is_err() {
            saw_error = true;
            break;
        }
    }
    assert!(saw_error, "expected the walker to surface a Fasta error");
    assert!(
        walker.next().is_none(),
        "walker must return None after its first error",
    );
}

// ---------------------------------------------------------------------
// PreparedRead::length — internal length-invariant checks
// ---------------------------------------------------------------------

#[test]
fn prepared_read_length_returns_validated_length_on_well_formed_read() {
    // 4-bp Match: seq.len = bq_baq.len = cigar-consumed = 4.
    let r = snp_read("r", 1, b"ACGT", &[30; 4]);
    assert_eq!(r.length(), Ok(4));
}

#[test]
fn prepared_read_length_with_insertion_includes_inserted_bases() {
    // CIGAR 2M1I2M consumes 5 read bases. seq and bq_baq must
    // match. Confirms the helper sums the I op (read-consuming)
    // and not just the M ops.
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 4,
        cigar: vec![CigarOp::Match(2), CigarOp::Insertion(1), CigarOp::Match(2)],
        seq: b"ACXGT".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    assert_eq!(r.length(), Ok(5));
}

#[test]
fn prepared_read_length_with_deletion_excludes_deleted_bases() {
    // CIGAR 2M2D2M consumes 4 read bases (D op is reference-only).
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(2), CigarOp::Deletion(2), CigarOp::Match(2)],
        seq: b"ACGT".to_vec(),
        bq_baq: vec![30; 4],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    assert_eq!(r.length(), Ok(4));
}

#[test]
fn prepared_read_length_seq_bq_mismatch_returns_typed_error() {
    let mut r = snp_read("r", 1, b"ACGT", &[30; 4]);
    r.bq_baq.pop(); // bq_baq now length 3, seq length 4.
    assert_eq!(
        r.length(),
        Err(super::ReadLengthError::SeqBqMismatch {
            seq_len: 4,
            bq_baq_len: 3,
        })
    );
}

#[test]
fn prepared_read_length_cigar_seq_mismatch_returns_typed_error() {
    // seq.len = bq_baq.len = 4 but cigar consumes 3 read bases.
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACGT".to_vec(),
        bq_baq: vec![30; 4],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    assert_eq!(
        r.length(),
        Err(super::ReadLengthError::CigarSeqMismatch {
            cigar_consumed: 3,
            seq_len: 4,
        })
    );
}

#[test]
fn prepared_read_length_checks_seq_bq_before_cigar() {
    // Both invariants violated; the seq/bq check runs first so
    // SeqBqMismatch wins. Locks in the deterministic ordering.
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACGT".to_vec(), // length 4
        bq_baq: vec![30; 3],   // length 3 ≠ seq.len → SeqBqMismatch
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),

        mapq: 60,
    };
    assert!(matches!(
        r.length(),
        Err(super::ReadLengthError::SeqBqMismatch { .. })
    ));
}

// ---------------------------------------------------------------------
// The spliced fixture (D6) — the failure this whole change exists for
// ---------------------------------------------------------------------

/// A 60-base reference, `ACGT` fifteen times. The pattern only has to make the reads'
/// matched bases distinguishable from each other; every read below agrees with it, so the
/// fixture's subject is the *shape* of a witness and never an allele.
const SPLICED_REF: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

/// The reference bases at 1-based inclusive `from..=to`.
fn spliced_ref_bases(from: usize, to: usize) -> Vec<u8> {
    SPLICED_REF.as_bytes()[from - 1..to].to_vec()
}

/// Two reads over `SPLICED_REF`, and the deletion's length is the knob.
///
/// - **the spliced read**, `3M 15N 3M` from 28: exon 1 at 28–30, a 15-base intron at 31–45,
///   exon 2 at 46–48. A `Skip` emits no event (`cigar_cursor.rs`), so the read witnesses
///   **six** positions in two runs and nothing between them;
/// - **the deletion read**, `3M <del>D 3M` from 26: matches at 26–28, then `del` deleted
///   positions from 29. It anchors a record at **28** — the base before the deletion — whose
///   footprint the deletion widens to `28 ..= 28 + del`.
///
/// The deletion read is why the spliced read's hole is *inside a record at all*: an intron
/// cannot widen a record on its own, so without an indel allele spanning it the two exons
/// would simply be separate records and there would be no hole to represent (spec §8).
fn spliced_and_deleting_reads(deletion_len: u32) -> Vec<PreparedRead> {
    let deletion_len_usize = deletion_len as usize;
    let spliced = PreparedRead {
        chrom_id: 0,
        alignment_start: 28,
        alignment_end: 48,
        cigar: vec![CigarOp::Match(3), CigarOp::Skip(15), CigarOp::Match(3)],
        seq: [spliced_ref_bases(28, 30), spliced_ref_bases(46, 48)].concat(),
        bq_baq: vec![30; 6],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("spliced"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),
        mapq: 60,
    };
    let deleting = PreparedRead {
        chrom_id: 0,
        alignment_start: 26,
        alignment_end: 31 + deletion_len,
        cigar: vec![
            CigarOp::Match(3),
            CigarOp::Deletion(deletion_len),
            CigarOp::Match(3),
        ],
        seq: [
            spliced_ref_bases(26, 28),
            spliced_ref_bases(29 + deletion_len_usize, 31 + deletion_len_usize),
        ]
        .concat(),
        bq_baq: vec![30; 6],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("deleting"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
        read_group: ReadGroupId(0),
        mapq: 60,
    };
    vec![deleting, spliced]
}

/// The record the deletion anchors at 28, and the spliced read's observation in it (if any).
fn spliced_reads_observation(deletion_len: u32) -> (Locus, Option<SequenceObservation>) {
    let records = drive_walker(
        spliced_and_deleting_reads(deletion_len),
        MockFasta::new(SPLICED_REF),
    );
    let record = records
        .iter()
        .find(|record| record.anchor() == 28)
        .expect("the deletion anchors a record at 28")
        .clone();
    // The spliced read's bases are its two exons concatenated; the deletion read's observation over
    // the same footprint is the anchor base alone, so the two cannot be confused.
    let observation = record
        .observations
        .iter()
        .find(|observation| observation.bases.len() > 1)
        .cloned();
    (record, observation)
}

/// **A read blind in the middle of a record is recorded, with both of its runs** — the
/// regression anchor for the whole change, and the only fixture here drawn from a real
/// failure rather than constructed to exercise a branch (spec §7, §8; arch §6).
///
/// A 20-base deletion widens the record at 28 to `28 ..= 48`, twenty-one positions. The
/// spliced read witnessed **six** of them — three in each exon — and, before C3, was
/// **absent from the record entirely**: `apply_events_into` answered `None` for a
/// non-contiguous witness and the drop path fired once per position the record was affected
/// at. The evidence was not merely mis-described, it was discarded.
///
/// The assertions are the two halves of that: the read is *there*, and its witness says
/// `[(0,3), (18,21)]` — two runs — rather than a span that swallows the fifteen positions it
/// never saw.
#[test]
fn a_spliced_read_across_a_widened_record_is_recorded_with_both_of_its_runs() {
    let (record, observation) = spliced_reads_observation(20);
    assert_eq!(
        record.footprint_len(),
        21,
        "28 ..= 48, the anchor plus 20 deleted"
    );

    let observation = observation.expect(
        "the spliced read's observation — before C3 the walk dropped it, and this fixture is \
         what says it no longer does",
    );
    assert_eq!(
        &*observation.bases,
        &[spliced_ref_bases(28, 30), spliced_ref_bases(46, 48)].concat()[..],
        "the bases are the two exons and nothing from the intron — the fold copies no \
         reference base the read did not sequence",
    );
    assert_eq!(
        observation.read_witness,
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 3), (18, 21)])
                .expect("two runs"),
        },
        "exon 1 at locus offsets 0..3, exon 2 at 18..21, and the intron's fifteen positions \
         are a hole rather than a stretch the read is credited with",
    );
    assert_eq!(
        record.reads_without_observation, 0,
        "a holed witness is evidence, not a read that witnessed nothing — the counter this \
         read used to land in now means what its name says",
    );
}

/// **One base of footprint width decides it**, which is what makes the fixture above a
/// knife-edge rather than a decoration (spec §8).
///
/// The hole exists only where the record's footprint reaches *past* the intron into exon 2:
///
/// - a **17**-base deletion widens the record to `28 ..= 45`, one position short of exon 2 at
///   46, so the spliced read's witness inside it is exon 1 alone — one run, and the old
///   representation described it perfectly;
/// - **18** reaches 46, the first base of exon 2, and the same read is holed at once.
///
/// So the change earns its keep on a single deleted base, and a test that only ever built the
/// 20-base case could not tell the two readings apart.
#[test]
fn one_more_deleted_base_is_what_turns_the_spliced_read_into_a_holed_one() {
    let (short, short_observation) = spliced_reads_observation(17);
    assert_eq!(short.footprint_len(), 18, "28 ..= 45");
    assert_eq!(
        short_observation
            .expect("exon 1 alone is still an observation")
            .read_witness,
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 3)]).expect("one run"),
        },
        "the footprint stops before exon 2, so the read witnessed one contiguous stretch",
    );

    let (wide, wide_observation) = spliced_reads_observation(18);
    assert_eq!(wide.footprint_len(), 19, "28 ..= 46");
    assert_eq!(
        wide_observation
            .expect("the holed observation")
            .read_witness,
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 3), (18, 19)])
                .expect("two runs"),
        },
        "one more deleted base reaches the first position of exon 2, and the witness is two \
         runs with a fifteen-position hole between them",
    );
}

/// **The walk's own counter sees the spliced read, which is what makes the number obtainable
/// from a real BAM** (owner, 2026-07-31).
///
/// The same two counts exist on the parity census, but that lives behind `#[cfg(test)]` and only
/// measures loci where *production's* walker also produced a record — so it can never answer "how
/// often does a read see a locus in two pieces on real RNA-seq", which is the one open
/// measurement this representation was built for (spec §8). These ride the walk's `RunSummary`,
/// so pointing the generic dump at any BAM reports them.
///
/// The numbers are the fixture's, and they are exact rather than "greater than zero": **one**
/// holed read, blind over **15** positions — the intron at 31–45, which is the gap between the
/// two exons the read did see.
#[test]
fn the_walks_summary_counts_the_spliced_read_and_the_positions_it_was_blind_over() {
    let (_records, summary) =
        drive_walker_with_summary(spliced_and_deleting_reads(20), MockFasta::new(SPLICED_REF));
    assert_eq!(
        (summary.reads_with_holed_witness, summary.hole_positions),
        (1, 15),
        "one read saw the locus in two pieces and was blind over the 15-base intron between \
         them; the deletion read beside it witnessed the whole footprint and is not holed",
    );

    // The other side of the knife-edge: at 17 deleted bases the footprint stops before exon 2,
    // so the same read witnesses one contiguous piece and the counter must stay at zero. Without
    // this half, a counter that simply counted every partial read would pass the assertion above.
    let (_records, unholed) =
        drive_walker_with_summary(spliced_and_deleting_reads(17), MockFasta::new(SPLICED_REF));
    assert_eq!(
        (unholed.reads_with_holed_witness, unholed.hole_positions),
        (0, 0),
        "a partial witness of one run is not a hole",
    );
}
