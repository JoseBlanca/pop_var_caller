//! **Byte-parity of ng's prepared read against production's — the port anchor.**
//!
//! `LeftAlignPreparer` is a port of production's `process_read` fold, minus the stages ng moved
//! elsewhere: `G2`/`F1` are step-1 filters and BAQ is deferred sine die, so what remains is `F3`,
//! the left-alignment. This test is what makes that claim checkable.
//!
//! # What parity means here, exactly
//!
//! The oracle is production's `process_read(read, None, …)` — its `--no-baq` arm — run in-process
//! with **`max_read_mismatch_fraction: None`**. That setting is not cosmetic: ng moved `F1` to step
//! 1, so leaving it on would make production *drop* reads ng keeps, and the two keep-sets would
//! diverge for reasons that have nothing to do with preparation.
//!
//! Parity is asserted **field by field** rather than with `assert_eq!`, because `PreparedRead`
//! derives no `PartialEq` — and because a field-by-field failure names the field that moved, which
//! is the difference between "parity broke" and "`alignment_end` broke".
//!
//! # The two references, and why the second must NOT match
//!
//! - On an **all-uppercase** reference the two must be byte-identical. That is the port anchor.
//! - On a **soft-masked** (lowercase) reference they must **differ**, and the count of differing
//!   reads is the measurement this whole divergence was argued for. Production hands `F3` raw
//!   reference bytes while read bases were uppercased at decode, so its shift loop compares `A`
//!   against `a`, finds them unequal, and never advances: **production's left-alignment does
//!   nothing at all inside a soft-masked region** — which is to say, inside repeats, which is where
//!   indels are ambiguous and left-alignment is the entire point. ng fetches the uppercased view
//!   and shifts normally (spec §6, §9).
//!
//! # Why this lives in `src/ng/` and reads production
//!
//! It is **ng's test**: its subject is ng's preparer and production is only the yardstick. It is
//! `#[cfg(test)]`, so shipping ng code carries no test-only dependency on production — the same
//! shape as [`delimit_parity`](crate::ng::alignment::delimit_parity).

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use noodles_fasta as fasta;
use tempfile::TempDir;

use super::ReadPreparer;
use super::left_align::LeftAlignPreparer;
use crate::bam::alignment_input::{MappedRead, cigar_ref_span, read_exceeds_mismatch_fraction};
use crate::fasta::{ContigEntry, ContigList};
use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::ref_seq::{RefSeq, ResidentRefSeq};
use crate::ng::types::ContigId;
use crate::ng::types::ReadGroupId;
use crate::pileup::per_sample::read_processor::{
    RawContigRefCache, ReadOutcome, ReadProcessingConfig, process_read,
};
use crate::pileup::walker::{CigarOp, PreparedRead};

/// One contig of a fixture reference: its name and its bases, verbatim.
struct FixtureContig {
    name: &'static str,
    bases: &'static str,
}

/// Write a FASTA and its `.fai` holding `contigs` exactly as given — **case preserved**, which is
/// the whole point of the soft-masked arm.
///
/// Each contig is one header line and one sequence line, so the `.fai` arithmetic stays trivial:
/// `line_bases` is the contig length and `line_width` is one more, for the newline.
fn write_fasta(contigs: &[FixtureContig]) -> (TempDir, PathBuf) {
    let dir = tempfile::tempdir().expect("a temp dir");
    let fasta_path = dir.path().join("ref.fa");
    let fai_path = dir.path().join("ref.fa.fai");

    let mut fa = File::create(&fasta_path).expect("create fasta");
    let mut fai = File::create(&fai_path).expect("create fai");
    let mut offset: u64 = 0;
    for contig in contigs {
        let header = format!(">{}\n", contig.name);
        fa.write_all(header.as_bytes()).expect("write header");
        offset += header.len() as u64;

        fa.write_all(contig.bases.as_bytes()).expect("write bases");
        fa.write_all(b"\n").expect("write newline");

        let length = contig.bases.len() as u64;
        writeln!(
            fai,
            "{}\t{}\t{}\t{}\t{}",
            contig.name,
            length,
            offset,
            length,
            length + 1
        )
        .expect("write fai entry");
        offset += length + 1;
    }
    (dir, fasta_path)
}

/// A `Repository` over a fixture FASTA — what both sides read their bases through.
fn repository(fasta_path: &Path) -> fasta::Repository {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .expect("the fixture fasta has a sibling .fai");
    let adapter = fasta::repository::adapters::IndexedReader::new(reader);
    fasta::Repository::new(adapter)
}

fn contig_list(contigs: &[FixtureContig]) -> ContigList {
    ContigList {
        entries: contigs
            .iter()
            .map(|contig| ContigEntry {
                name: contig.name.to_string(),
                length: contig.bases.len() as u64,
                md5: None,
            })
            .collect(),
    }
}

/// Production's answer: the `--no-baq` arm of `process_read`, with `F1` disabled.
fn production_prepared(
    read: MappedRead,
    fasta_path: &Path,
    contigs: &ContigList,
) -> Option<PreparedRead> {
    let mut raw_ref = RawContigRefCache::new(repository(fasta_path), contigs.clone());
    let config = ReadProcessingConfig {
        // Disabled deliberately: ng runs the mismatch-fraction filter in step 1, so leaving it on
        // here would drop reads ng keeps and the keep-sets would differ for unrelated reasons.
        max_read_mismatch_fraction: None,
        mismatch_bq_floor: 0,
    };
    match process_read(read, None, &mut raw_ref, &config) {
        ReadOutcome::Prepared(prepared) => Some(prepared),
        ReadOutcome::Dropped(_) => None,
    }
}

/// ng's answer.
///
/// The fixture is built as production's `MappedRead`, because the other half of
/// this comparison is production's own path and both sides must start from one
/// record. Converting here is the whole of the difference: ng's read is that
/// read plus its read group, and read preparation reads none of it.
fn ng_prepared(read: MappedRead, fasta_path: &Path, contigs: &ContigList) -> PreparedRead {
    let reference = ResidentRefSeq::new(repository(fasta_path), contigs.clone());
    let preparer = LeftAlignPreparer::with_default_normalizer(reference);
    let mut scratch = <LeftAlignPreparer<ResidentRefSeq> as ReadPreparer>::Scratch::default();
    preparer
        .prepare_read(to_aligned_read(read), &mut scratch)
        .expect("the fixture windows are all in range")
        .expect("v1 never declines a read")
}

/// The same record as ng's read type, carrying a placeholder read group that
/// nothing on this path looks at.
fn to_aligned_read(read: MappedRead) -> AlignedRead {
    AlignedRead {
        qname: read.qname,
        flag: read.flag,
        ref_id: read.ref_id,
        pos: read.pos,
        mapq: read.mapq,
        cigar: read.cigar,
        seq: read.seq,
        qual: read.qual,
        mate_ref_id: read.mate_ref_id,
        mate_pos: read.mate_pos,
        adaptor_boundary: read.adaptor_boundary,
        read_group: ReadGroupId(0),
    }
}

/// Compare every field of a `PreparedRead`, naming the one that moved.
///
/// `PreparedRead` derives no `PartialEq`, and writing one here would be worse than this: an
/// equality that silently ignored a field it did not know about is exactly the failure a parity
/// test exists to catch. Listing the fields means a field added to production forces this to be
/// updated rather than passing vacuously.
#[track_caller]
fn assert_same_prepared_read(ours: &PreparedRead, production: &PreparedRead, context: &str) {
    assert_eq!(ours.chrom_id, production.chrom_id, "chrom_id [{context}]");
    assert_eq!(
        ours.alignment_start, production.alignment_start,
        "alignment_start [{context}]"
    );
    assert_eq!(
        ours.alignment_end, production.alignment_end,
        "alignment_end [{context}]"
    );
    assert_eq!(ours.cigar, production.cigar, "cigar [{context}]");
    assert_eq!(ours.seq, production.seq, "seq [{context}]");
    assert_eq!(ours.bq_baq, production.bq_baq, "bq_baq [{context}]");
    assert_eq!(
        ours.mq_log_err, production.mq_log_err,
        "mq_log_err [{context}]"
    );
    assert_eq!(ours.mapq, production.mapq, "mapq [{context}]");
    assert_eq!(
        ours.is_reverse_strand, production.is_reverse_strand,
        "is_reverse_strand [{context}]"
    );
    assert_eq!(ours.qname, production.qname, "qname [{context}]");
    assert_eq!(
        ours.mate_role, production.mate_role,
        "mate_role [{context}]"
    );
    assert_eq!(
        ours.adaptor_boundary, production.adaptor_boundary,
        "adaptor_boundary [{context}]"
    );
}

/// A read at 1-based `pos` on contig 0.
fn read_at(qname: &str, pos: u64, cigar: Vec<CigarOp>, seq: &[u8]) -> MappedRead {
    MappedRead {
        qname: qname.as_bytes().to_vec(),
        flag: 0,
        ref_id: 0,
        pos,
        mapq: 60,
        cigar,
        seq: seq.to_vec(),
        qual: (0..seq.len()).map(|i| 20 + (i as u8 % 15)).collect(),
        mate_ref_id: None,
        mate_pos: None,
        adaptor_boundary: None,
        source_file_index: 0,
    }
}

/// The fixture reference, shared by every case; each read is placed so its window sits inside it.
///
/// ```text
///  pos  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
///  base C C T A A A A G G  T  T  T  T  T  A  C  G  C  A  T  G  C  G  C  G
///                \_____/       \________/
///                A-run 4–7       T-run 10–14
/// ```
///
/// The two homopolymer runs are where an indel has somewhere to slide; the unique sequence from 16
/// on is where it has nowhere. **Every read below was checked base by base** against this: its
/// read-consuming CIGAR ops sum to `seq.len()`, and its `M` blocks match the reference. That is not
/// pedantry — a CIGAR whose read length disagrees with `seq` makes `left_align_cigar` bail *safe*
/// and hand the CIGAR back untouched, so both sides would "agree" while normalizing nothing.
const REFERENCE: &[FixtureContig] = &[FixtureContig {
    name: "chr1",
    bases: "CCTAAAAGGTTTTTACGCATGCGCG",
}];

/// Every read the parity fixture runs, with the case each one covers and, where it shifts, the
/// spelling it should end up with.
fn fixture_reads() -> Vec<(&'static str, MappedRead)> {
    vec![
        (
            "no indel — the pass-through path, where nothing may move",
            // ref 1–10 = CCTAAAAGGT
            read_at("plain", 1, vec![CigarOp::Match(10)], b"CCTAAAAGGT"),
        ),
        (
            "deletion in the A-run, spelled at its right end — must move left (M4 D1 M5 → M1 D1 M8)",
            // ref 3–12 = TAAAAGGTTT, minus one A. M4=TAAA(3–6), D1=A(7), M5=GGTTT(8–12).
            read_at(
                "del_in_a_run",
                3,
                vec![CigarOp::Match(4), CigarOp::Deletion(1), CigarOp::Match(5)],
                b"TAAAGGTTT",
            ),
        ),
        (
            "insertion after the T-run — must move left (M6 I1 M1 → M1 I1 M6)",
            // ref 9–15 = GTTTTTA, plus one inserted T. M6=GTTTTT(9–14), I1=T, M1=A(15).
            read_at(
                "ins_in_t_run",
                9,
                vec![CigarOp::Match(6), CigarOp::Insertion(1), CigarOp::Match(1)],
                b"GTTTTTTA",
            ),
        ),
        (
            "the same deletion already leftmost — normalization must be a no-op",
            // M1=T(3), D1=A(4), M8=AAAGGTTT(5–12).
            read_at(
                "already_left",
                3,
                vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(8)],
                b"TAAAGGTTT",
            ),
        ),
        (
            "deletion in unique sequence — nowhere to slide",
            // ref 16–23 = CGCATGCG, minus the A at 19. M3=CGC(16–18), D1=A(19), M4=TGCG(20–23).
            read_at(
                "unique_context",
                16,
                vec![CigarOp::Match(3), CigarOp::Deletion(1), CigarOp::Match(4)],
                b"CGCTGCG",
            ),
        ),
        (
            "soft-clipped read with a shiftable deletion — clips consume read but not reference",
            // S2 (arbitrary bases) then the del_in_a_run case: 2 + 4 + 5 = 11 read bases.
            read_at(
                "soft_clipped",
                3,
                vec![
                    CigarOp::SoftClip(2),
                    CigarOp::Match(4),
                    CigarOp::Deletion(1),
                    CigarOp::Match(5),
                ],
                b"GGTAAAGGTTT",
            ),
        ),
        (
            "reverse strand — the flag-derived fields must agree too",
            {
                let mut read = read_at(
                    "reverse",
                    3,
                    vec![CigarOp::Match(4), CigarOp::Deletion(1), CigarOp::Match(5)],
                    b"TAAAGGTTT",
                );
                read.flag = 0x10;
                read
            },
        ),
        ("paired, second of pair — mate_role must agree", {
            let mut read = read_at("paired", 1, vec![CigarOp::Match(10)], b"CCTAAAAGGT");
            read.flag = 0x1 | 0x80;
            read
        }),
    ]
}

/// **The port anchor.** On an all-uppercase reference, ng's prepared read must be byte-identical to
/// production's `--no-baq` output, field for field, for every case in the fixture.
#[test]
fn ng_matches_production_on_an_uppercase_reference() {
    let (_dir, fasta_path) = write_fasta(REFERENCE);
    let contigs = contig_list(REFERENCE);

    for (case, read) in fixture_reads() {
        let ours = ng_prepared(read.clone(), &fasta_path, &contigs);
        let theirs = production_prepared(read, &fasta_path, &contigs)
            .expect("production keeps every fixture read with F1 disabled");
        assert_same_prepared_read(&ours, &theirs, case);
    }
}

/// **The fixture would prove nothing if left-alignment never fired on it.** Byte-parity between two
/// implementations that both correctly do nothing is not evidence that either one shifts an indel.
///
/// This caught three real defects in the first draft of the fixture above — a read whose deletion
/// fell on a `G` rather than in the homopolymer it was named for, one whose `M` block did not match
/// the reference at all, and one whose CIGAR consumed nine read bases while its `seq` had ten,
/// which makes `left_align_cigar` bail safe and normalize nothing. All three passed byte-parity.
///
/// So the assertion is by **name**, not by count: a floor of "at least N moved" would still have
/// let two of those three through.
#[test]
fn the_parity_fixture_actually_exercises_left_alignment() {
    let (_dir, fasta_path) = write_fasta(REFERENCE);
    let contigs = contig_list(REFERENCE);

    let mut rewritten: Vec<String> = fixture_reads()
        .into_iter()
        .filter(|(_, read)| {
            let before = read.cigar.clone();
            ng_prepared(read.clone(), &fasta_path, &contigs).cigar != before
        })
        .map(|(_, read)| String::from_utf8(read.qname).expect("fixture qnames are ASCII"))
        .collect();
    rewritten.sort();

    // Named, not counted. `reverse` carries the same shiftable spelling as `del_in_a_run` — it is
    // the same deletion on the other strand — so it moves too; a bare count hid that.
    assert_eq!(
        rewritten,
        vec!["del_in_a_run", "ins_in_t_run", "reverse", "soft_clipped"],
        "exactly the shiftable reads should move, and only those"
    );
}

/// The same reference, **soft-masked**: identical bases, lowercased. This is what a repeat-masked
/// assembly looks like, and every analysis-set human reference ships this way.
const SOFT_MASKED_REFERENCE: &[FixtureContig] = &[FixtureContig {
    name: "chr1",
    bases: "cctaaaaggtttttacgcatgcgcg",
}];

/// **The divergence, measured — and the production defect, demonstrated.**
///
/// Production hands `F3` *raw* reference bytes, while read bases were uppercased at decode. So on a
/// soft-masked reference its shift loop compares `A` against `a`, finds them unequal, and never
/// advances: left-alignment silently does nothing, in exactly the sequence — repeats and
/// low-complexity — where indel placement is ambiguous and left-alignment is the whole point.
///
/// ng fetches the uppercased view, so it shifts here exactly as it does on unmasked sequence.
///
/// This test asserts three things, and the middle one is the finding:
///
/// 1. ng's answer on masked sequence is **identical to its own answer on unmasked sequence** — the
///    masking makes no difference to it, which is the property the decision was made for;
/// 2. **production's answer is the mapper's original CIGAR, unshifted** — it normalized nothing;
/// 3. therefore the two disagree, on exactly the reads that have somewhere to shift.
///
/// Spec §6, §9, §11 record this as a known production defect ng fixes; production is frozen, so the
/// number this test reports is the whole of what we can say about its size.
#[test]
fn production_left_alignment_does_nothing_on_a_soft_masked_reference() {
    let (_upper_dir, upper_path) = write_fasta(REFERENCE);
    let upper_contigs = contig_list(REFERENCE);
    let (_masked_dir, masked_path) = write_fasta(SOFT_MASKED_REFERENCE);
    let masked_contigs = contig_list(SOFT_MASKED_REFERENCE);

    let mut disagreed: Vec<String> = Vec::new();
    for (case, read) in fixture_reads() {
        let mapper_cigar = read.cigar.clone();
        let qname = String::from_utf8(read.qname.clone()).expect("ASCII qname");

        let ours_masked = ng_prepared(read.clone(), &masked_path, &masked_contigs);
        let ours_unmasked = ng_prepared(read.clone(), &upper_path, &upper_contigs);
        let theirs_masked = production_prepared(read, &masked_path, &masked_contigs)
            .expect("production keeps every fixture read with F1 disabled");

        // 1. Masking is invisible to ng.
        assert_eq!(
            ours_masked.cigar, ours_unmasked.cigar,
            "ng must give the same answer masked or not [{case}]"
        );
        // 2. Production normalized nothing: its output is the mapper's own spelling.
        assert_eq!(
            theirs_masked.cigar, mapper_cigar,
            "production's left-alignment should be inert on masked sequence [{case}]"
        );

        if ours_masked.cigar != theirs_masked.cigar {
            disagreed.push(qname);
        }
    }
    disagreed.sort();

    // 3. The two disagree on exactly the reads with somewhere to shift — the same four the
    //    uppercase fixture proves ng moves. On an unmasked reference these are byte-identical;
    //    the masking alone produces the difference.
    assert_eq!(
        disagreed,
        vec!["del_in_a_run", "ins_in_t_run", "reverse", "soft_clipped"],
        "every read with a shiftable indel should be mis-placed by production here"
    );
}

/// **Discharges the spec's step-1/step-2 ordering argument (§5).**
///
/// Production runs the mismatch-fraction filter `F1` *after* left-alignment; ng runs it in step 1,
/// *before*. The spec argues the verdict cannot move, because a gap only shifts across bases the
/// read and the reference share, so the match/mismatch tally is unchanged — but that is an argument,
/// and the invariant production carries is a `debug_assert` that compiles out of release *and*
/// counts a different quantity than `F1` does (raw, case-sensitive mismatches, versus `F1`'s
/// BQ-floored count over ATGC-comparable columns).
///
/// So this checks the thing itself: for every fixture read, at every threshold, `F1`'s verdict is
/// the same before and after left-alignment.
#[test]
fn left_alignment_cannot_change_the_mismatch_fraction_verdict() {
    let (_dir, fasta_path) = write_fasta(REFERENCE);
    let contigs = contig_list(REFERENCE);

    for (case, read) in fixture_reads() {
        let mut window = Vec::new();
        let span = u64::from(cigar_ref_span(&read.cigar));
        if span == 0 {
            continue;
        }
        ResidentRefSeq::new(repository(&fasta_path), contigs.clone())
            .fetch_into(ContigId(0), read.pos, span, &mut window)
            .expect("the fixture windows are in range");

        let after = ng_prepared(read.clone(), &fasta_path, &contigs);

        // Sweep the whole threshold range: a verdict that agreed only at one threshold would say
        // nothing, since the counts could still differ and simply fall the same side of it.
        for threshold in [0.0_f32, 0.05, 0.10, 0.25, 0.5, 1.0] {
            let before_verdict = read_exceeds_mismatch_fraction(
                &read.cigar,
                &read.seq,
                &read.qual,
                &window,
                0,
                threshold,
            );
            let after_verdict = read_exceeds_mismatch_fraction(
                &after.cigar,
                &after.seq,
                &after.bq_baq,
                &window,
                0,
                threshold,
            );
            assert_eq!(
                before_verdict, after_verdict,
                "F1's verdict moved across left-alignment at threshold {threshold} [{case}]"
            );
        }
    }
}

/// The shifted spellings themselves, so a regression that moved an indel to the *wrong* place —
/// rather than not moving it — is caught here and not only by the parity comparison.
#[test]
fn the_shiftable_reads_land_on_their_leftmost_spelling() {
    let (_dir, fasta_path) = write_fasta(REFERENCE);
    let contigs = contig_list(REFERENCE);
    let by_name = |name: &str| {
        let (_, read) = fixture_reads()
            .into_iter()
            .find(|(_, read)| read.qname == name.as_bytes())
            .expect("a fixture read by that name");
        ng_prepared(read, &fasta_path, &contigs)
    };

    // The deletion slides from the right end of the A-run to its left end.
    assert_eq!(
        by_name("del_in_a_run").cigar,
        vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(8)]
    );
    // The insertion slides from the right end of the T-run to just after the G that precedes it.
    assert_eq!(
        by_name("ins_in_t_run").cigar,
        vec![CigarOp::Match(1), CigarOp::Insertion(1), CigarOp::Match(6)]
    );
    // The soft clip is carried through untouched while the deletion inside it still shifts.
    assert_eq!(
        by_name("soft_clipped").cigar,
        vec![
            CigarOp::SoftClip(2),
            CigarOp::Match(1),
            CigarOp::Deletion(1),
            CigarOp::Match(8)
        ]
    );
    // And the two that must not move.
    assert_eq!(
        by_name("already_left").cigar,
        vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(8)]
    );
    assert_eq!(
        by_name("unique_context").cigar,
        vec![CigarOp::Match(3), CigarOp::Deletion(1), CigarOp::Match(4)]
    );
}
