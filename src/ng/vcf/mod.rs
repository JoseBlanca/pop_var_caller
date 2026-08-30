//! **The one file a run writes** — SNPs, indels and repeat tracts in a single VCF, where
//! production writes two (`doc/devel/ng/spec/vcf_output.md`).
//!
//! This module holds what a record *is* once its locus is gone, and — as later steps land —
//! how it becomes bytes. It computes no quality, fits no frequency and decides nothing about
//! which loci deserve to exist: the numbers arrive finished from the calling loop and the
//! quality stage, and *what is dropped* belongs to the emission steps
//! (`doc/devel/ng/spec/ng_proposal.md` steps 11a and 11b). What this module owns is the
//! spelling.
//!
//! **Why a carrier type at all, rather than encoding straight from a called locus.** Several
//! of the file's columns are read counts — `DP`, `AD`, and the pooled mapping qualities behind
//! `MQREF`/`MQALT`/`MQDIFF` — and no type downstream of the merge holds them:
//! [`LocusInference`](crate::ng::calling::LocusInference) carries the calls and the qualities,
//! while the reads themselves live in the merge's `SampleSupport` and are released with the
//! locus. So the counts are summed **in the worker, while the evidence is still in hand**, and
//! travel here. It is the same shape, for the same reason, as the quality module's nine-number
//! artifact summary (`doc/devel/ng/spec/calling_quality.md` §3.3): everything cohort-shaped is
//! summed away where it lives, and only scalars and short per-allele vectors cross the
//! boundary. The reference base an empty allele has to be padded with travels the same way
//! (§5), because it lies *outside* the record's own span and nothing downstream holds the
//! reference.
//!
//! **What the encoder still derives, and it is a closed list.** The alternative alleles are the
//! allele table past the reference; `AC` and `AN` are counted from the called genotypes; the
//! record's `INFO/DP` is the sum of the samples'; `AF` is [`VcfRecord::expected_copies`]
//! normalised; `MQDIFF` is a subtraction of two means this record already holds; and `REPCN` is
//! [`TractAnnotation::repeat_copies_of`] over the alleles a sample was called on. **No quality,
//! frequency or count is re-derived from anything but what this record carries.**
//!
//! **One record shape covers both kinds of locus** (spec §3). A repeat-tract record differs
//! from a SNP/indel one by carrying a [`TractAnnotation`] — which is what writes the `STR`
//! flag, `RU` and `PERIOD` — and nothing else. That is possible because the two paths already
//! agree on what a locus is: a reference span, a table of sequence alleles, and per-sample
//! support against that table. The reference partition keeps them off each other's bases
//! (`doc/devel/ng/spec/typed_regions.md` §2.3), so no two records ever describe the same
//! position.

use crate::ng::calling::quality::artifact_correction::ArtifactPenalties;
use crate::ng::types::{GenomeRegion, Genotype, Motif, Phred, Position};

/// The first base of a contig, in the 1-based inclusive coordinates
/// [`GenomeRegion`] uses — the one position at which an empty allele must be padded from the
/// right, because there is no base to its left (spec §5).
const FIRST_POSITION_OF_A_CONTIG: Position = Position(1);

/// **One record of the file, complete, with its locus released.**
///
/// Everything here is finished: the site quality is the corrected one, the genotypes are the
/// loop's, the counts were summed where the reads were.
///
/// **The site quality is the corrected value and is never recomputed here.** Production
/// recomputes its correction at VCF-encode time, and for sixteen days its emission gate
/// compared the engine's baseline while the corrected number went into the `QUAL` column — 40
/// false positives at 30× on GIAB HG002 and 64 at 50×, against 14 and 14 once both read one
/// function (`doc/devel/ng/spec/calling_quality.md` §3.5). ng carries one quality, written
/// upstream; this type is a place it rests, not a place it is derived.
///
/// **Most fields are private, and for one reason:** they are parallel to, or indexed by, the
/// allele table, or they carry an invariant [`Self::new`] checks once. Public fields would let
/// a consumer push an allele afterwards — leaving every sample's `AD`, the expected copies, the
/// pooled mapping qualities and any called genotype describing a table that no longer exists —
/// or write a backwards span past the check that refused it. The two that stay public,
/// [`Self::site_quality`] and [`Self::artifact_penalties`], are valid at every value their own
/// types admit and bind nothing else.
#[derive(Clone, PartialEq, Debug)]
pub struct VcfRecord {
    /// **How unlikely it is that no sample here carries a non-reference allele**, on the Phred
    /// scale, already corrected for artifact shape. Written, never recomputed — see the type's
    /// own comment for the defect that makes this a rule.
    pub site_quality: Phred,
    /// **What the two artifact tests took off the site quality** — `None` where they did not
    /// run, which is every repeat-tract record until the tract quality document exists
    /// (spec §6), and every locus the correction had nothing to test.
    ///
    /// Carried so the uncorrected quality stays recoverable as the sum of the three, wherever
    /// the correction's floor at zero did not bite. That is the annotate-and-defer convention
    /// the rest of the field follows, and it costs two floats a record.
    pub artifact_penalties: Option<ArtifactPenalties>,
    /// The reference span this record covers, before any anchor base is added — the anchor
    /// rule is a presentation step at encode time, not a property of the locus. Read through
    /// [`Self::region`].
    region: GenomeRegion,
    /// **The distinct alleles, the reference first**, each the full sequence over
    /// [`Self::region`]. Read through [`Self::alleles`].
    ///
    /// A record with the reference alone is legal and is what a locus the caller looked at and
    /// could not call comes to: `ALT` is written `.`, every sample is a no-call, and the
    /// [`FilterVerdict`] says why (spec §8).
    ///
    /// **An *alternative* may be empty** — a full-tract deletion — and that is what
    /// [`Self::padding_base`] and the encoder's anchor rule exist for. **The reference may
    /// not**: it is the bases the span spells, and an empty one would reach the `REF` column
    /// as an unparseable record. `CandidateAlleles::new` refuses it upstream for exactly this
    /// reason, naming this column; refusing it here too keeps the writer from being the weaker
    /// gate.
    alleles: Vec<Box<[u8]>>,
    /// **The cohort's expected copies of each allele**, parallel to [`Self::alleles`] — what
    /// `AF` is written from, after normalising over the called allele count.
    ///
    /// **Carried rather than recomputed downstream, and that is a rule with a reason.** This
    /// is the calling loop's converged fit; `AF` derived instead from the *called genotypes*
    /// is a different number, because a call has already discarded the uncertainty these
    /// counts hold (`LocusInference::cohort_expected_copies`, whose own doc says so). Read
    /// through [`Self::expected_copies`].
    expected_copies: Vec<f64>,
    /// **One entry per sample of the run, in the run's sample order** — dense, including the
    /// samples that had no coverage, because a sample column exists in every record whether or
    /// not the sample had anything to say.
    ///
    /// **Dense, where the merge's same-named field is sparse.** `CohortObservation::per_sample`
    /// holds only the covering samples, each naming its own index; by the time a record is
    /// built, every sample of the run needs a column. Read through [`Self::sample_columns`].
    sample_columns: Vec<SampleColumn>,
    /// **The cohort-pooled mapping qualities, one entry per allele**, parallel to
    /// [`Self::alleles`]. Read through [`Self::allele_mapq`].
    allele_mapq: Vec<MapqPool>,
    /// **The reference base an empty allele is padded with** — `Some` exactly where some
    /// alternative is empty, `None` otherwise (spec §5).
    ///
    /// It lies outside [`Self::region`] on either side, so it cannot be recovered from this
    /// record's own alleles, and nothing downstream of the worker holds the reference to fetch
    /// it from: `ReferenceInfo` carries contig geometry and digests, not bases, and a run
    /// driven from a `.fai` alone has no sequence at all. Resolving it here is what lets the
    /// encoder stay a pure function of the record — and what keeps it from falling back on
    /// production's invented `N`, which spec §5 does not port.
    padding_base: Option<PaddingBase>,
    /// The record's one `FILTER` value. Read through [`Self::filter`].
    filter: FilterVerdict,
    /// **`Some` exactly at a repeat tract**, and what makes this a tract record: it writes the
    /// `STR` flag, `RU` and `PERIOD` together, so no consumer can meet a record that claims one
    /// and not the others. Read through [`Self::repeat_tract`].
    repeat_tract: Option<TractAnnotation>,
}

impl VcfRecord {
    /// Assemble a record, checking every parity the allele table binds and every invariant the
    /// format's own columns rest on.
    ///
    /// # Panics
    ///
    /// On a record that cannot describe a locus, or that the encoder could not write: a
    /// backwards span; an empty allele table or an empty reference allele; a reference whose
    /// length disagrees with the span it is meant to spell; a per-allele vector of the wrong
    /// length; a called genotype naming an allele the table does not hold; an empty cohort; a
    /// missing or wrong-sided padding base; or pooled mapping qualities counting a different
    /// set of reads than the samples' `AD` does. Every one of them is a wiring defect in the
    /// stage that built the record — the merge and the loop agree on the table by construction
    /// — so they are assertions rather than a `Result` nobody could act on, which is
    /// [`LocusInference::new`](crate::ng::calling::LocusInference)'s choice for the same
    /// reason.
    #[allow(
        clippy::too_many_arguments,
        reason = "the ten are the record's own fields and the list is the point: each is \
                  owned by a different part of the format — the calling loop's genotypes and \
                  fitted copies, the quality stage's two numbers, the merge's counts, the \
                  reference's padding base, the emission step's filter — and the constructor \
                  exists so that a record cannot be assembled with one of them forgotten. \
                  Grouping them into sub-structs would invent a shape the file does not have, \
                  and every parameter has a distinct type, so no two can be transposed \
                  silently. The same ruling, for the same reason, as LocusInference::new"
    )]
    pub fn new(
        region: GenomeRegion,
        alleles: Vec<Box<[u8]>>,
        expected_copies: Vec<f64>,
        sample_columns: Vec<SampleColumn>,
        allele_mapq: Vec<MapqPool>,
        padding_base: Option<PaddingBase>,
        site_quality: Phred,
        artifact_penalties: Option<ArtifactPenalties>,
        filter: FilterVerdict,
        repeat_tract: Option<TractAnnotation>,
    ) -> Self {
        assert!(
            region.start <= region.end,
            "a record covers a stretch of reference, so its region cannot run backwards: \
             {region}"
        );
        assert!(
            !alleles.is_empty(),
            "every locus holds at least the reference allele, so a record with an empty \
             allele table has lost the table rather than having none — even a locus the \
             caller refused carries its reference and writes `ALT .`"
        );
        assert!(
            !alleles[0].is_empty(),
            "a record's reference allele is the bases its span spells — one for a SNP, \
             several for an indel, the whole tract for a repeat — and an empty one reaches \
             the REF column as an unparseable record. An empty *alternative* is a deletion \
             and is legal"
        );
        assert_eq!(
            alleles[0].len() as u64,
            region.len(),
            "the reference allele is the bases of the whole span, so a REF of {} bases over \
             a {}-base region describes a different stretch of reference than the record's \
             own POS claims",
            alleles[0].len(),
            region.len()
        );
        assert_eq!(
            expected_copies.len(),
            alleles.len(),
            "the cohort's expected allele copies run parallel to the allele table: one entry \
             per allele, reference first, so that AF names the alleles this record holds"
        );
        assert!(
            !sample_columns.is_empty(),
            "a record carries one column per sample of the run and a run has at least one \
             sample, so a record naming no sample has lost them"
        );
        assert_eq!(
            allele_mapq.len(),
            alleles.len(),
            "the pooled mapping qualities run parallel to the allele table: one entry per \
             allele, reference first, so that MQREF and MQALT name the alleles this record \
             actually holds"
        );
        for (index, column) in sample_columns.iter().enumerate() {
            assert_eq!(
                column.read_counts.allele_reads().len(),
                alleles.len(),
                "sample {index}'s allele read counts are {} entries against an allele table \
                 of {}: AD is written one entry per allele, reference first, so a tally of a \
                 different width was pooled against a different table",
                column.read_counts.allele_reads().len(),
                alleles.len()
            );
            if let Some(genotype) = column.call.genotype() {
                assert!(
                    genotype
                        .alleles()
                        .iter()
                        .all(|id| usize::from(id.get()) < alleles.len()),
                    "sample {index}'s call names an allele this record does not hold: the \
                     table has {} alleles, and a call carried across a renumbering without \
                     being remapped is how a genotype comes to point past it",
                    alleles.len()
                );
            }
        }
        for (allele, pool) in allele_mapq.iter().enumerate() {
            let attributed: u64 = sample_columns
                .iter()
                .map(|column| u64::from(column.read_counts.allele_reads()[allele]))
                .sum();
            assert_eq!(
                pool.reads, attributed,
                "allele {allele}'s mapping qualities were pooled over {} reads while the \
                 samples' AD attributes {attributed} to it: MQREF and MQALT are means over \
                 the same reads AD counts, so two totals mean two different pools",
                pool.reads
            );
        }
        let some_allele_is_empty = alleles.iter().any(|allele| allele.is_empty());
        assert_eq!(
            padding_base.is_some(),
            some_allele_is_empty,
            "a record with an empty allele is written by padding every allele with a \
             reference base beside the span, and one without needs no such base — so a \
             padding base is carried exactly when some allele is empty, and this record {} \
             one while {} (spec §5)",
            if padding_base.is_some() {
                "carries"
            } else {
                "carries no"
            },
            if some_allele_is_empty {
                "an allele is empty"
            } else {
                "every allele spells bases"
            }
        );
        if let Some(base) = padding_base {
            let at_contig_start = region.start == FIRST_POSITION_OF_A_CONTIG;
            assert_eq!(
                matches!(base, PaddingBase::Right(_)),
                at_contig_start,
                "the padding base is the one to the left of the span, except where the span \
                 starts at the contig's first base and there is none — this record starts at \
                 {} and carries a {} base",
                region.start.get(),
                if matches!(base, PaddingBase::Right(_)) {
                    "right-hand"
                } else {
                    "left-hand"
                }
            );
        }
        Self {
            site_quality,
            artifact_penalties,
            region,
            alleles,
            expected_copies,
            sample_columns,
            allele_mapq,
            padding_base,
            filter,
            repeat_tract,
        }
    }

    /// The reference span, before any padding base is applied.
    #[inline]
    #[must_use]
    pub fn region(&self) -> GenomeRegion {
        self.region
    }

    /// **The allele table, the reference at index 0**, each entry the full sequence over
    /// [`Self::region`].
    #[inline]
    #[must_use]
    pub fn alleles(&self) -> &[Box<[u8]>] {
        &self.alleles
    }

    /// The reference allele's bases — the `REF` column, before any padding base.
    #[inline]
    #[must_use]
    pub fn reference(&self) -> &[u8] {
        // PANIC-FREE: `Self::new` refuses an empty allele table, so index 0 always exists.
        &self.alleles[0]
    }

    /// The alternative alleles, in table order — the `ALT` column. Empty at a record that
    /// established none, which is written `.`.
    #[inline]
    #[must_use]
    pub fn alternatives(&self) -> &[Box<[u8]>] {
        // PANIC-FREE: `Self::new` refuses an empty allele table, so this slice is in range.
        &self.alleles[1..]
    }

    /// The cohort's expected copies of each allele, parallel to [`Self::alleles`] — what `AF`
    /// is written from.
    #[inline]
    #[must_use]
    pub fn expected_copies(&self) -> &[f64] {
        &self.expected_copies
    }

    /// One column per sample of the run, in the run's sample order.
    #[inline]
    #[must_use]
    pub fn sample_columns(&self) -> &[SampleColumn] {
        &self.sample_columns
    }

    /// The cohort-pooled mapping qualities, parallel to [`Self::alleles`].
    #[inline]
    #[must_use]
    pub fn allele_mapq(&self) -> &[MapqPool] {
        &self.allele_mapq
    }

    /// The reference base an empty allele is padded with, and which side it came from.
    #[inline]
    #[must_use]
    pub fn padding_base(&self) -> Option<PaddingBase> {
        self.padding_base
    }

    /// The record's one `FILTER` value.
    #[inline]
    #[must_use]
    pub fn filter(&self) -> FilterVerdict {
        self.filter
    }

    /// **What marks this a repeat-tract record** — `None` at a SNP or indel.
    #[inline]
    #[must_use]
    pub fn repeat_tract(&self) -> Option<&TractAnnotation> {
        self.repeat_tract.as_ref()
    }

    /// Whether the `STR` flag is written. **Not a field**: the flag is the presence of the
    /// annotation, so it cannot disagree with the `RU` and `PERIOD` beside it (spec §3).
    #[inline]
    #[must_use]
    pub fn is_repeat_tract(&self) -> bool {
        self.repeat_tract.is_some()
    }
}

/// **The reference base an empty allele is padded with**, and the side it was taken from.
///
/// VCF cannot spell an empty allele, so a deletion that removes a whole span is written by
/// giving every allele of the record one flanking reference base. Ordinarily that is the base
/// to the **left**, and the record's `POS` moves one left with it. At the very first base of a
/// contig there is nothing to the left, and the base to the **right** is appended instead with
/// `POS` unmoved (spec §5).
///
/// **Production's repeat-tract writer invents the letter `N` in that second case** — a base the
/// reference does not contain, at an unshifted position. That is the one behaviour of the two
/// production writers this format deliberately does not port.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PaddingBase {
    /// The reference base immediately before the record's span. Every allele is prefixed with
    /// it and `POS` moves one base left.
    Left(u8),
    /// The reference base immediately after the record's span, for a span that starts at the
    /// contig's first base. Every allele is suffixed with it and `POS` does not move.
    Right(u8),
}

impl PaddingBase {
    /// The base itself, whichever side it came from.
    #[inline]
    #[must_use]
    pub fn base(self) -> u8 {
        match self {
            Self::Left(base) | Self::Right(base) => base,
        }
    }
}

/// **What one sample contributes to a record** — its call, and what its reads showed.
///
/// The two are independent, and that is the point. A sample can have reads and no call, and the
/// file says so: `./.` with the evidence still written beside it, because a sample no-called on
/// three reads and one no-called on none are different facts a filter can use (spec §7).
#[derive(Clone, PartialEq, Debug)]
pub struct SampleColumn {
    /// The genotype and its quality, or a no-call.
    pub call: SampleCall,
    /// What this sample's reads showed — `DP` and the per-allele `AD`.
    pub read_counts: SampleReadCounts,
}

/// **What the file says this sample is** — a genotype, or no call at all.
///
/// **This is the VCF's notion of a no-call, and it is deliberately wider than the calling
/// loop's.** [`SampleGenotypeCall::Missing`](crate::ng::calling::SampleGenotypeCall) means one
/// specific thing — candidate selection cut an allele this sample's own reads had earned, so
/// the locus is called over a set that cannot hold what the sample carries — and it arises on
/// the SNP/indel path only. The file has to spell `./.` in three situations (spec §7): that
/// one, a sample with no coverage at the locus, and (when emission adopts one) a sample below a
/// per-sample quality floor. A repeat-tract locus the caller refused writes every sample as a
/// no-call (spec §8), which the loop's enum refuses to represent at a tract at all.
///
/// So the two are not the same idea and one type cannot carry both without one of the two
/// documents being wrong. **The conversion is the worker-side assembly's** — notably what to do
/// with a sample the loop *called* from the prior alone, having seen no reads, which spec §7
/// says must not be force-called.
///
/// **A no-call carries no quality**, and an enum is what enforces it: there is no `GQ` to read,
/// so nothing can write a genotype quality beside a genotype that was never scored.
#[derive(Clone, PartialEq, Debug)]
pub enum SampleCall {
    /// The file names this sample's alleles.
    Called {
        /// Which alleles it carries, one per copy of its genome — indices into the record's
        /// allele table, ascending.
        genotype: Genotype,
        /// How sure the caller is of that genotype: the `GQ` column.
        genotype_quality: Phred,
    },
    /// The file writes `./.` — and `GQ` as missing, never as zero. The two mean different
    /// things, and conflating them turns an absent call into a merely uncertain one.
    NoCall,
}

impl SampleCall {
    /// The alleles called, or `None` at a no-call.
    #[inline]
    #[must_use]
    pub fn genotype(&self) -> Option<&Genotype> {
        match self {
            Self::Called { genotype, .. } => Some(genotype),
            Self::NoCall => None,
        }
    }

    /// The `GQ` value, or `None` at a no-call — which has none, rather than having zero.
    #[inline]
    #[must_use]
    pub fn genotype_quality(&self) -> Option<Phred> {
        match self {
            Self::Called {
                genotype_quality, ..
            } => Some(*genotype_quality),
            Self::NoCall => None,
        }
    }

    /// Whether this sample's `GT` is written `./.`.
    #[inline]
    #[must_use]
    pub fn is_no_call(&self) -> bool {
        matches!(self, Self::NoCall)
    }
}

/// **One sample's read counts at one record: the `DP` and `AD` columns.**
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SampleReadCounts {
    /// Reads whose observed sequence matched each allele, parallel to the record's allele
    /// table, reference first.
    allele_reads: Vec<u32>,
    /// Every read this sample observed at the locus, whether or not a written allele
    /// explains it.
    depth: u32,
}

impl SampleReadCounts {
    /// Build a sample's counts.
    ///
    /// # Panics
    ///
    /// When the depth is below the reads attributed to alleles. `DP` counts a superset of what
    /// `AD` splits up, so the difference is *how many of this sample's reads no written allele
    /// explains* — a negative difference is not a thin locus, it is two counts taken from
    /// different pools.
    //
    // No `#[must_use]`, matching this module's other two constructors and
    // `LocusInference::new`: a test that exercises the refusal calls it for the panic alone,
    // and the attribute would make that call a lint error rather than a passing test.
    pub fn new(allele_reads: Vec<u32>, depth: u32) -> Self {
        let attributed = Self::attributed_reads(&allele_reads);
        assert!(
            u64::from(depth) >= attributed,
            "a sample's depth is every read it showed at the locus and its allele counts \
             split up part of that, so a depth of {depth} under {attributed} attributed reads \
             means the two were pooled over different sets of reads"
        );
        Self {
            allele_reads,
            depth,
        }
    }

    /// The reads attributed to written alleles — the total `DP` has to cover.
    ///
    /// Summed in `u64` because it runs before [`Self::new`]'s check, where the slice is still
    /// unvalidated and its total may exceed a `u32`.
    fn attributed_reads(allele_reads: &[u32]) -> u64 {
        allele_reads.iter().map(|reads| u64::from(*reads)).sum()
    }

    /// The `AD` column: reads matching each allele, reference first.
    #[inline]
    #[must_use]
    pub fn allele_reads(&self) -> &[u32] {
        &self.allele_reads
    }

    /// The `DP` column: every read the sample observed here.
    #[inline]
    #[must_use]
    pub fn depth(&self) -> u32 {
        self.depth
    }

    /// **`DP` minus the sum of `AD`: the reads no written allele explains.** Stutter at a
    /// tract, a candidate selection dropped, or noise at a SNP — the per-sample artifact
    /// signal a downstream filter can read without re-running anything.
    #[inline]
    #[must_use]
    pub fn unexplained_reads(&self) -> u32 {
        // PANIC-FREE: `Self::new` refused a depth below this sum, so the sum fits a `u32` and
        // the subtraction cannot go below zero; `saturating_sub` makes that structural rather
        // than a claim a reader has to check.
        let attributed: u32 = self.allele_reads.iter().sum();
        self.depth.saturating_sub(attributed)
    }
}

/// **The cohort-pooled mapping qualities of one allele's reads** — what `MQREF` and `MQALT`
/// are means of.
///
/// **Pooled over every sample of the cohort that had a read on this allele, called or not.**
/// That is production's rule — its writer sums each sample's per-allele observation count and
/// mapping-quality total with no test on whether the sample was called
/// (`src/vcf/record_encode.rs:391-401`) — and it is the same principle spec §7 applies in
/// writing a no-called sample's `DP` and `AD` beside its `./.`: the evidence exists even where
/// the call does not.
///
/// A real variant's reads map about as well as the reference's; reads carried in from a
/// paralogous copy elsewhere in the genome map worse, and `MQDIFF` — this allele's mean less
/// the reference's — is where that shows.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct MapqPool {
    /// How many reads across the cohort supported this allele.
    ///
    /// **The same reads `AD` counts**, which is why [`VcfRecord::new`] asserts this equals the
    /// samples' `AD` for the allele: in production the two are literally one field, and a pool
    /// summed over a different set would give an `MQDIFF` computed over reads the `AD` beside
    /// it does not describe.
    pub reads: u64,
    /// Their mapping qualities, summed.
    pub mapq_sum: u64,
}

impl MapqPool {
    /// The mean mapping quality, or `None` where no read in the cohort reached this allele.
    ///
    /// **`None` is not zero and the file spells them differently:** an allele nobody's reads
    /// reached has no mean, and writing `0` would claim every read mapped as badly as possible.
    #[inline]
    #[must_use]
    pub fn mean(self) -> Option<f64> {
        (self.reads > 0).then(|| self.mapq_sum as f64 / self.reads as f64)
    }
}

/// **The one `FILTER` value a record carries.**
///
/// Five values in one namespace: production's four — one from its SNP/indel writer, three from
/// its repeat-tract writer — plus `PASS`. There is no missing filter; a record always says
/// which of these it is.
///
/// **A loop that ran out of passes reaches the file as [`Self::EmDidNotConverge`] rather than
/// as a separate flag.** The plan's step A1 lists `converged` among the record's contents; it is
/// carried here instead, because spec §8 gives a record exactly one filter value and a boolean
/// beside it would be a second spelling of the same fact — the shape that lets two fields
/// disagree. Nothing is lost: the fact still reaches the output, which is what
/// `calling_em_loop.md` §6 requires of it.
///
/// **Which refusals are written and which are dropped is not this module's decision.** The
/// emission steps own it, and production answers it both ways today: its SNP/indel writer drops
/// every rejection silently, while its tract writer emits the locus on its filter with every
/// sample no-called. What the format fixes is the vocabulary and the shape of a written
/// refusal (spec §8).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FilterVerdict {
    /// Every filter passed.
    Pass,
    /// **The calling loop ran out of passes before its frequencies stopped moving.** The record
    /// is written, not dropped: one hard locus must not kill a cohort run, and a genotype from
    /// a loop that did not settle is a weaker claim than one from a loop that did — which
    /// nothing downstream could otherwise tell (`calling_em_loop.md` §6).
    EmDidNotConverge,
    /// The tract's allele-length distribution is inconsistent with its motif period.
    NotPeriodic,
    /// More candidate alleles segregate at the tract than the caller admits.
    TooManyAlleles,
    /// Too little cohort depth to call the tract.
    LowDepth,
}

impl FilterVerdict {
    /// The value as it is written in the `FILTER` column.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "PASS",
            Self::EmDidNotConverge => "EMNoConv",
            Self::NotPeriodic => "notPeriodic",
            Self::TooManyAlleles => "tooManyAlleles",
            Self::LowDepth => "lowDepth",
        }
    }
}

/// **What makes a record a repeat-tract record**: the tract's motif, from which the `STR` flag,
/// `RU` and `PERIOD` are all written.
///
/// **`PERIOD` is not stored, it is the motif's length.** The two cannot disagree, which is the
/// whole reason to hold one field rather than two — production writes the period and never the
/// motif at all, so a consumer wanting the repeat unit has to go back to the reference for it.
///
/// **Holds the crate's [`Motif`]** rather than bytes of its own: a motif is `1..=MAX_MOTIF_LEN`
/// bases by construction, never heap-allocated, and hands out its period as a checked type. So
/// there is no empty motif to refuse here and no `PERIOD=0` to write — the state is
/// unrepresentable rather than merely rejected.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct TractAnnotation {
    motif: Motif,
}

impl TractAnnotation {
    /// Build the annotation from the tract's motif, on the reference strand.
    #[inline]
    pub fn new(motif: Motif) -> Self {
        Self { motif }
    }

    /// The `RU` field: the motif's bases.
    #[inline]
    #[must_use]
    pub fn motif(&self) -> &[u8] {
        self.motif.as_bytes()
    }

    /// The `PERIOD` field: the motif's length in bases.
    #[inline]
    #[must_use]
    pub fn period(&self) -> usize {
        self.motif.period()
    }

    /// **How many whole repeat units an allele's sequence holds** — one entry of `REPCN`.
    ///
    /// Derived from the allele rather than stored, so it cannot drift from the sequence the
    /// same record writes. Taken over the tract sequence as the record holds it, which is
    /// before the encoder's padding base is applied: that base is presentation and would
    /// otherwise inflate a count by a fraction of a unit.
    ///
    /// Truncating, as production's is: a tract carrying a partial final unit reports the whole
    /// ones it has.
    #[inline]
    #[must_use]
    pub fn repeat_copies_of(&self, allele: &[u8]) -> usize {
        // PANIC-FREE: a `Motif` is at least one base by construction, so the period is never
        // zero and this division cannot trap.
        allele.len() / self.motif.ssr_period().get() as usize
    }
}

#[cfg(test)]
mod tests;
