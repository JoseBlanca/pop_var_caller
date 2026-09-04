//! **What a run has to say about itself when it finishes** — what it wrote, what it refused,
//! what it could not speak for, and which of its numbers rest on a measurement.
//!
//! # Why a run has to say anything at all
//!
//! **A VCF cannot distinguish ground the caller examined and found nothing at from ground it
//! never spoke for.** That is `cohort_merge.md` §3.3's whole argument for the failed-locus count
//! being *counted* and not merely dropped — "in the VCF alone, that absence is indistinguishable
//! from *analysed and found nothing*" — and the same holds for every other kind of nothing this
//! caller produces: a repeat tract it has not built a generator for, a locus where the allele cap
//! left nobody callable, a sample whose reads were all filtered away. Each looks, in the file,
//! exactly like a quiet stretch of genome.
//!
//! So the report is the run saying which. It is a small document and its job is arithmetic a
//! reader would otherwise have to do: **the analysed ground partitions into what was called and
//! the two kinds of what was not**, and every count here is a share of a stated whole rather
//! than a bare number.
//!
//! # Contigs are named, and that is a change
//!
//! `GenomeRegion`'s own `Display` writes `contig 0:15-15`, because a region is a position with no
//! reference beside it and cannot spend a contig table it does not have. A *run* has one, so
//! everything here that a person reads names its chromosome — `chr1:15-15` — and the index never
//! reaches the page. That was recorded as the place it stops being acceptable, and this is it.
//!
//! # What it is not
//!
//! **Not a log, and not a place to put spans.** Where the refused loci themselves surface — a
//! sidecar, a BED — is `cohort_merge.md` §14's open question 5; this states the counts and the
//! ground they cover, which is what §3.3 fixes as reaching the user.

use std::fmt::Write as _;

use crate::fasta::ContigList;
use crate::ng::calling::parameters_file::ParametersFile;
use crate::ng::psp::ParameterValue;
use crate::ng::read::filtering::ReadFilterCounts;
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::region_typing::GenomeRegions;
use crate::ng::types::{GenomeRegion, ReadGroupId};

use super::callers::{CohortCallingTallies, CohortWalkTallies, WrittenCohort};
use super::psp_caller::StoredCohortTallies;

/// **The two bounds a run's refusals were measured against**, so the report can name each beside
/// the count it explains — and name it as the flag a person types, which is the whole of what
/// makes the advice actionable.
#[derive(Clone, Copy, Debug)]
pub struct BoundsTheRunCalledUnder {
    /// `--max-cohort-locus-span`, in reference bases.
    pub max_cohort_locus_span: u32,
    /// `--max-candidate-alleles`, the reference counted among them.
    pub max_candidate_alleles: u16,
}

/// **A run's own account of itself**, gathered when it finishes and rendered for a person.
///
/// Built from what the run produced and the two things that give its numbers meaning: the
/// contig table, so ground can be named rather than numbered, and the parameters file, so the
/// numbers behind the calls can be said to rest on a measurement or not.
pub struct RunReport<'a> {
    /// What the calling itself produced — the same in both modes, because none of it depends
    /// on where the observations came from.
    calling: &'a CohortCallingTallies,
    /// **What the run knows about each sample, which is the one part that differs by mode.**
    samples: WhatEachSampleDid<'a>,
    /// The run's contigs, so a span can be named.
    contigs: &'a ContigList,
    /// The run's read groups, so a library can be named where its reads were dropped.
    read_groups: &'a ReadGroups,
    /// The numbers the run scored with, as the file beside its output states them.
    parameters: &'a ParametersFile,
    /// **Which ground the run undertook to speak for**, in a phrase — the whole argument for a
    /// report is that a VCF cannot say this, and a report that only counted the bases could not
    /// say it either. Measured in review: two runs printed "analysed ground: 3000 bases" and
    /// nothing on either page named the chromosome.
    where_the_ground_is: String,
    /// How much genome the run was asked to call over, in bases.
    analysed_bases: u64,
    /// The widest locus this run undertook to assemble, so its refusals can name the bound they
    /// were measured against.
    max_cohort_locus_span: u32,
    /// The most alleles this run called a locus over, for the same reason.
    max_candidate_alleles: u16,
}

/// **How this run came by its observations, and therefore what it can say about each sample.**
///
/// A walk counts what it saw — the ground it handled, the reads its filters kept and dropped —
/// and hands those counts to the report. **A stored file carries none of them**: they belong to
/// the read cursor of a walk that ran in another process, and what that walk counted is in the
/// file's provenance rather than in its records. So a run over stored files reports what it
/// *drew* instead, and this enum is where the report stops pretending the two are the same.
///
/// The alternative — one shape with the walk's fields left at zero in psp mode — was rejected
/// for the reason the report exists at all: a zero reads as *measured and none*, where the
/// truth is *not measured here*.
enum WhatEachSampleDid<'a> {
    /// Direct mode: every sample's reads were walked by this run.
    TheirReadsWereWalked(&'a CohortWalkTallies),
    /// psp mode: every sample's stored file was read by this run.
    TheirStoredFilesWereRead(&'a StoredCohortTallies),
}

/// How many refused spans a report prints before it stops and says how many are left.
///
/// **A run whose bound is badly set refuses thousands**, and `cohort_merge.md` §3.3 says what a
/// reader does with a non-zero count: inspect the spans, and if their lengths cluster just above
/// the bound, raise it and call again. A handful is enough to see a cluster; a thousand lines is
/// a log, which this is not.
const SPANS_A_REPORT_SHOWS: usize = 5;

impl<'a> RunReport<'a> {
    /// Gather what a finished run has to say.
    #[must_use]
    pub fn of(
        written: &'a WrittenCohort,
        contigs: &'a ContigList,
        read_groups: &'a ReadGroups,
        parameters: &'a ParametersFile,
        analysed: &GenomeRegions,
        bounds: BoundsTheRunCalledUnder,
    ) -> Self {
        Self {
            calling: &written.calling,
            samples: WhatEachSampleDid::TheirReadsWereWalked(&written.walk),
            contigs,
            read_groups,
            parameters,
            where_the_ground_is: describe(analysed, contigs),
            analysed_bases: analysed.iter().map(|region| region.len()).sum(),
            max_cohort_locus_span: bounds.max_cohort_locus_span,
            max_candidate_alleles: bounds.max_candidate_alleles,
        }
    }

    /// Gather what a finished run over **stored files** has to say (spec §5.3).
    ///
    /// **The calling half is the same value direct mode's report is built from**, because
    /// calling does not know where its observations came from. What differs is the second
    /// argument: `stored` says what this run drew out of each psp, where direct mode hands over
    /// what each sample's walk saw — see [`WhatEachSampleDid`].
    #[must_use]
    pub fn of_a_stored_cohort(
        calling: &'a CohortCallingTallies,
        stored: &'a StoredCohortTallies,
        contigs: &'a ContigList,
        read_groups: &'a ReadGroups,
        parameters: &'a ParametersFile,
        analysed: &GenomeRegions,
        bounds: BoundsTheRunCalledUnder,
    ) -> Self {
        Self {
            calling,
            samples: WhatEachSampleDid::TheirStoredFilesWereRead(stored),
            contigs,
            read_groups,
            parameters,
            where_the_ground_is: describe(analysed, contigs),
            analysed_bases: analysed.iter().map(|region| region.len()).sum(),
            max_cohort_locus_span: bounds.max_cohort_locus_span,
            max_candidate_alleles: bounds.max_candidate_alleles,
        }
    }

    /// **A read group as its file named it** — the sample and the `@RG ID`, not an index.
    ///
    /// The same argument the contigs are named under: a run holds the table, so a person reading
    /// its report should not have to count from zero to find out which library dropped their
    /// reads. In a cohort with several libraries a sample the number cannot be mapped back at
    /// all.
    fn library(&self, read_group: Option<ReadGroupId>) -> String {
        let Some(id) = read_group else {
            return "reads declaring no read group".to_owned();
        };
        // **`get` panics on an id outside the table**, which cannot happen: these ids come from
        // the cursors this same table opened. The scan is over a handful of entries.
        match self.read_groups.iter().find(|(known, _)| *known == id) {
            Some((_, group)) => format!("library {} of {}", group.id, group.sample),
            None => format!("read group {}", id.0),
        }
    }

    /// **The whole report, as the lines a person reads.**
    ///
    /// Lines rather than printed output, so that what the report says is something a test can
    /// hold — the summary's own text was the one part of this command a mutation could change
    /// with the suite still green.
    #[must_use]
    pub fn lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        self.what_was_written(&mut lines);
        self.what_was_not_called(&mut lines);
        self.what_the_reads_did(&mut lines);
        self.what_the_numbers_rest_on(&mut lines);
        lines
    }

    /// What reached the file, and what was called but established nothing.
    fn what_was_written(&self, lines: &mut Vec<String>) {
        lines.push(format!("records written: {}", self.calling.records_written));
        lines.push(format!(
            "loci called: {} — {} written, {} establishing no variant and so left out",
            self.calling.loci_called(),
            self.calling.records_written,
            self.calling.loci_called_but_not_written,
        ));
    }

    /// **The ground this run could not speak for, and why** — three different kinds of nothing,
    /// kept apart because a reader acts on each differently.
    ///
    /// # The shares are of what the walk covered, not of what the BED asked for
    ///
    /// **A repeat tract is typed and walked whole even where a BED cuts one** (spec §4.2's
    /// emission rule: findings whole, generic clipped), so the bases the walk was handed can
    /// exceed the bases asked for — measured in review, a BED of 120 bases inside two tracts
    /// charged 240 to *not built yet*, and dividing by the 120 printed **200.0%**. The parts
    /// are right; the denominator was wrong. Every share here is of the three parts' own sum,
    /// so they add to a hundred by construction, and the two totals are printed side by side
    /// whenever they differ so nobody has to work out why.
    fn what_was_not_called(&self, lines: &mut Vec<String>) {
        self.the_ground_the_run_spoke_for(lines);
        self.what_the_loci_came_to(lines);
    }

    /// **How the analysed ground was accounted for, in bases** — printed only by a run that
    /// walked the reads itself.
    ///
    /// **A run over stored files cannot say this and does not guess.** The three parts are a
    /// walk's own region tally, and no psp records one: the ground a sample was walked over is
    /// in its header, but how much of it a generator spoke for is not. So psp mode prints the
    /// ground it called over and stops there, rather than printing three numbers it would have
    /// to invent — which is the same rule the per-sample section follows below.
    fn the_ground_the_run_spoke_for(&self, lines: &mut Vec<String>) {
        let WhatEachSampleDid::TheirReadsWereWalked(walk) = self.samples else {
            // **The ground still gets named**, because it is the one thing a VCF cannot say and
            // psp mode knows it exactly — every file's header records it and the cohort was
            // refused unless they agreed (spec §6.2).
            lines.push(format!(
                "analysed ground: {} — {} bases, as every file's header records them",
                self.where_the_ground_is, self.analysed_bases,
            ));
            return;
        };
        // **Every sample walks the same ground, so one sample's region tally is the run's.**
        // The loci each sample's walk emitted differ and are not this section's business.
        let Some(ground) = walk.per_sample.first().map(|walk| &walk.regions) else {
            return;
        };
        let covered = ground.regions_handled_bp
            + ground.unhandled_not_implemented_bp
            + ground.unhandled_out_of_scope_bp;
        lines.push(format!(
            "analysed ground: {} — {} bases asked for, in {} typed region{}",
            self.where_the_ground_is,
            self.analysed_bases,
            ground.regions_in,
            plural(ground.regions_in),
        ));
        if covered != self.analysed_bases {
            lines.push(format!(
                "  the walk covers a typed region whole, so it spoke for {covered} bases; the \
                 shares below are of that",
            ));
        }
        lines.push(format!(
            "  called: {} bases ({})",
            ground.regions_handled_bp,
            share_of(ground.regions_handled_bp, covered),
        ));
        lines.push(format!(
            "  not called — clusters of repeats too close together to have clean flanks: {} \
             bases ({})",
            ground.unhandled_not_implemented_bp,
            share_of(ground.unhandled_not_implemented_bp, covered),
        ));
        lines.push(format!(
            "  not called — tandem arrays longer than this run types as callable: {} bases ({})",
            ground.unhandled_out_of_scope_bp,
            share_of(ground.unhandled_out_of_scope_bp, covered),
        ));
    }

    /// **What became of the loci the run built** — which is a different fact from how much
    /// ground it covered, and one both modes can state, because it is the calling's.
    ///
    /// **Loci, not bases, and a different fact from the base lines above.** Those say what
    /// ground no generator looked at; these say what was looked at, built and merged across
    /// the cohort — and then called, filtered, or not scored at all. A run that printed only
    /// the base lines would say a refused tract's ground was *called*, which it was not.
    fn what_the_loci_came_to(&self, lines: &mut Vec<String>) {
        // **The two refusals are here because the file cannot say them.** A tract refused as not
        // periodic is called over the reference tract alone, so every sample is homozygous
        // reference and no record is written; in the file it is indistinguishable from a tract
        // nobody varied at (`doc/devel/ng/spec/vcf_output.md` §9).
        let tracts = &self.calling.tracts;
        if tracts.built() > 0 {
            lines.push(format!(
                "repeat tracts: {} built, of which {} called",
                tracts.built(),
                tracts.called,
            ));
            if tracts.not_periodic > 0 {
                lines.push(format!(
                    "  not called — the reads do not vary in whole motif units (notPeriodic): {}",
                    tracts.not_periodic,
                ));
            }
            if tracts.too_many_alleles > 0 {
                lines.push(format!(
                    "  called over fewer sequences than segregate there (tooManyAlleles): {}",
                    tracts.too_many_alleles,
                ));
            }
            if tracts.without_whole_repeats > 0 {
                lines.push(format!(
                    "  not called — a candidate shorter than one copy of the motif: {}",
                    tracts.without_whole_repeats,
                ));
            }
            if tracts.bundles_set_aside > 0 {
                lines.push(format!(
                    "  not called — clusters of repeats with no clean flanks, which nothing \
                     builds a caller for yet: {}",
                    tracts.bundles_set_aside,
                ));
            }
        }

        self.refused_loci(
            lines,
            "loci the merge declined to assemble for being too wide",
            &self.calling.loci_too_wide_to_assemble,
            &format!(
                "the bound is --max-cohort-locus-span {}; raise it and call again if the \
                 lengths above cluster just past it",
                self.max_cohort_locus_span,
            ),
        );
        self.refused_loci(
            lines,
            "loci where the allele cap left no sample callable",
            &self.calling.loci_with_nobody_to_call,
            &format!(
                "the cap is --max-candidate-alleles {}; raise it to keep more of what they vary \
                 over",
                self.max_candidate_alleles,
            ),
        );
    }

    /// One refusal's count, a handful of its spans by name, and what a reader does about it.
    ///
    /// **The advice is printed only where the count is non-zero**, because a line telling
    /// somebody how to fix a thing that did not happen is a line they have to read and discard.
    fn refused_loci(
        &self,
        lines: &mut Vec<String>,
        what: &str,
        spans: &[GenomeRegion],
        what_to_do: &str,
    ) {
        lines.push(format!("{what}: {}", spans.len()));
        if spans.is_empty() {
            return;
        }
        for span in spans.iter().take(SPANS_A_REPORT_SHOWS) {
            // **The length beside the span**, because the advice below asks a reader to compare
            // it against a bound and subtracting two coordinates by eye is what they would
            // otherwise do, once a line, on every line.
            lines.push(format!("  {} ({} bases)", self.named(*span), span.len()));
        }
        if spans.len() > SPANS_A_REPORT_SHOWS {
            lines.push(format!(
                "  … and {} more",
                spans.len() - SPANS_A_REPORT_SHOWS
            ));
        }
        lines.push(format!("  {what_to_do}"));
    }

    /// **What each sample's reads did**, and the two ways a sample can contribute nothing.
    ///
    /// **Three situations were written with one sentence and two of them were wrong** (review,
    /// 2026-09-01). A sample the filters emptied — every read a duplicate, every read below the
    /// mapping-quality floor — was printed as having *no reads* four lines above the line saying
    /// it had 720 of them; and so was a sample whose reads exist over the analysed ground but
    /// whose ground no generator walked. A geneticist does something different in each case:
    /// check the sample sheet, check the duplicate marking, check nothing at all.
    ///
    /// **And what the file says about such a sample is `./.`, not a genotype.** This line used
    /// to claim the sample "still carries a genotype, from the prior alone", which is what the
    /// *loop* does and not what the *file* writes: `vcf_output.md` §7.1 no-calls a sample whose
    /// likelihoods are flat, and F1 implemented it. The claim was reassurance that stopped a
    /// reader opening the file, and the file said the opposite.
    fn what_the_reads_did(&self, lines: &mut Vec<String>) {
        match self.samples {
            WhatEachSampleDid::TheirReadsWereWalked(walk) => self.what_each_walk_saw(lines, walk),
            WhatEachSampleDid::TheirStoredFilesWereRead(stored) => {
                what_each_stored_file_gave(lines, stored);
            }
        }
    }

    /// **What each sample's reads did**, for a run that walked them — the section above's body.
    fn what_each_walk_saw(&self, lines: &mut Vec<String>, walked: &CohortWalkTallies) {
        let mut spoke = 0_usize;
        let mut emptied_by_the_filters = Vec::new();
        let mut nothing_reached_the_caller = Vec::new();
        for walk in &walked.per_sample {
            let did: Vec<WhatTheFiltersDid> = walk
                .read_filters
                .iter()
                .map(|(_, counts)| what_the_filters_did(counts))
                .collect();
            let kept: u64 = did.iter().map(|did| did.kept).sum();
            let dropped: u64 = did.iter().map(|did| did.dropped).sum();
            match (kept, dropped) {
                (0, 0) => nothing_reached_the_caller.push(walk.sample_name.as_str()),
                (0, _) => emptied_by_the_filters.push(walk.sample_name.as_str()),
                _ => spoke += 1,
            }
        }
        lines.push(format!(
            "samples: {} — {spoke} whose reads the caller used, {} whose reads the filters took \
             all of, {} that contributed none",
            walked.per_sample.len(),
            emptied_by_the_filters.len(),
            nothing_reached_the_caller.len(),
        ));
        if !emptied_by_the_filters.is_empty() {
            lines.push(format!(
                "  every read filtered out: {} — written ./. wherever that left nothing to score",
                emptied_by_the_filters.join(", "),
            ));
        }
        if !nothing_reached_the_caller.is_empty() {
            lines.push(format!(
                "  no read reached the caller: {} — either the sample has none over this ground, \
                 or none of its ground was walked; written ./.",
                nothing_reached_the_caller.join(", "),
            ));
        }

        for walk in &walked.per_sample {
            let did: Vec<(Option<ReadGroupId>, WhatTheFiltersDid)> = walk
                .read_filters
                .iter()
                .map(|(read_group, counts)| (*read_group, what_the_filters_did(counts)))
                .collect();
            let kept: u64 = did.iter().map(|(_, did)| did.kept).sum();
            let dropped: u64 = did.iter().map(|(_, did)| did.dropped).sum();
            // **A sample that read nothing has nothing to say about its filters**, and the
            // line above already named it. Two zeroes here would be a line to discard.
            if kept == 0 && dropped == 0 {
                continue;
            }
            lines.push(format!(
                "  {}: {kept} reads kept, {dropped} dropped by the read filters",
                walk.sample_name,
            ));
            for (read_group, did) in &did {
                if did.why.is_empty() {
                    continue;
                }
                lines.push(format!("    {}: {}", self.library(*read_group), did.why));
            }
            // **Not a drop, and on its own line** — a file's read groups need not all belong to
            // the sample being opened, and a record of somebody else's says nothing about how
            // this one behaved.
            let elsewhere: u64 = did
                .iter()
                .map(|(_, did)| did.records_of_another_sample)
                .sum();
            if elsewhere > 0 {
                lines.push(format!(
                    "    {elsewhere} records in these files belong to another sample and were \
                     skipped, which is not a drop",
                ));
            }
        }
    }

    /// **Which of the run's numbers rest on a measurement of reads, and which are constants** —
    /// `parameters_file.md` §8's question, answered from the file the run wrote beside its VCF.
    ///
    /// **The groups that were *not* fitted are what the reader needs**, and they are named
    /// rather than counted: a run whose contamination and slippage are compiled-in constants is
    /// a different claim from one whose calibration is, and a count of five says neither.
    fn what_the_numbers_rest_on(&self, lines: &mut Vec<String>) {
        let fitted = self.parameters.what_the_run_fitted();
        // **Not "parameters", which is the line naming the file.** Two lines opening on one word
        // and saying different things is two lines a reader has to disambiguate.
        lines.push(format!(
            "numbers behind the calls: {} of {} groups the file says were fitted",
            fitted.fitted().len(),
            fitted.groups(),
        ));
        if !fitted.not_fitted().is_empty() {
            lines.push(format!(
                "  taken from constants or supplied, not measured here: {}",
                fitted
                    .not_fitted()
                    .iter()
                    .map(|group| group.in_the_readers_words())
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }
    }

    /// A span with its chromosome named — [`name_of`], with the run's own table.
    fn named(&self, region: GenomeRegion) -> String {
        name_of(region, self.contigs)
    }
}

/// **What each sample's stored file gave this run** — psp mode's counterpart of the section
/// above, and deliberately a different set of facts.
///
/// **Both numbers are this run's own measurements, not the file's claims.** A psp records no
/// count of reads kept or dropped, so nothing here restates a walk; what it states is how many
/// stored loci this run drew out of each file, and how many reads went into the comparison at
/// one of them on average — summed from the record head's `reads-compared-with-reference` as
/// the loci went past.
///
/// **That average is not the sample's depth, and the line does not call it depth.** The head's
/// count is the keep rule's own denominator, and
/// [its field](crate::ng::psp::RecordHead::reads_compared_with_reference) lists what it leaves
/// out: reads a filter turned away, reads the per-position cap discarded, reads that covered
/// the locus and produced no observation, and reads whose witness stopped inside it. At a
/// repeat tract that 40 reads cover but only 22 anchor both borders of, this says 22. It is the
/// number the admission rule was applied with, which is why it is worth printing; a reader who
/// wants coverage wants a different number and no psp carries one.
///
/// **A file that gave no loci is named, and that is the third of direct mode's three cases.**
/// A sample whose reads the filters emptied and one that had no reads at all are the same
/// thing here — an empty file — because the walk that could tell them apart is not this run.
///
/// **The read filters are compared, and this is the only place they are.** Every psp records
/// the filters its walk applied and spec §6.1 compares them against nothing, so without this a
/// cohort assembled from files walked under different filters would call without a word about
/// it. The line is printed only where the files disagree, which is never on a cohort walked by
/// one `generate-psps` invocation.
///
/// **What can differ today is the build, not a flag.** `generate-psps` exposes no read-filter
/// option and walks with the compiled-in policy, so two psps disagree here only when they were
/// written by ng builds whose defaults differ — a cohort walked over months, which is exactly
/// the cohort psp mode exists to make possible (spec §2). The moment those settings become
/// flags, which Milestone C left as the owner's call, this becomes the check on what a person
/// typed.
fn what_each_stored_file_gave(lines: &mut Vec<String>, stored: &StoredCohortTallies) {
    let empty: Vec<&str> = stored
        .per_sample
        .iter()
        .filter(|sample| sample.read.loci_read == 0)
        .map(|sample| sample.sample_name.as_str())
        .collect();
    lines.push(format!(
        "samples: {} — {} whose stored file gave this run loci, {} whose file held none over \
         this ground",
        stored.per_sample.len(),
        stored.per_sample.len() - empty.len(),
        empty.len(),
    ));
    if !empty.is_empty() {
        lines.push(format!(
            "  no locus over this ground: {} — written ./. wherever that left nothing to score",
            empty.join(", "),
        ));
    }
    for sample in &stored.per_sample {
        // **A sample whose file held nothing has nothing to say about its depth**, and the line
        // above already named it. `0 loci read, — reads a locus` would be a line to discard.
        let Some(depth) = sample.read.mean_reads_a_locus() else {
            continue;
        };
        lines.push(format!(
            "  {}: {} loci read, {depth:.1} reads a locus compared with the reference",
            sample.sample_name, sample.read.loci_read,
        ));
    }
    where_the_walks_disagreed(lines, stored);
}

/// **Where the cohort's files were not walked alike**, one line a setting that differs.
///
/// Each line names the setting as its key is spelled in the psp header, then every value the
/// cohort holds with the samples that hold it — so a reader sees at once which files are the
/// odd ones and by how much. Nothing is printed where every file agrees.
fn where_the_walks_disagreed(lines: &mut Vec<String>, stored: &StoredCohortTallies) {
    // **Keys from every file, not from the first**, so a setting one walk recorded and another
    // did not is a disagreement rather than an absence nobody looks at.
    let mut settings: Vec<&str> = stored
        .per_sample
        .iter()
        .flat_map(|sample| {
            sample
                .read_filters_the_walk_applied
                .keys()
                .map(String::as_str)
        })
        .collect();
    settings.sort_unstable();
    settings.dedup();

    let mut differing: Vec<String> = Vec::new();
    for setting in settings {
        // **Values in the order the samples were given**, and each carrying its samples, so the
        // line reads as *these files say this, those say that*.
        let mut values: Vec<(String, Vec<&str>)> = Vec::new();
        for sample in &stored.per_sample {
            let value = match sample.read_filters_the_walk_applied.get(setting) {
                Some(value) => as_a_person_reads_it(value),
                // **Absent is a value here.** A file whose walk did not record this setting
                // differs from one that did, and saying so is the point of the check.
                None => "not recorded".to_owned(),
            };
            match values.iter_mut().find(|(seen, _)| *seen == value) {
                Some((_, samples)) => samples.push(sample.sample_name.as_str()),
                None => values.push((value, vec![sample.sample_name.as_str()])),
            }
        }
        if values.len() < 2 {
            continue;
        }
        differing.push(format!(
            "    {setting}: {}",
            values
                .iter()
                .map(|(value, samples)| format!("{value} for {}", samples.join(", ")))
                .collect::<Vec<_>>()
                .join("; "),
        ));
    }
    if differing.is_empty() {
        return;
    }
    lines.push(
        "  not every file was walked under the same read filters, and nothing refuses that — \
         these samples were judged on different reads:"
            .to_owned(),
    );
    lines.extend(differing);
}

/// **A psp header's recorded setting, as the header itself spells it.**
///
/// A quoted string would put quotation marks round `off`, which is how a switched-off filter
/// is spelled in the header — so the value goes out bare, exactly as somebody reading the
/// file's own TOML would see it.
fn as_a_person_reads_it(value: &ParameterValue) -> String {
    match value {
        ParameterValue::Integer(number) => number.to_string(),
        ParameterValue::Float(number) => number.to_string(),
        ParameterValue::Boolean(flag) => flag.to_string(),
        ParameterValue::String(text) => text.clone(),
    }
}

/// **Which ground a run undertook to speak for**, short enough for one line.
///
/// One interval is named whole; several are the first and the last with the count between them,
/// which is what a person needs to tell one BED from another at a glance. A run over every base
/// of every contig says so rather than listing them.
pub(crate) fn describe(analysed: &GenomeRegions, contigs: &ContigList) -> String {
    let mut spans = analysed.iter();
    let Some(first) = spans.next() else {
        return "nothing".to_owned();
    };
    let rest = spans.count();
    let named = |region: GenomeRegion| name_of(region, contigs);
    if rest == 0 {
        return named(first);
    }
    let last = analysed.iter().last().unwrap_or(first);
    format!("{} intervals, {} … {}", rest + 1, named(first), named(last))
}

/// A span with its chromosome named, 1-based and inclusive — `chr1:15-15`.
///
/// **Falls back to the index only where the table has no such contig**, which is a table and a
/// span built over different references and cannot happen in a run; saying `contig 7` there is
/// more use than a panic in a report.
fn name_of(region: GenomeRegion, contigs: &ContigList) -> String {
    let mut out = String::new();
    match contigs.entries.get(region.contig.get() as usize) {
        Some(entry) => out.push_str(&entry.name),
        None => {
            let _ = write!(out, "contig {}", region.contig.get());
        }
    }
    let _ = write!(out, ":{}-{}", region.start.get(), region.end.get());
    out
}

/// The `s` that makes a count read, or nothing where the count is one.
pub(crate) fn plural(n: u64) -> &'static str {
    if n == 1 { "" } else { "s" }
}

/// `n` as a share of `whole`, or `—` where the whole is zero.
///
/// **A percentage and not a natural frequency, and the reason is the range.** These shares run
/// from a few bases in a hundred million to nearly all of the ground, and *6.0%* reads where
/// *6 bases in every 100* does not once the numerator is 40 million.
pub(crate) fn share_of(n: u64, whole: u64) -> String {
    if whole == 0 {
        return "—".to_owned();
    }
    format!("{:.1}%", 100.0 * n as f64 / whole as f64)
}

/// **How many of a read group's reads the filters turned down, and why** — both answers from
/// one exhaustive destructuring, so a counter added to [`ReadFilterCounts`] later must be
/// routed here or this stops compiling.
///
/// **That guarantee is the reason for the shape**, and it is `ReadFilterCounts::add`'s own: a
/// tally that silently stopped counting one reason "would under-report drops without failing
/// anything". Naming the fields one by one gives it up, which is what the first version of this
/// did.
///
/// **`other_sample` is not a drop and is not counted as one.** Its own field says so — *"the
/// rest of this struct answers 'how did this read group behave?', and a read belonging to
/// someone else says nothing about that: counting it as a drop would make a shared file look
/// like a low-quality one"*. On a cohort sharing one multi-sample BAM it would dominate every
/// sample's drop count. It is reported on its own terms instead.
fn what_the_filters_did(counts: &ReadFilterCounts) -> WhatTheFiltersDid {
    let ReadFilterCounts {
        kept,
        duplicate,
        low_mapq,
        supplementary,
        secondary,
        unmapped,
        qc_fail,
        too_short,
        high_mismatch_fraction,
        bad_cigar,
        other_sample,
    } = counts;

    // One entry a `DropReason`, and there are nine of them.
    let by_reason = [
        ("duplicate", *duplicate),
        ("mapping quality too low", *low_mapq),
        ("supplementary", *supplementary),
        ("secondary", *secondary),
        ("unmapped", *unmapped),
        ("failed the sequencer's own QC", *qc_fail),
        ("too short", *too_short),
        ("too many mismatches", *high_mismatch_fraction),
        ("unusable CIGAR", *bad_cigar),
    ];
    WhatTheFiltersDid {
        kept: *kept,
        dropped: by_reason.iter().map(|(_, n)| n).sum(),
        // **Only the reasons that fired.** Nine reasons and a run trips two or three; the other
        // six as zeros are six numbers a reader scans past to find the one that matters.
        why: by_reason
            .iter()
            .filter(|(_, n)| *n > 0)
            .map(|(why, n)| format!("{n} {why}"))
            .collect::<Vec<_>>()
            .join(", "),
        records_of_another_sample: *other_sample,
    }
}

/// What one read group's filters did, as the report states it.
struct WhatTheFiltersDid {
    /// Reads that came through.
    kept: u64,
    /// Reads a filter turned down, over the nine reasons there are.
    dropped: u64,
    /// Those reasons that fired, in the reader's words — empty where none did.
    why: String,
    /// **Records this open skipped as another sample's**, which is not a drop: a file's read
    /// groups need not all belong to the sample being opened.
    records_of_another_sample: u64,
}

#[cfg(test)]
mod tests;
