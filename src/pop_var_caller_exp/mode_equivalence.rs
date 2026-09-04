//! **The oracle that justifies psp mode**: one cohort, one set of parameters, called both ways,
//! compared byte for byte (spec [`run_streaming.md`](../../../doc/devel/ng/spec/run_streaming.md)
//! §12.3).
//!
//! # Why this is a module of its own
//!
//! It belongs to neither command. `call-from-alignments` and `generate-psps` + `call-from-psps`
//! are two routes to one answer, and the claim under test is about the *pair* — spec §1.1's goal
//! 1, that the calling function cannot tell where its observations came from, and the
//! sufficiency test for the format: anything a psp fails to carry surfaces here, where a
//! write-read round trip would pass.
//!
//! # What is compared: the whole file
//!
//! Every byte — the `##contig` lines with their checksums, the `##INFO` and `##FORMAT`
//! declarations, the sample columns in their order, and every genotype, likelihood and quality
//! of every record. **Including `##commandline`**, which records the command a person typed and
//! which is here the *test binary's* command line, the same on both routes because both routes
//! run inside one process.
//!
//! **At a shell it is not the same, and that is the one exemption the real-data oracle takes.**
//! `call-from-alignments` and `call-from-psps` are typed differently, so the line differs by
//! construction rather than by anything the calls did, and
//! `scripts/ng_mode_equivalence_oracle.sh` drops it. That is the same kind of exemption
//! spec §12.1 already grants the psp header's timestamp: a field a comparison deliberately
//! skips, named, rather than a comparison quietly weakened. Everything the calls themselves
//! decide is compared here without any exemption at all.
//!
//! **Both routes write the same `--output` path, and the file is removed between them.** The
//! header names the parameters file beside the VCF, so two different output names would make two
//! different headers for identical calls and the comparison would fail on a filename. Reusing
//! one path costs something else, and it was measured: with the second route writing its VCF
//! anywhere but that path, the comparison re-read the first route's file and compared it with
//! itself — a green test making no claim at all. So the path is deleted before the second route
//! runs, and a run that did not recreate it fails on the read.
//!
//! # The run-level invariances live here too
//!
//! Spec §12's other properties of a psp-mode run — that the order the files are named in does
//! not change the calls (§12.6), that a cohort each of whose samples was walked in its own
//! invocation calls (§12.7), and that ground a sample analysed and found nothing in reads back
//! as analysed and empty rather than as never looked at (§12.9) — are all comparisons between
//! whole runs, so they are here rather than in either command's own tests. Each is a run
//! against a run, which is what this module is for.
//!
//! # The real-data half is a script, and it is where the claim has weight
//!
//! What runs here is a fixture: one contig, two samples, three variants, twenty reads a
//! position. The claim that psp mode calls what direct mode calls is made on real reads by
//! `scripts/ng_mode_equivalence_oracle.sh` — measured 2026-09-04 on six tomato accessions over
//! the first two 100 kb intervals of `benchmarks/tomato1/regions.bed`: **599 records,
//! byte-identical apart from that one line**, with the parameters file beside each identical
//! too. The fixture is what keeps it from silently stopping being true between such runs.
//!
//! # Where this oracle stops, and it is worth knowing before trusting it
//!
//! **It compares VCFs, so a stored locus that produces no record is not compared.** 578 of the
//! fixture's 581 stored loci per sample carry no variant, and a psp route that dropped one of
//! them would pass — the file has no record there either way. What holds the psp equal to the
//! walk *field for field* is Milestone B's own oracle
//! (`examples/ng_psp_gather_oracle.rs`), not this one.
//!
//! **A defaults run cannot see which read group a read belongs to**, because every group is
//! then scored with the same numbers — so under `--defaults` alone a psp route that dropped the
//! walk-local-to-run-wide renumbering entirely passes, which was measured. The second test here
//! calls both routes under a parameters file whose three read groups carry *different*
//! base-quality multipliers, which makes the identity observable in every genotype quality; it
//! is also the only place psp mode's supplied-parameters path meets direct mode's.
//!
//! **What no run here exercises** is what only a parameters *fit* reads: the stored sum of
//! squared mapping qualities and the count of reads that covered a locus without producing an
//! observation. Destroying either on write leaves both comparisons green. That is a property of
//! a run that does not fit rather than of the fixture, and the real-data script shares it.

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use crate::ng::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES;
    use crate::ng::calling::parameters_file::beside_the_vcf;
    use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
    use crate::ng::region_typing::segment_criteria::{
        DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
    };
    use crate::ng::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN;
    use crate::pop_var_caller_exp::call_from_alignments::{
        CallFromAlignmentsArgs, run_call_from_alignments,
    };
    use crate::pop_var_caller_exp::call_from_psps::{CallFromPspsArgs, run_call_from_psps};
    use crate::pop_var_caller_exp::generate_psps::{
        GeneratePspsArgs, psp_path_for, run_generate_psps,
    };
    use crate::pop_var_caller_exp::test_fixtures::{
        AVaryingCohort, FIRST_SAMPLES_SUBSTITUTION, SECOND_SAMPLES_SUBSTITUTION, TRACT,
        a_varying_cohort_on_disk,
    };

    /// A VCF's whole contents, as the lines a comparison holds and a failure prints.
    ///
    /// **Lines rather than bytes**, only so that a failure names which line differs instead of
    /// saying two files are unequal. Nothing is filtered out.
    fn comparable(vcf: &Path) -> Vec<String> {
        std::fs::read_to_string(vcf)
            .expect("the VCF reads")
            .lines()
            .map(str::to_owned)
            .collect()
    }

    /// How many records a VCF holds — the lines that are not header.
    fn records(lines: &[String]) -> usize {
        lines.iter().filter(|line| !line.starts_with('#')).count()
    }

    /// Call the cohort by walking its alignment files, and hand back the VCF's comparable lines.
    ///
    /// **`output` is the same path both routes write**, and each route's caller renames it
    /// afterwards — see the module note.
    fn called_from_alignments(cohort: &AVaryingCohort, output: &Path) -> Vec<String> {
        run_call_from_alignments(&alignments_args(cohort, output))
            .expect("the cohort calls from its alignment files");
        comparable(output)
    }

    /// What `call-from-alignments` is handed for this cohort — the defaults run's arguments,
    /// which the supplied-parameters run below then edits.
    fn alignments_args(cohort: &AVaryingCohort, output: &Path) -> CallFromAlignmentsArgs {
        CallFromAlignmentsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output: output.to_path_buf(),
            regions: None,
            parameters: None,
            defaults: true,
            ploidy: None,
            build_index_if_missing: false,
            max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
            max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
            cohort_locus_builder_regions_len: None,
            threads: 0,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        }
    }

    /// Walk the cohort into psps and call those, and hand back the VCF's comparable lines.
    ///
    /// **The psps are named one by one in the order the alignment files were**, so both routes
    /// put the sample columns in one order. Whether *any* order gives the same calls is spec
    /// §12.6's question and step F3's, not this one's.
    fn called_from_psps(cohort: &AVaryingCohort, output: &Path, samples: &[&str]) -> Vec<String> {
        // **The first route's file goes before the second route runs.** Without this the
        // comparison below can read one file twice and pass with the psp route writing nothing
        // — measured: pointing `call-from-psps` at another path left the test green.
        std::fs::remove_file(output).expect("the first route wrote the file this one replaces");
        let psps = walked_into_psps(cohort);
        run_call_from_psps(&psps_args(cohort, output, &psps, samples))
            .expect("the stored cohort calls");
        assert!(
            output.is_file(),
            "the stored cohort's run must have written {}, and the comparison would otherwise \
             read the other route's file twice",
            output.display(),
        );
        comparable(output)
    }

    /// Walk the cohort into psps and hand back the directory they are in.
    fn walked_into_psps(cohort: &AVaryingCohort) -> PathBuf {
        let psps = cohort.directory.path().join("psps");
        run_generate_psps(&GeneratePspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output_dir: psps.clone(),
            regions: None,
            force: true,
            build_index_if_missing: false,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        })
        .expect("the cohort walks into psps");
        psps
    }

    /// What `call-from-psps` is handed for this cohort.
    ///
    /// **The psps are named one by one in the order the alignment files were**, so both routes
    /// put the sample columns in one order. Whether *any* order gives the same calls is spec
    /// §12.6's question and step F3's, not this one's.
    fn psps_args(
        cohort: &AVaryingCohort,
        output: &Path,
        psps: &Path,
        samples: &[&str],
    ) -> CallFromPspsArgs {
        CallFromPspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            psps: samples
                .iter()
                .map(|sample| psp_path_for(psps, sample))
                .collect::<Vec<PathBuf>>(),
            output: output.to_path_buf(),
            parameters: None,
            defaults: true,
            ploidy: None,
            max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
            max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
            cohort_locus_builder_regions_len: None,
            threads: 0,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        }
    }

    /// **The oracle: one cohort called both ways gives one VCF** (spec §12.3).
    ///
    /// **The record count is asserted before the comparison**, and that assertion is what stops
    /// this test passing for the wrong reason. Two empty files are equal; a fixture that stopped
    /// producing calls — a reference the catalog began routing to the repeat path, a read the
    /// filters began dropping — would leave a green test making no claim at all. The fixture is
    /// built to vary at three positions and the cohort is diploid, so the run writes at least
    /// those three records.
    #[test]
    fn one_cohort_called_both_ways_gives_one_vcf() {
        let cohort = a_varying_cohort_on_disk();
        // **One path, written twice and read between**, because the header names the parameters
        // file beside the VCF: two output names would give two different headers for identical
        // calls and the comparison would fail on a filename.
        let output = cohort.directory.path().join("calls.vcf");

        let from_alignments = called_from_alignments(&cohort, &output);
        let from_psps = called_from_psps(&cohort, &output, &["one", "two"]);

        assert!(
            records(&from_alignments) >= WHAT_THE_FIXTURE_VARIES_AT,
            "the fixture must produce calls or this comparison claims nothing; got {} records",
            records(&from_alignments),
        );
        // **What the records are, not merely how many.** A run writing three records somewhere
        // else entirely would satisfy the count and would be testing something other than what
        // the fixture says it varies at.
        the_records_are_where_the_fixture_varies(&from_alignments);
        // **The two samples' columns must differ somewhere**, or a defect that gave one
        // sample's observations to the other could not be seen: every record would carry the
        // same string twice whichever sample it came from.
        assert!(
            from_alignments
                .iter()
                .filter(|line| !line.starts_with('#'))
                .any(|line| {
                    let columns: Vec<&str> = line.split('\t').collect();
                    columns.len() > 10 && columns[9] != columns[10]
                }),
            "the two samples must be called differently somewhere:\n{}",
            from_alignments.join("\n"),
        );
        assert_eq!(
            from_alignments, from_psps,
            "the same cohort and the same parameters must give the same VCF, whole",
        );
    }

    /// **How many records the fixture must at least produce**: one at each sample's own
    /// substitution, and one at the repeat tract the first sample shortens.
    const WHAT_THE_FIXTURE_VARIES_AT: usize = 3;

    /// **Every record sits where the fixture varies, and every place it varies has a record.**
    ///
    /// The positions are 1-based in a VCF and 0-based in the fixture, and the tract's record
    /// starts one base before the tract so that a deletion has an anchoring base — so this
    /// compares against the three the fixture declares, allowing that one base of slack.
    fn the_records_are_where_the_fixture_varies(lines: &[String]) {
        let mut at: Vec<u64> = lines
            .iter()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                line.split('\t')
                    .nth(1)
                    .expect("a record names its position")
                    .parse()
                    .expect("a VCF position is a number")
            })
            .collect();
        at.sort_unstable();
        let varies_at = [
            FIRST_SAMPLES_SUBSTITUTION as u64 + 1,
            TRACT.0 as u64,
            SECOND_SAMPLES_SUBSTITUTION as u64 + 1,
        ];
        assert_eq!(
            at.len(),
            varies_at.len(),
            "one record a place the fixture varies, and got {at:?} against {varies_at:?}",
        );
        for (found, expected) in at.iter().zip(varies_at) {
            assert!(
                found.abs_diff(expected) <= 1,
                "a record at {found} where the fixture varies at {expected}; all of them: \
                 {at:?} against {varies_at:?}",
            );
        }
    }

    /// **The fixture's ground is typed the way its doc says it is**, and this assertion is what
    /// keeps that true rather than a reader's word.
    ///
    /// Three regions: ordinary sequence, the deliberate repeat tract, ordinary sequence. Two
    /// things could break it silently, and both would leave the oracle green while it stopped
    /// testing what it claims — the tract falling below its floor and the whole contig becoming
    /// generic again, or the generator's own sequence drifting a base and some other stretch
    /// crossing a floor. Measured on the shipped sequence: the longest accidental tandem stretch
    /// is 5.5 copies of a period-2 motif, one copy short of its six-copy floor.
    #[test]
    fn the_fixtures_ground_is_typed_as_its_doc_says() {
        use crate::ng::region_typing::RegionKind;
        use crate::pop_var_caller_exp::run_ground::{self, GroundRequest, RepeatRouting};

        let cohort = a_varying_cohort_on_disk();
        let request = GroundRequest {
            reference: &cohort.reference,
            catalog: Some(&cohort.catalog),
            regions: None,
            routing: RepeatRouting {
                min_copies: MinCopies::default(),
                min_period: DEFAULT_MIN_PERIOD,
                max_period: DEFAULT_MAX_PERIOD,
                max_str_len: DEFAULT_MAX_STR_LEN,
                min_purity: DEFAULT_MIN_PURITY,
            },
        };
        let (_info, with_checksums) = the_references_two_views(&cohort.reference);
        let analysed = run_ground::analysed_regions(&request, &with_checksums.contig_list())
            .expect("the whole reference");
        let segmentation =
            run_ground::segments_over(&request, &analysed, &with_checksums).expect("it types");

        let kinds: Vec<(u64, u64, bool)> = segmentation
            .segments()
            .iter()
            .map(|region| {
                (
                    region.region.start.get(),
                    region.region.end.get(),
                    matches!(region.kind, RegionKind::SsrSegment(_)),
                )
            })
            .collect();
        assert_eq!(
            kinds.len(),
            3,
            "ordinary sequence, the tract, ordinary sequence — and got {kinds:?}",
        );
        assert!(
            !kinds[0].2 && kinds[1].2 && !kinds[2].2,
            "only the middle region is a repeat tract, and got {kinds:?}",
        );
        assert_eq!(
            (kinds[1].0, kinds[1].1),
            (TRACT.0 as u64 + 1, TRACT.1 as u64),
            "the tract is where the fixture put it, 1-based and inclusive: {kinds:?}",
        );
    }

    /// The reference as a run holds it: opened for bases, and read to the end for its per-contig
    /// checksums. The catalog check needs the second.
    fn the_references_two_views(
        fasta: &Path,
    ) -> (
        std::sync::Arc<crate::ng::reference_info::ReferenceInfo>,
        std::sync::Arc<crate::ng::reference_info::ReferenceInfo>,
    ) {
        use crate::ng::reference_info::{
            ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
        };
        let cache = std::sync::Arc::new(ReferenceInfoCache::new());
        let (info, verify) = read_reference_verifying_or_creating_fai(
            &cache,
            fasta.to_path_buf(),
            ReferenceCheck::VerifyAgainstIndex,
        )
        .expect("the fixture's reference reads");
        let with_checksums = match verify {
            Some(handle) => handle.join().expect("it verifies"),
            None => std::sync::Arc::clone(&info),
        };
        (info, with_checksums)
    }

    /// **The same cohort under a parameters file whose read groups differ**, called both ways,
    /// gives one VCF — and unlike the defaults run above, this one can see which read group a
    /// read belongs to.
    ///
    /// **Why it exists**: each sample's psp numbers its own read groups from zero, and the
    /// calling run merges those tables into one run-wide numbering (spec §6.2). Under
    /// `--defaults` every group is scored with the same numbers, so dropping that renumbering
    /// changes no genotype and the comparison above stays green — measured. Give the three
    /// groups different base-quality multipliers and the identity reaches every likelihood, so a
    /// read scored against another group's calibration moves the record it is in.
    ///
    /// It is also the only place psp mode's supplied-parameters path is compared against direct
    /// mode's: the by-name match of samples and read groups (spec §6.2, §12.5), the reference
    /// binding, and the numbers reaching the calling loop.
    #[test]
    fn one_cohort_called_both_ways_under_read_groups_that_score_differently_gives_one_vcf() {
        let cohort = a_varying_cohort_on_disk();
        let output = cohort.directory.path().join("calls.vcf");

        // A defaults run writes the file this test then edits — the shape spec §7 invites, and
        // the only way to get a parameters file naming exactly this cohort's read groups.
        called_from_alignments(&cohort, &output);
        let parameters = cohort.directory.path().join("differing.parameters.toml");
        std::fs::rename(beside_the_vcf(&output), &parameters).expect("the file moves");
        make_the_read_groups_score_differently(&parameters);

        let from_alignments = called_from_alignments_under(&cohort, &output, &parameters);
        let from_psps = called_from_psps_under(&cohort, &output, &["one", "two"], &parameters);

        the_records_are_where_the_fixture_varies(&from_alignments);
        assert_eq!(
            from_alignments, from_psps,
            "the same cohort under the same parameters must give the same VCF, whole",
        );
    }

    /// Give each read group its own base-quality multiplier, so that which group a read belongs
    /// to reaches the genotype it is scored into.
    ///
    /// **Three different numbers and none of them one**, because one is the identity: a file
    /// that left any group at 1.0 would score that group's reads exactly as the defaults run
    /// did, and a mis-attribution into it would be invisible again.
    fn make_the_read_groups_score_differently(parameters: &Path) {
        let text = std::fs::read_to_string(parameters).expect("the parameters file reads");
        let mut rewritten = String::with_capacity(text.len());
        for line in text.lines() {
            let multiplier =
                ["0.25", "2.5", "4.0"]
                    .iter()
                    .enumerate()
                    .find_map(|(group, value)| {
                        line.contains(&format!(
                            "{{ read_group = {group}, error_probability_multiplier"
                        ))
                        .then_some(*value)
                    });
            match multiplier {
                Some(value) => {
                    rewritten.push_str(&line.replace("value = 1.0", &format!("value = {value}")))
                }
                None => rewritten.push_str(line),
            }
            rewritten.push('\n');
        }
        assert!(
            rewritten.contains("value = 0.25")
                && rewritten.contains("value = 2.5")
                && rewritten.contains("value = 4.0"),
            "all three read groups must have been given their own multiplier:\n{rewritten}",
        );
        std::fs::write(parameters, rewritten).expect("the parameters file rewrites");
    }

    /// [`called_from_alignments`], under a supplied parameters file rather than the defaults.
    fn called_from_alignments_under(
        cohort: &AVaryingCohort,
        output: &Path,
        parameters: &Path,
    ) -> Vec<String> {
        let mut args = alignments_args(cohort, output);
        args.parameters = Some(parameters.to_path_buf());
        args.defaults = false;
        run_call_from_alignments(&args).expect("the cohort calls under the supplied numbers");
        comparable(output)
    }

    /// [`called_from_psps`], under a supplied parameters file rather than the defaults.
    fn called_from_psps_under(
        cohort: &AVaryingCohort,
        output: &Path,
        samples: &[&str],
        parameters: &Path,
    ) -> Vec<String> {
        std::fs::remove_file(output).expect("the first route wrote the file this one replaces");
        let psps = walked_into_psps(cohort);
        let mut args = psps_args(cohort, output, &psps, samples);
        args.parameters = Some(parameters.to_path_buf());
        args.defaults = false;
        run_call_from_psps(&args).expect("the stored cohort calls under the supplied numbers");
        assert!(
            output.is_file(),
            "the stored cohort's run must have written it"
        );
        comparable(output)
    }

    // -----------------------------------------------------------------
    // The run-level invariances (spec §12.6, §12.7, §12.9)
    // -----------------------------------------------------------------

    /// **The order the psps are named in does not change the calls** (spec §12.6).
    ///
    /// **What it proves is that nothing joins by position.** Every per-sample number a run
    /// carries — the inbreeding coefficient, the read-group calibration, the column a genotype
    /// is written into — is keyed by the sample's name, and a run that joined any of them by the
    /// order its arguments arrived in would give a reordered cohort another sample's numbers
    /// with nothing failing. That is the failure this test exists for, and it is invisible in a
    /// VCF: every field would be in range and every column would have a genotype.
    ///
    /// **Sample for sample, not line for line.** The columns follow the order given, so the two
    /// files' text differs by construction; what must match is each sample's own fields at each
    /// locus.
    #[test]
    fn the_order_the_psps_are_named_in_does_not_change_the_calls() {
        let cohort = a_varying_cohort_on_disk();
        let output = cohort.directory.path().join("calls.vcf");

        called_from_alignments(&cohort, &output);
        let one_way = called_from_psps(&cohort, &output, &["one", "two"]);
        std::fs::remove_file(&output).expect("the file is there to remove");
        let psps = cohort.directory.path().join("psps");
        run_call_from_psps(&psps_args(&cohort, &output, &psps, &["two", "one"]))
            .expect("the same cohort in the other order calls");
        let the_other_way = comparable(&output);

        assert_ne!(
            one_way, the_other_way,
            "the two orders put the sample columns differently, so the files must differ — if \
             they did not, this test would be comparing one order with itself",
        );
        assert_eq!(
            genotypes_by_sample(&one_way),
            genotypes_by_sample(&the_other_way),
            "every sample must be called the same at every locus whichever order its file was \
             named in",
        );
    }

    /// Each locus's per-sample fields, keyed by the sample's name — what `--psp`'s order must
    /// not change.
    fn genotypes_by_sample(lines: &[String]) -> Vec<(String, Vec<(String, String)>)> {
        let samples: Vec<&str> = lines
            .iter()
            .find(|line| line.starts_with("#CHROM"))
            .expect("the VCF names its columns")
            .split('\t')
            .skip(9)
            .collect();
        lines
            .iter()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                let at = format!("{}:{}", columns[0], columns[1]);
                let mut called: Vec<(String, String)> = samples
                    .iter()
                    .zip(columns.iter().skip(9))
                    .map(|(sample, call)| ((*sample).to_owned(), (*call).to_owned()))
                    .collect();
                called.sort();
                (at, called)
            })
            .collect()
    }

    /// **A cohort each of whose samples was walked in its own invocation calls, and calls the
    /// same** (spec §12.7).
    ///
    /// **Without this the ordinary psp-mode run does not work at all.** A gatherer sees one
    /// sample's files, so it numbers that sample's read groups from zero and every sample's
    /// first group comes back as identifier `0`; the calling run has to merge those tables into
    /// one run-wide numbering (spec §6.1, §6.2). The failure without it is a refusal at the
    /// parameters fit rather than anything visible in the walk — and this cohort makes it
    /// reachable, because its first sample declares two read groups and its second one, so the
    /// separately-walked tables collide on `0` and disagree about what `1` means.
    ///
    /// **Compared against the one-invocation cohort**, not merely asserted to run: what must
    /// hold is that walking a cohort in pieces is the same run as walking it whole, which is the
    /// property psp mode exists for (spec §2).
    #[test]
    fn a_cohort_walked_one_sample_at_a_time_calls_what_one_invocation_calls() {
        let cohort = a_varying_cohort_on_disk();
        let output = cohort.directory.path().join("calls.vcf");

        called_from_alignments(&cohort, &output);
        let together = called_from_psps(&cohort, &output, &["one", "two"]);

        // One invocation a sample, each into a directory of its own — the way this command's own
        // documentation tells a person to spread a cohort over cores or machines.
        let mut apart = Vec::new();
        for (sample, alignment) in ["one", "two"].iter().zip(&cohort.alignments) {
            let directory = cohort.directory.path().join(format!("walked-{sample}"));
            run_generate_psps(&GeneratePspsArgs {
                reference: cohort.reference.clone(),
                catalog: Some(cohort.catalog.clone()),
                alignments: vec![alignment.clone()],
                output_dir: directory.clone(),
                regions: None,
                force: true,
                build_index_if_missing: false,
                min_copies: MinCopies::default(),
                min_period: DEFAULT_MIN_PERIOD,
                max_period: DEFAULT_MAX_PERIOD,
                max_str_len: DEFAULT_MAX_STR_LEN,
                min_purity: DEFAULT_MIN_PURITY,
            })
            .expect("one sample walks on its own");
            apart.push(psp_path_for(&directory, sample));
        }

        std::fs::remove_file(&output).expect("the file is there to remove");
        let mut args = psps_args(&cohort, &output, &cohort.directory.path().join("psps"), &[]);
        args.psps = apart;
        run_call_from_psps(&args).expect("separately-walked samples make a cohort");

        assert_eq!(
            comparable(&output),
            together,
            "a cohort walked one sample at a time must call what one invocation calls",
        );
    }

    /// **Ground a sample analysed and found nothing in reads back as analysed and empty**, not
    /// as ground it never looked at (spec §12.9).
    ///
    /// **A VCF cannot say which**, which is the whole reason a psp records its analysed regions:
    /// a sample with no reads over a stretch and a sample that never looked at it produce the
    /// same absence of records. So what this asserts is on the file rather than in it — the psp
    /// of a sample with no reads holds no record and its header still claims the whole ground,
    /// and the run calls that cohort rather than refusing it or dropping the sample.
    #[test]
    fn a_sample_that_analysed_ground_and_found_nothing_is_not_a_sample_that_never_looked() {
        use crate::ng::psp::PspReader;
        use crate::pop_var_caller_exp::test_fixtures::a_cohort_on_disk;

        // `zeta` carries three reads and `alpha` none, over ground both were walked across.
        let cohort = a_cohort_on_disk();
        let psps = cohort.directory.path().join("psps");
        run_generate_psps(&GeneratePspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output_dir: psps.clone(),
            regions: None,
            force: true,
            build_index_if_missing: false,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        })
        .expect("the cohort walks into psps");

        let mut empty = PspReader::open(&psp_path_for(&psps, "alpha")).expect("alpha's psp opens");
        let ground: Vec<_> = empty
            .header()
            .segmentation_inputs
            .analysed_regions
            .iter()
            .collect();
        assert!(
            !ground.is_empty(),
            "the file claims the ground it was walked over, which is what tells this apart from \
             a sample that never looked",
        );
        assert_eq!(
            empty.records().expect("the walk starts").count(),
            0,
            "and it holds no record, because the sample had no read there",
        );

        let output = cohort.directory.path().join("calls.vcf");
        let mut args = psps_args_for(&cohort.reference, &cohort.catalog, &output);
        args.psps = vec![psp_path_for(&psps, "zeta"), psp_path_for(&psps, "alpha")];
        run_call_from_psps(&args).expect("a cohort holding an analysed-but-empty sample calls");

        let called = comparable(&output);
        assert!(
            called
                .iter()
                .find(|line| line.starts_with("#CHROM"))
                .is_some_and(|line| line.contains("zeta") && line.contains("alpha")),
            "both samples keep their column, the empty one included:\n{}",
            called.join("\n"),
        );
    }

    /// [`psps_args`]'s shape for a cohort that is not [`AVaryingCohort`] — the psps are filled in
    /// by the caller.
    fn psps_args_for(reference: &Path, catalog: &Path, output: &Path) -> CallFromPspsArgs {
        CallFromPspsArgs {
            reference: reference.to_path_buf(),
            catalog: Some(catalog.to_path_buf()),
            psps: Vec::new(),
            output: output.to_path_buf(),
            parameters: None,
            defaults: true,
            ploidy: None,
            max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
            max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
            cohort_locus_builder_regions_len: None,
            threads: 0,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        }
    }
}
