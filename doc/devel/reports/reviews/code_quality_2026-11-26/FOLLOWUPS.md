
- CRAM major-version check missing on the live reader (found by D1; deleted cram_input had it). Add to read/input/open_bam.rs with a test: refuse non-3.x CRAM. Autonomous (unsupported input only). **FIXED (F8, branch cq-f8):** `read_header` in read/input/open_bam.rs refuses a CRAM whose major version is not 3 with `AlignmentFileError::UnsupportedCramVersion`; test `a_cram_with_major_version_4_is_refused_at_the_header_read`.
- RepeatLadder computes modal_repeat_count on every locus; only tests read it (found by D2). Candidate speed cleanup, gate timing.
- Release-only test failure: alignment::left_align_structured::tests::an_offset_past_the_reference_panics expects a debug_assert message (pre-existing; owner U8).
- Two rustdoc "redundant explicit link target" warnings: stratum_fits.rs:64, cohort_merge/close.rs:18.
