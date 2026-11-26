
- CRAM major-version check missing on the live reader (found by D1; deleted cram_input had it). Add to read/input/open_bam.rs with a test: refuse non-3.x CRAM. Autonomous (unsupported input only).
