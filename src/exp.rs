//! **Research switches for the GIAB alpha study (branch `research/giab-alpha`), not for main.**
//!
//! Each switch is an environment variable read once per process. Unset, every switch leaves the
//! caller exactly as main builds it, so a run without them is the baseline.
//!
//! **None is defined at present.** The first, `PVC_EXP_OFF_GRID_SHARE`, set the repeat-tract
//! periodicity check's share; it went with the check itself on 2026-10-22
//! (`doc/devel/ng/spec/candidate_alleles_ssr.md` §7). The module stays for the switches to come.
