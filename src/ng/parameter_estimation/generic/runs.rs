//! The inbreeding coefficient: a two-state hidden Markov model over windows, deciding
//! which of them lie in a **run of homozygosity** — a stretch of genome where both
//! copies descend from one recent ancestor, and which therefore carries almost no
//! heterozygotes. The rest of the genome carries them at the sample's own rate.
//!
//! The chain walks the windows of a contig deciding which state each is in, and the
//! coefficient is the share of the analysable genome the inside state claims —
//! weighted by reference positions rather than by loci, so a window dense in widened
//! indel loci is not under-weighted (`spec/parameter_prepass_generic.md` §6).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §5.3. Implemented in
//! Milestone E.
