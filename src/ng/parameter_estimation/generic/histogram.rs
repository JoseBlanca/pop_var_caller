//! The cell table: a tally of what the sites looked like.
//!
//! Each covered position reduces to a small key — how many reads covered it, how many
//! of those showed something other than the reference, and, when a sample has more
//! than one read group and few alternative reads, which library those reads came
//! from. The table counts how many positions showed each key. Eight hundred million
//! positions become a few hundred counters, and the fits lose nothing, because a site
//! enters the likelihood only through that key
//! (`spec/parameter_prepass_generic.md` §4): the positions that all looked alike are
//! scored once and multiplied by how many there were.
//!
//! **The library attribution is the part that is easy to drop and cannot be.** With
//! it forgotten, a key of total depth and total alternative count sees only the
//! share-weighted mean error rate and nothing else about the individual libraries —
//! the likelihood is exactly flat along every combination holding that mean fixed, so
//! no amount of genome separates them. Keeping which library each of the first few
//! alternative reads came from is what breaks that flatness
//! (`arch/parameter_prepass_generic.md` §2.2).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §2.2. Implemented in
//! Milestone B; the depth ladder it is binned by lands in Milestone A.
