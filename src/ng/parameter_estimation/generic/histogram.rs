//! The cell table: a tally of what the sites looked like.
//!
//! Each covered position reduces to two numbers — how many reads covered it, and how
//! many of those showed something other than the reference — and the table counts how
//! many positions showed each combination. Eight hundred million positions become a
//! few hundred counters, and nothing is lost, because a site enters the likelihood
//! only through those two numbers: five sites that all looked like twelve reads with
//! one alternative are scored once and multiplied by five
//! (`spec/parameter_prepass_generic.md` §4).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §2.2. Implemented in
//! Milestone B; the depth ladder it is binned by lands in Milestone A.
