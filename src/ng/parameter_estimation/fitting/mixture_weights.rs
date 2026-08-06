//! The concave climb: given how likely each genotype makes each cell, find the
//! genotype frequencies that best explain the whole table.
//!
//! Shared by two of the four fits — the error-rate scan, which climbs to the best
//! frequencies at every rung of its ladder, and the sample's own rates, climbed once
//! on the whole-sample table. The surface is concave, so the climb cannot get stuck
//! on a false summit and a failure to converge is a bug rather than a data condition
//! (`spec/parameter_prepass.md` §3.1).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §4.1. Implemented in
//! Milestone D.
