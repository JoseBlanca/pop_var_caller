//! One locus reduced to one cell key: how many reads covered the site, and how many
//! of those showed something other than the reference.
//!
//! **The only place that decides what counts as an alternative read**, which is why
//! it is its own file rather than a method on the locus type — a locus cannot know a
//! model's answer to that question
//! (`arch/parameter_prepass_generic.md` §2.3).
//!
//! Implemented in Milestone C.
