//! The comparator: Hardy–Weinberg at one estimated allele frequency, substituted in as
//! though it were the truth.
//!
//! **Empty until plan step F** (`doc/devel/ng/impl_plan/calling_prior.md`). It is the route
//! this caller does *not* take, kept behind the same seam only so the change the default
//! makes stays measurable — never a shipping default.
//!
//! Estimating the frequency and squaring the estimate undercounts homozygotes by exactly
//! the variance of that frequency, and that variance is how badly the frequency is pinned
//! down — so the gap is negligible with a thousand samples and largest at one sample and
//! low depth, the corner this caller commits to supporting. Measured on the GIAB trio,
//! each sample called on its own at 5×, swapping this route for the default took SNP
//! genotype accuracy at true variants from **83.6% to 94.6%**, and the sites where a
//! sample carrying two copies of the variant was called heterozygous from **214 to 8**,
//! with the emitted variant set byte-identical (`doc/devel/ng/spec/calling_priors.md`
//! §2.2). That is one corner — one sample at a time, high-quality human data, 5× — and it
//! is the corner the change was aimed at.
