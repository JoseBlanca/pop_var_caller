//! What a read at a repeat tract does to the length it shows, and the search that fits it.
//!
//! Three numbers describe slippage and they answer different questions: how often a read
//! shows a length other than its allele's (the level), which way it moves when it does
//! (strongly asymmetric — a read at a tomato dinucleotide is 4.9 times as likely to have
//! lost a repeat as to have gained one), and how far it moves (usually one copy, sometimes
//! two, with the same fall-off in both directions). The fourth number this path fits, the
//! per-base substitution rate, is not slippage at all and is a division rather than a
//! search (`spec/parameter_prepass_ssr.md` §3, §4.1).
//!
//! Empty until Milestone A4; the mathematics that uses these lands in Milestone D.
