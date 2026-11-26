//! **How many files this process may hold open, and what to tell a person whose run needs more.**
//!
//! A run that opens every sample's inputs at once and finds the limit too low fails part-way, with
//! "Too many open files", after the work before it is lost. On rick, 157 CRAMs failed so under
//! Linux's default soft limit of 1,024 (2026-11-25), where `ulimit -n 65536` let the same run
//! finish. So each command that holds one or more files a sample counts what it will hold, and
//! refuses at start-up, before reading any data, if that is above the soft limit. **It never
//! raises the limit itself**: the operator decides, and the message says how.
//!
//! The limit is read with `getrlimit(RLIMIT_NOFILE)` through `rustix`, which needs no `unsafe`
//! here. *Soft* is what the process may open now; *hard* is the most the soft limit may be raised
//! to without an administrator. Either may be unlimited.
//!
//! The counts themselves live with the commands: `src/run/callers.rs` for alignment files,
//! [`DESCRIPTORS_A_PSP_RUN_NEEDS_BESIDES_ITS_PSPS`] for psp files.

/// **What a psp-reading command holds besides its psps**: one descriptor a psp for the whole run
/// (measured, `call-from-psps` and `estimate-parameters` on 9 and 29 tomato1 psps at 1 and 8
/// threads: exactly one a psp), and 5 to 7 others — the three standard streams, the log, the
/// reference, the scratch directory's spill. The rest of the 32 is slack, for the psp prefetcher's
/// short-lived opens, the output and its index, the catalog and whatever the runtime holds.
pub const DESCRIPTORS_A_PSP_RUN_NEEDS_BESIDES_ITS_PSPS: u64 = 32;

/// **This process's limit on open files**: `None` in a field means unlimited.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DescriptorLimit {
    /// What the process may open now (`ulimit -n`, the soft limit).
    pub soft: Option<u64>,
    /// The most the soft limit may be raised to without an administrator (`ulimit -Hn`).
    pub hard: Option<u64>,
}

impl DescriptorLimit {
    /// The limit this process runs under.
    #[must_use]
    pub fn of_this_process() -> Self {
        let limit = rustix::process::getrlimit(rustix::process::Resource::Nofile);
        Self {
            soft: limit.current,
            hard: limit.maximum,
        }
    }
}

/// **A run needs more open files than its soft limit allows**: how many, the limit, and what to do,
/// for the command's error to carry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DescriptorShortfall {
    /// What the run needs at once.
    pub needed: u64,
    /// The soft limit it is above.
    pub soft: u64,
    /// The hard limit; `None` for unlimited.
    pub hard: Option<u64>,
    /// The samples that fit under the hard limit, where the hard limit is below the need too.
    pub samples_that_fit: Option<usize>,
}

impl DescriptorShortfall {
    /// **Whether a run needing `needed` fits `limit`**, and if not, the shortfall. `samples_at`
    /// gives what a run of that many of the cohort's samples would need; it finds how many fit
    /// under the hard limit.
    #[must_use]
    pub fn judge(
        needed: u64,
        limit: DescriptorLimit,
        samples: usize,
        samples_at: impl Fn(usize) -> u64,
    ) -> Option<Self> {
        let soft = limit.soft?;
        if needed <= soft {
            return None;
        }
        let samples_that_fit = limit.hard.filter(|&hard| hard < needed).map(|hard| {
            (0..samples)
                .rev()
                .find(|&fewer| samples_at(fewer) <= hard)
                .unwrap_or(0)
        });
        Some(Self {
            needed,
            soft,
            hard: limit.hard,
            samples_that_fit,
        })
    }

    /// **The soft limit to suggest**: a quarter above the need, rounded up to a power of two, and
    /// no higher than the hard limit. `None` where the hard limit is below the need.
    #[must_use]
    pub fn suggested_soft_limit(&self) -> Option<u64> {
        let roomy = (self.needed + self.needed / 4)
            .next_power_of_two()
            .max(1_024);
        match self.hard {
            None => Some(roomy),
            Some(hard) if hard >= roomy => Some(roomy),
            Some(hard) if hard >= self.needed => Some(hard),
            Some(_) => None,
        }
    }

    /// **What the person can do, in words**: the command that raises the soft limit, or, where the
    /// hard limit is too low as well, the administrator or a smaller cohort.
    #[must_use]
    pub fn remedy(&self) -> String {
        let hard = self
            .hard
            .map_or_else(|| "unlimited".to_owned(), |hard| hard.to_string());
        match self.suggested_soft_limit() {
            Some(raise) => format!(
                "Raise the soft limit in the shell that starts pvc, then run it again: \
                 `ulimit -n {raise}` (the hard limit, {hard}, allows it). In a systemd unit, set \
                 LimitNOFILE={raise}."
            ),
            None => {
                let fewer = match self.samples_that_fit {
                    Some(0) | None => String::new(),
                    Some(fit) => format!(", or call at most {fit} samples at once"),
                };
                format!(
                    "The hard limit, {hard}, is below the need too, so the soft limit cannot be \
                     raised far enough without the administrator: ask them to raise the hard \
                     limit to at least {} (`nofile` in /etc/security/limits.conf, or LimitNOFILE= \
                     in a systemd unit){fewer}.",
                    self.needed
                )
            }
        }
    }
}

/// **A psp-reading run needs more open files than its soft limit allows.** Each psp stays open
/// for the whole run, one descriptor each, whatever the thread count.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error(
    "this run reads {psps} psp files and needs {} files open at once, one a psp and {} besides, \
     and this process may open {} (its soft limit, `ulimit -n`). Every psp stays open for the \
     whole run, so the thread count does not change this. {}",
    shortfall.needed,
    DESCRIPTORS_A_PSP_RUN_NEEDS_BESIDES_ITS_PSPS,
    shortfall.soft,
    shortfall.remedy()
)]
pub struct PspsOverTheLimit {
    /// The psps the run reads.
    pub psps: usize,
    /// What the run needs, and the limits.
    pub shortfall: DescriptorShortfall,
}

/// **Whether a run reading `psps` psp files fits `limit`**, checked before any is opened.
///
/// # Errors
///
/// [`PspsOverTheLimit`] where the need is above the soft limit.
pub fn refuse_more_psps_than_the_limit(
    psps: usize,
    limit: DescriptorLimit,
) -> Result<(), PspsOverTheLimit> {
    let at = |psps: usize| psps as u64 + DESCRIPTORS_A_PSP_RUN_NEEDS_BESIDES_ITS_PSPS;
    match DescriptorShortfall::judge(at(psps), limit, psps, at) {
        None => Ok(()),
        Some(shortfall) => Err(PspsOverTheLimit { psps, shortfall }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two descriptors a sample and 32 besides.
    fn cohort(samples: usize) -> u64 {
        samples as u64 * 2 + 32
    }

    /// **Within the soft limit, or with no limit, a run goes on.**
    #[test]
    fn a_run_within_the_soft_limit_goes_on() {
        let limit = DescriptorLimit {
            soft: Some(1_024),
            hard: Some(4_096),
        };
        assert_eq!(DescriptorShortfall::judge(1_024, limit, 496, cohort), None);
        let unlimited = DescriptorLimit {
            soft: None,
            hard: None,
        };
        assert_eq!(
            DescriptorShortfall::judge(u64::MAX, unlimited, 1, cohort),
            None
        );
    }

    /// **rick's failure**: 157 CRAMs at 15 descriptors each and 32 besides, 2,387, under a soft
    /// limit of 1,024 and a hard one of 1,048,576. The suggestion is a quarter above the need,
    /// rounded up to a power of two: 4,096.
    #[test]
    fn above_the_soft_limit_the_suggestion_is_a_roomy_power_of_two() {
        let needed = 157 * 15 + 32;
        let shortfall = DescriptorShortfall::judge(
            needed,
            DescriptorLimit {
                soft: Some(1_024),
                hard: Some(1_048_576),
            },
            157,
            |samples| samples as u64 * 15 + 32,
        )
        .expect("2,387 is over 1,024");
        assert_eq!(shortfall.needed, 2_387);
        assert_eq!(shortfall.samples_that_fit, None);
        assert_eq!(shortfall.suggested_soft_limit(), Some(4_096));
        assert_eq!(
            shortfall.remedy(),
            "Raise the soft limit in the shell that starts pvc, then run it again: `ulimit -n \
             4096` (the hard limit, 1048576, allows it). In a systemd unit, set LimitNOFILE=4096."
        );
    }

    /// **A hard limit between the need and the roomy value**: the hard limit is suggested.
    #[test]
    fn a_hard_limit_just_above_the_need_is_suggested_itself() {
        let shortfall = DescriptorShortfall::judge(
            2_387,
            DescriptorLimit {
                soft: Some(1_024),
                hard: Some(3_000),
            },
            157,
            |samples| samples as u64 * 15 + 32,
        )
        .expect("over");
        assert_eq!(shortfall.suggested_soft_limit(), Some(3_000));
    }

    /// **A hard limit below the need**: the administrator, or fewer samples — the most that fit
    /// under the hard limit, (4,096 − 32) / 15 = 270 of 400.
    #[test]
    fn a_hard_limit_below_the_need_names_the_administrator_and_the_samples_that_fit() {
        let per_sample = |samples: usize| samples as u64 * 15 + 32;
        let shortfall = DescriptorShortfall::judge(
            per_sample(400),
            DescriptorLimit {
                soft: Some(1_024),
                hard: Some(4_096),
            },
            400,
            per_sample,
        )
        .expect("6,032 is over 1,024");
        assert_eq!(shortfall.samples_that_fit, Some(270));
        assert_eq!(shortfall.suggested_soft_limit(), None);
        let remedy = shortfall.remedy();
        for expected in [
            "The hard limit, 4096, is below the need too",
            "raise the hard limit to at least 6032",
            "or call at most 270 samples at once.",
        ] {
            assert!(remedy.contains(expected), "{expected}: {remedy}");
        }
    }

    /// **A psp run needs one descriptor a psp and 32 besides**: 2,176 psps under a soft limit of
    /// 1,024 are refused, with the arithmetic and the command; 992 fit exactly.
    #[test]
    fn a_psp_run_is_refused_above_the_soft_limit_naming_the_arithmetic() {
        let limit = DescriptorLimit {
            soft: Some(1_024),
            hard: Some(524_288),
        };
        refuse_more_psps_than_the_limit(992, limit).expect("992 + 32 is 1,024");
        let refused = refuse_more_psps_than_the_limit(2_176, limit).expect_err("2,208 is over");
        assert_eq!(
            refused.to_string(),
            "this run reads 2176 psp files and needs 2208 files open at once, one a psp and 32 \
             besides, and this process may open 1024 (its soft limit, `ulimit -n`). Every psp \
             stays open for the whole run, so the thread count does not change this. Raise the \
             soft limit in the shell that starts pvc, then run it again: `ulimit -n 4096` (the \
             hard limit, 524288, allows it). In a systemd unit, set LimitNOFILE=4096."
        );
    }
}
