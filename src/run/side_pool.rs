//! **A small pool of threads beside the run's own**, for work that runs in bursts and must not
//! queue behind the long tasks on the run's pool: the hand-over's (`crate::run::psp_caller`) and
//! the reading of a round's ground (`crate::run::cohort_merge::stream`). Each caller says how many
//! threads and why; this holds what the two have in common.

/// A pool of its own, or none — its work then runs where it was asked for.
pub(crate) struct SidePool(Option<rayon::ThreadPool>);

impl SidePool {
    /// A pool of `threads` threads named `name-0`, `name-1`, …; none when `threads` is zero or the
    /// threads cannot be started.
    pub(crate) fn of(threads: usize, name: &'static str) -> Self {
        if threads == 0 {
            return Self(None);
        }
        Self(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .thread_name(move |at| format!("{name}-{at}"))
                .build()
                .ok(),
        )
    }

    /// Run `work` on the pool, or where it was called when there is none.
    pub(crate) fn run<R: Send>(&self, work: impl FnOnce() -> R + Send) -> R {
        match &self.0 {
            Some(pool) => pool.install(work),
            None => work(),
        }
    }
}
