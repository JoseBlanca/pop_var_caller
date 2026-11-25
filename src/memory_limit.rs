//! **The most memory this process may hold**: the smaller of the machine's memory and, on Linux,
//! the limit of the control group the process runs in.
//!
//! Physical memory alone is not the limit a run meets. A job started under a systemd unit with
//! `MemoryMax=`, or in a container, is killed when its control group reaches that limit, however
//! much memory the machine has. So on Linux both are read, and the smaller one is the limit:
//!
//! - **the machine's memory**: `MemTotal` in `/proc/meminfo` on Linux; `hw.memsize` on macOS, read
//!   by running `sysctl -n hw.memsize`, since asking the kernel directly needs `unsafe`, which this
//!   crate forbids;
//! - **the control group's limit** (Linux): cgroup v2's `memory.max`, or cgroup v1's
//!   `memory.limit_in_bytes`, for the group `/proc/self/cgroup` names and every group above it,
//!   since a parent's limit binds its children too. `max`, and v1's page-rounded `i64::MAX`, mean
//!   no limit.
//!
//! What cannot be read is left out. When neither can be read — another operating system, or a
//! `/proc` that is not mounted — there is no limit to check against, and the caller says so in
//! the log and goes on.

use std::path::{Path, PathBuf};

/// **The most memory this process may hold, and where that number came from.**
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MemoryLimit {
    /// The limit, in bytes.
    pub bytes: u64,
    /// Where it came from.
    pub source: LimitSource,
}

/// **Which limit was the smaller one.**
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LimitSource {
    /// The machine's physical memory.
    MachineMemory,
    /// The control group's `memory.max` (cgroup v2), such as a systemd unit's `MemoryMax=`.
    ControlGroupV2,
    /// The control group's `memory.limit_in_bytes` (cgroup v1).
    ControlGroupV1,
}

impl std::fmt::Display for LimitSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::MachineMemory => "the machine's memory",
            Self::ControlGroupV2 => "the memory.max of the control group pvc runs in",
            Self::ControlGroupV1 => "the memory.limit_in_bytes of the control group pvc runs in",
        })
    }
}

/// **This process's memory limit**: the smaller of the machine's memory and its control group's
/// limit, or `None` when neither can be read.
#[must_use]
pub fn of_this_process() -> Option<MemoryLimit> {
    let machine = machine_memory().map(|bytes| MemoryLimit {
        bytes,
        source: LimitSource::MachineMemory,
    });
    let group = control_group_limit(&|path| std::fs::read_to_string(path).ok());
    if let Some((limit, file)) = &group {
        tracing::debug!(
            "memory limit: the control group's, {} bytes, from {}",
            limit.bytes,
            file.display()
        );
    }
    smaller(machine, group.map(|(limit, _)| limit))
}

/// The smaller of two limits, either of which may be missing.
fn smaller(a: Option<MemoryLimit>, b: Option<MemoryLimit>) -> Option<MemoryLimit> {
    match (a, b) {
        (Some(a), Some(b)) => Some(if b.bytes < a.bytes { b } else { a }),
        (a, None) => a,
        (None, b) => b,
    }
}

/// The machine's physical memory, in bytes.
#[cfg(target_os = "linux")]
fn machine_memory() -> Option<u64> {
    total_from_meminfo(&std::fs::read_to_string("/proc/meminfo").ok()?)
}

/// The machine's physical memory, in bytes: `sysctl -n hw.memsize`.
#[cfg(target_os = "macos")]
fn machine_memory() -> Option<u64> {
    let output = std::process::Command::new("/usr/sbin/sysctl")
        .args(["-n", "hw.memsize"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    String::from_utf8(output.stdout).ok()?.trim().parse().ok()
}

/// Elsewhere the machine's memory is not read.
#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn machine_memory() -> Option<u64> {
    None
}

/// `MemTotal` from the text of `/proc/meminfo`, in bytes.
#[cfg_attr(not(target_os = "linux"), allow(dead_code))]
fn total_from_meminfo(text: &str) -> Option<u64> {
    let kib: u64 = text
        .lines()
        .find_map(|line| line.strip_prefix("MemTotal:"))?
        .trim()
        .strip_suffix("kB")?
        .trim()
        .parse()
        .ok()?;
    kib.checked_mul(1024)
}

/// **The tightest control-group limit on this process**, and the file it was read from: every
/// group from the one `/proc/self/cgroup` names up to the root, cgroup v2 and v1 alike. `read`
/// returns a file's text, or `None` where it cannot be read.
pub(crate) fn control_group_limit(
    read: &dyn Fn(&Path) -> Option<String>,
) -> Option<(MemoryLimit, PathBuf)> {
    let membership = read(Path::new("/proc/self/cgroup"))?;
    let mut tightest: Option<(MemoryLimit, PathBuf)> = None;
    for line in membership.lines() {
        // `hierarchy-ID:controller-list:cgroup-path`; the path may itself hold a `:`.
        let mut fields = line.splitn(3, ':');
        let (Some(_), Some(controllers), Some(group)) =
            (fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        let (mount, file, source) = if controllers.is_empty() {
            ("/sys/fs/cgroup", "memory.max", LimitSource::ControlGroupV2)
        } else if controllers
            .split(',')
            .any(|controller| controller == "memory")
        {
            (
                "/sys/fs/cgroup/memory",
                "memory.limit_in_bytes",
                LimitSource::ControlGroupV1,
            )
        } else {
            continue;
        };
        let mut directory = PathBuf::from(mount).join(group.trim_start_matches('/'));
        loop {
            let path = directory.join(file);
            let bytes = read(&path).and_then(|text| match source {
                LimitSource::ControlGroupV1 => limit_from_v1(&text),
                _ => limit_from_v2(&text),
            });
            if let Some(bytes) = bytes
                && tightest
                    .as_ref()
                    .is_none_or(|(limit, _)| bytes < limit.bytes)
            {
                tightest = Some((MemoryLimit { bytes, source }, path));
            }
            if directory == Path::new(mount) || !directory.pop() {
                break;
            }
        }
    }
    tightest
}

/// A cgroup v2 `memory.max`: a number of bytes, or `max` for none.
fn limit_from_v2(text: &str) -> Option<u64> {
    match text.trim() {
        "max" => None,
        number => number.parse().ok(),
    }
}

/// **A cgroup v1 `memory.limit_in_bytes`.** No limit is written as the largest page-rounded
/// `i64`, 9223372036854771712 on 4-kB pages; anything from 2⁶⁰ bytes up is read as none.
fn limit_from_v1(text: &str) -> Option<u64> {
    let bytes: u64 = text.trim().parse().ok()?;
    (bytes < 1 << 60).then_some(bytes)
}

#[cfg(test)]
mod tests;
