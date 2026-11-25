use std::collections::HashMap;
use std::path::Path;

use super::*;

/// A file system of the given files, every other path unreadable.
fn files(entries: &[(&str, &str)]) -> impl Fn(&Path) -> Option<String> + use<> {
    let map: HashMap<String, String> = entries
        .iter()
        .map(|(path, text)| ((*path).to_owned(), (*text).to_owned()))
        .collect();
    move |path: &Path| map.get(path.to_str()?).cloned()
}

/// **A systemd unit's `MemoryMax=`** under cgroup v2: the unit's own `memory.max`, its slice
/// unlimited.
#[test]
fn a_systemd_units_memory_max_is_read_under_cgroup_v2() {
    let read = files(&[
        ("/proc/self/cgroup", "0::/system.slice/run-u42.service\n"),
        (
            "/sys/fs/cgroup/system.slice/run-u42.service/memory.max",
            "126701535232\n",
        ),
        ("/sys/fs/cgroup/system.slice/memory.max", "max\n"),
    ]);
    let (limit, file) = control_group_limit(&read).expect("a limit");
    assert_eq!(limit.bytes, 126_701_535_232);
    assert_eq!(limit.source, LimitSource::ControlGroupV2);
    assert_eq!(
        file,
        Path::new("/sys/fs/cgroup/system.slice/run-u42.service/memory.max")
    );
}

/// **A parent's limit binds its children**: the tighter of the two is the limit, wherever it is.
#[test]
fn the_tightest_limit_above_the_group_is_the_limit() {
    let read = files(&[
        ("/proc/self/cgroup", "0::/user.slice/job.scope\n"),
        ("/sys/fs/cgroup/user.slice/job.scope/memory.max", "max\n"),
        ("/sys/fs/cgroup/user.slice/memory.max", "8000000000\n"),
    ]);
    let (limit, _) = control_group_limit(&read).expect("a limit");
    assert_eq!(limit.bytes, 8_000_000_000);
}

/// **`max` everywhere is no limit**, and so is a group whose files cannot be read.
#[test]
fn no_limit_where_every_group_says_max_or_nothing_can_be_read() {
    let unlimited = files(&[
        ("/proc/self/cgroup", "0::/a/b\n"),
        ("/sys/fs/cgroup/a/b/memory.max", "max\n"),
        ("/sys/fs/cgroup/a/memory.max", "max\n"),
    ]);
    assert_eq!(control_group_limit(&unlimited), None);
    assert_eq!(control_group_limit(&files(&[])), None);
    let unreadable = files(&[("/proc/self/cgroup", "0::/a/b\n")]);
    assert_eq!(control_group_limit(&unreadable), None);
}

/// **A container with its own cgroup namespace** sees its group as `/`: the mount's own
/// `memory.max` is its limit.
#[test]
fn a_containers_root_group_limit_is_read() {
    let read = files(&[
        ("/proc/self/cgroup", "0::/\n"),
        ("/sys/fs/cgroup/memory.max", "17179869184\n"),
    ]);
    let (limit, _) = control_group_limit(&read).expect("a limit");
    assert_eq!(limit.bytes, 17_179_869_184);
}

/// **cgroup v1**: the memory controller's line, its `memory.limit_in_bytes`; the page-rounded
/// `i64::MAX` the kernel writes for no limit is none.
#[test]
fn cgroup_v1_limits_are_read_and_its_no_limit_value_is_none() {
    let membership = "12:pids:/user.slice\n\
                      4:cpu,cpuacct:/user.slice\n\
                      7:memory:/user.slice/job\n\
                      0::/user.slice/job\n";
    let read = files(&[
        ("/proc/self/cgroup", membership),
        (
            "/sys/fs/cgroup/memory/user.slice/job/memory.limit_in_bytes",
            "64000000000\n",
        ),
        (
            "/sys/fs/cgroup/memory/user.slice/memory.limit_in_bytes",
            "9223372036854771712\n",
        ),
    ]);
    let (limit, file) = control_group_limit(&read).expect("a limit");
    assert_eq!(limit.bytes, 64_000_000_000);
    assert_eq!(limit.source, LimitSource::ControlGroupV1);
    assert_eq!(
        file,
        Path::new("/sys/fs/cgroup/memory/user.slice/job/memory.limit_in_bytes")
    );
    let unlimited = files(&[
        ("/proc/self/cgroup", "7:memory:/\n"),
        (
            "/sys/fs/cgroup/memory/memory.limit_in_bytes",
            "9223372036854771712\n",
        ),
    ]);
    assert_eq!(control_group_limit(&unlimited), None);
}

/// `MemTotal` is read in kB and given in bytes; a file without it gives nothing.
#[test]
fn meminfo_total_is_read_in_bytes() {
    let text = "MemTotal:       131072000 kB\nMemFree:         1000 kB\n";
    assert_eq!(total_from_meminfo(text), Some(131_072_000 * 1024));
    assert_eq!(total_from_meminfo("MemFree: 1 kB\n"), None);
}

/// **The smaller limit wins**, whichever it is, and either may be missing.
#[test]
fn the_smaller_of_the_machine_and_the_group_is_the_limit() {
    let machine = MemoryLimit {
        bytes: 128_000_000_000,
        source: LimitSource::MachineMemory,
    };
    let group = MemoryLimit {
        bytes: 118_000_000_000,
        source: LimitSource::ControlGroupV2,
    };
    assert_eq!(smaller(Some(machine), Some(group)), Some(group));
    assert_eq!(smaller(Some(group), Some(machine)), Some(group));
    assert_eq!(smaller(Some(machine), None), Some(machine));
    assert_eq!(smaller(None, None), None);
}

/// **On the platforms pvc runs on, the limit is read**: a laptop or a server has a machine memory.
#[cfg(any(target_os = "linux", target_os = "macos"))]
#[test]
fn this_machines_limit_is_read() {
    let limit = of_this_process().expect("the machine's memory is readable here");
    assert!(limit.bytes > 100_000_000, "{limit:?}");
}
