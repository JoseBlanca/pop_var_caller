# Floating-point maths that rounds alike on macOS and Linux

*Status: draft, 2026-09-13. Branch `portable-float`, worktree `../pop_var_caller-portable-float`.*

This plan turns the owner's brief of 2026-09-13 into build order. It is **not** a place for new
design. The brief settles the goal and a candidate fix, and asks for the fix's two costs to be
measured before it is adopted; Milestone A is that measurement, and Milestones B and C run only if
the owner approves at the checkpoint before each.

## 1. The problem and the candidate fix

Rust's `f64::ln`, `exp`, `powf`, `log10`, `ln_1p`, `exp_m1` and the other transcendental functions
call the operating system's maths library, and IEEE 754 does not require those to round correctly.
glibc on Linux and Apple's libm on macOS can return values one or two units apart in the last
binary place. On 2026-09-13 that difference flipped a tie in the repeat-tract delimiter: at Phred
quality 4 the aligner's match score was `0xbfe03ee1794fdd29` on macOS and `0xbfe03ee1794fdd27` on
glibc, and one generated read's tract measured a base longer on macOS (test
`alignment::delimit_parity`). Commit `75722336` fixed that table alone, by writing its 256 values
out as the bits glibc produced (`src/alignment/emission.rs`).

**The candidate fix.** Route every such call through one module, `crate::float`, whose functions
delegate to the `libm` crate — a pure-Rust port of musl's maths library, already a dependency for
`lgamma` (`src/genetics.rs`). A clippy `disallowed_methods` rule stops new code calling the std
methods directly. Because the same Rust source runs on both platforms, the same inputs give the
same bits.

**Two costs are unmeasured**, and measuring them is Milestone A:

- **Speed.** glibc's `ln` and `exp` are heavily optimised; libm's may be slower, and the caller
  calls them in its read-likelihood, frequency and fitting loops.
- **Output.** libm and glibc do not round alike either, so Linux's own results move once, in the
  last bits, and a tie can break the other way.

## 2. Scope

**In**

- Milestone A: an inventory of every transcendental call in `src/`; a measurement of std against
  libm per function, for speed and for how often the bits differ, on both platforms; and the
  caller's own speed and memory before any change. Nothing shipped changes.
- Milestone B (if approved): the `float` module, the call-site conversion, the clippy ban, and the
  same measurements repeated.
- Milestone C (if approved): a macOS CI job, a cross-platform digest test, a new identity baseline,
  and a PROJECT_STATUS entry.

**Out**

- Examples (`examples/`, 343 calls) are research tools, not the caller. Whether the clippy ban
  reaches them is decided at Checkpoint A with the inventory in hand; converting them is not in
  this plan either way.
- Floating-point *formatting and parsing* (`{}` in the VCF writer, `str::parse::<f64>`) and the
  basic operations `+ − × ÷` and `sqrt` are IEEE-exact or implemented in Rust's own core library,
  so they are already identical across platforms. The inventory confirms nothing else in the
  arithmetic depends on the platform (fused multiply-add, x87 extended precision).

## 3. Principles that fixed the order

1. **Measure before changing anything.** The owner will decide from Milestone A whether to convert
   everything, only the calls outside hot loops, or nothing. No step of A edits `src/`.
2. **Baselines come from the tree the change will be compared against.** The caller's benches and
   runs are taken on this worktree at `75722336`, on the same machines and with the same commands
   Milestone B will repeat, and each is repeated enough times to state its noise.
3. **An oracle confirmed before it is trusted.** The identity oracle's digests are reproduced on
   the untouched worktree before any later difference is read as caused by a change.
4. **A step that moves output silently is its own commit.** In Milestone B each conversion group
   lands alone, with the oracle run before and after, so a moved genotype can be bisected to one
   group.
5. **Reuse over rewrite.** The digest oracle is `scripts/promote_ng_oracle.sh`, unchanged; the
   benches are the five existing criterion benches, unchanged.

## 4. Preconditions (checked before step A1)

- The worktree is at `75722336` with a clean tree.
- A release build in the container reproduces the five identity-oracle digests the brief lists
  (`b74f3edf…` twice, `5e6d874d…`, `25d388c8…`, `259baf04…`).
- `cargo` works natively on the macOS host, and `./scripts/dev.sh` reaches the container.

## 5. The steps

The machines: **macOS** is the host (Apple M5 Pro, 18 cores, 64 GB), running natively; **Linux**
is the Apple `container` VM on the same host (arm64 glibc, 8 CPUs, 16 GB). Linux on x86_64 — the
CI runner and the Linux dev box — is not measured here; its glibc picks a different `exp`/`log`
build at load time (with or without fused multiply-add), which the report must say it did not
cover.

### Milestone A — measure; change nothing shipped

- ✅ **A1. The call-site inventory.** A report listing every transcendental call in `src/` outside
  tests: function, file and line, the function it sits in, and how often it runs — per base, per
  read, per locus, per fitting iteration, once per run, once per table — and whether its argument
  is a constant computed from literals (writable as bits at no speed cost, as `75722336` did). It
  also covers `powi`, `sqrt`, the trigonometric functions, `libm::lgamma`, the `wide` SIMD types,
  and any `mul_add`. Hot-loop classification is checked against a sampling profile of one run of the
  caller where the code reading is not decisive.
  *Depends:* preconditions. *Source:* brief, "Milestone A", first bullet.
- ☐ **A2. std against libm, per function.** An example program, `examples/float_libm_vs_std.rs`,
  that for each function the inventory found evaluates std and libm over the input ranges the
  inventory names, (a) timing both in repeated trials and reporting the median and spread, (b)
  counting inputs where the two return different bits, and (c) writing each side's outputs to a
  file so the macOS and Linux runs can be compared bit for bit — showing how often the two
  platforms' std disagree today, and that libm's outputs match across them. Run on both machines.
  *Depends:* A1. *Source:* brief, "Milestone A", second bullet.
- ☐ **A3. The caller's speed baseline.** The five criterion benches (`ng_site_quality_perf`,
  `ng_ssr_delimiter_perf`, `ng_generic_pileup_perf`, `ng_psp_perf`, `ng_joint_fit_perf` with
  `--features bench-fixtures`) on both machines, each at its native parallelism and run twice to
  see run-to-run noise; and wall time and peak memory of `call-from-alignments`, `generate-psps`
  and `call-from-psps` on the four tomato accessions, repeated, on both machines. A script,
  `scripts/portable_float_baseline.sh`, holds the exact commands so Milestone B repeats them.
  *Depends:* preconditions. *Source:* brief, "Milestone A", third bullet.

> **Checkpoint A: the costs are known.** Report the per-function slowdown on each platform, where
> the hot calls are, and a recommendation: convert everything, convert only the cold calls and pin
> the constants, or drop the idea. Pause for the owner's decision.

### Milestone B — convert (only if approved at A)

Provisional; its grouping is revised at Checkpoint A from A1's inventory.

- ☐ **B1. The `float` module.** `src/float.rs`: one function per transcendental the caller uses,
  each delegating to `libm`, plus `powi` as an explicit multiplication loop if A2 shows std's `powi`
  is not bit-stable; unit tests pin a handful of outputs as bits.
  *Depends:* A. *Source:* brief, fix step 1.
- ☐ **B2…Bn. Convert call sites, one group per commit.** Groups follow the module tree (alignment;
  calling likelihoods; calling inference, prior and quality; parameter estimation; the rest), with
  constants computed from literals written as bits where A1 marked them. **Own commit each, do not
  bundle**; the identity oracle runs before and after each, and each commit message records whether
  the digests moved.
  *Depends:* B1. *Source:* brief, fix step 1.
- ☐ **B(n+1). The clippy ban.** `clippy.toml` `disallowed-methods` on the std methods, with the
  scope over examples as decided at Checkpoint A.
  *Depends:* all conversion steps. *Source:* brief, fix step 2.
- ☐ **B(n+2). Measure the change.** A3's benches and runs repeated on the same machines; the identity
  oracle's VCFs diffed against A's, counting changed records, genotypes and repeat-tract calls and
  explaining each; the GIAB benchmark if the output moved; the full test suite on both platforms.
  *Depends:* B(n+1). *Source:* brief, "Milestone B".

> **Checkpoint B: the cost is paid or refused.** Report the speed cost, the output change and a
> recommendation. Pause for the owner's decision.

### Milestone C — lock it in (only if approved at B)

- ☐ **C1.** A macOS job in `.github/workflows/ci.yml` running the test suite. *Source:* fix step 3.
- ☐ **C2.** A cross-platform end-to-end test: a synthetic cohort called, the VCF's digest asserted.
  *Source:* fix step 4.
- ☐ **C3.** A new identity baseline and a PROJECT_STATUS entry. *Source:* brief, "Milestone C".

> **Checkpoint C: done.** Pause for review.

## 6. Verification summary

| milestone | proven by |
|---|---|
| A | the identity oracle reproducing the brief's five digests on the untouched tree; every figure in the reports traceable to a log under `tmp/` and the command that made it |
| B | the test suite green on Linux and macOS; A2's libm output files identical across the two machines; the oracle's VCF diff explained call by call; benches against A3 with noise stated |
| C | the macOS CI job green; the digest test passing on both platforms |

## 7. Gates for every commit

In the container: `cargo fmt --all -- --check`; `cargo clippy --all-targets --all-features -- -D
warnings`; `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features`; `cargo test --lib
--tests --all-features` (4,787 passed, 0 failed, 4 ignored at `75722336`); the release build; the
identity oracle. On macOS: `cargo test --lib --tests --all-features` (4,786 passed). `cargo bench
--no-run` takes `-j 2` in the container.

## 8. Deviations

None yet.
