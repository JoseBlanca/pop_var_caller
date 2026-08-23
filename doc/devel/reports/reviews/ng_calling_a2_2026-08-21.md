# Code review — ng calling foundations, step A2 (`Genotype`)

*2026-08-21. Branch `ng-calling-foundations`, reviewed at `2cf8be6e` (step A1) with the step's
uncommitted working-tree diff applied. Three category sub-agents, each in its own git worktree.
Per-category audit trail in the gitignored `tmp/review_2026-08-21_ng-calling-a2/`.*

## 1. Scope

**What was reviewed:** the working-tree diff of step A2 of
[`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md) — a single-file, purely
additive change to [`src/ng/types.rs`](../../../../src/ng/types.rs), 128 insertions and 0
deletions at the time of review. One type, `Genotype`, in a new section at the end of the file,
plus three tests.

**Deliberately out of scope:** the rest of the repository, and `AlleleId` and `Phred`, reviewed
and committed as step A1.

**Categories dispatched.** Three agents for 128 additive lines defining one type with no callers:

| category | reason |
|---|---|
| `reliability` | always; and the only category that mutation-tests. Also asked directly whether `new` should refuse an empty vector |
| `naming` | always; the diff is vocabulary, and its doc comment argues two deliberate omissions |
| `idiomatic` + `smells` + `defaults` | **one agent for three**, since a 128-line diff adding one type does not warrant three fan-out slots. Findings kept under separate headings |

Not dispatched: `errors` (the type is infallible and adds no error path — the one error-shaped
question, whether an empty multiset should be a `Result`, was routed to `reliability` instead
because it turns on this file's constructor conventions), `module_structure` (one file, nothing
moved), `unsafe_concurrency` (none present; the crate `forbid`s `unsafe_code`), `tooling`
(`Cargo.toml` untouched), `extras`, `refactor_safety` (purely additive, zero callers).

## 2. Verdict

**Approve-with-changes.** No Blockers. Two Majors, both from `reliability` and both applied. The
mutation pass is again the reason this is not an Approve: two mutations survived all three
submitted tests, and both broke the exact property the type exists to provide.

## 3. Execution status

Orchestrator runs, passed verbatim into every sub-agent's prompt:

| command | exit | result |
|---|---|---|
| `cargo fmt` | 0 | clean |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 4.46s` |
| `cargo test --lib ng::types` | 0 | `44 passed; 0 failed; 0 ignored; 0 measured; 3915 filtered out` |

`reliability` built and mutated in its own worktree; the other two found no finding that needed a
build and said so rather than spending the time.

Findings labelled "Needs verification": **zero**.

## 4. Open questions and assumptions

1. **Should `Genotype::new` refuse an empty vector, and at what price?** Two agents reached
   opposite recommendations, and the disagreement is resolved in §6 under M2. `reliability`'s
   assumption, checked: `grep -rn "Genotype::new" src/ | grep -v src/ng/types.rs` returns nothing,
   so outside its own tests the constructor has no callers — the minter is Milestone C.
2. **Is an id at or above 256 reachable?** No, through today's caps: `AlleleId`'s own doc records
   that production keeps 6 alleles per record and refuses above 16. The agent flagged this against
   itself, which is why the wide-id mutation is the *weaker* half of M1's evidence and the
   three-copy one the stronger.

## 5. Top 3 priorities

1. **M1 — a wrong sort survives all three tests from three copies up.** Every fixture used at most
   four copies and ids `0..=3`; an ordinary "already sorted?" fast path reproduces them all and
   silently makes two spellings of one genotype compare unequal.
2. **M2 — `new(vec![])` succeeded**, while `Ploidy::try_new(0)` in the same file refuses, and the
   type's doc says the genotype's length *is* the ploidy.
3. **Mi1 — a test named `..._at_any_ploidy` that exercises two ploidies**, a wording this file has
   already been bitten by once and records against itself.

## 6. Findings

### Major

**M1: `src/ng/types.rs` (`Genotype::new`) — the algebraic law the type rests on had no property
test, and two realistic wrong sorts survived every fixture.** **Category:** reliability.
**Confidence:** High — both survivors were executed and their behaviour change printed.

`Genotype::new` is a pure function whose correctness *is* a law: `new(v) == new(any permutation of
v)`. Every fixture in the three submitted tests used ids in `0..=3` and at most four copies, and
two mutations reproduce all of them:

| mutation | why it is realistic | what it breaks |
|---|---|---|
| `if alleles.first() > alleles.last() { sort_unstable() }` | an ordinary "is it already sorted?" fast path | leaves interior disorder from **three** copies up — and triploid and tetraploid regions are in scope by this project's own commitment, which `Ploidy`'s doc states |
| `sort_unstable_by_key(\|a\| a.0 as u8)` | a narrowed sort key | misorders any id at or above 256 — a hazard against the type's stated contract, not against today's data |

Measured by probe, under the fast path at three copies:

```
PROBE3 a: [AlleleId(0), AlleleId(2), AlleleId(1)] b: [AlleleId(0), AlleleId(1), AlleleId(2)] equal=false
```

Unmutated, `equal=true`. So under either survivor two spellings of one multiset compare **unequal**
— nothing panics, and the wrong count reaches the cohort allele frequencies, exactly as the
submitted test's own doc comment says ("one heterozygote in a cohort would count as two").

The agent also found a **free second fix**: the tetraploid fixture `[3, 0, 2, 0]` begins with its
largest id and ends with its smallest, which is the one arrangement that trips a "reversed?"
heuristic into sorting anyway — which is precisely why it survived that fixture. Swapping two
entries so the fixture starts and ends in order kills the mutation with no new test at all.

**M2: `src/ng/types.rs` (`Genotype::new`) — `new(vec![])` succeeded, and the invariant that
forbids it was prose with nothing behind it.** **Category:** reliability (raised as a Major),
`idiomatic` (raised independently as a Major), `naming` (raised as a cross-category note).
**Confidence:** High for the behaviour, executed:

```
PROBE-EMPTY constructed=Genotype([]) len=0 alleles=[]
PROBE-EMPTY equals-another-empty=true ploidy_try_new(0)=Err(Ploidy(0))
```

Two representations of one quantity with two different domains. The type's own doc says "how many
alleles it holds **is** the ploidy at that region", `Ploidy::try_new(0)` refuses because a genome
with no copies is not a genome, and `Genotype` accepted it in silence. The submitted doc comment's
defence — "a failure the only minter cannot produce" — is a claim about a component that has not
been written yet, made about a constructor that is `pub` on ng's surface.

The concrete failure, at Milestone C: expanding a genotype-table row means walking a flat
row-major counts row and pushing one `AlleleId` per count. If a shape is ever built with an allele
count of zero — an empty candidate table, an off-by-one on a row stride, a prune that emptied a
locus — the walk pushes nothing, `new` accepts it, and the empty genotype travels through
`SampleGenotypeCall` to `to_vcf` and out as a `GT` field naming no allele: a sample declared to
have no genome at that locus, written to the VCF, with nothing between the table and the writer
objecting and nothing panicking.

**The two agents recommended opposite fixes, and the disagreement is the useful part.**
`idiomatic` wanted `try_new -> Result` with an `EmptyGenotype` variant, on the grounds that the
file rules ploidy zero illegal in one type and legal in another. `reliability` argued against a
`Result` from this file's own convention: `types.rs` draws a deliberate line between *constrained*
newtypes wrapping untrusted external scalars (`Phred`, `Ploidy`, `ErrorRate`, `GenotypeFrequency`,
`InbreedingF` — all `try_new`, all with a `DomainError` variant) and *unconstrained* ones wrapping
internal indices (`ContigId`, `Position`, `AlleleId`), and `AlleleId`'s doc states that policy
outright. A `Genotype` is a container of internal indices, not a parsed external value; a `Result`
on a path that cannot fail is discharged with `.expect()` at the one call site, which relocates
the panic without adding a guarantee, and costs a `DomainError` variant whose message nothing will
read.

`reliability`'s middle answer — `assert!` plus a `#[should_panic]` test — is what was applied. It
makes the state unrepresentable, which was `idiomatic`'s actual complaint, at the price of one
length comparison per sample per locus on the last pass only, which is not the loop's hot path.

### Minor

- **Mi1 (reliability): `a_genotype_holds_one_allele_per_genome_copy_at_any_ploidy` promised a
  universal two fixtures cannot deliver**, and this file already records the same lesson against
  itself for `Ploidy` ("Named `every` and checking three of 255 was the gap"). The difference is
  that `Ploidy`'s domain is a finite `u8` and can be enumerated, while a genotype's length domain
  cannot — so **no fixture can establish "no ceiling"**, and that limit belongs stated on the test.
- **Mi2 (smells): the constructor's stated warrant for taking a `Vec` was not quite true.**
  "`into_boxed_slice` on an exactly-sized vector hands the buffer straight over: no allele is
  copied twice" holds only for callers who size exactly; `into_boxed_slice` shrinks first, so the
  ordinary `Vec::new()` plus one `push` per genome copy (length 2, capacity 4 for a diploid)
  reallocates and copies. The cost is a handful of `u16`s; the misleading guarantee is the defect.
- **Mi3 (smells): the derived `Ord` had no stated meaning**, against this file's own habit —
  `Ploidy` and `DomainError` both explain their derives — and its behaviour at mixed ploidy
  (`[0]` before `[0, 0]`) is not what a reader predicts from the word "multiset".

### Nits

The `HashSet` and `cmp` assertions in the multiset test cannot fail once the `alleles()`
comparison above them has passed, because `Genotype` has one field and all three impls are derived
— worth keeping as guards against a future second field, worth a half-line saying they are not
coverage today. The `HashSet` was named `hashes` while holding genotypes. *Pass* was undefined in
the sentence justifying the whole shape of the type. The "transposition hole" sentence named no
pair at risk, and `Ploidy`'s own doc gives division-by-zero rather than transposition as why it
exists. `std::cmp::Ordering::Equal` was fully qualified where the file's other ordering test
imports it.

Argued down and not filed, each with the reason recorded: `new` is the right name (the file
reserves `try_new` for constructors that *reject*, and `Motif::new` is already a plain — and
fallible — `new` here); the sort cannot surprise a caller, because `new` takes the `Vec` by value
so nobody can observe their own vector reordered, and `from_unsorted` would imply a `from_sorted`
sibling that will never exist; `alleles()` matches `Motif::as_bytes()`, since the file's `get()` is
for `Copy` scalars returning a primitive by value; the banner is form-identical to the motif
section's; `sort_unstable` is sound and `sort` would allocate; the derive set is the file's rule
minus the `Copy` a heap-owning type cannot have; nothing is dead code; the doc's length is one
paragraph per decision; not deriving `Default` is correct and must stay so, since it would mint
exactly the empty genotype M2 is about; and `impl IntoIterator` over `Vec` buys generality but no
allocation (`&[AlleleId]` would be the worst of the three).

## 7. Out of scope observations

- **`src/ng/mod.rs` has a curated `pub use types::{…}` list that names none of A1's or A2's three
  new types.** It also omits `LogProb`, `Position`, `GenomeRegion` and `ReadGroupId`, so it is a
  partial convenience list rather than a contract, and adding to it is step B1's business when it
  wires `calling/` in. Recorded, not changed.
- **`arch/ng_step_interfaces.md` §5 still records the decision that a `Genotype` is "reached via
  `alleles()` / `ploidy()` / `is_homozygous()`".** Now stale, since §2's own sketch is superseded by
  `arch/calling_priors.md` §3.2 on homozygosity. A design-doc edit, which this loop does not make.
- **The `#[inline]` asymmetry** — `alleles()` carries one, `new` does not. Deliberate: `new` sorts
  and allocates, and is called once per sample per locus.
- **The three red aggregate gates on `main`** (`clippy --all-targets`, `test --all-targets`,
  `doc --no-deps`) are unchanged by this step; recorded in
  [A1's review](ng_calling_a1_2026-08-21.md) §7.

## 8. Missing tests to add now

All three supplied as runnable code by `reliability`, all three applied:
`a_genotype_sorts_its_alleles_and_keeps_the_multiset` (the property test, in the file's existing
`proptest!` block), `a_genotype_cannot_be_built_from_no_alleles_at_all`, and
`a_genotype_new_is_idempotent_on_its_own_alleles` — the last catching a canonical form that is not
a fixed point, which the multiset test cannot see because a sort-then-rotate rule still makes two
spellings agree with each other.

## 9. What's good

- **The mutation pass reported four numbers, not two** — 8 run, 5 killed, 2 survived, 1 changed no
  behaviour — and the no-op (`sort` for `sort_unstable`) is argued rather than asserted: `AlleleId`
  wraps a plain `u16` with a derived total order, so equal entries are the same bit pattern and no
  fixture could separate the two. The agent explicitly says no test should be written against it.
- **The agent flagged the weaker half of its own evidence.** It could have presented both survivors
  as equal; instead it recorded that ids at or above 256 are unreachable through today's caps, so
  the three-copy mutation carries the finding.
- **It found a fix that costs nothing** — the two-entry swap in the tetraploid fixture — and
  verified it separately from the proptest, so the cheap guard and the thorough one are both on
  the table.
- **Two agents disagreed on M2's remedy and both arguments were made from this file's
  conventions**, which is what made the disagreement decidable rather than a matter of taste.
- **`naming` argued four rules down** rather than filing findings for the appearance of a
  violation, and cited the specific line of the specific sibling for each.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::types
./scripts/dev.sh cargo test --all-targets --all-features   # expect the pre-existing bench panic
```
