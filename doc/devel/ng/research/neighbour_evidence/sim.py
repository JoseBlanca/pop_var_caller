"""Genotype-likelihood-level cohort simulator with known truth.

Each chromosome copy is a mosaic of founder haplotypes (segments of mean `segment` bases);
inside autozygous runs (mean `run` bases, genome share F) the second copy equals the first.
Reads: Poisson(depth) per site per sample, each from a random copy, 1% error. The likelihood row
is computed from the read counts as pvc's per-read model would for a SNP at one base quality.

usage: uv run --with numpy sim.py <population> <n_samples> <depth> [seed]
"""

import sys
import time

import numpy as np

import hmm

POPULATIONS = {
    # tomato-like: few haplotypes, long shared segments, mostly autozygous
    "inbred": dict(founders=8, segment=300_000, inbreeding=0.8, run=2_000_000),
    # outbred, many haplotypes, short segments, almost no autozygosity (an outbred human)
    "outbred": dict(founders=400, segment=30_000, inbreeding=0.01, run=300_000),
    # outbred but a single sample from a consanguineous family (F = 1/16)
    "consanguineous": dict(founders=400, segment=30_000, inbreeding=0.0625, run=3_000_000),
}
LENGTH = 4_000_000
SITE_SPACING = 1_000
ERROR = 0.01
SEED_TOTAL, SEED_ALT = 10.0, 0.05  # single-sample concentration, alt share 0.5%


def mosaic(rng, length_sites, positions, founders, segment):
    out = np.empty(length_sites, dtype=np.int64)
    i = 0
    while i < length_sites:
        f = rng.integers(founders)
        end_pos = positions[i] + rng.exponential(segment)
        j = np.searchsorted(positions, end_pos)
        out[i:max(j, i + 1)] = f
        i = max(j, i + 1)
    return out


def autozygous_track(rng, positions, f, run):
    run_not = run * (1 - f) / f
    track = np.empty(len(positions), dtype=bool)
    state = rng.random() < f
    next_change = positions[0] + rng.exponential(run if state else run_not)
    for i, p in enumerate(positions):
        while p > next_change:
            state = not state
            next_change += rng.exponential(run if state else run_not)
        track[i] = state
    return track


def simulate(population, n_samples, depth, seed):
    rng = np.random.default_rng(seed)
    pop = POPULATIONS[population]
    positions = np.sort(rng.choice(LENGTH, LENGTH // SITE_SPACING, replace=False))
    n_sites = len(positions)
    # founder alleles: derived allele frequency among founders from a 1/x spectrum
    freq = np.exp(rng.uniform(np.log(1 / pop["founders"]), np.log(0.5), n_sites))
    founder_alleles = (rng.random((pop["founders"], n_sites)) < freq).astype(np.int64)
    truth = np.empty((n_samples, n_sites), dtype=np.int64)  # alt copies 0, 1, 2
    autoz = np.empty((n_samples, n_sites), dtype=bool)
    counts = np.empty((n_samples, n_sites, 2), dtype=np.int64)
    for s in range(n_samples):
        c1 = mosaic(rng, n_sites, positions, pop["founders"], pop["segment"])
        c2 = mosaic(rng, n_sites, positions, pop["founders"], pop["segment"])
        track = autozygous_track(rng, positions, pop["inbreeding"], pop["run"])
        c2 = np.where(track, c1, c2)
        a1 = founder_alleles[c1, np.arange(n_sites)]
        a2 = founder_alleles[c2, np.arange(n_sites)]
        truth[s] = a1 + a2
        autoz[s] = track
        reads = rng.poisson(depth, n_sites)
        alt_share = truth[s] / 2 * (1 - ERROR) + (1 - truth[s] / 2) * ERROR
        alt = rng.binomial(reads, alt_share)
        counts[s, :, 0], counts[s, :, 1] = reads - alt, alt
    return positions, truth, autoz, counts


def likelihood_rows(counts):
    ref, alt = counts[..., 0], counts[..., 1]
    log_l = np.stack([
        ref * np.log(1 - ERROR) + alt * np.log(ERROR),
        (ref + alt) * np.log(0.5),
        ref * np.log(ERROR) + alt * np.log(1 - ERROR),
    ], axis=-1)
    log_l -= log_l.max(axis=-1, keepdims=True)
    return np.exp(log_l)


def biallelic_priors(alpha):
    """Vectorised hmm.prior_rows for two alleles: alpha (..., 2) -> (independent, descent)."""
    r, a = alpha[..., 0], alpha[..., 1]
    t = r + a
    ind = np.stack([r * (r + 1), 2 * r * a, a * (a + 1)], axis=-1) / (t * (t + 1))[..., None]
    ibd = np.stack([r / t, np.zeros_like(t), a / t], axis=-1)
    return ind, ibd


def cohort_alphas(lik, inbreeding):
    """Leave-one-out concentrations from a few rounds of the cohort's expected allele copies,
    the way pvc's frequency loop builds them."""
    n, m, _ = lik.shape
    seed = np.array([SEED_TOTAL - SEED_ALT, SEED_ALT])
    alt_copies = np.zeros((n, m))
    for _ in range(4):
        total_alt = alt_copies.sum(axis=0)
        loo_alt = total_alt[None, :] - alt_copies
        alpha = np.stack([seed[0] + 2 * (n - 1) - loo_alt, seed[1] + loo_alt], axis=-1)
        ind, ibd = biallelic_priors(alpha)
        f = inbreeding[:, None, None]
        post = lik * ((1 - f) * ind + f * ibd)
        post /= post.sum(axis=-1, keepdims=True)
        alt_copies = post[..., 1] + 2 * post[..., 2]
    total_alt = alt_copies.sum(axis=0)
    alphas = np.empty((n, m, 2))
    for s in range(n):
        loo_alt = total_alt - alt_copies[s]
        alphas[s, :, 0] = seed[0] + 2 * (n - 1) - loo_alt
        alphas[s, :, 1] = seed[1] + loo_alt
    return alphas


def genome_wide_inbreeding(lik, alphas):
    """pvc fits one F per sample from its heterozygote deficit; here: the F that maximises the
    sample's likelihood over a grid, which is what that fit estimates."""
    grid = np.linspace(0, 0.99, 34)
    ind, ibd = biallelic_priors(alphas)
    li, lb = (lik * ind).sum(-1), (lik * ibd).sum(-1)
    scores = np.stack([np.log((1 - f) * li + f * lb).sum(axis=1) for f in grid], axis=1)
    return grid[np.argmax(scores, axis=1)]


def score(posteriors, truth):
    called = np.array([np.argmax(p) for p in posteriors])
    nonref = truth > 0
    return dict(
        errors=int((called != truth).sum()),
        nonref_wrong=int((called[nonref] != truth[nonref]).sum()),
        nonref=int(nonref.sum()),
        het_wrong=int((called[truth == 1] != 1).sum()),
        het=int((truth == 1).sum()),
    )


def run(population, n_samples, depth, seed=1, scored=4, templates_k=16):
    positions, truth, autoz, counts = simulate(population, n_samples, depth, seed)
    lik = likelihood_rows(counts)
    if n_samples == 1:
        alphas = np.broadcast_to(np.array([SEED_TOTAL - SEED_ALT, SEED_ALT]), lik.shape[:2] + (2,))
    else:
        alphas = cohort_alphas(lik, np.full(n_samples, POPULATIONS[population]["inbreeding"]))
    f_hat = genome_wide_inbreeding(lik, alphas)
    ind, ibd = biallelic_priors(np.asarray(alphas, dtype=float))
    rows = [[(lik[s, i], ind[s, i], ibd[s, i]) for i in range(lik.shape[1])]
            for s in range(n_samples)]
    base = [hmm.baseline(rows[s], f_hat[s]) for s in range(n_samples)]
    results = {"baseline": [], "autozygosity": [], "haplotypes": []}
    timing = {"autozygosity": 0.0, "haplotypes": 0.0}
    dosage = np.array([[p[1] / 2 + p[2] for p in b] for b in base])  # each sample's alt share
    for s in range(min(scored, n_samples)):
        results["baseline"].append(score(base[s], truth[s]))
        t = time.time()
        (az, _, _), _ = hmm.fit_autozygosity(rows[s], positions, f_hat[s])
        timing["autozygosity"] += time.time() - t
        results["autozygosity"].append(score(az, truth[s]))
        if n_samples > 1:
            t = time.time()
            others = [j for j in range(n_samples) if j != s]
            # the K others whose calls agree best with this sample's over the region
            agree = [-np.abs(dosage[j] - dosage[s]).sum() for j in others]
            chosen = [others[j] for j in np.argsort(agree)[::-1][:templates_k]]
            tmpl = [np.stack([[1 - dosage[j, i], dosage[j, i]] for j in chosen]) for i in range(len(positions))]
            best = None
            for segment in (3e4, 3e5):
                post, ll = hmm.haplotypes(rows[s], positions, tmpl, f_hat[s],
                                          POPULATIONS[population]["run"], segment)
                if best is None or ll > best[1]:
                    best = (post, ll)
            timing["haplotypes"] += time.time() - t
            results["haplotypes"].append(score(best[0], truth[s]))
    summary = {}
    for name, rs in results.items():
        if not rs:
            continue
        summary[name] = {k: sum(r[k] for r in rs) for k in rs[0]}
    return summary, timing, f_hat[:scored], autoz[:scored].mean(axis=1)


if __name__ == "__main__":
    population, n, depth = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    summary, timing, f_hat, f_true = run(population, n, depth, seed)
    for name, r in summary.items():
        print(f"{population}\t{n}\t{depth}\t{name}\t{r['errors']}\t{r['nonref_wrong']}\t{r['nonref']}"
              f"\t{r['het_wrong']}\t{r['het']}\t{timing.get(name, 0):.1f}")
    print(f"# F fitted {np.round(f_hat, 2)} true share autozygous {np.round(f_true, 2)}", file=sys.stderr)
