"""Prototype neighbour-evidence models over pvc's per-sample genotype likelihoods.

Diploid only. A record has A alleles; genotypes are in VCF order (0/0, 0/1, 1/1, 0/2, ...).
Every per-site input is a likelihood row normalised to max 1 (scale does not matter).

Three ways to call a sample's genotype at each record:
  baseline      pvc's prior: (1-F) * independent draws + F * identical by descent, F genome-wide
  autozygosity  the same two branches, mixed by a two-state HMM's posterior at that site
  haplotypes    each copy copies one of K templates (other samples' allele-probability vectors);
                the autozygous state makes both copies copy one template
"""

import numpy as np


def genotype_pairs(n_alleles):
    """(j, k) with j <= k for each genotype index in VCF order."""
    return [(j, k) for k in range(n_alleles) for j in range(k + 1)]


def prior_rows(alpha):
    """The two branches of pvc's prior for one sample at one record, each a probability row.

    Independent draws: Dirichlet-multinomial with concentration alpha, ploidy 2.
    Identical by descent: homozygote a with probability alpha_a / sum(alpha).
    """
    alpha = np.asarray(alpha, dtype=float)
    total = alpha.sum()
    pairs = genotype_pairs(len(alpha))
    independent = np.empty(len(pairs))
    descent = np.zeros(len(pairs))
    for g, (j, k) in enumerate(pairs):
        if j == k:
            independent[g] = alpha[j] * (alpha[j] + 1) / (total * (total + 1))
            descent[g] = alpha[j] / total
        else:
            independent[g] = 2 * alpha[j] * alpha[k] / (total * (total + 1))
    return independent / independent.sum(), descent / descent.sum()


def posterior(lik, prior):
    p = lik * prior
    return p / p.sum()


def baseline(sites, inbreeding):
    """sites: list of (lik, independent, descent). Returns genotype posteriors."""
    return [posterior(l, (1 - inbreeding) * ind + inbreeding * ibd) for l, ind, ibd in sites]


def _two_state_forward_backward(emit, stay):
    """emit: (n, 2) emission per site; stay: (n, 2) probability each state persists from the
    previous site (stay[0] unused); stationary start. Returns (n, 2) state posteriors and log
    likelihood."""
    n = emit.shape[0]
    fwd = np.empty((n, 2))
    scale = np.empty(n)
    start = stay_stationary(stay[1] if n > 1 else np.array([0.5, 0.5]))
    a = start * emit[0]
    scale[0] = a.sum()
    fwd[0] = a / scale[0]
    for i in range(1, n):
        s0, s1 = stay[i]
        prev = fwd[i - 1]
        a = np.array([prev[0] * s0 + prev[1] * (1 - s1), prev[0] * (1 - s0) + prev[1] * s1]) * emit[i]
        scale[i] = a.sum()
        fwd[i] = a / scale[i]
    bwd = np.ones((n, 2))
    for i in range(n - 2, -1, -1):
        s0, s1 = stay[i + 1]
        b = bwd[i + 1] * emit[i + 1]
        bwd[i] = np.array([s0 * b[0] + (1 - s0) * b[1], (1 - s1) * b[0] + s1 * b[1]])
        bwd[i] /= bwd[i].sum()
    post = fwd * bwd
    post /= post.sum(axis=1, keepdims=True)
    return post, np.log(scale).sum()


def stay_stationary(stay_pair):
    s0, s1 = stay_pair
    leave0, leave1 = 1 - s0, 1 - s1
    if leave0 + leave1 <= 0:
        return np.array([0.5, 0.5])
    return np.array([leave1, leave0]) / (leave0 + leave1)


def autozygosity(sites, positions, inbreeding, run_bases):
    """Two-state HMM per sample. State 1 = autozygous. run_bases = mean autozygous run length;
    the non-autozygous run length follows so the stationary share is `inbreeding`.

    Returns (genotype posteriors, per-site P(autozygous), log likelihood)."""
    f = min(max(inbreeding, 1e-4), 1 - 1e-4)
    run_ibd = run_bases
    run_not = run_bases * (1 - f) / f
    pos = np.asarray(positions, dtype=float)
    gaps = np.diff(pos, prepend=pos[0])
    gaps[gaps < 0] = 1e9  # a new contig: forget
    stay = np.stack([np.exp(-gaps / run_not), np.exp(-gaps / run_ibd)], axis=1)
    emit = np.array([[l @ ind, l @ ibd] for l, ind, ibd in sites])
    emit /= emit.max(axis=1, keepdims=True)
    state, loglik = _two_state_forward_backward(emit, stay)
    out = []
    for (l, ind, ibd), (p_not, p_ibd) in zip(sites, state):
        out.append(p_not * posterior(l, ind) + p_ibd * posterior(l, ibd))
    return out, state[:, 1], loglik


def fit_autozygosity(sites, positions, inbreeding, run_grid=(1e4, 3e4, 1e5, 3e5, 1e6, 3e6)):
    """Pick the run length by likelihood, re-estimate F as the mean autozygous posterior
    (two rounds). Returns the result of `autozygosity` and the chosen (F, run)."""
    f = inbreeding
    best = None
    for _ in range(2):
        for run in run_grid:
            res = autozygosity(sites, positions, f, run)
            if best is None or res[2] > best[0][2]:
                best = (res, (f, run))
        f = float(np.clip(best[0][1].mean(), 0.01, 0.99))
    return best


def haplotypes(sites, positions, templates, inbreeding, run_bases, segment_bases, mismatch=0.01):
    """Copying model for one target sample.

    templates: per site an array (K, A) — each template's probability of carrying each allele
    (another sample's expected allele share). States: autozygous with template k (K states), or
    not autozygous with templates (k1, k2) (K*K states). Each copy switches template with
    probability 1 - exp(-gap/segment_bases) between sites; the autozygous state persists as in
    `autozygosity`. Emission: the target's genotype likelihood under the templates' alleles,
    blurred by `mismatch`. Returns (genotype posteriors, log likelihood)."""
    f = min(max(inbreeding, 1e-4), 1 - 1e-4)
    run_ibd, run_not = run_bases, run_bases * (1 - f) / f
    pos = np.asarray(positions, dtype=float)
    n = len(sites)
    K = templates[0].shape[0]
    fwd_d, fwd_p, scales = [], [], []
    emits = []

    def emissions(i):
        lik, ind, ibd = sites[i]
        q = templates[i]
        A = q.shape[1]
        q = q * (1 - mismatch) + mismatch / A
        pairs = genotype_pairs(A)
        M = np.zeros((A, A))
        for g, (j, k) in enumerate(pairs):
            M[j, k] = M[k, j] = lik[g]
        e_pair = q @ M @ q.T  # (K, K): copy one from k1, copy two from k2
        e_desc = np.einsum("ka,a->k", q, np.diag(M))  # both copies the same template allele
        return e_desc, e_pair, M, q

    prev_d = prev_p = None
    for i in range(n):
        e_desc, e_pair, M, q = emissions(i)
        emits.append((M, q, e_desc, e_pair))
        if i == 0 or pos[i] < pos[i - 1]:
            d = np.full(K, f / K)
            p = np.full((K, K), (1 - f) / (K * K))
        else:
            gap = pos[i] - pos[i - 1]
            s_d, s_n = np.exp(-gap / run_ibd), np.exp(-gap / run_not)
            keep = np.exp(-gap / segment_bases)
            mass_d, mass_p = prev_d.sum(), prev_p.sum()
            # within the autozygous class: one template, switches like one copy
            d_in = keep * prev_d + (1 - keep) * mass_d / K
            # within the other class: each copy switches on its own
            rows, cols = prev_p.sum(axis=1), prev_p.sum(axis=0)
            p_in = (keep * keep * prev_p + keep * (1 - keep) * (rows[:, None] + cols[None, :]) / K
                    + (1 - keep) ** 2 * mass_p / (K * K))
            # class changes: redraw templates uniformly
            d = s_d * d_in + (1 - s_n) * mass_p / K
            p = s_n * p_in + (1 - s_d) * mass_d / (K * K)
        d = d * e_desc
        p = p * e_pair
        s = d.sum() + p.sum()
        d, p = d / s, p / s
        fwd_d.append(d); fwd_p.append(p); scales.append(s)
        prev_d, prev_p = d, p

    # backward pass, same transitions transposed (they are symmetric within a class)
    bwd_d = [None] * n
    bwd_p = [None] * n
    bwd_d[-1] = np.ones(K); bwd_p[-1] = np.ones((K, K))
    for i in range(n - 2, -1, -1):
        _, _, e_desc, e_pair = emits[i + 1]
        bd = bwd_d[i + 1] * e_desc
        bp = bwd_p[i + 1] * e_pair
        if pos[i + 1] < pos[i]:
            bwd_d[i] = np.full(K, f / K * bd.sum() + (1 - f) / (K * K) * bp.sum())
            bwd_p[i] = np.full((K, K), bwd_d[i][0])
        else:
            gap = pos[i + 1] - pos[i]
            s_d, s_n = np.exp(-gap / run_ibd), np.exp(-gap / run_not)
            keep = np.exp(-gap / segment_bases)
            d_in = keep * bd + (1 - keep) * bd.sum() / K
            rows, cols = bp.sum(axis=1), bp.sum(axis=0)
            p_in = (keep * keep * bp + keep * (1 - keep) * (rows[:, None] + cols[None, :]) / K
                    + (1 - keep) ** 2 * bp.sum() / (K * K))
            bwd_d[i] = s_d * d_in + (1 - s_d) * bp.sum() / (K * K)
            bwd_p[i] = s_n * p_in + (1 - s_n) * bd.sum() / K
        t = bwd_d[i].sum() + bwd_p[i].sum()
        bwd_d[i] /= t; bwd_p[i] /= t

    out = []
    for i in range(n):
        M, q, e_desc, e_pair = emits[i]
        A = q.shape[1]
        # state posterior with this site's own emission divided back out, so the target's
        # likelihood enters once, below
        wd = fwd_d[i] * bwd_d[i] / e_desc
        wp = fwd_p[i] * bwd_p[i] / e_pair
        # allele-pair weights implied by the template posterior, times the target's likelihood
        pair_w = q.T @ wp @ q  # (A, A) ordered
        desc_w = wd @ q  # (A,)
        lik = sites[i][0]
        pairs = genotype_pairs(A)
        g = np.empty(len(pairs))
        for gi, (j, k) in enumerate(pairs):
            w = pair_w[j, k] + (pair_w[k, j] if j != k else 0.0)
            if j == k:
                w += desc_w[j]
            g[gi] = w * lik[gi]
        out.append(g / g.sum())
    return out, float(np.log(scales).sum())
