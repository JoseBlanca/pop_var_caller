"""Allele-length frequencies of each repeat class, estimated from the genotype rows without the
average over tract frequencies, and few-number shapes fitted the same way.

For a class, every (tract, sample) row carries the sample's genotype likelihoods with the stutter
held. Whatever an individual tract's frequencies, a sample's alleles at a random tract of the class
are drawn from the class's shares p, so the composite log-likelihood

    sum over rows of log sum_g prior(g | p, F) L(g),
    prior(a, a) = F p_a + (1 - F) p_a^2,  prior(a, b) = 2 (1 - F) p_a p_b,

is maximised by the class's shares. The free fit (13 shares) is by expectation maximisation; its
errors come from resampling whole tracts (rows of one tract are not independent). Each shape is
fitted by maximising the same composite likelihood over its few numbers.

Usage: rows_shapes.py ROWS_DIR [min_tracts] [bootstraps]
"""

import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

SPAN = 6  # replaced by the span genotypes.tsv records, in main()
OFFSETS = np.arange(-SPAN, SPAN + 1)
K = len(OFFSETS)
MAX_TRACTS = 1500


def load(path):
    with open(path, "rb") as handle:
        rows, width = np.frombuffer(handle.read(16), dtype="<u8")
        rows, width = int(rows), int(width)
        tract = np.frombuffer(handle.read(4 * rows), dtype="<u4")
        excess = np.frombuffer(handle.read(8 * rows), dtype="<f8")
        values = np.frombuffer(handle.read(8 * rows * width), dtype="<f8").reshape(rows, width)
    return tract, excess, values


def genotype_order(directory):
    """The genotype pairs' allele classes, and setting the span from them."""
    global SPAN, OFFSETS, K, DISTANCE
    pairs = np.loadtxt(os.path.join(directory, "genotypes.tsv"), skiprows=1, dtype=int)
    SPAN = int(pairs.max())
    OFFSETS = np.arange(-SPAN, SPAN + 1)
    K = len(OFFSETS)
    DISTANCE = np.arange(1, SPAN + 1, dtype=float)
    return pairs[:, 0] + SPAN, pairs[:, 1] + SPAN


FIRST, SECOND = None, None


def priors(p, excess):
    """Per row, per genotype, prior(g | p, F)."""
    hom = FIRST == SECOND
    independent = np.where(hom, p[FIRST] ** 2, 2 * p[FIRST] * p[SECOND])
    by_descent = np.where(hom, p[FIRST], 0.0)
    return excess[:, None] * by_descent[None, :] + (1 - excess)[:, None] * independent[None, :]


def composite(p, excess, values, weights=None):
    totals = (priors(p, excess) * values).sum(axis=1)
    logs = np.log(np.maximum(totals, 1e-300))
    return float(logs.sum() if weights is None else (weights * logs).sum())


def em(excess, values, weights=None, start=None, iterations=2000, tolerance=1e-10):
    """Free shares by expectation maximisation, the number of allele draws as the latent count."""
    p = np.full(K, 1.0 / K) if start is None else start.copy()
    w = np.ones(len(excess)) if weights is None else weights
    hom = FIRST == SECOND
    for _ in range(iterations):
        joint_ibd = excess[:, None] * np.where(hom, p[FIRST], 0.0)[None, :] * values
        joint_ind = (1 - excess)[:, None] * np.where(hom, p[FIRST] ** 2, 2 * p[FIRST] * p[SECOND])[None, :] * values
        total = (joint_ibd + joint_ind).sum(axis=1)
        total = np.maximum(total, 1e-300)
        ibd = (joint_ibd / total[:, None]) * w[:, None]
        ind = (joint_ind / total[:, None]) * w[:, None]
        draws = np.zeros(K)
        np.add.at(draws, FIRST, ibd.sum(axis=0))
        np.add.at(draws, FIRST, ind.sum(axis=0))
        np.add.at(draws, SECOND, ind.sum(axis=0))
        new = draws / draws.sum()
        if np.abs(new - p).max() < tolerance:
            p = new
            break
        p = new
    return p


def safe_exp(x):
    return math.exp(max(min(x, 50.0), -50.0))


def expit(x):
    x = max(min(x, 700.0), -700.0)
    return 1.0 / (1.0 + math.exp(-x))


def with_reference(a, rest):
    rest = rest.copy()
    rest[SPAN] = 0.0
    rest = rest / rest.sum()
    out = (1 - a) * rest
    out[SPAN] = a
    return out


def side_tails(below_share, below, above):
    w = np.zeros(K)
    w[OFFSETS < 0] = below / below.sum() * below_share
    w[OFFSETS > 0] = above / above.sum() * (1 - below_share)
    return w


DISTANCE = np.arange(1, SPAN + 1, dtype=float)


def geometric(r):
    return r ** (DISTANCE - 1)


def power(s):
    return DISTANCE ** (-s)


def flip(v):
    return v[::-1]


SHAPES = {
    # the stepwise mutation model: one fall-off, symmetric
    "smm": (1, lambda x: (lambda r: (lambda w: w / w.sum())(r ** np.abs(OFFSETS).astype(float)))(expit(x[0]))),
    "ref_geo4": (4, lambda x: with_reference(expit(x[0]), side_tails(expit(x[1]), flip(geometric(expit(x[2]))), geometric(expit(x[3]))))),
    "ref_power3": (3, lambda x: with_reference(expit(x[0]), side_tails(expit(x[1]), flip(power(safe_exp(x[2]))), power(safe_exp(x[2]))))),
    "ref_power4": (4, lambda x: with_reference(expit(x[0]), side_tails(expit(x[1]), flip(power(safe_exp(x[2]))), power(safe_exp(x[3]))))),
    # a power tail whose first step is free: the one-repeat share set apart on each side
    "ref_step_power6": (
        6,
        lambda x: with_reference(
            expit(x[0]),
            side_tails(
                expit(x[1]),
                flip(np.concatenate(([1.0], expit(x[2]) * power(safe_exp(x[3]))[1:] / power(safe_exp(x[3]))[1]))),
                np.concatenate(([1.0], expit(x[4]) * power(safe_exp(x[5]))[1:] / power(safe_exp(x[5]))[1])),
            ),
        ),
    ),
}


def fit_shape(shape, excess, values, free):
    count, build = SHAPES[shape]

    def negative(x):
        q = build(x)
        if not np.all(np.isfinite(q)):
            return 1e300
        return -composite(np.maximum(q, 1e-300), excess, values)

    rng = np.random.default_rng(7)
    best = None
    for attempt in range(8):
        start = rng.normal(0.0, 1.0, size=count)
        if attempt == 0:
            start[0] = math.log(free[SPAN] / max(1 - free[SPAN], 1e-9)) if count > 1 else 0.0
        result = minimize(negative, start, method="Nelder-Mead", options={"maxiter": 3000 * count, "xatol": 1e-7, "fatol": 1e-9})
        if best is None or result.fun < best.fun:
            best = result
    return build(best.x), -best.fun


def main():
    global FIRST, SECOND
    directory = sys.argv[1]
    min_tracts = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    bootstraps = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    FIRST, SECOND = genotype_order(directory)
    names = sorted(
        (f for f in os.listdir(directory) if f.endswith(".bin")),
        key=lambda f: tuple(int(x) for x in f[:-4].split("_")),
    )
    rng = np.random.default_rng(11)
    totals = {shape: {"within": 0, "shares": 0, "loss": [], "changed": 0, "rows": 0, "changed_by_class": []} for shape in SHAPES}
    print("class\ttracts\trows\tshares -6..+6 (free fit; bootstrap error in brackets for those above 0.005)")
    details = []
    for name in names:
        tract, excess, values = load(os.path.join(directory, name))
        tracts = int(tract.max()) + 1 if len(tract) else 0
        if tracts > MAX_TRACTS:
            kept = np.sort(rng.choice(tracts, size=MAX_TRACTS, replace=False))
            keep = np.isin(tract, kept)
            renumber = np.full(tracts, -1)
            renumber[kept] = np.arange(MAX_TRACTS)
            tract, excess, values = renumber[tract[keep]], excess[keep], values[keep]
            tracts = MAX_TRACTS
        if tracts < min_tracts:
            continue
        free = em(excess, values)
        free_ll = composite(free, excess, values)
        # Bootstrap over tracts: each resample weights a row by how many times its tract was drawn.
        samples = []
        for _ in range(bootstraps):
            draw = np.bincount(rng.integers(0, tracts, size=tracts), minlength=tracts)
            samples.append(em(excess, values, weights=draw[tract].astype(float), start=free, iterations=300, tolerance=1e-7))
        errors = np.std(samples, axis=0, ddof=1)
        shares = " ".join(f"{s:.3f}" + (f"({e:.3f})" if s > 0.005 else "") for s, e in zip(free, errors))
        print(f"{name[:-4]}\t{tracts}\t{len(tract)}\t{shares}")
        row = {"class": name[:-4], "tracts": tracts}
        for shape in SHAPES:
            q, ll = fit_shape(shape, excess, values, free)
            misses = np.abs(q - free) / np.maximum(errors, 1e-6)
            totals[shape]["within"] += int((misses <= 2).sum())
            totals[shape]["shares"] += K
            totals[shape]["loss"].append(free_ll - ll)
            # The one-sample call: each row's most likely genotype under its own reads and each prior.
            map_free = (priors(free, excess) * values).argmax(axis=1)
            map_shape = (priors(np.maximum(q, 1e-300), excess) * values).argmax(axis=1)
            changed = int((map_free != map_shape).sum())
            totals[shape]["changed"] += changed
            totals[shape]["rows"] += len(map_free)
            totals[shape]["changed_by_class"].append(changed / len(map_free))
            row[shape] = (free_ll - ll, float(np.abs(q - free).max()), changed)
        details.append(row)

    print()
    print("shape\tnumbers\tshares within 2 bootstrap errors\tmedian log-lik loss\tworst log-lik loss\tone-sample calls changed, per 10,000 rows\tworst class, per 10,000")
    for shape, (count, _) in SHAPES.items():
        t = totals[shape]
        print(f"{shape}\t{count}\t{t['within']} of {t['shares']}\t{np.median(t['loss']):.1f}\t{np.max(t['loss']):.1f}\t{1e4 * t['changed'] / t['rows']:.1f}\t{1e4 * max(t['changed_by_class']):.0f}")
    print()
    print("per class: for each shape, composite log-lik loss against the free fit / largest share miss / rows whose one-sample call changes")
    print("class\ttracts\t" + "\t".join(SHAPES))
    for row in details:
        print(f"{row['class']}\t{row['tracts']}\t" + "\t".join(f"{row[s][0]:.1f}/{row[s][1]:.3f}/{row[s][2]}" for s in SHAPES))


if __name__ == "__main__":
    main()
