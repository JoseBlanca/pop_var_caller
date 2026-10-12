"""Allele-length frequencies of each repeat class, estimated from the genotype rows without the
average over tract frequencies, and few-number shapes fitted the same way.

For a class, every (tract, sample) row carries the sample's genotype likelihoods with the stutter
held. Whatever an individual tract's frequencies, a sample's alleles at a random tract of the class
are drawn from the class's shares p, so the composite log-likelihood

    sum over rows of log sum_g prior(g | p, F) L(g),
    prior(a, a) = F p_a + (1 - F) p_a^2,  prior(a, b) = 2 (1 - F) p_a p_b,

is maximised by the class's shares. The free fit (one share per allele length) is by expectation
maximisation; its errors come from resampling whole tracts (rows of one tract are not independent).
Each shape is fitted by maximising the same composite likelihood over its few numbers.

**The two outermost allele lengths are catch-alls.** The window runs -span..+span (span 9 since
2026-10-12) while reads are recorded only to one repeat less, so the class at -span holds every
allele span or more repeats shorter, and +span every one span or more longer. A tail that falls
with distance cannot follow mass piled up there. Each shape therefore comes in two variants:

- the plain shape (as before 2026-10-12): its tail runs to the edge, which it treats as one more
  length;
- the `_edges` variant: the shape covers -(span-1)..+(span-1), and the two edge shares are two free
  numbers of their own, outside the tail.

Two measures compare a shape with the free fit, unchanged by the edges:

- **composite log-likelihood loss**: the free fit's composite log-likelihood minus the shape's,
  over the same rows. The edge allele is scored at exactly +-span, as the genotype rows give it.
- **one-sample calls changed**: rows whose most likely genotype (prior times the row's own
  likelihoods) differs between the free shares and the shape's. A genotype with an edge allele is
  one call, whatever distance beyond the edge the allele really sits at.

Usage:
    str_allele_length_shapes.py ROWS_DIR [min_tracts] [bootstraps]
        [--shapes all|plain|edges|NAME,NAME...] [--max-tracts N] [--jobs N] [--numbers FILE]

`--numbers FILE` writes each class's fitted numbers for every shape run, in the units a reader
can interpret (shares and powers), one row per class and shape.
"""

import argparse
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

SPAN = 6  # replaced by the span genotypes.tsv records, in set_span()
OFFSETS = np.arange(-SPAN, SPAN + 1)
K = len(OFFSETS)
FIRST, SECOND = None, None


def load(path):
    with open(path, "rb") as handle:
        rows, width = np.frombuffer(handle.read(16), dtype="<u8")
        rows, width = int(rows), int(width)
        tract = np.frombuffer(handle.read(4 * rows), dtype="<u4")
        excess = np.frombuffer(handle.read(8 * rows), dtype="<f8")
        values = np.frombuffer(handle.read(8 * rows * width), dtype="<f8").reshape(rows, width)
    return tract, excess, values


def set_span(directory):
    """Read the genotype pairs' allele classes, and set the span from them."""
    global SPAN, OFFSETS, K, FIRST, SECOND
    pairs = np.loadtxt(os.path.join(directory, "genotypes.tsv"), skiprows=1, dtype=int)
    SPAN = int(pairs.max())
    OFFSETS = np.arange(-SPAN, SPAN + 1)
    K = len(OFFSETS)
    FIRST, SECOND = pairs[:, 0] + SPAN, pairs[:, 1] + SPAN


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


# ---------------------------------------------------------------------
# The shapes. Each builder takes its numbers and the half-width it covers, and returns the shares
# over -span..+span.
# ---------------------------------------------------------------------


def with_reference(a, rest, span):
    rest = rest.copy()
    rest[span] = 0.0
    rest = rest / rest.sum()
    out = (1 - a) * rest
    out[span] = a
    return out


def side_tails(below_share, below, above, span):
    """Shorter alleles weighted by `below` (nearest last), longer by `above` (nearest first)."""
    w = np.zeros(2 * span + 1)
    w[:span] = below / below.sum() * below_share
    w[span + 1:] = above / above.sum() * (1 - below_share)
    return w


def distance(span):
    return np.arange(1, span + 1, dtype=float)


def geometric(r, span):
    return r ** (distance(span) - 1)


def power(s, span):
    return distance(span) ** (-s)


def step_power(second, s, span):
    """One-repeat weight 1; then `second` at two repeats, falling as a power of distance."""
    d = distance(span)
    return np.concatenate(([1.0], second * (d[1:] / 2.0) ** (-s)))


def smm(x, span):
    r = expit(x[0])
    w = r ** np.abs(np.arange(-span, span + 1)).astype(float)
    return w / w.sum()


PLAIN = {
    # the stepwise mutation model: one fall-off, symmetric
    "smm": (1, smm),
    "ref_geo4": (4, lambda x, n: with_reference(expit(x[0]), side_tails(expit(x[1]), geometric(expit(x[2]), n)[::-1], geometric(expit(x[3]), n), n), n)),
    "ref_power3": (3, lambda x, n: with_reference(expit(x[0]), side_tails(expit(x[1]), power(safe_exp(x[2]), n)[::-1], power(safe_exp(x[2]), n), n), n)),
    "ref_power4": (4, lambda x, n: with_reference(expit(x[0]), side_tails(expit(x[1]), power(safe_exp(x[2]), n)[::-1], power(safe_exp(x[3]), n), n), n)),
    # a power tail whose first step is free: the one-repeat share set apart on each side
    "ref_step_power6": (
        6,
        lambda x, n: with_reference(
            expit(x[0]),
            side_tails(expit(x[1]), step_power(expit(x[2]), safe_exp(x[3]), n)[::-1], step_power(expit(x[4]), safe_exp(x[5]), n), n),
            n,
        ),
    ),
}

# Which of each plain shape's numbers are powers (reported as exp(x)); the rest are shares
# (reported as expit(x)).
POWER_NUMBERS = {"smm": [], "ref_geo4": [], "ref_power3": [2], "ref_power4": [2, 3], "ref_step_power6": [3, 5]}
NUMBER_NAMES = {
    "smm": ["fall_off"],
    "ref_geo4": ["ref_share", "shorter_split", "fall_off_below", "fall_off_above"],
    "ref_power3": ["ref_share", "shorter_split", "power"],
    "ref_power4": ["ref_share", "shorter_split", "power_below", "power_above"],
    "ref_step_power6": ["ref_share", "shorter_split", "two_step_weight_below", "power_below", "two_step_weight_above", "power_above"],
}


def with_edges(base):
    """The base shape over -(span-1)..+(span-1), and the two edge shares as free numbers."""
    count, build = base

    def edged(x, span):
        inner = build(x[:count], span - 1)
        below, above = safe_exp(x[count]), safe_exp(x[count + 1])
        total = 1.0 + below + above
        return np.concatenate(([below / total], inner / total, [above / total]))

    return count + 2, edged


SHAPES = dict(PLAIN)
SHAPES.update({f"{name}_edges": with_edges(shape) for name, shape in PLAIN.items()})


def natural_numbers(shape, x, q):
    """The shape's numbers as shares and powers, plus figures read off its shares."""
    base = shape[: -len("_edges")] if shape.endswith("_edges") else shape
    out = {}
    for i, name in enumerate(NUMBER_NAMES[base]):
        out[name] = safe_exp(x[i]) if i in POWER_NUMBERS[base] else expit(x[i])
    out["edge_below"] = float(q[0])
    out["edge_above"] = float(q[-1])
    # Read off the shares, inside the edges, so the plain and edged shapes compare.
    inner = q[1:-1] / q[1:-1].sum()
    centre = SPAN - 1
    shorter, longer = inner[:centre].sum(), inner[centre + 1:].sum()
    out["inner_ref_share"] = float(inner[centre])
    out["inner_shorter_split"] = float(shorter / max(shorter + longer, 1e-300))
    out["one_shorter_of_shorter"] = float(inner[centre - 1] / max(shorter, 1e-300))
    out["one_longer_of_longer"] = float(inner[centre + 1] / max(longer, 1e-300))
    return out


def fit_shape(shape, excess, values, free):
    count, build = SHAPES[shape]

    def negative(x):
        q = build(x, SPAN)
        if not np.all(np.isfinite(q)):
            return 1e300
        return -composite(np.maximum(q, 1e-300), excess, values)

    rng = np.random.default_rng(7)
    best = None
    for attempt in range(8):
        start = rng.normal(0.0, 1.0, size=count)
        if attempt == 0:
            start[0] = math.log(free[SPAN] / max(1 - free[SPAN], 1e-9)) if count > 1 else 0.0
            if shape.endswith("_edges"):
                # The edge numbers start at the free fit's edge shares.
                rest = max(1.0 - free[0] - free[-1], 1e-9)
                start[-2] = math.log(max(free[0], 1e-9) / rest)
                start[-1] = math.log(max(free[-1], 1e-9) / rest)
        result = minimize(negative, start, method="Nelder-Mead", options={"maxiter": 3000 * count, "xatol": 1e-7, "fatol": 1e-9})
        if best is None or result.fun < best.fun:
            best = result
    return build(best.x, SPAN), -best.fun, best.x


def fit_class(job):
    """One class: the free fit, its bootstrap errors, and every shape against it."""
    directory, name, shapes, max_tracts, bootstraps, min_tracts = job
    period, repeats = (int(v) for v in name[:-4].split("_"))
    rng = np.random.default_rng([11, period, repeats])
    tract, excess, values = load(os.path.join(directory, name))
    tracts = int(tract.max()) + 1 if len(tract) else 0
    if tracts > max_tracts:
        kept = np.sort(rng.choice(tracts, size=max_tracts, replace=False))
        keep = np.isin(tract, kept)
        renumber = np.full(tracts, -1)
        renumber[kept] = np.arange(max_tracts)
        tract, excess, values = renumber[tract[keep]], excess[keep], values[keep]
        tracts = max_tracts
    if tracts < min_tracts or tracts == 0:
        return None
    started = time.time()
    free = em(excess, values)
    free_ll = composite(free, excess, values)
    # Bootstrap over tracts: each resample weights a row by how many times its tract was drawn.
    samples = []
    for _ in range(bootstraps):
        draw = np.bincount(rng.integers(0, tracts, size=tracts), minlength=tracts)
        samples.append(em(excess, values, weights=draw[tract].astype(float), start=free, iterations=300, tolerance=1e-7))
    errors = np.std(samples, axis=0, ddof=1) if bootstraps > 1 else np.full(K, np.nan)
    map_free = (priors(free, excess) * values).argmax(axis=1)
    # The free fit's own noise: the calls that change when it is swapped for one of its refits.
    noise = [int((map_free != (priors(np.maximum(b, 1e-300), excess) * values).argmax(axis=1)).sum()) for b in samples]
    out = {
        "class": name[:-4], "tracts": tracts, "rows": len(tract), "free": free, "errors": errors,
        "noise": float(np.mean(noise)) if noise else float("nan"), "shapes": {},
    }
    for shape in shapes:
        q, ll, x = fit_shape(shape, excess, values, free)
        map_shape = (priors(np.maximum(q, 1e-300), excess) * values).argmax(axis=1)
        out["shapes"][shape] = {
            "loss": free_ll - ll,
            "within": int((np.abs(q - free) / np.maximum(errors, 1e-6) <= 2).sum()),
            "miss": float(np.abs(q - free).max()),
            "changed": int((map_free != map_shape).sum()),
            "numbers": natural_numbers(shape, x, q),
        }
    print(f"{name[:-4]}: {tracts} tracts, {len(tract)} rows, {time.time() - started:.0f} s", file=sys.stderr, flush=True)
    return out


def chosen_shapes(spec):
    if spec == "none":
        return []
    if spec == "all":
        return list(SHAPES)
    if spec == "plain":
        return list(PLAIN)
    if spec == "edges":
        return [name for name in SHAPES if name.endswith("_edges")]
    names = spec.split(",")
    for name in names:
        if name not in SHAPES:
            raise SystemExit(f"unknown shape {name}; known: {', '.join(SHAPES)}")
    return names


def init_worker(directory):
    set_span(directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rows_dir")
    parser.add_argument("min_tracts", nargs="?", type=int, default=0)
    parser.add_argument("bootstraps", nargs="?", type=int, default=100)
    parser.add_argument("--shapes", default="all", help="all, plain (the shapes as before 2026-10-12), edges, none (the free fit and its noise only), or a comma list")
    parser.add_argument("--max-tracts", type=int, default=1500, help="tracts drawn per class at most (default 1500)")
    parser.add_argument("--jobs", type=int, default=1, help="classes fitted at once, one process each")
    parser.add_argument("--numbers", help="write each class's fitted numbers per shape to this TSV")
    args = parser.parse_args()

    set_span(args.rows_dir)
    shapes = chosen_shapes(args.shapes)
    names = sorted(
        (f for f in os.listdir(args.rows_dir) if f.endswith(".bin")),
        key=lambda f: tuple(int(x) for x in f[:-4].split("_")),
    )
    jobs = [(args.rows_dir, name, shapes, args.max_tracts, args.bootstraps, args.min_tracts) for name in names]
    if args.jobs > 1:
        import multiprocessing

        with multiprocessing.get_context("spawn").Pool(args.jobs, initializer=init_worker, initargs=(args.rows_dir,)) as pool:
            results = pool.map(fit_class, jobs, chunksize=1)
    else:
        results = [fit_class(job) for job in jobs]
    results = [r for r in results if r is not None]

    print(f"class\ttracts\trows\tshares -{SPAN}..+{SPAN} (free fit; bootstrap error in brackets for those above 0.005; the two ends are 'this far or further')")
    for r in results:
        shares = " ".join(f"{s:.3f}" + (f"({e:.3f})" if s > 0.005 else "") for s, e in zip(r["free"], r["errors"]))
        print(f"{r['class']}\t{r['tracts']}\t{r['rows']}\t{shares}")

    print()
    print("shape\tnumbers\tshares within 2 bootstrap errors\tmedian log-lik loss\tworst log-lik loss\tone-sample calls changed, per 10,000 rows\tworst class, per 10,000")
    for shape in shapes:
        fits = [r["shapes"][shape] for r in results]
        rows = sum(r["rows"] for r in results)
        changed = sum(f["changed"] for f in fits)
        worst = max(f["changed"] / r["rows"] for f, r in zip(fits, results))
        within = sum(f["within"] for f in fits)
        losses = [f["loss"] for f in fits]
        print(f"{shape}\t{SHAPES[shape][0]}\t{within} of {K * len(results)}\t{np.median(losses):.1f}\t{np.max(losses):.1f}\t{1e4 * changed / rows:.1f}\t{1e4 * worst:.0f}")
    rows = sum(r["rows"] for r in results)
    noise = sum(r["noise"] for r in results)
    worst = max(r["noise"] / r["rows"] for r in results)
    print(f"free fit vs its own refits (noise)\t-\t-\t-\t-\t{1e4 * noise / rows:.1f}\t{1e4 * worst:.0f}")
    print()
    print("per class: for each shape, composite log-lik loss against the free fit / largest share miss / rows whose one-sample call changes")
    print("class\ttracts\trows\tnoise (mean calls changed by a refit)\t" + "\t".join(shapes))
    for r in results:
        print(f"{r['class']}\t{r['tracts']}\t{r['rows']}\t{r['noise']:.1f}\t" + "\t".join(f"{r['shapes'][s]['loss']:.1f}/{r['shapes'][s]['miss']:.3f}/{r['shapes'][s]['changed']}" for s in shapes))

    if args.numbers:
        keys = sorted({key for r in results for s in shapes for key in r["shapes"][s]["numbers"]})
        with open(args.numbers, "w") as handle:
            handle.write("class\ttracts\trows\tshape\tloss\tchanged\t" + "\t".join(keys) + "\tfree_edge_below\tfree_edge_above\n")
            for r in results:
                for s in shapes:
                    f = r["shapes"][s]
                    cells = [f"{f['numbers'][k]:.4f}" if k in f["numbers"] else "NA" for k in keys]
                    handle.write(f"{r['class']}\t{r['tracts']}\t{r['rows']}\t{s}\t{f['loss']:.1f}\t{f['changed']}\t" + "\t".join(cells) + f"\t{r['free'][0]:.4f}\t{r['free'][-1]:.4f}\n")


if __name__ == "__main__":
    main()
