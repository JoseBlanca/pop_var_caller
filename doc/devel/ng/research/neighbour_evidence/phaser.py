"""A prototype phaser for the neighbour-evidence pass: rounds of the pass's own copying walk, with
each sample's two haplotypes as its templates' contribution, scored on founder_grid.py cells
against the simulator's truth.

usage: uv run --with numpy python phaser.py [--rounds R] [--haplotypes 16 [--one-per-sample]]
                                           [--hide N] [--dump DIR] <cell dir> [...]

--dump DIR writes each round's haplotypes, with the truth, to DIR/<cell>.round<r>.pickle, for
phase_accuracy.py.

Each sample carries two haplotypes, soft: at each record, each copy's probability of carrying the
alternative. Round 1 starts every template sample's two haplotypes at its GP allele share, so the
walk is the pass as now, with each template counted twice. In every round, for each sample:

  1. its 16 template samples — chosen once, as pvc chooses them — give 32 haplotype templates, plus
     the background template (the cohort's frequencies);
  2. the walk gives its genotype posterior at every record (what the pass would write);
  3. one path through the states is drawn from the walk's forward probabilities, backwards — which
     template each copy follows at each record, or the one both copies follow where autozygous;
  4. along that path, each copy's allele at each record is re-estimated from the followed
     templates' haplotypes and the sample's own genotype likelihood: P(copy one a, copy two b) ∝
     q̃_k1(a) · q̃_k2(b) · L(a, b), so a copy stays uncertain where its template and the reads are.

With --haplotypes 16 (spec neighbour_evidence_phasing.md §10 Q2) the templates are the 16 best
haplotypes of all other samples' instead of both haplotypes of the 16 best samples, so the walk has
as many templates as the pass has now. Round 1 is then the pass as now exactly (16 samples, their GP
shares, each once). From round 2 a haplotype's score is the pessimistic rate at which it carries an
allele the sample's genotype lacks — x·P(0/0) + (1 − x)·P(1/1) per record, x its alternative
probability, the genotype the previous round's posterior — over the sample's informative records
(at most 500, spread evenly) where the haplotype's sample is informative too, the 16 lowest kept;
with --one-per-sample, at most one haplotype of each sample.

With --hide N (the fill check without the echo, owner 2026-11-02): before round 1, N records where
each sample's reads are informative are drawn at random, and at them the sample is treated as
having no reads in every round — its walk and its path read a flat likelihood there. Its averaged
neighbours-only prediction there therefore holds none of its own reads at that record, through any
other sample's haplotype either, which is the situation of a real ./. . The fill check compares that
prediction, where its GQ reaches 10, with the genotype the sample's own evidence settles (posterior
≥ 0.999 from its reads, its inbreeding and the cohort frequency, as the pass's settled_genotype). The
genotype scored there is the prediction times the sample's likelihood, renormalised.

All samples are updated together at the end of a round, from one seed, so the result does not
depend on the order or the threads. A drawn path, not the posterior's marginals, because the two
copies are symmetric in the walk: averaged, each copy's marginal is the same and the phase is lost.

Scored as phasing_ceiling.py: biallelic SNP records matching a truth record, the sample's genotype
from the round's posterior; the walk's lengths fixed at the simulator's, and the pass's rules for
leaving a sample ./. not applied. The target's own reads phase its templates here, as they would
in the pass — the comparison with Beagle, run without the target, is slightly in the phaser's favour.
"""
import math, multiprocessing, os, sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phasing_ceiling as pc  # noqa: E402
from hmm import background_shares, haplotypes, stays  # noqa: E402

MISMATCH = 0.01
SEED = 20261102


def lik_matrix(lik):
    l0, l1, l2 = lik
    return np.array([[l0, l1], [l1, l2]])


def draw_path(rows, target, q_sites, inbreeding, rng, liks):
    """One state path, drawn backwards from the forward probabilities; then each copy's
    probability of carrying the alternative along it."""
    f = min(max(inbreeding, 1e-4), 1 - 1e-4)
    n = len(rows)
    T = q_sites[0].shape[0]
    pos = [r['pos'] for r in rows]
    fwd_d, fwd_p, Ms = [], [], []
    prev_d = prev_p = None
    for i, r in enumerate(rows):
        q = q_sites[i]
        M = lik_matrix(liks[i])
        Ms.append(M)
        e_pair = q @ M @ q.T
        e_desc = q @ np.diag(M)
        if prev_d is None:
            d = np.full(T, f / T)
            p = np.full((T, T), (1 - f) / (T * T))
        else:
            gap = pos[i] - pos[i - 1]
            s_d, s_n = stays(gap, f, pc.RUN_BASES)
            keep = math.exp(-gap / pc.SEGMENT_BASES)
            mass_d, mass_p = prev_d.sum(), prev_p.sum()
            d_in = keep * prev_d + (1 - keep) * mass_d / T
            rows_, cols = prev_p.sum(axis=1), prev_p.sum(axis=0)
            p_in = (keep * keep * prev_p + keep * (1 - keep) * (rows_[:, None] + cols[None, :]) / T
                    + (1 - keep) ** 2 * mass_p / (T * T))
            d = s_d * d_in + (1 - s_n) * mass_p / T
            p = s_n * p_in + (1 - s_d) * mass_d / (T * T)
        d = d * e_desc
        p = p * e_pair
        s = d.sum() + p.sum()
        d, p = d / s, p / s
        fwd_d.append(d); fwd_p.append(p)
        prev_d, prev_p = d, p

    def pick(w_d, w_p):
        flat = np.concatenate([w_d, w_p.ravel()])
        flat = flat / flat.sum()
        index = rng.choice(len(flat), p=flat)
        return ('d', index) if index < T else ('p', divmod(index - T, T))

    states = [None] * n
    states[-1] = pick(fwd_d[-1], fwd_p[-1])
    for i in range(n - 2, -1, -1):
        gap = pos[i + 1] - pos[i]
        s_d, s_n = stays(gap, f, pc.RUN_BASES)
        keep = math.exp(-gap / pc.SEGMENT_BASES)
        kind, k = states[i + 1]
        if kind == 'p':
            k1, k2 = k
            w_p = np.full((T, T), (1 - keep) ** 2 / (T * T))
            w_p[k1, :] += keep * (1 - keep) / T
            w_p[:, k2] += keep * (1 - keep) / T
            w_p[k1, k2] += keep * keep
            w_p = fwd_p[i] * s_n * w_p
            w_d = fwd_d[i] * (1 - s_d) / (T * T)
        else:
            w_d = np.full(T, (1 - keep) / T)
            w_d[k] += keep
            w_d = fwd_d[i] * s_d * w_d
            w_p = fwd_p[i] * (1 - s_n) / T
        states[i] = pick(w_d, w_p)

    first = np.empty(n)
    second = np.empty(n)
    for i, (kind, k) in enumerate(states):
        q, M = q_sites[i], Ms[i]
        if kind == 'p':
            joint = np.outer(q[k[0]], q[k[1]]) * M
        else:
            joint = np.diag(q[k] * np.diag(M))
        joint = joint / joint.sum()
        first[i] = joint[1, :].sum()
        second[i] = joint[:, 1].sum()
    return first, second


def template_q(rows, pairs, haps, contig):
    """Per record, the blurred allele shares of the haplotype templates `pairs` — (sample, copy) —
    and the background, last."""
    out = []
    for i, r in enumerate(rows):
        alt = [haps[t][contig][c][i] for t, c in pairs]
        q = np.array([[1 - a, a] for a in alt] + [list(background_shares([r['af']]))])
        out.append(q * (1 - MISMATCH) + MISMATCH / 2)
    return out


def choose_haplotypes(cell, target, previous):
    """The 16 haplotypes (sample, copy) whose alleles the sample's genotype lacks least."""
    samples, _, _ = pc.DATA[cell]
    n = len(samples)
    picked = [(contig, i) for contig, rows in pc.CONTIGS[cell].items()
              for i, r in enumerate(rows) if pc.informative(r['lik'][target])]
    if len(picked) > pc.COMPARED:
        picked = [picked[j * len(picked) // pc.COMPARED] for j in range(pc.COMPARED)]
    d = np.zeros((n, 2))
    shared = np.zeros(n)
    for contig, i in picked:
        r = pc.CONTIGS[cell][contig][i]
        p00, _, p11 = previous[contig][i]
        for o in range(n):
            if o == target or not pc.informative(r['lik'][o]):
                continue
            shared[o] += 1
            for c in (0, 1):
                x = HAPS[cell][o][contig][c][i]
                d[o, c] += x * p00 + (1 - x) * p11
    scored = []
    for o in range(n):
        if o == target:
            continue
        for c in (0, 1):
            mean = (d[o, c] + 1) / (shared[o] + 2)
            scored.append((mean + 2 * math.sqrt(mean * (1 - mean) / (shared[o] + 3)), o, c))
    ranked = sorted(scored)
    if ONE_PER_SAMPLE:
        seen, kept = set(), []
        for _, o, c in ranked:
            if o not in seen:
                seen.add(o)
                kept.append((o, c))
        return kept[:pc.K]
    return [(o, c) for _, o, c in ranked[:pc.K]]


def one_sample(args):
    cell, target, round_index = args
    samples, records, inbreeding = pc.DATA[cell]
    templates = pc.choose(records, target, len(samples))
    if HAPLOTYPE_TEMPLATES is None:
        pairs = [(t, c) for t in templates for c in (0, 1)]
    elif round_index == 1:
        pairs = [(t, 0) for t in templates]
    else:
        pairs = choose_haplotypes(cell, target, PREVIOUS[target])
    rng = np.random.default_rng([SEED, round_index, target])
    new, posts = {}, {}
    for contig, rows in pc.CONTIGS[cell].items():
        q_sites = template_q(rows, pairs, HAPS[cell], contig)
        # hmm.haplotypes blurs the shares itself, so it is given them unblurred
        raw = [(q - MISMATCH / 2) / (1 - MISMATCH) for q in q_sites]
        hidden = HIDDEN.get(cell, {}).get(target, set())
        liks = [np.ones(3) if (contig, i) in hidden else r['lik'][target] for i, r in enumerate(rows)]
        posteriors, _ = haplotypes([(lik, None, None) for lik in liks],
                                   [r['pos'] for r in rows], raw, inbreeding[target],
                                   pc.RUN_BASES, pc.SEGMENT_BASES)
        posts[contig] = np.array(posteriors)
        new[contig] = draw_path(rows, target, q_sites, inbreeding[target], rng, liks)
    return target, posts, new


def score(cell, posts):
    """Counts over every sample, from each sample's genotype posteriors `posts[target][contig]`."""
    counts = dict(reads=0, calling_wrong=0, pass_wrong=0, missing=0, filled=0, filled_wrong=0,
                  checked=0, disagreed=0)
    samples, _, inbreeding = pc.DATA[cell]
    for target, by_contig in posts.items():
        hidden = HIDDEN.get(cell, {}).get(target, set())
        for contig, rows in pc.CONTIGS[cell].items():
            for i, (r, post) in enumerate(zip(rows, by_contig[contig])):
                if (contig, i) in hidden:
                    # the fill check: the clean prediction against the settled genotype
                    settled = settled_genotype(r['lik'][target], r['af'], inbreeding[target])
                    best_prediction = int(np.argmax(post))
                    if (settled is not None
                            and -10 * math.log10(max(1 - post[best_prediction], 1e-10)) >= 10):
                        counts['checked'] += 1
                        counts['disagreed'] += best_prediction != settled
                    post = post * r['lik'][target]
                    post = post / post.sum()
                truth = tuple(sorted(r['truth'][target]))
                best = int(np.argmax(post))
                genotype = [(0, 0), (0, 1), (1, 1)][best]
                if r['called'][target] is not None:
                    counts['reads'] += 1
                    counts['calling_wrong'] += r['called'][target] != truth
                    counts['pass_wrong'] += genotype != truth
                else:
                    counts['missing'] += 1
                    if -10 * math.log10(max(1 - post[best], 1e-10)) >= 10:
                        counts['filled'] += 1
                        counts['filled_wrong'] += genotype != truth
    return counts


def settled_genotype(lik, af, inbreeding):
    """The genotype index the sample's own evidence settles (posterior ≥ 0.999), or None."""
    p, f = af, min(max(inbreeding, 0.0), 1.0)
    prior = np.array([f * (1 - p) + (1 - f) * (1 - p) ** 2, 2 * (1 - f) * p * (1 - p),
                      f * p + (1 - f) * p * p])
    joint = prior * lik
    total = joint.sum()
    best = int(np.argmax(joint))
    return best if total > 0 and joint[best] >= 0.999 * total else None


def row(name, label, c):
    share = f" ({100 * c['filled_wrong'] / c['filled']:.1f}%)" if c['filled'] else ''
    check = (f" | {c['disagreed']} of {c['checked']} ({100 * c['disagreed'] / c['checked']:.1f}%)"
             if c.get('checked') else ' | ')
    return (f"| {name} | {label} | {c['calling_wrong']:,} → {c['pass_wrong']:,} of {c['reads']:,} | "
            f"{c['missing']:,} | {c['filled']:,} | {c['filled_wrong']:,}{share}{check} |")


HAPS = {}
HIDDEN = {}
PREVIOUS = {}
HAPLOTYPE_TEMPLATES = None
ONE_PER_SAMPLE = False


def main():
    global HAPLOTYPE_TEMPLATES, PREVIOUS, ONE_PER_SAMPLE
    args = sys.argv[1:]
    rounds = 4
    hide = 0
    dump = None
    while args[0].startswith('--'):
        if args[0] == '--rounds':
            rounds = int(args[1])
        elif args[0] == '--dump':
            dump = args[1]
        elif args[0] == '--hide':
            hide = int(args[1])
        elif args[0] == '--haplotypes':
            HAPLOTYPE_TEMPLATES = int(args[1])
        elif args[0] == '--one-per-sample':
            ONE_PER_SAMPLE = True
            args = args[1:]
            continue
        args = args[2:]
    print('| cell | round | wrong where samples have reads, calling → pass | missing | filled (GQ ≥ 10) | of them wrong | fill check: wrong of checked |')
    print('|---|---|---|---|---|---|---|')
    for cell in args:
        pc.DATA[cell] = pc.load(cell)
        samples, records, _ = pc.DATA[cell]
        pc.CONTIGS[cell] = {}
        for r in records:
            pc.CONTIGS[cell].setdefault(r['contig'], []).append(r)
        pc.RANKS[id(records)] = [pc.ranking(records, t, len(samples)) for t in range(len(samples))]
        HAPS[cell] = {t: {contig: (np.array([r['shares'][t][1] for r in rows]),) * 2
                          for contig, rows in pc.CONTIGS[cell].items()}
                      for t in range(len(samples))}
        name = os.path.basename(os.path.normpath(cell))
        HIDDEN[cell] = {}
        if hide:
            for t in range(len(samples)):
                candidates = [(contig, i) for contig, rows in pc.CONTIGS[cell].items()
                              for i, r in enumerate(rows) if pc.informative(r['lik'][t])]
                rng = np.random.default_rng([SEED, 999, t])
                chosen = rng.choice(len(candidates), size=min(hide, len(candidates)), replace=False)
                HIDDEN[cell][t] = {candidates[j] for j in chosen}
        summed = None
        for round_index in range(1, rounds + 1):
            posts, new_haps = {}, {}
            with ProcessPoolExecutor(6, mp_context=multiprocessing.get_context('fork')) as pool:
                for target, post, new in pool.map(one_sample, [(cell, t, round_index)
                                                               for t in range(len(samples))]):
                    new_haps[target] = new
                    posts[target] = post
            HAPS[cell] = new_haps
            if dump:
                import pickle
                os.makedirs(dump, exist_ok=True)
                name_ = os.path.basename(os.path.normpath(cell))
                with open(os.path.join(dump, f'{name_}.round{round_index}.pickle'), 'wb') as fh:
                    pickle.dump({'haplotypes': new_haps,
                                 'records': [[(r['pos'], r['truth']) for r in rows]
                                             for rows in pc.CONTIGS[cell].values()],
                                 'contigs': list(pc.CONTIGS[cell].keys())}, fh)
            PREVIOUS = posts
            print(row(name, str(round_index), score(cell, posts)), flush=True)
            # the average over rounds from the second on: the first is the pass as now
            if round_index >= 2:
                if summed is None:
                    summed = {t: {c: p.copy() for c, p in by.items()} for t, by in posts.items()}
                else:
                    for t, by in posts.items():
                        for c, p in by.items():
                            summed[t][c] += p
                averaged = {t: {c: p / (round_index - 1) for c, p in by.items()}
                            for t, by in summed.items()}
                print(row(name, f'mean of 2–{round_index}', score(cell, averaged)), flush=True)


if __name__ == '__main__':
    main()
