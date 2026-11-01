"""Write the Rust model's reference walks: random windows run through hmm.haplotypes.

usage: uv run --with numpy python make_reference_walks.py <output file>

Format, whitespace-separated, one case after another:
  case <templates> <sites> <segment_bases> <run_bases> <autozygous_share> <mismatch>
  site <position> <alleles> <likelihoods...> <template shares, template by template...>
  expect <log_likelihood>
  posterior <genotype posteriors...>      (one line a site)
"""
import sys
import numpy as np
import hmm

rng = np.random.default_rng(7)
out = open(sys.argv[1], 'w')
cases = [
    # templates, sites, allele counts drawn from, segment, run, share, mismatch, gap scale
    (3, 8, [2], 3e4, 1e5, 0.8, 0.01, 2e3),
    (4, 12, [2, 3, 4], 3e5, 2e6, 0.5, 0.01, 5e4),
    (2, 6, [2, 5], 1e4, 3e4, 0.2, 0.05, 1e3),
]
for templates, n, allele_counts, segment, run, share, mismatch, gap in cases:
    positions = np.cumsum(1 + rng.exponential(gap, n)).astype(np.int64)
    sites, tmpls, lines = [], [], []
    for position in positions:
        a = int(rng.choice(allele_counts))
        g = a * (a + 1) // 2
        lik = np.exp(rng.normal(0, 3, g))
        lik /= lik.max()
        shares = rng.dirichlet(np.full(a, 0.3), templates)
        sites.append((lik, None, None))
        tmpls.append(shares)
        lines.append('site %d %d %s %s' % (position, a, ' '.join('%.17g' % x for x in lik),
                                            ' '.join('%.17g' % x for x in shares.ravel())))
    posts, ll = hmm.haplotypes(sites, positions.astype(float), tmpls, share, run, segment, mismatch)
    out.write('case %d %d %.17g %.17g %.17g %.17g\n' % (templates, n, segment, run, share, mismatch))
    out.write('\n'.join(lines) + '\n')
    out.write('expect %.17g\n' % ll)
    for p in posts:
        out.write('posterior %s\n' % ' '.join('%.17g' % x for x in p))
