"""Write the Rust model's reference walks: random windows run through hmm.haplotypes.

usage: uv run --with numpy python make_reference_walks.py <output file>

Format, whitespace-separated, one case after another:
  case <templates> <sites> <segment_bases> <run_bases> <autozygous_share> <mismatch> <background>
  site <position> <alleles> <likelihoods...> <template shares, template by template...> [<alt frequencies...>]
  expect <log_likelihood>
  posterior <genotype posteriors...>      (one line a site)

<templates> counts the sampled templates only. With <background> 1 each site line ends with the
record's alternative-allele frequencies, and the walk ran with one more template, the background,
whose shares hmm.background_shares makes from them (spec §3.3).
"""
import sys
import numpy as np
import hmm

rng = np.random.default_rng(7)
out = open(sys.argv[1], 'w')
cases = [
    # templates, sites, allele counts drawn from, segment, run, share, mismatch, gap scale,
    # background
    (3, 8, [2], 3e4, 1e5, 0.8, 0.01, 2e3, False),
    (4, 12, [2, 3, 4], 3e5, 2e6, 0.5, 0.01, 5e4, False),
    (2, 6, [2, 5], 1e4, 3e4, 0.2, 0.05, 1e3, False),
    # The background template, at sites where the sample's reads disagree with both sampled
    # templates, so its rows carry real weight.
    (2, 10, [2, 3, 4], 3e4, 3e5, 0.9, 0.01, 5e3, True),
]
for templates, n, allele_counts, segment, run, share, mismatch, gap, background in cases:
    positions = np.cumsum(1 + rng.exponential(gap, n)).astype(np.int64)
    sites, tmpls, lines = [], [], []
    for position in positions:
        a = int(rng.choice(allele_counts))
        g = a * (a + 1) // 2
        lik = np.exp(rng.normal(0, 3, g))
        lik /= lik.max()
        shares = rng.dirichlet(np.full(a, 0.3), templates)
        line = 'site %d %d %s %s' % (position, a, ' '.join('%.17g' % x for x in lik),
                                     ' '.join('%.17g' % x for x in shares.ravel()))
        if background:
            alt_frequencies = rng.dirichlet(np.ones(a))[1:]
            shares = np.vstack([shares, hmm.background_shares(alt_frequencies)])
            line += ' ' + ' '.join('%.17g' % x for x in alt_frequencies)
        sites.append((lik, None, None))
        tmpls.append(shares)
        lines.append(line)
    posts, ll = hmm.haplotypes(sites, positions.astype(float), tmpls, share, run, segment, mismatch)
    out.write('case %d %d %.17g %.17g %.17g %.17g %d\n'
              % (templates, n, segment, run, share, mismatch, int(background)))
    out.write('\n'.join(lines) + '\n')
    out.write('expect %.17g\n' % ll)
    for p in posts:
        out.write('posterior %s\n' % ' '.join('%.17g' % x for x in p))
