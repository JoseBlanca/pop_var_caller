"""The benchmarks, in the order a full run takes them."""

from . import (external_chr22, giab_chr20, giab_genome, giab_panel, giab_trio, hg002_5mb, hg002_tr,
               inbreeding_sim, refine_sim, smoke, speed, stutter_sim, tomato_agreement)

# Accuracy first, so a run cut short has it: the largest human set, the rest of GIAB, the
# simulations. Then speed, about 14 h on its own, and tomato_agreement, which scores the calls speed
# makes and so has to follow it.
BENCHMARKS = {
    "giab_chr20": giab_chr20,
    "giab_panel": giab_panel,
    "giab_trio": giab_trio,
    "hg002_5mb": hg002_5mb,
    "hg002_tr": hg002_tr,
    "stutter_sim": stutter_sim,
    "inbreeding_sim": inbreeding_sim,
    "refine_sim": refine_sim,
    "speed": speed,
    "tomato_agreement": tomato_agreement,
    "external_chr22": external_chr22,
    # Days long, so out of a full run (IN_FULL_RUN): run it with --only giab_genome.
    "giab_genome": giab_genome,
    "smoke": smoke,
}
