"""The benchmarks, in the order a full run takes them."""

from . import external_chr22, giab_chr20, giab_panel, giab_trio, hg002_5mb, hg002_tr, smoke, speed, stutter_sim, tomato_agreement

# The most quoted figures first, so a run cut short has them: speed, then the largest human set.
BENCHMARKS = {
    "speed": speed,
    "tomato_agreement": tomato_agreement,
    "giab_chr20": giab_chr20,
    "giab_panel": giab_panel,
    "giab_trio": giab_trio,
    "hg002_5mb": hg002_5mb,
    "hg002_tr": hg002_tr,
    "stutter_sim": stutter_sim,
    "external_chr22": external_chr22,
    "smoke": smoke,
}
