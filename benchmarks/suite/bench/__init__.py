"""The benchmarks, in the order a full run takes them."""

from . import giab_panel, giab_trio, hg002_5mb, hg002_tr, smoke

BENCHMARKS = {
    "giab_panel": giab_panel,
    "giab_trio": giab_trio,
    "hg002_5mb": hg002_5mb,
    "hg002_tr": hg002_tr,
    "smoke": smoke,
}
