"""The benchmarks, in the order a full run takes them."""

from . import giab_panel, smoke

BENCHMARKS = {
    "giab_panel": giab_panel,
    "smoke": smoke,
}
