#!/usr/bin/env bash
# Sweep the simulator: population x cohort size x depth. One process, one thread.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
for pop in inbred consanguineous outbred; do
  for depth in 1 3 10 30; do
    for seed in 1 2 3 4 5 6 7 8; do uv run --with numpy python sim.py $pop 1 $depth $seed; done
    for n in 10 63; do uv run --with numpy python sim.py $pop $n $depth 1; done
  done
done
