"""Score exp implementations against correctly rounded values.

Reads the directories `examples/float_exp_check.rs measure` writes. For every `stride`-th argument
of each range, taken in order of argument value so the sample spans the whole range, it computes exp at 200 bits with mpmath, rounds to the nearest double, and counts,
for each implementation, the arguments whose output is not that double and the largest distance
in units in the last place (steps between adjacent doubles). Given two directories it also counts
the arguments whose outputs differ between them, per implementation.

Usage:
  uv run --no-project --with mpmath python scripts/float_exp_accuracy.py <dir> [stride] [<other-dir>]
"""

import struct
import sys
from pathlib import Path

from mpmath import exp, mp, mpf

mp.prec = 200
SIDES = ("platform", "libm", "table")


def read(path: Path) -> list[int]:
    data = path.read_bytes()
    return list(struct.unpack(f"<{len(data) // 8}Q", data))


def as_double(bits: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", bits))[0]


def correctly_rounded_bits(x: float) -> int:
    # exp at 200 bits, as an exact integer mantissa and binary exponent; Python's int / int is
    # correctly rounded to the nearest double, subnormals included, which mpmath's own conversion
    # is not guaranteed to be below 2^-1022.
    mantissa, exponent = exp(mpf(x)).man_exp
    if exponent >= 0:
        rounded = float(int(mantissa) << exponent)
    else:
        rounded = int(mantissa) / (1 << -exponent)
    return struct.unpack("<Q", struct.pack("<d", rounded))[0]


def steps_apart(a: int, b: int) -> int:
    # Same sign for exp's outputs; subnormals and normals share one ordering in the bits.
    return abs(a - b)


def main() -> None:
    directory = Path(sys.argv[1])
    stride = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    other = Path(sys.argv[3]) if len(sys.argv) > 3 else None
    print("range\tscored\t" + "\t".join(f"{s} wrong\t{s} max steps" for s in SIDES))
    index = 0
    while (directory / f"exp.{index}.arguments.f64").exists():
        arguments = read(directory / f"exp.{index}.arguments.f64")
        outputs = {side: read(directory / f"exp.{index}.{side}.f64") for side in SIDES}
        wrong = {side: 0 for side in SIDES}
        worst = {side: 0 for side in SIDES}
        scored = 0
        order = sorted(range(len(arguments)), key=lambda i: as_double(arguments[i]))
        for position in order[::stride]:
            truth = correctly_rounded_bits(as_double(arguments[position]))
            scored += 1
            for side in SIDES:
                got = outputs[side][position]
                if got != truth:
                    wrong[side] += 1
                    worst[side] = max(worst[side], steps_apart(got, truth))
        print(
            f"{index}\t{scored}\t"
            + "\t".join(f"{wrong[s]}\t{worst[s]}" for s in SIDES)
        )
        if other is not None:
            differing = {
                side: sum(
                    1
                    for a, b in zip(outputs[side], read(other / f"exp.{index}.{side}.f64"))
                    if a != b
                )
                for side in SIDES
            }
            print(f"{index}\tagainst {other}\t" + "\t".join(f"{s} {differing[s]}" for s in SIDES))
        index += 1


main()
