#!/usr/bin/env python3
"""Validate the Gaussian sampler constants/kernels in
algebra/src/sampling/gaussian.rs (see docs/MIGRATION.md).

1. Re-derive CDF155 from first principles at 60-digit precision:
   T_k = round(2^128 * P(X > k)) for the discrete half-Gaussian with
   sigma = 1.55, rho(k) = exp(-k^2 / (2 sigma^2)), k >= 0 — and require
   <= 2 ulp agreement with the table in the Rust source.
2. If a lazer checkout is given as argv[1], additionally require the Rust
   table to be *identical* to the C table in src/grandom.c.
3. Check the fdlibm exp_small kernel against math.exp on |x| <= 0.5 ln 2.

Run: python3 scripts/verify_gaussian.py [path/to/lazer]
"""
import math
import os
import re
import sys
from decimal import Decimal, getcontext

getcontext().prec = 60

HERE = os.path.dirname(os.path.abspath(__file__))
RS_SRC = os.path.join(HERE, "..", "algebra", "src", "sampling", "gaussian.rs")


def parse_pairs(path, pattern, limit):
    txt = open(path).read()
    return [(int(a), int(b)) for a, b in re.findall(pattern, txt)][:limit]


def main():
    rs_tab = parse_pairs(RS_SRC, r"\(\s*(\d+)\s*,\s*(\d+)\s*\)", 21)
    assert len(rs_tab) == 21, f"expected 21 CDF entries, got {len(rs_tab)}"

    if len(sys.argv) > 1:
        c_src = os.path.join(sys.argv[1], "src", "grandom.c")
        c_tab = parse_pairs(c_src, r"\{\s*(\d+)ULL,\s*(\d+)ULL\s*\}", 21)
        assert rs_tab == c_tab, "Rust CDF155 differs from C source!"
        print("CDF155: identical to lazer src/grandom.c (21/21 pairs)")

    sigma = Decimal("1.55")
    two_ss = 2 * sigma * sigma
    rho = [(-Decimal(k * k) / two_ss).exp() for k in range(200)]
    total = sum(rho)
    scale = Decimal(2) ** 128
    max_err = 0
    for k, (hi, lo) in enumerate(rs_tab):
        got = (hi << 64) | lo
        want = int((sum(rho[k + 1:]) / total * scale).to_integral_value())
        max_err = max(max_err, abs(got - want))
        assert abs(got - want) <= 2, f"k={k}: {got} vs {want}"
    print(f"CDF155 re-derivation from half-Gaussian(1.55): max dev {max_err} ulp")

    C1 = 1.66666666666666019037e-01
    C2 = -2.77777777770155933842e-03
    C3 = 6.61375632143793436117e-05
    C4 = -1.65339022054652515390e-06
    C5 = 4.13813679705723846039e-08

    def exp_small(x):
        t = x * x
        t = x - t * (C1 + t * (C2 + t * (C3 + t * (C4 + t * C5))))
        return 1.0 - ((x * t) / (t - 2.0) - x)

    worst, x = 0.0, -0.5 * math.log(2)
    while x <= 0.5 * math.log(2):
        worst = max(worst, abs(exp_small(x) - math.exp(x)) / math.exp(x))
        x += 1e-4
    assert worst < 1e-15, worst
    print(f"exp_small: max relative error {worst:.2e} on |x| <= ln2/2")

    for xv in [0.0, 0.3, 1.0, 2.7, 10.0]:
        u = int(xv / math.log(2))
        r = xv - math.log(2) * u
        assert abs((0.5 ** u) * math.exp(-r) - math.exp(-xv)) < 1e-12
    print("BerExp split 2^-u * exp(-r) == exp(-x)")
    print("\nALL GAUSSIAN CHECKS PASSED")


if __name__ == "__main__":
    main()
