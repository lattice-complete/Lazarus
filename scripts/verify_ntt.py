#!/usr/bin/env python3
"""Reference-model cross-validation for algebra::ntt (see docs/MIGRATION.md).

This is a 1:1 Python port of the Rust NTT engine — same deterministic prime
search, same runtime root derivation, same Longa–Naehrig bit-reversed
negacyclic butterflies, same Garner CRT — validated against an exact
integer schoolbook negacyclic product on full-range inputs.

Run: python3 scripts/verify_ntt.py
"""
import random

Q = (1 << 32) - 99


# ---- algebra/src/ntt/prime.rs ----
def is_prime_u64(n):
    if n < 2:
        return False
    for sp in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n == sp:
            return True
        if n % sp == 0:
            return False
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    for w in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        x = pow(w, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def find_root_2adic(p, a):
    odd = (p - 1) >> a
    x = 2
    while True:
        y = pow(x, odd, p)
        if pow(y, 1 << (a - 1), p) == p - 1:
            return y
        x += 1


def search_below_2_62():
    m = (1 << 42) - 1
    while True:
        p = (m << 20) + 1
        if is_prime_u64(p):
            a = ((p - 1) & -(p - 1)).bit_length() - 1
            return p, a
        m -= 1


# ---- algebra/src/ntt/mod.rs ----
def bitrev(x, bits):
    return int(format(x, f"0{bits}b")[::-1], 2) if bits else 0


class Table:
    def __init__(self, d, p, two_adicity, root):
        log2d = d.bit_length() - 1
        k = (2 * d).bit_length() - 1
        assert k <= two_adicity
        psi = pow(root, 1 << (two_adicity - k), p)
        assert pow(psi, 2 * d, p) == 1 and pow(psi, d, p) == p - 1
        ipsi = pow(psi, p - 2, p)
        pows = [pow(psi, i, p) for i in range(d)]
        ipows = [pow(ipsi, i, p) for i in range(d)]
        self.psi_rev = [pows[bitrev(i, log2d)] for i in range(d)]
        self.ipsi_rev = [ipows[bitrev(i, log2d)] for i in range(d)]
        self.d_inv = pow(d, p - 2, p)
        self.d, self.p = d, p

    def forward(self, a):
        d, p = self.d, self.p
        t, m = d, 1
        while m < d:
            t >>= 1
            for i in range(m):
                j1 = 2 * i * t
                s = self.psi_rev[m + i]
                for j in range(j1, j1 + t):
                    u, v = a[j], a[j + t] * s % p
                    a[j] = (u + v) % p
                    a[j + t] = (u - v) % p
            m <<= 1

    def inverse(self, a):
        d, p = self.d, self.p
        t, m = 1, d
        while m > 1:
            h = m >> 1
            j1 = 0
            for i in range(h):
                s = self.ipsi_rev[h + i]
                for j in range(j1, j1 + t):
                    u, v = a[j], a[j + t]
                    a[j] = (u + v) % p
                    a[j + t] = (u - v) * s % p
                j1 += 2 * t
            t <<= 1
            m = h
        for j in range(d):
            a[j] = a[j] * self.d_inv % p


def negacyclic_mul_centered(a, b, tables, p0, p1):
    res = []
    for (p, tb) in tables:
        fa = [x % p for x in a]
        fb = [x % p for x in b]
        tb.forward(fa)
        tb.forward(fb)
        fc = [x * y % p for x, y in zip(fa, fb)]
        tb.inverse(fc)
        res.append(fc)
    p0_inv_p1 = pow(p0 % p1, p1 - 2, p1)
    P = p0 * p1
    out = []
    for r0, r1 in zip(res[0], res[1]):
        diff = (r1 + p1 - (r0 % p1)) % p1
        t = diff * p0_inv_p1 % p1
        x = r0 + p0 * t
        out.append(x - P if x > P // 2 else x)
    return out


def schoolbook(a, b):
    d = len(a)
    out = [0] * d
    for i in range(d):
        for j in range(d):
            pr = a[i] * b[j]
            if i + j < d:
                out[i + j] += pr
            else:
                out[i + j - d] -= pr
    return out


def main():
    p0, a0 = 0xFFFFFFFF00000001, 32
    assert is_prime_u64(p0) and ((p0 - 1) & ((1 << a0) - 1)) == 0
    p1, a1 = search_below_2_62()
    print(f"p0 = {p0} (Goldilocks, 2-adicity {a0})")
    print(f"p1 = {p1} = m*2^20+1 (2-adicity {a1})  [deterministic search]")
    assert a1 >= 13, "must support degree 4096"
    r0 = find_root_2adic(p0, a0)
    r1 = find_root_2adic(p1, a1)
    print(f"roots: {r0}, {r1} (orders verified)")

    bound = 4096 * (Q // 2 + 1) ** 2
    assert bound < p0 * p1 // 2, "CRT bound violated"
    print(f"bound ok: d*(q/2)^2 = 2^{bound.bit_length()} < P/2 = 2^{(p0 * p1 // 2).bit_length()}")

    random.seed(1)
    for d in [2, 64, 128, 512]:
        tabs = [(p0, Table(d, p0, a0, r0)), (p1, Table(d, p1, a1, r1))]
        for trial in range(8 if d <= 128 else 2):
            a = [random.randint(-(1 << 31), 1 << 31) for _ in range(d)]
            b = [random.randint(-(1 << 31), 1 << 31) for _ in range(d)]
            got = negacyclic_mul_centered(a, b, tabs, p0, p1)
            assert got == schoolbook(a, b), f"MISMATCH d={d} trial={trial}"
        print(f"d={d}: NTT == schoolbook on full-range random inputs OK")

    tabs2 = [(p0, Table(2, p0, a0, r0)), (p1, Table(2, p1, a1, r1))]
    assert negacyclic_mul_centered([1, 2], [3, 4], tabs2, p0, p1) == [-5, 10]
    print("(1+2X)(3+4X) mod X^2+1 = -5+10X OK")
    print("\nALL NTT/CRT CHECKS PASSED")


if __name__ == "__main__":
    main()
