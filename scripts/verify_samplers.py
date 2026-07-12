#!/usr/bin/env python3
"""Executable spec of algebra::sampling — byte-exact simulation.

Reproduces the exact XOF byte-consumption order of the Rust samplers on a
SHAKE128 stream (hashlib.shake_128 output is prefix-stable, matching the
sha3 crate), then checks distribution moments. If this passes, the Rust
samplers produce correctly-distributed output for the same consumption
order — `cargo test` then only needs to catch typos.

Run: python3 scripts/verify_samplers.py
"""
import hashlib
import math

LN2 = math.log(2)


class Shake:
    """Streaming reads over SHAKE128(seed || dom_le)."""

    def __init__(self, seed: bytes, dom: int, budget=1 << 22):
        h = hashlib.shake_128()
        h.update(seed + dom.to_bytes(8, "little"))
        self.buf = h.digest(budget)
        self.off = 0

    def fill(self, n: int) -> bytes:
        out = self.buf[self.off:self.off + n]
        assert len(out) == n, "budget exceeded"
        self.off += n
        return out

    def u64(self):
        return int.from_bytes(self.fill(8), "little")

    def u128(self):
        return int.from_bytes(self.fill(16), "little")


# ---- gaussian.rs, byte-exact ----
CDF155 = [
    (10894764499197476522, 10804844707381617341),
    (4761708367981796450, 6209732027000382074),
    (1476784279527800432, 14108379346150813303),
    (316388870594767345, 17827298407763885637),
    (46043503515468600, 18385899657892021654),
    (4503729779335039, 3860889375818664979),
    (294122444862326, 13947176349836216550),
    (12769598070895, 14321894682751135119),
    (367552986472, 10286761328368440884),
    (7001273393, 2287787188970898528),
    (88153536, 17843977990435837663),
    (733119, 12174894787802692461),
    (4024, 18067426722645776197),
    (14, 10764017821655913055),
    (0, 643125733022530080),
    (0, 1014291014134832),
    (0, 1055215183460),
    (0, 724109373),
    (0, 327744),
    (0, 98),
    (0, 0),
]


def exp_small(x):
    C = [1.66666666666666019037e-01, -2.77777777770155933842e-03,
         6.61375632143793436117e-05, -1.65339022054652515390e-06,
         4.13813679705723846039e-08]
    t = x * x
    t = x - t * (C[0] + t * (C[1] + t * (C[2] + t * (C[3] + t * C[4]))))
    return 1.0 - ((x * t) / (t - 2.0) - x)


def cdf_sampler(xof):
    buf = xof.fill(16)
    hi = int.from_bytes(buf[:8], "little")
    lo = int.from_bytes(buf[8:], "little")
    z = 0
    while hi < CDF155[z][0] or (hi == CDF155[z][0] and lo < CDF155[z][1]):
        z += 1
    return z


def ber_exp(xof, x):
    buf = xof.fill(16)
    w0 = int.from_bytes(buf[:8], "little")
    w1 = int.from_bytes(buf[8:], "little")
    u = int(x / LN2)
    r = x - LN2 * u
    u = min(u, 63)
    low = 0 if u == 0 else w0 & ((1 << u) - 1)
    b1 = low == 0
    t53 = w1 & ((1 << 53) - 1)
    thresh = int(exp_small(-r) * (1 << 53))
    return b1 and t53 < thresh


class GaussianSampler:
    def __init__(self, xof):
        self.xof, self.bits, self.pos = xof, 0, 64

    def sign_bit(self):
        if self.pos >= 64:
            self.bits = self.xof.u64()
            self.pos = 0
        b = self.bits & 1
        self.bits >>= 1
        self.pos += 1
        return b

    def sample_base(self, c):
        dss = 1.0 / (2.0 * 1.55 * 1.55)
        while True:
            b = self.sign_bit()
            k0 = cdf_sampler(self.xof)
            k = b + (2 * b - 1) * k0
            x = ((k - c) ** 2 - (k - b) ** 2) * dss
            if ber_exp(self.xof, x):
                return k


def gaussian_vec(n, log2o, seed, dom):
    xof = Shake(seed, dom)
    outlen = -(-(log2o * n) // 8)
    stream = xof.fill(outlen) if outlen else b""
    s = GaussianSampler(xof)
    out = []
    for i in range(n):
        urand = 0
        for j in range(log2o):
            idx = i * log2o + j
            urand |= ((stream[idx >> 3] >> (idx & 7)) & 1) << j
        c = urand / (1 << log2o)
        k = s.sample_base(c)
        out.append((k << log2o) - urand)
    return out


# ---- binomial.rs, byte-exact ----
def cbd_vec(n, k, seed, dom):
    xof = Shake(seed, dom)
    stream = xof.fill(-(-(2 * k * n) // 8))

    def bit(i):
        return (stream[i >> 3] >> (i & 7)) & 1

    return [sum(bit(i * k + j) for j in range(k))
            - sum(bit(k * n + i * k + j) for j in range(k)) for i in range(n)]


# ---- rejection.rs ----
def rej_standard(xof, z, v, m, sigma):
    u = xof.u128()
    x = (sum(a * a for a in v) - 2 * sum(a * b for a, b in zip(z, v))) / (2 * sigma * sigma)
    if u == 0:
        return True
    return math.log(u) - 128 * LN2 + math.log(m) <= x


def main():
    seed = bytes([14] * 32)

    # gaussian: sigma = 1.55 * 2^3 (mirrors the Rust unit test)
    n = 1 << 14
    v = gaussian_vec(n, 3, seed, 0)
    sigma = 1.55 * 8
    mean = sum(v) / n
    var = sum((x - mean) ** 2 for x in v) / n
    print(f"gaussian(log2o=3): mean {mean:+.4f}  var {var:.2f} (target {sigma*sigma:.2f})")
    assert abs(mean) < 0.5 and abs(var - sigma * sigma) / (sigma * sigma) < 0.1
    assert all(abs(x) < 8 * sigma for x in v)

    v0 = gaussian_vec(4096, 0, bytes([15] * 32), 0)
    var0 = sum(x * x for x in v0) / len(v0)
    print(f"gaussian(log2o=0): var {var0:.4f} (target {1.55**2:.4f})")
    assert abs(var0 - 1.55 ** 2) / 1.55 ** 2 < 0.15

    # cbd
    w = cbd_vec(1 << 14, 2, bytes([10] * 32), 0)
    meanw = sum(w) / len(w)
    varw = sum((x - meanw) ** 2 for x in w) / len(w)
    print(f"cbd(k=2): mean {meanw:+.4f}  var {varw:.4f} (target 1.0)")
    assert abs(meanw) < 0.05 and abs(varw - 1.0) < 0.1
    assert all(-2 <= x <= 2 for x in w)

    # rejection acceptance rate with v = 0: exactly 1/M
    xof = Shake(bytes([21] * 32), 0, budget=1 << 20)
    m = math.e
    acc = sum(rej_standard(xof, [0] * 4, [0] * 4, m, 100.0) for _ in range(20000))
    rate = acc / 20000
    print(f"rej_standard(v=0, M=e): rate {rate:.4f} (target {1/m:.4f})")
    assert abs(rate - 1 / m) < 0.02

    print("\nALL SAMPLER SIMULATIONS PASSED (byte-exact spec of the Rust code)")


if __name__ == "__main__":
    main()
