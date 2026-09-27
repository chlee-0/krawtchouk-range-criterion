#!/usr/bin/env python3
"""Self-contained exact verifier for the simplified certificate in case (4, 7, 2, 7).

This file is intentionally human-readable: the Gaussian-integer matrices B_i
and the Gaussian-rational factor L are listed below as sparse coordinate data,
rather than hidden in a compressed/base64 payload.

The notation matches the manuscript:
  A(x) = sum_i x_i A_i, with the band offsets (0, 1, 3, 5);
  mathcal B = vertical_stack(B_1, B_2, B_3, B_4);
  R = L L*;
  delta = 1/5000.

Sparse B tuples are (i, row, col, real, imag), using 1-based indices, and mean
(B_i)[row,col] = real + imag*j.  Sparse L tuples are
(row, col, real, imag), also 1-based, for the numerator of L; the common
denominator is L_DEN.

All rank and positivity checks use exact integer arithmetic (Bareiss/Sylvester).
Standard library only.
"""
from __future__ import annotations

import hashlib
import json
import time
from math import lcm

if not __debug__:
    raise RuntimeError("Run without python -O")

M, N_SECOND, K, ELL = (4, 7, 2, 7)
assert M == 4
OFFSETS = (0, 1, 3, 5)
DELTA = (1, 5000)
L_DEN = 500
L_COLS = 9

# Human-readable exact certificate data.
B_ENTRIES = [
    (1, 1, 7, -1, 0),
    (1, 2, 4, 0, 1),
    (1, 2, 7, 0, 1),
    (1, 3, 2, 1, 0),
    (1, 3, 3, 0, 1),
    (1, 3, 4, -1, 0),
    (1, 4, 3, 0, 1),
    (1, 4, 5, -1, 0),
    (1, 4, 6, 0, 1),
    (1, 4, 7, 0, -1),
    (1, 5, 1, -1, 0),
    (1, 5, 4, 0, 1),
    (1, 5, 6, 0, 1),
    (1, 6, 6, 0, -1),
    (1, 6, 7, 1, -1),
    (1, 7, 2, -1, 0),
    (1, 7, 5, 1, 0),
    (2, 1, 4, 0, -1),
    (2, 1, 5, 0, -1),
    (2, 1, 7, 0, -1),
    (2, 2, 2, -1, 0),
    (2, 2, 4, 1, 0),
    (2, 2, 5, -1, 0),
    (2, 4, 1, 0, -1),
    (2, 4, 3, -1, 0),
    (2, 5, 2, -1, 1),
    (2, 5, 4, 0, -1),
    (2, 5, 5, 1, 0),
    (2, 5, 6, 0, 1),
    (2, 6, 3, 0, -1),
    (2, 6, 4, -1, 0),
    (2, 7, 2, 0, -1),
    (2, 7, 5, 0, -1),
    (2, 7, 6, 1, 0),
    (2, 7, 7, -1, 0),
    (3, 1, 5, 1, 0),
    (3, 1, 6, 0, -1),
    (3, 2, 1, 1, 1),
    (3, 2, 4, 0, -1),
    (3, 2, 6, 0, -1),
    (3, 3, 2, 1, 0),
    (3, 3, 3, -1, 0),
    (3, 3, 5, -1, 0),
    (3, 6, 1, 0, 1),
    (3, 6, 3, -1, -1),
    (3, 7, 1, 0, -1),
    (3, 7, 4, 1, 0),
    (4, 1, 2, 1, 0),
    (4, 1, 6, 0, 1),
    (4, 1, 7, 0, 1),
    (4, 2, 2, 1, 0),
    (4, 2, 3, 1, 0),
    (4, 2, 5, -1, 0),
    (4, 3, 2, 0, 1),
    (4, 3, 5, 0, 1),
    (4, 3, 7, 1, 0),
    (4, 4, 1, 0, -1),
    (4, 4, 3, 1, 0),
    (4, 4, 5, 1, 0),
    (4, 5, 1, 0, 1),
    (4, 5, 4, -1, 0),
]

L_NUM_ENTRIES = [
    (1, 1, 344, -34),
    (1, 2, 20, 60),
    (1, 3, 129, 108),
    (1, 4, -17, 29),
    (1, 5, -151, -1),
    (1, 6, -40, -72),
    (1, 7, -7, 109),
    (1, 8, 5, -31),
    (1, 9, 23, -1),
    (2, 1, -77, -17),
    (2, 2, 486, 0),
    (3, 1, 13, 5),
    (3, 2, 5, -21),
    (3, 3, 11, -23),
    (3, 4, -1, 8),
    (3, 5, -7, 18),
    (3, 6, -71, -70),
    (3, 7, -62, -88),
    (3, 8, -4, 57),
    (3, 9, 35, -5),
    (4, 1, -2, 41),
    (4, 2, -68, 8),
    (4, 3, -86, 18),
    (4, 4, 37, 18),
    (4, 5, 57, 111),
    (4, 6, -19, 29),
    (4, 7, -63, 9),
    (4, 8, 64, -18),
    (4, 9, -121, 0),
    (5, 1, -21, -17),
    (5, 2, -13, 27),
    (5, 3, -36, -20),
    (5, 4, 118, -31),
    (5, 5, -8, 41),
    (5, 6, -7, 45),
    (5, 7, -63, 12),
    (5, 8, 53, 46),
    (5, 9, -67, 18),
    (6, 1, -16, 44),
    (6, 2, -94, -4),
    (6, 3, -118, 93),
    (6, 4, 30, 45),
    (6, 5, 40, 68),
    (6, 6, 4, -5),
    (6, 7, -1, 51),
    (6, 8, 62, 38),
    (6, 9, 70, -21),
    (7, 1, -22, 22),
    (7, 2, -4, 3),
    (7, 3, 3, 0),
    (7, 4, 66, -27),
    (7, 5, -1, 10),
    (7, 6, 0, 4),
    (7, 7, 43, -33),
    (7, 8, 50, 18),
    (7, 9, 112, -3),
    (8, 1, 1, 16),
    (8, 2, -13, 32),
    (8, 3, -3, 7),
    (8, 4, 54, 58),
    (8, 5, 83, -32),
    (8, 6, 4, 55),
    (8, 7, -36, -10),
    (8, 8, -115, -26),
    (8, 9, 22, -52),
    (9, 1, 523, 0),
    (10, 1, -30, -15),
    (10, 2, 380, -23),
    (10, 3, -160, 74),
    (10, 4, -13, 82),
    (10, 5, 88, -11),
    (10, 6, -125, -16),
    (10, 7, -12, 70),
    (10, 8, -2, -59),
    (10, 9, 8, 71),
    (11, 1, -23, 17),
    (11, 2, 39, 20),
    (11, 3, 114, -15),
    (11, 4, 18, 47),
    (11, 5, -73, -15),
    (11, 6, -18, 12),
    (11, 7, -9, 13),
    (11, 8, -57, 29),
    (11, 9, -48, 53),
    (12, 1, 19, -59),
    (12, 2, -37, 35),
    (12, 3, -5, 47),
    (12, 4, -25, 44),
    (12, 5, 235, 0),
    (13, 1, -35, -27),
    (13, 2, -9, 42),
    (13, 3, 56, 104),
    (13, 4, -53, 47),
    (13, 5, 20, -48),
    (13, 6, 16, -8),
    (13, 7, -37, 88),
    (13, 8, 1, 6),
    (13, 9, -30, -31),
    (14, 1, 9, 78),
    (14, 2, -1, 43),
    (14, 3, 44, 4),
    (14, 4, 248, 0),
    (15, 1, 40, -16),
    (15, 2, -18, 26),
    (15, 3, -50, 16),
    (15, 4, -24, -58),
    (15, 5, 2, 25),
    (15, 6, -35, -79),
    (15, 7, -4, 34),
    (15, 8, 29, 101),
    (15, 9, -23, -18),
    (16, 1, -3, 34),
    (16, 2, -25, 72),
    (16, 3, 26, 81),
    (16, 4, 21, -24),
    (16, 5, -69, 38),
    (16, 6, -121, -5),
    (16, 7, -4, 12),
    (16, 8, 50, 75),
    (16, 9, 37, 35),
    (17, 1, 3, -9),
    (17, 2, 2, -1),
    (17, 3, 20, -3),
    (17, 4, -25, -3),
    (17, 5, -55, -40),
    (17, 6, 9, 30),
    (17, 7, -85, 29),
    (17, 8, -67, -82),
    (17, 9, -52, -56),
    (18, 1, 398, -23),
    (18, 2, 58, -10),
    (18, 3, 102, 36),
    (18, 4, 102, 87),
    (18, 5, -34, 12),
    (18, 6, 5, -69),
    (18, 7, -87, 97),
    (18, 8, 17, -49),
    (18, 9, -23, 48),
    (19, 1, 45, -17),
    (19, 2, 383, 1),
    (19, 3, -286, 0),
    (20, 1, -10, -77),
    (20, 2, -34, 15),
    (20, 3, -23, 74),
    (20, 4, -70, 25),
    (20, 5, -26, -34),
    (20, 6, 15, 18),
    (20, 7, -197, 0),
    (21, 1, -2, 17),
    (21, 2, -6, 25),
    (21, 3, -4, 6),
    (21, 4, -12, -13),
    (21, 5, -22, 16),
    (21, 6, 16, -9),
    (21, 7, 56, 127),
    (21, 8, -36, -64),
    (21, 9, -58, -88),
    (22, 1, 46, -8),
    (22, 2, -33, -55),
    (22, 3, -52, -42),
    (22, 4, 58, -5),
    (22, 5, -51, 61),
    (22, 6, 200, 0),
    (23, 1, 12, 0),
    (23, 2, -17, 18),
    (23, 3, -41, 9),
    (23, 4, -37, -20),
    (23, 5, -86, 28),
    (23, 6, 22, -58),
    (23, 7, -5, -12),
    (23, 8, 141, 0),
    (24, 1, -29, -11),
    (24, 2, -23, 38),
    (24, 3, -19, 16),
    (24, 4, 3, 28),
    (24, 5, -13, 18),
    (24, 6, -133, 117),
    (24, 7, 1, -10),
    (24, 8, 37, 62),
    (24, 9, -43, -37),
    (25, 1, -38, -23),
    (25, 2, 53, -7),
    (25, 3, 71, -2),
    (25, 4, -22, -10),
    (25, 5, -53, 31),
    (25, 6, 58, 33),
    (25, 7, 37, -14),
    (25, 8, 36, 66),
    (25, 9, 13, -38),
    (26, 1, -5, 16),
    (26, 2, 0, 10),
    (26, 3, -3, -16),
    (26, 4, 57, 3),
    (26, 5, 80, 8),
    (26, 6, -74, 45),
    (26, 7, -67, 19),
    (26, 8, -90, -47),
    (26, 9, -4, 47),
    (27, 1, 451, -71),
    (27, 2, 40, -28),
    (27, 3, 8, -19),
    (27, 4, 19, 100),
    (27, 5, 50, -3),
    (27, 6, -2, 9),
    (27, 7, -26, 52),
    (27, 8, 38, 58),
    (27, 9, -4, 20),
    (28, 1, -38, -15),
    (28, 2, 408, 19),
    (28, 3, -163, 73),
    (28, 4, 98, 1),
    (28, 5, 20, 56),
    (28, 6, 57, -51),
    (28, 7, -47, 33),
    (28, 8, 3, 73),
    (28, 9, -7, -58),
]

def build_band_A(n, k, offsets):
    """Return the four coefficient matrices A_i from the manuscript band pencil."""
    out = []
    for shift in offsets:
        a = [[0] * n for _ in range(k)]
        for row in range(k):
            col = row + shift
            assert 0 <= col < n
            a[row][col] = 1
        out.append(a)
    return out

def build_B(n, ell, entries):
    """Return real/imaginary parts of vertical_stack(B_1,...,B_4)."""
    N = 4 * n
    re = [[0] * ell for _ in range(N)]
    im = [[0] * ell for _ in range(N)]
    for block, row, col, a, b in entries:
        i = (block - 1) * n + (row - 1)
        j = col - 1
        assert 1 <= block <= 4 and 1 <= row <= n and 1 <= col <= ell
        assert re[i][j] == 0 and im[i][j] == 0
        re[i][j], im[i][j] = a, b
    return re, im

def build_L(n, cols, entries):
    """Return integer numerator real/imaginary parts of L."""
    N = 4 * n
    re = [[0] * cols for _ in range(N)]
    im = [[0] * cols for _ in range(N)]
    for row, col, a, b in entries:
        i, j = row - 1, col - 1
        assert 1 <= row <= N and 1 <= col <= cols
        assert re[i][j] == 0 and im[i][j] == 0
        re[i][j], im[i][j] = a, b
    return re, im

A_RE = build_band_A(N_SECOND, K, OFFSETS)
A_IM = [[[0] * N_SECOND for _ in range(K)] for _ in range(4)]
B_RE, B_IM = build_B(N_SECOND, ELL, B_ENTRIES)
L_RE, L_IM = build_L(N_SECOND, L_COLS, L_NUM_ENTRIES)

DATA = {
    "n": N_SECOND,
    "k": K,
    "ell": ELL,
    "A_re": A_RE,
    "A_im": A_IM,
    "B_re": B_RE,
    "B_im": B_IM,
    "L_re": L_RE,
    "L_im": L_IM,
    "L_den": L_DEN,
    "delta": list(DELTA),
}


def transpose(x):
    return [list(t) for t in zip(*x)]

def product(ar, ai, br, bi):
    n, k, m = len(ar), len(br), len(br[0])
    cr = [[sum(ar[i][t] * br[t][j] - ai[i][t] * bi[t][j] for t in range(k))
           for j in range(m)] for i in range(n)]
    ci = [[sum(ar[i][t] * bi[t][j] + ai[i][t] * br[t][j] for t in range(k))
           for j in range(m)] for i in range(n)]
    return cr, ci

def adjoint(ar, ai):
    return transpose(ar), [[-v for v in row] for row in transpose(ai)]

def gram_rows(ar, ai):
    return product(ar, ai, *adjoint(ar, ai))

def gram_cols(ar, ai):
    return product(*adjoint(ar, ai), ar, ai)

def pt(x, n):
    N = 4 * n
    return [[x[(j // n) * n + i % n][(i // n) * n + j % n]
             for j in range(N)] for i in range(N)]

def realify(ar, ai):
    n = len(ar)
    for i in range(n):
        for j in range(n):
            assert ar[i][j] == ar[j][i] and ai[i][j] == -ai[j][i], "not Hermitian"
    return [ar[i] + [-v for v in ai[i]] for i in range(n)] +            [ai[i] + ar[i] for i in range(n)]

def bareiss_pd(ar, ai):
    a = realify(ar, ai)
    n = len(a)
    prev = 1
    h = hashlib.sha256()
    maxbits = 0
    for k in range(n):
        p = a[k][k]
        if p <= 0:
            raise AssertionError(f"nonpositive leading principal determinant at {k + 1}")
        maxbits = max(maxbits, p.bit_length())
        h.update((str(p) + "\n").encode())
        if k == n - 1:
            break
        col = [a[i][k] for i in range(n)]
        for i in range(k + 1, n):
            for j in range(i, n):
                val, rem = divmod(p * a[i][j] - col[i] * col[j], prev)
                assert rem == 0, "Bareiss division not exact"
                a[i][j] = a[j][i] = val
        for i in range(k + 1, n):
            a[i][k] = a[k][i] = 0
        prev = p
    return {
        "real_size": n,
        "positive_leading_minors": n,
        "largest_minor_bits": maxbits,
        "minors_sha256": h.hexdigest(),
    }

def check_matrix(re, im, rows, cols):
    assert len(re) == len(im) == rows
    assert all(len(row) == cols for row in re + im)
    assert all(type(v) is int for row in re + im for v in row)

def verify_data(d):
    start = time.monotonic()
    n, k, ell = d["n"], d["k"], d["ell"]
    N = 4 * n
    ar, ai = d["A_re"], d["A_im"]
    br, bi = d["B_re"], d["B_im"]
    lr, li = d["L_re"], d["L_im"]
    q = d["L_den"]
    dn, dd = d["delta"]
    r = len(lr[0])

    for re_part, im_part in zip(ar, ai):
        check_matrix(re_part, im_part, k, n)
    check_matrix(br, bi, N, ell)
    check_matrix(lr, li, N, r)

    # The stacked matrix mathcal A = [A_1 A_2 A_3 A_4].
    aar = [[v for t in range(4) for v in ar[t][i]] for i in range(k)]
    aai = [[v for t in range(4) for v in ai[t][i]] for i in range(k)]

    agr, agi = gram_cols(aar, aai)
    bgr, bgi = gram_rows(br, bi)
    rr, ri = gram_rows(lr, li)

    fg, fi = pt(agr, n), pt(agi, n)
    rg, rgi = pt(rr, n), pt(ri, n)
    scale = lcm(q * q, dd)

    # scale * ((mathcal A*mathcal A)^Gamma + mathcal B mathcal B*
    #          - R^Gamma - delta I), with R = L L* and L = L_num/q.
    sr = [[
        scale * (fg[i][j] + bgr[i][j])
        - (scale // (q * q)) * rg[i][j]
        - (dn * (scale // dd) if i == j else 0)
        for j in range(N)
    ] for i in range(N)]
    si = [[
        scale * (fi[i][j] + bgi[i][j])
        - (scale // (q * q)) * rgi[i][j]
        for j in range(N)
    ] for i in range(N)]

    out = {
        "case": [4, n, k, ell],
        "B_nnz": sum(br[i][j] != 0 or bi[i][j] != 0
                     for i in range(N) for j in range(ell)),
        "L_shape": [N, r],
        "L_nnz": sum(lr[i][j] != 0 or li[i][j] != 0
                     for i in range(N) for j in range(r)),
        "L_den": q,
        "L_numerator_height": max(abs(v) for row in lr + li for v in row),
        "delta": [dn, dd],
    }
    out["A_full_row_rank"] = bareiss_pd(*gram_rows(aar, aai))
    out["B_full_column_rank"] = bareiss_pd(*gram_cols(br, bi))
    out["R_positive_semidefinite"] = "R = L L*, structural exact certificate"
    out["SOS_slack_positive_definite"] = bareiss_pd(sr, si)
    out["exact_verified"] = True
    out["seconds"] = time.monotonic() - start
    return out

def main():
    # These assertions make the manuscript/code correspondence explicit.
    assert (4, DATA["n"], DATA["k"], DATA["ell"]) == (4, 7, 2, 7)
    assert OFFSETS == (0, 1, 3, 5)
    assert DATA["delta"] == [1, 5000]
    assert L_DEN == 500
    result = verify_data(DATA)
    print(json.dumps(result, indent=2))
    print("EXACT PASS: 4727")

if __name__ == "__main__":
    main()
