#!/usr/bin/env python3
"""
Self-contained exact verification for the exceptional case (4,15,5,12), using the height-one Gaussian-integer witness.

The script prints the explicit data A_1,...,A_4, B_1,...,B_4,
L, R = L L^*, and delta, and then checks exactly the hypotheses
of the matrix sum-of-squares criterion used in the paper.

Standard library only. No floating point arithmetic or eigensolver is used.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd
import base64
import gzip
import json

CERTIFICATE_GZIP_BASE64 = """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"""

DATA = json.loads(gzip.decompress(base64.b64decode(CERTIFICATE_GZIP_BASE64)))
# Normalize the historical certificate field names to the notation used below.
if "Bstack_re" not in DATA:
    DATA["Bstack_re"] = DATA["U_re"]
    DATA["Bstack_im"] = DATA["U_im"]

n, k, ell = DATA["n"], DATA["k"], DATA["ell"]
N = 4 * n
d = DATA["denominator"]
dA = DATA["A_denominator"]
delta = Fraction(DATA["delta_num"], DATA["delta_den"])

A_re, A_im = DATA["A_re"], DATA["A_im"]
Bstack_re, Bstack_im = DATA["Bstack_re"], DATA["Bstack_im"]
L_re, L_im = DATA["L_re"], DATA["L_im"]

# Reduce the common denominator of the B_i for human-readable output.
gB = 0
for row in Bstack_re + Bstack_im:
    for z in row:
        gB = gcd(gB, abs(z))
gB = gcd(gB, d)
B_den = d // gB
Bstack_re_red = [[z // gB for z in row] for row in Bstack_re]
Bstack_im_red = [[z // gB for z in row] for row in Bstack_im]
B_re = [Bstack_re_red[i*n:(i+1)*n] for i in range(4)]
B_im = [Bstack_im_red[i*n:(i+1)*n] for i in range(4)]

def gram(re, im):
    """Return integer numerator matrices (Re, Im) of M M^*."""
    r, c = len(re), len(re[0])
    gre = [[sum(re[i][t]*re[j][t] + im[i][t]*im[j][t] for t in range(c))
            for j in range(r)] for i in range(r)]
    gim = [[sum(im[i][t]*re[j][t] - re[i][t]*im[j][t] for t in range(c))
            for j in range(r)] for i in range(r)]
    return gre, gim

R_re_num, R_im_num = gram(L_re, L_im)
R_den = d * d

def block_swap(X):
    return [[X[(j//n)*n + i % n][(i//n)*n + j % n] for j in range(N)]
            for i in range(N)]

def realify(re, im):
    m = len(re)
    for i in range(m):
        for j in range(m):
            if re[i][j] != re[j][i] or im[i][j] != -im[j][i]:
                raise AssertionError("matrix is not Hermitian")
    return ([re[i] + [-z for z in im[i]] for i in range(m)]
            + [im[i] + re[i] for i in range(m)])

def rank_over_q(rows):
    a = [row.copy() for row in rows]
    nr, nc = len(a), len(a[0])
    rank = row = 0
    for col in range(nc):
        piv = next((i for i in range(row, nr) if a[i][col] != 0), None)
        if piv is None:
            continue
        a[row], a[piv] = a[piv], a[row]
        for i in range(row + 1, nr):
            if a[i][col]:
                f, g = a[row][col], a[i][col]
                a[i] = [f*a[i][j] - g*a[row][j] for j in range(nc)]
        rank += 1
        row += 1
        if row == nr:
            break
    return rank

def complex_rank(re, im):
    # Realification [[Re,-Im],[Im,Re]] has twice the complex rank.
    M = [r + [-z for z in s] for r, s in zip(re, im)]
    M += [s + r for r, s in zip(re, im)]
    rr = rank_over_q(M)
    if rr % 2:
        raise AssertionError("realification has odd rank")
    return rr // 2

def bareiss_positive(a):
    """Check positive definiteness by Sylvester + fraction-free Bareiss."""
    a = [row.copy() for row in a]
    m = len(a)
    prev = 1
    pivots = []
    for q in range(m):
        p = a[q][q]
        if p <= 0:
            raise AssertionError(f"nonpositive leading principal minor at {q+1}")
        pivots.append(p)
        for i in range(q + 1, m):
            aiq = a[i][q]
            for j in range(i, m):
                num = p*a[i][j] - aiq*a[q][j]
                val, rem = divmod(num, prev)
                if rem:
                    raise AssertionError("Bareiss division not exact")
                a[i][j] = val
                a[j][i] = val
        prev = p
    return pivots

def fmt_gaussian(a, b):
    if b == 0:
        return str(a)
    if a == 0:
        if b == 1: return "i"
        if b == -1: return "-i"
        return f"{b}i"
    sign = "+" if b > 0 else "-"
    bb = abs(b)
    return f"{a}{sign}{'' if bb == 1 else bb}i"

def print_gaussian_matrix(name, re, im, den=1):
    prefix = f"{name} =" if den == 1 else f"{name} = (1/{den}) *"
    print(prefix)
    for rr, ii in zip(re, im):
        print("  [" + ", ".join(fmt_gaussian(a,b) for a,b in zip(rr,ii)) + "]")
    print()

def print_data():
    print("=" * 78)
    print("EXPLICIT CERTIFICATE DATA")
    print("=" * 78)
    print(f"(m,n,k,ell) = (4,{n},{k},{ell})")
    print(f"delta = {delta}\n")

    for t in range(4):
        print_gaussian_matrix(f"A_{t+1}", A_re[t], A_im[t], dA)

    for t in range(4):
        print_gaussian_matrix(f"B_{t+1}", B_re[t], B_im[t], B_den)

    print_gaussian_matrix("L", L_re, L_im, d)

    print(f"R = L L^* = (R_re + i R_im)/{R_den}")
    print("R_re =")
    for row in R_re_num:
        print("  [" + ", ".join(map(str,row)) + "]")
    print("R_im =")
    for row in R_im_num:
        print("  [" + ", ".join(map(str,row)) + "]")
    print()

def verify():
    checks = []

    # calA = [A_1 A_2 A_3 A_4].
    calA_re = [[A_re[t][a][c] for t in range(4) for c in range(n)]
               for a in range(k)]
    calA_im = [[A_im[t][a][c] for t in range(4) for c in range(n)]
               for a in range(k)]
    rankA = complex_rank(calA_re, calA_im)
    checks.append(("rank(calA) = k", rankA == k, f"rank = {rankA}, k = {k}"))

    # calB is the vertical stack of B_1,...,B_4.
    rankB = complex_rank(Bstack_re, Bstack_im)
    checks.append(("rank(calB) = ell", rankB == ell,
                   f"rank = {rankB}, ell = {ell}"))

    hermR = all(R_re_num[i][j] == R_re_num[j][i]
                and R_im_num[i][j] == -R_im_num[j][i]
                for i in range(N) for j in range(N))
    checks.append(("R is Hermitian", hermR, "R = L L^*"))
    checks.append(("R is positive semidefinite", True,
                   "structurally certified by R = L L^*"))
    checks.append(("delta > 0", delta > 0, f"delta = {delta}"))

    # F = calA^* calA, with denominator dA^2.
    F_re = [[sum(calA_re[a][i]*calA_re[a][j] + calA_im[a][i]*calA_im[a][j]
                 for a in range(k)) for j in range(N)] for i in range(N)]
    F_im = [[sum(calA_re[a][i]*calA_im[a][j] - calA_im[a][i]*calA_re[a][j]
                 for a in range(k)) for j in range(N)] for i in range(N)]
    Fg_re, Fg_im = block_swap(F_re), block_swap(F_im)

    BB_re, BB_im = gram(Bstack_re, Bstack_im)
    Rg_re, Rg_im = block_swap(R_re_num), block_swap(R_im_num)

    # Clear denominators in
    # (calA^*calA)^Gamma + calB calB^* - R^Gamma - delta I_N.
    q, p = DATA["delta_den"], DATA["delta_num"]
    dd, aa = d*d, dA*dA
    S_re = [[q*dd*Fg_re[i][j] + q*aa*BB_re[i][j] - q*aa*Rg_re[i][j]
             - (p*dd*aa if i == j else 0)
             for j in range(N)] for i in range(N)]
    S_im = [[q*dd*Fg_im[i][j] + q*aa*BB_im[i][j] - q*aa*Rg_im[i][j]
             for j in range(N)] for i in range(N)]

    pivots = bareiss_positive(realify(S_re, S_im))
    checks.append(("SOS inequality",
                   len(pivots) == 2*N,
                   f"all {len(pivots)} leading principal minors of the "
                   f"{2*N}x{2*N} realification are positive"))

    print("=" * 78)
    print("EXACT VERIFICATION")
    print("=" * 78)
    ok = True
    for label, passed, detail in checks:
        ok = ok and passed
        print(f"[{'PASS' if passed else 'FAIL'}] {label}: {detail}")
    print()
    print("FINAL RESULT:", "PASS" if ok else "FAIL")
    if not ok:
        raise SystemExit(1)

if __name__ == "__main__":
    print_data()
    verify()
