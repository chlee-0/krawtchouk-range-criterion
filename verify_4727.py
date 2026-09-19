#!/usr/bin/env python3
"""
Self-contained exact verification for the exceptional case (4,7,2,7).

The script prints the explicit data A_1,...,A_4, B_1,...,B_4,
L, R = L L^*, and delta, and then checks exactly the hypotheses
of the matrix sum-of-squares criterion used in the paper.

Standard library only. No floating point arithmetic or eigensolver is used.
"""
from __future__ import annotations
from fractions import Fraction
from math import gcd
import json

DATA = json.loads(r'''{"n":7,"k":2,"ell":7,"denominator":10000,"A_denominator":1,"delta_num":1,"delta_den":104,"A_re":[[[1,0,0,0,0,0,0],[0,1,0,0,0,0,0]],[[0,1,0,0,0,0,0],[0,0,1,0,0,0,0]],[[0,0,0,1,0,0,0],[0,0,0,0,1,0,0]],[[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]]],"A_im":[[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]],"Bstack_re":[[0,0,0,0,0,0,-10000],[0,0,0,0,0,0,0],[0,10000,0,-10000,0,0,0],[10000,0,0,0,-10000,0,0],[-10000,0,0,0,-10000,-10000,0],[0,0,0,0,0,0,10000],[0,-10000,0,0,10000,0,0],[0,0,0,0,0,0,0],[0,-10000,0,10000,-10000,0,0],[0,0,0,0,0,0,0],[-10000,0,-10000,0,0,0,0],[0,-10000,0,0,10000,0,0],[0,0,0,-10000,0,0,0],[0,-10000,0,0,0,10000,-10000],[0,0,0,0,10000,0,0],[10000,0,0,0,0,0,0],[0,10000,-10000,0,-10000,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,-10000,0,0,0,0],[0,0,0,10000,0,0,0],[0,10000,0,0,0,0,0],[0,10000,10000,0,-10000,0,0],[0,0,0,0,0,0,10000],[0,0,10000,0,10000,0,0],[0,0,0,-10000,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],"Bstack_im":[[0,0,0,0,0,0,0],[0,0,0,10000,0,0,10000],[0,0,10000,0,0,0,0],[0,0,10000,0,0,10000,-10000],[0,0,0,10000,0,10000,0],[0,0,0,0,0,-10000,-10000],[0,0,0,10000,0,0,0],[0,0,0,-10000,-10000,0,-10000],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[-10000,0,10000,0,0,0,0],[0,10000,0,-10000,0,10000,0],[0,0,-10000,0,0,0,0],[0,-10000,0,0,-10000,0,0],[0,0,0,0,0,-10000,0],[10000,0,0,-10000,0,-10000,0],[-10000,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,-10000,0],[10000,0,-10000,0,0,0,0],[-10000,0,0,0,0,0,0],[0,0,0,0,0,10000,10000],[0,0,0,0,0,0,0],[0,10000,0,0,10000,0,10000],[-10000,0,0,0,0,0,0],[10000,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],"L_re":[[7240,-2314,4711,-3874,946,-571,1210,-1131,647,-774,594,-880,169,-15],[-2280,-354,3272,-868,-1427,-2612,-130,590,167,-438,687,393,24,-2],[-141,1416,1302,601,1518,-1380,-1042,-245,-1434,270,632,664,548,3],[1229,1978,-4197,-912,-1622,-2215,1351,1917,1053,1269,-63,-128,0,0],[-687,3177,615,-818,-443,-2712,1471,-1644,292,1622,-1005,-51,179,0],[-1007,-154,-3212,-289,2493,1035,-71,-919,594,-635,575,-49,-8,0],[-669,-176,-68,-1004,-1908,-251,-1683,278,-1047,-1349,1771,597,205,-5],[106,640,46,3580,-388,-375,1093,-495,-402,-1079,85,1221,-74,-8],[9698,-3440,-2437,1018,139,-378,-1584,1644,564,-186,-1208,301,-322,-4],[-1685,-2476,611,1741,1124,-212,-642,1328,-1238,1112,-793,-129,81,-9],[309,102,4571,-375,-481,-862,1640,940,947,486,4,-118,-669,4],[1518,307,-778,2555,-95,-3787,381,212,1063,6,1593,-188,309,7],[-251,617,1508,-1031,753,638,987,528,1381,421,812,191,-130,-6],[953,2320,-604,2378,-2030,-765,177,711,-566,-778,-430,468,-15,-8],[1204,310,-447,-2752,1171,-275,-656,879,1156,-1291,-671,439,74,13],[481,1517,1272,-3636,2450,-446,-2883,1158,-670,-1380,920,137,67,5],[399,-177,730,-67,543,-60,1643,-1879,240,1091,-2699,-683,518,5],[5867,-4390,620,1343,591,1400,2513,-795,-2490,-545,1175,370,-69,8],[-1619,162,-3263,594,1398,1397,-1250,15,1036,-1012,-1011,-1092,364,-4],[-99,902,-392,-1458,2467,-1566,1986,-623,-1506,1169,314,-192,-115,2],[921,-231,-68,-900,-782,908,976,-544,1156,-330,-1081,-224,192,-1],[497,-226,-2355,-409,-1303,1215,1670,-949,2104,-603,547,-146,136,4],[497,891,-1457,-3518,178,2119,275,-270,-37,801,-23,472,56,-1],[-311,942,1210,-1797,1199,-2742,-1598,861,-137,236,-737,-710,62,-4],[-997,308,28,-1413,191,960,201,-154,-726,633,381,-643,-283,10],[-778,-47,1383,1927,-18,-2418,69,-2404,5,-1002,180,-232,-414,-4],[6474,-3703,-1926,1604,2656,-222,-90,43,787,1447,678,-183,-18,0],[-1805,712,-2185,52,-444,1066,1645,-491,482,-725,1505,-404,-10,0]],"L_im":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-3734,-6745,-1108,-864,315,-947,3520,1416,-696,-1456,-157,8,343,-3],[-346,-308,-143,-37,2003,-1911,-559,-807,-245,-1104,-1453,900,-358,10],[1141,-855,1574,2743,-2024,781,-749,-1829,1263,202,340,-619,-119,8],[616,-210,-2039,-2377,-3520,811,401,-1376,-1959,933,55,-242,-285,-3],[3075,6,2104,572,-1533,827,203,9,-2355,758,-23,-723,537,7],[617,1451,-1335,-1132,-978,-370,597,-232,-943,728,858,-399,257,-8],[73,846,-45,561,-74,2287,613,1623,790,898,257,605,11,-1],[444,-429,-414,-1033,421,924,1572,357,-87,-45,-364,1139,168,1],[-2955,-6968,1651,1986,-1322,133,-1274,1061,285,-543,-113,-804,662,2],[-201,-1138,-3148,110,556,-1386,-46,-2235,-693,1353,1054,617,303,4],[-2135,-153,3497,883,-399,1764,1662,1710,-1069,272,-160,399,-17,6],[144,-668,3281,-14,-870,950,1325,576,1615,162,-779,-155,-39,3],[2903,-654,-420,-2602,-4849,-1344,-1194,-791,-780,73,-1063,973,240,3],[-403,-1505,419,-2034,612,598,-1468,1494,-699,585,485,-818,197,5],[958,-932,684,-1094,-1678,1822,-2433,-1349,1274,-341,46,1177,219,4],[-428,-175,151,259,-244,-536,514,411,1879,435,-452,702,508,-4],[-224,-2208,-702,12,-2167,-1411,-1602,-371,1077,146,1079,-618,-159,3],[-3753,-6771,-1340,-1206,818,-402,-601,-1830,200,788,-32,415,-212,-2],[-907,85,1781,-1233,881,669,1366,-1209,2078,366,1336,1001,180,-1],[801,689,402,425,-287,301,-1226,1544,525,1182,-27,64,33,-3],[1560,767,-866,-1836,-351,-1809,2213,-487,438,-2624,-76,-472,262,2],[498,-581,-28,398,516,-16,103,-788,-918,-1849,-1551,986,-57,-1],[203,340,-469,-236,-1879,1003,117,279,175,-1404,240,-1245,-293,3],[467,-64,-2348,-58,-920,-2067,1663,3239,-44,-412,173,-252,-111,8],[-585,-639,816,897,516,-638,-1427,-2523,201,-291,-133,-984,-124,0],[-2219,-929,-2103,-1655,-547,593,530,-301,-1267,-786,-878,-663,-377,-6],[-2683,-7193,-65,-2208,-1006,-1505,-365,387,811,1375,-121,647,-418,0]]}''')

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
