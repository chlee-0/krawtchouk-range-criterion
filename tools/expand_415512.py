#!/usr/bin/env python3
"""One-shot converter from the legacy compressed certificate to readable sparse data."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from pprint import pformat

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "verify_415512.py"

spec = importlib.util.spec_from_file_location("legacy_415512", TARGET)
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)
d = mod.DATA["415512"]

n, k, ell = d["n"], d["k"], d["ell"]
assert (n, k, ell) == (15, 5, 12)
offsets = (0, 1, 6, 10)
assert d["delta"] == [1, 20000]
assert d["L_den"] == 2000

def band_A(n, k, offsets):
    out = []
    for shift in offsets:
        a = [[0] * n for _ in range(k)]
        for row in range(k):
            a[row][row + shift] = 1
        out.append(a)
    return out

assert d["A_re"] == band_A(n, k, offsets)
assert all(v == 0 for mat in d["A_im"] for row in mat for v in row)

B_entries = []
for stacked_row, (rr, ii) in enumerate(zip(d["B_re"], d["B_im"])):
    block, row0 = divmod(stacked_row, n)
    for col0, (a, b) in enumerate(zip(rr, ii)):
        if a or b:
            B_entries.append((block + 1, row0 + 1, col0 + 1, a, b))

L_entries = []
for row0, (rr, ii) in enumerate(zip(d["L_re"], d["L_im"])):
    for col0, (a, b) in enumerate(zip(rr, ii)):
        if a or b:
            L_entries.append((row0 + 1, col0 + 1, a, b))

assert len(B_entries) == 446
assert len(L_entries) == 1200

source = TARGET.read_text()
tail_start = source.index("def transpose")
main_start = source.index("def main():", tail_start)
common_tail = source[tail_start:main_start]

header = """#!/usr/bin/env python3
\"\"\"Self-contained exact verifier for the simplified certificate in case (4, 15, 5, 12).

This file is intentionally human-readable: the Gaussian-integer matrices B_i
and the Gaussian-rational factor L are listed below as sparse coordinate data,
rather than hidden in a compressed/base64 payload.

The notation matches the manuscript:
  A(x) = sum_i x_i A_i, with the band offsets (0, 1, 6, 10);
  mathcal B = vertical_stack(B_1, B_2, B_3, B_4);
  R = L L*;
  delta = 1/20000.

Sparse B tuples are (i, row, col, real, imag), using 1-based indices, and mean
(B_i)[row,col] = real + imag*j. Sparse L tuples are
(row, col, real, imag), also 1-based, for the numerator of L; the common
denominator is L_DEN.

All rank and positivity checks use exact integer arithmetic (Bareiss/Sylvester).
Standard library only.
\"\"\"
from __future__ import annotations

import hashlib
import json
import time
from math import lcm

if not __debug__:
    raise RuntimeError("Run without python -O")

M, N_SECOND, K, ELL = (4, 15, 5, 12)
assert M == 4
OFFSETS = (0, 1, 6, 10)
DELTA = (1, 20000)
L_DEN = 2000
L_COLS = 25

# Human-readable exact certificate data.
"""

builders = """
def build_band_A(n, k, offsets):
    \"\"\"Return the four coefficient matrices A_i from the manuscript band pencil.\"\"\"
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
    \"\"\"Return real/imaginary parts of vertical_stack(B_1,...,B_4).\"\"\"
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
    \"\"\"Return integer numerator real/imaginary parts of L.\"\"\"
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

"""

footer = """
def main():
    # These assertions make the manuscript/code correspondence explicit.
    assert (4, DATA["n"], DATA["k"], DATA["ell"]) == (4, 15, 5, 12)
    assert OFFSETS == (0, 1, 6, 10)
    assert DATA["delta"] == [1, 20000]
    assert L_DEN == 2000
    result = verify_data(DATA)
    print(json.dumps(result, indent=2))
    print("EXACT PASS: 415512")

if __name__ == "__main__":
    main()
"""

text = (
    header
    + "B_ENTRIES = " + pformat(B_entries, width=96, sort_dicts=False) + "\n\n"
    + "L_NUM_ENTRIES = " + pformat(L_entries, width=96, sort_dicts=False) + "\n\n"
    + builders
    + common_tail
    + footer
)
TARGET.write_text(text)
print(f"wrote {TARGET} ({len(text)} chars); B nnz={len(B_entries)}, L nnz={len(L_entries)}")
