# Independent exact check of upstream lrsnash (gmp):
#  1. every reported equilibrium is an exact Nash equilibrium (Fractions)
#  2. for nondegenerate games, the reported set equals the set found by
#     independent support enumeration.
import random, subprocess, sys, itertools
from fractions import Fraction as F

def solve(M, rhs):
    # exact Gaussian elimination; returns None if singular
    n = len(M); A = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return None
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]

def support_enum(A, B):
    m, n = len(A), len(A[0]); out = set()
    for k in range(1, min(m, n) + 1):
        for I in itertools.combinations(range(m), k):
            for J in itertools.combinations(range(n), k):
                # y on J makes rows in I indifferent: sum_j A[i][j] y_j - u = 0, sum y = 1
                My = [[F(A[i][j]) for j in J] + [F(-1)] for i in I] + [[F(1)] * k + [F(0)]]
                ys = solve(My, [F(0)] * k + [F(1)])
                Mx = [[F(B[i][j]) for i in I] + [F(-1)] for j in J] + [[F(1)] * k + [F(0)]]
                xs = solve(Mx, [F(0)] * k + [F(1)])
                if ys is None or xs is None: continue
                if any(v <= 0 for v in ys[:k]) or any(v <= 0 for v in xs[:k]): continue
                x = [F(0)] * m; y = [F(0)] * n
                for a, i in enumerate(I): x[i] = xs[a]
                for a, j in enumerate(J): y[j] = ys[a]
                if is_ne(A, B, x, y): out.add((tuple(x), tuple(y)))
    return out

def is_ne(A, B, x, y):
    m, n = len(A), len(A[0])
    if any(v < 0 for v in x) or any(v < 0 for v in y) or sum(x) != 1 or sum(y) != 1: return False
    r = [sum(A[i][j] * y[j] for j in range(n)) for i in range(m)]
    c = [sum(B[i][j] * x[i] for i in range(m)) for j in range(n)]
    v1 = sum(x[i] * r[i] for i in range(m)); v2 = sum(y[j] * c[j] for j in range(n))
    return max(r) == v1 and max(c) == v2

def run_lrs(binary, A, B, path):
    m, n = len(A), len(A[0])
    with open(path, 'w') as f:
        f.write(f"{m} {n}\n\n" + "\n".join(" ".join(map(str, r)) for r in A) + "\n\n" + "\n".join(" ".join(map(str, r)) for r in B) + "\n")
    out = subprocess.run([binary, path], capture_output=True, text=True).stdout
    eqs, ys = [], []
    for line in out.splitlines():
        t = line.split()
        if not t: continue
        if t[0] == '2' and len(t) == n + 2: ys.append(([F(v) for v in t[1:n+1]], F(t[n+1])))
        elif t[0] == '1' and len(t) == m + 2:
            x = [F(v) for v in t[1:m+1]]; v2 = F(t[m+1])
            for y, v1 in ys: eqs.append((x, y, v1, v2))
            ys = []
    return eqs

binary, trials, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed); bad_ne = bad_set = dup = checked_set = total = 0
for t in range(trials):
    m, n = rng.randint(1, 6), rng.randint(1, 6)
    degenerate = t % 2 == 0
    hi = 3 if degenerate else 10**6
    A = [[rng.randint(-hi, hi) for _ in range(n)] for _ in range(m)]
    B = [[rng.randint(-hi, hi) for _ in range(n)] for _ in range(m)]
    if degenerate and m > 1 and rng.random() < .3: A[1] = A[0][:]; B[1] = B[0][:]
    eqs = run_lrs(binary, A, B, f"v{seed}.txt"); total += len(eqs)
    for x, y, v1, v2 in eqs:
        r = sum(A[i][j] * x[i] * y[j] for i in range(m) for j in range(n))
        c = sum(B[i][j] * x[i] * y[j] for i in range(m) for j in range(n))
        if not is_ne(A, B, x, y) or r != v1 or c != v2:
            bad_ne += 1; print("NOT NE", A, B, x, y, v1, v2)
    keys = [(tuple(x), tuple(y)) for x, y, _, _ in eqs]
    if len(set(keys)) != len(keys): dup += 1
    if not eqs: bad_set += 1; print("NO EQUILIBRIUM REPORTED", A, B)
    if not degenerate:
        checked_set += 1
        if set(keys) != support_enum(A, B):
            bad_set += 1; print("SET MISMATCH", A, B, sorted(set(keys)), sorted(support_enum(A, B)))
print(f"trials={trials} equilibria={total} not_NE={bad_ne} set_mismatch_or_empty={bad_set} (set-checked={checked_set}) games_with_duplicates={dup}")
