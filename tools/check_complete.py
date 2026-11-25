#!/usr/bin/env python3
"""Completeness check for the extreme-equilibrium enumeration (nash.hpp).

For random small games, most of them degenerate, the equilibria that
nash.hpp reports (through tools/dump_eq) must be exactly the completely
labelled vertex pairs of the two best-response polyhedra. This script
computes those independently, by brute force in exact rational arithmetic.
It also checks that no pair is reported twice.

Polyhedra and labels, as in nash.hpp (labels 1..rows are row strategies,
rows+1..rows+cols column strategies):
  P1, variables (x_1..x_rows, v):  x_i >= 0 (label i),
      v - sum_t B[t][j] x_t >= 0 (label rows+j),  sum x = 1.
  P2, variables (y_1..y_cols, u):  u - sum_t A[i][t] y_t >= 0 (label i),
      y_j >= 0 (label rows+j),  sum y = 1.
A pair of vertices (x, v), (y, u) is an extreme equilibrium if and only if
every label is tight at x or at y.

  g++ -O2 -std=c++17 -Isrc tools/dump_eq.cc -lgmp -o dump_eq
  python3 -I tools/check_complete.py ./dump_eq trials maxdim seed
"""
import itertools
import random
import subprocess
import sys
from fractions import Fraction


def solve(M, rhs):
    """Exact solution of the square system M z = rhs, or None if singular."""
    n = len(M)
    a = [[Fraction(v) for v in row] + [Fraction(r)] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None:
            return None
        a[c], a[p] = a[p], a[c]
        piv = a[c][c]
        a[c] = [v / piv for v in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c]
                a[r] = [v - f * w for v, w in zip(a[r], a[c])]
    return [a[r][n] for r in range(n)]


def vertices(G, k):
    """Vertices of {z in Q^k : g.z >= 0 for every g in G, z_1+...+z_{k-1} = 1}.
    Constraint G[l] has label l + 1. Returns {vertex: frozenset of the labels
    tight there}. A vertex is the unique solution of the equation and k - 1
    tight constraints, so trying every (k-1)-subset finds them all."""
    eq = [1] * (k - 1) + [0]
    out = {}
    for S in itertools.combinations(range(len(G)), k - 1):
        z = solve([G[l] for l in S] + [eq], [0] * (k - 1) + [1])
        if z is None:
            continue
        vals = [sum(c * zi for c, zi in zip(g, z)) for g in G]
        if min(vals) < 0:
            continue
        out[tuple(z)] = frozenset(l + 1 for l, v in enumerate(vals) if v == 0)
    return out


def brute_force(m, n, A, B):
    P1 = ([[1 if t == i else 0 for t in range(m)] + [0] for i in range(m)] +
          [[-B[t][j] for t in range(m)] + [1] for j in range(n)])
    P2 = ([[-A[i][t] for t in range(n)] + [1] for i in range(m)] +
          [[1 if t == j else 0 for t in range(n)] + [0] for j in range(n)])
    V1, V2 = vertices(P1, m + 1), vertices(P2, n + 1)
    labels = frozenset(range(1, m + n + 1))
    eqs = {(x, y) for x, lx in V1.items() for y, ly in V2.items()
           if lx | ly == labels}
    # degenerate: some vertex has more tight labels than its dimension
    degenerate = (any(len(l) > m for l in V1.values()) or
                  any(len(l) > n for l in V2.values()))
    return eqs, degenerate


def game(rng, maxdim):
    m, n = rng.randint(1, maxdim), rng.randint(1, maxdim)
    kind = rng.randrange(6)
    hi = [1, 2, 3, 9, 100, 2][kind]
    A = [[rng.randint(0, hi) for _ in range(n)] for _ in range(m)]
    B = [[rng.randint(0, hi) for _ in range(n)] for _ in range(m)]
    if kind == 5:  # constant-sum
        B = [[hi - a for a in row] for row in A]
    if m > 1 and rng.random() < 0.3:  # duplicate row
        A[1], B[1] = A[0][:], B[0][:]
    if n > 1 and rng.random() < 0.3:  # duplicate column
        for i in range(m):
            A[i][1], B[i][1] = A[i][0], B[i][0]
    if rng.random() < 0.05:  # every profile is an equilibrium
        A = [[0] * n for _ in range(m)]
        B = [[0] * n for _ in range(m)]
    return m, n, A, B


def parse(text):
    games, cur = [], None
    for line in text.splitlines():
        t = line.split()
        if t[0] == 'G':
            cur = []
        elif t[0] == 'E':
            bar = t.index('|')
            xs = [int(v) for v in t[1:bar]]
            ys = [int(v) for v in t[bar + 1:]]
            cur.append((tuple(Fraction(v, xs[0]) for v in xs[1:]),
                        tuple(Fraction(v, ys[0]) for v in ys[1:])))
        elif t[0] == 'END':
            games.append(cur)
        else:
            raise SystemExit('unexpected dump_eq output: ' + line)
    return games


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    exe = sys.argv[1]
    trials, maxdim, seed = (int(v) for v in sys.argv[2:])
    rng = random.Random(seed)
    games = [game(rng, maxdim) for _ in range(trials)]
    lines = [str(trials)]
    for m, n, A, B in games:
        lines += [f'{m} {n}', ' '.join(str(v) for r in A for v in r),
                  ' '.join(str(v) for r in B for v in r)]
    run = subprocess.run([exe], input='\n'.join(lines) + '\n',
                         capture_output=True, text=True, check=True)
    reported = parse(run.stdout)
    if len(reported) != trials:
        raise SystemExit(f'dump_eq reported {len(reported)} of {trials} games')
    bad = dup = degenerate = total = 0
    for (m, n, A, B), got in zip(games, reported):
        want, degen = brute_force(m, n, A, B)
        degenerate += degen
        total += len(want)
        if len(got) != len(set(got)):
            dup += 1
        if set(got) != want:
            bad += 1
            if bad <= 3:
                print(f'MISMATCH {m}x{n} A={A} B={B}')
                print(f'  missing: {sorted(want - set(got))}')
                print(f'  extra:   {sorted(set(got) - want)}')
    print(f'games={trials} maxdim={maxdim} degenerate_games={degenerate} '
          f'extreme_equilibria={total} mismatches={bad} duplicates={dup}')
    return 1 if bad or dup else 0


if __name__ == '__main__':
    sys.exit(main())
