// Header-only port of the lrsnash driver (src/lib.cc, upstream lrsnashlib.c,
// derived from nash.c by Terje Lensberg). GPL v2 or later, see COPYING.
//
// Enumerates the vertices x of player 1's best-response polyhedron with lrs.
// For each x it enumerates the vertices y of the face of player 2's
// polyhedron given by the labels missing at x. Every (x, y) found is an
// extreme equilibrium and is passed to a visitor. The pivot sequence is the
// same as in upstream lrsnash, so the equilibria arrive in the same order.
//
// Differences from lib.cc that do not change which equilibria are produced,
// or their order:
//  * Vertex outputs are not gcd-reduced (reducearray). Upstream reduces
//    every player-1 vertex, but only uses it when an equilibrium is printed.
//    The visitor gets unreduced numerators and the common denominator.
//  * FirstTime, P2orig and linindex are members, not globals.
//  * The visitor can stop the enumeration early.
#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "lrs.hpp"

namespace LRSNash {
namespace detail {

template <class Ar> struct Nash {
  using T_ar = Ar;
  using T = typename Ar::T;
  using D = Dict<Ar>;
  using P = Problem<Ar>;

  Ar ar;
  P Q1, Q2;
  D *P1 = nullptr, *P2orig = nullptr;
  std::vector<long> linindex;
  std::vector<T> output1, output2;
  std::vector<long> num;
  std::vector<char> mark;
  std::vector<long> undo;  // warm-restart pivots (basis pos, cobasis pos)
  bool FirstTime = true;
  bool stop = false;

  // Visitor interface:
  //   bool on_eq2(const T *y, long n2)   y[0] = denominator, y[1..cols] = y,
  //                                      y[cols+1] = player 1's payoff.
  //                                      Return false to stop.
  //   bool on_eq1(const T *x, long n1)   x of the P1 vertex that produced at
  //                                      least one equilibrium (called after
  //                                      its on_eq2 calls), x[rows+1] =
  //                                      player 2's payoff. Return false to
  //                                      stop.
  // A(i, j) / B(i, j): payoffs of players 1 and 2 (0-based, i row, j col).
  template <class PayA, class PayB, class Visit>
  void solve(long rows, long cols, PayA A, PayB B, Visit &visit) {
    FirstTime = true;
    stop = false;
    num.resize(std::max(rows, cols) + 2);

    Q1.init(&ar, rows + cols + 1, rows + 2);
    P1 = Q1.lrs_alloc_dic();
    // BuildRepP1Is1: x_i >= 0, then v - (B^T x)_j >= 0, then sum x = 1
    FillNonnegativityRows(Q1, P1, 1, rows);
    for (long row = rows + 1; row < rows + 1 + cols; row++) {
      long s = row - (rows + 1);
      num[0] = 0;
      for (long t = 0; t < rows; t++) num[t + 1] = -(long)B(t, s);
      num[rows + 1] = 1;
      Q1.set_row(P1, row, num.data(), false);
    }
    FillLinearityRow(Q1, P1);
    FillFirstRow(Q1, P1);

    Q2.init(&ar, rows + cols + 1, cols + 2);
    P2orig = Q2.lrs_alloc_dic();
    // BuildRepP1Is0: u - (A y)_i >= 0, then y_j >= 0, then sum y = 1
    for (long row = 1; row < 1 + rows; row++) {
      long s = row - 1;
      num[0] = 0;
      for (long t = 0; t < cols; t++) num[t + 1] = -(long)A(s, t);
      num[cols + 1] = 1;
      Q2.set_row(P2orig, row, num.data(), false);
    }
    FillNonnegativityRows(Q2, P2orig, rows + 1, rows + cols);
    FillLinearityRow(Q2, P2orig);
    FillFirstRow(Q2, P2orig);

    linindex.assign(P2orig->m + P2orig->d + 2, 0);
    output1.resize(Q1.n + Q1.m);
    output2.resize(Q1.n + Q1.m);

    if (!Q1.lrs_getfirstbasis(P1))
      throw std::logic_error("lrs_getfirstbasis failed");

    long col = Q1.nredundcol;
    do {
      if (Q1.lrs_getsolution(P1, output1, col, false)) {
        long found = 0;
        nash2_main(visit, found);
        if (found > 0) {
          if (!origin(Q1, output1) && !visit.on_eq1(output1.data(), Q1.n))
            stop = true;
        }
        if (stop) break;
      }
    } while (Q1.lrs_getnextbasis(P1));
    Q1.free_ring(P1);
  }

  static bool origin(const P &Q, const std::vector<T> &out) {
    for (long i = 1; i < Q.n; i++)
      if (!Ar::zero(out[i])) return false;
    return true;
  }

  // ------------------------------------------------------------ input rows
  void FillNonnegativityRows(P &Q, D *Pd, long firstRow, long lastRow) {
    long n = Q.n;
    for (long row = firstRow; row <= lastRow; row++) {
      num[0] = 0;
      for (long col = 1; col < n; col++)
        num[col] = (row - firstRow + 1 == col) ? 1 : 0;
      Q.set_row(Pd, row, num.data(), false);
    }
  }
  void FillLinearityRow(P &Q, D *Pd) {
    long n = Q.n;
    num[0] = -1;
    for (long i = 1; i < n - 1; i++) num[i] = 1;
    num[n - 1] = 0;
    Q.set_row(Pd, Q.m, num.data(), true);
  }
  void FillFirstRow(P &Q, D *Pd) {
    for (long i = 0; i < Q.n; i++) num[i] = 1;
    Q.set_row(Pd, 0, num.data(), false);
  }

  // ------------------------------------------------------------ player 2
  //
  // Same decisions as lib.cc nash2_main, organised to avoid work in the
  // common case. After the first successful getabasis2 (FirstTime false),
  // upstream does: P2 = copy(P2orig); warm-restart pivots on P2;
  // P2orig = copy(P2); remove the linearity columns of P2; primalfeasible(P2).
  // In about 98% of calls (random 9x9 games) primalfeasible fails at its
  // first step, with no pivot. Here the warm-restart pivots are applied to
  // P2orig in place (undone if getabasis2 would fail, which leaves P2orig
  // unchanged as upstream does), and primalfeasible's first step is evaluated
  // read-only on P2orig with the linearity columns skipped. P2 is copied and
  // reduced only if the face may be nonempty. The dictionaries, and so every
  // later decision, are identical to upstream's.
  template <class Visit> void nash2_main(Visit &visit, long &numequilib) {
    compute_linearities();
    D *P2;
    if (FirstTime) {
      P2 = Q2.lrs_getdic();
      Q2.copy_dict(P2, P2orig);
      if (!lrs_getfirstbasis2(P2)) {
        Q2.free_ring(P2);
        return;
      }
    } else {
      if (!warm_restart_inplace()) return;    // getabasis2 FALSE
      if (precheck_infeasible()) return;      // primalfeasible FALSE, no pivot
      P2 = Q2.lrs_getdic();
      Q2.copy_dict(P2, P2orig);
      remove_linearities(P2);
      if (!finish_firstbasis2(P2)) {
        Q2.free_ring(P2);
        return;
      }
    }
    do {
      if (Q2.lrs_getsolution(P2, output2, 0, false)) {
        if (!origin(Q2, output2)) {
          numequilib++;
          if (!visit.on_eq2(output2.data(), Q2.n)) {
            stop = true;
            break;
          }
        }
      }
    } while (Q2.lrs_getnextbasis(P2));
    Q2.free_ring(P2);
  }

  // !FirstTime branch of getabasis2, applied to P2orig in place. Returns
  // false (with P2orig restored exactly) where upstream returns FALSE. The
  // integer pivot is an involution: pivoting again on the same entering /
  // leaving pair restores A, det, B, C, Row and Col exactly.
  bool warm_restart_inplace() {
    D *Pd = P2orig;
    long *B = Pd->B.data(), *C = Pd->C.data(), *Row = Pd->Row.data(),
         *Col = Pd->Col.data();
    const long *linearity = Q2.linearity.data();
    const long m = Pd->m, d = Pd->d, nlinearity = Q2.nlinearity;
    // upstream zeroes row 0 of the copy before getabasis2; after the first
    // successful call P2orig's row 0 is already zero, so this is a no-op
    for (long j = 0; j <= d; j++) ar.set(Pd->a(0, j), 0);
    for (long i = 1; i <= m + d; i++) linindex[i] = 0;
    for (long i = 0; i < nlinearity; i++) linindex[d + linearity[i]] = 1;
    undo.clear();
    for (long i = 1; i <= m; i++) {
      if (!linindex[B[i]]) continue;
      long k = 0;
      while (k < d && (linindex[C[k]] || ar.zero(Pd->a(Row[i], Col[k])))) k++;
      if (k < d) {
        long j = i;
        if (C[k] > B[j]) i--;
        Q2.pivot(Pd, j, k);
        Q2.update(Pd, &j, &k);
        undo.push_back(j);
        undo.push_back(k);
      } else if (!ar.zero(Pd->a(Row[i], 0))) {
        for (long u = (long)undo.size() - 2; u >= 0; u -= 2) {
          long uj = undo[u], uk = undo[u + 1];
          Q2.pivot(Pd, uj, uk);
          Q2.update(Pd, &uj, &uk);
        }
        return false;
      }
    }
    return true;
  }

  // primalfeasible's first step on the dictionary that getabasis2 would
  // produce from P2orig (cobasic linearity columns removed): true iff it
  // would return FALSE without pivoting. Removing columns does not change
  // the basis order, Row, or the order of the remaining cobasic columns.
  bool precheck_infeasible() {
    D *Pd = P2orig;
    const long m = Pd->m, d = Pd->d;
    const long *Row = Pd->Row.data(), *Col = Pd->Col.data(),
               *C = Pd->C.data();
    long i = Q2.lastdv + 1;
    while (i <= m && !ar.negative(Pd->a(Row[i], 0))) i++;
    if (i > m) return false;  // already primal feasible
    for (long k = 0; k < d; k++)
      if (!linindex[C[k]] && ar.positive(Pd->a(Row[i], Col[k])))
        return false;  // primalfeasible would pivot
    return true;
  }

  // tail of getabasis2: remove the cobasic linearity columns
  void remove_linearities(D *Pd) {
    const long *C = Pd->C.data();
    const long *linearity = Q2.linearity.data();
    long d = Pd->d;
    for (long i = 0; i < Q2.nlinearity; i++) {
      long k = 0;
      while (k < d && C[k] != linearity[i] + d) k++;
      if (k < d) {
        Q2.removecobasicindex(Pd, k);
        d = Pd->d;
      }
    }
  }

  // lrs_getfirstbasis2 after getabasis2. The renumbering of Q2.inequality is
  // skipped: Q2.inequality is only read in the FirstTime branch, where it is
  // rebuilt first.
  bool finish_firstbasis2(D *&Dp) {
    long i, j;
    const long m = Dp->m;
    const long nredundcol = Q2.nredundcol;
    long d = Dp->d;
    if (nredundcol > 0) {
      std::vector<T> lin(Q2.n);
      for (i = 0; i < nredundcol; i++) {
        Q2.lrs_getray(Dp, Dp->Col[0], Dp->C[0] + i, lin);
        if (!Q2.removecobasicindex(Dp, 0L)) return false;
      }
      d = Dp->d;
    }
    if (!Q2.primalfeasible(Dp)) return false;
    for (j = 1; j <= d; j++) {
      ar.copy(Dp->a(0, j), Dp->det);
      ar.storesign(Dp->a(0, j), -1);
    }
    ar.set(Dp->a(0, 0), 0);
    Dp->bound = ~(std::uint64_t)0;
    long *B = Dp->B.data(), *C = Dp->C.data();
    while (C[0] <= m) {
      i = C[0];
      C[0] = B[i];
      B[i] = i;
      P::reorder1(C, Dp->Col.data(), 0, d);
    }
    if (Q2.n - 1 > Dp->d) Dp = Q2.resize(Dp);
    return true;
  }

  void compute_linearities() {
    long *linearity = Q2.linearity.data();
    long nlinearity = 0;
    for (long i = Q1.lastdv + 1; i <= P1->m; i++) {
      if (!ar.zero(P1->a(P1->Row[i], 0))) {
        long j = Q1.inequality[P1->B[i] - Q1.lastdv];
        if (Q1.nlinearity == 0 || j < Q1.linearity[0])
          linearity[nlinearity++] = j;
      }
    }
    if (Q1.nlinearity > 0) linearity[nlinearity++] = Q1.linearity[0];
    // upstream calls reorder() nlinearity-1 times: a full sort of distinct
    // row indices
    std::sort(linearity, linearity + nlinearity);
    Q2.nlinearity = nlinearity;
  }

  // lib.cc lrs_getfirstbasis2 (Lin is not needed and dropped)
  bool lrs_getfirstbasis2(D *&Dp) {
    long i, j, k;
    long m = Dp->m, d = Dp->d;
    long nlinearity = Q2.nlinearity;
    long *linearity = Q2.linearity.data();
    long *inequality = Q2.inequality.data();
    long nredundcol;

    // linearities first, then the other rows m..1 (upstream does this with
    // a quadratic search)
    mark.assign(m + 1, 0);
    for (i = 0; i < nlinearity; i++) {
      inequality[i] = linearity[i];
      mark[linearity[i]] = 1;
    }
    k = nlinearity;
    for (i = m; i >= 1; i--)
      if (!mark[i]) inequality[k++] = i;
    for (j = 0; j <= d; j++) ar.set(Dp->a(0, j), 0);

    if (!getabasis2(Dp, inequality)) return false;

    nredundcol = Q2.nredundcol;
    d = Dp->d;

    for (i = 1; i <= m; i++) inequality[i] = i;
    if (nlinearity > 0) {
      for (i = 0; i < nlinearity; i++) inequality[linearity[i]] = 0;
      k = 1;
      for (i = 1; i <= m - nlinearity; i++) {
        while (k <= m && inequality[k] == 0) k++;
        inequality[i] = inequality[k++];
      }
    }
    if (nredundcol > 0) {
      std::vector<T> lin(Q2.n);
      for (i = 0; i < nredundcol; i++) {
        Q2.lrs_getray(Dp, Dp->Col[0], Dp->C[0] + i, lin);
        if (!Q2.removecobasicindex(Dp, 0L)) return false;
      }
    }
    if (!Q2.primalfeasible(Dp)) return false;

    for (j = 1; j <= d; j++) {
      ar.copy(Dp->a(0, j), Dp->det);
      ar.storesign(Dp->a(0, j), -1);
    }
    ar.set(Dp->a(0, 0), 0);
    Dp->bound = ~(std::uint64_t)0;

    long *B = Dp->B.data(), *C = Dp->C.data();
    while (C[0] <= m) {
      i = C[0];
      C[0] = B[i];
      B[i] = i;
      P::reorder1(C, Dp->Col.data(), 0, d);
    }
    if (Q2.n - 1 > Dp->d) Dp = Q2.resize(Dp);
    return true;
  }

  // lib.cc getabasis2
  bool getabasis2(D *Pd, long order[]) {
    long i, j, k;
    long *B = Pd->B.data(), *C = Pd->C.data(), *Row = Pd->Row.data(),
         *Col = Pd->Col.data();
    long *linearity = Q2.linearity.data();
    long *redundcol = Q2.redundcol.data();
    long m = Pd->m, d = Pd->d;
    long nlinearity = Q2.nlinearity;
    long nredundcol = 0;

    if (!FirstTime) {
      for (i = 1; i <= m + d; i++) linindex[i] = 0;
      for (i = 0; i < nlinearity; i++) linindex[d + linearity[i]] = 1;

      for (i = 1; i <= m; i++) {
        if (linindex[B[i]]) {
          k = 0;
          while (k < d && (linindex[C[k]] || ar.zero(Pd->a(Row[i], Col[k]))))
            k++;
          if (k < d) {
            j = i;
            if (C[k] > B[j]) i--;
            Q2.pivot(Pd, j, k);
            Q2.update(Pd, &j, &k);
          } else {
            if (!ar.zero(Pd->a(Row[i], 0))) return false;
          }
        }
      }
    } else {
      for (j = 0; j < m; j++) {
        i = 0;
        while (i <= m && B[i] != d + order[j]) i++;
        if (j < nlinearity && i > m) return false;
        if (i <= m) {
          k = 0;
          while (C[k] <= d && ar.zero(Pd->a(Row[i], Col[k]))) k++;
          if (C[k] <= d) {
            Q2.pivot(Pd, i, k);
            Q2.update(Pd, &i, &k);
          } else if (j < nlinearity) {
            if (ar.zero(Pd->a(Row[i], 0)))
              linearity[j] = 0;
            else
              return false;
          }
        }
      }
      i = 0;
      k = 0;
      while (k < nlinearity) {
        while (k < nlinearity && linearity[k] == 0) k++;
        if (k < nlinearity) linearity[i++] = linearity[k++];
      }
      nlinearity = i;
      Q2.nlinearity = i;
      k = 0;
      while (k < d && C[k] <= d) {
        if (C[k] <= d) redundcol[nredundcol++] = C[k];
        k++;
      }
      Q2.nredundcol = nredundcol;
      Q2.lastdv = d - nredundcol;
      FirstTime = false;
    }

    Q2.copy_dict(P2orig, Pd);

    for (i = 0; i < nlinearity; i++) {
      k = 0;
      while (k < d && C[k] != linearity[i] + d) k++;
      if (k < d) {
        Q2.removecobasicindex(Pd, k);
        d = Pd->d;
      }
    }
    return true;
  }
};

}  // namespace detail
}  // namespace LRSNash
