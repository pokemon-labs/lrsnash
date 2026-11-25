// Header-only port of the lrsnash driver (src/lib.cc, upstream lrsnashlib.c,
// derived from nash.c by Terje Lensberg). GPL v2 or later, see COPYING.
//
// For a bimatrix game with payoffs A (player 1, rows) and B (player 2,
// columns), every extreme equilibrium is a pair of vertices of the two
// best-response polyhedra (Avis, Rosenberg, Savani, von Stengel 2010):
//   P1:  x >= 0,  v - (B^T x)_j >= 0 (j = 1..cols),  sum x = 1   (rows 1..m)
//   P2:  u - (A y)_i >= 0 (i = 1..rows),  y >= 0,  sum y = 1
// with the rows numbered so that row k of P1 and row k of P2 carry the same
// label: rows 1..rows are row strategy i (x_i = 0 / i is a best response),
// rows rows+1..rows+cols column strategy j, and row m = rows+cols+1 is the
// equation sum = 1. lrs enumerates the vertices x of P1. For each x, the
// labels x does not have (its slack rows with nonzero value) must be tight in
// P2, so the vertices y of that face of P2 are enumerated with those rows as
// extra equations. Every (x, y) found is an extreme equilibrium and is passed
// to a visitor. The pivot sequence is the same as in upstream lrsnash, so the
// equilibria arrive in the same order.
//
// Mapping to lib.cc: BuildRepP1Is1/BuildRepP1Is0 -> build_p1/build_p2,
// nash2_main -> enumerate_face, the linearity loop at its top ->
// face_linearities, lrs_getfirstbasis2 + getabasis2 -> p2_first_basis (first
// call) or warm_restart + precheck_infeasible + finish_basis (later calls),
// P2orig -> P2base, FirstTime -> !have_p2_base.
//
// Differences from lib.cc that do not change which equilibria are produced,
// or their order:
//  * Vertex outputs are not gcd-reduced (reducearray). Upstream reduces
//    every player-1 vertex, but only uses it when an equilibrium is printed.
//    The visitor gets unreduced numerators and the common denominator.
//  * The origin checks are dropped: sum x = 1 (sum y = 1) is an equation of
//    the system, so no vertex is the origin.
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
  using T_ar = Ar;  // used by solve.hpp
  using T = typename Ar::T;
  using D = Dict<Ar>;
  using P = Problem<Ar>;

  Ar ar;
  P Q1, Q2;
  D *P1 = nullptr;
  // P2's dictionary after the first successful basis search, with the
  // linearity columns still present; the starting point of every later
  // face (upstream P2orig).
  D *P2base = nullptr;
  bool have_p2_base = false;  // upstream !FirstTime
  std::vector<char> is_lin;   // is_lin[k]: slack k is a face equation
  std::vector<T> out1, out2;
  std::vector<long> undo;     // warm-restart pivots (basis pos, cobasis pos)
  bool stop = false;

  // Visitor interface:
  //   bool on_eq2(const T *y, long n2)   y[0] = denominator, y[1..cols] = y,
  //                                      y[cols+1] = player 1's payoff.
  //                                      Return false to stop.
  //   bool on_eq1(const T *x, long n1)   x of the P1 vertex that produced at
  //                                      least one equilibrium (called after
  //                                      its on_eq2 calls, also when one of
  //                                      them stopped), x[rows+1] = player
  //                                      2's payoff. Return false to stop.
  // A(i, j) / B(i, j): payoffs of players 1 and 2 (0-based, i row, j col).
  template <class PayA, class PayB, class Visit>
  void solve(long rows, long cols, PayA A, PayB B, Visit &visit) {
    have_p2_base = false;
    stop = false;
    build_p1(rows, cols, B);
    build_p2(rows, cols, A);
    is_lin.assign(Q2.m + Q2.nvars + 2, 0);
    out1.resize(Q1.n);
    out2.resize(Q2.n);

    if (!Q1.first_basis(P1)) throw std::logic_error("lrsnash: player 1 has no first basis");
    do {
      if (!Q1.get_vertex(P1, out1)) continue;
      if (enumerate_face(visit) > 0 && !visit.on_eq1(out1.data(), Q1.n)) stop = true;
      if (stop) break;
    } while (Q1.next_basis(P1));
    Q1.release_ring(P1);
  }

  // Player 1's current vertex x, in on_eq1's layout. It is valid during the
  // on_eq2 calls of its face too, as it is computed before the face.
  const T *vertex1() const { return out1.data(); }

  // pivots done by the last solve (both players), for tests and benchmarks
  long pivot_count() const { return Q1.pivots + Q2.pivots; }
  // dictionaries allocated so far (cumulative). A repeat solve of the same
  // shape reuses the pool and allocates none.
  std::size_t dict_allocations() const { return Q1.dicts_allocated + Q2.dicts_allocated; }

  // ------------------------------------------------------------ input rows
  // upstream BuildRepP1Is1: variables x_1..x_rows, v (column rows+1)
  template <class PayB> void build_p1(long rows, long cols, PayB &B) {
    const long m = rows + cols + 1;
    Q1.init(&ar, m, rows + 2);
    P1 = Q1.initial_dict();
    for (long r = 1; r <= rows; r++)  // x_r >= 0
      Q1.set_row(P1, r, [r](long c) { return c == r ? 1L : 0L; }, false);
    for (long s = 0; s < cols; s++)  // v - sum_t B(t, s) x_t >= 0
      Q1.set_row(P1, rows + 1 + s, [&, s](long c) {
        return c == 0 ? 0L : c <= rows ? -(long)B(c - 1, s) : 1L;
      }, false);
    set_sum_row(Q1, P1);
  }

  // upstream BuildRepP1Is0: variables y_1..y_cols, u (column cols+1)
  template <class PayA> void build_p2(long rows, long cols, PayA &A) {
    const long m = rows + cols + 1;
    Q2.init(&ar, m, cols + 2);
    P2base = Q2.initial_dict();
    for (long s = 0; s < rows; s++)  // u - sum_t A(s, t) y_t >= 0
      Q2.set_row(P2base, 1 + s, [&, s](long c) {
        return c == 0 ? 0L : c <= cols ? -(long)A(s, c - 1) : 1L;
      }, false);
    for (long r = 1; r <= cols; r++)  // y_r >= 0
      Q2.set_row(P2base, rows + r, [r](long c) { return c == r ? 1L : 0L; }, false);
    set_sum_row(Q2, P2base);
  }

  // row m: -1 + sum of the strategy variables = 0
  static void set_sum_row(P &Q, D *Pd) {
    const long last = Q.nvars;  // the payoff variable
    Q.set_row(Pd, Q.m, [last](long c) { return c == 0 ? -1L : c < last ? 1L : 0L; }, true);
  }

  // ------------------------------------------------------------ player 2
  // The face equations for the current P1 vertex: every P1 slack that is
  // basic with a nonzero value (a label x does not have), plus sum y = 1.
  // Upstream sorts them with nlinearity-1 passes of reorder(), which is a
  // full sort of these distinct row numbers.
  void face_linearities() {
    std::vector<long> &lin = Q2.linearity;
    lin.clear();
    const long *B = P1->B.data(), *Row = P1->Row.data();
    const long first_lin = Q1.linearity.empty() ? 0 : Q1.linearity[0];
    for (long i = Q1.nvars + 1; i <= P1->m; i++) {
      if (ar.zero(P1->a(Row[i], 0))) continue;
      const long r = Q1.slack_row[B[i] - Q1.nvars];
      if (Q1.linearity.empty() || r < first_lin) lin.push_back(r);
    }
    if (!Q1.linearity.empty()) lin.push_back(first_lin);
    std::sort(lin.begin(), lin.end());
  }

  // upstream nash2_main: enumerate the vertices y of the face of P2 for the
  // current P1 vertex; returns how many were reported.
  //
  // Upstream, after the first successful basis search, does for every face:
  // P2 = copy(P2orig); warm-restart pivots on P2; P2orig = copy(P2); remove
  // the linearity columns of P2; primalfeasible(P2). In about 98% of calls
  // (random 9x9 games) primalfeasible fails at its first step, with no pivot.
  // Here the warm-restart pivots are applied to P2base in place (undone if
  // upstream's getabasis2 would fail, which leaves P2orig unchanged as
  // upstream does), and primalfeasible's first step is evaluated read-only
  // on P2base with the linearity columns skipped. P2 is copied and reduced
  // only if the face may be nonempty. The dictionaries, and so every later
  // decision, are identical to upstream's.
  template <class Visit> long enumerate_face(Visit &visit) {
    face_linearities();
    D *P2;
    if (!have_p2_base) {
      P2 = Q2.new_ring();
      Q2.copy_dict(P2, P2base);
      if (!p2_first_basis(P2)) {
        Q2.release_ring(P2);
        return 0;
      }
    } else {
      if (!warm_restart()) return 0;          // upstream getabasis2 FALSE
      if (precheck_infeasible()) return 0;    // primalfeasible FALSE, no pivot
      P2 = Q2.new_ring();
      Q2.copy_dict(P2, P2base);
      Q2.remove_linearity_columns(P2);
      if (!Q2.finish_basis(P2, false)) {
        Q2.release_ring(P2);
        return 0;
      }
    }
    long found = 0;
    do {
      if (!Q2.get_vertex(P2, out2)) continue;
      ++found;
      if (!visit.on_eq2(out2.data(), Q2.n)) {
        stop = true;
        break;
      }
    } while (Q2.next_basis(P2));
    Q2.release_ring(P2);
    return found;
  }

  // lib.cc lrs_getfirstbasis2 with getabasis2's FirstTime branch: as
  // Problem::first_basis, but the dictionary is saved as P2base before the
  // linearity columns are removed (even if primalfeasible then fails). Q2's
  // slack_row is not maintained: only player 1's is read.
  bool p2_first_basis(D *&Dp) {
    Q2.start_order(Q2.slack_row.data());
    for (long j = 0; j <= Dp->d; j++) ar.set(Dp->a(0, j), 0);
    if (!Q2.find_basis(Dp, Q2.slack_row.data())) return false;
    have_p2_base = true;
    Q2.copy_dict(P2base, Dp);
    Q2.remove_linearity_columns(Dp);
    return Q2.finish_basis(Dp, false);
  }

  // getabasis2's !FirstTime branch, applied to P2base in place: pivot each
  // basic face-equation slack out, if a non-equation cobasic column allows
  // it. Returns false (with P2base restored exactly) where upstream returns
  // FALSE: a basic equation slack that cannot leave and is nonzero. The
  // integer pivot is an involution: pivoting again on the same entering /
  // leaving pair restores A, det, B, C, Row and Col exactly.
  //
  // P2base's row 0 is zero (it was saved before finish_basis wrote the
  // objective, and pivots keep a zero row zero), so upstream's zeroing of
  // row 0 here is a no-op and is skipped.
  bool warm_restart() {
    D *Pd = P2base;
    const long *B = Pd->B.data(), *C = Pd->C.data(), *Row = Pd->Row.data(),
               *Col = Pd->Col.data();
    const long m = Pd->m, d = Pd->d;  // d = nvars: no column was removed
    std::fill(is_lin.begin(), is_lin.end(), 0);
    for (long r : Q2.linearity) is_lin[d + r] = 1;
    undo.clear();
    for (long i = 1; i <= m; i++) {
      if (!is_lin[B[i]]) continue;
      long k = 0;
      while (k < d && (is_lin[C[k]] || ar.zero(Pd->a(Row[i], Col[k])))) k++;
      if (k < d) {
        long j = i;
        if (C[k] > B[j]) i--;  // the entering variable sorts below: revisit i
        Q2.pivot(Pd, j, k);
        undo.push_back(j);
        undo.push_back(k);
      } else if (!ar.zero(Pd->a(Row[i], 0))) {
        for (long u = (long)undo.size() - 2; u >= 0; u -= 2) {
          long uj = undo[u], uk = undo[u + 1];
          Q2.pivot(Pd, uj, uk);
        }
        return false;
      }
    }
    return true;
  }

  // primalfeasible's first step on the dictionary that remove_linearity_
  // columns would produce from P2base: true iff it would return false
  // without pivoting. Removing columns does not change the basis order,
  // Row, or the order of the remaining cobasic columns.
  bool precheck_infeasible() {
    const D *Pd = P2base;
    const long m = Pd->m, d = Pd->d;
    const long *Row = Pd->Row.data(), *Col = Pd->Col.data(), *C = Pd->C.data();
    long i = Q2.nvars + 1;
    while (i <= m && !ar.negative(Pd->a(Row[i], 0))) i++;
    if (i > m) return false;  // already primal feasible
    for (long k = 0; k < d; k++)
      if (!is_lin[C[k]] && ar.positive(Pd->a(Row[i], Col[k])))
        return false;  // primalfeasible would pivot
    return true;
  }
};

}  // namespace detail
}  // namespace LRSNash
