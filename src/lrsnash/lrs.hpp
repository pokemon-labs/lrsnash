// Header-only, global-free port of the parts of lrslib 7.2 used by lrsnash:
// lexicographic reverse search over the vertices of one polyhedron.
//
// Derived from lrslib.c / lrslib.h, Copyright David Avis 2005-2022 and
// contributors (Bremner's dictionary cache, Lensberg's nash driver).
// GPL v2 or later, see COPYING.
//
// ----------------------------------------------------------------- model
// A Problem is the system  b_r + sum_c a_rc x_c >= 0  (r = 1..m), with
// variables x_1..x_nvars, where nvars = n - 1, plus the equations listed in
// `linearity`. Each row r has a slack variable with index nvars + r. In a
// dictionary (Dict) every variable is either basic or cobasic:
//   B[0..m]   basic variable indices, sorted. B[0] = 0 is the objective row.
//   Row[i]    the row of A that holds the basic variable B[i].
//   C[0..d-1] cobasic variable indices, sorted; C[d] is a sentinel
//             (m + nvars + 1) that stops searches.
//   Col[j]    the column of A that holds the cobasic variable C[j].
//   A         (m + 1) x stride integers, row major. Row Row[i] is the
//             equation of B[i] in the cobasic variables; column 0 is the
//             right-hand side. Row 0 is the objective.
//   det       the common denominator: every entry of A is a numerator
//             over det (integer pivoting keeps A integral).
// d starts at nvars and drops by one for every column removed (cobasic
// linearity slacks). After the first basis, x_1..x_nvars are basic at
// positions 1..nvars (they are the smallest indices), so the current vertex
// is  x_k = A[Row[k]][0] / det.
//
// ------------------------------------------------- mapping to upstream
//   lrs_dic / lrs_dat             Dict / Problem (all file-level globals,
//                                 e.g. dict_count, live inside Problem)
//   lrs_alloc_dat, lrs_alloc_dic  init, initial_dict
//   lrs_getdic                    new_ring
//   copy_dict, pushQ, cache_dict  copy_dict, cache_push, cache_dict
//   check_cache                   restore_cached
//   resize                        reset_ring
//   lrs_free_dic                  release_ring
//   lrs_set_row_mp                set_row
//   pivot + update                pivot (always called as a pair upstream)
//   reorder1                      resort
//   lrs_ratio                     ratio_test
//   lexmin(P, 0)                  lexmin
//   reverse                       is_reverse_pivot
//   selectpivot                   select_pivot
//   primalfeasible                make_primal_feasible
//   removecobasicindex            remove_cobasic
//   getabasis                     find_basis + remove_linearity_columns
//   lrs_getfirstbasis             first_basis (= start_order, find_basis,
//                                 remove_linearity_columns, finish_basis)
//   lrs_getnextbasis              next_basis
//   lrs_getsolution(P, out, 0)    get_vertex (lrs_getvertex)
//   inequality                    slack_row
//   count[3]                      pivots
//
// ------------------------------------------- what is left out, and why
// Only the code reachable from lrsnash is kept: no hull, voronoi,
// nonnegative, maximize/minimize, bound, restart, givenstart or extract
// options, and nothing is printed. The following upstream code is dead for
// the lrsnash systems (see nash.hpp for their shape) and is removed:
//  * Redundant columns and rays (lrs_getray, the ismin branch of lexmin,
//    lrs_getsolution with col != 0). Both players' systems contain an
//    identity block (x_i >= 0, y_j >= 0), so the columns have full rank and
//    getabasis finds no redundant column; lastdv is always nvars.
//    find_basis checks this and throws std::logic_error if it ever fails.
//  * The gcd/lcm reduction of input rows (lrs_set_row_mp). Every lrsnash
//    row has an entry of +-1, so the row gcd is 1 and nothing is divided.
//  * The output gcd reduction (reducearray): lrsnash wants unreduced
//    numerators over det.
//  * depth and the MAXD depth limits (9.2e18, never reached), every count[]
//    except the pivot count, minratio[m], objnum/objden, m_A (always m) and
//    d_orig.
// Dictionaries come from a per-Problem pool, all at full size, and are
// reused across solves of the same shape.
#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "arith.hpp"

namespace LRSNash {
namespace detail {

template <class Ar> struct Dict {
  using T = typename Ar::T;
  static constexpr std::uint64_t kUnknownBound = ~(std::uint64_t)0;

  std::vector<T> A;
  long stride = 0;  // = nvars + 1, fixed for the Problem
  long m = 0, d = 0;
  T det;
  // magnitude hint for the fixed-width pivot kernels (see arith.hpp).
  // Reset to kUnknownBound after any direct write to A.
  std::uint64_t bound = kUnknownBound;
  bool lexmin = true;   // the basis is lexicographically minimal
  std::vector<long> B, Row, C, Col;
  long cache_r = 0, cache_s = 0;  // the pivot that left this dictionary
  Dict *prev = this, *next = this;  // cache ring

  Dict(long m_, long nvars)
      : A((m_ + 1) * (nvars + 1)), stride(nvars + 1), m(m_), d(nvars),
        B(m_ + 1), Row(m_ + 1), C(nvars + 1), Col(nvars + 1) {}

  T &a(long r, long c) { return A[r * stride + c]; }
  const T &a(long r, long c) const { return A[r * stride + c]; }
};

template <class Ar> struct Problem {
  using T = typename Ar::T;
  using D = Dict<Ar>;
  static constexpr unsigned long kCacheSize = 50;  // upstream dict_limit

  Ar *ar = nullptr;
  long m = 0;       // constraint rows
  long n = 0;       // columns of the input: right-hand side + nvars
  long nvars = 0;   // = n - 1 = upstream lastdv
  std::vector<long> linearity;  // rows that are equations, sorted
  // slack_row[k]: input row of the k-th slack after first_basis relabels the
  // slacks (upstream inequality). Only kept up to date by first_basis.
  std::vector<long> slack_row;
  long pivots = 0;

  // ---------------------------------------------------- dictionary storage
  // Cache ring (upstream Qhead/Qtail and the global dict_count): the
  // dictionaries of the last pivots, for backtracking without a reverse
  // pivot. Qtail is the current dictionary.
  D *Qhead = nullptr, *Qtail = nullptr;
  unsigned long dict_count = 0;
  std::vector<std::unique_ptr<D>> pool;
  std::vector<D *> free_list;
  std::size_t dicts_allocated = 0;  // cumulative, never reset (for tests)
  std::vector<long> ratio_buf;   // upstream minratio
  std::vector<char> mark;

  D *new_dict() {
    D *p;
    if (!free_list.empty()) {
      p = free_list.back();
      free_list.pop_back();
      p->d = nvars;
    } else {
      pool.push_back(std::make_unique<D>(m, nvars));
      ++dicts_allocated;
      p = pool.back().get();
    }
    p->prev = p->next = p;
    return p;
  }

  // Returns every dictionary of P's ring to the pool (upstream lrs_free_dic
  // stops at Qhead and can leak part of the ring).
  void release_ring(D *P) {
    if (!P) return;
    D *p = P;
    do {
      D *nx = p->next;
      free_list.push_back(p);
      p = nx;
    } while (p != P);
    Qhead = Qtail = nullptr;
  }

  // upstream lrs_alloc_dat + the problem part of lrs_alloc_dic
  void init(Ar *a, long m_, long n_) {
    ar = a;
    m = m_;
    n = n_;
    nvars = n - 1;
    pivots = 0;
    linearity.clear();
    linearity.reserve(m + 1);
    slack_row.assign(m + nvars + 1, 0);
    ratio_buf.assign(m + 1, 0);
    // Keep the pool across solves of the same shape. (The shape is m and
    // stride; d varies within a solve and is reset by new_dict.)
    if (!pool.empty() && (pool[0]->m != m || pool[0]->stride != n)) pool.clear();
    free_list.clear();
    for (auto &p : pool) free_list.push_back(p.get());
    Qhead = Qtail = nullptr;
  }

  // upstream lrs_alloc_dic: a zero dictionary with every slack basic, which
  // becomes the cache ring.
  D *initial_dict() {
    D *p = new_ring();
    ar->set(p->det, 1);
    p->lexmin = true;
    p->bound = D::kUnknownBound;
    for (long i = 0; i <= m; i++)
      for (long j = 0; j <= nvars; j++) ar->set(p->a(i, j), 0);
    for (long i = 0; i <= m; i++) {
      p->B[i] = (i == 0) ? 0 : nvars + i;
      p->Row[i] = i;
    }
    for (long j = 0; j < nvars; j++) {
      p->C[j] = j + 1;
      p->Col[j] = j + 1;
    }
    p->C[nvars] = m + nvars + 1;  // sentinel
    p->Col[nvars] = 0;
    return p;
  }

  // upstream lrs_getdic: a fresh dictionary that becomes the cache ring.
  // Upstream never resets the (global, shared with player 1) dict_count
  // here, which silently disables caching for player 2 once 50 dictionaries
  // exist. Every ring gets its own budget here; this only changes how
  // backtracking is done (restore instead of reverse pivot), not the result.
  D *new_ring() {
    D *p = new_dict();
    Qhead = Qtail = p;
    dict_count = 1;
    return p;
  }

  void copy_dict(D *dest, const D *src) {
    assert(dest != src);
    const long m_ = src->m, d = src->d;
    if constexpr (std::is_trivially_copyable_v<T>) {
      // one contiguous block; every dictionary of a Problem has the same stride
      std::memcpy(dest->A.data(), src->A.data(), sizeof(T) * (m_ + 1) * src->stride);
    } else {
      for (long r = 0; r <= m_; r++)
        for (long s = 0; s <= d; s++) ar->copy(dest->a(r, s), src->a(r, s));
    }
    dest->cache_r = src->cache_r;
    dest->cache_s = src->cache_s;
    dest->m = m_;
    dest->d = d;
    dest->lexmin = src->lexmin;
    dest->bound = src->bound;
    ar->copy(dest->det, src->det);
    std::memcpy(dest->B.data(), src->B.data(), sizeof(long) * (m_ + 1));
    std::memcpy(dest->C.data(), src->C.data(), sizeof(long) * (d + 1));
    std::memcpy(dest->Row.data(), src->Row.data(), sizeof(long) * (m_ + 1));
    std::memcpy(dest->Col.data(), src->Col.data(), sizeof(long) * (d + 1));
  }

  // upstream pushQ: advance Qtail, growing the ring up to kCacheSize and
  // then overwriting the oldest entry.
  void cache_push() {
    if (Qtail->next == Qhead) {
      if (dict_count < kCacheSize) {
        D *p = new_dict();
        p->next = Qtail->next;
        Qtail->next->prev = p;
        Qtail->next = p;
        p->prev = Qtail;
        dict_count++;
        Qtail = p;
      } else {
        Qhead = Qhead->next;
        Qtail = Qtail->next;
      }
    } else {
      Qtail = Qtail->next;
    }
  }

  // Save Dp (and the pivot r, s about to be done on it) and continue on a
  // copy.
  void cache_dict(D *&Dp, long r, long s) {
    Dp->cache_r = r;
    Dp->cache_s = s;
    cache_push();
    if (Qtail == Dp) return;
    copy_dict(Qtail, Dp);
    Dp = Qtail;
  }

  // Backtrack by restoring the previous dictionary, if it is cached.
  bool restore_cached(D *&Dp, long &r, long &s) {
    if (Qtail == Qhead) return false;
    Qtail = Qtail->prev;
    Dp = Qtail;
    r = Qtail->cache_r;
    s = Qtail->cache_s;
    return true;
  }

  // upstream resize(): every dictionary is full size, so only the cache ring
  // is reset to P.
  D *reset_ring(D *P) {
    for (D *p = P->next; p != P;) {
      D *nx = p->next;
      free_list.push_back(p);
      p = nx;
    }
    Qhead = Qtail = P;
    P->next = P->prev = P;
    return P;
  }

  // ------------------------------------------------------------ input rows
  // upstream lrs_set_row_mp for integer rows whose gcd is at most 1, which
  // upstream stores as given (it divides a row only by a gcd > 1; every
  // lrsnash row has an entry +-1): A[row][j] = num(j). eq: the row is an
  // equation.
  template <class Num> void set_row(D *P, long row, Num num, bool eq) {
#ifndef NDEBUG
    long g = 0;
    for (long j = 0; j <= nvars; j++) g = std::gcd(g, num(j));
    assert(g <= 1 && "set_row: upstream would divide this row by its gcd");
#endif
    P->bound = D::kUnknownBound;
    for (long j = 0; j <= nvars; j++) ar->set(P->a(row, j), num(j));
    if (eq) linearity.push_back(row);
  }

  // ------------------------------------------------------- basic steps
  // upstream reorder1: a[pos] has changed; move it to its sorted place in
  // a[0..len-1], moving b[] along with it. Returns that place.
  static long resort(long a[], long b[], long pos, long len) {
    while (pos > 0 && a[pos] < a[pos - 1]) {
      std::swap(a[pos], a[pos - 1]);
      std::swap(b[pos], b[pos - 1]);
      --pos;
    }
    while (pos < len - 1 && a[pos] > a[pos + 1]) {
      std::swap(a[pos], a[pos + 1]);
      std::swap(b[pos], b[pos + 1]);
      ++pos;
    }
    return pos;
  }

  // upstream pivot() + update(): exchange basic B[r] with cobasic C[s]. The
  // integer pivot is the arithmetic policy's kernel (exactly upstream's
  // result). On return r and s are the new positions of the entering and
  // leaving variables in B and C.
  void pivot(D *P, long &r, long &s) {
    ++pivots;
    ar->pivot(P->A.data(), P->stride, P->m + 1, P->d + 1, P->Row[r], P->Col[s], P->det,
              P->bound);
    long *B = P->B.data(), *C = P->C.data();
    const long leave = B[r], enter = C[s];
    B[r] = enter;
    r = resort(B, P->Row.data(), r, P->m + 1);  // >= 1: B[0] = 0 < enter
    C[s] = leave;
    s = resort(C, P->Col.data(), s, P->d);
  }

  // upstream lrs_ratio: the lexicographic minimum ratio test for entering
  // column col. Returns the basis position of the leaving variable, or 0 if
  // the column is unbounded.
  long ratio_test(const D *P, long col) {
    const long *B = P->B.data(), *Row = P->Row.data(), *Col = P->Col.data();
    long *cand = ratio_buf.data();  // candidate basis positions
    const long d = P->d;
    long ncand = 0;
    for (long j = nvars + 1; j <= P->m; j++) {  // branch-free append
      cand[ncand] = j;
      ncand += ar->negative(P->a(Row[j], col));
    }
    if (ncand == 0) return 0;
    // Break ties column by column: first the right-hand side, then the
    // cobasic columns in index order, skipping columns of basic variables.
    // Nmin/Dmin point into A, which is not changed here; each pass sets them
    // at its first candidate before reading them.
    const T *Nmin = nullptr, *Dmin = nullptr;
    long ratiocol = 0, start = 0, nstart = 0, nnext = 0;
    long bindex = d + 1, cindex = 0, basicindex = d;
    while (ncand > 1) {
      if (B[bindex] == basicindex) {
        if (cand[start] == bindex) {
          start++;
          ncand--;
        }
        bindex++;
      } else {
        if (basicindex != d) ratiocol = Col[cindex++];
        for (long j = start; j < start + ncand; j++) {
          const long i = Row[cand[j]];
          long comp = 1;
          if (j != start) {
            const T &Air = P->a(i, ratiocol);
            if constexpr (Ar::fixed_width) {
              // the sign tests below only skip this product, and give its
              // sign (both denominators are negative); for fixed-width
              // types the product is cheaper than the branches
              comp = ar->comprod(*Nmin, P->a(i, col), Air, *Dmin);
            } else if (ar->positive(*Nmin) || ar->negative(Air)) {
              if (ar->negative(*Nmin) || ar->positive(Air))
                comp = ar->comprod(*Nmin, P->a(i, col), Air, *Dmin);
              else
                comp = -1;
            } else if (ar->zero(*Nmin) && ar->zero(Air)) {
              comp = 0;
            }
            if (ratiocol == 0) comp = -comp;
          }
          if (comp == 1) {
            nstart = j;
            Nmin = &P->a(i, ratiocol);
            Dmin = &P->a(i, col);
            nnext = 1;
          } else if (comp == 0) {
            cand[nstart + nnext++] = cand[j];
          }
        }
        ncand = nnext;
        start = nstart;
      }
      basicindex++;
    }
    return cand[start];
  }

  // upstream lexmin(P, 0): is the basis lexicographically minimal for the
  // current vertex? (No degenerate row may allow a smaller index in.)
  bool lexmin(const D *P) const {
    const long *B = P->B.data(), *C = P->C.data(), *Row = P->Row.data(),
               *Col = P->Col.data();
    for (long i = nvars + 1; i <= P->m; i++) {
      const long r = Row[i];
      if (!ar->zero(P->a(r, 0))) continue;
      for (long j = 0; j < P->d; j++)
        if (B[i] > C[j] && !ar->zero(P->a(r, Col[j]))) return false;
    }
    return true;
  }

  // upstream reverse(): is pivoting C[s] in (leaving B[r]) the reverse of
  // the simplex step that leads back to this dictionary? Sets r.
  bool is_reverse_pivot(const D *P, long &r, long s) {
    const long *B = P->B.data(), *C = P->C.data(), *Row = P->Row.data(),
               *Col = P->Col.data();
    const long col = Col[s];
    if (!ar->negative(P->a(0, col))) return false;
    r = ratio_test(P, col);
    if (r == 0) return false;
    const long row = Row[r];
    for (long i = 0; i < P->d && C[i] < B[r]; i++) {
      if (i == s) continue;
      const long j = Col[i];
      const T &A0j = P->a(0, j), &Arj = P->a(row, j);
      if constexpr (Ar::fixed_width) {
        // the same test without short-circuit branches (comprod is cheap
        // and cannot throw for fixed-width types)
        const bool pos0 = ar->positive(A0j), neg0 = ar->negative(A0j);
        const bool posr = ar->positive(Arj), negr = ar->negative(Arj);
        const bool lt = ar->comprod(A0j, P->a(row, col), P->a(0, col), Arj) == -1;
        if ((pos0 | negr) & ((!neg0 & !posr) | lt)) return false;
      } else if (ar->positive(A0j) || ar->negative(Arj)) {
        if ((!ar->negative(A0j) && !ar->positive(Arj)) ||
            ar->comprod(A0j, P->a(row, col), P->a(0, col), Arj) == -1)
          return false;
      }
    }
    return true;
  }

  // upstream selectpivot: Bland's rule (smallest index with a positive cost).
  bool select_pivot(const D *P, long &r, long &s) {
    const long *Col = P->Col.data();
    const long d = P->d;
    r = 0;
    s = d;
    long j = 0;
    while (j < d && !ar->positive(P->a(0, Col[j]))) j++;
    if (j < d) {
      s = j;
      r = ratio_test(P, Col[j]);
      if (r != 0) return true;
    }
    return false;
  }

  // upstream primalfeasible: pivot until every slack's right-hand side is
  // nonnegative. False if some row cannot be repaired (infeasible).
  bool make_primal_feasible(D *P) {
    const long *Row = P->Row.data(), *Col = P->Col.data();
    for (;;) {
      long i = nvars + 1;
      while (i <= P->m && !ar->negative(P->a(Row[i], 0))) i++;
      if (i > P->m) return true;
      long j = 0;
      while (j < P->d && !ar->positive(P->a(Row[i], Col[j]))) j++;
      if (j >= P->d) return false;
      pivot(P, i, j);
    }
  }

  // upstream removecobasicindex: delete cobasic variable C[k] and its column
  // (a cobasic equation slack is 0 at every vertex). Higher indices shift
  // down by one; the last column of A moves into the freed one.
  void remove_cobasic(D *P, long k) {
    long *B = P->B.data(), *C = P->C.data(), *Col = P->Col.data();
    const long m_ = P->m, d = P->d;
    const long cindex = C[k], deloc = Col[k];
    for (long i = 1; i <= m_; i++)
      if (B[i] > cindex) B[i]--;
    for (long j = k; j < d; j++) {
      C[j] = C[j + 1] - 1;
      Col[j] = Col[j + 1];
    }
    if (deloc != d) {
      for (long i = 0; i <= m_; i++) ar->copy(P->a(i, deloc), P->a(i, d));
      long j = 0;
      while (Col[j] != d) j++;
      Col[j] = deloc;
    }
    P->d--;
  }

  // ------------------------------------------------------ the first basis
  // The order in which find_basis tries the slacks: the linearities, then
  // the other rows from m down to 1 (upstream lrs_getfirstbasis).
  void start_order(long *order) {
    mark.assign(m + 1, 0);
    long k = 0;
    for (long r : linearity) {
      order[k++] = r;
      mark[r] = 1;
    }
    for (long i = m; i >= 1; i--)
      if (!mark[i]) order[k++] = i;
  }

  // The first part of upstream getabasis: pivot every decision variable into
  // the basis, slacks in the given order, and drop linearities whose row is
  // 0 = 0. False if a linearity is inconsistent. Afterwards the remaining
  // linearities are all cobasic: pivots only bring decision variables in.
  bool find_basis(D *P, const long *order) {
    const long *B = P->B.data(), *C = P->C.data(), *Row = P->Row.data(),
               *Col = P->Col.data();
    const long m_ = P->m, d = P->d;
    const long nlin = (long)linearity.size();
    for (long j = 0; j < m_; j++) {
      long i = 0;
      while (i <= m_ && B[i] != d + order[j]) i++;
      if (j < nlin && i > m_) return false;
      if (i > m_) continue;
      long k = 0;
      while (C[k] <= d && ar->zero(P->a(Row[i], Col[k]))) k++;
      if (C[k] <= d) {
        pivot(P, i, k);
      } else if (j < nlin) {
        if (!ar->zero(P->a(Row[i], 0))) return false;
        linearity[j] = 0;  // redundant equation
      }
    }
    long kept = 0;
    for (long r : linearity)
      if (r != 0) linearity[kept++] = r;
    linearity.resize(kept);
    // a decision variable left cobasic would be a redundant column
    if (C[0] <= d)
      throw std::logic_error("lrsnash: rank-deficient system (redundant column)");
    return true;
  }

  // The tail of upstream getabasis: remove the columns of cobasic linearity
  // slacks (they are 0 at every vertex of the face). Slack r has index
  // nvars + r before any removal; each removal shifts the higher indices
  // down by one and lowers d by one, and the linearities are sorted, so the
  // slack of r is looked up as r + d with the current d (as upstream).
  void remove_linearity_columns(D *P) {
    for (long r : linearity) {
      const long *C = P->C.data();
      long k = 0;
      while (k < P->d && C[k] != r + P->d) k++;
      if (k < P->d) remove_cobasic(P, k);
    }
  }

  // The rest of upstream lrs_getfirstbasis: reach a feasible basis, set the
  // objective to minus the sum of the cobasic variables, and relabel the
  // cobasic slacks so that they have the smallest slack indices (this makes
  // the starting basis the root of the reverse search tree). If track is
  // true, slack_row is set up and follows the relabelling.
  bool finish_basis(D *&Dp, bool track) {
    const long d = Dp->d;
    long *ineq = slack_row.data();
    if (track) {
      for (long i = 1; i <= m; i++) ineq[i] = i;
      if (!linearity.empty()) {
        for (long r : linearity) ineq[r] = 0;
        long k = 1;
        for (long i = 1; i <= m - (long)linearity.size(); i++) {
          while (k <= m && ineq[k] == 0) k++;
          ineq[i] = ineq[k++];
        }
      }
    }
    if (!make_primal_feasible(Dp)) return false;
    for (long j = 1; j <= d; j++) {
      ar->copy(Dp->a(0, j), Dp->det);
      ar->storesign(Dp->a(0, j), -1);
    }
    ar->set(Dp->a(0, 0), 0);
    Dp->bound = D::kUnknownBound;
    long *B = Dp->B.data(), *C = Dp->C.data();
    while (C[0] <= m) {
      const long i = C[0];
      if (track) std::swap(ineq[B[i] - nvars], ineq[C[0] - nvars]);
      C[0] = B[i];
      B[i] = i;
      resort(C, Dp->Col.data(), 0, d);
    }
    if (nvars > Dp->d) Dp = reset_ring(Dp);
    return true;
  }

  // upstream lrs_getfirstbasis: from the initial dictionary to the root of
  // the reverse search. False if the system is infeasible.
  bool first_basis(D *&Dp) {
    start_order(slack_row.data());
    for (long j = 0; j <= Dp->d; j++) ar->set(Dp->a(0, j), 0);
    if (!find_basis(Dp, slack_row.data())) return false;
    remove_linearity_columns(Dp);
    return finish_basis(Dp, true);
  }

  // ------------------------------------------------------- reverse search
  // upstream lrs_getnextbasis: the next dictionary in reverse search order.
  // Backtracks by restoring a cached dictionary or, if it has been evicted,
  // by the forward simplex pivot.
  bool next_basis(D *&Dp) {
    long i = 0, j = 0;
    const long m_ = Dp->m, d = Dp->d;
    bool backtrack = false;
    while (j < d || Dp->B[m_] != m_) {
      if (backtrack) {
        backtrack = false;
        if (!restore_cached(Dp, i, j)) {
          select_pivot(Dp, i, j);
          pivot(Dp, i, j);
        }
        j++;
      }
      while (j < d && !is_reverse_pivot(Dp, i, j)) j++;
      if (j == d) {
        backtrack = true;
      } else {
        cache_dict(Dp, i, j);
        pivot(Dp, i, j);
        Dp->lexmin = lexmin(Dp);
        return true;
      }
    }
    return false;
  }

  // upstream lrs_getvertex without the gcd reduction: if the basis is
  // lexicographically minimal (each vertex is reported once), out[0] = det
  // and out[k] = numerator of x_k, k = 1..nvars.
  bool get_vertex(const D *P, std::vector<T> &out) const {
    if (!P->lexmin) return false;
    ar->copy(out[0], P->det);
    for (long k = 1; k <= nvars; k++) ar->copy(out[k], P->a(P->Row[k], 0));
    return true;
  }
};

}  // namespace detail
}  // namespace LRSNash
