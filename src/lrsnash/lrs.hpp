// Header-only, global-free port of the parts of lrslib 7.2 used by lrsnash.
//
// Derived from lrslib.c / lrslib.h, Copyright David Avis 2005-2022 and
// contributors (Bremner's dictionary cache, Lensberg's nash driver).
// GPL v2 or later, see COPYING.
//
// Mapping to upstream:
//   lrs_dic  -> Dict<Ar>       (one dictionary; A stored row-major)
//   lrs_dat  -> Problem<Ar>    (problem data + dictionary cache + counters)
//   all file-level globals (dict_count, dict_limit, lrs_global_list, ...)
//   live inside Problem.
//
// Only the code paths reachable from lrsnash are kept. That means no hull,
// voronoi, nonnegative, maximize/minimize, bound, restart, givenstart or
// extract options, and nothing is printed. Function names follow upstream so
// that the two can be compared side by side.
//
// Notable simplifications that do not change results:
//  * objnum/objden are not maintained by pivot(). Upstream only prints them.
//  * Dictionaries come from a per-Problem pool and are never freed during a
//    solve. resize() is a no-op because every dictionary is allocated at the
//    full size.
#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>
#include <type_traits>
#include <vector>

#include "arith.hpp"

namespace LRSNash {
namespace detail {

template <class Ar> struct Dict {
  using T = typename Ar::T;
  std::vector<T> A;  // (m_A+1) x (stride) row major
  long stride = 0;
  long m = 0, m_A = 0, d = 0, d_orig = 0;
  long lexflag = 1, depth = 0, i = 0, j = 0;
  T det;
  // magnitude hint for the fixed-width pivot kernels (see arith.hpp);
  // ~0 means unknown. Set to ~0 after any direct write to A.
  std::uint64_t bound = ~(std::uint64_t)0;
  std::vector<long> B, Row, C, Col;
  Dict *prev = this, *next = this;

  Dict(long m_, long d_, long m_A_)
      : A((m_A_ + 1) * (d_ + 1)), stride(d_ + 1), m(m_), m_A(m_A_), d(d_),
        d_orig(d_), B(m_ + 1), Row(m_ + 1), C(d_ + 1), Col(d_ + 1) {}

  T &a(long r, long c) { return A[r * stride + c]; }
  const T &a(long r, long c) const { return A[r * stride + c]; }
  T *row(long r) { return &A[r * stride]; }
};

template <class Ar> struct Problem {
  using T = typename Ar::T;
  using D = Dict<Ar>;
  static constexpr long MAXD_DEPTH = 9223372036854775807L;  // upstream MAXD

  Ar *ar = nullptr;

  std::vector<T> Gcd, Lcm;
  std::vector<long> redundcol, inequality, linearity, minratio, temparray;
  long m = 0, n = 0, lastdv = 0;
  long count[10] = {0};
  long nredundcol = 0, nlinearity = 0;

  // dictionary cache (was Qhead/Qtail + file globals dict_count/dict_limit)
  D *Qhead = nullptr, *Qtail = nullptr;
  unsigned long dict_count = 0, dict_limit = 50;

  // storage
  std::vector<std::unique_ptr<D>> pool;
  std::vector<D *> free_list;

  // ---------------------------------------------------------------- memory
  D *new_lrs_dic() {
    if (!free_list.empty()) {
      D *p = free_list.back();
      free_list.pop_back();
      p->d = p->d_orig = n - 1;
      p->prev = p->next = p;
      return p;
    }
    pool.emplace_back(new D(m, n - 1, m));
    return pool.back().get();
  }

  // upstream lrs_free_dic: returns every dictionary of P's cache ring to the
  // pool (upstream stops at Qhead and can leak part of the ring).
  void free_ring(D *P) {
    if (!P) return;
    D *p = P;
    do {
      D *nx = p->next;
      free_list.push_back(p);
      p = nx;
    } while (p != P);
    Qhead = Qtail = nullptr;
  }

  // upstream lrs_alloc_dat + the problem-data part of lrs_alloc_dic.
  void init(Ar *a, long m_, long n_) {
    ar = a;
    m = m_;
    n = n_;
    nlinearity = 0;
    nredundcol = 0;
    for (long &c : count) c = 0;
    count[2] = 1;
    long d = n - 1;
    linearity.assign(m + d + 1, 0);
    inequality.assign(m + d + 1, 0);
    redundcol.assign(m + d + 1, 0);
    minratio.assign(m + d + 1, 0);
    temparray.assign(m + d + 1, 0);
    Gcd.resize(m + 1);
    Lcm.resize(m + 1);
    lastdv = d;
    // drop all dictionaries from previous solves of a different shape
    if (!pool.empty() && (pool[0]->m != m || pool[0]->d_orig != n - 1 ||
                          pool[0]->stride != n)) {
      pool.clear();
      free_list.clear();
    } else {
      free_list.clear();
      for (auto &p : pool) free_list.push_back(p.get());
    }
    Qhead = Qtail = nullptr;
  }

  // upstream lrs_alloc_dic (dictionary part). Creates a dictionary and makes
  // it the cache ring.
  D *lrs_alloc_dic() {
    D *p = new_lrs_dic();
    p->next = p->prev = p;
    Qhead = Qtail = p;
    dict_count = 1;
    dict_limit = 50;
    long d = n - 1;
    p->d = p->d_orig = d;
    p->m = m;
    p->m_A = m;
    p->depth = 0;
    p->lexflag = 1;
    ar->set(p->det, 1);
    p->bound = ~(std::uint64_t)0;
    for (long i = 0; i <= m; i++)
      for (long j = 0; j <= d; j++) ar->set(p->a(i, j), 0);
    for (long i = 0; i <= m; i++) {
      p->B[i] = (i == 0) ? 0 : d + i;
      p->Row[i] = i;
    }
    for (long j = 0; j < d; j++) {
      p->C[j] = j + 1;
      p->Col[j] = j + 1;
    }
    p->C[d] = m + d + 1;
    p->Col[d] = 0;
    return p;
  }

  // upstream lrs_getdic: fresh dictionary that becomes the cache ring.
  // Upstream never resets the (global, shared with player 1) dict_count
  // here, which silently disables caching for player 2 once 50 dictionaries
  // exist. We give every enumeration its own budget; this only changes how
  // backtracking is done (restore vs. reverse pivot), not the result.
  D *lrs_getdic() {
    D *p = new_lrs_dic();
    p->next = p->prev = p;
    Qhead = Qtail = p;
    dict_count = 1;
    return p;
  }

  void copy_dict(D *dest, D *src) {
    assert(dest != src);
    long m_ = src->m, m_A = src->m_A, d = src->d;
    if constexpr (std::is_trivially_copyable_v<T>) {
      // one contiguous block; every dictionary of a Problem has the same stride
      std::memcpy(dest->A.data(), src->A.data(),
                  sizeof(T) * (m_A + 1) * src->stride);
    } else {
      for (long r = 0; r <= m_A; r++)
        for (long s = 0; s <= d; s++) ar->copy(dest->a(r, s), src->a(r, s));
    }
    dest->i = src->i;
    dest->j = src->j;
    dest->m = m_;
    dest->d = d;
    dest->d_orig = src->d_orig;
    dest->m_A = src->m_A;
    dest->depth = src->depth;
    dest->lexflag = src->lexflag;
    dest->bound = src->bound;
    ar->copy(dest->det, src->det);
    std::memcpy(dest->B.data(), src->B.data(), sizeof(long) * (m_ + 1));
    std::memcpy(dest->C.data(), src->C.data(), sizeof(long) * (d + 1));
    std::memcpy(dest->Row.data(), src->Row.data(), sizeof(long) * (m_ + 1));
    std::memcpy(dest->Col.data(), src->Col.data(), sizeof(long) * (d + 1));
  }

  void pushQ() {
    if (Qtail->next == Qhead) {
      if (dict_count < dict_limit) {
        D *p = new_lrs_dic();
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

  void cache_dict(D *&Dp, long i, long j) {
    if (dict_limit > 1) {
      Dp->i = i;
      Dp->j = j;
      pushQ();
      if (Qtail == Dp) return;
      copy_dict(Qtail, Dp);
    }
    Dp = Qtail;
  }

  bool check_cache(D *&Dp, long *i_p, long *j_p) {
    if (Qtail == Qhead) return false;
    Qtail = Qtail->prev;
    Dp = Qtail;
    *i_p = Qtail->i;
    *j_p = Qtail->j;
    return true;
  }

  // upstream resize(): all dictionaries are full size, so only the cache
  // ring is reset to P (as upstream does after copying).
  D *resize(D *P) {
    if (Qhead != P || Qtail != P || P->next != P) {
      // return the other ring members to the pool
      D *p = P->next;
      while (p != P) {
        D *nx = p->next;
        free_list.push_back(p);
        p = nx;
      }
    }
    P->d_orig = P->d;
    Qhead = Qtail = P;
    P->next = P->prev = P;
    return P;
  }

  // ------------------------------------------------------- input rows
  // upstream lrs_set_row_mp specialised to integer input (all den == 1)
  void set_row(D *P, long row, const long *num, bool eq) {
    P->bound = ~(std::uint64_t)0;
    long d = P->d;
    long i = row;
    T &G = Gcd[i], &L = Lcm[i];
    T temp, mpone;
    ar->set(mpone, 1);
    ar->set(L, 1);
    ar->set(G, 0);
    for (long j = 0; j <= d; j++) {
      ar->set(P->a(i, j), num[j]);
      ar->copy(temp, P->a(i, j));
      ar->gcd(G, temp);
    }
    ar->storesign(G, 1);
    ar->storesign(L, 1);
    if (ar->greater(G, mpone) || ar->greater(L, mpone)) {
      T one_;
      ar->set(one_, 1);
      for (long j = 0; j <= d; j++) {
        ar->exactdiv(temp, P->a(i, j), G);
        ar->mul(temp, L, temp);
        ar->exactdiv(P->a(i, j), temp, one_);
      }
    }
    if (eq) {
      linearity[nlinearity] = row;
      nlinearity++;
    }
  }

  // --------------------------------------------------- core lrs routines
  static void reorder(long a[], long range) {
    long i, temp;
    for (i = 0; i < range - 1; i++)
      if (a[i] > a[i + 1]) {
        temp = a[i];
        a[i] = a[i + 1];
        a[i + 1] = temp;
      }
    for (i = range - 2; i >= 0; i--)
      if (a[i] > a[i + 1]) {
        temp = a[i];
        a[i] = a[i + 1];
        a[i + 1] = temp;
      }
  }

  static void reorder1(long a[], long b[], long newone, long range) {
    long temp;
    while (newone > 0 && a[newone] < a[newone - 1]) {
      temp = a[newone];
      a[newone] = a[newone - 1];
      a[newone - 1] = temp;
      temp = b[newone];
      b[newone] = b[newone - 1];
      b[--newone] = temp;
    }
    while (newone < range - 1 && a[newone] > a[newone + 1]) {
      temp = a[newone];
      a[newone] = a[newone + 1];
      a[newone + 1] = temp;
      temp = b[newone];
      b[newone] = b[newone + 1];
      b[++newone] = temp;
    }
  }

  // upstream pivot(): the integer pivot itself is the arithmetic policy's
  // dense kernel (Ar::pivot), which gives exactly the upstream result.
  void pivot(D *P, long bas, long cob) {
    count[3]++;
    ar->pivot(P->A.data(), P->stride, P->m_A + 1, P->d + 1, P->Row[bas],
              P->Col[cob], P->det, P->bound);
  }

  void update(D *P, long *i, long *j) {
    long *B = P->B.data(), *C = P->C.data();
    long m_ = P->m, d = P->d;
    long leave = B[*i], enter = C[*j];
    B[*i] = enter;
    reorder1(B, P->Row.data(), *i, m_ + 1);
    C[*j] = leave;
    reorder1(C, P->Col.data(), *j, d);
    for (*i = 1; B[*i] != enter; (*i)++)
      ;
    for (*j = 0; C[*j] != leave; (*j)++)
      ;
  }

  long lrs_ratio(D *P, long col) {
    long i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
    long firstime;
    long degencount, ndegencount;
    long *B = P->B.data(), *Row = P->Row.data(), *Col = P->Col.data();
    long *mr = minratio.data();
    long m_ = P->m, d = P->d;
    nstart = 0;
    ndegencount = 0;
    degencount = 0;
    mr[P->m] = 1;
    for (j = lastdv + 1; j <= m_; j++) {
      if (ar->negative(P->a(Row[j], col))) {
        mr[degencount++] = j;
        if (ar->zero(P->a(Row[j], 0))) mr[P->m] = 0;
      }
    }
    if (degencount == 0) return degencount;
    T *Nmin = nullptr, *Dmin = nullptr;  // point into A: values are copied
    T NminV, DminV;                      // upstream copies them
    ratiocol = 0;
    start = 0;
    bindex = d + 1;
    cindex = 0;
    basicindex = d;
    while (degencount > 1) {
      if (B[bindex] == basicindex) {
        if (mr[start] == bindex) {
          start++;
          degencount--;
        }
        bindex++;
      } else {
        firstime = 1;
        if (basicindex != d) ratiocol = Col[cindex++];
        for (j = start; j < start + degencount; j++) {
          i = Row[mr[j]];
          comp = 1;
          if (firstime)
            firstime = 0;
          else {
            const T &Air = P->a(i, ratiocol);
            if (ar->positive(*Nmin) || ar->negative(Air)) {
              if (ar->negative(*Nmin) || ar->positive(Air))
                comp = ar->comprod(*Nmin, P->a(i, col), Air, *Dmin);
              else
                comp = -1;
            } else if (ar->zero(*Nmin) && ar->zero(Air))
              comp = 0;
            if (ratiocol == 0) comp = -comp;
          }
          if (comp == 1) {
            nstart = j;
            ar->copy(NminV, P->a(i, ratiocol));
            ar->copy(DminV, P->a(i, col));
            Nmin = &NminV;
            Dmin = &DminV;
            ndegencount = 1;
          } else if (comp == 0)
            mr[nstart + ndegencount++] = mr[j];
        }
        degencount = ndegencount;
        start = nstart;
      }
      basicindex++;
    }
    return mr[start];
  }

  bool ismin(D *P, long r, long s) {
    long m_A = P->m_A;
    for (long i = 1; i <= m_A; i++)
      if ((i != r) && ar->negative(P->a(i, s)) &&
          ar->comprod(P->a(i, 0), P->a(r, s), P->a(i, s), P->a(r, 0)))
        return false;
    return true;
  }

  bool lexmin(D *P, long col) {
    long r, s, i, j;
    long *B = P->B.data(), *C = P->C.data(), *Row = P->Row.data(),
         *Col = P->Col.data();
    long m_ = P->m, d = P->d;
    for (i = lastdv + 1; i <= m_; i++) {
      r = Row[i];
      if (ar->zero(P->a(r, col)))
        for (j = 0; j < d; j++) {
          s = Col[j];
          if (B[i] > C[j]) {
            if (ar->zero(P->a(r, 0))) {
              if (!ar->zero(P->a(r, s))) return false;
            } else if (ar->negative(P->a(r, s)) && ismin(P, r, s)) {
              return false;
            }
          }
        }
    }
    return true;
  }

  bool reverse(D *P, long *r, long s) {
    long i, j, row, col;
    long *B = P->B.data(), *C = P->C.data(), *Row = P->Row.data(),
         *Col = P->Col.data();
    long d = P->d;
    col = Col[s];
    if (!ar->negative(P->a(0, col))) {
      minratio[P->m] = 0;
      return false;
    }
    *r = lrs_ratio(P, col);
    if (*r == 0) {
      minratio[P->m] = 0;
      return false;
    }
    row = Row[*r];
    for (i = 0; i < d && C[i] < B[*r]; i++)
      if (i != s) {
        j = Col[i];
        const T &A0j = P->a(0, j), &Arj = P->a(row, j);
        if (ar->positive(A0j) || ar->negative(Arj))
          if ((!ar->negative(A0j) && !ar->positive(Arj)) ||
              ar->comprod(A0j, P->a(row, col), P->a(0, col), Arj) == -1) {
            minratio[P->m] = 0;
            return false;
          }
      }
    return true;
  }

  bool selectpivot(D *P, long *r, long *s) {
    long j, col;
    long *Col = P->Col.data();
    long d = P->d;
    *r = 0;
    *s = d;
    j = 0;
    while ((j < d) && (!ar->positive(P->a(0, Col[j])))) j++;
    if (j < d) {
      *s = j;
      col = Col[j];
      *r = lrs_ratio(P, col);
      if (*r != 0) return true;
    }
    return false;
  }

  bool primalfeasible(D *P) {
    long i, j;
    long *Row = P->Row.data(), *Col = P->Col.data();
    long m_ = P->m, d = P->d;
    for (;;) {
      i = lastdv + 1;
      while (i <= m_ && !ar->negative(P->a(Row[i], 0))) i++;
      if (i <= m_) {
        j = 0;
        while (j < d && !ar->positive(P->a(Row[i], Col[j]))) j++;
        if (j >= d) return false;
        pivot(P, i, j);
        update(P, &i, &j);
      } else
        return true;
    }
  }

  bool removecobasicindex(D *P, long k) {
    long i, j, cindex, deloc;
    long *B = P->B.data(), *C = P->C.data(), *Col = P->Col.data();
    long m_ = P->m, d = P->d;
    cindex = C[k];
    deloc = Col[k];
    for (i = 1; i <= m_; i++)
      if (B[i] > cindex) B[i]--;
    for (j = k; j < d; j++) {
      C[j] = C[j + 1] - 1;
      Col[j] = Col[j + 1];
    }
    if (deloc != d) {
      for (i = 0; i <= m_; i++) ar->copy(P->a(i, deloc), P->a(i, d));
      j = 0;
      while (Col[j] != d) j++;
      Col[j] = deloc;
    }
    P->d--;
    return true;
  }

  bool getabasis(D *P, long order[]) {
    long i, j, k;
    long *B = P->B.data(), *C = P->C.data(), *Row = P->Row.data(),
         *Col = P->Col.data();
    long m_ = P->m, d = P->d, nlin = nlinearity;
    long nred = 0;
    for (j = 0; j < m_; j++) {
      i = 0;
      while (i <= m_ && B[i] != d + order[j]) i++;
      if (j < nlin && i > m_) return false;
      if (i <= m_) {
        k = 0;
        while (C[k] <= d && ar->zero(P->a(Row[i], Col[k]))) k++;
        if (C[k] <= d) {
          pivot(P, i, k);
          update(P, &i, &k);
        } else if (j < nlin) {
          if (ar->zero(P->a(Row[i], 0)))
            linearity[j] = 0;
          else
            return false;
        }
      }
    }
    i = 0;
    k = 0;
    while (k < nlin) {
      while (k < nlin && linearity[k] == 0) k++;
      if (k < nlin) linearity[i++] = linearity[k++];
    }
    nlin = i;
    nlinearity = i;
    k = 0;
    while (k < d && C[k] <= d) {
      if (C[k] <= d) redundcol[nred++] = C[k];
      k++;
    }
    nredundcol = nred;
    lastdv = d - nred;
    for (i = 0; i < nlin; i++) {
      k = 0;
      while (k < d && C[k] != linearity[i] + d) k++;
      if (k >= d) return false;
      if (!removecobasicindex(P, k)) return false;
      d = P->d;
    }
    return true;
  }

  void getnextoutput(D *P, long i, long col, T &out) {
    ar->copy(out, P->a(P->Row[i], col));
  }

  void reducearray(std::vector<T> &p, long n_) {
    T divisor, Temp;
    long i = 0;
    while ((i < n_) && ar->zero(p[i])) i++;
    if (i == n_) return;
    ar->copy(divisor, p[i]);
    ar->storesign(divisor, 1);
    i++;
    while (i < n_) {
      if (!ar->zero(p[i])) {
        ar->copy(Temp, p[i]);
        ar->storesign(Temp, 1);
        ar->gcd(divisor, Temp);
      }
      i++;
    }
    for (i = 0; i < n_; i++)
      if (!ar->zero(p[i])) {
        ar->copy(Temp, p[i]);
        ar->exactdiv(p[i], Temp, divisor);
      }
  }

  bool lrs_getvertex(D *P, std::vector<T> &output, bool reduce = true) {
    long i, ind, ired;
    long lexflag = P->lexflag;
    if (lexflag) {
      ++count[1];
      if (P->depth > count[8]) count[8] = P->depth;
    }
    if (!lexflag) return false;
    i = 1;
    ired = 0;
    ar->copy(output[0], P->det);
    for (ind = 1; ind < n; ind++)
      if ((ired < nredundcol) && (redundcol[ired] == ind)) {
        ar->set(output[ind], 0);
        ired++;
      } else {
        getnextoutput(P, i, 0, output[ind]);
        i++;
      }
    if (reduce) reducearray(output, n);
    if (lexflag && ar->one(output[0])) ++count[4];
    return true;
  }

  bool lrs_getray(D *P, long col, long redcol, std::vector<T> &output) {
    long i, ind, ired;
    if (redcol == n) ++count[0];
    i = 1;
    ired = 0;
    for (ind = 0; ind < n; ind++) {
      if (ind == 0)
        ar->set(output[0], 0);
      else if ((ired < nredundcol) && (redundcol[ired] == ind)) {
        if (redcol == ind)
          ar->copy(output[ind], P->det);
        else
          ar->set(output[ind], 0);
        ired++;
      } else {
        getnextoutput(P, i, col, output[ind]);
        i++;
      }
    }
    reducearray(output, n);
    return true;
  }

  bool lrs_getsolution(D *P, std::vector<T> &output, long col,
                       bool reduce = true) {
    long j;
    long *Row = P->Row.data();
    if (col == 0) return lrs_getvertex(P, output, reduce);
    if (!ar->negative(P->a(0, col))) return false;
    j = lastdv + 1;
    while (j <= P->m && !ar->negative(P->a(Row[j], col))) j++;
    if (j <= P->m) return false;
    if (lexmin(P, col)) {
      if (P->depth > count[8]) count[8] = P->depth;
      return lrs_getray(P, col, n, output);
    }
    return false;
  }

  // upstream lrs_getfirstbasis (Lin is not needed by nash and is dropped)
  bool lrs_getfirstbasis(D *&Dp) {
    long i, j, k;
    long m_ = Dp->m, d = Dp->d;
    long nlin = nlinearity;
    long *ineq = inequality.data();
    for (i = 0; i < nlin; i++) ineq[i] = linearity[i];
    k = nlin;
    for (i = m_; i >= 1; i--) {
      j = 0;
      while (j < k && ineq[j] != i) j++;
      if (j == k) ineq[k++] = i;
    }
    for (j = 0; j <= d; j++) ar->set(Dp->a(0, j), 0);
    if (!getabasis(Dp, ineq)) return false;
    nlin = nlinearity;
    long nred = nredundcol, ldv = lastdv;
    d = Dp->d;
    for (i = 1; i <= m_; i++) ineq[i] = i;
    if (nlin > 0) {
      for (i = 0; i < nlin; i++) ineq[linearity[i]] = 0;
      k = 1;
      for (i = 1; i <= m_ - nlin; i++) {
        while (k <= m_ && ineq[k] == 0) k++;
        ineq[i] = ineq[k++];
      }
    }
    if (nred > 0) {
      std::vector<T> lin(n);
      for (i = 0; i < nred; i++) {
        lrs_getray(Dp, Dp->Col[0], Dp->C[0] + i, lin);
        if (!removecobasicindex(Dp, 0)) return false;
      }
    }
    if (!primalfeasible(Dp)) return false;
    for (j = 1; j <= d; j++) {
      ar->copy(Dp->a(0, j), Dp->det);
      ar->storesign(Dp->a(0, j), -1);
    }
    ar->set(Dp->a(0, 0), 0);
    Dp->bound = ~(std::uint64_t)0;
    long *B = Dp->B.data(), *C = Dp->C.data();
    while (C[0] <= m_) {
      i = C[0];
      j = ineq[B[i] - ldv];
      ineq[B[i] - ldv] = ineq[C[0] - ldv];
      ineq[C[0] - ldv] = j;
      C[0] = B[i];
      B[i] = i;
      reorder1(C, Dp->Col.data(), 0, d);
    }
    if (n - 1 > Dp->d) Dp = resize(Dp);
    return true;
  }

  bool lrs_getnextbasis(D *&Dp) {
    long i = 0, j = 0;
    long m_ = Dp->m, d = Dp->d;
    bool backtrack = false;
    while ((j < d) || (Dp->B[m_] != m_)) {
      if (Dp->depth >= MAXD_DEPTH) backtrack = true;
      if (backtrack) {
        backtrack = false;
        if (!check_cache(Dp, &i, &j)) {
          Dp->depth--;
          selectpivot(Dp, &i, &j);
          pivot(Dp, i, j);
          update(Dp, &i, &j);
        }
        j++;
      }
      if (Dp->depth < -MAXD_DEPTH) break;
      while ((j < d) && (!reverse(Dp, &i, j))) j++;
      if (j == d)
        backtrack = true;
      else {
        cache_dict(Dp, i, j);
        Dp->depth++;
        pivot(Dp, i, j);
        update(Dp, &i, &j);
        Dp->lexflag = lexmin(Dp, 0);
        count[2]++;
        return true;
      }
    }
    return false;
  }
};

}  // namespace detail
}  // namespace LRSNash
