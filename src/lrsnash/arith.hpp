// Arithmetic policies for the header-only lrsnash port.
//
// Derived from lrslib 7.2 (lrslong.h/.c, lrsmp.h/.c, lrsgmp.h/.c),
// Copyright David Avis and contributors, GPL v2 or later (see COPYING).
//
// Every policy is a small object that owns its scratch space, so there is
// no global or thread_local state. A policy exposes:
//
//   using T;                       storage type of one integer
//   set(T&, long)                  a = v                      (itomp)
//   copy(T&, const T&)             a = b
//   zero/positive/negative/one     predicates
//   sign(const T&)                 -1 if a < 0 else +1  (lrs: zero is POS)
//   storesign(T&, long s)          a = |a| * s
//   changesign(T&)                 a = -a
//   greater(const T&, const T&)    a > b
//   mul(T& c, const T& a, const T& b)        c = a*b
//   qpiv(T& a, const T& b, const T& c, const T& d, const T& e)
//                                  a = (a*b - c*d) / e   (exact division)
//   comprod(a, b, c, d)            sign(a*b - c*d)
//   exactdiv(T& c, T& a, T& b)     c = a/b exactly (a may be clobbered)
//   gcd(T& u, T& v)                u = gcd(u, v)   (v may be clobbered)
//   lcm(T& a, const T& b)          a = lcm(a, b)
//   to_double(const T&)
//   pivot(T* A, stride, rows, cols, r, s, T& det, uint64_t& bound)
//                                  the whole lrs integer pivot on a dense
//                                  row-major block (see dense_pivot below).
//                                  `bound` is a per-dictionary magnitude hint
//                                  for the fixed-width kernels (~0: unknown);
//                                  the others ignore it.
//
// Fast exact division. In lrs's integer pivot every division
// (a*b - c*d)/det is exact (Bareiss / Sylvester identity: every dictionary
// entry is a minor of the input). Exact division by an odd e equals
// multiplication by e^-1 mod 2^w. With det = 2^k * e', a pivot computes
// e'^-1 once and then each entry costs two wrapping multiplies, a subtract
// and an arithmetic shift, with no hardware division and no libgcc __divti3.
// The wrapping result is the true quotient only if that quotient is known to
// fit. The I64/I128 kernels prove this per row with a cheap magnitude bound,
// and fall back to an exact double-width computation (I64: 128-bit, I128:
// 256-bit) for rows where the bound fails. So an overflow is still reported
// if and only if a true result does not fit. Stored values are kept in the
// symmetric range [-MAX, MAX] so that negation can never overflow.
//
// The fixed-width policies (I64, I128) check for overflow *exactly*: an
// operation throws only if its true result does not fit. Upstream lrs's SAFE
// mode instead throws when an operand exceeds a conservative bound, for
// example |x| > 2^31 for a 64-bit qpiv, which overflowed on roughly 13% of
// 9x9 games. Both give the same results whenever neither throws.
#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#if __has_include(<gmp.h>)
#include <gmp.h>
#define LRSNASH_HAVE_GMP 1
#endif

namespace LRSNash {

struct overflow_error : std::overflow_error {
  using std::overflow_error::overflow_error;
};

namespace arith {

[[noreturn]] inline void overflow(const char *what) {
  throw overflow_error(what);
}

// Inverse of an odd a modulo 2^64 / 2^128 (Newton iteration; (3a)^2 is
// correct to 5 bits and every step doubles the number of correct bits).
inline std::uint64_t inv64(std::uint64_t a) {
  std::uint64_t x = (3 * a) ^ 2;
  x *= 2 - a * x;
  x *= 2 - a * x;
  x *= 2 - a * x;
  x *= 2 - a * x;
  return x;
}
inline unsigned __int128 inv128(unsigned __int128 a) {
  unsigned __int128 x = inv64((std::uint64_t)a);
  x *= 2 - a * x;
  return x;
}
inline unsigned ctz128(unsigned __int128 a) {
  std::uint64_t lo = (std::uint64_t)a;
  return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((std::uint64_t)(a >> 64));
}
inline double u128_to_double(unsigned __int128 a) {
  return (double)(std::uint64_t)(a >> 64) * 0x1p64 + (double)(std::uint64_t)a;
}

// Reference full-tableau pivot (see I64::tableau_pivot). Requires Ars > 0.
template <class Ar>
void dense_tableau_pivot(Ar &ar, typename Ar::T *A, long stride, long rows,
                         long cols, long r, long s, typename Ar::T &det) {
  using T = typename Ar::T;
  T *Arow = A + r * stride;
  T Ars;
  ar.copy(Ars, Arow[s]);
  for (long i = 0; i < rows; i++) {
    if (i == r) continue;
    T *Ai = A + i * stride;
    for (long j = 0; j < cols; j++)
      if (j != s) ar.qpiv(Ai[j], Ars, Ai[s], Arow[j], det);
    ar.set(Ai[s], 0);
  }
  ar.copy(det, Ars);
}

// Reference implementation of the lrs integer pivot (upstream pivot(),
// without the objnum/objden bookkeeping) on a dense row-major block.
//   rows 0..rows-1, columns 0..cols-1, row stride `stride`,
//   (r, s) the pivot position, det > 0 the current determinant.
// For i != r, j != s:  A[i][j] = (A[i][j]*Ars - A[i][s]*A[r][j]) / (sg*det)
// with sg = sign(Ars); then column s (Ars < 0) or row r (Ars > 0) is
// negated, A[r][s] = sg*det, det = |Ars|.
template <class Ar>
void dense_pivot(Ar &ar, typename Ar::T *A, long stride, long rows, long cols,
                 long r, long s, typename Ar::T &det) {
  using T = typename Ar::T;
  T *Arow = A + r * stride;
  T Ars;
  ar.copy(Ars, Arow[s]);
  ar.storesign(det, ar.sign(Ars));
  for (long i = 0; i < rows; i++) {
    if (i == r) continue;
    T *Ai = A + i * stride;
    for (long j = 0; j < cols; j++)
      if (j != s) ar.qpiv(Ai[j], Ars, Ai[s], Arow[j], det);
  }
  if (ar.sign(Ars) == 1) {
    for (long j = 0; j < cols; j++)
      if (!ar.zero(Arow[j])) ar.changesign(Arow[j]);
  } else {
    for (long i = 0; i < rows; i++)
      if (!ar.zero(A[i * stride + s])) ar.changesign(A[i * stride + s]);
  }
  ar.copy(Arow[s], det);
  ar.copy(det, Ars);
  ar.storesign(det, 1);
}

// ---------------------------------------------------------------------------
// 64-bit integers, exact overflow detection via 128-bit intermediates.
// ---------------------------------------------------------------------------
struct I64 {
  using T = std::int64_t;
  static constexpr const char *name = "int64";

  static void set(T &a, long v) { a = v; }
  static void copy(T &a, const T &b) { a = b; }
  static bool zero(const T &a) { return a == 0; }
  static bool positive(const T &a) { return a > 0; }
  static bool negative(const T &a) { return a < 0; }
  static bool one(const T &a) { return a == 1; }
  static long sign(const T &a) { return a < 0 ? -1 : 1; }
  static void changesign(T &a) {
    if (a == std::numeric_limits<T>::min()) overflow("I64 changesign");
    a = -a;
  }
  static void storesign(T &a, long s) {
    if ((s < 0) != (a < 0) && a != 0) changesign(a);
  }
  static bool greater(const T &a, const T &b) { return a > b; }
  static void mul(T &c, const T &a, const T &b) {
    if (__builtin_mul_overflow(a, b, &c)) overflow("I64 mul");
  }
  static void qpiv(T &a, const T &b, const T &c, const T &d, const T &e) {
    __int128 t = (__int128)a * b - (__int128)c * d;  // cannot overflow
    t /= e;
    if (t > std::numeric_limits<T>::max() || t < std::numeric_limits<T>::min())
      overflow("I64 qpiv");
    a = (T)t;
  }
  static long comprod(const T &a, const T &b, const T &c, const T &d) {
    __int128 x = (__int128)a * b, y = (__int128)c * d;
    return x > y ? 1 : (x < y ? -1 : 0);
  }
  static void exactdiv(T &c, T &a, T &b) { c = a / b; }
  static void gcd(T &u, T &v) {
    std::uint64_t ul = u < 0 ? 0 - (std::uint64_t)u : (std::uint64_t)u;
    std::uint64_t vl = v < 0 ? 0 - (std::uint64_t)v : (std::uint64_t)v;
    if (ul == 0) {
      ul = vl;
    } else {
      while (vl != 0) {
        std::uint64_t r = ul % vl;
        ul = vl;
        vl = r;
      }
    }
    if (ul > (std::uint64_t)std::numeric_limits<T>::max()) overflow("I64 gcd");
    u = (T)ul;
  }
  void lcm(T &a, const T &b) {
    T u = a, v = b, q;
    gcd(u, v);
    exactdiv(q, a, u);
    mul(a, q, b);
  }
  static double to_double(const T &a) { return (double)a; }

  // Dense lrs pivot, same result as dense_pivot(). With sg = sign(Ars),
  // e = det = 2^k e' and t = a*|Ars| - (sg*b)*c (= sg*(a*Ars - b*c)), the new
  // entry is q = t/e. Fast path: q*2^k == t * e'^-1 (mod 2^64), so
  // q = (int64)(a*R - (sg*b)*C_j) >> k with R = |Ars| e'^-1 and
  // C_j = c_j e'^-1 precomputed. That is exact whenever |q| < 2^(63-k).
  // `bound` is a magnitude bound kept with the dictionary: every entry
  // satisfies |a| <= bound + 1 (bound == ~0: unknown). If
  //   (bound+1) (|Ars| + |row r|max) / e <= 2^(62-k)
  // every row takes the fast path. Otherwise each row is checked with its own
  // maximum, and rows that fail use exact 128-bit arithmetic. The new bound is
  // accumulated in the main loop.
  static std::uint64_t mag(T a) { return (std::uint64_t)(a ^ (a >> 63)); }
  void pivot(T *A, long stride, long rows, long cols, long r, long s, T &det,
             std::uint64_t &bound) {
    pivot_core<true>(A, stride, rows, cols, r, s, det, bound);
  }
  // `bound` for a dense block (OR of mag over all entries)
  static std::uint64_t magnitude_bound(const T *A, long stride, long rows,
                                       long cols) {
    std::uint64_t b = 0;
    for (long i = 0; i < rows; i++)
      for (long j = 0; j < cols; j++) b |= mag(A[i * stride + j]);
    return b;
  }
  // true if a pivot at (r, s) provably cannot overflow (every row takes the
  // fast path). Conservative: false does not mean it would overflow.
  static bool pivot_safe(const T *A, long stride, long cols, long r, long s,
                         const T &det, std::uint64_t bound) {
    if (bound == ~(std::uint64_t)0) return false;
    const T *Arow = A + r * stride;
    const std::uint64_t aars = (std::uint64_t)(Arow[s] < 0 ? -Arow[s] : Arow[s]);
    const std::uint64_t e = (std::uint64_t)det;
    const unsigned k = __builtin_ctzll(e);
    std::uint64_t Mr = 0;
    for (long j = 0; j < cols; j++) Mr |= mag(Arow[j]);
    return ((double)bound + 1.0) * ((double)aars + (double)Mr + 1.0) <=
           (double)e * (double)(1ULL << (62 - k));
  }
  static void from_i64(T &a, std::int64_t v) { a = v; }
  // Full-tableau (Edmonds/Bareiss) pivot, used by Lemke-Howson: for i != r,
  // A[i][j] = (A[i][j]*Ars - A[i][s]*A[r][j]) / det for all j (so column s
  // becomes 0), row r unchanged, det = Ars. Requires Ars > 0.
  void tableau_pivot(T *A, long stride, long rows, long cols, long r, long s,
                     T &det, std::uint64_t &bound) {
    pivot_core<false>(A, stride, rows, cols, r, s, det, bound);
  }
  template <bool LRS>
  void pivot_core(T *A, long stride, long rows, long cols, long r, long s,
                  T &det, std::uint64_t &bound) {
    using u64 = std::uint64_t;
    using i128 = __int128;
    using u128 = unsigned __int128;
    T *Arow = A + r * stride;
    const T ars = Arow[s];
    const bool neg = ars < 0;
    const u64 aars = neg ? (u64)0 - (u64)ars : (u64)ars;
    const u64 e = (u64)det;
    const unsigned k = __builtin_ctzll(e);
    const u64 inv = inv64(e >> k);
    const u64 R = aars * inv;
    if ((long)C.size() < cols) C.resize(cols);
    u64 *Cj = C.data();
    u64 Mr = 0;
    for (long j = 0; j < cols; j++) {
      Cj[j] = (u64)Arow[j] * inv;
      Mr |= mag(Arow[j]);
    }
    const double lim = (double)e * (double)(1ULL << (62 - k));
    const double daars = (double)aars, dMr = (double)Mr + 1.0;
    const bool all_fast =
        bound != ~(u64)0 && ((double)bound + 1.0) * (daars + dMr) <= lim;
    u64 nb = Mr | e | aars;
    for (long i = 0; i < rows; i++) {
      if (i == r) continue;
      T *Ai = A + i * stride;
      const T b = Ai[s];
      const T bs = neg ? -b : b;
      bool fast = all_fast;
      if (!fast) {
        u64 Mi = 0;
        for (long j = 0; j < cols; j++) Mi |= mag(Ai[j]);
        fast = ((double)Mi + 1.0) * daars + (double)(bs < 0 ? -bs : bs) * dMr <=
               lim;
      }
      u64 acc = mag(bs);
      if (fast) {
        const u64 ubs = (u64)bs;
        for (long j = 0; j < cols; j++) {
          const T q = (T)((u64)Ai[j] * R - ubs * Cj[j]) >> k;
          Ai[j] = q;
          acc |= mag(q);
        }
      } else {
        const i128 D = (i128)e * (i128)std::numeric_limits<T>::max();
        const u128 W = (u128)D * 2;
        for (long j = 0; j < cols; j++) {
          i128 t = (i128)Ai[j] * (i128)aars - (i128)bs * (i128)Arow[j];
          if ((u128)t + (u128)D > W) overflow("I64 pivot");
          const T q = (T)((u64)(t >> k) * inv);
          Ai[j] = q;
          acc |= mag(q);
        }
      }
      if (LRS) Ai[s] = bs;
      nb |= acc;
    }
    if (LRS) {
      if (!neg)
        for (long j = 0; j < cols; j++) Arow[j] = -Arow[j];
      Arow[s] = neg ? -(T)e : (T)e;
    }
    det = (T)aars;
    bound = nb;
  }

 private:
  std::vector<std::uint64_t> C;  // pivot scratch
};

// ---------------------------------------------------------------------------
// 128-bit integers. Fast path uses builtin overflow checks; when an
// intermediate product overflows we fall back to exact 256-bit arithmetic so
// that only a result that truly does not fit throws.
// ---------------------------------------------------------------------------
struct I128 {
  using T = __int128;
  using U = unsigned __int128;
  static constexpr const char *name = "int128";
  static constexpr T MAX = (T)(~(U)0 >> 1);
  static constexpr T MIN = -MAX - 1;

  static void set(T &a, long v) { a = v; }
  static void copy(T &a, const T &b) { a = b; }
  static bool zero(const T &a) { return a == 0; }
  static bool positive(const T &a) { return a > 0; }
  static bool negative(const T &a) { return a < 0; }
  static bool one(const T &a) { return a == 1; }
  static long sign(const T &a) { return a < 0 ? -1 : 1; }
  static void changesign(T &a) {
    if (a == MIN) overflow("I128 changesign");
    a = -a;
  }
  static void storesign(T &a, long s) {
    if ((s < 0) != (a < 0) && a != 0) changesign(a);
  }
  static bool greater(const T &a, const T &b) { return a > b; }
  static void mul(T &c, const T &a, const T &b) {
    if (__builtin_mul_overflow(a, b, &c)) overflow("I128 mul");
  }

  // --- exact 256-bit helpers (magnitudes, hi:lo) ---
  struct W256 {
    U hi, lo;
  };
  static U uabs(T a) { return a < 0 ? (U)0 - (U)a : (U)a; }
  static W256 umul(U a, U b) {
    const U M = ~(std::uint64_t)0;
    U a0 = a & M, a1 = a >> 64, b0 = b & M, b1 = b >> 64;
    U p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    U mid = (p00 >> 64) + (p01 & M) + (p10 & M);
    W256 r;
    r.lo = (p00 & M) | (mid << 64);
    r.hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
    return r;
  }
  static int ucmp(const W256 &x, const W256 &y) {
    if (x.hi != y.hi) return x.hi > y.hi ? 1 : -1;
    if (x.lo != y.lo) return x.lo > y.lo ? 1 : -1;
    return 0;
  }
  static W256 uadd(const W256 &x, const W256 &y) {
    W256 r;
    r.lo = x.lo + y.lo;
    r.hi = x.hi + y.hi + (r.lo < x.lo);
    return r;
  }
  static W256 usub(const W256 &x, const W256 &y) {  // requires x >= y
    W256 r;
    r.lo = x.lo - y.lo;
    r.hi = x.hi - y.hi - (x.lo < y.lo);
    return r;
  }

  static long comprod(const T &a, const T &b, const T &c, const T &d) {
    T x, y;
    if (!__builtin_mul_overflow(a, b, &x) && !__builtin_mul_overflow(c, d, &y))
      return x > y ? 1 : (x < y ? -1 : 0);
    // signs: s = -1, 0, +1
    int sx = (a == 0 || b == 0) ? 0 : (((a < 0) != (b < 0)) ? -1 : 1);
    int sy = (c == 0 || d == 0) ? 0 : (((c < 0) != (d < 0)) ? -1 : 1);
    if (sx != sy) return sx > sy ? 1 : -1;
    if (sx == 0) return 0;
    int m = ucmp(umul(uabs(a), uabs(b)), umul(uabs(c), uabs(d)));
    return sx > 0 ? m : -m;
  }

  static void qpiv(T &a, const T &b, const T &c, const T &d, const T &e) {
    T x, y, z;
    if (!__builtin_mul_overflow(a, b, &x) && !__builtin_mul_overflow(c, d, &y) &&
        !__builtin_sub_overflow(x, y, &z)) {
      a = z / e;
      return;
    }
    qpiv_slow(a, b, c, d, e);
  }
  static void qpiv_slow(T &a, const T &b, const T &c, const T &d, const T &e) {
    bool nx = (a != 0 && b != 0) && ((a < 0) != (b < 0));
    bool ny = (c != 0 && d != 0) && ((c < 0) != (d < 0));
    W256 X = umul(uabs(a), uabs(b)), Y = umul(uabs(c), uabs(d));
    // z = X*(nx?-1:1) - Y*(ny?-1:1)
    W256 Z;
    bool nz;
    if (nx != ny) {  // same effective sign: |z| = X + Y, sign of x
      Z = uadd(X, Y);
      nz = nx;
    } else if (ucmp(X, Y) >= 0) {
      Z = usub(X, Y);
      nz = nx;
    } else {
      Z = usub(Y, X);
      nz = !nx;
    }
    U ue = uabs(e);
    if (Z.hi >= ue) overflow("I128 qpiv");  // quotient >= 2^128
    U rem = Z.hi, q = 0;
    for (int i = 127; i >= 0; --i) {
      bool top = rem >> 127;  // rem < ue <= 2^127, so top is always 0
      (void)top;
      rem = (rem << 1) | ((Z.lo >> i) & 1);
      q <<= 1;
      if (rem >= ue) {
        rem -= ue;
        q |= 1;
      }
    }
    bool neg = nz != (e < 0);
    if (q > (U)MAX) overflow("I128 qpiv");
    a = neg ? -(T)q : (T)q;
  }
  static void exactdiv(T &c, T &a, T &b) { c = a / b; }
  static void gcd(T &u, T &v) {
    U ul = uabs(u), vl = uabs(v);
    if (ul == 0) {
      ul = vl;
    } else {
      while (vl != 0) {
        U r = ul % vl;
        ul = vl;
        vl = r;
      }
    }
    if (ul > (U)MAX) overflow("I128 gcd");
    u = (T)ul;
  }
  void lcm(T &a, const T &b) {
    T u = a, v = b, q;
    gcd(u, v);
    exactdiv(q, a, u);
    mul(a, q, b);
  }
  static double to_double(const T &a) { return (double)a; }

  // Dense lrs pivot, same result as dense_pivot(); see I64::pivot. Here the
  // wrapping arithmetic is mod 2^128, the fast-path condition is
  // |q| < 2^(127-k), and the fallback is exact 256-bit arithmetic. Bounds use
  // only the high 64 bits of each entry: |a| <= (hmag(a) + 1) * 2^64, and
  // `bound` bounds hmag over the dictionary (~0: unknown).
  static std::uint64_t hmag(T a) {
    const std::int64_t h = (std::int64_t)(a >> 64);
    return (std::uint64_t)(h ^ (h >> 63));
  }
  static std::uint64_t magnitude_bound(const T *A, long stride, long rows,
                                       long cols) {
    std::uint64_t b = 0;
    for (long i = 0; i < rows; i++)
      for (long j = 0; j < cols; j++) b |= hmag(A[i * stride + j]);
    return b;
  }
  static bool pivot_safe(const T *A, long stride, long cols, long r, long s,
                         const T &det, std::uint64_t bound) {
    if (bound == ~(std::uint64_t)0) return false;
    const T *Arow = A + r * stride;
    const U aars = uabs(Arow[s]);
    const U e = (U)det;
    const unsigned k = ctz128(e);
    std::uint64_t Mr = 0;
    for (long j = 0; j < cols; j++) Mr |= hmag(Arow[j]);
    return ((double)bound + 1.0) * 0x1p64 *
               (u128_to_double(aars) + ((double)Mr + 1.0) * 0x1p64) <=
           u128_to_double(e) * __builtin_ldexp(1.0, 126 - (int)k);
  }
  static void from_i64(T &a, std::int64_t v) { a = v; }
  static void from_i128(T &a, __int128 v) { a = v; }
  void pivot(T *A, long stride, long rows, long cols, long r, long s, T &det,
             std::uint64_t &bound) {
    pivot_core<true>(A, stride, rows, cols, r, s, det, bound);
  }
  // Full-tableau (Edmonds/Bareiss) pivot, used by Lemke-Howson: for i != r,
  // A[i][j] = (A[i][j]*Ars - A[i][s]*A[r][j]) / det for all j (so column s
  // becomes 0), row r unchanged, det = Ars. Requires Ars > 0.
  void tableau_pivot(T *A, long stride, long rows, long cols, long r, long s,
                     T &det, std::uint64_t &bound) {
    pivot_core<false>(A, stride, rows, cols, r, s, det, bound);
  }
  template <bool LRS>
  void pivot_core(T *A, long stride, long rows, long cols, long r, long s,
                  T &det, std::uint64_t &bound) {
    using u64 = std::uint64_t;
    T *Arow = A + r * stride;
    const T ars = Arow[s];
    const bool neg = ars < 0;
    const U aars = uabs(ars);
    const U e = (U)det;
    const unsigned k = ctz128(e);
    const U inv = inv128(e >> k);
    const U R = aars * inv;
    if ((long)C.size() < cols) C.resize(cols);
    U *Cj = C.data();
    u64 Mr = 0;
    for (long j = 0; j < cols; j++) {
      Cj[j] = (U)Arow[j] * inv;
      Mr |= hmag(Arow[j]);
    }
    const double lim = u128_to_double(e) * __builtin_ldexp(1.0, 126 - (int)k);
    const double daars = u128_to_double(aars), dMr = ((double)Mr + 1.0) * 0x1p64;
    const bool all_fast = bound != ~(u64)0 &&
                          ((double)bound + 1.0) * 0x1p64 * (daars + dMr) <= lim;
    u64 nb = Mr | (u64)(e >> 64) | (u64)(aars >> 64);
    for (long i = 0; i < rows; i++) {
      if (i == r) continue;
      T *Ai = A + i * stride;
      const T b = Ai[s];
      const T bs = neg ? -b : b;
      bool fast = all_fast;
      if (!fast) {
        u64 Mi = 0;
        for (long j = 0; j < cols; j++) Mi |= hmag(Ai[j]);
        fast = ((double)Mi + 1.0) * 0x1p64 * daars +
                   u128_to_double(uabs(bs)) * dMr <=
               lim;
      }
      u64 acc = hmag(bs);
      if (fast) {
        const U ubs = (U)bs;
        for (long j = 0; j < cols; j++) {
          const T q = (T)((U)Ai[j] * R - ubs * Cj[j]) >> k;
          Ai[j] = q;
          acc |= hmag(q);
        }
      } else {
        pivot_row_exact(Ai, Arow, cols, aars, bs, e, k, inv);
        for (long j = 0; j < cols; j++) acc |= hmag(Ai[j]);
      }
      if (LRS) Ai[s] = bs;
      nb |= acc;
    }
    if (LRS) {
      if (!neg)
        for (long j = 0; j < cols; j++) Arow[j] = -Arow[j];
      Arow[s] = neg ? -(T)e : (T)e;
    }
    det = (T)aars;
    bound = nb;
  }

  // Ai[j] = (Ai[j]*aars - bs*Arow[j]) / e exactly, in 256-bit arithmetic.
  static void pivot_row_exact(T *Ai, const T *Arow, long cols, U aars, T bs,
                              U e, unsigned k, U inv) {
    const W256 E = umul(e, (U)MAX);  // largest allowed |t|
    const bool nbs = bs < 0;
    const U abs_bs = uabs(bs);
    for (long j = 0; j < cols; j++) {
      const T a = Ai[j], c = Arow[j];
      const W256 P1 = umul(uabs(a), aars), P2 = umul(abs_bs, uabs(c));
      const bool n1 = a < 0, n2 = nbs != (c < 0);
      W256 Z;
      bool nz;
      if (n1 != n2) {
        Z = uadd(P1, P2);
        nz = n1;
      } else if (ucmp(P1, P2) >= 0) {
        Z = usub(P1, P2);
        nz = n1;
      } else {
        Z = usub(P2, P1);
        nz = !n1;
      }
      if (ucmp(Z, E) > 0) overflow("I128 pivot");
      const U zk = k == 0 ? Z.lo : (Z.lo >> k) | (Z.hi << (128 - k));
      const U q = zk * inv;  // |t| / e, exact and <= MAX
      Ai[j] = nz ? -(T)q : (T)q;
    }
  }

 private:
  std::vector<U> C;  // pivot scratch
};

// ---------------------------------------------------------------------------
// lrsmp: lrs's own multiple precision, base 10^9 digits in long long words.
// Layout of a number a: a[0] = +-length (length counts a[0] itself, so zero
// has length 2 with a[1] == 0); sign of a[0] is the sign of the number.
// W is the maximum number of digit words (upstream lrs_digits; the default
// DEC2DIG(1000) == 112 matches upstream lrs_mp_init(0)).
// ---------------------------------------------------------------------------
template <long W = 112> struct MP {
  static constexpr long DIGITS = W;
  static constexpr long long BASE = 1000000000LL;
  static constexpr long long MAXD = 9223372036854775807LL;
  static constexpr const char *name = "lrsmp";
  // Two words of headroom: normalize/divint may touch one word past length.
  struct T {
    long long v[W + 3];
    T() {
      v[0] = 2;
      v[1] = 0;
    }
    long long &operator[](long i) { return v[i]; }
    const long long &operator[](long i) const { return v[i]; }
  };

  static long length(const T &a) { return a[0] > 0 ? a[0] : -a[0]; }
  static void storelength(T &a, long la) { a[0] = a[0] > 0 ? la : -la; }
  static long sign(const T &a) { return a[0] < 0 ? -1 : 1; }
  static void storesign(T &a, long sa) { a[0] = a[0] > 0 ? sa * a[0] : -sa * a[0]; }
  static void changesign(T &a) { a[0] = -a[0]; }
  static bool positive(const T &a) { return !(a[0] < 2 || (a[0] == 2 && a[1] == 0)); }
  static bool negative(const T &a) { return !(a[0] > -2 || (a[0] == -2 && a[1] == 0)); }
  static bool zero(const T &a) { return (a[0] == 2 || a[0] == -2) && a[1] == 0; }
  static bool one(const T &a) { return a[0] == 2 && a[1] == 1; }

  static void copy(T &a, const T &b) {
    long lb = length(b);
    for (long i = 0; i <= lb; i++) a[i] = b[i];
  }

  static void set(T &a, long in) {
    a[0] = 2;
    for (long i = 1; i < W; i++) a[i] = 0;
    long long x = in;
    if (x < 0) {
      storesign(a, -1);
      x = -x;
    }
    long i = 0;
    while (x != 0) {
      i++;
      a[i] = x - BASE * (x / BASE);
      x = x / BASE;
      storelength(a, i + 1);
    }
  }

  static bool greater(const T &a, const T &b) {
    if (a[0] > b[0]) return true;
    if (a[0] < b[0]) return false;
    for (long i = length(a) - 1; i >= 1; i--) {
      if (a[i] < b[i]) return sign(a) != 1;
      if (a[i] > b[i]) return sign(a) != -1;
    }
    return false;
  }

  static void normalize(T &a) {
    long long cy, i, la;
    la = length(a);
  start:
    cy = 0;
    for (i = 1; i < la; i++) {
      cy = (a[i] += cy) / BASE;
      a[i] -= cy * BASE;
      if (a[i] < 0) {
        a[i] += BASE;
        cy--;
      }
    }
    while (cy > 0) {
      if (i > W + 1) overflow("lrsmp normalize");
      a[i++] = cy % BASE;
      cy /= BASE;
    }
    if (cy < 0) {
      a[la - 1] += cy * BASE;
      for (i = 1; i < la; i++) a[i] = (-a[i]);
      storesign(a, sign(a) == 1 ? -1 : 1);
      goto start;
    }
    while (a[i - 1] == 0 && i > 2) i--;
    if (i > W) overflow("lrsmp digits");
    storelength(a, i);
    if (i == 2 && a[1] == 0) storesign(a, 1);
  }

  // c = a*b ; b and c may coincide
  static void mulint(const T &a, const T &b, T &c) {
    long nlength, i, j, la, lb;
    la = length(a);
    lb = length(b);
    nlength = la + lb - 2;
    if (nlength > W) overflow("lrsmp mulint");
    for (i = 0; i < la - 2; i++) c[lb + i] = 0;
    for (i = lb - 1; i > 0; i--) {
      for (j = 2; j < la; j++)
        if ((c[i + j - 1] += b[i] * a[j]) >
            MAXD - (BASE - 1) * (BASE - 1) - MAXD / BASE) {
          c[i + j - 1] -= (MAXD / BASE) * BASE;
          c[i + j] += MAXD / BASE;
        }
      c[i] = b[i] * a[1];
    }
    long sa = sign(a), sb = sign(b);
    storelength(c, nlength);
    storesign(c, sa == sb ? 1 : -1);
    normalize(c);
  }
  void mul(T &c, const T &a, const T &b) {
    if (&c == &a) {  // upstream mulint allows only b == c aliasing
      t1 = a;
      mulint(t1, b, c);
    } else {
      mulint(a, b, c);
    }
  }

  // a = a*ka + b*kb
  static void linint(T &a, long ka, const T &b, long kb) {
    long i, la, lb;
    la = length(a);
    lb = length(b);
    for (i = 1; i < la; i++) a[i] *= ka;
    if (sign(a) != sign(b)) kb = (-kb);
    if (lb > la) {
      storelength(a, lb);
      for (i = la; i < lb; i++) a[i] = 0;
    }
    for (i = 1; i < lb; i++) a[i] += kb * b[i];
    normalize(a);
  }

  // c = a/b; a contains remainder on return; b is restored. (Knuth D)
  static void divint(T &a, T &b, T &c) {
    long long cy, la, lb, lc, d1, s, t, sig;
    long long i, j, qh;
    sig = sign(a) * sign(b);
    la = length(a);
    lb = length(b);
    lc = la - lb + 2;
    if (la < lb) {
      c[0] = 2;  // storelength(c,TWO); storesign(c,POS)
      c[1] = 0;
      normalize(c);
      return;
    }
    for (i = 1; i < lc; i++) c[i] = 0;
    storelength(c, lc);
    storesign(c, (sign(a) == sign(b)) ? 1 : -1);
    if (lb == 2) {
      cy = 0;
      t = b[1];
      for (i = la - 1; i > 0; i--) {
        cy = cy * BASE + a[i];
        a[i] = 0;
        cy -= (c[i] = cy / t) * t;
      }
      a[1] = cy;
      storesign(a, (cy == 0) ? 1 : sign(a));
      storelength(a, 2);
      storesign(c, sig);
      normalize(c);
      return;
    }
    d1 = BASE / (b[lb - 1] + 1);
    if (d1 > 1) {
      cy = 0;
      for (i = 1; i < la; i++) {
        cy = (a[i] = a[i] * d1 + cy) / BASE;
        a[i] %= BASE;
      }
      a[i] = cy;
      cy = 0;
      for (i = 1; i < lb; i++) {
        cy = (b[i] = b[i] * d1 + cy) / BASE;
        b[i] %= BASE;
      }
      b[i] = cy;
    } else {
      a[la] = 0;
      b[lb] = 0;
    }
    for (j = 0; j <= la - lb; j++) {
      if (a[la - j] == b[lb - 1])
        qh = BASE - 1;
      else {
        s = (a[la - j] * BASE + a[la - j - 1]);
        qh = s / b[lb - 1];
        while (qh * b[lb - 2] > (s - qh * b[lb - 1]) * BASE + a[la - j - 2]) qh--;
      }
      cy = 0;
      for (i = 1; i <= lb; i++) {
        s = qh * b[i] + cy;
        a[la - j - lb + i] -= s % BASE;
        cy = s / BASE;
        if (a[la - j - lb + i] < 0) {
          a[la - j - lb + i] += BASE;
          cy++;
        }
      }
      if (cy) {
        qh--;
        cy = 0;
        for (i = 1; i <= lb; i++) {
          a[la - j - lb + i] += b[i] + cy;
          cy = a[la - j - lb + i] / BASE;
          a[la - j - lb + i] %= BASE;
        }
      }
      c[la - lb - j + 1] = qh;
    }
    for (i = lc; c[i - 1] == 0 && i > 2; i--)
      ;
    storelength(c, i);
    if (i == 2 && c[1] == 0) storesign(c, 1);
    cy = 0;
    for (i = lb - 1; i >= 1; i--) {
      cy = (a[i] += cy * BASE) % d1;
      a[i] /= d1;
    }
    for (i = la; a[i - 1] == 0 && i > 2; i--)
      ;
    storelength(a, i);
    if (i == 2 && a[1] == 0) storesign(a, 1);
    if (cy) throw std::logic_error("lrsmp divide error");
    for (i = lb - 1; i >= 1; i--) {
      cy = (b[i] += cy * BASE) % d1;
      b[i] /= d1;
    }
  }

  static void exactdiv(T &c, T &a, T &b) { divint(a, b, c); }

  // u = gcd(u, v), v destroyed
  void gcd(T &u, T &v) {
    T &r = t2;
    unsigned long long ul, vl;
    long i;
    constexpr unsigned long long maxspval = 9;  // MAXD / BASE^2
    constexpr long maxsplen = 4;
    if (greater(v, u)) goto bigv;
  bigu:
    if (zero(v)) return;
    if ((i = length(u)) < maxsplen || (i == maxsplen && (unsigned long long)u[maxsplen - 1] < maxspval))
      goto quickfinish;
    divint(u, v, r);
    normalize(u);
  bigv:
    if (zero(u)) {
      copy(u, v);
      return;
    }
    if ((i = length(v)) < maxsplen || (i == maxsplen && (unsigned long long)v[maxsplen - 1] < maxspval))
      goto quickfinish;
    divint(v, u, r);
    normalize(v);
    goto bigu;
  quickfinish:
    ul = vl = 0;
    for (i = length(u) - 1; i > 0; i--) ul = BASE * ul + u[i];
    for (i = length(v) - 1; i > 0; i--) vl = BASE * vl + v[i];
    if (ul > vl) goto qv;
  qu:
    if (!vl) {
      for (i = 1; ul; i++) {
        u[i] = ul % BASE;
        ul = ul / BASE;
      }
      storelength(u, i);
      return;
    }
    ul %= vl;
  qv:
    if (!ul) {
      for (i = 1; vl; i++) {
        u[i] = vl % BASE;
        vl = vl / BASE;
      }
      storelength(u, i);
      return;
    }
    vl %= ul;
    goto qu;
  }

  void lcm(T &a, const T &b) {
    copy(t3, a);
    copy(t4, b);
    gcd(t3, t4);
    exactdiv(t4, a, t3);  // t4 = a/u
    mulint(t4, b, a);
  }

  long comprod(const T &a, const T &b, const T &c, const T &d) {
    mulint(a, b, t1);
    mulint(c, d, t2);
    linint(t1, 1, t2, -1);
    if (positive(t1)) return 1;
    if (negative(t1)) return -1;
    return 0;
  }

  // a = (a*b - c*d)/e, exactly as upstream pivot() for non-LRSLONG builds
  void qpiv(T &a, const T &b, const T &c, const T &d, T &e) {
    mulint(a, b, t1);
    mulint(c, d, t2);
    linint(t1, 1, t2, -1);  // decint
    divint(t1, e, a);
  }

  void pivot(T *A, long stride, long rows, long cols, long r, long s, T &det,
             std::uint64_t &) {
    dense_pivot(*this, A, stride, rows, cols, r, s, det);
  }
  void tableau_pivot(T *A, long stride, long rows, long cols, long r, long s,
                     T &det, std::uint64_t &) {
    dense_tableau_pivot(*this, A, stride, rows, cols, r, s, det);
  }
  static std::uint64_t magnitude_bound(const T *, long, long, long) {
    return ~(std::uint64_t)0;
  }
  static bool pivot_safe(const T *, long, long, long, long, const T &,
                         std::uint64_t) {
    return true;  // overflow is checked exactly inside the arithmetic
  }
  static void from_i64(T &a, std::int64_t v) { from_i128(a, v); }
  static void from_i128(T &a, __int128 v) {
    unsigned __int128 u = v < 0 ? (unsigned __int128)0 - (unsigned __int128)v
                                : (unsigned __int128)v;
    long i = 1;
    a[1] = 0;
    while (u != 0) {
      a[i++] = (long long)(u % (unsigned __int128)BASE);
      u /= (unsigned __int128)BASE;
    }
    a[0] = i < 2 ? 2 : i;
    if (v < 0) a[0] = -a[0];
  }

  static double to_double(const T &a) {
    double x = 0, y = 1.0;
    long la = length(a);
    for (long i = 1; i < la; i++) {
      x = x + y * a[i];
      y = y * BASE;
    }
    return negative(a) ? -x : x;
  }

 private:
  T t1, t2, t3, t4;  // scratch
};

#ifdef LRSNASH_HAVE_GMP
// ---------------------------------------------------------------------------
// GMP
// ---------------------------------------------------------------------------
struct GMP {
  static constexpr const char *name = "gmp";
  struct T {
    mpz_t z;
    T() { mpz_init(z); }
    T(const T &o) { mpz_init_set(z, o.z); }
    T &operator=(const T &o) {
      mpz_set(z, o.z);
      return *this;
    }
    ~T() { mpz_clear(z); }
  };
  static void set(T &a, long v) { mpz_set_si(a.z, v); }
  static void copy(T &a, const T &b) { mpz_set(a.z, b.z); }
  static bool zero(const T &a) { return mpz_sgn(a.z) == 0; }
  static bool positive(const T &a) { return mpz_sgn(a.z) > 0; }
  static bool negative(const T &a) { return mpz_sgn(a.z) < 0; }
  static bool one(const T &a) { return mpz_cmp_si(a.z, 1) == 0; }
  static long sign(const T &a) { return mpz_sgn(a.z) < 0 ? -1 : 1; }
  static void changesign(T &a) { mpz_neg(a.z, a.z); }
  static void storesign(T &a, long s) {
    if (s * sign(a) < 0) mpz_neg(a.z, a.z);
  }
  static bool greater(const T &a, const T &b) { return mpz_cmp(a.z, b.z) > 0; }
  static void mul(T &c, const T &a, const T &b) { mpz_mul(c.z, a.z, b.z); }
  void qpiv(T &a, const T &b, const T &c, const T &d, const T &e) {
    mpz_mul(t1.z, a.z, b.z);
    mpz_submul(t1.z, c.z, d.z);
    mpz_divexact(a.z, t1.z, e.z);
  }
  long comprod(const T &a, const T &b, const T &c, const T &d) {
    mpz_mul(t1.z, a.z, b.z);
    mpz_mul(t2.z, c.z, d.z);
    int i = mpz_cmp(t1.z, t2.z);
    return i > 0 ? 1 : (i < 0 ? -1 : 0);
  }
  static void exactdiv(T &c, T &a, T &b) { mpz_divexact(c.z, a.z, b.z); }
  static void gcd(T &u, T &v) { mpz_gcd(u.z, u.z, v.z); }
  void lcm(T &a, const T &b) {
    mpz_set(t1.z, a.z);
    mpz_gcd(t1.z, t1.z, b.z);
    mpz_divexact(t2.z, a.z, t1.z);
    mpz_mul(a.z, t2.z, b.z);
  }
  static double to_double(const T &a) { return mpz_get_d(a.z); }

  void pivot(T *A, long stride, long rows, long cols, long r, long s, T &det,
             std::uint64_t &) {
    dense_pivot(*this, A, stride, rows, cols, r, s, det);
  }
  void tableau_pivot(T *A, long stride, long rows, long cols, long r, long s,
                     T &det, std::uint64_t &) {
    dense_tableau_pivot(*this, A, stride, rows, cols, r, s, det);
  }
  static std::uint64_t magnitude_bound(const T *, long, long, long) {
    return ~(std::uint64_t)0;
  }
  static bool pivot_safe(const T *, long, long, long, long, const T &,
                         std::uint64_t) {
    return true;
  }
  static void from_i64(T &a, std::int64_t v) { mpz_set_si(a.z, v); }
  static void from_i128(T &a, __int128 v) {
    unsigned __int128 u = v < 0 ? (unsigned __int128)0 - (unsigned __int128)v
                                : (unsigned __int128)v;
    mpz_set_ui(a.z, (unsigned long)(u >> 64));
    mpz_mul_2exp(a.z, a.z, 64);
    mpz_add_ui(a.z, a.z, (unsigned long)(std::uint64_t)u);
    if (v < 0) mpz_neg(a.z, a.z);
  }

 private:
  T t1, t2;
};
#endif

}  // namespace arith
}  // namespace LRSNash
