// Contract test for the fixed-width pivot kernels in src/lrsnash/arith.hpp
// (I64::pivot, I128::pivot), against exact GMP arithmetic.
//
// A random integer matrix (det = 1) is pivoted along a random sequence of
// nonzero pivot elements. Every step is done exactly in GMP with the
// reference kernel (dense_pivot), so every division
// is exact (Bareiss). Before each step the fixed-width kernel gets the same
// input, and the test checks that:
//  * it throws overflow_error if and only if some true result entry (or the
//    new det) does not fit in T, as for the exact reference;
//  * otherwise every entry and det are bit-identical to GMP;
//  * the returned bound is sound: unless it is ~0 (unknown), every entry a
//    satisfies |a| <= bound + 1 (I64) or |a| <= (bound + 1) * 2^64 (I128),
//    which is what the fast-path test of the next pivot relies on.
// The fast path is exercised with an unknown bound, an exact bound and an
// overestimated bound. Inputs are shaped to reach every path: small and
// large magnitudes, pivots with det = 2^k e' for small and large k
// (k >= 64 for I128), rows on both sides of the fast-path limit, and
// results just inside and just outside the range of T. I64 runs once with
// its plain kernel and, where the CPU has AVX2, once with the AVX2 kernel.
//   g++ -O2 -std=c++17 -Isrc tools/test_kernel.cc -lgmp -o test_kernel
//   ./test_kernel trials seed
#include "lrsnash/arith.hpp"
#include <gmp.h>
#include <cstdio>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

using namespace LRSNash;
using namespace LRSNash::arith;

static void to_mpz(mpz_t z, std::int64_t a) { mpz_set_si(z, a); }
static void to_mpz(mpz_t z, __int128 a) {
  unsigned __int128 u = a < 0 ? -(unsigned __int128)a : (unsigned __int128)a;
  mpz_set_ui(z, (std::uint64_t)(u >> 64));
  mpz_mul_2exp(z, z, 64);
  mpz_add_ui(z, z, (std::uint64_t)u);
  if (a < 0) mpz_neg(z, z);
}
// z fits in T (symmetric range [-MAX, MAX], as the kernels keep it)
template <class T> static bool fits(const mpz_t z);
template <> bool fits<std::int64_t>(const mpz_t z) {
  return mpz_sizeinbase(z, 2) <= 63;  // |z| < 2^63
}
template <> bool fits<__int128>(const mpz_t z) { return mpz_sizeinbase(z, 2) <= 127; }
template <class T> static T from_mpz(const mpz_t z);
template <> std::int64_t from_mpz<std::int64_t>(const mpz_t z) { return mpz_get_si(z); }
template <> __int128 from_mpz<__int128>(const mpz_t z) {
  mpz_t a, lo;
  mpz_inits(a, lo, nullptr);
  mpz_abs(a, z);
  mpz_fdiv_r_2exp(lo, a, 64);
  unsigned __int128 u = mpz_get_ui(lo);
  mpz_fdiv_q_2exp(a, a, 64);
  u |= (unsigned __int128)mpz_get_ui(a) << 64;
  mpz_clears(a, lo, nullptr);
  return mpz_sgn(z) < 0 ? -(__int128)u : (__int128)u;
}

// I64's kernel choice for this run (see I64::avx2)
static bool g_use_avx2 = false;
template <class Ar> static void choose_kernel(Ar &ar) {
#ifdef LRSNASH_AVX2
  if constexpr (std::is_same_v<Ar, I64>) ar.avx2 = g_use_avx2;
#endif
  (void)ar;
}

struct Stats {
  long steps = 0, ovf = 0, bad = 0, bad_ovf = 0, bad_bound = 0, k0 = 0, ksmall = 0, kbig = 0;
  long fast_known = 0;
};

// One pivot sequence; pivot elements of both signs.
template <class Ar>
static void sequence(std::mt19937_64 &g, Stats &st, long rows, long cols, int mode) {
  using T = typename Ar::T;
  GMP gm;
  Ar ar;
  choose_kernel(ar);
  const long stride = cols + (long)(g() % 3);  // stride > cols sometimes
  std::vector<GMP::T> M(rows * stride);
  GMP::T det;
  gm.set(det, 1);
  // initial entries: mode 0 small, 1 medium, 2 near the limit of T,
  // 3 multiples of a large power of two (large k later)
  const int tbits = sizeof(T) == 8 ? 63 : 127;
  for (long i = 0; i < rows; i++)
    for (long j = 0; j < cols; j++) {
      mpz_t &z = M[i * stride + j].z;
      int bits = mode == 0 ? 3 : mode == 1 ? 1 + (int)(g() % 20) : mode == 2 ? (int)(g() % (tbits / 2 + 4)) : 1 + (int)(g() % 6);
      mpz_set_ui(z, g() >> (64 - std::max(1, std::min(bits, 63))));
      if (bits > 63) mpz_mul_2exp(z, z, bits - 63);
      if (mode == 3 && g() % 2) mpz_mul_2exp(z, z, sizeof(T) == 8 ? 20 + g() % 25 : 40 + g() % 50);
      if (g() % 2) mpz_neg(z, z);
      if (g() % 5 == 0) mpz_set_ui(z, 0);
    }
  mpz_t z;
  mpz_init(z);
  for (int step = 0; step < 12; step++) {
    // pick a nonzero pivot element
    long r = -1, s = -1;
    for (int tries = 0; tries < 50 && r < 0; tries++) {
      long rr = (long)(g() % rows), ss = (long)(g() % cols);
      if (mpz_sgn(M[rr * stride + ss].z) != 0) r = rr, s = ss;
    }
    if (r < 0) break;
    // the fixed-width input, if the current matrix fits
    bool in_fits = fits<T>(det.z);
    for (long i = 0; i < rows && in_fits; i++)
      for (long j = 0; j < cols && in_fits; j++) in_fits = fits<T>(M[i * stride + j].z);
    if (!in_fits) break;
    std::vector<T> F(rows * stride, (T)0x5a5a);  // padding sentinels
    for (long i = 0; i < rows; i++)
      for (long j = 0; j < cols; j++) F[i * stride + j] = from_mpz<T>(M[i * stride + j].z);
    T fdet = from_mpz<T>(det.z);
    std::uint64_t bound;
    const int bmode = (int)(g() % 3);
    if (bmode == 0) bound = ~(std::uint64_t)0;
    else bound = Ar::magnitude_bound(F.data(), stride, rows, cols);
    if (bmode == 2 && bound < (1ULL << 40)) bound = bound * 3 + (g() % 1000);
    if (bound != ~(std::uint64_t)0) ++st.fast_known;
    {  // k of this pivot
      const unsigned k = mpz_scan1(det.z, 0);
      if (k == 0) ++st.k0; else if (k < 64) ++st.ksmall; else ++st.kbig;
    }
    // exact step
    dense_pivot(gm, M.data(), stride, rows, cols, r, s, det);
    bool out_fits = fits<T>(det.z);
    for (long i = 0; i < rows; i++)
      for (long j = 0; j < cols; j++) out_fits = out_fits && fits<T>(M[i * stride + j].z);
    // fixed-width step
    bool threw = false;
    try {
      ar.pivot(F.data(), stride, rows, cols, r, s, fdet, bound);
    } catch (overflow_error &) {
      threw = true;
    }
    ++st.steps;
    if (threw != !out_fits) {
      if (++st.bad_ovf <= 3) printf("OVERFLOW MISMATCH %s %ldx%ld step %d: threw=%d fits=%d\n", Ar::name, rows, cols, step, threw, out_fits);
      break;
    }
    if (threw) {
      ++st.ovf;
      break;
    }
    bool same = true, sound = true;
    to_mpz(z, fdet);
    same &= mpz_cmp(z, det.z) == 0;
    mpz_t lim;  // |a| <= lim for every entry if bound is known
    mpz_init_set_ui(lim, bound);
    mpz_add_ui(lim, lim, 1);
    if (sizeof(T) == 16) mpz_mul_2exp(lim, lim, 64);
    for (long i = 0; i < rows; i++)
      for (long j = 0; j < stride; j++) {
        const T v = F[i * stride + j];
        if (j >= cols) {
          same &= v == (T)0x5a5a;  // padding untouched
          continue;
        }
        to_mpz(z, v);
        same &= mpz_cmp(z, M[i * stride + j].z) == 0;
        if (bound != ~(std::uint64_t)0 && mpz_cmpabs(z, lim) > 0) sound = false;
      }
    mpz_clear(lim);
    if (!same && ++st.bad <= 3) printf("VALUE MISMATCH %s %ldx%ld step %d\n", Ar::name, rows, cols, step);
    if (!sound && ++st.bad_bound <= 3)
      printf("UNSOUND BOUND %s %ldx%ld step %d: bound %llu\n", Ar::name, rows, cols, step, (unsigned long long)bound);
    if (!same) break;
  }
  mpz_clear(z);
}

// One pivot on a constructed input whose results straddle the fast-path
// limit |q| < 2^(w-k) (w = 63 for I64, 127 for I128), where the wrapping
// arithmetic stops being exact. Any input works as long as every division is
// exact: with det e = 2^k e', pivot Ars = +-e*m, pivot-row entries
// c_j = e*g_j, the results are q_j = (a_j*Ars - b*c_j)/e = +-a_j*m - b*g_j.
// a_j is drawn so that |q_j| is spread over [2^(w-k-2), 2^(w-k+2)].
template <class Ar>
static void boundary(std::mt19937_64 &g, Stats &st) {
  using T = typename Ar::T;
  constexpr int w = sizeof(T) == 8 ? 63 : 127;
  GMP gm;
  Ar ar;
  choose_kernel(ar);
  // cols up to 13: the AVX2 kernel only runs for rows of 8 or more entries
  const long rows = 2 + (long)(g() % 4), cols = 2 + (long)(g() % 12), stride = cols;
  const long r = (long)(g() % rows), s = (long)(g() % cols);
  const int k = (int)(g() % (w - 3));  // det = 2^k e' must fit with room for m
  mpz_t e, m, t, z;
  mpz_inits(e, m, t, z, nullptr);
  mpz_set_ui(e, 1 + 2 * (g() % 50));  // e' odd
  mpz_mul_2exp(e, e, k);
  mpz_set_ui(m, 1 + g() % 7);
  std::vector<GMP::T> M(rows * stride);
  bool ok = true;
  const bool bzero = g() % 2;  // b = 0 makes the bound estimate tight
  for (long j = 0; j < cols; j++) {
    mpz_t &c = M[r * stride + j].z;
    if (j == s) {
      mpz_mul(c, e, m);
      if (g() % 2) mpz_neg(c, c);
    } else {
      mpz_set_si(c, (long)(g() % 7) - 3);  // g_j
      mpz_mul(c, c, e);
    }
    ok &= fits<T>(c);
  }
  for (long i = 0; i < rows; i++) {
    if (i == r) continue;
    for (long j = 0; j < cols; j++) {
      mpz_t &a = M[i * stride + j].z;
      if (j == s) {
        mpz_set_si(a, bzero ? 0 : (long)(g() % 9) - 4);
      } else {
        // |a| m ~ 2^(w-k) * 2^u, u in [-2, 2]
        const int bits = w - k - 2 + (int)(g() % 5);
        mpz_set_ui(a, 1);
        mpz_mul_2exp(a, a, bits);
        mpz_set_ui(t, g());
        mpz_mul_2exp(t, t, bits > 64 ? bits - 64 : 0);
        if (bits < 64) mpz_fdiv_q_2exp(t, t, 64 - bits);
        mpz_add(a, a, t);  // a in [2^bits, 2^(bits+1))
        mpz_fdiv_q(a, a, m);
        if (g() % 2) mpz_neg(a, a);
        if (g() % 6 == 0) mpz_set_si(a, (long)(g() % 5) - 2);
      }
      ok &= fits<T>(a);
    }
  }
  ok &= fits<T>(e);
  if (ok) {
    std::vector<T> F(rows * stride);
    for (size_t q = 0; q < M.size(); q++) F[q] = from_mpz<T>(M[q].z);
    T fdet = from_mpz<T>(e);
    std::uint64_t bound = g() % 2 ? ~(std::uint64_t)0 : Ar::magnitude_bound(F.data(), stride, rows, cols);
    GMP::T det;
    mpz_set(det.z, e);
    dense_pivot(gm, M.data(), stride, rows, cols, r, s, det);
    bool out_fits = fits<T>(det.z);
    for (auto &v : M) out_fits = out_fits && fits<T>(v.z);
    bool threw = false;
    try {
      ar.pivot(F.data(), stride, rows, cols, r, s, fdet, bound);
    } catch (overflow_error &) {
      threw = true;
    }
    ++st.steps;
    if (k == 0) ++st.k0; else if (k < 64) ++st.ksmall; else ++st.kbig;
    if (threw != !out_fits) {
      if (++st.bad_ovf <= 3) printf("BOUNDARY OVERFLOW MISMATCH %s k=%d: threw=%d fits=%d\n", Ar::name, k, threw, out_fits);
    } else if (threw) {
      ++st.ovf;
    } else {
      bool same = true;
      to_mpz(z, fdet);
      same &= mpz_cmp(z, det.z) == 0;
      for (size_t q = 0; q < M.size(); q++) {
        to_mpz(z, F[q]);
        same &= mpz_cmp(z, M[q].z) == 0;
      }
      if (!same && ++st.bad <= 3) printf("BOUNDARY VALUE MISMATCH %s k=%d\n", Ar::name, k);
    }
  }
  mpz_clears(e, m, t, z, nullptr);
}

template <class Ar> static long run(const char *kernel, int trials, std::mt19937_64 &g) {
  Stats st, sb;
  for (int t = 0; t < trials; t++) {
    const long rows = 1 + (long)(g() % 12), cols = 1 + (long)(g() % 14);
    sequence<Ar>(g, st, rows, cols, t % 4);
    boundary<Ar>(g, sb);
  }
  printf("  %-6s %-6s steps=%ld overflow=%ld k=0:%ld 0<k<64:%ld k>=64:%ld known_bound=%ld  value_mismatch=%ld overflow_mismatch=%ld unsound_bound=%ld\n",
         Ar::name, kernel, st.steps, st.ovf, st.k0, st.ksmall, st.kbig, st.fast_known, st.bad, st.bad_ovf, st.bad_bound);
  printf("  %-6s %-6s boundary: steps=%ld overflow=%ld k=0:%ld 0<k<64:%ld k>=64:%ld  value_mismatch=%ld overflow_mismatch=%ld\n",
         Ar::name, kernel, sb.steps, sb.ovf, sb.k0, sb.ksmall, sb.kbig, sb.bad, sb.bad_ovf);
  return st.bad + st.bad_ovf + st.bad_bound + sb.bad + sb.bad_ovf;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s trials seed\n", argv[0]);
    return 2;
  }
  const int trials = atoi(argv[1]);
  std::mt19937_64 g(atoll(argv[2]));
  long fail = 0;
  printf("trials=%d\n", trials);
  g_use_avx2 = false;
  fail += run<I64>("plain", trials, g);
#ifdef LRSNASH_AVX2
  if (__builtin_cpu_supports("avx2")) {
    g_use_avx2 = true;
    fail += run<I64>("avx2", trials, g);
  } else {
    printf("  (no AVX2 on this CPU: AVX2 kernel not tested)\n");
  }
#endif
  fail += run<I128>("plain", trials, g);
  return fail ? 1 : 0;
}
