// Prints, exactly, every extreme equilibrium that the lrsnash enumeration
// (src/lrsnash/nash.hpp, GMP arithmetic) reports for games read from stdin.
// Used by tools/check_complete.py.
//
// Input:  number of games, then for each game: rows cols, A (rows*cols ints,
//         row-major, player 1), B (rows*cols ints, row-major, player 2).
// Output: for each game "G <index>", then one line per reported equilibrium
//           E <d1> <x_1..x_rows> <v> | <d2> <y_1..y_cols> <u>
//         (numerators over d1 / d2; v is player 2's payoff, u player 1's),
//         then "END".
//   g++ -O2 -std=c++17 -Isrc tools/dump_eq.cc -lgmp -o dump_eq
#include "lrsnash/nash.hpp"
#include <gmp.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using LRSNash::arith::GMP;

static std::string fmt(const GMP::T *o, long n) {
  std::string s;
  for (long k = 0; k < n; k++) {
    char *c = mpz_get_str(nullptr, 10, o[k].z);
    if (k) s += ' ';
    s += c;
    free(c);
  }
  return s;
}

// on_eq2 reports the y's of the current x's face, then on_eq1 reports x.
struct Dump {
  std::vector<std::string> ys;
  bool on_eq2(const GMP::T *o, long n) {
    ys.push_back(fmt(o, n));
    return true;
  }
  bool on_eq1(const GMP::T *o, long n) {
    const std::string x = fmt(o, n);
    for (auto &y : ys) printf("E %s | %s\n", x.c_str(), y.c_str());
    ys.clear();
    return true;
  }
};

int main() {
  long games;
  if (scanf("%ld", &games) != 1) return 2;
  LRSNash::detail::Nash<GMP> en;
  for (long g = 0; g < games; g++) {
    long m, n;
    if (scanf("%ld %ld", &m, &n) != 2) return 2;
    std::vector<int> A(m * n), B(m * n);
    for (auto &v : A)
      if (scanf("%d", &v) != 1) return 2;
    for (auto &v : B)
      if (scanf("%d", &v) != 1) return 2;
    printf("G %ld\n", g);
    Dump d;
    en.solve(m, n, [&](long i, long j) { return (long)A[i * n + j]; },
             [&](long i, long j) { return (long)B[i * n + j]; }, d);
    if (!d.ys.empty()) printf("PENDING %zu\n", d.ys.size());  // must not happen
    printf("END\n");
  }
}
