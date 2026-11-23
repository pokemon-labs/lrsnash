// Benchmark for the header-only solvers (src/lrsnash/solve.hpp).
// Same game generator as benchmark.cc / the legacy bench: mt19937(seed),
// payoffs uniform in [0, den+1]; for "full" the column player's matrix is an
// independent draw.
//
//   g++ -O2 -std=c++17 -Isrc tools/bench.cc -lgmp -o bench
//   ./bench rows cols den trials seed fast|full|all|first|any auto|i64|i128|mp|gmp
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>

#include "lrsnash/solve.hpp"

using namespace LRSNash;

template <class S>
static void run(S &s, int rows, int cols, int den, int trials, unsigned seed,
                const std::string &mode) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dis(0, den + 1);
  std::vector<int> a(rows * cols), b(rows * cols);
  std::vector<float> rs(rows), cs(cols);
  double chk = 0;
  long ovf = 0, neq = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int t = 0; t < trials; ++t) {
    for (auto &x : a) x = dis(gen);
    if (mode != "fast")
      for (auto &x : b) x = dis(gen);
    FloatOneSumOutput out{rs.data(), cs.data(), 0.f};
    try {
      if (mode == "fast") {
        FastInput in{rows, cols, a.data(), den};
        s.solve_fast(&in, &out);
      } else if (mode == "all") {
        neq += (long)s.solve_all(rows, cols, a.data(), b.data()).size();
      } else {
        Input in{rows, cols, den, a.data(), b.data()};
        s.solve_full(&in, &out,
                     mode == "first" ? Mode::First
                                     : mode == "any" ? Mode::Any : Mode::Last);
      }
    } catch (overflow_error &) {
      ++ovf;
      continue;
    }
    chk += out.value;
  }
  double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  printf("%dx%d den=%d %-5s trials=%d: %10.2f us/game  overflow=%ld%s chk=%.6f\n",
         rows, cols, den, mode.c_str(), trials, 1e6 * sec / trials, ovf,
         mode == "all" ? (" eq/game=" + std::to_string((double)neq / trials)).c_str() : "", chk);
}

int main(int argc, char **argv) {
  if (argc < 8) {
    fprintf(stderr, "usage: %s rows cols den trials seed fast|full|all|first|any auto|i64|i128|mp|gmp\n", argv[0]);
    return 1;
  }
  int rows = atoi(argv[1]), cols = atoi(argv[2]), den = atoi(argv[3]), trials = atoi(argv[4]);
  unsigned seed = atoi(argv[5]);
  std::string mode = argv[6], ar = argv[7];
  printf("[%-4s] ", ar.c_str());
  if (ar == "i64") { Solver64 s; run(s, rows, cols, den, trials, seed, mode); }
  else if (ar == "i128") { Solver128 s; run(s, rows, cols, den, trials, seed, mode); }
  else if (ar == "mp") { SolverMP s; run(s, rows, cols, den, trials, seed, mode); }
#ifdef LRSNASH_HAVE_GMP
  else if (ar == "gmp") { SolverGMP s; run(s, rows, cols, den, trials, seed, mode); }
#endif
  else { Solver s; run(s, rows, cols, den, trials, seed, mode); }
}
