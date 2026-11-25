Notice: This software contains modifications to the original source code of the program 'lrslib-072', which is licensed under the GNU General Public License (GPL) Version 2, June 1991. The modifications made were to remove code that was not necessary to compile 'lrsnash1' and to make the code available as a library for C++ projects.

The modifications made are released under the same GPL license terms as the original software. The original copyright and license notice is preserved.

## Usage

Header-only C++17. Add `src/` to the include path (or link the CMake target
`lrsnash`), and link GMP if it is installed:

```cpp
#include "lib.h"

LRSNash::Solver64 solver;  // keep one per thread and reuse it
LRSNash::ConstantSumGame g{rows, cols, data, den};
LRSNash::FloatEquilibrium out{row_strategy, col_strategy, 0.f};
solver.solve(&g, &out);
```

`src/lib.h` gives an overview, and `src/lrsnash/solve.hpp` documents the whole
API: general-sum and rational games, every extreme equilibrium (also exactly),
and the choice of arithmetic (int64, int128, lrsmp, GMP).
