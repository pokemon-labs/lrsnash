#pragma once

namespace LRSNash {

struct FastInput {
  long rows;
  long cols;
  int *data;
  int den;
};

struct Input {
  long rows;
  long cols;
  int den;
  int *row_data;
  int *col_data;
};

struct FloatOneSumOutput {
  float *row_strategy;
  float *col_strategy;
  float value;
};

void solve_fast(const FastInput *g, FloatOneSumOutput *gg);
void solve_full(const Input *g, FloatOneSumOutput *gg);

}