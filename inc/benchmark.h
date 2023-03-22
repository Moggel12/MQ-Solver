#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <stdio.h>
#include <time.h>

#include "mq_config.h"

#define VAR(TYPE, Z) TYPE g_iLine_##Z

#define BEGIN_BENCH(ID) VAR(clock_t, bench_##ID) = clock();

#define END_BENCH(ID) ID += (clock() - g_iLine_bench_##ID);

#define READOUT_BENCH(ID)                        \
  size_t msec_##ID = ID * 1000 / CLOCKS_PER_SEC; \
  printf(#ID " timings: %zus, %zums\n", msec_##ID / 1000, msec_##ID % 1000);

extern size_t g_solve_time;
extern size_t g_recover_time;
extern size_t g_recover_eval_time;
extern size_t g_recover_interp_time;
// extern size_t g_output_time;
extern size_t g_fes_time;
extern size_t g_ek_time;
extern size_t g_matrix_time;
extern size_t g_eval_time;
extern size_t g_hist_time;
extern size_t g_recover_eval;
extern size_t g_recover_interp;

void readout_benchmarks();

void e2e_benchmark(size_t rounds, poly_t *systems[], size_t n, size_t m);

#endif  // !BENCHMARK_H