#include "benchmark.h"

#include "mq.h"

size_t g_solve_time = 0;
size_t g_recover_time = 0;
size_t g_recover_eval_time = 0;
size_t g_recover_interp_time = 0;
// size_t g_output_time = 0;
size_t g_fes_time = 0;
size_t g_ek_time = 0;
size_t g_matrix_time = 0;
size_t g_eval_time = 0;
size_t g_hist_time = 0;
size_t g_recover_eval = 0;
size_t g_recover_interp = 0;

double g_recover_eval_avg = 0;
double g_recover_interp_avg = 0;

void readout_benchmarks()
{
  READOUT_BENCH(g_solve_time)
  READOUT_BENCH(g_recover_time)
  READOUT_BENCH(g_recover_eval_time)
  READOUT_BENCH(g_recover_interp_time)
  // READOUT_BENCH(g_output_time)
  READOUT_BENCH(g_fes_time)
  READOUT_BENCH(g_ek_time)
  READOUT_BENCH(g_matrix_time)
  READOUT_BENCH(g_eval_time)
  READOUT_BENCH(g_hist_time)
  printf("FES recovery evaluations: %f\n", g_recover_eval_avg);
  printf("FES recovery interpolations: %f\n", g_recover_interp_avg);
  printf("FES recovery ratio (eval/interp): %f\n",
         g_recover_eval_avg / g_recover_interp_avg);
}

void e2e_benchmark(size_t rounds, poly_t *systems[], size_t n, size_t m)
{
  size_t succeeded_r = rounds;
  for (size_t r = 0; r < rounds; r++)
  {
    poly_t sol;
    uint8_t error = solve(systems[r], n, m, &sol);
    if (error)
    {
      succeeded_r--;
      printf("Failed rounds %zu\n", r);
    }
  }

  g_solve_time /= succeeded_r;
  g_recover_time /= succeeded_r;
  g_recover_eval_time /= succeeded_r;
  g_recover_interp_time /= succeeded_r;
  // g_output_time /= succeeded_r;
  g_fes_time /= succeeded_r;
  g_ek_time /= succeeded_r;
  g_matrix_time /= succeeded_r;
  g_eval_time /= succeeded_r;
  g_hist_time /= succeeded_r;
  g_recover_eval_avg = g_recover_eval / (double)succeeded_r;
  g_recover_interp_avg = g_recover_interp / (double)succeeded_r;

  readout_benchmarks();
}