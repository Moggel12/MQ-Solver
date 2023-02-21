#include "benchmark.h"

#include "mq.h"

size_t g_solve_time = 0;
size_t g_recover_time = 0;
size_t g_output_time = 0;
size_t g_fes_time = 0;
size_t g_ek_time = 0;
size_t g_matrix_time = 0;

void readout_benchmarks()
{
  READOUT_BENCH(g_solve_time)
  READOUT_BENCH(g_recover_time)
  READOUT_BENCH(g_output_time)
  READOUT_BENCH(g_fes_time)
  READOUT_BENCH(g_ek_time)
  READOUT_BENCH(g_matrix_time)
}

void e2e_benchmark(size_t rounds, poly_t *systems[], size_t n, size_t m)
{
  size_t succeeded_r = rounds;
  for (size_t r = 0; r < rounds; r++)
  {
    unsigned int sol;
    uint8_t error = solve(systems[r], n, m, &sol);
    if (error)
    {
      succeeded_r--;
      printf("Failed rounds %zu\n", r);
    }
  }

  g_solve_time = g_solve_time / succeeded_r;
  g_recover_time = g_recover_time / succeeded_r;
  g_output_time = g_output_time / succeeded_r;
  g_fes_time = g_fes_time / succeeded_r;
  g_ek_time = g_ek_time / succeeded_r;
  g_matrix_time = g_matrix_time / succeeded_r;

  readout_benchmarks();
}