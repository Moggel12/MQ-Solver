#include "benchmark.h"

#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "mq.h"
#include "fes.h"

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

  int pid = getpid();
  char cmd[100];

  sprintf(cmd, "cat /proc/%d/status | grep Vm > procinfo_%d_%zu_%zu\n", pid, pid, n, m);
  system(cmd);
}

#if !defined(REG256) && !defined(REG128)

void fes_benchmark(size_t rounds, poly_t *systems[], size_t n, size_t m)
{
  clock_t fes_time;
  size_t succeeded_r = rounds;

  g_solve_time = 0;

  poly_t *solutions = malloc((1 << n) * sizeof(poly_t));

  for (int r = 0; r < rounds; r++)
  {
    clock_t current_time = clock();

    fes(systems[r], n, m, solutions);

    fes_time += (current_time - fes_time);
  }

  g_solve_time = fes_time / succeeded_r;
  
  size_t msec = (g_solve_time) * 1000 / CLOCKS_PER_SEC;
  printf("FES solve time: %zus, %zums\n", msec / 1000, msec % 1000);

  int pid = getpid();
  char cmd[100];

  sprintf(cmd, "cat /proc/%d/status | grep Vm > procinfo_fes_%d_%zu_%zu\n", pid, pid, n, m);
  system(cmd);
}

#endif