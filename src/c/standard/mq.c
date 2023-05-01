#include "mq.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "common_utils.h"
#include "fes.h"
#include "mq_config.h"
#include "standard_utils.h"

typedef struct SolutionsStruct
{
  PotentialSolution *solutions;
  size_t amount;
} SolutionsStruct;

unsigned int compute_e_k(container_t *mat, container_t *new_sys,
                         container_t *old_sys, int l, int n)
{
  unsigned int deg = 0;
  for (int s = 0; s < l; s++)
  {
    int s_poly_deg = 0;

    new_sys[0] = GF2_ADD(new_sys[0], parity(GF2_MUL(old_sys[0], mat[s])) << s);

    int par;

    for (int i = 1; i <= n; i++)
    {
      par = parity(GF2_MUL(old_sys[i], mat[s]));

      new_sys[i] = GF2_ADD(new_sys[i], par << s);

      if (s_poly_deg == 0)
      {
        s_poly_deg = par;
      }
    }
    unsigned int idx = lex_idx(0, 1, n);
    for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
      {
        par = parity(GF2_MUL(old_sys[idx], mat[s]));

        new_sys[idx] = GF2_ADD(new_sys[idx], par << s);

        idx++;
        if ((s_poly_deg < 2) && (par == 1)) s_poly_deg = 2;
      }
    }

    deg += s_poly_deg;
  }
  return deg;
}

uint8_t solve(container_t *system, unsigned int n, unsigned int m,
              container_t *sol)
{
  BEGIN_BENCH(g_solve_time)

  srand(RSEED);  // Seeding rand

  unsigned int amnt_sys_vars = (n_choose_k(n, 2) + n + 1);

  unsigned int n1 = (unsigned int)ceil(n / 5.4);
  unsigned int l = n1 + 1;
  unsigned int k = 0;

  SolutionsStruct *potential_solutions[MAX_HISTORY] = {0};
  size_t hist_progress[MAX_HISTORY] = {0};

  container_t *rand_sys = malloc(amnt_sys_vars * sizeof(container_t));

  container_t *rand_mat = malloc(l * sizeof(container_t));

  for (; k < MAX_HISTORY; k++)
  {
#ifdef _DEBUG
    printf("# Commencing round %u\n", k);
#endif

    memset(rand_sys, 0, amnt_sys_vars * sizeof(container_t));

    BEGIN_BENCH(g_matrix_time)

    unsigned int error = gen_matrix(rand_mat, l, m);

    END_BENCH(g_matrix_time)

    PotentialSolution *curr_potentials = calloc(
        1 << (n - n1),
        sizeof(PotentialSolution));  // TODO: Change this to a suitable size
    if (error || !curr_potentials)
    {
      k--;
      continue;
    }

    BEGIN_BENCH(g_ek_time)

    unsigned int w = compute_e_k(rand_mat, rand_sys, system, l, n) - n1;

    END_BENCH(g_ek_time)

    size_t len_out = 0;

    BEGIN_BENCH(g_recover_time);

    error = fes_recover(rand_sys, n, n1, w + 1, curr_potentials, &len_out);

    END_BENCH(g_recover_time);

    if (error)
    {
      free(rand_sys);
      free(rand_mat);
      return 1;
    }

    SolutionsStruct *s = malloc(sizeof(SolutionsStruct));
    if (!s)
    {
      for (unsigned int i = 0; i <= k; i++)
      {
        free(potential_solutions[l]->solutions);
        free(potential_solutions[l]);
      }
      free(rand_sys);
      free(rand_mat);
      return 1;
    }

    s->solutions = curr_potentials, s->amount = len_out;
    potential_solutions[k] = s;

    for (size_t idx = 0; idx < len_out; idx++)
    {
      PotentialSolution idx_solution = curr_potentials[idx];

      BEGIN_BENCH(g_hist_time)

      for (unsigned int k1 = 0; k1 < k; k1++)
      {
        for (; hist_progress[k1] < potential_solutions[k1]->amount;
             hist_progress[k1]++)
        {
          PotentialSolution k1_solution =
              potential_solutions[k1]->solutions[hist_progress[k1]];

          if (k1_solution.y_idx > idx_solution.y_idx)
          {
            k1++;
            break;
          }
          else if (INT_EQ(idx_solution.y_idx, k1_solution.y_idx) &&
                   INT_EQ(idx_solution.z_bits, k1_solution.z_bits))
          {
            container_t gray_y = idx_solution.y_idx ^ (idx_solution.y_idx >> 1);
            container_t solution =
                GF2_ADD(gray_y, INT_LSHIFT(idx_solution.z_bits, (n - n1)));
            if (!eval(system, n, solution))
            {
              *sol = solution;

              for (unsigned int i = 0; i <= k; i++)
              {
                free(potential_solutions[i]->solutions);
                free(potential_solutions[i]);
              }

              free(rand_sys);
              free(rand_mat);

              END_BENCH(g_solve_time)

              return 0;
            }
          }
        }
      }

      END_BENCH(g_hist_time)
    }

    memset(hist_progress, 0, sizeof(size_t) * MAX_HISTORY);
  }

  for (unsigned int i = 0; i < MAX_HISTORY; i++)
  {
    free(potential_solutions[i]->solutions);
    free(potential_solutions[i]);
  }
  free(rand_sys);
  free(rand_mat);

  return 1;
}
