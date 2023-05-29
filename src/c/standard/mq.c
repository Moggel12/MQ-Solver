#include "mq.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "fes.h"
#include "mq_config.h"
#include "mq_uni.h"
#include "utils.h"

#if defined(_DEBUG)
size_t solver_rounds = 0;
#endif

typedef struct SolutionsStruct
{
  poly_t *solutions;
  size_t amount;
} SolutionsStruct;

poly_t parity(poly_t bits) { return __builtin_parity(bits); }

unsigned int compute_p_k(poly_t *mat, poly_t *new_sys, poly_t *old_sys, int l,
                         int n)
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

uint8_t solve(poly_t *system, unsigned int n, unsigned int m, poly_t *sol)
{
  BEGIN_BENCH(g_solve_time)

  srand(RSEED);  // Seeding rand

  unsigned int amnt_sys_vars = (n_choose_k(n, 2) + n + 1);

  unsigned int n1 = (unsigned int)ceil(n / 5.4);
  unsigned int l = n1 + 1;
  unsigned int k = 0;

  SolutionsStruct *potential_solutions[MAX_HISTORY] = {0};
  size_t hist_progress[MAX_HISTORY] = {0};

  poly_t *rand_sys = malloc(amnt_sys_vars * sizeof(poly_t));

  poly_t *rand_mat = malloc(l * sizeof(poly_t));

  for (; k < MAX_HISTORY; k++)
  {
    printf("# Commencing round %u\n", k);

    memset(rand_sys, 0, amnt_sys_vars * sizeof(poly_t));

    BEGIN_BENCH(g_matrix_time)

    unsigned int error = gen_matrix(rand_mat, l, m);

    END_BENCH(g_matrix_time)

    poly_t *curr_potentials = calloc(
        1 << (n - n1),
        sizeof(poly_t));
    if (error || !curr_potentials)
    {
      k--;
      continue;
    }

    BEGIN_BENCH(g_ek_time)

    unsigned int w = compute_p_k(rand_mat, rand_sys, system, l, n) - n1;

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

    poly_t *tmp = curr_potentials;
    curr_potentials =
        realloc(curr_potentials, len_out * sizeof(poly_t));
    if (!curr_potentials) curr_potentials = tmp;

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
      poly_t idx_solution = curr_potentials[idx];

      BEGIN_BENCH(g_hist_time)

      for (unsigned int k1 = 0; k1 < k; k1++)
      {
        for (; hist_progress[k1] < potential_solutions[k1]->amount;
             hist_progress[k1]++)
        {
          poly_t k1_solution =
              potential_solutions[k1]->solutions[hist_progress[k1]];

          if (GF2_MUL(k1_solution, INT_MASK((n - n1))) > GF2_MUL(idx_solution, INT_MASK((n - n1))))
          {
            break;
          }
          else if (INT_EQ(k1_solution, idx_solution))
          {
            poly_t gray_y = GRAY(GF2_MUL(idx_solution, INT_MASK((n - n1))));
            poly_t solution =
                GF2_ADD(gray_y, GF2_MUL(idx_solution, INT_LSHIFT(INT_MASK(n1), (n - n1))));

            if (!eval(system, n, solution))
            {
              *sol = solution;

              for (unsigned int i = 0; i <= k; i++)
              {
                g_stored_solutions += potential_solutions[i]->amount;
                free(potential_solutions[i]->solutions);
                free(potential_solutions[i]);
              }

              free(rand_sys);
              free(rand_mat);

              END_BENCH(g_solve_time)

#if defined(_DEBUG)
              solver_rounds = k + 1;
#endif

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

#if defined(_DEBUG)
  solver_rounds = MAX_HISTORY;
#endif

  return 1;
}
