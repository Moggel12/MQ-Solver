#include "mq.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "fes.h"
#include "mq_config.h"
#include "utils.h"

typedef struct SolutionsStruct
{
  vars_t *solutions;
  size_t amount;
} SolutionsStruct;

unsigned int compute_e_k(poly_t *mat, poly_t *new_sys, poly_t *old_sys, int l,
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
        if ((s_poly_deg < 2) && ((par << 1) == 2)) s_poly_deg = par << 1;
      }
    }

    deg += s_poly_deg;
  }
  return deg;
}

uint8_t output_potentials(poly_t *system, unsigned int n, unsigned int n1,
                          unsigned int w, vars_t *out, size_t *out_size)
{
  vars_t *evals = calloc(1 << (n - n1), sizeof(vars_t));
  if (!evals) return 1;

  BEGIN_BENCH(g_recover_time)

  uint8_t error = fes_recover(system, n, n1, w + 1, evals);

  if (error) return error;

  END_BENCH(g_recover_time)

  *out_size = 0;

  for (int y_hat = 0; y_hat < (1 << (n - n1)); y_hat++)
  {
    if (!VARS_IS_ZERO(VARS_IDX(evals[y_hat], 0)))
    {
      out[*out_size] = y_hat;
      vars_t z_bits =
          VARS_RSHIFT(GF2_ADD(evals[y_hat], VARS_MASK((n1 + 1))), 1);

      out[(*out_size)] =
          GF2_ADD(out[(*out_size)], VARS_LSHIFT(z_bits, (n - n1)));
      (*out_size)++;
    }
  }
  free(evals);
  return 0;
}

// TODO: Reconsider error handling and it's performance impact.
uint8_t solve(poly_t *system, unsigned int n, unsigned int m, vars_t *sol)
{
  BEGIN_BENCH(g_solve_time)

  srand(RSEED);  // Seeding rand

  unsigned int amnt_sys_vars = (n_choose_k(n, 2) + n + 1);

  unsigned int n1 = (unsigned int)ceil(n / 5.4);
  unsigned int l = n1 + 1;
  unsigned int k = 0;

  SolutionsStruct *potential_solutions[MAX_HISTORY] = {0};

  // TODO: poly_t is not guaranteed large enough for l bits.
  poly_t *rand_sys = malloc(amnt_sys_vars * sizeof(poly_t));

  poly_t *rand_mat = malloc(l * sizeof(poly_t));

  for (; k < MAX_HISTORY; k++)
  {
#ifdef _DEBUG
    printf("# Commencing round %u\n", k);
#endif

    memset(rand_sys, 0, amnt_sys_vars * sizeof(poly_t));

    BEGIN_BENCH(g_matrix_time)

    unsigned int error = gen_matrix(rand_mat, l, m);

    END_BENCH(g_matrix_time)

    vars_t *curr_potentials = calloc(
        1 << (n - n1), sizeof(vars_t));  // TODO: Change this to a suitable size
    if (error || !curr_potentials)
    {
      k--;
      continue;
    }

    BEGIN_BENCH(g_ek_time)

    unsigned int w = compute_e_k(rand_mat, rand_sys, system, l, n) - n1;

    END_BENCH(g_ek_time)

    BEGIN_BENCH(g_output_time)

    size_t len_out = 0;
    error = output_potentials(rand_sys, n, n1, w, curr_potentials, &len_out);
    if (error)
    {
      free(rand_sys);
      free(rand_mat);
      return 1;
    }
    END_BENCH(g_output_time)

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
      vars_t y = GF2_MUL(curr_potentials[idx], VARS_MASK((n - n1)));

      for (unsigned int k1 = 0; k1 < k; k1++)
      {
        for (size_t old_idx = 0; old_idx < potential_solutions[k1]->amount;
             old_idx++)
        {
          vars_t old_y = GF2_MUL(potential_solutions[k1]->solutions[old_idx],
                                 VARS_MASK((n - n1)));

          if (old_y > y)
          {
            break;
          }
          else if (VARS_EQ(curr_potentials[idx],
                           potential_solutions[k1]->solutions[old_idx]))
          {
            if (!eval(system, n, curr_potentials[idx]))
            {
              *sol = curr_potentials[idx];

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
    }
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
