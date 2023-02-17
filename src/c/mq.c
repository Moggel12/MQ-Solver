#include "mq.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fes.h"
#include "utils.h"

typedef struct SolutionsStruct
{
  vars_t *solutions;
  size_t amount;
} SolutionsStruct;

void print_bits(unsigned int len, vars_t bits)
{
  for (unsigned int i = 0; i < len; i++)
  {
    printf("%u", (bits >> i) & 1);
  }
}

unsigned int compute_e_k(poly_t *mat, poly_t *new_sys, poly_t *old_sys, int l,
                         int n)
{
  unsigned int deg = 0;
  for (int s = 0; s < l; s++)
  {
    int s_poly_deg = 0;
    new_sys[0] =
        GF2_ADD(new_sys[0], __builtin_parity(GF2_MUL(old_sys[0], mat[s])) << s);

    int parity;

    for (int i = 1; i <= n; i++)
    {
      parity = __builtin_parity(GF2_MUL(old_sys[i], mat[s]));

      new_sys[i] = GF2_ADD(new_sys[i], parity << s);

      s_poly_deg = parity;
    }
    unsigned int idx = lex_idx(0, 1, n);
    for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
      {
        parity = __builtin_parity(GF2_MUL(old_sys[idx], mat[s]));

        new_sys[idx] = GF2_ADD(new_sys[idx], parity << s);

        idx++;

        s_poly_deg = parity << 1;
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

  uint8_t error = fes_recover(system, n, n1, w + 1, evals);
  if (error) return error;

  *out_size = 0;

  for (int y_hat = 0; y_hat < (1 << (n - n1)); y_hat++)
  {
    // printf("%zu\n", *out_size);
    for (unsigned int i = 0; i < n1 + 1; i++)
      //   printf("%u ", (evals[y_hat] >> i) & 1);
      // printf("\n");
      if ((evals[y_hat] & 1) == 1)
      {
        out[*out_size] = y_hat;
        vars_t z_bits = ~(evals[y_hat] >> 1);

        out[(*out_size)] = GF2_ADD(out[(*out_size)], (z_bits << (n - n1)));
        (*out_size)++;
      }
  }

  return 0;
}

// TODO: Reconsider error handling and it's performance impact.
uint8_t solve(poly_t *system, unsigned int n, unsigned int m, vars_t *sol)
{
  srand(RSEED);  // Seeding rand

  unsigned int amnt_sys_vars = (n_choose_k(n, 2) + n + 1);

  unsigned int n1 = (unsigned int)ceil(n / 5.4);
  unsigned int l = n1 + 1;
  unsigned int k = 0;

  printf("%u, %u, %u, %u\n", n, m, n1, l);

  SolutionsStruct *potential_solutions[MAX_HISTORY] = {0};

  // printf("Initialized list!\n");

  // TODO: poly_t is not guaranteed large enough for l bits.
  poly_t *rand_sys = malloc(amnt_sys_vars * sizeof(poly_t));

  poly_t *rand_mat = malloc(l * sizeof(poly_t));

  for (; k < MAX_HISTORY; k++)
  {
    memset(rand_sys, 0, amnt_sys_vars * sizeof(poly_t));

    printf("Round (%u)\n", k);

    unsigned int error = gen_matrix(rand_mat, l, m);
    // printf("Generated matrix\n");

    for (unsigned int i = 0; i < l; i++)
    {
      for (unsigned int j = 0; j < m; j++)
      {
        printf("%u ", (rand_mat[i] >> j) & 1);
      }
      printf("\n");
    }

    vars_t *curr_potentials = calloc(
        1 << (n - n1), sizeof(vars_t));  // TODO: Change this to a suitable size
    if (error || !curr_potentials)
    {
      k--;
      continue;
    }
    // printf("Allocated curr_potentials!\n");

    // printf("Memory before compute\n");
    // for (unsigned int i = 0; i < n_choose_k(n, 2) + n + 1; i++)
    // {
    //   printf("%u ", rand_sys[i]);
    // }
    // printf("\n");
    unsigned int w = compute_e_k(rand_mat, rand_sys, system, l, n);
    // printf("compute_e_k\n");

    // printf("Memory after compute\n");
    for (unsigned int i = 0; i < n_choose_k(n, 2) + n + 1; i++)
    {
      printf("%u ", rand_sys[i]);
    }
    printf("\n");
    // printf("Computed new system!\n");

    // append_list(potential_solutions, out, 1 << (n - n1));

    size_t len_out = 0;
    error = output_potentials(rand_sys, n, n1, w, curr_potentials, &len_out);
    if (error)
    {
      free(rand_sys);
      free(rand_mat);
      return 1;
    }
    // printf("ran output_solutions (Got %zu solutions)\n", len_out);

    // for (unsigned int i = 0; i < len_out; i++)
    // {
    //   printf("%u\n", curr_potentials[i]);
    // }

    // printf("Successfully ran output_potentials (Got %zu solutions)\n",
    // len_out);

    // potential_solutions[k]->solutions = curr_potentials;
    // potential_solutions[k]->amount = len_out;
    SolutionsStruct *s = malloc(sizeof(SolutionsStruct));
    if (!s)
    {
      for (unsigned int l = 0; l < k; l++)
      {
        free(potential_solutions[l]->solutions);
        free(potential_solutions[l]);
      }
      free(rand_sys);
      free(rand_mat);
      return 1;
    }
    // printf("allocated struct\n");

    s->solutions = curr_potentials, s->amount = len_out;
    potential_solutions[k] = s;

    // printf("Appended list!\n");
    // printf("Added solutions\n");

    // TODO: Sequentially go through stored solutions and compare, instead of
    // using lookup

    for (size_t idx = 0; idx < len_out; idx++)
    {
      size_t bin_y = gray_to_bin(curr_potentials[idx] & ((1 << (n - n1)) - 1));
      for (unsigned int k1 = 0; k1 < k; k1++)
      {
        for (size_t old_idx = 0; old_idx < potential_solutions[k1]->amount;
             old_idx++)
        {
          size_t bin_k1 =
              gray_to_bin(potential_solutions[k1]->solutions[old_idx] &
                          ((1 << (n - n1)) - 1));
          // printf("%zu > %zu?\n", bin_k1, bin_y);
          if (bin_k1 > bin_y)
          {
            break;
          }
          else if (curr_potentials[idx] ==
                   potential_solutions[k1]->solutions[old_idx])
          {
            if (!eval(system, n, curr_potentials[idx]))
            {
              // printf("Found something\n");
              for (unsigned int l = 0; l < k; l++)
              {
                free(potential_solutions[l]->solutions);
                free(potential_solutions[l]);
              }
              free(rand_sys);
              free(rand_mat);
              *sol = curr_potentials[idx];
              return 0;
            }
          }
        }
      }
    }
  }
  // printf("Found nothing!\n");
  for (unsigned int l = 0; l < k; l++)
  {
    free(potential_solutions[l]->solutions);
    free(potential_solutions[l]);
  }
  free(rand_sys);
  free(rand_mat);
  return 1;
}