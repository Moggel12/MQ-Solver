#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "common_utils.h"
#include "fes.h"
#include "mq.h"
#include "mq_config.h"
#include "mq_vec.h"
#include "vector_utils.h"

// TODO: Remove
void print_bits(container_t val, int amount)
{
  for (int i = 0; i < amount; i++)
  {
    printf("%u ", (val >> i) & 1);
  }
  printf("\n");
}

// TODO: Change container_t to a "sub_poly_t" and create "poly_t".
container_vec_t compute_e_k(container_t *mat, container_vec_t *new_sys,
                            container_t *old_sys, size_t sys_len, int l, int n)
{
  container_vec_t deg = VEC_0;
  container_t mon = 0;
  container_t lin_deg = 0, quad_deg = 0;

  for (int v = 0; v < (1 << FIXED_VARS) * 4; v++)
  {
    int v_sys_deg = 0;

    for (int s = 0; s < l; s++)
      mon = GF2_ADD(
          mon, parity(GF2_MUL(old_sys[(v / 4) * sys_len], mat[(v * l) + s]))
                   << s);

    new_sys[0] = VEC_INSERT(new_sys[0], mon, v);
    mon = 0;

    int par;

    for (int i = 1; i <= n; i++)
    {
      for (int s = 0; s < l; s++)
      {
        par = parity(GF2_MUL(old_sys[(v / 4) * sys_len + i], mat[(v * l) + s]));
        mon = GF2_ADD(mon, par << s);
      }

      lin_deg |= mon;

      new_sys[i] = VEC_INSERT(new_sys[i], mon, v);
      mon = 0;
    }
    unsigned int idx = lex_idx(0, 1, n);
    for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
      {
        for (int s = 0; s < l; s++)
        {
          par = parity(
              GF2_MUL(old_sys[(v / 4) * sys_len + idx], mat[(v * l) + s]));
          mon = GF2_ADD(mon, par << s);
        }

        quad_deg |= mon;

        new_sys[idx] = VEC_INSERT(new_sys[idx], mon, v);
        mon = 0;

        idx++;
      }
    }

    lin_deg ^= (lin_deg & quad_deg);
    v_sys_deg = __builtin_popcount(lin_deg) + 2 * __builtin_popcount(quad_deg);
    deg = VEC_INSERT(deg, v_sys_deg, v);
  }

  return deg;  // TODO: Add degree computation
}

void fix_poly(container_t *system, container_t *fixed_system,
              container_t *assignment, int new_n, int n)
{
  fixed_system[0] = system[0];
  for (int i = 1; i <= n; i++)
  {
    if (i <= FIXED_VARS)
    {
      fixed_system[0] ^= (-assignment[i - 1]) & system[i];
    }
    else
    {
      fixed_system[i - FIXED_VARS] = system[i];
    }
  }

  int idx = lex_idx(0, 1, n);
  int new_idx = lex_idx(0, 1, new_n);
  for (int i = 0; i < n; i++)
  {
    for (int j = (i + 1); j < n; j++)
    {
      if (i < FIXED_VARS)
      {
        if (j < FIXED_VARS)
        {
          fixed_system[0] ^=
              ((-assignment[i]) & (-assignment[j])) & system[idx];
        }
        else
        {
          fixed_system[j - FIXED_VARS + 1] ^= (-assignment[i]) & system[idx];
        }
      }
      else
      {
        fixed_system[new_idx] = system[idx];
        new_idx++;
      }
      idx++;
    }
  }
}

uint8_t solve(container_t *system, unsigned int n, unsigned int m,
              container_t *sol)
{
  srand(RSEED);  // Seeding rand

  unsigned int new_n = n - FIXED_VARS;

  size_t sys_len = (n_choose_k(new_n, 2) + new_n + 1);

  container_t *fixed_system =
      calloc((1 << FIXED_VARS) * sys_len, sizeof(container_t));

  container_t assignment[FIXED_VARS] = {0};
  for (container_t i = 0; i < (1 << FIXED_VARS); i++)
  {
    for (int j = 0; j < FIXED_VARS; j++) assignment[j] = (i >> j) & 1;
    fix_poly(system, fixed_system + i * sys_len, assignment, new_n, n);
  }

  unsigned int n1 = (unsigned int)ceil(new_n / 5.4);
  unsigned int l = n1 + 1;
  unsigned int k = 0;

  container_vec_t *rand_sys =
      aligned_alloc(ALIGNMENT, sys_len * sizeof(container_vec_t));

  container_t *rand_mat =
      aligned_alloc(ALIGNMENT, l * sizeof(container_t) * (1 << FIXED_VARS) * 4);

  memset(rand_mat, 0, l * sizeof(container_t) * (1 << FIXED_VARS) * 4);

  for (; k < MAX_HISTORY; k++)
  {
#ifdef _DEBUG
    printf("# Commencing round %u\n", k);
#endif

    memset(rand_sys, 0, sys_len * sizeof(container_vec_t));
    uint8_t exit_code;

    for (int i = 0; i < (1 << FIXED_VARS) * 4; i++)
    {
      exit_code =
          gen_matrix(rand_mat + i * l, l,
                     m);  // TODO: Change to generate non vectorized matrices.
      if (exit_code == 1)
      {
        free(rand_sys);
        free(rand_mat);
        free(fixed_system);
        printf("Matrix error\n");
        return 1;
      }
    }

    container_t curr_potentials = 0;

    container_vec_t w = VEC_SUB(
        compute_e_k(rand_mat, rand_sys, fixed_system, sys_len, l, new_n),
        VEC_ASSIGN_ONE(n1));  // TODO: Change to take in arrays and compute
                              // vector (do not vectorize parities).

    exit_code = fes_recover_vectorized(system, rand_sys, new_n, n1,
                                       VEC_ADD(w, VEC_1), &curr_potentials);
    switch (exit_code)
    {
      case 0:
        free(rand_sys);
        free(rand_mat);
        free(fixed_system);

        *sol = curr_potentials;

        return 0;

      case 1:
        free(rand_sys);
        free(rand_mat);
        free(fixed_system);
        printf("fes_recover mem error\n");
        return 1;

      default:
        break;
    }
  }

  free(rand_sys);
  free(rand_mat);
  free(fixed_system);

  return 1;
}
