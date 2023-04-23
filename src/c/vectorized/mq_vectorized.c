#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "fes.h"
#include "mq.h"
#include "mq_config.h"
#include "mq_vec.h"
#include "utils.h"

container_vec_t compute_e_k(container_vec_t *mat, container_vec_t *new_sys,
                            container_t *old_sys, int l, int n)
{
  container_vec_t deg = VEC_0;

  for (int s = 0; s < l; s++)
  {
    container_vec_t s_poly_deg = VEC_0;

    new_sys[0] = VEC_GF2_ADD(
        new_sys[0],
        VEC_LSHIFT(parity(VEC_GF2_MUL(VEC_ASSIGN_ONE(old_sys[0]), mat[s])), s));

    container_vec_t par;

    for (int i = 1; i <= n; i++)
    {
      par = parity(VEC_GF2_MUL(VEC_ASSIGN_ONE(old_sys[i]), mat[s]));

      new_sys[i] = GF2_ADD(new_sys[i], VEC_LSHIFT(par, s));
      // printf("#\n");
      // print_register(s_poly_deg);
      // print_register(VEC_EQ(s_poly_deg, VEC_0));
      // printf("#\n");
      s_poly_deg = VEC_BLEND(s_poly_deg, par, VEC_EQ(s_poly_deg, VEC_0));
    }
    unsigned int idx = lex_idx(0, 1, n);
    for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
      {
        par = parity(VEC_GF2_MUL(VEC_ASSIGN_ONE(old_sys[idx]), mat[s]));

        new_sys[idx] = VEC_GF2_ADD(new_sys[idx], VEC_LSHIFT(par, s));

        idx++;

        container_vec_t cond =
            VEC_GF2_MUL(VEC_LT(s_poly_deg, VEC_ASSIGN_ONE(2)),
                        VEC_EQ(par, VEC_ASSIGN_ONE(1)));
        // print_register(s_poly_deg);
        // printf("#\n");
        // print_register(s_poly_deg);
        // print_register(cond);
        // printf("#\n");
        s_poly_deg = VEC_BLEND(s_poly_deg, VEC_ASSIGN_ONE(2), cond);
      }
    }
    // printf("===\n");
    deg = VEC_ADD(deg, s_poly_deg);
  }
  return deg;
}

uint8_t solve(container_t *system, unsigned int n, unsigned int m,
              container_t *sol)
{
  srand(RSEED);  // Seeding rand

  size_t amnt_sys_vars = (n_choose_k(n, 2) + n + 1);

  unsigned int n1 = (unsigned int)ceil(n / 5.4);
  unsigned int l = n1 + 1;
  unsigned int k = 0;

  container_vec_t *rand_sys =
      aligned_alloc(ALIGNMENT, amnt_sys_vars * sizeof(container_vec_t));

  container_vec_t *rand_mat =
      aligned_alloc(ALIGNMENT, l * sizeof(container_vec_t));

  // memset(rand_sys, 0, amnt_sys_vars * sizeof(container_vec_t));
  memset(rand_mat, 0, l * sizeof(container_vec_t));

  for (; k < MAX_HISTORY; k++)
  {
#ifdef _DEBUG
    printf("# Commencing round %u\n", k);
#endif

    memset(rand_sys, 0, amnt_sys_vars * sizeof(container_vec_t));

    uint8_t exit_code = gen_matrix(rand_mat, l, m);
    if (exit_code == 1)
    {
      free(rand_sys);
      free(rand_mat);

      return 1;
    }

    container_t curr_potentials = 0;

    container_vec_t w = VEC_SUB(compute_e_k(rand_mat, rand_sys, system, l, n),
                                VEC_ASSIGN_ONE(n1));
    // print_register(w);
    // if (k == 0)
    // {
    //   for (int i = 0; i < amnt_sys_vars; i++)
    //   {
    //     print_register(rand_sys[i]);
    //   }
    // }
    // printf("#\n");

    exit_code = fes_recover_vectorized(system, rand_sys, n, n1,
                                       VEC_ADD(w, VEC_1), &curr_potentials);
    switch (exit_code)
    {
      case 0:
        free(rand_sys);
        free(rand_mat);

        *sol = curr_potentials;

        return 0;

      case 1:
        free(rand_sys);
        free(rand_mat);

        return 1;

      default:
        break;
    }
  }

  free(rand_sys);
  free(rand_mat);

  return 1;
}
