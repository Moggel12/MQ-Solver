#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "benchmark.h"
#include "binom.h"
#include "fes.h"
#include "mq.h"

state *init_state(unsigned int n, unsigned int n1, uint8_t *prefix)
{
  state *s = aligned_alloc(ALIGNMENT, sizeof(state));
  if (!s) return NULL;

  s->d1 = aligned_alloc(ALIGNMENT, n1 * sizeof(poly_vec_t));

  if (!(s->d1))
  {
    free(s);

    return NULL;
  }

  s->prefix = calloc(n - n1, sizeof(uint8_t));

  if (!(s->prefix))
  {
    free(s->d1);
    free(s);

    return NULL;
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    s->prefix[i] = prefix[i];
  }

  s->d2 = aligned_alloc(ALIGNMENT, n1 * n1 * sizeof(poly_vec_t));

  if (!(s->d2))
  {
    if (!prefix) free(s->prefix);

    free(s->d1);
    free(s);

    return NULL;
  }
  memset(s->d1, 0, n1 * sizeof(poly_vec_t));
  memset(s->d2, 0, n1 * n1 * sizeof(poly_vec_t));

  s->i = 0;
  s->y = VEC_0;

  return s;
}

void destroy_state(state *s)
{
  if (!s) return;
  free(s->d1);
  free(s->d2);
  free(s->prefix);
  free(s);
}

unsigned int bit1(poly_t i) { return trailing_zeros(i); }

unsigned int bit2(poly_t i) { return bit1(i ^ (i & -i)); }

// Assumes arr has been allocated with arr_len bits.
poly_t bits(poly_t i, unsigned int *arr, poly_t arr_len)
{
  if (i == 0)
  {
    if (arr_len > 0) arr[0] = (-1u);
    return 0;
  }

  unsigned int sum = 0;

  for (unsigned int j = 0; j < arr_len; j++)
  {
    arr[j] = bit1(i);

    if (arr[j] == (-1u)) break;

    sum++;

    i = GF2_ADD(i, (i & -i));
  }

  return sum;
}


unsigned int monomial_to_index(size_t mon, unsigned int n,
                               unsigned int boundary)
{
  unsigned int i;
  int d = 0;
  unsigned int index = 0;
  unsigned int index_d = 0;
  for (i = 0; i < n; i++)
    if (!(((mon >> i) & 1) == 0) && (i <= boundary))
    {
      d++;
      index += lk_binom[i * BINOM_DIM2 + d];
      index_d += lk_binom[n * BINOM_DIM2 + d];
    }
  index = index_d - index;

  return index;
}

state *init(state *s, poly_vec_t *systems, unsigned int n, unsigned int n1,
            uint8_t *prefix)
{
  s = init_state(n, n1, prefix);

  if (!s)
  {
    return NULL;
  }

  s->y = systems[0];

  for (unsigned int k = 0; k < n1; k++)
  {
    for (unsigned int j = 0; j < k; j++)
    {
      s->d2[k * n1 + j] = systems[lex_idx(j + (n - n1), k + (n - n1), n)];
    }
  }

  s->d1[0] = systems[1 + n - n1];

  for (unsigned int k = 1; k < n1; k++)
  {
    s->d1[k] = VEC_GF2_ADD(s->d2[k * n1 + (k - 1)], systems[1 + k + (n - n1)]);
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (prefix[i] == 0) continue;

    for (unsigned int k = 0; k < n1; k++)
    {
      s->d1[k] = VEC_GF2_ADD(s->d1[k], systems[lex_idx(i, k + (n - n1), n)]);
    }

    s->y = VEC_GF2_ADD(s->y, systems[i + 1]);
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (prefix[i] == 0) continue;

    for (unsigned int j = i + 1; j < (n - n1); j++)
    {
      if (prefix[j] == 0) continue;

      s->y = VEC_GF2_ADD(s->y, systems[lex_idx(i, j, n)]);
    }
  }

  return s;
}

state *update(state *s, poly_vec_t *systems, unsigned int n, unsigned int n1,
              uint8_t *prefix)
{
  if (!s)
  {
    return init(s, systems, n, n1, prefix);
  }

  uint8_t *on = malloc(n - n1);

  if (!on)
  {
    return NULL;
  }

  uint8_t *off = malloc(n - n1);

  if (!off)
  {
    free(on);
    return NULL;
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    off[i] = (prefix[i] == 0) && (s->prefix[i] == 1) ? 1 : 0;
    on[i] = (s->prefix[i] == 0) && (prefix[i] == 1) ? 1 : 0;
  }

  for (unsigned int idx = 0; idx < (n - n1); idx++)
  {
    if (off[idx] == 0) continue;

    for (unsigned int k = 0; k < n1; k++)
    {
      unsigned int small_idx = idx;
      unsigned int large_idx = k + (n - n1);

      if (idx > k + (n - n1))
      {
        small_idx = k + (n - n1);
        large_idx = idx;
      }

      s->d1[k] =
          VEC_GF2_ADD(s->d1[k], systems[lex_idx(small_idx, large_idx, n)]);
    }

    s->y = VEC_GF2_ADD(s->y, systems[idx + 1]);
  }

  for (unsigned int i = 0; i < (n - n1); ++i)
  {
    if (off[i] == 0) continue;

    for (unsigned int j = 0; j < (n - n1); j++)
    {
      if (!((s->prefix[j] == 1) && (off[j] == 0))) continue;

      unsigned int small_idx = i;
      unsigned int large_idx = j;

      if (i > j)
      {
        small_idx = j;
        large_idx = i;
      }

      s->y = VEC_GF2_ADD(s->y, systems[lex_idx(small_idx, large_idx, n)]);
    }
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (off[i] == 0) continue;

    for (unsigned int j = i + 1; j < (n - n1); j++)
    {
      if (off[j] == 0) continue;

      s->y = VEC_GF2_ADD(s->y, systems[lex_idx(i, j, n)]);
    }
  }

  // Turn new variables on
  for (unsigned int idx = 0; idx < (n - n1); idx++)
  {
    if (on[idx] == 0) continue;

    for (unsigned int k = 0; k < n1; k++)
    {
      unsigned int small_idx = idx;
      unsigned int large_idx = k + (n - n1);

      if (idx > k + (n - n1))
      {
        small_idx = k + (n - n1);
        large_idx = idx;
      }

      s->d1[k] =
          VEC_GF2_ADD(s->d1[k], systems[lex_idx(small_idx, large_idx, n)]);
    }

    s->y = VEC_GF2_ADD(s->y, systems[idx + 1]);
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (on[i] == 0) continue;

    for (unsigned int j = 0; j < (n - n1); j++)
    {
      if (!((prefix[j] == 1) && (on[j] == 0))) continue;

      unsigned int small_idx = i;
      unsigned int large_idx = j;

      if (i > j)
      {
        small_idx = j;
        large_idx = i;
      }

      s->y = VEC_GF2_ADD(s->y, systems[lex_idx(small_idx, large_idx, n)]);
    }
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (on[i] == 0) continue;

    for (unsigned int j = i + 1; j < (n - n1); j++)
    {
      if (on[j] == 0) continue;

      s->y = VEC_GF2_ADD(s->y, systems[lex_idx(i, j, n)]);
    }
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    s->prefix[i] = prefix[i];
  }
  free(off);
  free(on);

  return s;
}

static inline void step(state *s, unsigned int n1)
{
  s->i = INC(s->i);

  unsigned int alpha1 = bit1(s->i);
  unsigned int alpha2 = bit2(s->i);

  if (alpha2 < (-1u))
  {
    s->d1[alpha1] = VEC_GF2_ADD(s->d1[alpha1], s->d2[alpha2 * n1 + alpha1]);
  }

  s->y = VEC_GF2_ADD(s->y, s->d1[alpha1]);
}

state *fes_eval_parity(poly_vec_t *systems, unsigned int n, unsigned int n1,
                       uint8_t *prefix, state *s, poly_vec_t *parities)
{
  if (!s)
  {
    s = init(s, systems, n, n1, prefix);

    if (!s)
    {
      return NULL;
    }
  }
  int pre_x = 0;
  for (unsigned int pos = 0; pos < (n - n1); pos++)
  {
    pre_x += prefix[pos] * (1 << pos);
  }

  poly_vec_t zero_mask = VEC_IS_ZERO(s->y);

  poly_vec_t added = VEC_GF2_ADD(*parities, VEC_MASK((n1 + 1)));

  *parities = VEC_BLEND(*parities, added, zero_mask);

  while (s->i < ((((size_t) 1) << n1) - 1))
  {
    step(s, n1);

    sub_poly_t z = GRAY(s->i);

    zero_mask = VEC_IS_ZERO(s->y);

    added = VEC_SETBIT(*parities, 0, 1);

    *parities = VEC_BLEND(*parities, added, zero_mask);

    for (unsigned int pos = 0; pos < n1; pos++)
    {
      poly_vec_t z_mask =
          VEC_AND(zero_mask, VEC_ASSIGN_ONE(-INT_IS_ZERO(INT_IDX(z, pos))));

      added = VEC_SETBIT(*parities, (pos + 1), 1);

      *parities = VEC_BLEND(*parities, added, z_mask);
    }
  }

  s->i = 0;
  for (unsigned int i = 0; i < (n1 - 1); i++)
  {
    s->d1[i] = VEC_GF2_ADD(s->d1[i], s->d2[(n1 - 1) * n1 + i]);
  }
  s->y = VEC_GF2_ADD(
      s->y,
      VEC_GF2_ADD(s->d1[n1 - 1],
                  s->d2[(n1 - 1) * n1 + ((n1 - 2) > (n1 - 1) ? 0 : (n1 - 2))]));

  return s;
}

state *part_eval(poly_vec_t *systems, uint8_t *prefix, unsigned int n,
                 unsigned int n1, poly_vec_t *parities, state *s)
{
  s = update(s, systems, n, n1, prefix);

  s = fes_eval_parity(systems, n, n1, prefix, s, parities);

  if (!s)
  {
    destroy_state(s);
    return NULL;
  }

  return s;
}

uint8_t fes_recover_vectorized(poly_t *system, poly_vec_t *e_k_systems,
                               unsigned int n, unsigned int n1, poly_vec_t deg,
                               poly_t *result)
{
  state *s = NULL;

  sub_poly_t max_deg = VEC_MAX(deg);

  uint8_t *prefix = calloc((n - n1), sizeof(uint8_t));
  if (!prefix) return 1;

  size_t d_size = 0;
  for (unsigned int i = 0; i <= max_deg; i++)
  {
    d_size += lk_binom[(n - n1) * BINOM_DIM2 + i];
  }

  poly_vec_t *d = aligned_alloc(ALIGNMENT, d_size * sizeof(poly_vec_t));
  if (!d) return 1;
  memset(d, 0, d_size * sizeof(poly_vec_t));

  unsigned int *alpha = calloc(max_deg, sizeof(unsigned int));
  if (!alpha) return 1;

  poly_vec_t parities = VEC_0;

  s = part_eval(e_k_systems, prefix, n, n1, &parities, s);

  if (!s)
  {
    free(alpha);
    free(d);
    free(prefix);
    return 1;
  }

  if (_avx_sol_overlap(parities))
  {
    poly_t solution, fixed_solution;
    PotentialSolution potential_solutions[(1 << FIXED_VARS) * 6] = {0};

    int amount = _avx_extract_sol(parities, potential_solutions);

    for (int count = 0; count < amount; count++)
    {
      fixed_solution = GF2_ADD(
          0,
          INT_LSHIFT(GF2_ADD(potential_solutions[count].solution, INT_MASK(n1)),
                     (n - n1)));

      solution = GF2_ADD(potential_solutions[count].fixed_var,
                         INT_LSHIFT(fixed_solution, FIXED_VARS));

      if (!eval(system, n + FIXED_VARS, solution))
      {
        *result = solution;

        destroy_state(s);
        free(prefix);
        free(alpha);
        free(d);

        return 0;
      }
    }
  }

  d[0] = parities;

  parities = VEC_0;

  poly_vec_t mask;

  for (unsigned int si = 1; si < (1u << (n - n1)); si++)
  {
    if (hamming_weight(si) > max_deg)
    {
      g_recover_eval++;
      BEGIN_BENCH(g_recover_eval_time)

      poly_t len_alpha = bits(si, alpha, max_deg);

      for (unsigned int j = len_alpha; j-- > 0;)
      {
        mask = ~VEC_GT(VEC_ASSIGN_ONE(j), deg);

        unsigned int idx =
            (j == 0) ? 0 : monomial_to_index(si, n - n1, alpha[j - 1]);

        poly_vec_t added =
            VEC_GF2_ADD(d[idx], d[monomial_to_index(si, n - n1, alpha[j])]);

        d[idx] = VEC_BLEND(d[idx], added, mask);
      }

      END_BENCH(g_recover_eval_time)
    }
    else
    {
      g_recover_interp++;
      BEGIN_BENCH(g_recover_interp_time)

      poly_t len_alpha = bits(si, alpha, max_deg);

      unsigned int gray_si = GRAY(si);

      for (unsigned int pos = 0; pos < (n - n1); pos++)
      {
        prefix[pos] = (1 & (gray_si >> pos));
      }

      BEGIN_BENCH(g_fes_time)

      s = part_eval(e_k_systems, prefix, n, n1, &parities, s);

      END_BENCH(g_fes_time)

      if (!s)
      {
        free(alpha);
        free(d);
        free(prefix);
        return 1;
      }

      poly_vec_t prev = d[0];
      d[0] = parities;

      parities = VEC_0;

      for (unsigned int j = 1; j <= len_alpha; j++)
      {
        mask = ~VEC_GT(VEC_ASSIGN_ONE(j), deg);

        poly_vec_t tmp = VEC_0;

        if (j < n)
        {
          tmp =
              VEC_BLEND(tmp, d[monomial_to_index(si, n - n1, alpha[j - 1])], mask);
        }

        unsigned int idx =
            (j == 1) ? 0 : monomial_to_index(si, n - n1, alpha[j - 2]);

        poly_vec_t added = VEC_GF2_ADD(d[idx], prev);

        d[monomial_to_index(si, n - n1, alpha[j - 1])] =
            VEC_BLEND(d[monomial_to_index(si, n - n1, alpha[j - 1])], added, mask);

        if (j < n)
        {
          prev = VEC_BLEND(prev, tmp, mask);
        }
      }

      END_BENCH(g_recover_interp_time)
    }

    if (_avx_sol_overlap(d[0]))
    {
      poly_t solution, fixed_solution;
      PotentialSolution potential_solutions[6 * (1 << FIXED_VARS)] = {0};

      int amount = _avx_extract_sol(d[0], potential_solutions);

      for (int count = 0; count < amount; count++)
      {
        fixed_solution = GF2_ADD(
            GRAY(si), INT_LSHIFT(GF2_ADD(potential_solutions[count].solution,
                                         INT_MASK(n1)),
                                 (n - n1)));

        solution = GF2_ADD(potential_solutions[count].fixed_var,
                           INT_LSHIFT(fixed_solution, FIXED_VARS));

        if (!eval(system, n + FIXED_VARS, solution))
        {
          *result = solution;

          destroy_state(s);
          free(prefix);
          free(alpha);
          free(d);

          return 0;
        }
      }
    }
  }
  destroy_state(s);
  free(prefix);
  free(alpha);
  free(d);

  return 2;
}