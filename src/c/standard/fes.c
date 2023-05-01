#include "fes.h"

#include <stdio.h>
#include <time.h>

#include "benchmark.h"
#include "binom.h"

state *init_state(unsigned int n, unsigned int n1, uint8_t *prefix)
{
  state *s = malloc(sizeof(state));
  if (!s) return NULL;

  s->d1 = calloc(n1, sizeof(container_t));

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

  s->d2 = calloc(n1 * n1, sizeof(container_t));

  if (!(s->d2))
  {
    if (!prefix) free(s->prefix);

    free(s->d1);
    free(s);

    return NULL;
  }

  s->i = 0;
  s->y = 0;

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

unsigned int bit1(container_t i) { return trailing_zeros(i); }

unsigned int bit2(container_t i) { return bit1(GF2_ADD(i, INT_LSB(i))); }

// Assumes arr has been allocated with arr_len bits.
unsigned int bits(container_t i, unsigned int *arr, unsigned int arr_len)
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

    i = GF2_ADD(i, INT_LSB(i));
  }

  return sum;
}

unsigned int monomial_to_index(container_t mon, unsigned int n,
                               unsigned int boundary)
{
  unsigned int i;
  int d = 0;
  unsigned int index = 0;
  unsigned int index_d = 0;
  for (i = 0; i < n; i++)
    if (!INT_IS_ZERO(INT_IDX(mon, i)) && (i <= boundary))
    {
      d++;

      index += lk_binom[i * BINOM_DIM2 + d];
      index_d += lk_binom[n * BINOM_DIM2 + d];
    }
  index = index_d - index;

  return index;
}

state *init(state *s, container_t *system, unsigned int n, unsigned int n1,
            uint8_t *prefix)
{
  s = init_state(n, n1, prefix);

  if (!s)
  {
    return NULL;
  }

  s->y = system[0];

  for (unsigned int k = 0; k < n1; k++)
  {
    for (unsigned int j = 0; j < k; j++)
    {
      s->d2[k * n1 + j] = system[lex_idx(j + (n - n1), k + (n - n1), n)];
    }
  }

  s->d1[0] = system[1 + n - n1];

  for (unsigned int k = 1; k < n1; k++)
  {
    s->d1[k] = GF2_ADD(s->d2[k * n1 + (k - 1)], system[1 + k + (n - n1)]);
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (prefix[i] == 0) continue;

    for (unsigned int k = 0; k < n1; k++)
    {
      s->d1[k] = GF2_ADD(s->d1[k], system[lex_idx(i, k + (n - n1), n)]);
    }

    s->y = GF2_ADD(s->y, system[i + 1]);
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (prefix[i] == 0) continue;

    for (unsigned int j = i + 1; j < (n - n1); j++)
    {
      if (prefix[j] == 0) continue;

      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
    }
  }

  return s;
}

state *update(state *s, container_t *system, unsigned int n, unsigned int n1,
              uint8_t *prefix)
{
  if (!s)
  {
    return init(s, system, n, n1, prefix);
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

      s->d1[k] = GF2_ADD(s->d1[k], system[lex_idx(small_idx, large_idx, n)]);
    }

    s->y = GF2_ADD(s->y, system[idx + 1]);
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

      s->y = GF2_ADD(s->y, system[lex_idx(small_idx, large_idx, n)]);
    }
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (off[i] == 0) continue;

    for (unsigned int j = i + 1; j < (n - n1); j++)
    {
      if (off[j] == 0) continue;

      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
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

      s->d1[k] = GF2_ADD(s->d1[k], system[lex_idx(small_idx, large_idx, n)]);
    }

    s->y = GF2_ADD(s->y, system[idx + 1]);
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

      s->y = GF2_ADD(s->y, system[lex_idx(small_idx, large_idx, n)]);
    }
  }

  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (on[i] == 0) continue;

    for (unsigned int j = i + 1; j < (n - n1); j++)
    {
      if (on[j] == 0) continue;

      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
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

  unsigned int k1 = bit1(s->i);
  unsigned int k2 = bit2(s->i);

  if (k2 < (-1u))
  {
    s->d1[k1] = GF2_ADD(s->d1[k1], s->d2[k2 * n1 + k1]);
  }

  s->y = GF2_ADD(s->y, s->d1[k1]);
}

state *fes_eval_parity(container_t *system, unsigned int n, unsigned int n1,
                       uint8_t *prefix, state *s, container_t *parities)
{
  if (!s)
  {
    s = init(s, system, n, n1, prefix);

    if (!s)
    {
      return NULL;
    }
  }

  if (s->y == 0)
  {
    *parities = GF2_ADD(*parities, INT_MASK((n1 + 1)));
  }

  while (s->i < ((1 << n1) - 1))
  {
    step(s, n1);

    container_t z = GF2_ADD(s->i, INT_RSHIFT(s->i, 1));

    if (INT_IS_ZERO(s->y))
    {
      *parities = INT_SETBIT(*parities, 0, 1);

      for (unsigned int pos = 0; pos < n1; pos++)
      {
        if (INT_IS_ZERO(INT_IDX(z, pos)))
        {
          *parities = INT_SETBIT(*parities, (pos + 1), 1);
        }
      }
    }
  }

  s->i = 0;
  for (unsigned int i = 0; i < (n1 - 1); i++)
  {
    s->d1[i] = GF2_ADD(s->d1[i], s->d2[(n1 - 1) * n1 + i]);
  }

  s->y = GF2_ADD(
      s->y,
      GF2_ADD(s->d1[n1 - 1],
              s->d2[(n1 - 1) * n1 + ((n1 - 2) > (n1 - 1) ? 0 : (n1 - 2))]));

  return s;
}

state *fes_eval_solutions(container_t *system, unsigned int n, unsigned int n1,
                          uint8_t *prefix, state *s, container_t *solutions,
                          unsigned int *sol_amount)
{
  if (!s)
  {
    s = init(s, system, n, n1, prefix);

    if (!s)
    {
      return NULL;
    }
  }

  uint64_t pre_x = 0;
  for (unsigned int i = 0; i < (n - n1); i++)
  {
    if (prefix[i] == 0) continue;

    pre_x += (1 << i);
  }

  if (s->y == 0)
  {
    solutions[(*sol_amount)++] = ((s->i ^ (s->i >> 1)) << (n - n1) | pre_x);
  }

  while (s->i < ((1 << n1) - 1))
  {
    step(s, n1);

    if (s->y == 0)
    {
      solutions[(*sol_amount)++] = ((s->i ^ (s->i >> 1)) << (n - n1) | pre_x);
    }
  }

  s->i = 0;
  for (unsigned int i = 0; i < (n1 - 1); i++)
  {
    s->d1[i] = GF2_ADD(s->d1[i], s->d2[(n1 - 1) * n1 + i]);
  }

  s->y = GF2_ADD(
      s->y,
      GF2_ADD(s->d1[n1 - 1],
              s->d2[(n1 - 1) * n1 + ((n1 - 2) > (n1 - 1) ? 0 : (n1 - 2))]));

  return s;
}

state *part_eval(container_t *system, uint8_t *prefix, unsigned int n,
                 unsigned int n1, container_t *parities, state *s)
{
  s = update(s, system, n, n1, prefix);

  s = fes_eval_parity(system, n, n1, prefix, s, parities);

  if (!s)
  {
    destroy_state(s);
    return NULL;
  }

  return s;
}

uint8_t fes_recover(container_t *system, unsigned int n, unsigned int n1,
                    unsigned int deg, PotentialSolution *results,
                    size_t *res_size)
{
  state *s = NULL;

  uint8_t *prefix = calloc((n - n1), sizeof(uint8_t));
  if (!prefix) return 1;

  size_t d_size = 0;
  for (unsigned int i = 0; i <= deg; i++)
  {
    d_size += lk_binom[(n - n1) * BINOM_DIM2 + i];
  }
  container_t *d =
      calloc(d_size,
             sizeof(container_t));  // TODO: Find suitable datastructure for d
                                    // and initialize.
  if (!d) return 1;

  unsigned int *k = calloc(deg, sizeof(unsigned int));
  if (!k) return 1;

  container_t parities = INT_0;

  s = part_eval(system, prefix, n, n1, &parities, s);

  if (!s)
  {
    free(k);
    free(d);
    free(prefix);
    return 1;
  }

  if (!INT_IS_ZERO(INT_IDX(parities, 0)))
  {
    results[0].y_idx = 0;
    results[0].z_bits = INT_RSHIFT(GF2_ADD(parities, INT_MASK((n1 + 1))), 1);
    (*res_size)++;
  }
  d[0] = parities;

  parities = 0;

  for (unsigned int si = 1; si < (1u << (n - n1)); si++)
  {
    if (hamming_weight(si) > deg)
    {
      g_recover_eval++;

      BEGIN_BENCH(g_recover_eval_time)

      unsigned int len_k = bits(si, k, deg);

      for (unsigned int j = len_k; j-- > 0;)
      {
        unsigned int idx =
            (j == 0) ? 0 : monomial_to_index(si, n - n1, k[j - 1]);
        if (GRAY(si) == 2912)
          // printf("%u %u\n", d[idx], d[monomial_to_index(si, n - n1, k[j])]);
          d[idx] = GF2_ADD(d[idx], d[monomial_to_index(si, n - n1, k[j])]);
      }

      END_BENCH(g_recover_eval_time)
    }
    else
    {
      g_recover_interp++;

      BEGIN_BENCH(g_recover_interp_time)

      unsigned int len_k = bits(si, k, deg);
      unsigned int gray_si = si ^ (si >> 1);
      for (unsigned int pos = 0; pos < (n - n1); pos++)
      {
        prefix[pos] = (1 & (gray_si >> pos));
      }

      BEGIN_BENCH(g_fes_time)

      s = part_eval(system, prefix, n, n1, &parities, s);

      END_BENCH(g_fes_time)

      if (!s)
      {
        free(k);
        free(d);
        free(prefix);
        return 1;
      }

      unsigned int prev = d[0];
      d[0] = parities;

      parities = INT_0;

      for (unsigned int j = 1; j <= len_k; j++)
      {
        unsigned int tmp;

        if (j < n)
        {
          tmp = d[monomial_to_index(si, n - n1, k[j - 1])];
        }

        unsigned int idx =
            (j == 1) ? 0 : monomial_to_index(si, n - n1, k[j - 2]);

        d[monomial_to_index(si, n - n1, k[j - 1])] = GF2_ADD(d[idx], prev);

        if (j < n)
        {
          prev = tmp;
        }
      }
      END_BENCH(g_recover_interp_time)
    }
    if (!INT_IS_ZERO(INT_IDX(d[0], 0)))
    {
      results[*res_size].y_idx = si;
      results[*res_size].z_bits =
          INT_RSHIFT(GF2_ADD(d[0], INT_MASK((n1 + 1))), 1);
      (*res_size)++;
    }
  }
  destroy_state(s);
  free(prefix);
  free(k);
  free(d);

  return 0;
}

unsigned int bruteforce(container_t *system, unsigned int n, unsigned int n1,
                        unsigned int d, container_t *solutions)
{
  unsigned int sol_amount = 0;

  uint8_t *prefix = malloc(n - n1);

  if (!prefix)
  {
    free(solutions);
    return -1;
  }

  state *s = NULL;

  for (unsigned int i = 0; i < (1u << (n - n1)); i++)
  {
    if (hamming_weight(i) <= d)
    {
      for (unsigned int pos = 0; pos < (n - n1); pos++)
      {
        prefix[pos] = (1 & (i >> pos));
      }

      s = update(s, system, n, n1, prefix);
      s = fes_eval_solutions(system, n, n1, prefix, s, solutions, &sol_amount);
    }
  }

  free(prefix);
  destroy_state(s);

  return sol_amount;
}

unsigned int fes(container_t *system, unsigned int n, unsigned int m,
                 container_t *solutions)
{
  unsigned int sol_amount = 0;

  clock_t fes_time = clock();
  state *s =
      fes_eval_solutions(system, n, n, NULL, NULL, solutions, &sol_amount);
  size_t msec = (clock() - fes_time) * 1000 / CLOCKS_PER_SEC;
  printf("FES solve time: %zus, %zums\n", msec / 1000, msec % 1000);

  if (!s) return 0;

  return sol_amount;
}