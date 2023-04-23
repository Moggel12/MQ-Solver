#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "benchmark.h"
#include "binom.h"
#include "fes.h"

state *init_state(unsigned int n, unsigned int n1, uint8_t *prefix)
{
  state *s = malloc(sizeof(state));
  if (!s) return NULL;

  s->d1 = aligned_alloc(ALIGNMENT, n1 * sizeof(container_vec_t));

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

  s->d2 = aligned_alloc(ALIGNMENT, n1 * n1 * sizeof(container_vec_t));

  if (!(s->d2))
  {
    if (!prefix) free(s->prefix);

    free(s->d1);
    free(s);

    return NULL;
  }

  memset(s->d1, 0, n1 * sizeof(container_vec_t));
  memset(s->d2, 0, n1 * n1 * sizeof(container_vec_t));

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

unsigned int bit1(container_t i) { return trailing_zeros(i); }

unsigned int bit2(container_t i) { return bit1(i ^ (i & -i)); }

// Assumes arr has been allocated with arr_len bits.
container_t bits(container_t i, unsigned int *arr, container_t arr_len)
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

// TODO: Check if vector is needed here
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

state *init(state *s, container_vec_t *systems, unsigned int n, unsigned int n1,
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

state *update(state *s, container_vec_t *systems, unsigned int n,
              unsigned int n1, uint8_t *prefix)
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

  unsigned int k1 = bit1(s->i);
  unsigned int k2 = bit2(s->i);

  if (k2 < (-1u))
  {
    s->d1[k1] = VEC_GF2_ADD(s->d1[k1], s->d2[k2 * n1 + k1]);
  }

  s->y = VEC_GF2_ADD(s->y, s->d1[k1]);
}

state *fes_eval_parity(container_vec_t *systems, unsigned int n,
                       unsigned int n1, uint8_t *prefix, state *s,
                       container_vec_t *parities)
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

  container_vec_t zero_mask = VEC_IS_ZERO(s->y);
  // if (pre_x == 209)
  // {
  //   printf("# Constant\n");
  //   print_register(*parities);
  //   print_register(zero_mask);
  // }
  container_vec_t added = VEC_GF2_ADD(*parities, VEC_MASK((n1 + 1)));

  *parities = VEC_BLEND(*parities, added, zero_mask);
  // if (pre_x == 209) print_register(*parities);
  // if (pre_x == 209) printf("#\n");

  while (s->i < ((1u << n1) - 1))
  {
    step(s, n1);

    container_t z = GRAY(s->i);

    zero_mask = VEC_IS_ZERO(s->y);

    added = VEC_SETBIT(*parities, 0, 1);
    // if (pre_x == 209)
    // {
    //   printf("# U0\n");
    //   print_register(*parities);
    //   print_register(zero_mask);
    // }
    *parities = VEC_BLEND(*parities, added, zero_mask);
    // if (pre_x == 209) print_register(*parities);
    // if (pre_x == 209) printf("#\n");

    for (unsigned int pos = 0; pos < n1; pos++)
    {
      container_vec_t z_mask =
          VEC_AND(zero_mask, VEC_ASSIGN_ONE(-INT_IS_ZERO(INT_IDX(z, pos))));

      added = VEC_SETBIT(*parities, (pos + 1), 1);
      // if (pre_x == 209)
      // {
      //   printf("# U_%u\n", pos + 1);
      //   printf("%u\n", z);
      //   print_register(*parities);
      //   print_register(z_mask);
      // }
      *parities = VEC_BLEND(*parities, added, z_mask);
      // if (pre_x == 209) print_register(*parities);
      // if (pre_x == 209) printf("#\n");
    }
  }

  s->i = 0;
  for (unsigned int i = 0; i < (n1 - 1); i++)
  {
    s->d1[i] = VEC_GF2_ADD(s->d1[i], s->d2[(n1 - 1) * n1 + i]);
  }
  // if (pre_x == 209) printf("--\n");
  s->y = VEC_GF2_ADD(
      s->y,
      VEC_GF2_ADD(s->d1[n1 - 1],
                  s->d2[(n1 - 1) * n1 + ((n1 - 2) > (n1 - 1) ? 0 : (n1 - 2))]));

  return s;
}

state *part_eval(container_vec_t *systems, uint8_t *prefix, unsigned int n,
                 unsigned int n1, container_vec_t *parities, state *s)
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

uint8_t fes_recover_vectorized(container_t *system,
                               container_vec_t *e_k_systems, unsigned int n,
                               unsigned int n1, container_vec_t deg,
                               container_t *result)
{
  state *s = NULL;
  // print_register(deg);
  container_t max_deg = VEC_MAX(deg);

  uint8_t *prefix = calloc((n - n1), sizeof(uint8_t));
  if (!prefix) return 1;

  size_t d_size = 0;
  for (unsigned int i = 0; i <= max_deg; i++)
  {
    d_size += lk_binom[(n - n1) * BINOM_DIM2 + i];
  }

  container_vec_t *d =
      aligned_alloc(ALIGNMENT, d_size * sizeof(container_vec_t));
  if (!d) return 1;
  memset(d, 0, d_size * sizeof(container_vec_t));

  unsigned int *k = calloc(max_deg, sizeof(unsigned int));
  if (!k) return 1;

  container_vec_t parities = VEC_0;

  s = part_eval(e_k_systems, prefix, n, n1, &parities, s);

  if (!s)
  {
    free(k);
    free(d);
    free(prefix);
    return 1;
  }

  if (VEC_SOL_OVERLAP(parities))
  {
    container_t solution = GF2_ADD(
        0,
        INT_LSHIFT(GF2_ADD(VEC_EXTRACT_SOL(parities), INT_MASK(n1)), (n - n1)));

    if (!eval(system, n, solution))
    {
      *result = solution;

      destroy_state(s);
      free(prefix);
      free(k);
      free(d);

      return 0;
    }
  }
  // results[0] = parities;
  d[0] = parities;

  parities = VEC_0;

  container_vec_t mask;

  for (unsigned int si = 1; si < (1u << (n - n1)); si++)
  {
    if (hamming_weight(si) > max_deg)
    {
      container_t len_k = bits(si, k, max_deg);

      // if (GRAY(si) == 2912)
      // {
      //   for (int i = 0; i < len_k; i++)
      //   {
      //     printf("%u ", k[i]);
      //   }
      //   printf("\n");
      // }

      for (unsigned int j = len_k; j-- > 0;)
      {
        mask = ~VEC_GT(VEC_ASSIGN_ONE(j), deg);
        // if (GRAY(si) == 2912) printf("%u: ", j);
        // if (GRAY(si) == 2912) print_register(mask);

        unsigned int idx =
            (j == 0) ? 0 : monomial_to_index(si, n - n1, k[j - 1]);
        // if (GRAY(si) == 2912) printf("VAL:\n");
        // if (GRAY(si) == 2912)
        //   printf("IDX: %u %u\n", idx, monomial_to_index(si, n - n1, k[j]));
        // if (GRAY(si) == 2912) print_register(d[idx]);
        // if (GRAY(si) == 2912)
        // print_register(d[monomial_to_index(si, n - n1, k[j])]);
        // if (GRAY(si) == 2912) printf("\n");
        // printf("%u %u\n", idx, monomial_to_index(si, n - n1, k[j]));
        // printf("%u\n", j);
        // print_register(d[idx]);
        // print_register(d[monomial_to_index(si, n - n1, k[j])]);

        container_vec_t added =
            VEC_GF2_ADD(d[idx], d[monomial_to_index(si, n - n1, k[j])]);

        d[idx] = VEC_BLEND(d[idx], added, mask);
        // if (GRAY(si) == 209) print_register(d[idx]);
      }
      // printf("===\n");

      // if (GRAY(si) == 2912) printf("TEST\n");
    }
    else
    {
      container_t len_k = bits(si, k, max_deg);

      unsigned int gray_si = GRAY(si);

      for (unsigned int pos = 0; pos < (n - n1); pos++)
      {
        prefix[pos] = (1 & (gray_si >> pos));
      }

      s = part_eval(e_k_systems, prefix, n, n1, &parities, s);

      if (!s)
      {
        free(k);
        free(d);
        free(prefix);
        return 1;
      }

      container_vec_t prev = d[0];
      d[0] = parities;

      parities = VEC_0;
      // print_register(d[0]);

      for (unsigned int j = 1; j <= len_k; j++)
      {
        mask = ~VEC_GT(VEC_ASSIGN_ONE(j), deg);

        container_vec_t tmp;

        if (j < n)
        {
          // printf("Setting tmp\n");
          // tmp = prev;
          tmp =
              VEC_BLEND(tmp, d[monomial_to_index(si, n - n1, k[j - 1])], mask);
        }

        unsigned int idx =
            (j == 1) ? 0 : monomial_to_index(si, n - n1, k[j - 2]);
        // printf("%u\n", j);
        // print_register(d[idx]);
        // print_register(prev);
        container_vec_t added = VEC_GF2_ADD(d[idx], prev);

        d[monomial_to_index(si, n - n1, k[j - 1])] =
            VEC_BLEND(d[monomial_to_index(si, n - n1, k[j - 1])], added, mask);

        if (j < n)
        {
          // printf("Setting prev\n");
          // prev = tmp;
          prev = VEC_BLEND(prev, tmp, mask);
        }
      }
      // printf("===\n");
    }
    // if (GRAY(si) == 2912) print_register(d[0]);
    if (VEC_SOL_OVERLAP(d[0]))
    {
      // if (GRAY(si) == 2438) print_register(d[0]);
      // if (GRAY(si) == 2438) printf("%u\n", VEC_EXTRACT_SOL(d[0]));
      // if (GRAY(si) == 2438)
      //   printf("2438 overlaps! (%u)\n",
      //          GF2_ADD(VEC_EXTRACT_SOL(d[0]), INT_MASK(n1)));
      // if (GRAY(si) == 19) printf("Solutions overlap!\n");
      container_t solution = GF2_ADD(
          GRAY(si),
          INT_LSHIFT(GF2_ADD(VEC_EXTRACT_SOL(d[0]), INT_MASK(n1)), (n - n1)));
      if (!eval(system, n, solution))
      {
        *result = solution;

        destroy_state(s);
        free(prefix);
        free(k);
        free(d);

        return 0;
      }
    }
    // if (GRAY(si) == 19) printf("===\n");
  }
  destroy_state(s);
  free(prefix);
  free(k);
  free(d);

  return 2;
}