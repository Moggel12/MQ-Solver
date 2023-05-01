#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>

#include "benchmark.h"
#include "common_utils.h"
#include "vector_utils.h"

unsigned int hamming_weight(unsigned int x)
{
  return __builtin_popcount(x);
}  // TODO: common

unsigned int trailing_zeros(unsigned int v)  // TODO: common
{
  unsigned int c;
  if (v)
  {
    v = (v ^ (v - 1)) >> 1;
    for (c = 0; v; c++)
    {
      v >>= 1;
    }
  }
  else
  {
    c = (-1u);
  }
  return c;
}

container_t parity(container_t bits)
{
  return __builtin_parity(bits);
}  // TODO: common

container_t gen_row(unsigned int m)
{
  return (rand() & ((1 << m) - 1));
}  // TODO: common

int lex_idx(unsigned int i, unsigned int j, unsigned int n)  // TODO: common
{
  int sum = 0;
  for (unsigned int k = 1; k < i + 2; k++)
  {
    sum += (n - k);
  }
  return n + sum - (n - j - 1);
}

int n_choose_k(int n, int k)  // TODO: common
{
  if (k > n) return 0;
  if (k == n) return 1;
  if (k > (n - k)) k = n - k;
  double c = 1;
  for (long i = 1; i <= k; i++)
  {
    c *= n--;
    c /= i;
  }
  return (int)c;
}

unsigned int gen_matrix(container_t *mat, unsigned int n_rows,
                        unsigned int n_columns)  // TODO: common
{
  container_t *mat_copy = malloc(n_rows * sizeof(container_t));
  if (!mat_copy) return 1;

  unsigned int rank;

  do
  {
    rank = 0;
    for (unsigned int i = 0; i < n_rows; i++)
    {
      mat_copy[i] = mat[i] = gen_row(n_columns);
    }
    for (unsigned int i = 0; i < n_rows; i++)
    {
      if (!INT_IS_ZERO(mat_copy[i]))
      {
        rank++;
        unsigned int pivot_elm = INT_LSB(mat_copy[i]);
        for (unsigned int j = i + 1; j < n_rows; j++)
        {
          if (!INT_IS_ZERO(GF2_MUL(mat_copy[j], pivot_elm)))
          {
            mat_copy[j] = GF2_ADD(mat_copy[j], mat_copy[i]);
          }
        }
      }
    }

  } while (rank != n_rows);

  free(mat_copy);

  return rank != n_rows;
}

container_t eval(container_t *system, size_t n,
                 container_t var_values)  // TODO: common
{
  BEGIN_BENCH(g_eval_time)

  container_t table[2] = {INT_0, INT_FF};

  container_t res = system[0];

  for (unsigned int i = 0; i < n; i++)
  {
    res = GF2_ADD(res, GF2_MUL(table[INT_IDX(var_values, i)], system[i + 1]));
  }

  int idx = lex_idx(0, 1, n);
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = i + 1; j < n; j++)
    {
      container_t monomial =
          GF2_MUL(table[INT_IDX(var_values, i)], table[INT_IDX(var_values, j)]);
      res = GF2_ADD(res, GF2_MUL(monomial, system[idx]));
      idx++;
    }
  }

  END_BENCH(g_eval_time)

  return res;
}

void print_register(__m128i a)  // TODO: Remove
{
  printf("%u ", _mm_extract_epi16(a, 0));
  printf("%u ", _mm_extract_epi16(a, 1));
  printf("%u ", _mm_extract_epi16(a, 2));
  printf("%u ", _mm_extract_epi16(a, 3));
  printf("%u ", _mm_extract_epi16(a, 4));
  printf("%u ", _mm_extract_epi16(a, 5));
  printf("%u ", _mm_extract_epi16(a, 6));
  printf("%u ", _mm_extract_epi16(a, 7));
  printf("\n");
}

int _avx_sol_overlap(container_vec_t reg)
{
  for (int i = 0; i < (1 << FIXED_VARS); i++)
  {
    container_vec_t reg_alt = _avx_broadcast(_avx_rotate(reg, 4 * i));
    container_vec_t perm = reg_alt, mask = _avx_zero();

    for (int j = 0; j < (VECTORIZED_ROUNDS - 1); j++)
    {
      perm = _mm_alignr_epi8(perm, perm, 2);
      mask = _avx_or(mask, _avx_eq(reg_alt, perm));
    }
    mask =
        _avx_and(mask, _avx_eq(_avx_and(reg_alt, _avx_set1(1)), _avx_set1(1)));

    if (!_avx_testz(mask, mask)) return 1;
  }

  return 0;
}

int _avx_extract_sol(container_vec_t reg, PotentialSolution *solutions)
{
  int sol_idx = 0;

  for (int i = 0; i < (1 << FIXED_VARS); i++)
  {
    container_vec_t reg_alt = _avx_broadcast(_avx_rotate(reg, 4 * i));
    container_vec_t perm = reg_alt;
    container_vec_t mask_eq = _avx_zero();
    container_vec_t mask_zbit =
        _avx_eq(_avx_and(reg_alt, _avx_set1(1)), _avx_set1(1));

    for (int j = 0; j < (VECTORIZED_ROUNDS - 1); j++)
    {
      perm = _mm_alignr_epi8(perm, perm, 2);
      mask_eq = _avx_and(mask_zbit, _mm_cmpeq_epi16(reg_alt, perm));

      if (!_avx_testz(mask_eq, mask_eq))
      {
        unsigned int m_mask1 = _avx_m_mask(mask_eq) & 0xF;
        unsigned int m_mask2 = m_mask1 & ~1;
        unsigned int idx = trailing_zeros(m_mask1);

        solutions[sol_idx].fixed_var = i;
        solutions[sol_idx++].solution = _avx_extract(reg_alt, idx) >> 1;

        if (m_mask2 > 0)
        {
          solutions[sol_idx].fixed_var = i;
          solutions[sol_idx++].solution =
              _avx_extract(reg_alt, trailing_zeros(m_mask2)) >> 1;
        }
      }
    }
  }

  for (int i = 0; i < sol_idx; i++)
  {
    for (int j = i + 1; j < sol_idx; j++)
    {
      if ((solutions[i].solution == solutions[j].solution) &&
          (solutions[i].fixed_var == solutions[j].fixed_var))
      {
        for (int k = j; k < sol_idx - 1; k++)
        {
          solutions[k] = solutions[k + 1];
        }
        sol_idx--;

        j--;
      }
    }
  }

  return sol_idx;
}

// TOOD: Add _avx_max here