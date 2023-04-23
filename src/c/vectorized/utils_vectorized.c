#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>

#include "benchmark.h"
#include "utils.h"

unsigned int hamming_weight(unsigned int x) { return __builtin_popcount(x); }

unsigned int trailing_zeros(unsigned int v)
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

container_vec_t parity(container_vec_t bits)
{
  return VEC_ASSIGN(__builtin_parity(VEC_EXTRACT(bits, 3)),
                    __builtin_parity(VEC_EXTRACT(bits, 2)),
                    __builtin_parity(VEC_EXTRACT(bits, 1)),
                    __builtin_parity(VEC_EXTRACT(bits, 0)));
}

container_t gen_row(unsigned int m) { return (rand() & ((1 << m) - 1)); }

int lex_idx(unsigned int i, unsigned int j, unsigned int n)
{
  int sum = 0;
  for (unsigned int k = 1; k < i + 2; k++)
  {
    sum += (n - k);
  }
  return n + sum - (n - j - 1);
}

int n_choose_k(int n, int k)
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

uint8_t gen_matrix(container_vec_t *mat, unsigned int n_rows,
                   unsigned int n_columns)
{
  container_t *mat_copy = malloc(n_rows * sizeof(container_t));
  if (!mat_copy) return 1;

  unsigned int rank;

  for (int idx = 4; idx-- > 0;)
  {
    do
    {
      rank = 0;
      for (unsigned int i = 0; i < n_rows; i++)
      {
        mat_copy[i] = gen_row(n_columns);
        mat[i] = VEC_INSERT(mat[i], mat_copy[i], idx);
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
  }

  free(mat_copy);

  return 0;
}

container_t eval(container_t *system, size_t n, container_t var_values)
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

#if defined(REG128)
///////////////// AVX 128

// TODO: Remove
void print_register(__m128i a)
{
  printf("%u ", _mm_extract_epi32(a, 0));
  printf("%u ", _mm_extract_epi32(a, 1));
  printf("%u ", _mm_extract_epi32(a, 2));
  printf("%u ", _mm_extract_epi32(a, 3));
  printf("\n");
}

__m128i _avx128_insert(__m128i reg, uint32_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm_insert_epi32(reg, val, 0);
    case 1:
      return _mm_insert_epi32(reg, val, 1);
    case 2:
      return _mm_insert_epi32(reg, val, 2);
    case 3:
      return _mm_insert_epi32(reg, val, 3);
    default:
      return _mm_setzero_si128();
  }
}

uint32_t _avx128_max(__m128i a)
{
  __m128i tmp = a;
  __m128i perm, mask;

  for (int i = 0; i < 4; i++)
  {
    perm = (__m128i)_mm_permute_ps((__m128)tmp, _MM_SHUFFLE(0, 3, 2, 1));
    mask = _mm_cmpgt_epi32(tmp, perm);
    tmp = (__m128i)_mm_blendv_ps((__m128)perm, (__m128)tmp, (__m128)mask);
  }

  return _mm_extract_epi32(tmp, 0);
}

int _avx128_sol_overlap(__m128i a)
{
  __m128i perm = a, mask = _mm_setzero_si128();

  for (int i = 0; i < 3; i++)
  {
    perm = (__m128i)_mm_permute_ps((__m128)perm, _MM_SHUFFLE(0, 3, 2, 1));
    mask = _mm_or_si128(mask, _mm_cmpeq_epi32(a, perm));
  }
  mask = _mm_and_si128(
      mask,
      _mm_cmpeq_epi32(_mm_and_si128(a, _mm_set1_epi32(1)), _mm_set1_epi32(1)));

  return !_mm_testz_ps((__m128)mask, (__m128)mask);
}

uint32_t _avx128_extract_sol(__m128i a)
{
  __m128i perm = a;
  __m128i mask_eq = _mm_setzero_si128();
  __m128i mask_zbit =
      _mm_cmpeq_epi32(_mm_and_si128(a, _mm_set1_epi32(1)), _mm_set1_epi32(1));

  for (int i = 0; i < 3; i++)
  {
    perm = (__m128i)_mm_permute_ps((__m128)perm, _MM_SHUFFLE(0, 3, 2, 1));
    mask_eq = _mm_and_si128(mask_zbit, _mm_cmpeq_epi32(a, perm));

    if (!_mm_testz_ps((__m128)mask_eq, (__m128)mask_eq))
    {
      break;
    }
  }

  int m_mask = _mm_movemask_ps((__m128)mask_eq);
  int idx = trailing_zeros(m_mask);

  switch (idx)
  {
    case 0:
      return _mm_extract_epi32(a, 0) >> 1;
    case 1:
      return _mm_extract_epi32(a, 1) >> 1;
    case 2:
      return _mm_extract_epi32(a, 2) >> 1;
    case 3:
      return _mm_extract_epi32(a, 3) >> 1;
  }

  return 0;
}

#elif defined(REG256)
//////////////// AVX 256

// TODO: Remove
void print_register(__m256i a)
{
  printf("%u ", _mm256_extract_epi64(a, 0));
  printf("%u ", _mm256_extract_epi64(a, 1));
  printf("%u ", _mm256_extract_epi64(a, 2));
  printf("%u ", _mm256_extract_epi64(a, 3));
  printf("\n");
}

__m256i _avx256_insert(__m256i reg, uint64_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm256_insert_epi64(reg, val, 0);
    case 1:
      return _mm256_insert_epi64(reg, val, 1);
    case 2:
      return _mm256_insert_epi64(reg, val, 2);
    case 3:
      return _mm256_insert_epi64(reg, val, 3);
    default:
      return _mm256_setzero_si256();
  }
}

uint64_t _avx256_max(__m256i a)
{
  __m256i tmp = a;
  __m256i perm, mask;

  for (int i = 0; i < 4; i++)
  {
    perm = (__m256i)_mm_permute_pd((__m256d)tmp, _MM_SHUFFLE(0, 3, 2, 1));
    mask = _mm256_cmpgt_epi64(tmp, perm);
    tmp = (__m256i)_mm256_blendv_pd((__m256d)perm, (__m256d)tmp, (__m256d)mask);
  }

  return _mm256_extract_epi64(tmp, 0);
}

int _avx256_sol_overlap(__m128i a)
{
  __m256i perm = a;
  __m256d = mask = _mm256_setzero_si256();

  for (int i = 0; i < 3; i++)
  {
    perm = _mm256_permute_pd((__m256d)perm, _MM_SHUFFLE(0, 3, 2, 1));
    mask = _mm256_or_si256(mask, _mm256_cmpeq_epi64(a, perm));
  }
  mask = _mm256_and_si256(
      mask, _mm256_cmpeq_epi64(_mm256_and_si256(a, _mm256_set1_epi64x(1)),
                               _mm256_set1_epi64x(1)));

  return !_mm256_testz_pd((__m256d)mask, (__m256d)mask);
}

uint64_t _avx256_extract_sol(__m256i a)
{
  __m256i perm = a;
  __m256i mask_eq;
  __m256i mask_zbit = _mm256_cmpeq_epi64(
      _mm256_and_si256(a, _mm256_set1_epi64x(1)), _mm256_set1_epi64x(1));

  int found = 0;

  for (int i = 0; i < 3; i++)
  {
    perm = (__m256i)_mm256_permute_pd((__m256)perm, _MM_SHUFFLE(0, 3, 2, 1));
    mask_eq = _mm256_and_si256(mask_zbit, _mm256_cmpeq_epi64(a, perm));

    if (!_mm256_testz_pd((__m256d)mask_eq, (__m256d)mask_eq))
    {
      found = 1;
      break;
    }
  }

  int m_mask = _mm256_movemask_pd((__m256d)mask_eq);
  int idx = trailing_zeros(m_mask);

  switch (idx)
  {
    case 0:
      return _mm256_extract_epi64(a, 0) >> 1;
    case 1:
      return _mm256_extract_epi64(a, 1) >> 1;
    case 2:
      return _mm256_extract_epi64(a, 2) >> 1;
    case 3:
      return _mm256_extract_epi64(a, 3) >> 1;
  }

  return 0;
}
#endif