#include "vector_utils.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>

#include "benchmark.h"
#include "utils.h"

sub_poly_t parity(poly_t bits)
{
  return __builtin_parity(bits);
}  // TODO: common

int _avx_sol_overlap(poly_vec_t reg)
{
  for (int i = 0; i < (1 << FIXED_VARS); i++)
  {
    poly_vec_t reg_alt = _avx_broadcast(_avx_rotate(reg, 4 * i));
    poly_vec_t perm = reg_alt, mask = _avx_zero();

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

int _avx_extract_sol(poly_vec_t reg, PotentialSolution *solutions)
{
  int sol_idx = 0;

  for (int i = 0; i < (1 << FIXED_VARS); i++)
  {
    poly_vec_t reg_alt = _avx_broadcast(_avx_rotate(reg, 4 * i));
    poly_vec_t perm = reg_alt;
    poly_vec_t mask_eq = _avx_zero();
    poly_vec_t mask_zbit =
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

sub_poly_t _avx_max(poly_vec_t a)
{
  poly_vec_t tmp = a;
  poly_vec_t perm, mask;

  for (int i = 0; i < (VECTOR_SIZE - 1); i++)
  {
    perm = _avx_rotate(tmp, 1);
    mask = _avx_gt(tmp, perm);
    tmp = (poly_vec_t)_avx_blend(perm, tmp, mask);
  }

  return _avx_extract(tmp, 0);
}

// TOOD: Add _avx_max here