#ifndef VEC128_CONFIG_H
#define VEC128_CONFIG_H

#include <immintrin.h>

#include "common_utils.h"

#if defined(INT8)  ////////////////////////////////////// EPI8

#define VECTOR_SIZE 16

#define FIXED_VARS 2

#define INT_TYPE uint8_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu8, &p)
#endif

#define _avx_sll(i, w) _avx128_sll_8(i, w)
#define _avx_srl(i, w) _avx128_srl_8(i, w)
#define _avx_eq(a, b) _mm_cmpeq_epi8(a, b)
#define _avx_lt(a, b) _mm_cmplt_epi8(a, b)
#define _avx_gt(a, b) _mm_cmpgt_epi8(a, b)
#define _avx_set1(val) _mm_set1_epi8(val)
#define _avx_insert(a, val, idx) _avx_insert(a, val, idx)
#define _avx_extract(a, idx) _avx_extract(a, idx)
#define _avx_m_mask(a) _mm_movemask_epi8(a)
#define _avx_add(a, b) _mm_add_epi8(a, b)
#define _avx_sub(a, b) _mm_sub_epi8(a, b)
#define _avx_rotate(a, w) _avx_rotate(a, w)
#define _avx_broadcast(a) _mm_set1_epi32(_mm_extract_epi32(a, 0))

static inline __m128i _avx_insert(__m128i reg, uint8_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm_insert_epi8(reg, val, 0);
    case 1:
      return _mm_insert_epi8(reg, val, 1);
    case 2:
      return _mm_insert_epi8(reg, val, 2);
    case 3:
      return _mm_insert_epi8(reg, val, 3);
    case 4:
      return _mm_insert_epi8(reg, val, 4);
    case 5:
      return _mm_insert_epi8(reg, val, 5);
    case 6:
      return _mm_insert_epi8(reg, val, 6);
    case 7:
      return _mm_insert_epi8(reg, val, 7);
    case 8:
      return _mm_insert_epi8(reg, val, 8);
    case 9:
      return _mm_insert_epi8(reg, val, 9);
    case 10:
      return _mm_insert_epi8(reg, val, 10);
    case 11:
      return _mm_insert_epi8(reg, val, 11);
    case 12:
      return _mm_insert_epi8(reg, val, 12);
    case 13:
      return _mm_insert_epi8(reg, val, 13);
    case 14:
      return _mm_insert_epi8(reg, val, 14);
    case 15:
      return _mm_insert_epi8(reg, val, 15);
    default:
      return reg;
  }
}

static inline uint8_t _avx_extract(__m128i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm_extract_epi8(reg, 0);
    case 1:
      return _mm_extract_epi8(reg, 1);
    case 2:
      return _mm_extract_epi8(reg, 2);
    case 3:
      return _mm_extract_epi8(reg, 3);
    case 4:
      return _mm_extract_epi8(reg, 4);
    case 5:
      return _mm_extract_epi8(reg, 5);
    case 6:
      return _mm_extract_epi8(reg, 6);
    case 7:
      return _mm_extract_epi8(reg, 7);
    case 8:
      return _mm_extract_epi8(reg, 8);
    case 9:
      return _mm_extract_epi8(reg, 9);
    case 10:
      return _mm_extract_epi8(reg, 10);
    case 11:
      return _mm_extract_epi8(reg, 11);
    case 12:
      return _mm_extract_epi8(reg, 12);
    case 13:
      return _mm_extract_epi8(reg, 13);
    case 14:
      return _mm_extract_epi8(reg, 14);
    case 15:
      return _mm_extract_epi8(reg, 15);
    default:
      return 0;
  }
}

static inline __m128i _avx_rotate(__m128i reg, int width)
{
  width = width % 16;  // 128 bit integers
  switch (width)
  {
    case 0:
      return reg;
    case 1:
      return _mm_alignr_epi8(reg, reg, 1);
    case 2:
      return _mm_alignr_epi8(reg, reg, 2);
    case 3:
      return _mm_alignr_epi8(reg, reg, 3);
    case 4:
      return _mm_alignr_epi8(reg, reg, 4);
    case 5:
      return _mm_alignr_epi8(reg, reg, 5);
    case 6:
      return _mm_alignr_epi8(reg, reg, 6);
    case 7:
      return _mm_alignr_epi8(reg, reg, 7);
    case 8:
      return _mm_alignr_epi8(reg, reg, 8);
    case 9:
      return _mm_alignr_epi8(reg, reg, 9);
    case 10:
      return _mm_alignr_epi8(reg, reg, 10);
    case 11:
      return _mm_alignr_epi8(reg, reg, 11);
    case 12:
      return _mm_alignr_epi8(reg, reg, 12);
    case 13:
      return _mm_alignr_epi8(reg, reg, 13);
    case 14:
      return _mm_alignr_epi8(reg, reg, 14);
    case 15:
      return _mm_alignr_epi8(reg, reg, 15);
  }
  return _mm_setzero_si128();
}

static inline __m128i _avx128_sll_8(__m128i reg, int w)
{
  int bits = 8 - w;
  __m128i lo = _mm_and_si128(reg, _mm_set1_epi16(0x00FF));
  __m128i hi = _mm_and_si128(reg, _mm_set1_epi16(0xFF00));
  if (bits <= 0)
  {
    return _mm_setzero_si128();
  }
  else
  {
    __m128i mask = _mm_set1_epi16((1 << bits) - 1);
    return _mm_xor_si128(
        _mm_sll_epi16(hi, _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w)),
        _mm_sll_epi16(_mm_and_si128(lo, mask),
                      _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w)));
  }
}

static inline __m128i _avx128_srl_8(__m128i reg, int w)
{
  int bits = 8 - w;
  __m128i lo = _mm_and_si128(reg, _mm_set1_epi16(0x00FF));
  __m128i hi = _mm_and_si128(reg, _mm_set1_epi16(0xFF00));
  if (bits <= 0)
  {
    return _mm_setzero_si128();
  }
  else
  {
    __m128i mask = _mm_sll_epi16(_mm_set1_epi16((1 << bits) - 1),
                                 _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w + 8));

    return _mm_xor_si128(
        _mm_srl_epi16(_mm_and_si128(hi, mask),
                      _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w)),
        _mm_srl_epi16(lo, _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w)));
  }
}

static inline uint8_t _avx_max(__m128i a)
{
  __m128i tmp = a;
  __m128i perm, mask;

  for (int i = 0; i < 15; i++)
  {
    perm = _mm_alignr_epi8(tmp, tmp, 1);
    mask = _mm_cmpgt_epi8(tmp, perm);
    tmp = (__m128i)_mm_blendv_ps((__m128)perm, (__m128)tmp, (__m128)mask);
  }

  return _mm_extract_epi8(tmp, 0);
}

#elif defined(INT16)  ////////////////////////////////////// EPI16

#define VECTOR_SIZE 8

#define FIXED_VARS 1

#define INT_TYPE uint16_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu16, &p)
#endif

#define _avx_sll(i, w) _mm_sll_epi16(i, _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w))
#define _avx_srl(i, w) _mm_srl_epi16(i, _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, w))
#define _avx_eq _mm_cmpeq_epi16
#define _avx_lt _mm_cmplt_epi16
#define _avx_gt _mm_cmpgt_epi16
#define _avx_set1 _mm_set1_epi16
#define _avx_insert _avx128_insert_16
#define _avx_extract _avx128_extract_16
#define _avx_m_mask _avx128_movemask_16
#define _avx_add _mm_add_epi16
#define _avx_sub _mm_sub_epi16
#define _avx_rotate(a, w) _avx_rotate(a, w)
#define _avx_broadcast(a) _mm_set1_epi64x(_mm_extract_epi64(a, 0))

static inline __m128i _avx128_insert_16(__m128i reg, uint16_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm_insert_epi16(reg, val, 0);
    case 1:
      return _mm_insert_epi16(reg, val, 1);
    case 2:
      return _mm_insert_epi16(reg, val, 2);
    case 3:
      return _mm_insert_epi16(reg, val, 3);
    case 4:
      return _mm_insert_epi16(reg, val, 4);
    case 5:
      return _mm_insert_epi16(reg, val, 5);
    case 6:
      return _mm_insert_epi16(reg, val, 6);
    case 7:
      return _mm_insert_epi16(reg, val, 7);
    default:
      return reg;
  }
}

static inline uint16_t _avx128_extract_16(__m128i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm_extract_epi16(reg, 0);
    case 1:
      return _mm_extract_epi16(reg, 1);
    case 2:
      return _mm_extract_epi16(reg, 2);
    case 3:
      return _mm_extract_epi16(reg, 3);
    case 4:
      return _mm_extract_epi16(reg, 4);
    case 5:
      return _mm_extract_epi16(reg, 5);
    case 6:
      return _mm_extract_epi16(reg, 6);
    case 7:
      return _mm_extract_epi16(reg, 7);
    default:
      return 0;
  }
}

static inline __m128i _avx_rotate(__m128i reg, int width)
{
  width = width % 8;  // 128 bit integers
  switch (width)
  {
    case 0:
      return reg;
    case 1:
      return _mm_alignr_epi8(reg, reg, 2);
    case 2:
      return _mm_alignr_epi8(reg, reg, 2 * 2);
    case 3:
      return _mm_alignr_epi8(reg, reg, 2 * 3);
    case 4:
      return _mm_alignr_epi8(reg, reg, 2 * 4);
    case 5:
      return _mm_alignr_epi8(reg, reg, 2 * 5);
    case 6:
      return _mm_alignr_epi8(reg, reg, 2 * 6);
    case 7:
      return _mm_alignr_epi8(reg, reg, 2 * 7);
  }
  return _mm_setzero_si128();
}

static inline int _avx128_movemask_16(__m128i reg)
{
  int mask = _mm_movemask_epi8(reg);
  int mask16 = 0;
  for (int i = 0; i < 16; i += 2)
  {
    if ((mask >> (i + 1)) & 1) mask16 ^= (1 << (i / 2));
  }
  return mask16;
}

static inline uint16_t _avx_max(__m128i reg)
{
  __m128i tmp = reg;
  __m128i perm, mask;
  for (int i = 0; i < 7; i++)
  {
    perm = _mm_alignr_epi8(tmp, tmp, 2);
    mask = _mm_cmpgt_epi16(tmp, perm);
    tmp = _mm_blendv_epi8(perm, tmp, mask);
  }

  return _mm_extract_epi32(tmp, 0);
}

#elif defined(INT32)  ////////////////////////////////////// EPI32

#define VECTOR_SIZE 4

#define FIXED_VARS 0

#define INT_TYPE uint32_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif

#define _avx_sll(reg, w) _mm_sll_epi32(reg, _mm_set_epi32(0, 0, 0, w))
#define _avx_srl(reg, w) _mm_srl_epi32(reg, _mm_set_epi32(0, 0, 0, w))
#define _avx_eq _mm_cmpeq_epi32
#define _avx_lt _mm_cmplt_epi32
#define _avx_gt _mm_cmpgt_epi32
#define _avx_set1 _mm_set1_epi32
#define _avx_insert _avx128_insert_32
#define _avx_extract _avx128_extract_32
#define _avx_m_mask(a) _mm_movemask_ps((__m128)a)
#define _avx_add _mm_add_epi32
#define _avx_sub _mm_sub_epi32
#define _avx_rotate(a) _mm_alignr_epi8(a, a, 4)

static inline __m128i _avx128_insert_32(__m128i reg, uint32_t val, int index)
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

static inline uint32_t _avx128_extract_32(__m128i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm_extract_epi32(reg, 0);
    case 1:
      return _mm_extract_epi32(reg, 1);
    case 2:
      return _mm_extract_epi32(reg, 2);
    case 3:
      return _mm_extract_epi32(reg, 3);
    default:
      return 0;
  }
}

static inline uint32_t _avx_max(__m128i reg)
{
  __m128i tmp = reg;
  __m128i perm, mask;

  for (int i = 0; i < 4; i++)
  {
    perm = (__m128i)_mm_permute_ps((__m128)tmp, _MM_SHUFFLE(0, 3, 2, 1));
    mask = _mm_cmpgt_epi32(tmp, perm);
    tmp = (__m128i)_mm_blendv_ps((__m128)perm, (__m128)tmp, (__m128)mask);
  }

  return _mm_extract_epi32(tmp, 0);
}

#endif

#define _avx_xor(a, b) _mm_xor_si128(a, b)
#define _avx_and(a, b) _mm_and_si128(a, b)
#define _avx_or(a, b) _mm_or_si128(a, b)
#define _avx_blend(a, b, mask) _mm_blendv_epi8(a, b, mask)
#define _avx_testz(a, b) _mm_testz_si128(a, b)
#define _avx_zero() _mm_setzero_si128()

#endif  // !VEC128_CONFIG_H