#ifndef VEC256_CONFIG_H
#define VEC256_CONFIG_H

#include <immintrin.h>

// TODO: Add support for computing parity across a vector (both 256 and 128, no
// matter which integer size).

#define _avx_xor(a, b) _mm256_xor_si256(a, b)
#define _avx_and(a, b) _mm256_and_si256(a, b)
#define _avx_or(a, b) _mm256_or_si256()
#define _avx_blend(a, b, mask) _mm256_blendv_epi8(a, b, mask)
#define _avx_testz(a, b) _mm256_testz_si256(a, b)
#define _avx_zero _mm256_setzero_si256()

// TODO: Fix naming and stuff.

#if defined(INT8)  ////////////////////////////////////// EPI8

#define VECTOR_SIZE 32

#define SUB_POLY_TYPE uint8_t

#define FIXED_VARS 3

#if defined(_DEBUG)
#define _DEBUG_READ_SUB(p) scanf("%" SCNu8, &p)
#endif

#define _avx_sll(i, w) _avx256_sll_8(i, w)
#define _avx_srl(i, w) _avx256_srl_8(i, w)
#define _avx_eq _mm256_cmpeq_epi8
#define _avx_lt _mm256_cmplt_epi8
#define _avx_gt _mm256_cmpgt_epi8
#define _avx_set1 _mm256_set1_epi8
#define _avx_insert _avx256_insert_8
#define _avx_extract _avx256_extract_8
#define _avx_m_mask _mm256_movemask_epi8
#define _avx_add _mm256_add_epi8
#define _avx_sub _mm256_sub_epi8
#define _avx_rotate(a, w) _avx_rotate(a, w)
#define _avx_broadcast(a) _mm256_broadcastd_epi32(a)

static inline __m256i _avx256_insert_8(__m256i reg, uint8_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm256_insert_epi8(reg, val, 0);
    case 1:
      return _mm256_insert_epi8(reg, val, 1);
    case 2:
      return _mm256_insert_epi8(reg, val, 2);
    case 3:
      return _mm256_insert_epi8(reg, val, 3);
    case 4:
      return _mm256_insert_epi8(reg, val, 4);
    case 5:
      return _mm256_insert_epi8(reg, val, 5);
    case 6:
      return _mm256_insert_epi8(reg, val, 6);
    case 7:
      return _mm256_insert_epi8(reg, val, 7);
    case 8:
      return _mm256_insert_epi8(reg, val, 8);
    case 9:
      return _mm256_insert_epi8(reg, val, 9);
    case 10:
      return _mm256_insert_epi8(reg, val, 10);
    case 11:
      return _mm256_insert_epi8(reg, val, 11);
    case 12:
      return _mm256_insert_epi8(reg, val, 12);
    case 13:
      return _mm256_insert_epi8(reg, val, 13);
    case 14:
      return _mm256_insert_epi8(reg, val, 14);
    case 15:
      return _mm256_insert_epi8(reg, val, 15);
    case 16:
      return _mm256_insert_epi8(reg, val, 16);
    case 17:
      return _mm256_insert_epi8(reg, val, 17);
    case 18:
      return _mm256_insert_epi8(reg, val, 18);
    case 19:
      return _mm256_insert_epi8(reg, val, 19);
    case 20:
      return _mm256_insert_epi8(reg, val, 20);
    case 21:
      return _mm256_insert_epi8(reg, val, 21);
    case 22:
      return _mm256_insert_epi8(reg, val, 22);
    case 23:
      return _mm256_insert_epi8(reg, val, 23);
    case 24:
      return _mm256_insert_epi8(reg, val, 24);
    case 25:
      return _mm256_insert_epi8(reg, val, 25);
    case 26:
      return _mm256_insert_epi8(reg, val, 26);
    case 27:
      return _mm256_insert_epi8(reg, val, 27);
    case 28:
      return _mm256_insert_epi8(reg, val, 28);
    case 29:
      return _mm256_insert_epi8(reg, val, 29);
    case 30:
      return _mm256_insert_epi8(reg, val, 30);
    case 31:
      return _mm256_insert_epi8(reg, val, 31);
    default:
      return reg;
  }
}

static inline uint8_t _avx256_extract_8(__m256i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm256_extract_epi8(reg, 0);
    case 1:
      return _mm256_extract_epi8(reg, 1);
    case 2:
      return _mm256_extract_epi8(reg, 2);
    case 3:
      return _mm256_extract_epi8(reg, 3);
    case 4:
      return _mm256_extract_epi8(reg, 4);
    case 5:
      return _mm256_extract_epi8(reg, 5);
    case 6:
      return _mm256_extract_epi8(reg, 6);
    case 7:
      return _mm256_extract_epi8(reg, 7);
    case 8:
      return _mm256_extract_epi8(reg, 8);
    case 9:
      return _mm256_extract_epi8(reg, 9);
    case 10:
      return _mm256_extract_epi8(reg, 10);
    case 11:
      return _mm256_extract_epi8(reg, 11);
    case 12:
      return _mm256_extract_epi8(reg, 12);
    case 13:
      return _mm256_extract_epi8(reg, 13);
    case 14:
      return _mm256_extract_epi8(reg, 14);
    case 15:
      return _mm256_extract_epi8(reg, 15);
    case 16:
      return _mm256_extract_epi8(reg, 16);
    case 17:
      return _mm256_extract_epi8(reg, 17);
    case 18:
      return _mm256_extract_epi8(reg, 18);
    case 19:
      return _mm256_extract_epi8(reg, 19);
    case 20:
      return _mm256_extract_epi8(reg, 20);
    case 21:
      return _mm256_extract_epi8(reg, 21);
    case 22:
      return _mm256_extract_epi8(reg, 22);
    case 23:
      return _mm256_extract_epi8(reg, 23);
    case 24:
      return _mm256_extract_epi8(reg, 24);
    case 25:
      return _mm256_extract_epi8(reg, 25);
    case 26:
      return _mm256_extract_epi8(reg, 26);
    case 27:
      return _mm256_extract_epi8(reg, 27);
    case 28:
      return _mm256_extract_epi8(reg, 28);
    case 29:
      return _mm256_extract_epi8(reg, 29);
    case 30:
      return _mm256_extract_epi8(reg, 30);
    case 31:
      return _mm256_extract_epi8(reg, 31);
    default:
      return 0;
  }
}

static inline __m256i _avx_rotate(__m256i reg, int width)
{
  width = width % 32;
  switch (width)
  {
    case 0:
      return _mm256_alignr_epi8(reg, reg, 0);
    case 1:
      return _mm256_alignr_epi8(reg, reg, 1);
    case 2:
      return _mm256_alignr_epi8(reg, reg, 2);
    case 3:
      return _mm256_alignr_epi8(reg, reg, 3);
    case 4:
      return _mm256_alignr_epi8(reg, reg, 4);
    case 5:
      return _mm256_alignr_epi8(reg, reg, 5);
    case 6:
      return _mm256_alignr_epi8(reg, reg, 6);
    case 7:
      return _mm256_alignr_epi8(reg, reg, 7);
    case 8:
      return _mm256_alignr_epi8(reg, reg, 8);
    case 9:
      return _mm256_alignr_epi8(reg, reg, 9);
    case 10:
      return _mm256_alignr_epi8(reg, reg, 10);
    case 11:
      return _mm256_alignr_epi8(reg, reg, 11);
    case 12:
      return _mm256_alignr_epi8(reg, reg, 12);
    case 13:
      return _mm256_alignr_epi8(reg, reg, 13);
    case 14:
      return _mm256_alignr_epi8(reg, reg, 14);
    case 15:
      return _mm256_alignr_epi8(reg, reg, 15);
    case 16:
      return _mm256_alignr_epi8(reg, reg, 16);
    case 17:
      return _mm256_alignr_epi8(reg, reg, 17);
    case 18:
      return _mm256_alignr_epi8(reg, reg, 18);
    case 19:
      return _mm256_alignr_epi8(reg, reg, 19);
    case 20:
      return _mm256_alignr_epi8(reg, reg, 20);
    case 21:
      return _mm256_alignr_epi8(reg, reg, 21);
    case 22:
      return _mm256_alignr_epi8(reg, reg, 22);
    case 23:
      return _mm256_alignr_epi8(reg, reg, 23);
    case 24:
      return _mm256_alignr_epi8(reg, reg, 24);
    case 25:
      return _mm256_alignr_epi8(reg, reg, 25);
    case 26:
      return _mm256_alignr_epi8(reg, reg, 26);
    case 27:
      return _mm256_alignr_epi8(reg, reg, 27);
    case 28:
      return _mm256_alignr_epi8(reg, reg, 28);
    case 29:
      return _mm256_alignr_epi8(reg, reg, 29);
    case 30:
      return _mm256_alignr_epi8(reg, reg, 30);
    case 31:
      return _mm256_alignr_epi8(reg, reg, 31);
  }
  return reg;
}

static inline __m256i _avx256_sll_8(__m256i reg, int w)
{
  int bits = 8 - w;
  __m256i lo = _mm256_and_si256(reg, _mm256_set1_epi16(0x00FF));
  __m256i hi = _mm256_and_si256(reg, _mm256_set1_epi16(0xFF00));
  if (bits <= 0)
  {
    return _mm256_setzero_si256();
  }
  else
  {
    __m256i mask = _mm256_set1_epi16((1 << bits) - 1);
    return _mm256_xor_si256(
        _mm256_sll_epi16(hi, _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w)),
        _mm256_sll_epi16(_mm256_and_si256(lo, mask),
                         _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w)));
  }
}

static inline __m256i _avx256_srl_8(__m256i reg, int w)
{
  int bits = 8 - w;
  __m256i lo = _mm256_and_si256(reg, _mm256_set1_epi16(0x00FF));
  __m256i hi = _mm256_and_si256(reg, _mm256_set1_epi16(0xFF00));
  if (bits <= 0)
  {
    return _mm256_setzero_si256();
  }
  else
  {
    __m256i mask =
        _mm256_sll_epi16(_mm256_set1_epi16((1 << bits) - 1),
                         _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w + 8));

    return _mm256_xor_si256(
        _mm256_srl_epi16(_mm256_and_si256(hi, mask),
                         _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w)),
        _mm256_srl_epi16(lo, _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w)));
  }
}

#elif defined(INT16)  ////////////////////////////////////// EPI16

#define VECTOR_SIZE 16

#define SUB_POLY_TYPE uint16_t

#define FIXED_VARS 2

#if defined(_DEBUG)
#define _DEBUG_READ_SUB(p) scanf("%" SCNu16, &p)
#endif

#define _avx_sll(i, w) \
  _mm256_sll_epi16(i, _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w))
#define _avx_srl(i, w) \
  _mm256_srl_epi16(i, _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, w))
#define _avx_eq _mm256_cmpeq_epi16
#define _avx_lt _mm256_cmplt_epi16
#define _avx_gt _mm256_cmpgt_epi16
#define _avx_set1 _mm256_set1_epi16
#define _avx_insert _avx256_insert_16
#define _avx_extract _avx256_extract_16
#define _avx_m_mask _avx256_movemask_16
#define _avx_add _mm256_add_epi16
#define _avx_sub _mm256_sub_epi16
#define _avx_rotate(a, w) _avx_rotate(a, w)
#define _avx_broadcast(a) _mm256_broadcastq_epi64(a)

static inline __m256i _avx256_insert_16(__m256i reg, uint32_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm256_insert_epi32(reg, val, 0);
    case 1:
      return _mm256_insert_epi32(reg, val, 1);
    case 2:
      return _mm256_insert_epi32(reg, val, 2);
    case 3:
      return _mm256_insert_epi32(reg, val, 3);
    case 4:
      return _mm256_insert_epi16(reg, val, 4);
    case 5:
      return _mm256_insert_epi16(reg, val, 5);
    case 6:
      return _mm256_insert_epi16(reg, val, 6);
    case 7:
      return _mm256_insert_epi16(reg, val, 7);
    case 8:
      return _mm256_insert_epi32(reg, val, 8);
    case 9:
      return _mm256_insert_epi32(reg, val, 9);
    case 10:
      return _mm256_insert_epi32(reg, val, 10);
    case 11:
      return _mm256_insert_epi32(reg, val, 11);
    case 12:
      return _mm256_insert_epi16(reg, val, 12);
    case 13:
      return _mm256_insert_epi16(reg, val, 13);
    case 14:
      return _mm256_insert_epi16(reg, val, 14);
    case 15:
      return _mm256_insert_epi16(reg, val, 15);
    default:
      return _mm256_setzero_si256();
  }
}

static inline uint16_t _avx256_extract_16(__m256i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm256_extract_epi16(reg, 0);
    case 1:
      return _mm256_extract_epi16(reg, 1);
    case 2:
      return _mm256_extract_epi16(reg, 2);
    case 3:
      return _mm256_extract_epi16(reg, 3);
    case 4:
      return _mm256_extract_epi16(reg, 4);
    case 5:
      return _mm256_extract_epi16(reg, 5);
    case 6:
      return _mm256_extract_epi16(reg, 6);
    case 7:
      return _mm256_extract_epi16(reg, 7);
    case 8:
      return _mm256_extract_epi16(reg, 8);
    case 9:
      return _mm256_extract_epi16(reg, 9);
    case 10:
      return _mm256_extract_epi16(reg, 10);
    case 11:
      return _mm256_extract_epi16(reg, 11);
    case 12:
      return _mm256_extract_epi16(reg, 12);
    case 13:
      return _mm256_extract_epi16(reg, 13);
    case 14:
      return _mm256_extract_epi16(reg, 14);
    case 15:
      return _mm256_extract_epi16(reg, 15);
    default:
      return 0;
  }
}

static inline __m256i _avx_rotate(__m256i reg, int width)
{
  width = width % 16;
  switch (width)
  {
    case 0:
      return reg;
    case 1:
      return _mm256_alignr_epi8(reg, reg, 2);
    case 2:
      return _mm256_alignr_epi8(reg, reg, 2 * 2);
    case 3:
      return _mm256_alignr_epi8(reg, reg, 2 * 3);
    case 4:
      return _mm256_alignr_epi8(reg, reg, 2 * 4);
    case 5:
      return _mm256_alignr_epi8(reg, reg, 2 * 5);
    case 6:
      return _mm256_alignr_epi8(reg, reg, 2 * 6);
    case 7:
      return _mm256_alignr_epi8(reg, reg, 2 * 7);
    case 8:
      return _mm256_alignr_epi8(reg, reg, 2 * 8);
    case 9:
      return _mm256_alignr_epi8(reg, reg, 2 * 9);
    case 10:
      return _mm256_alignr_epi8(reg, reg, 2 * 10);
    case 11:
      return _mm256_alignr_epi8(reg, reg, 2 * 11);
    case 12:
      return _mm256_alignr_epi8(reg, reg, 2 * 12);
    case 13:
      return _mm256_alignr_epi8(reg, reg, 2 * 13);
    case 14:
      return _mm256_alignr_epi8(reg, reg, 2 * 14);
    case 15:
      return _mm256_alignr_epi8(reg, reg, 2 * 15);
  }
  return reg;
}

static inline static inline int _avx256_movemask_16(__m256i reg)
{
  int mask = _mm256_movemask_epi8(reg);
  int mask16 = 0;
  for (int i = 0; i < 16; i += 2)
  {
    if ((mask >> (i + 1)) & 1) mask16 ^= (1 << (i / 2));
  }
  return mask16;
}

#elif defined(INT32)  ////////////////////////////////////// EPI32

#define VECTOR_SIZE 8

#define SUB_POLY_TYPE uint32_t

#define FIXED_VARS 1

#if defined(_DEBUG)
#define _DEBUG_READ_SUB(p) scanf("%" SCNu32, &p)
#endif

#define _avx_sll _mm256_sll_epi32
#define _avx_srl _mm256_srl_epi32
#define _avx_eq _mm256_cmpeq_epi32
#define _avx_lt _mm256_cmplt_epi32
#define _avx_gt _mm256_cmpgt_epi32
#define _avx_set1 _mm256_set1_epi32
#define _avx_insert _avx256_insert_32
#define _avx_extract _avx256_extract_32
#define _avx_m_mask(a) _mm256_movemask_ps((__m256)a)
#define _avx_add _mm256_add_epi32
#define _avx_sub _mm256_sub_epi32
#define _avx_rotate(a, w) _avx_rotate(a, w)
#define _avx_broadcast(a) _mm256_broadcastsi128_si256(_mm256_castsi256_si128(a))

static inline __m256i _avx256_insert_32(__m256i reg, uint32_t val, int index)
{
  switch (index)
  {
    case 0:
      return _mm256_insert_epi32(reg, val, 0);
    case 1:
      return _mm256_insert_epi32(reg, val, 1);
    case 2:
      return _mm256_insert_epi32(reg, val, 2);
    case 3:
      return _mm256_insert_epi32(reg, val, 3);
    case 4:
      return _mm256_insert_epi32(reg, val, 4);
    case 5:
      return _mm256_insert_epi32(reg, val, 5);
    case 6:
      return _mm256_insert_epi32(reg, val, 6);
    case 7:
      return _mm256_insert_epi32(reg, val, 7);
    default:
      return _mm256_setzero_si256();
  }
}

static inline uint32_t _avx256_extract_32(__m256i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm256_extract_epi32(reg, 0);
    case 1:
      return _mm256_extract_epi32(reg, 1);
    case 2:
      return _mm256_extract_epi32(reg, 2);
    case 3:
      return _mm256_extract_epi32(reg, 3);
    case 4:
      return _mm256_extract_epi32(reg, 4);
    case 5:
      return _mm256_extract_epi32(reg, 5);
    case 6:
      return _mm256_extract_epi32(reg, 6);
    case 7:
      return _mm256_extract_epi32(reg, 7);
    default:
      return 0;
  }
}

static inline __m256i _avx_rotate(__m256i reg, int width)
{
  width = width % 8;
  switch (width)
  {
    case 0:
      return reg;
    case 1:
      return _mm256_alignr_epi8(reg, reg, 4);
    case 2:
      return _mm256_alignr_epi8(reg, reg, 4 * 2);
    case 3:
      return _mm256_alignr_epi8(reg, reg, 4 * 3);
    case 4:
      return _mm256_alignr_epi8(reg, reg, 4 * 4);
    case 5:
      return _mm256_alignr_epi8(reg, reg, 4 * 5);
    case 6:
      return _mm256_alignr_epi8(reg, reg, 4 * 6);
    case 7:
      return _mm256_alignr_epi8(reg, reg, 4 * 7);
  }
  return reg;
}

#elif defined(INT64)  ////////////////////////////////////// EPI64

#define VECTOR_SIZE 4

#define SUB_POLY_TYPE uint64_t

#define FIXED_VARS 0

#if defined(_DEBUG)
#define _DEBUG_READ_SUB(p) scanf("%" SCNu64, &p)
#endif

#define _avx_sll _mm256_sll_epi64
#define _avx_srl _mm256_srl_epi64
#define _avx_eq _mm256_cmpeq_epi64
#define _avx_lt _mm256_cmplt_epi64
#define _avx_gt _mm256_cmpgt_epi64
#define _avx_set1 _mm256_set1_epi64
#define _avx_insert _avx256_insert_64
#define _avx_extract _avx256_extract_64
#define _avx_m_mask(a) _mm256_movemask_pd((__m256d)a)
#define _avx_add _mm256_add_epi64
#define _avx_sub _mm256_sub_epi64
#define _avx_rotate(a, w) _avx_rotate(a, w)
#define _avx_testz(a, b) _mm_testz_si128(a, a)

static inline __m256i _avx256_insert_64(__m256i reg, uint64_t val, int idx)
{
  switch (idx)
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

static inline uint64_t _avx256_extract_64(__m256i reg, int idx)
{
  switch (idx)
  {
    case 0:
      return _mm256_exctract_epi64(reg, 0);
    case 1:
      return _mm256_exctract_epi64(reg, 1);
    case 2:
      return _mm256_exctract_epi64(reg, 2);
    case 3:
      return _mm256_exctract_epi64(reg, 3);
    default:
      return _mm256_setzero_si256();
  }
}

#endif

#endif  // !VEC256_CONFIG_H