#ifndef MQ_VEC_H
#define MQ_VEC_H
// TODO: Make checks to see if API below is available compile-time and change to
// use this API instead of VEC_

/*
 * These need to be available either through macros or as C functions (for
 * reference; see vec128_config.h and vec256_config.h)
 *
 * _avx_sll(a, w) : Left-shift the bits of each element of a by w bits.
 * _avx_srl(a, w) : Right-shifth the bits of each element of a by w bits.
 * _avx_eq(a, b) : Returns a vector mask containing the conditional a == b
 *                 (follows AVX mask structure).
 * _avx_lt(a, b) : Returns a vector mask containing the conditional a < b
 *                 (follows AVX mask structure).
 * _avx_gt(a, b) : Returns a vector mask containing the conditional a > b
 *                 (follows AVX mask structure).
 * _avx_set1(val) : Return a vector with all elements set to val.
 * _avx_insert(a, val, idx) : Insert the value val at idx into a.
 * _avx_extract(a, idx) : Extract the element at idx of the vector.
 * _avx_m_mask(a) : Movemask for the given vector-packing.
 * _avx_add(a, b) : Add entry-wise.
 * _avx_sub(a, b) : Subtract entry-wise.
 * _avx_rotate(a, w) : Rotate the elements of a by w "placements".
 * _avx_xor(a, b) : Bitwise xor of a and b.
 * _avx_and(a, b) : Bitwise and of a and b.
 * _avx_or(a, b) : Bitwise or of a and b.
 * _avx_blend(a, b, mask) : Blend the value of a and b using mask (see AVX blend
 *                          instructions).
 * _avx_testz(a, b) : Test whether the bitwise and of a and b returns zero or
 *                    not (see testz in AVX).
 * _avx_zero() : A vector register of all zeroes.
 * _avx_broadcast(a) : Broadcast The lower four elements onto an entire vector.
 */

#if defined(REG128)  ////////////// 128 bits, split register into multiple
                     /// systems

#include "vec128_config.h"

#define VEC_GF2_ADD(a, b) _mm_xor_si128(a, b)
#define VEC_GF2_MUL(a, b) _mm_and_si128(a, b)
#define VEC_AND(a, b) _mm_and_si128(a, b)

#define VEC_BLEND(a, b, mask) _mm_blendv_epi8(a, b, mask)
#define VEC_0 _mm_setzero_si128()

#elif defined( \
    REG256)  ////////////// 256 bits, split register into multiple systems

#include "vec256_config.h"

#define VEC_GF2_ADD(a, b) _mm256_xor_si256(a, b)
#define VEC_GF2_MUL(a, b) _mm256_and_si256(a, b)
#define VEC_AND(a, b) _mm256_and_si256(a, b)

#define VEC_0 _mm256_setzero_si256()
#define VEC_BLEND(a, b, mask) _mm256_blendv_epi8(a, b, mask)

#endif

#define VEC_LSHIFT(i, w) _avx_sll(i, w)
#define VEC_RSHIFT(i, w) _avx_srl(i, w)
#define VEC_SETBIT(p, i, b) _avx_xor(p, _avx_set1(b << i))
#define VEC_1 _avx_set1(1)
#define VEC_IS_ZERO(p) _avx_eq(p, VEC_0)
#define VEC_MASK(b) _avx_set1((1 << b) - 1)
#define VEC_EQ(a, b) _avx_eq(a, b)
#define VEC_LT(a, b) _avx_lt(a, b)
#define VEC_GT(a, b) _avx_gt(a, b)
#define VEC_ASSIGN_ONE(a) _avx_set1(a)
#define VEC_MAX(a) _avx_max(a)
#define VEC_M_MASK(a) _avx_m_mask(a)
#define VEC_INSERT(reg, val, i) _avx_insert(reg, val, i)

#define VEC_ADD(a, b) _avx_add(a, b)
#define VEC_SUB(a, b) _avx_sub(a, b)

#define VEC_EXTRACT(a, i) _avx_extract(a, i)

#define VECTORIZED_ROUNDS 4

#endif  // !MQ_VEC_H