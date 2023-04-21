#ifndef MQ_VEC_H
#define MQ_VEC_H

#if defined(REG128)  ////////////// 128 bits, split register into multiple
                     /// systems

#define VEC_GF2_ADD(a, b) _mm_xor_si128(a, b)
#define VEC_GF2_MUL(a, b) _mm_and_si128(a, b)

#define VEC_LSHIFT(i, w) _mm_sll_epi32(i, _mm_set_epi32(0, 0, 0, w))
#define VEC_RSHIFT(i, w) _mm_srl_epi32(i, _mm_set_epi32(0, 0, 0, w))
#define VEC_IDX(p, i) _mm_and_si128(POLY_RSHIFT(p, i), _mm_set1_epi32(1))
#define VEC_SETBIT(p, i, b) _mm_xor_si128(p, _mm_set1_epi32(b << i))
#define VEC_FF _mm_set1_epi32(-1u)
#define VEC_0 _mm_setzero_si128()
#define VEC_1 _mm_set1_epi32(1)
#define VEC_LSB(i) _mm_and_si128(i, -i)
#define VEC_IS_ZERO(p) _mm_cmpeq_epi32(p, VEC_0)
#define VEC_MASK(b) _mm_set1_epi32((1 << b) - 1)
#define VEC_EQ(a, b) _mm_cmpeq_epi32(a, b)
#define VEC_LT(a, b) _mm_cmplt_epi32(a, b)
#define VEC_GT(a, b) _mm_cmpgt_epi32(a, b)
#define VEC_ASSIGN(a, b, c, d) _mm_set_epi32(a, b, c, d)
#define VEC_COND_LOAD(a, cond) _mm_maskload_epi32(a, cond)
#define VEC_ASSIGN_ONE(a) _mm_set1_epi32(a)
#define VEC_BLEND(a, b, mask) _mm_blendv_epi8(a, b, mask)
#define VEC_EXTRACT_SOL _avx128_extract_sol
#define VEC_SOL_OVERLAP _avx128_sol_overlap
#define VEC_AND(a, b) ((__m128i)_mm_and_ps((__m128)a, (__m128)b))
#define VEC_BITCMP(a, b) (!_mm_testz_si128(a, b))
#define VEC_MAX(a) _avx128_max(a)
#define VEC_H_OR(a) (!_mm_testz_ps((__m128)a, (__m128)a))
#define VEC_M_MASK(a) _mm_movemask_ps((__m128)a)
#define C_TO_VEC(a) _mm_cvtsi32_si128(a)
#define VEC_INSERT(reg, val, i) _avx128_insert(reg, val, i)

#define VEC_ADD(a, b) _mm_add_epi32(a, b)
#define VEC_SUB(a, b) _mm_sub_epi32(a, b)

#define VEC_EXTRACT(a, i) _mm_extract_epi32(a, i)

#elif defined( \
    REG256)  ////////////// 256 bits, split register into multiple systems

#define VEC_GF2_ADD(a, b) _mm256_xor_si256(a, b)
#define VEC_GF2_MUL(a, b) _mm256_and_si256(a, b)

#define VEC_LSHIFT(i, w) _mm256_sll_epi64(i, _mm256_set_epi64x(0, 0, 0, w))
#define VEC_RSHIFT(i, w) _mm256_srl_epi64(i, _mm256_set_epi64x(0, 0, 0, w))
#define VEC_IDX(v, i) _mm256_and_si256(POLY_RSHIFT(v, i), _mm256_set1_epi64x(1))
#define VEC_SETBIT(p, i, b) _mm256_xor_si256(p, _mm256_set1_epi64x(b << i))
#define VEC_FF _mm256_set1_epi64x((-1ull))
#define VEC_0 _mm256_setzero_si256()
#define VEC_1 _mm256_set1_epi64x(1)
#define VEC_LSB(i) _mm256_and_si256(i, -i)
#define VEC_IS_ZERO(v) (v == 0)
#define VEC_MASK(b) (_mm256_set1_epi64x((1 << b) - 1))
#define VEC_EQ(a, b) _mm256_cmpeq_epi64(a, b)
#define VEC_LT(a, b) _mm256_cmplt_epi64(a, b)
#define VEC_GT(a, b) _mm256_cmpgt_epi64(a, b)
#define VEC_ASSIGN(a, b, c, d) _mm256_set_epi64x(a, b, c, d)
#define VEC_ASSIGN_ONE(a) _mm256_set1_epi64x(a)
#define VEC_BLEND(a, b, mask) _mm256_blendv_epi8(a, b, mask)
#define VEC_EXTRACT_SOL _avx256_extract_sol
#define VEC_SOL_OVERLAP _avx256_sol_overlap
#define VEC_AND(a, b) ((__m256i)_mm256_and_pd((__m256d)a, (__m256d)b))
#define VEC_BITCMP(a, b) (!_mm256_testz_si256(a, b))
#define VEC_MAX(a) _avx256_max(a)
#define VEC_H_OR(a) (!_mm256_testz_pd((__m256d)a, (__m256d)a))
#define VEC_M_MASK(a) _mm256_movemask_pd((__m256d)a)
#define C_TO_VEC(a) _mm256_set_epi64x(0, 0, 0, a)
#define VEC_INSERT(reg, val, i) _avx256_insert(reg, val, i)

#define VEC_BROAD_ADD(a, b) _mm256_add_epi64(a, b)
#define VEC_BROAD_SUB(a, b) _mm256_sub_epi64(a, b)

#define VEC_EXTRACT(a, i) (container_t) _mm256_extract_epi64(a, i)

#endif

#endif  // !MQ_VEC_H