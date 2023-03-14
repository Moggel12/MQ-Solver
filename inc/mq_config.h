#ifndef MQ_CONFIG_H
#define MQ_CONFIG_H

// TODO: Possibly change n and m to size_t, just looks a bit nicer.
// TODO: Make use of Macros for indexing, etc. (makes things more generic)

#include <immintrin.h>
#include <stdint.h>

#define MAX_HISTORY 30

#if defined(_DEBUG)
#define RSEED 42
#else
#define RSEED time(NULL)
// #define RSEED 42
#endif

// TODO: Change to function-like macro.
//////////////// POLYNOMIAL SIZE DEFINITIONS ////////////////
#if defined(P8)  ////////////// <= 8 polynomials

#define POLY_TYPE uint8_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu8, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P16)  ////////////// <= 16 polynomials

#define POLY_TYPE uint16_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu16, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P32)  ////////////// <= 32 polynomials

#define POLY_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P64)  ////////////// <= 64 polynomials

#define POLY_TYPE uint64_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu64, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P128)  ////////////// <= 128 polynomials

#define POLY_TYPE __m128i
#if defined(_DEBUG)
#define _DEBUG_READ_P(p)    \
  uint64_t high;            \
  uint64_t low;             \
  scanf("%" SCNu64, &high); \
  scanf("%" SCNu64, &low);  \
  p = _mm_set_epi64(high, low);
#endif

#define POLY_LSHIFT(i, w) i << w
#define POLY_RSHIFT(i, w) i >> w
#define POLY_IDX(p, i) ((p >> i) & 1)
#define POLY_SETBIT(p, i, b) (p ^ (b << i))
#define POLY_FF (-1u)
#define POLY_0
#define POLY_1 1
#define POLY_LSB(i) (i & -i)
#define POLY_IS_ZERO(p) (p == 0)
#define POLY_MASK(b) ((1 << b) - 1)

#elif defined(P256)  ////////////// <= 256 polynomials

#define POLY_TYPE __m256i
#if defined(_DEBUG)
#define _DEBUG_READ_P(p)        \
  uint64_t high;                \
  uint64_t high_mid;            \
  uint64_t low_mid;             \
  uint64_t low;                 \
  scanf("%" SCNu64, &high);     \
  scanf("%" SCNu64, &high_mid); \
  scanf("%" SCNu64, &low_mid);  \
  scanf("%" SCNu64, &low);      \
  p = _mm256_set_epi64x(high, high_mid, low_mid, low);
#endif

#define POLY_LSHIFT(i, w) i << w
#define POLY_RSHIFT(i, w) i >> w
#define POLY_IDX(p, i) ((p >> i) & 1)
#define POLY_SETBIT(p, i, b) (p ^ (b << i))
#define POLY_FF (-1u)
#define POLY_0
#define POLY_1 1
#define POLY_LSB(i) (i & -i)
#define POLY_IS_ZERO(p) (p == 0)
#define POLY_MASK(b) ((1 << b) - 1)

#else  ////////////// Default case <= 32 polynomials

#define POLY_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif
#define IS_INTEGER_REPR

#endif

// TODO: Change to function-like macro.
//////////////// VARIABLE SIZE DEFINITIONS ////////////////
#if defined(V8)

#define VARS_TYPE uint8_t  ////////////// <= 8 variables
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu8, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V16)  ////////////// <= 16 variables

#define VARS_TYPE uint16_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu16, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V32)  ////////////// <= 32 variables

#define VARS_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V64)  ////////////// <= 64 variables

#define VARS_TYPE uint64_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu64, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V128)  ////////////// <= 128 variables

#define VARS_TYPE __m128i
#if defined(_DEBUG)
#define _DEBUG_READ_V(v)    \
  uint64_t high;            \
  uint64_t low;             \
  scanf("%" SCNu64, &high); \
  scanf("%" SCNu64, &low);  \
  v = _mm_set_epi64(high, low);
#endif

#define VARS_LSHIFT(i, w) POLY_LSHIFT(i, w)
#define VARS_RSHIFT(i, w) POLY_RSHIFT(i, w)
#define VARS_IDX(v, i) POLY_IDX(v, i)
#define VARS_SETBIT(v, i, b) POLY_SETBIT(v, i, b)
#define VARS_FF POLY_FF
#define VARS_0 POLY_0
#define VARS_1 POLY_1
#define VARS_LSB(i) POLY_LSB(i)
#define VARS_IS_ZERO(v) POLY_IS_ZERO(v)
#define VARS_MASK(b) ((1 << b) - 1)
#define VARS_EQ(a, b) a == b

#elif defined(V256)  ////////////// <= 256 variables

#define VARS_TYPE __m256i
#if defined(_DEBUG)
#define _DEBUG_READ_V(v)        \
  uint64_t high;                \
  uint64_t high_mid;            \
  uint64_t low_mid;             \
  uint64_t low;                 \
  scanf("%" SCNu64, &high);     \
  scanf("%" SCNu64, &high_mid); \
  scanf("%" SCNu64, &low_mid);  \
  scanf("%" SCNu64, &low);      \
  v = _mm256_set_epi64x(high, high_mid, low_mid, low);
#endif

#define VARS_LSHIFT(i, w) POLY_LSHIFT(i, w)
#define VARS_RSHIFT(i, w) POLY_RSHIFT(i, w)
#define VARS_IDX(v, i) POLY_IDX(v, i)
#define VARS_SETBIT(v, i, b) POLY_SETBIT(v, i, b)
#define VARS_FF POLY_FF
#define VARS_0 POLY_0
#define VARS_1 POLY_1
#define VARS_LSB(i) POLY_LSB(i)
#define VARS_IS_ZERO(v) POLY_IS_ZERO(v)
#define VARS_MASK(b) ((1 << b) - 1)
#define VARS_EQ(a, b) a == b

#else  ////////////// Default case <= 32 variables

#define VARS_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#endif

//////////////// MACROS FOR OPERATIONS ON INTEGERS ////////////////

#ifdef IS_INTEGER_REPR

#define GF2_ADD(a, b) (a ^ b)
#define GF2_MUL(a, b) (a & b)
#define BIT_OR(a, b) (a | b)
#define INC(i) i + 1

#define POLY_LSHIFT(i, w) i << w
#define POLY_RSHIFT(i, w) i >> w
#define POLY_IDX(p, i) ((p >> i) & 1)
#define POLY_SETBIT(p, i, b) (p ^ (b << i))
#define POLY_FF (-1u)
#define POLY_0 (0)
#define POLY_1 1
#define POLY_LSB(i) (i & -i)
#define POLY_IS_ZERO(p) (p == 0)
#define POLY_MASK(b) ((1 << b) - 1)

#define VARS_LSHIFT(i, w) POLY_LSHIFT(i, w)
#define VARS_RSHIFT(i, w) POLY_RSHIFT(i, w)
#define VARS_IDX(v, i) POLY_IDX(v, i)
#define VARS_SETBIT(v, i, b) POLY_SETBIT(v, i, b)
#define VARS_FF POLY_FF
#define VARS_0 POLY_0
#define VARS_1 POLY_1
#define VARS_LSB(i) POLY_LSB(i)
#define VARS_IS_ZERO(v) POLY_IS_ZERO(v)
#define VARS_MASK(b) ((1 << b) - 1)
#define VARS_EQ(a, b) a == b

#endif

// Assumed used for bitsliced polynomials.
typedef POLY_TYPE poly_t;
typedef VARS_TYPE vars_t;

#endif  // !MQ_CONFIG_H

/*
#define GF2_ADD(a, b) (a ^ b)
#define GF2_MUL(a, b) (a & b)
#define BIT_OR(a, b) (a | b)
#define INC(i) i + 1



*/