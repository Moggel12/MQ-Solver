#ifndef MQ_CONFIG_H
#define MQ_CONFIG_H

// TODO: Possibly change n and m to size_t, just looks a bit nicer.
// TODO: Make use of Macros for indexing, etc. (makes things more generic)

#include <stdint.h>

#define MAX_HISTORY 30

#if defined(_DEBUG)
#define RSEED 42
#else
// #define RSEED time(NULL)
#define RSEED 42
#endif

// TODO: Change to function-like macro.
//////////////// POLYNOMIAL SIZE DEFINITIONS ////////////////
#if defined(P8)

#define POLY_TYPE uint8_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu8, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P16)

#define POLY_TYPE uint16_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu16, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P32)

#define POLY_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P64)

#define POLY_TYPE uint64_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu64, &p)
#endif
#define IS_INTEGER_REPR

#elif defined(P128)

#elif defined(P256)

#else  // Default case

#define POLY_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif
#define IS_INTEGER_REPR

#endif

// TODO: Change to function-like macro.
//////////////// VARIABLE SIZE DEFINITIONS ////////////////
#if defined(V8)

#define VARS_TYPE uint8_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu8, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V16)

#define VARS_TYPE uint16_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu16, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V32)

#define VARS_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V64)

#define VARS_TYPE uint64_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu64, &v)
#endif
#define IS_INTEGER_REPR

#elif defined(V128)

#elif defined(V256)

#else  // Default case

#define VARS_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#endif

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
