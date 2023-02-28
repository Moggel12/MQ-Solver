#ifndef MQ_CONFIG_H
#define MQ_CONFIG_H

// TODO: Possibly change n and m to size_t, just looks a bit nicer.
// TODO: Make use of Macros for indexing, etc. (makes things more generic)

#include <stdint.h>

#define MAX_HISTORY 500

#if defined(_DEBUG)
#define RSEED 42
#else
// #define RSEED time(NULL)
#define RSEED 42
#endif

// TODO: Change to function-like macro.
//////////////// POLYNOMIAL SIZE DEFINITIONS ////////////////
#ifdef P8

#define VARS_TYPE uint8_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu8, &v)
#endif
#define IS_INTEGER_REPR

#elifdef P16

#define VARS_TYPE uint16_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu16, &v)
#endif
#define IS_INTEGER_REPR

#elifdef P32

#define VARS_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#elifdef P64

#define VARS_TYPE uint64_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu64, &v)
#endif
#define IS_INTEGER_REPR

#elifdef P128

#elifdef P256

#else  // Default case

#define POLY_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#endif

// TODO: Change to function-like macro.
//////////////// VARIABLE SIZE DEFINITIONS ////////////////
#ifdef V8
#define VARS_TYPE uint8_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu8, &v)
#endif
#define IS_INTEGER_REPR

#elifdef V16

#define VARS_TYPE uint16_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu16, &v)
#endif
#define IS_INTEGER_REPR

#elifdef V32

#define VARS_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif
#define IS_INTEGER_REPR

#elifdef V64

#define VARS_TYPE uint64_t
#if defined(_DEBUG)
#define _DEBUG_READ_V(v) scanf("%" SCNu64, &v)
#endif
#define IS_INTEGER_REPR

#elifdef V128

#elifdef V256

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
#define POLY_FF -1u
#define POLY_0 0

#define VARS_LSHIFT(i, w) POLY_LSHIFT(i, w)
#define VARS_RSHIFT(i, w) POLY_RSHIFT(i, w)
#define VARS_IDX(p, i) POLY_IDX(p, i)
#define VARS_SETBIT(p, i, b) POLY_SETBIT(p, i, b)
#define VARS_FF -1u
#define VARS_0 0

#endif

// Assumed used for bitsliced polynomials.
typedef POLY_TYPE poly_t;
typedef VARS_TYPE vars_t;

#endif  // !MQ_CONFIG_H
