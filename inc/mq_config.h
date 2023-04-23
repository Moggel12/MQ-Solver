#ifndef MQ_CONFIG_H
#define MQ_CONFIG_H

#include <immintrin.h>
#include <stdint.h>

#include "mq_uni.h"

#if defined(REG128) || defined(REG256)
#include "mq_vec.h"
#endif

#define MAX_HISTORY 4  // TODO, change back

#if defined(_DEBUG)
#define RSEED 42
#else
#define RSEED time(NULL)
// #define RSEED 42
#endif

// TODO: Change to function-like macro.
//////////////// POLYNOMIAL SIZE DEFINITIONS ////////////////
#if defined(REG8)  ////////////// <= 8 polynomials

#define INT_TYPE uint8_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu8, &p)
#endif

#elif defined(REG16)  ////////////// <= 16 polynomials

#define INT_TYPE uint16_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu16, &p)
#endif

#elif defined(REG32)  ////////////// <= 32 polynomials

#define INT_TYPE uint32_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif

#elif defined(REG64)  ////////////// <= 64 polynomials

#define INT_TYPE uint64_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu64, &p)
#endif

#elif defined(REG128)  ////////////// vector of 64 bit integers

#define INT_TYPE uint32_t

typedef __m128i container_vec_t;

#define ALIGNMENT 16

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif

#elif defined(REG256)  ////////////// vector of 64 bit integers

#define INT_TYPE uint64_t

typedef __m256i container_vec_t;

#define ALIGNMENT 32

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu64, &p)
#endif

#else  ////////////// Default case <= 32 polynomials

#define INT_TYPE uint32_t
#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#endif

#endif

// Assumed used for bitsliced polynomials.
typedef INT_TYPE container_t;

#endif  // !MQ_CONFIG_H