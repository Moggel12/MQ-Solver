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

// Remove this definition if polynomials cannot be represented with ordinary
// integer registers and operations.
#define IS_INTEGER_REPR

// The operations below are assumed used on bitsliced polynomials only.
#define POLY_TYPE uint32_t
#define VARS_TYPE uint32_t

#if defined(_DEBUG)
#define _DEBUG_READ_P(p) scanf("%" SCNu32, &p)
#define _DEBUG_READ_V(v) scanf("%" SCNu32, &v)
#endif

#define GF2_ADD(a, b) (a ^ b)
#define GF2_MUL(a, b) (a & b)
#define POLY_IDX(p, i) (p & (1 << i))
#define INC(i) i + 1

// Assumed used for bitsliced polynomials.
typedef POLY_TYPE poly_t;
typedef VARS_TYPE vars_t;

#endif  // !MQ_CONFIG_H
