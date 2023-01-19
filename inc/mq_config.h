//
// Created by mikkelvestergaard on 1/12/23.
//

#ifndef MQ_SOLVER_MQ_CONFIG_H
#define MQ_SOLVER_MQ_CONFIG_H

#include "stdint.h"

// Remove this definition if polynomials cannot be represented with ordinary integer registers and operations.
#define IS_INTEGER_REPR

// The operations below are assumed used on bitsliced polynomials only.
#define POLY_AMOUNT uint8_t
#define POLY_VARS uint8_t
#define GF2_ADD(a, b) (a ^ b)
#define GF2_MUL(a, b) (a & b)
#define POLY_IDX(p, i) (p & (1 << i))
#define INC(i) i + 1

// Assumed used for bitsliced polynomials.
typedef POLY_AMOUNT poly_t;
typedef POLY_VARS vars_t;

#endif //MQ_SOLVER_MQ_CONFIG_H
