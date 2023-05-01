#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <stddef.h>
#include <stdint.h>

#include "mq_config.h"

typedef struct PotentialSolution
{
  uint8_t fixed_var;
  sub_poly_t solution;
} PotentialSolution;

sub_poly_t parity(poly_t bits);

int _avx_sol_overlap(poly_vec_t reg);

int _avx_extract_sol(poly_vec_t reg, PotentialSolution *solutions);

sub_poly_t _avx_max(poly_vec_t a);

#endif  // !VECTOR_UTILS_H
