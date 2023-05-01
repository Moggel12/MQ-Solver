#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <stddef.h>
#include <stdint.h>

#include "mq_config.h"

typedef struct PotentialSolution
{
  uint8_t fixed_var;
  container_t solution;
} PotentialSolution;

unsigned int gen_matrix(container_t *mat, unsigned int n_rows,
                        unsigned int n_columns);

container_t parity(container_t bits);

void print_register(__m128i a);  // TODO: Remove

int _avx_sol_overlap(container_vec_t reg);

int _avx_extract_sol(container_vec_t reg, PotentialSolution *solutions);

#endif  // !VECTOR_UTILS_H
