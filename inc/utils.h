#ifndef UTILS_H
#define UTILS_H
#include <stddef.h>
#include <stdint.h>

#include "mq_config.h"

int lex_idx(unsigned int i, unsigned int j, unsigned int n);

unsigned int hamming_weight(unsigned int x);

int gray_code(int i);

poly_t eval(poly_t *system, size_t n, vars_t var_values);

int n_choose_k(int n, int k);

unsigned int trailing_zeros(unsigned int v);

unsigned int gen_matrix(poly_t *mat, unsigned int n_rows,
                        unsigned int n_columns);

size_t gray_to_bin(size_t i);

#endif  // !UTILS_H
