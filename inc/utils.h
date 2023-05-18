#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include "mq_config.h"

unsigned int hamming_weight(unsigned int x);

unsigned int trailing_zeros(unsigned int v);

int lex_idx(unsigned int i, unsigned int j, unsigned int n);

int n_choose_k(int n, int k);

poly_t eval(poly_t *system, size_t n, poly_t var_values);

unsigned int gen_matrix(poly_t *mat, unsigned int n_rows,
                        unsigned int n_columns);

#if defined(REG256) || defined(REG128) || defined(REG64)

unsigned long long llrand(void);

#define _RAND() llrand()

#else 

#define _RAND() rand()

#endif

#endif  // !COMMON_UTILS_H