#ifndef UTILS_H
#define UTILS_H
#include <stddef.h>
#include <stdint.h>

#include "mq_config.h"

#if defined(REG128) || defined(REG256)

uint8_t gen_matrix(container_vec_t *mat, unsigned int n_rows,
                   unsigned int n_columns);

container_vec_t parity(container_vec_t bits);

#if defined(REG128)  ///////////////// AVX 128

void print_register(__m128i a);  // TODO: Remove

__m128i _avx128_insert(__m128i reg, uint32_t val, int index);

uint32_t _avx128_max();

uint32_t _avx128_extract_sol();

int _avx128_sol_overlap();

#endif  ///////////////// AVX 128

#if defined(REG256)  //////////////// AVX 256

__m256i _avx256_insert(__m256i reg, uint64_t val, int index);

uint64_t _avx256_max();

int _avx256_sol_overlap();

uint64_t _avx256_extract_sol();

#endif  //////////////// AVX 256

#else

unsigned int gen_matrix(container_t *mat, unsigned int n_rows,
                        unsigned int n_columns);

container_t parity(container_t bits);

#endif

unsigned int hamming_weight(unsigned int x);

unsigned int trailing_zeros(unsigned int v);

int lex_idx(unsigned int i, unsigned int j, unsigned int n);

int n_choose_k(int n, int k);

container_t eval(container_t *system, size_t n, container_t var_values);

#endif  // !UTILS_H
