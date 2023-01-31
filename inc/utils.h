#ifndef UTILS_H
#define UTILS_H
#include <stdint.h>
#include "mq_config.h"

int lex_idx(unsigned int i, unsigned int j, unsigned int n);

int hamming_weight(int x);

int gray_code(int i);

uint64_t *eval(uint64_t *system, int n, int m, uint64_t *values);

int n_choose_k(int n, int k);

unsigned int trailing_zeros(unsigned int v);

#endif // !UTILS_H
