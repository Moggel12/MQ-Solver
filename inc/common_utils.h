#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

unsigned int hamming_weight(unsigned int x);

unsigned int trailing_zeros(unsigned int v);

int lex_idx(unsigned int i, unsigned int j, unsigned int n);

int n_choose_k(int n, int k);

#endif  // !COMMON_UTILS_H