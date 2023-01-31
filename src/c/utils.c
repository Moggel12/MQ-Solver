#include <limits.h>
#include "utils.h"

int lex_idx(unsigned int i, unsigned int j, unsigned int n) {
  int sum = 0;
  for (int k = 1; k < i + 2; k++) {
    sum += (n - k);
  }
  return n + sum - (n - j - 1);
}

#ifdef IS_INTEGER_REPR
int hamming_weight(int x) {
  x -= (x >> 1) & 0x5555555555555555;
  x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
  x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
  return ((x * 0x0101010101010101) & 0xffffffffffffffff ) >> 56;
}

int gray_code(int i) { return i ^ (i >> 1); }

// Bit twiddling hacks
unsigned int trailing_zeros(unsigned int v) {
  unsigned int c;
  if (v) {
    v = (v ^ (v - 1)) >> 1;  // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++) {
      v >>= 1;
    }
  } else {
    c = CHAR_BIT * sizeof(v);
  }
  return c;
}
#else
#error "Remove this error whenever a representation without ordinary c integers is used, and suitable operations implemented"
#endif

int n_choose_k(int n, int k) {
  if (k > n) return 0;
  if (k == n) return 1;
  if (k > (n - k)) k = n - k;
  double c = 1;
  for (long i = 1; i <= k; i++) { // Multiplicative form
    c *= n--;
    c /= i;
  }
  return (int) c;
}

uint64_t *eval(uint64_t *system, int n, int m, uint64_t *values) {
  // TODO
}
