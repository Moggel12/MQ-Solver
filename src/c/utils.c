#include "utils.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int lex_idx(unsigned int i, unsigned int j, unsigned int n)
{
  int sum = 0;
  for (unsigned int k = 1; k < i + 2; k++)
  {
    sum += (n - k);
  }
  return n + sum - (n - j - 1);
}

#ifdef IS_INTEGER_REPR
int hamming_weight(int x)
{
  x -= (x >> 1) & 0x5555555555555555;
  x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
  x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
  return ((x * 0x0101010101010101) & 0xffffffffffffffff) >> 56;
}

int gray_code(int i) { return i ^ (i >> 1); }

// Bit twiddling hacks
unsigned int trailing_zeros(unsigned int v)
{
  unsigned int c;
  if (v)
  {
    v = (v ^ (v - 1)) >> 1;  // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++)
    {
      v >>= 1;
    }
  }
  else
  {
    c = (-1u);
  }
  return c;
}
#else
#error \
    "Remove this error whenever a representation without ordinary c integers is used, and suitable operations implemented"
#endif

int n_choose_k(int n, int k)
{
  if (k > n) return 0;
  if (k == n) return 1;
  if (k > (n - k)) k = n - k;
  double c = 1;
  for (long i = 1; i <= k; i++)
  {  // Multiplicative form
    c *= n--;
    c /= i;
  }
  return (int)c;
}

unsigned int gen_matrix(poly_t *mat, unsigned int n_rows,
                        unsigned int n_columns)
{
  poly_t *mat_copy = malloc(n_rows * sizeof(poly_t));
  if (!mat_copy) return 1;

  // TODO: Change to srand
  srand(RSEED);

  unsigned int rank;

  do
  {
    rank = 0;
    for (unsigned int i = 0; i < n_rows; i++)
    {
      mat_copy[i] = mat[i] = (rand() & ((1 << n_columns) - 1));
    }

    for (unsigned int i = 0; i < n_rows; i++)
    {
      if (mat_copy[i])
      {
        rank++;
        unsigned int pivot_elm = mat_copy[i] & -mat_copy[i];
        for (unsigned int j = i + 1; j < n_rows; j++)
        {
          if (mat_copy[j] & pivot_elm)
          {
            mat_copy[j] ^= mat_copy[i];
          }
        }
      }
    }

  } while (rank != n_rows);

  free(mat_copy);

  return rank != n_rows;
}

poly_t eval(poly_t *system, size_t n, vars_t var_values)
{
  poly_t table[2] = {0, (-1)};

  poly_t res = system[0];

  for (unsigned int i = 0; i < n; i++)
  {
    res = GF2_ADD(res, GF2_MUL(table[(var_values >> i) & 1], system[i + 1]));
  }

  int idx = lex_idx(0, 1, n);
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = i + 1; j < n; j++)
    {
      poly_t monomial =
          GF2_MUL(table[(var_values >> i) & 1], table[(var_values >> j) & 1]);
      res = GF2_ADD(res, GF2_MUL(monomial, system[idx]));
      idx++;
    }
  }

  return res;
}

size_t gray_to_bin(size_t i)
{
  int mask = i;
  while (mask)
  {
    mask >>= 1;
    i ^= mask;
  }
  return i;
}