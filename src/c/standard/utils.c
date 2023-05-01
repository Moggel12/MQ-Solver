#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "benchmark.h"
#include "common_utils.h"
#include "standard_utils.h"

unsigned int hamming_weight(unsigned int x) { return __builtin_popcount(x); }

int gray_code(int i) { return i ^ (i >> 1); }  // TODO: Alter and use in solve()

unsigned int trailing_zeros(unsigned int v)
{
  unsigned int c;
  if (v)
  {
    v = (v ^ (v - 1)) >> 1;
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

container_t parity(container_t bits) { return __builtin_parity(bits); }

int lex_idx(unsigned int i, unsigned int j, unsigned int n)
{
  int sum = 0;
  for (unsigned int k = 1; k < i + 2; k++)
  {
    sum += (n - k);
  }
  return n + sum - (n - j - 1);
}

int n_choose_k(int n, int k)
{
  if (k > n) return 0;
  if (k == n) return 1;
  if (k > (n - k)) k = n - k;
  double c = 1;
  for (long i = 1; i <= k; i++)
  {
    c *= n--;
    c /= i;
  }
  return (int)c;
}

container_t gen_row(unsigned int m) { return (rand() & ((1 << m) - 1)); }

unsigned int gen_matrix(container_t *mat, unsigned int n_rows,
                        unsigned int n_columns)
{
  container_t *mat_copy = malloc(n_rows * sizeof(container_t));
  if (!mat_copy) return 1;

  unsigned int rank;

  do
  {
    rank = 0;
    for (unsigned int i = 0; i < n_rows; i++)
    {
      mat_copy[i] = mat[i] = gen_row(n_columns);
    }
    for (unsigned int i = 0; i < n_rows; i++)
    {
      if (!INT_IS_ZERO(mat_copy[i]))
      {
        rank++;
        unsigned int pivot_elm = INT_LSB(mat_copy[i]);
        for (unsigned int j = i + 1; j < n_rows; j++)
        {
          if (!INT_IS_ZERO(GF2_MUL(mat_copy[j], pivot_elm)))
          {
            mat_copy[j] = GF2_ADD(mat_copy[j], mat_copy[i]);
          }
        }
      }
    }

  } while (rank != n_rows);

  free(mat_copy);

  return rank != n_rows;
}

container_t eval(container_t *system, size_t n, container_t var_values)
{
  BEGIN_BENCH(g_eval_time)

  container_t table[2] = {INT_0, INT_FF};

  container_t res = system[0];

  for (unsigned int i = 0; i < n; i++)
  {
    res = GF2_ADD(res, GF2_MUL(table[INT_IDX(var_values, i)], system[i + 1]));
  }

  int idx = lex_idx(0, 1, n);
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = i + 1; j < n; j++)
    {
      container_t monomial =
          GF2_MUL(table[INT_IDX(var_values, i)], table[INT_IDX(var_values, j)]);
      res = GF2_ADD(res, GF2_MUL(monomial, system[idx]));
      idx++;
    }
  }

  END_BENCH(g_eval_time)

  return res;
}