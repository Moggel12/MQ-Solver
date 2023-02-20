#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mq.h"
#include "mq_config.h"
#include "utils.h"

poly_t *read_poly_t_array(size_t len)
{
  poly_t *arr = malloc(len * sizeof(poly_t));
  if (!arr) return NULL;

  for (size_t i = 0; i < len; i++)
  {
    _DEBUG_READ_P(arr[i]);
  }

  return arr;
}

vars_t *read_vars_t_array(size_t len)
{
  vars_t *arr = malloc(len * sizeof(vars_t));
  if (!arr) return NULL;

  for (size_t i = 0; i < len; i++)
  {
    _DEBUG_READ_V(arr[i]);
  }

  return arr;
}

size_t read_size_t()
{
  size_t val;
  scanf("%zu\n", &val);
  return val;
}

unsigned int read_uint()
{
  unsigned int val;
  scanf("%u\n", &val);
  return val;
}

void test_compute_e_k(void)
{
  size_t l = read_size_t();
  size_t n = read_size_t();
  size_t sys_len = read_size_t();

  poly_t *mat = read_poly_t_array(l);
  if (!mat)
  {
    printf("Memory error\n");
    return;
  }

  poly_t *old_sys = read_poly_t_array(sys_len);
  if (!old_sys)
  {
    printf("Memory error\n");
    free(mat);
    return;
  }

  poly_t *new_sys = calloc(sys_len, sizeof(poly_t));
  if (!new_sys)
  {
    printf("Memory error\n");
    free(old_sys);
    free(mat);
    return;
  }

  unsigned int deg = read_size_t();

  poly_t *correct_new_sys = read_poly_t_array(sys_len);
  if (!correct_new_sys)
  {
    printf("Memory error\n");
    free(new_sys);
    free(old_sys);
    free(mat);
    return;
  }

  unsigned int res = compute_e_k(mat, new_sys, old_sys, l, n);

  uint8_t acc = 1;
  for (unsigned int i = 0; i < n_choose_k(n, 2) + n + 1; i++)
  {
    acc = acc && (new_sys[i] == correct_new_sys[i]);
    if (!acc)
    {
      printf("%u != %u\n", new_sys[i], correct_new_sys[i]);
      break;
    }
  }

  if (deg != res)
  {
    printf("Degrees do not correspond: %u != %u\n", deg, res);
    acc = 0;
  }

  if (acc)
  {
    printf("No errors found\n");
  }

  free(mat);
  free(old_sys);
  free(new_sys);
  free(correct_new_sys);
}

void test_eval(void)
{
  size_t len_system = read_size_t();
  size_t n = read_uint();

  poly_t *system = read_poly_t_array(len_system);
  if (!system) return;

  poly_t *solutions = read_poly_t_array(1 << n);
  if (!solutions) return;

  uint8_t acc = 1;
  for (size_t i = 0; i < (1u << n); i++)
  {
    poly_t val = eval(system, n, i);

    acc = acc && (solutions[i] == val);
    if (!acc)
    {
      printf("Found differing values: %u != %u\n", val, solutions[i]);
      break;
    }
  }

  if (acc)
  {
    printf("No errors found!\n");
  }
  free(system);
  free(solutions);
}

void test_gen_matrix(void)
{
  unsigned int l = read_uint();
  unsigned int m = read_uint();

  printf("%u %u\n", l, m);
  poly_t *sol_mat = read_poly_t_array(l);
  if (!sol_mat) return;

  poly_t *rand_mat = malloc(l * sizeof(poly_t));
  if (!rand_mat)
  {
    free(sol_mat);
    return;
  }

  unsigned int error = gen_matrix(rand_mat, l, m);
  if (error)
  {
    printf("Did not return rank %u matrix\n", l);
    free(rand_mat);
    free(sol_mat);
    return;
  }

  uint8_t acc = 1;
  for (unsigned int i = 0; i < l; i++)
  {
    acc = acc && (sol_mat[i] == rand_mat[i]);
    if (!acc)
    {
      printf("%u != %u", sol_mat[i], rand_mat[i]);
      break;
    }
  }

  free(rand_mat);
  free(sol_mat);

  if (acc)
  {
    printf("Found no errors.\n");
  }
}

void test_memory();

int main(void)
{
  // test_gen_matrix();
  poly_t sys[16] = {14, 14, 21, 19, 26, 8, 29, 26, 29, 8, 31, 16, 9, 24, 12, 1};
  unsigned int n = 5;
  unsigned int m = 5;
  vars_t sol;
  uint8_t err = solve(sys, n, m, &sol);
  printf("Error code %u\n", err);

  return 0;
}
