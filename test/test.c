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

int test_compute_e_k(void)
{
  size_t l = read_size_t();
  size_t n = read_size_t();
  size_t sys_len = read_size_t();

  poly_t *mat = read_poly_t_array(l);
  if (!mat)
  {
    printf("Memory error\n");
    return 1;
  }

  poly_t *old_sys = read_poly_t_array(sys_len);
  if (!old_sys)
  {
    printf("Memory error\n");
    free(mat);
    return 1;
  }

  poly_t *new_sys = calloc(sys_len, sizeof(poly_t));
  if (!new_sys)
  {
    printf("Memory error\n");
    free(old_sys);
    free(mat);
    return 1;
  }

  unsigned int deg = read_size_t();

  poly_t *correct_new_sys = read_poly_t_array(sys_len);
  if (!correct_new_sys)
  {
    printf("Memory error\n");
    free(new_sys);
    free(old_sys);
    free(mat);
    return 1;
  }

  unsigned int res = compute_e_k(mat, new_sys, old_sys, l, n);
  srand(RSEED);
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

  free(mat);
  free(old_sys);
  free(new_sys);
  free(correct_new_sys);

  if (!acc)
  {
    return 1;
  }
  return 0;
}

int test_eval(void)
{
  size_t len_system = read_size_t();
  size_t n = read_uint();

  poly_t *system = read_poly_t_array(len_system);
  if (!system) return 1;

  poly_t *solutions = read_poly_t_array(1 << n);
  if (!solutions) return 1;

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

  free(system);
  free(solutions);

  if (!acc)
  {
    return 1;
  }
  return 0;
}

int test_gen_matrix(void)
{
  srand(RSEED);
  unsigned int l = read_uint();
  unsigned int m = read_uint();

  printf("%u %u\n", l, m);
  poly_t *sol_mat = read_poly_t_array(l);
  if (!sol_mat) return 1;

  poly_t *rand_mat = malloc(l * sizeof(poly_t));
  if (!rand_mat)
  {
    free(sol_mat);
    return 1;
  }

  unsigned int error = gen_matrix(rand_mat, l, m);
  if (error)
  {
    printf("Did not return rank %u matrix\n", l);
    free(rand_mat);
    free(sol_mat);
    return 1;
  }

  uint8_t acc = 1;
  for (unsigned int i = 0; i < l; i++)
  {
    acc = acc && (sol_mat[i] == rand_mat[i]);
    if (!acc)
    {
      printf("%u != %u\n", sol_mat[i], rand_mat[i]);
      break;
    }
  }

  free(rand_mat);
  free(sol_mat);

  if (!acc)
  {
    return 1;
  }
  return 0;
}

int test_solve_sanitized()
{
  unsigned int n = read_uint();
  unsigned int m = read_uint();
  size_t sys_len = read_size_t();

  vars_t *py_sol = read_vars_t_array(1);

  if (!py_sol) return 1;

  poly_t *system = read_poly_t_array(sys_len);

  if (!system)
  {
    free(py_sol);
    return 1;
  }

  vars_t c_sol = 0;

  solve(system, n, m, &c_sol);

  free(system);

  if (*py_sol != c_sol)
  {
    free(py_sol);
    return 1;
  }

  free(py_sol);

  return 0;
}

int main(void)
{
  unsigned int choice = read_uint();
  switch (choice)
  {
    case 0:
      return test_compute_e_k();
    case 1:
      return test_eval();
    case 2:
      return test_gen_matrix();
    case 3:
      return test_solve_sanitized();
    default:
      printf("Invalid function chosen: %u\n", choice);
  }
  return 1;
}
