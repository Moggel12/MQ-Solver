#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "common_utils.h"
#include "mq.h"
#include "mq_config.h"
#if defined(REG128) || defined(REG256)
#include "vector_utils.h"
#else
#include "standard_utils.h"
#endif

container_t *read_container_t_array(size_t len)
{
  container_t *arr = malloc(len * sizeof(container_t));
  if (!arr) return NULL;

  for (size_t i = 0; i < len; i++)
  {
    _DEBUG_READ_P(arr[i]);
  }

  return arr;
}

container_t *read_container_t_array_static(container_t *arr, size_t len)
{
  for (size_t i = 0; i < len; i++)
  {
    _DEBUG_READ_P(arr[i]);
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

#if defined(REG128) || defined(REG256)  // TODO: Fix tests

int test_compute_e_k(void)
{
  size_t l = read_size_t();
  size_t n = read_size_t();
  size_t sys_len = read_size_t();

  container_t *mat = calloc(l * 4 * (1 << FIXED_VARS), sizeof(container_t));
  if (!mat)
  {
    printf("Memory error (allocating matrix)!\n");
    return 1;
  }

  for (int i = 0; i < (1 << FIXED_VARS) * 4; i++)
  {
    read_container_t_array_static(mat + i * l, l);
  }

  container_t *old_sys =
      calloc((1 << FIXED_VARS) * sys_len, sizeof(container_t));
  if (!old_sys)
  {
    printf("Memory error (allocating vectorized old sys)!\n");
    free(mat);
    return 1;
  }

  for (int i = 0; i < (1 << FIXED_VARS); i++)
  {
    read_container_t_array_static(old_sys + i * sys_len, sys_len);
  }

  container_vec_t *new_sys =
      aligned_alloc(ALIGNMENT, sizeof(container_vec_t) * sys_len);
  if (!new_sys)
  {
    printf("Memory error (allocating new sys)!\n");
    free(mat);
    free(old_sys);
    return 1;
  }

  memset(new_sys, 0, sizeof(container_vec_t) * sys_len);

  size_t correct_deg[VECTOR_SIZE];
  for (int i = 0; i < VECTOR_SIZE; i++)
  {
    correct_deg[i] = read_size_t();
  }

  container_t *correct_new_sys_arr[VECTOR_SIZE] = {0};

  for (int i = 0; i < VECTOR_SIZE; i++)
  {
    correct_new_sys_arr[i] = read_container_t_array(sys_len);
    if (!correct_new_sys_arr[i])
    {
      printf("Memory error (allocating correct sys)!\n");
      free(mat);
      free(old_sys);
      free(new_sys);
      for (int j = 0; j < i; j++)
      {
        free(correct_new_sys_arr[j]);
      }
      return 1;
    }
  }

  container_vec_t *correct_new_sys =
      aligned_alloc(ALIGNMENT, sys_len * sizeof(container_vec_t));
  if (!correct_new_sys)
  {
    printf("Memory error (allocating vectorized correct sys)!\n");
    free(new_sys);
    free(old_sys);
    for (int i = 0; i < VECTOR_SIZE; i++)
    {
      free(correct_new_sys_arr[i]);
    }
    free(mat);
    return 1;
  }

  for (size_t i = 0; i < sys_len; i++)
  {
    for (int j = 0; j < VECTOR_SIZE; j++)
    {
      correct_new_sys[i] =
          VEC_INSERT(correct_new_sys[i], correct_new_sys_arr[j][i], j);
    }
  }

  srand(RSEED);

  container_vec_t degrees = compute_e_k(mat, new_sys, old_sys, sys_len, l, n);

  uint8_t acc = 1;
  for (unsigned int i = 0; i < n_choose_k(n, 2) + n + 1; i++)
  {
    acc = acc && (VEC_M_MASK(VEC_EQ(new_sys[i], correct_new_sys[i])) ==
                  INT_MASK(VECTOR_SIZE));
    if (!acc)
    {
      print_register(new_sys[i]);
      print_register(correct_new_sys[i]);
      break;
    }
  }

  for (int i = 0; i < VECTOR_SIZE; i++)
  {
    if (correct_deg[i] != VEC_EXTRACT(degrees, VECTOR_SIZE - (i + 1)))
    {
      printf("Degrees do not correspond: %zu != %u\n", correct_deg[i],
             VEC_EXTRACT(degrees, VECTOR_SIZE - (i + 1)));
      acc = 0;
      break;
    }
  }

  free(mat);
  free(old_sys);
  free(new_sys);
  free(correct_new_sys);
  for (int i = 0; i < VECTOR_SIZE; i++)
  {
    free(correct_new_sys_arr[i]);
  }

  if (!acc)
  {
    return 1;
  }

  return 0;
}

int test_fix_poly(void)
{
  size_t n = read_size_t();
  size_t sys_len = read_size_t();
  size_t fixed_sys_len = read_size_t();
  size_t new_n = n - FIXED_VARS;

  container_t *old_sys = read_container_t_array(sys_len);
  if (!old_sys)
  {
    printf("Error allocating old system!\n");
    return 1;
  }

  container_t *correct_fixed_systems =
      calloc((1 << FIXED_VARS) * fixed_sys_len, sizeof(container_t));
  if (!correct_fixed_systems)
  {
    printf("Error allocating correct fixed systems!\n");
    free(old_sys);
    return 1;
  }
  for (int i = 0; i < (1 << FIXED_VARS); i++)
  {
    read_container_t_array_static(correct_fixed_systems + i * fixed_sys_len,
                                  fixed_sys_len);
  }

  container_t *fixed_systems =
      calloc((1 << FIXED_VARS) * fixed_sys_len, sizeof(container_t));
  if (!fixed_sys_len)
  {
    printf("Error allocating fixed systems!\n");
    free(correct_fixed_systems);
    free(old_sys);
    return 1;
  }

  container_t assignment[FIXED_VARS] = {0};
  for (container_t i = 0; i < (1 << FIXED_VARS); i++)
  {
    for (int j = 0; j < FIXED_VARS; j++) assignment[j] = (i >> j) & 1;
    fix_poly(old_sys, fixed_systems + i * fixed_sys_len, assignment, new_n, n);
  }

  uint8_t acc = 1;
  for (unsigned int i = 0; i < (1 << FIXED_VARS) * fixed_sys_len; i++)
  {
    acc = acc && (correct_fixed_systems[i] == fixed_systems[i]);
    if (!acc)
    {
      printf("System fixing not correct: %u != %u (index %u)\n",
             correct_fixed_systems[i], fixed_systems[i], i);
      free(old_sys);
      free(correct_fixed_systems);
      free(fixed_systems);
      return 1;
    }
  }

  free(old_sys);
  free(correct_fixed_systems);
  free(fixed_systems);

  return 0;
}

int test_solve_sanitized()
{
  unsigned int n = read_uint();
  unsigned int m = read_uint();
  size_t sys_len = read_size_t();

  container_t *system = read_container_t_array(sys_len);

  if (!system)
  {
    printf("Error receiving system\n");
    return 1;
  }

  container_t c_sol = 0;

  uint8_t error = solve(system, n, m, &c_sol);

  if (eval(system, n, c_sol) || error)
  {
    printf("Solution is invalid (error code: %u)\n", error);
    free(system);
    return 1;
  }

  free(system);

  return 0;
}

#else

int test_compute_e_k(void)
{
  size_t l = read_size_t();
  size_t n = read_size_t();
  size_t sys_len = read_size_t();

  container_t *mat = read_container_t_array(l);
  if (!mat)
  {
    printf("Memory error (allocating matrix)\n");
    return 1;
  }

  container_t *old_sys = read_container_t_array(sys_len);
  if (!old_sys)
  {
    printf("Memory error (allocating old sys)\n");
    free(mat);
    return 1;
  }

  container_t *new_sys = calloc(sys_len, sizeof(container_t));
  if (!new_sys)
  {
    printf("Memory error (allocating new sys)\n");
    free(old_sys);
    free(mat);
    return 1;
  }

  unsigned int deg = read_size_t();

  container_t *correct_new_sys = read_container_t_array(sys_len);
  if (!correct_new_sys)
  {
    printf("Memory error (allocating correct sys)\n");
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

int test_solve_sanitized()
{
  unsigned int n = read_uint();
  unsigned int m = read_uint();
  size_t sys_len = read_size_t();

  container_t *py_sol = read_container_t_array(1);

  if (!py_sol)
  {
    printf("Error receiving python solution\n");
    return 1;
  };

  container_t *system = read_container_t_array(sys_len);

  if (!system)
  {
    printf("Error receiving system\n");
    free(py_sol);
    return 1;
  }

  container_t c_sol = 0;

  solve(system, n, m, &c_sol);

  free(system);

  if (*py_sol != c_sol)
  {
    printf("Solutions differ\n");
    free(py_sol);
    return 1;
  }

  free(py_sol);

  return 0;
}

#endif

int test_gen_matrix(void)
{
  srand(RSEED);
  unsigned int l = read_uint();
  unsigned int m = read_uint();

  printf("%u %u\n", l, m);
  container_t *sol_mat = read_container_t_array(l);
  if (!sol_mat) return 1;

  container_t *rand_mat = malloc(l * sizeof(container_t));
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

int test_eval(void)
{
  size_t len_system = read_size_t();
  size_t n = read_uint();

  container_t *system = read_container_t_array(len_system);
  if (!system) return 1;

  container_t *solutions = read_container_t_array(1 << n);
  if (!solutions) return 1;

  uint8_t acc = 1;
  for (size_t i = 0; i < (1u << n); i++)
  {
    container_t val = eval(system, n, i);

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
#if defined(REG128) || defined(REG256)
    case 4:
      return test_fix_poly();
#endif
    default:
      printf("Invalid function chosen: %u\n", choice);
  }
  return 1;
}
