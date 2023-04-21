#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mq.h"
#include "mq_config.h"
#include "utils.h"

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

#if defined(REG128) || defined(REG256)

int test_compute_e_k(void)
{
  size_t l = read_size_t();
  size_t n = read_size_t();
  size_t sys_len = read_size_t();
  container_t *mat_arr[4];

  for (int i = 0; i < 4; i++)
  {
    mat_arr[i] = read_container_t_array(l);
    if (!mat_arr[i])
    {
      printf("Memory error\n");
      return 1;
    }
  }

  container_vec_t *mat = aligned_alloc(ALIGNMENT, l * sizeof(container_vec_t));
  for (unsigned int i = 0; i < l; i++)
  {
    mat[i] =
        VEC_ASSIGN(mat_arr[0][i], mat_arr[1][i], mat_arr[2][i], mat_arr[3][i]);
  }

  container_t *old_sys = read_container_t_array(sys_len);
  if (!old_sys)
  {
    printf("Memory error\n");
    for (int i = 0; i < 4; i++) free(mat_arr[i]);
    return 1;
  }

  container_vec_t *new_sys =
      aligned_alloc(ALIGNMENT, sys_len * sizeof(container_vec_t));
  if (!new_sys)
  {
    printf("Memory error\n");
    free(old_sys);
    for (int i = 0; i < 4; i++) free(mat_arr[i]);
    free(mat);
    return 1;
  }
  memset(new_sys, 0, sys_len * sizeof(container_vec_t));

  size_t deg[4];
  for (int i = 0; i < 4; i++)
  {
    deg[i] = read_size_t();
  }

  container_t *correct_new_sys_arr[4];

  for (int i = 0; i < 4; i++)
  {
    correct_new_sys_arr[i] = read_container_t_array(sys_len);
  }

  container_vec_t *correct_new_sys =
      aligned_alloc(ALIGNMENT, sys_len * sizeof(container_vec_t));
  if (!correct_new_sys)
  {
    printf("Memory error\n");
    free(new_sys);
    free(old_sys);
    for (int i = 0; i < 4; i++)
    {
      free(mat_arr[i]);
      free(correct_new_sys_arr[i]);
    }
    free(mat);
    return 1;
  }

  for (size_t i = 0; i < sys_len; i++)
  {
    correct_new_sys[i] =
        VEC_ASSIGN(correct_new_sys_arr[0][i], correct_new_sys_arr[1][i],
                   correct_new_sys_arr[2][i], correct_new_sys_arr[3][i]);
  }

  srand(RSEED);

  container_vec_t res = compute_e_k(mat, new_sys, old_sys, l, n);

  uint8_t acc = 1;
  for (unsigned int i = 0; i < n_choose_k(n, 2) + n + 1; i++)
  {
    acc = acc && (VEC_M_MASK(VEC_EQ(new_sys[i], correct_new_sys[i])) == 0xF);
    if (!acc)
    {
      print_register(new_sys[i]);
      print_register(correct_new_sys[i]);
      break;
    }
  }

  container_t res_arr[4] = {VEC_EXTRACT(res, 3), VEC_EXTRACT(res, 2),
                            VEC_EXTRACT(res, 1), VEC_EXTRACT(res, 0)};
  for (int i = 0; i < 4; i++)
  {
    if (deg[i] != res_arr[i])
    {
      printf("Degrees do not correspond: %zu != %u\n", deg[i], res_arr[i]);
      acc = 0;
      break;
    }
  }

  free(mat);
  free(old_sys);
  free(new_sys);
  free(correct_new_sys);
  for (int i = 0; i < 4; i++)
  {
    free(mat_arr[i]);
    free(correct_new_sys_arr[i]);
  }

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

  container_t *sol_mat_arr[4] = {0};
  for (int i = 0; i < 4; i++)
  {
    sol_mat_arr[i] = read_container_t_array(l);
    if (!sol_mat_arr[i])
    {
      for (int j = 0; j < i; j++) free(sol_mat_arr[j]);
    }
  }

  container_vec_t *sol_mat =
      aligned_alloc(ALIGNMENT, l * sizeof(container_vec_t));
  if (!sol_mat)
  {
    for (int i = 0; i < 4; i++) free(sol_mat_arr[i]);
  }
  for (unsigned int i = 0; i < l; i++)
  {
    sol_mat[i] = VEC_ASSIGN(sol_mat_arr[0][i], sol_mat_arr[1][i],
                            sol_mat_arr[2][i], sol_mat_arr[3][i]);
  }

  container_vec_t *rand_mat = malloc(l * sizeof(container_vec_t));
  if (!rand_mat)
  {
    free(sol_mat);
    for (int i = 0; i < 4; i++) free(sol_mat_arr[i]);
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
    acc = acc && (VEC_M_MASK(VEC_EQ(sol_mat[i], rand_mat[i])) == 0xF);
    if (!acc)
    {
      print_register(sol_mat[i]);
      print_register(rand_mat[i]);
      break;
    }
  }

  for (int i = 0; i < 4; i++) free(sol_mat_arr[i]);
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

  container_t *system = read_container_t_array(sys_len);
  if (!system)
  {
    printf("Error receiving system\n");
    return 1;
  }

  container_t c_sol = 0;

  solve(system, n, m, &c_sol);

  if (eval(system, n, c_sol))
  {
    printf("Invalid solution found!\n");
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
    printf("Memory error\n");
    return 1;
  }

  container_t *old_sys = read_container_t_array(sys_len);
  if (!old_sys)
  {
    printf("Memory error\n");
    free(mat);
    return 1;
  }

  container_t *new_sys = calloc(sys_len, sizeof(container_t));
  if (!new_sys)
  {
    printf("Memory error\n");
    free(old_sys);
    free(mat);
    return 1;
  }

  unsigned int deg = read_size_t();

  container_t *correct_new_sys = read_container_t_array(sys_len);
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
    default:
      printf("Invalid function chosen: %u\n", choice);
  }
  return 1;
}
