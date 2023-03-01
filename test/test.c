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
  // poly_t sys[16] = {14, 14, 21, 19, 26, 8, 29, 26, 29, 8, 31, 16, 9, 24, 12,
  poly_t sys[] = {
      30283809, 14054022, 10022038, 20297807, 11209775, 7463323,  20648542,
      25492091, 12618052, 14565323, 23165076, 11718925, 23605779, 22923,
      25203230, 13042352, 8600174,  15232996, 1138999,  1611790,  15902015,
      3861400,  23852149, 33069226, 21621085, 5025516,  1769262,  19895433,
      21970887, 29319362, 7460719,  7149483,  20451334, 17443565, 3588897,
      21223311, 2043094,  7158547,  6689853,  32935382, 19333217, 24725927,
      33519861, 26713032, 33448385, 16010160, 20898162, 18814196, 1336503,
      4242620,  22182197, 7673223,  29204136, 30778583, 30403735, 28894258,
      4574397,  15142109, 13616636, 3202675,  27744209, 19636469, 29248634,
      14097025, 26731896, 504038,   2976335,  17209208, 25246510, 22670619,
      16015943, 6479456,  3116172,  1947285,  24079648, 4302171,  29036394,
      9007627,  19502605, 7176490,  23632073, 26731492, 8693975,  6873107,
      31028925, 27358382, 2319066,  24983373, 25735833, 30748,    7236011,
      14876646, 22527673, 9139284,  7116978,  32564999, 26003065, 29033011,
      20063006, 3805357,  13784021, 5332386,  4728007,  14945059, 14830132,
      26435456, 18736952, 15876644, 7502162,  12483983, 10233263, 22827759,
      12565820, 2736387,  3055366,  4689119,  13257065, 31948591, 30780869,
      27002583, 10154899, 21436095, 30815190, 2438216,  26948378, 25156389,
      18019907, 28290414, 21191232, 9639429,  4427007,  13409208, 7682349,
      10030645, 17267100, 12117617, 33288742, 7335989,  2743390,  19235762,
      26230961, 19028154, 12431018, 7180667,  10536801, 13483406, 25264408,
      11039793, 15394328, 9068638,  29054183, 28170861, 26776393, 21961968,
      8473488,  15976635, 25286433, 15417351, 21816202, 24498943, 23748634,
      21185443, 20078938, 32397909, 29680205, 28837504, 25319360, 1146278,
      33058732, 19746744, 14571116, 11024499, 17662707, 8983580,  8569881,
      14727245, 1845031,  19059306, 14306004, 3276825,  26264472, 20330556,
      8331077,  7774600,  28859716, 12923144, 33232819, 3720374,  6233727,
      2250895,  4546142,  14487069, 24187827, 27316078, 30631120, 32413931,
      9810124,  25200375, 15440176, 29779923, 3251512,  27240219, 22261090,
      9852264,  6381240,  32651661, 25601372, 8207003,  6774573,  17506005,
      4015742,  11857647, 32989120, 8882692,  3542200,  13271931, 28576525,
      22357685, 30169265, 17484649, 30848404, 10002567, 27174515, 6317373,
      12503626, 4730745,  25510322, 8123624,  16596420, 4044519,  30105825,
      12243095, 1300859,  32886497, 32708542, 21628063, 32237608, 20241908,
      22872972, 6795381,  9635666,  16622345, 32310650, 14602813, 11556692,
      14819158, 25004052, 16865472, 3541780,  2303394,  16162727, 6609522,
      11369886, 16436323, 16751074, 1403054,  15912334, 18259666, 29636348,
      31265160, 17097697, 4689007,  2990226,  1503438,  11269605, 22763412,
      3175719,  20839332, 20356901, 262235,   17807800, 17547325, 19062028,
      23629414, 16898899, 9746407,  20343645, 21125983, 31253202, 6640004,
      26846228, 32514761, 3966969,  389678,   19066486, 26950228, 2451590,
      26228519, 21698850, 20976676, 10543147, 21243358, 18197999, 27122912,
      24122227, 2548736,  31416582, 8512557,  9060790,  8620488,  27460883,
      26230116, 27886161, 2968726,  27411632, 13703953, 26069705, 31323913,
      22721487, 3900332,  19322470, 9388851,  12803600, 13647240, 21680315,
      4441548,  12806510, 6838856,  16669484, 19143560, 20216714, 11057374,
      21816039, 5443736,  11276120, 315475,
  };
  unsigned int n = 25;
  unsigned int m = 25;
  vars_t sol;
  printf("test?\n");
  uint8_t err = solve(sys, n, m, &sol);
  printf("Error code %u\n", err);

  return 0;
}
