#include "utils.h"
#include "fes.h"
#include <stdlib.h>
#include <stdio.h>

union tmp {
  uint32_t system[8];
  uint8_t system_8[32];
};

int main(void) {
  // printf("Test");
  printf("Test0\n");
  union tmp t;
  uint32_t system[8] = { 0x00007FFF, 0, 0x00007FFF, 0, 0x00007FFF, 0, 0x00007FFF, 0 };
  for (int i = 0; i < 8; i++) {
    t.system[i] = system[i];
  }
  printf(""BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(t.system_8[0]));
  printf(""BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(t.system_8[1]));
  printf(""BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(t.system_8[2]));
  printf(""BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(t.system_8[3]));
  int n = 4;
  int m = 8;
  printf("Test1\n");
  uint8_t *sliced = slice(t.system_8, n, m);
  if (sliced == NULL) {
    printf("Error\n");
    return 1;
  } 
  printf("Test2\n");
  int coeffs = n_choose_k(n + 1, 2) + n + 1; 
  printf("%d\n", coeffs);
  for (int i = 0; i < coeffs; i++) {
    printf("0b"BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(sliced[i]));
  }
  // free(sliced);
  return 0;
}
