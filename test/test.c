#include "mq_config.h"
#include "fes.h"

int main(void) {
  uint8_t system[] = {38, 124, 50, 76, 20, 5, 97, 118, 69, 49, 49, 63, 127, 61, 5, 20, 80, 11, 24, 62, 70, 60, 11, 75, 61, 49, 64, 126, 82};
  unsigned int n = 7;
  unsigned int n1 = 2;
  unsigned int d = 1;
  vars_t results[1 << 5] = {0};
  fes_recover(system, n, n1, d, results);
  return 0;
}
