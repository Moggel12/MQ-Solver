#ifndef STANDARD_UTILS_H
#define STANDARD_UTILS_H

#include <stddef.h>
#include <stdint.h>

#include "mq_config.h"

unsigned int gen_matrix(container_t *mat, unsigned int n_rows,
                        unsigned int n_columns);

container_t parity(container_t bits);

container_t eval(container_t *system, size_t n, container_t var_values);

#endif  // !STANDARD_UTILS_H