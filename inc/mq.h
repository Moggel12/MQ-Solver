#ifndef MQ_H
#define MQ_H

#include <stdlib.h>

#include "mq_config.h"

#if defined(_DEBUG)
unsigned int compute_e_k(poly_t *mat, poly_t *new_sys, poly_t *old_sys, int l,
                         int n);
#endif

uint8_t solve(poly_t *system, unsigned int n, unsigned int m, vars_t *sol);

#endif  // !MQ_H