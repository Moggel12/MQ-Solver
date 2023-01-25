#ifndef FES_H
#define FES_H

#include "utils.h"
#include "mq_config.h"
#include <stdint.h>
#include <stdlib.h>


typedef struct state {
  vars_t i;
  poly_t y;
  poly_t *d1;
  poly_t *d2;
  uint8_t *prefix;
} state;

/*!
 *
 * @param system
 * @param n
 * @param n1
 * @param d
 * @param m
 * @return
 */
unsigned int bruteforce(poly_t *system, unsigned int n, unsigned int n1, unsigned int d, vars_t *solutions);

/*!
 *
 * @param system
 * @param n
 * @param n1
 * @param prefix
 * @param s
 * @param solutions
 * @param sol_amount
 */
void fes_eval_solutions(poly_t *system, unsigned int n, unsigned int n1, uint8_t *prefix, state *s, vars_t *solutions, unsigned int *sol_amount);

/*!
 *
 * @param system
 * @param n
 * @param n1
 * @param prefix
 * @param s
 * @return
 */
vars_t fes_eval_parity(poly_t *system, unsigned int n, unsigned int n1, uint8_t *prefix, state *s);

#endif // FES_H
