#ifndef MQ_H
#define MQ_H

#include <stdlib.h>

#include "mq_config.h"

#if defined(_DEBUG)

#if defined(REG128) || defined(REG256)
poly_vec_t compute_p_k(poly_t *mat, poly_vec_t *new_sys, poly_t *old_sys,
                       size_t sys_len, int l, int n);

void fix_poly(poly_t *system, poly_t *fixed_system, poly_t *assignment,
              int new_n, int n);
#else
unsigned int compute_p_k(poly_t *mat, poly_t *new_sys, poly_t *old_sys, int l,
                         int n);

extern size_t solver_rounds;
#endif

#endif

/*!
 * An implementation of Dinur's polynomial method algorithm that uses an
 * inverse FES procedure to interpolate polynomials, instead of using the
 * Möbius transform.
 *
 * @param system A bitsliced representation of the system (assumed to be
 * graded lexicographic order without *pure* quadratic terms).
 * @param n An int defining the amount of variables in the system.
 * @param m An int defining the amount of polynomials in the system.
 * @param sol A pointer to where the function should store the solution
 * found (if no errors occurred).
 * @return a byte indicating whether or not the value in the sol parameter
 * is a valid solution. Returns 0 if no error occurred (solution is valid)
 * and 1 if an error occurred.
 */
uint8_t solve(poly_t *system, unsigned int n, unsigned int m, poly_t *sol);

poly_t eval(poly_t *system, size_t n, poly_t var_values);

#endif  // !MQ_H