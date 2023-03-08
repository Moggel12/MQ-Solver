#ifndef MQ_H
#define MQ_H

#include <stdlib.h>

#include "mq_config.h"

#if defined(_DEBUG)
unsigned int compute_e_k(poly_t *mat, poly_t *new_sys, poly_t *old_sys, int l,
                         int n);
#endif

/*!
 * An implementation of Dinur's polynomial method algorithm that uses an inverse
 * FES procedure to interpolate polynomials, instead of using the Möbius
 * transform.
 *
 * @param system A bitsliced representation of the system (assumed to be graded
 * lexicographic order without *pure* quadratic terms).
 * @param n An int defining the amount of variables in the system.
 * @param m An int defining the amount of polynomials in the system.
 * @param sol A pointer to where the function should store the solution found
 * (if no errors occurred).
 * @return a byte indicating whether or not the value in the sol parameter is a
 * valid solution. Returns 0 if no error occurred (solution is valid) and 1 if
 * an error occurred.
 */
uint8_t solve(poly_t *system, unsigned int n, unsigned int m, vars_t *sol);

#endif  // !MQ_H