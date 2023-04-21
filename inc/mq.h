#ifndef MQ_H
#define MQ_H

#include <stdlib.h>

#include "mq_config.h"

#if defined(_DEBUG)

#if defined(REG128) || defined(REG256)
container_vec_t compute_e_k(container_vec_t *mat, container_vec_t *new_sys,
                            container_t *old_sys, int l, int n);
#else
unsigned int compute_e_k(container_t *mat, container_t *new_sys,
                         container_t *old_sys, int l, int n);
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
uint8_t solve(container_t *system, unsigned int n, unsigned int m,
              container_t *sol);

#endif  // !MQ_H