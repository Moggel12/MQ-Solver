#ifndef FES_H
#define FES_H

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mq_config.h"
#include "utils.h"

#define INC(i) i + 1

#if defined(REG128) || defined(REG256)
/*! A struct identifying the state of of a partial evaluation of a polynomial
 * system. */
typedef struct state
{
  container_t i; /*! The value of the input variables assigned. */  // TODO
  container_vec_t y; /*! The evaluation of the system on s->i. */
  container_vec_t
      *d1; /*! The first derivatives of the system evaluated on s->i. */
  container_vec_t *d2; /*! The second derivative of the system (constants). */
  uint8_t *prefix;     /*! The prefix used for the partial evaluation. */
} state;

// TODO: Update documentation of fes_recover
/*!
 * Computes the evaluations of all 2^(n - n1) inputs for the U_i polynomials in
 * Dinur's polynomial-method algorithm. The resulting evaluations are in a
bitsliced format.
 *
 * @param system The system related to the U_i polynomials in bitsliced format.
 * @param n The amount of variables in the system.
 * @param n1 The parameter denoting the (n - n1) y-bits and (n1) z-bits.
 * @param deg The degree of the identifying polynomial of *system*.
 * @param results The table of all evaluations and their inputs.
 * @return Returns 1 if an error occurred, else 0. If 1 is returned, the values
of *resuslts* are invalid.
 */
uint8_t fes_recover_vectorized(container_t *system,
                               container_vec_t *e_k_systems, unsigned int n,
                               unsigned int n1, container_vec_t deg,
                               container_t *result);

#else

typedef struct PotentialSolution
{
  container_t y_idx;
  container_t z_bits;
} PotentialSolution;

/*! A struct identifying the state of of a partial evaluation of a polynomial
 * system. */
typedef struct state
{
  container_t i; /*! The value of the input variables assigned. */  // TODO
  container_t y;   /*! The evaluation of the system on s->i. */
  container_t *d1; /*! The first derivatives of the system evaluated on s->i. */
  container_t *d2; /*! The second derivative of the system (constants). */
  uint8_t *prefix; /*! The prefix used for the partial evaluation. */
} state;

/*!
 * This does not (yet) conform to the interface given through mq_config.h, and
 * therefore does not support non integer polynomials.
 *
 * The bruteforce function splits the inputs to the system into a prefix of (n -
 * n1) bits and filters out solutions where the hamming weight is >= d.
 *
 * @param system The system to bruteforce in a bitsliced format.
 * @param n The amount of variables in *system*.
 * @param n1 The value used to denote where the input for *system* is split into
 * prefix and remainder.
 * @param d The value used to filter which prefixes should be included.
 * @param m The amount of polynomials in the system.
 * @return Returns the amount of solutions of 0xFF..FF if an error occurred.
 */
unsigned int bruteforce(container_t *system, unsigned int n, unsigned int n1,
                        unsigned int d, container_t *solutions);

/*!
 * Barebones FES implementation for comparing against Dinur's algorith (this is
 * essentially used internally on Dinur's algorithm).
 *
 * @param system The system (in a bitsliced format) to find a solution for.
 * @param n The amount of variables in *system*.
 * @param m The amount of polynomials in the system.
 * @param solutions An array of which the procedure should store the solutions.
 * The procedure expects there to be enough room for all solutions.
 * @return Returns the amount of solutions found.
 */
unsigned int fes(container_t *system, unsigned int n, unsigned int m,
                 container_t *solutions);

// TODO: Update documentation of fes_recover
/*!
 * Computes the evaluations of all 2^(n - n1) inputs for the U_i polynomials in
 * Dinur's polynomial-method algorithm. The resulting evaluations are in a
bitsliced format.
 *
 * @param system The system related to the U_i polynomials in bitsliced format.
 * @param n The amount of variables in the system.
 * @param n1 The parameter denoting the (n - n1) y-bits and (n1) z-bits.
 * @param deg The degree of the identifying polynomial of *system*.
 * @param results The table of all evaluations and their inputs.
 * @return Returns 1 if an error occurred, else 0. If 1 is returned, the values
of *resuslts* are invalid.
 */
uint8_t fes_recover(container_t *system, unsigned int n, unsigned int n1,
                    unsigned int deg, PotentialSolution *results,
                    size_t *res_size);

#endif

#endif  // FES_H
