#ifndef FES_H
#define FES_H

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mq_config.h"
#include "utils.h"

/*! A struct identifying the state of of a partial evaluation of a polynomial
 * system. */
typedef struct state
{
  vars_t i;        /*! The value of the input variables assigned. */
  poly_t y;        /*! The evaluation of the system on s->i. */
  poly_t *d1;      /*! The first derivatives of the system evaluated on s->i. */
  poly_t *d2;      /*! The second derivative of the system (constants). */
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
unsigned int bruteforce(poly_t *system, unsigned int n, unsigned int n1,
                        unsigned int d, vars_t *solutions);

/*!
 * This does not (yet) conform to the interface given through mq_config.h, and
 * therefore does not support non integer polynomials.
 *
 * @param system The system to exhaustively search in bitsliced format.
 * @param n The amount of input variables to *system*.
 * @param n1 The parameter denoting where the prefix and remainder of the input
 * is split.
 * @param prefix The prefix of variables for which *system* should be partially
 * evaluated.
 * @param s The current *state*, i.e. the old prefix, etc.
 * @param solutions A list of all solutions found. The solutions are stored from
 * solutions[0] to solutions[sol_amount - 1].
 * @param sol_amount The amount of solutions found.
 */
void fes_eval_solutions(poly_t *system, unsigned int n, unsigned int n1,
                        uint8_t *prefix, state *s, vars_t *solutions,
                        unsigned int *sol_amount);

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
uint8_t fes_recover(poly_t *system, unsigned int n, unsigned int n1,
                    unsigned int deg, vars_t *results);

#endif  // FES_H
