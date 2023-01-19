#include "fes.h"
#include "utils.h"
#include <string.h>
#include <limits.h>

#include <stdio.h>

state *init_state(unsigned int n1, uint8_t *prefix) {
  state *s = malloc(sizeof(state));

  if (!s) return NULL;

  s->d1 = calloc(n1, sizeof(vars_t));

  if (!(s->d1)) {
    free(s);

    return NULL;
  }

  if (!prefix) {
    s->prefix = calloc(n1, sizeof(s->prefix));

    if (!(s->prefix)) {
      free(s->d1);
      free(s);

      return NULL;
    }

  } else {
    s->prefix = prefix;
  }

  s->d2 = calloc(n1*n1, sizeof(vars_t));

  if (!(s->d2)) {
    if (!prefix) free(s->prefix);

    free(s->d1);
    free(s);

    return NULL;
  }

  s->i = 0;
  s->y = 0;

  return s;
}

void destroy_state(state *s) {
  free(s->d1);
  free(s->d2);
  free(s->prefix);
  free(s);
}

unsigned int bit1(vars_t i) { return trailing_zeros(i); }

unsigned int bit2(vars_t i) { return bit1(GF2_ADD(i, (i & -i))); }

state *init(state *s, poly_t *system, unsigned int n, unsigned int n1, uint8_t *prefix) {

  s = init_state(n1, prefix);

  if (!s) {
    return NULL;
  }

  s->y = system[0];

  for (unsigned int k = 0; k < n1; k++) {

    for (unsigned int j = 0; j < k; j++) {
      s->d2[k * n1 + j] = system[lex_idx(j + (n - n1), k + (n - n1), n)]; // TODO
    }

  }

  s->d1[0] = system[1 + n - n1];


  for (unsigned int k = 1; k < n1; k++) {
    s->d1[k] = GF2_ADD(s->d2[k * n1 + (k-1)], system[1 + k + (n - n1)]);
  }

  for (unsigned int i = 0; i < n1; i++) {

    if (prefix[i] == 0) continue; // TODO: Rethink how prefix is represented

    for (unsigned int k = 0; k < n1; k++) {
      s->d1[k] = GF2_ADD(s->d1[k], system[lex_idx(i, k + (n - n1), n)]);
    }

    s->y = GF2_ADD(s->y, system[i + 1]);
  }

  // 2-combs
  for (unsigned int i = 0; i < n1; i++) {

    if (prefix[i] == 0) continue; // TODO: Rethink how prefix is represented

    for (unsigned int j = i + 1; j < n1; j++) {
      if (prefix[j] == 0) continue; // TODO: Rethink how prefix is represented

      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
    }
  }
  return s;
}

// FREES prefix FROM STATE.
state *update(state *s, poly_t *system, unsigned int n, unsigned int n1, uint8_t *prefix) {
  uint8_t alloc_state = 0;

  if (!s) {
    s = init(s, system, n, n1, prefix);
    if (!s) {
      return NULL;
    }
    return s;
  }

  uint8_t *on = malloc(n1);
  if (!on) {
    return NULL;
  }
  uint8_t *off = malloc(n1);
  if (!off) {
    free(on);
    return NULL;
  }

  for (unsigned int i = 0; i < n1; i++) {
    off[i] = (prefix[i] == 0) && (s->prefix[i] == 1) ? 1 : 0;
    on[i] = (s->prefix[i] == 0) != (prefix[i] == 1) ? 1 : 0;
  }

  // Turn off variables not assigned in new prefix
  for (unsigned int idx = 0; idx < n1; idx++) {
    if (off[idx] == 0) continue; // TODO: Rethink how off is represented

    for (unsigned int k = 0; k < n1; k++) {
      s->d1[k] = GF2_ADD(s->d1[k], system[lex_idx(idx, k + (n - n1), n)]);
    }

    s->y = GF2_ADD(s->y, system[idx + 1]);
  }

  for (int i = 0; i < n1; ++i) {
    if (off[i] == 0) continue; // TODO: Rethink how off is represented

    for (unsigned int j = 0; j < n1; j++) {
      if ((s->prefix[j] == 1) && (off[j] == 1)) continue; // TODO: Rethink how prefix and off is represented

      if (i > j) { // Switch order for lex_idx
        int tmp = i;
        i = j;
        j = tmp;
      }
      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
    }
  }

  // 2-combs
  for (unsigned int i = 0; i < n1; i++) {
    if (off[i] == 0) continue; // TODO: Rethink how off is represented

    for (unsigned int j = i + 1; j < n1; j++) {
      if (off[j] == 0) continue; // TODO: Rethink how off is represented

      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
    }
  }

  // Turn new variables on
  for (unsigned int idx = 0; idx < n1; idx++) {
    if (on[idx] == 0) continue; // TODO: Rethink how on is represented

    for (unsigned int k = 0; k < n1; idx++) {
      s->d1[k] = GF2_ADD(s->d1[k], system[lex_idx(idx, k + (n - n1), n)]);
    }
  }

  for (unsigned int i = 0; i < n1; i++) {
    if (on[i] == 0) continue; // TODO: Rethink how on is represented

    for (unsigned int j = 0;j < n1; j++) {
      if ((prefix[j] == 1) && (on[j] == 1)) continue; // TODO: Rethink how prefix and on is represented

      if (i > j) { // Switch order for lex_idx
        int tmp = i;
        i = j;
        j = tmp;
      }
      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
    }
  }

  // 2-combs
  for (unsigned int i = 0; i < n1; i++) {
    if (on[i] == 0) continue; // TODO: Rethink how on is represented

    for (unsigned int j = i + 1; j < n1; j++) {
      if (on[j] == 0) continue; // TODO: Rethink how on is represented

      s->y = GF2_ADD(s->y, system[lex_idx(i, j, n)]);
    }
  }

  free(s->prefix);
  s->prefix = prefix;
  free(off);
  free(on);

  return s;
}

static void step(state *s, unsigned int n1) {
  s->i = INC(s->i);

  unsigned int k1 = bit1(s->i);
  unsigned int k2 = bit2(s->i);

  if (k2 < CHAR_BIT * sizeof(s->i)) {
    s->d1[k1] = GF2_ADD(s->d1[k1], s->d2[k2 * n1 + k1]);
  }

  s->y = GF2_ADD(s->y, s->d1[k1]);

}

vars_t fes_eval_parity(poly_t *system, unsigned int n, unsigned int n1, uint8_t *prefix, state *s) {

  if (!s) {
    init(s, system, n, n1, prefix);
    if (!s) {
      return -1; // 0xFF..FF
    }
  }

  uint64_t pre_x = 0;
  for (unsigned int i = 0; i < (n - n1); i++) {
    if (prefix[i] == 0) continue; // TODO: Change representation of prefixes

    pre_x += (1 << i);
  }

  vars_t parities = 0;

  if (s->y == 0) {
    parities = GF2_ADD(s->y, ((1 << n1) -1));
  }

  while (s->i < ((1 << n1) - 1)) {
    step(s, n1);

    if (s->y == 0) {
      parities = GF2_ADD(parities, 1);

      for (unsigned int pos = 0; pos < n1; pos++) {
        parities = GF2_ADD(parities, (1 << (pos + 1)));
      }
    }
  }

  s->i = 0;
  for (unsigned int i = 0; i < (n1 - 1); i++) {
    s->d1[i] = GF2_ADD(s->d1[i], s->d2[(n1 - 1) * n1 + i]);
  }

  s->y = GF2_ADD(s->y, GF2_ADD(s->d1[n1 - 1], s->d2[(n1 - 1) * n1 + n1 - 2]));

  return parities;
}

void fes_eval_solutions(poly_t *system, unsigned int n, unsigned int n1, uint8_t *prefix, state *s, vars_t *solutions, unsigned int *sol_amount) {
  if (!s) {
    s = init(s, system, n, n1, prefix);

    if (!s) {
      return;
    }
  }
  *sol_amount = 0;

  uint64_t pre_x = 0;
  for (unsigned int i = 0; i < n1; i++) {
    if (prefix[i] == 0) continue;

    pre_x += (1 << i);
  }

  if (s->y == 0) {
    solutions[(*sol_amount)++] = ((s->i ^ (s->i >> 1)) << (n - n1) | pre_x);
  }

  while (s->i < ((1 << n1) - 1)) {
    step(s, n1);

    if (s->y == 0) {
      solutions[(*sol_amount)++] = ((s->i ^ (s->i >> 1)) << (n - n1) | pre_x);
    }
  }

  s->i = 0;
  for (unsigned int i = 0; i < (n1 - 1); i++) {
    s->d1[i] = GF2_ADD(s->d1[i], s->d2[(n1 - 1) * n1 + i]);
  }

  s->y = GF2_ADD(s->y, GF2_ADD(s->d1[n1 - 1], s->d2[(n1 - 1) * n1 + (n1 - 2)]));
}

// Expects system pre-sliced
unsigned int bruteforce(poly_t *system, vars_t n, unsigned int n1, unsigned int d, unsigned int m, vars_t *solutions) {
  unsigned int sol_amount = 0;

  uint8_t *prefix = malloc(n1);
  if (!prefix) {
    free(solutions);
    return -1; // 0xFF..FF
  }

  state *s = NULL;

  for (unsigned int i = 0; i < (1 << (n - n1)); i++) {
    if (hamming_weight(i) <= d) {
      for (unsigned int pos = 0; pos < n1; pos++) {
        prefix[pos] = (1 & (i >> pos));
      }

      s = update(s, system, n, n1, prefix);
      fes_eval_solutions(system, n, n1, prefix, s, solutions, &sol_amount);
    }
  }

  destroy_state(s);
  return sol_amount;
}

void test(void) {
  poly_t polynomials[16] = { 5, 5, 25, 1, 7, 12, 17, 1, 28, 29, 2, 16, 21, 15, 21, 30};

  unsigned int n = 5;
  unsigned int sol_amount = 0;
  vars_t solutions[32] = {0};

  state *s = NULL;

  unsigned int n1 = 3;
  uint8_t prefix[3] = {1, 0, 0};

  printf("Initialized values...\n");

  s = update(s, polynomials, n, n1, prefix);
  fes_eval_solutions(polynomials, n, n1, prefix, s, solutions, &sol_amount);

  printf("Finished solver.\n");
  printf("Amount of solutions: %u\n", sol_amount);

  if (s) {
    for (unsigned int i = 0; i < sol_amount; i++) {
      printf("Solution %u\n", solutions[i]);
    }
  }
}
