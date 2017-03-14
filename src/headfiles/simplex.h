#ifndef __simplex_h
#define __simplex_h

#include "matrix.h"

/*
 *
 *  P: Permutation matrix
 *  L: Lower triangular matrix
 *  U: Upper triangular matrix
 *
 *  PA = UL
 *  PA = Pb
 *  UL = Pb
 *
 *  y = Ux
 *
 *  Solve Ly = Pb by forward substitution
 *  Solve Ux = y  by backward elimination
 *
 *
 * Variables:
 *  int *splitted_var := takes the unbounded variables, to be splitted
 *
 * */

typedef struct simplex{
  MAT *matrix;
  MAT *inv;
  MAT *UL;

  VECTOR *y;
  VECTOR *x;
  VECTOR *b;
  VECTOR *cost;

  int *basic_index;
  int *non_basic_index;
  int *splitted_var;

  int basis_size;
  int new_cols;

  double actual_cost;
} SIMPLEX;


double reduced_cost(SIMPLEX *model, VECTOR *p, int i);

int get_reduced_cost_index(SIMPLEX *model, VECTOR *p);
int simplex(SIMPLEX *model);

SIMPLEX *begin_simplex(char *nome);

void destroy_model(SIMPLEX *model);

#endif
