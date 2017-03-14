#ifndef __dual_h__
#define __dual_h__

#include "../headfiles/matrix.h"

// Assumes p = inv * cost
// Returns the reduced cost of the ith variable
double get_reduced_cost(VECTOR *cost, VECTOR *P,  MAT **A, int i);


int dual_entering_index(MAT **B, MAT** A, VECTOR *cost, VECTOR *p, int *index);

// Receives a feasible model for a linear programming problem matrix, basic and non-basic indices
// basic inverse matrix, x vector, as well as cost vector and b vector all set accordingly
// Runs the dual simplex method
// Returns 1 if an optimal solution was found, -1 when an "infinite optimization" is found
int dual(SIMPLEX *model);

#endif
