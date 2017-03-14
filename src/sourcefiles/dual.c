#include "../headfiles/interface.h"
#include "../headfiles/matrix.h"
/*#include "dev.c"*/

#define abs(a)\
   a < 0 ? -a : a

// Assumes p = inv * cost
// Returns the reduced cost of the ith variable
/*double reduced_cost(SIMPLEX *model, VECTOR *p, int i){*/
  /*return access_vector(model->cost, model->non_basic_index[i]) - mat_inner_product(model->matrix, p, model->non_basic_index[i]);*/
/*}*/

int dual_entering_index(SIMPLEX *model, VECTOR *p, int l, int *index){
  int i, j, k, r, m = -1, e = -1;
  double c;

  for(i = 0; i < model->inv->col_qtd; i++){
    /*for(j = 0; j < A->row_qtd; j++){*/
      /*r += access_matrix(j, i, A) * access_matrix(i, j, B);*/
    /*}*/

    r = abs(ith_line_lth_column_mut(model->inv, model->matrix, l, model->non_basic_index[i]));

    if(r < 0){
      c = reduced_cost(model, p, i);
      if(m > c/r || m == -1){
        m = c/r;
        e = i;
      }
    }
  }
  *index = e;

  return m;
}


// Receives a feasible model for a linear programming problem matrix, basic and non-basic indices
// basic inverse matrix, x vector, as well as cost vector and b vector all set accordingly
// Runs the dual simplex method
// Returns 1 if an optimal solution was found, -1 when an "infinite optimization" is found
int dual(SIMPLEX *model){

  int reduced_cost_index = 1, i, k;
  double min, j;

  VECTOR *p, *u;
  p = make_vector(model->basis_size);
  u = make_vector(model->matrix->col_qtd);

  p->index = model->inv->col_index;
  u->index = model->matrix->col_index;

  initialize_vector(p, nil_value);
  initialize_vector(u, nil_value);

  while(reduced_cost_index != -1){
    // If all the xs are non-negative, we have an optimal solutino
    for(k = 0; k < model->x->size; k++) if(access_vector(model->x, k) < 0) break;
    if(k == model->x->size - 1 && access_vector(model->x, k) >= 0) return 1;

    // Multiplies the inverse by the cost vector on the basic index and puts the result on p
    left_vet_mut(model->inv, model->cost, p, model->basic_index, NULL);
    min = dual_entering_index(model, p, k, &reduced_cost_index);

    if(reduced_cost_index == -1){
        destroy_vector(p);
        destroy_vector(u);
        return -1;
    } else {

      // Exchanges the one basic variable leaving for the one entering
      // But, not before setting apropriate variable to zero
      access_vector(model->x, model->basic_index[k]) = 0;
      swap_i(&model->basic_index[k], &model->non_basic_index[reduced_cost_index]);
      access_vector(model->x, model->basic_index[k]) = min;

      mat_right_vet_mut(model->inv, model->matrix, u, model->non_basic_index[reduced_cost_index]);

      // Update basic variables values
      for(i = 0; i < k;  i++){
        access_vector(model->x, model->basic_index[i]) = access_vector(model->x, model->basic_index[i]) - min * u->vector[i];
      }
      for(i++; i < model->basis_size;  i++){
        access_vector(model->x, model->basic_index[i]) = access_vector(model->x, model->basic_index[i]) - min * u->vector[i];
      }

      // Form new inverse matrix
      // sum -(row k of the inverse / u[k] * u[i]) to rows i
      for(i = 0; i != k && i < model->inv->row_qtd; i++){
        row_operation(model->inv, k, i, -(u->vector[i]/u->vector[k]));
      }
      for(i++; i != k && i < model->inv->row_qtd; i++){
        row_operation(model->inv, k, i, -(u->vector[i]/u->vector[k]));
      }
      for(i = 0; i < model->basis_size; i++){
        access_matrix(k, i, model->inv) /= u->vector[k];
      }
    }
  }
  destroy_vector(p);
  destroy_vector(u);
  return 0;
}
