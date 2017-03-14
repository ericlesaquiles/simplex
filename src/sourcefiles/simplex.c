#include "../headfiles/simplex.h"
#include "../headfiles/interface.h"
#include "../headfiles/matrix.h"

double reduced_cost(SIMPLEX *model, VECTOR *p, int i){
  return access_vector(model->cost, model->non_basic_index[i]) - mat_inner_product(model->matrix, p, model->non_basic_index[i]);

}

int get_reduced_cost_index(SIMPLEX *model, VECTOR *p){
  double c;
  int i;

  for(i = 0; i < model->matrix->col_qtd - model->basis_size; i++){
    /*c = access_vector(model->cost, model->non_basic_index[i]) - mat_inner_product(model->matrix, p, model->non_basic_index[i]);*/
    c = reduced_cost(model, p, i);
    if(c < 0){
      return i;
    }
  }
  return -1;
}

// Receives a feasible model for a linear programming problem matrix, basic and non-basic indices
// basic inverse matrix, x vector, as well as cost vector and b vector all set accordingly
// Runs the primal simplex method
// Returns 1 if an optimal solution was found, -1 when an "infinite optimization" is found
int simplex(SIMPLEX *model){

  int reduced_cost_index = 1, i, k;
  double min, j;

  VECTOR *p, *u;
  p = make_vector(model->basis_size);
  u = make_vector(model->inv->row_qtd);

  p->index = model->inv->col_index;
  u->index = model->inv->row_index;

  initialize_vector(p, nil_value);
  initialize_vector(u, nil_value);

  while(reduced_cost_index != -1){

    // Multiplies the inverse by the cost vector on the basic index and puts the result on p (that is: p = cost * inv)
    left_vet_mut(model->inv, model->cost, p, model->basic_index, NULL);

    reduced_cost_index = get_reduced_cost_index(model, p);
    model->actual_cost = inner_product(model->cost, model->x);

    if(reduced_cost_index == -1){
        destroy_vector(p);
        destroy_vector(u);
        return 1;
    } else {
      mat_right_vet_mut(model->inv, model->matrix, u, model->non_basic_index[reduced_cost_index]);

      for(i = 0; i < u->size; i++){
        if(u->vector[i] > 0) break;
      }
      // infinite optimization
      if(i == u->size){
        destroy_vector(p);
        destroy_vector(u);
        return -1;
      }

      // Finds minimum of ( Xb(i)/U(i) )
      // First, initializes min
      min = access_vector(model->x, model->basic_index[i])/u->vector[i];
      k = i;
      for(i = 0, j = min; i < model->basis_size; ++i){
        if(u->vector[i] > 0) j = access_vector(model->x, model->basic_index[i])/u->vector[i];

        if(j < min){
          min = j;
          k = i;
        }
      }

      // Exchanges the one basic variable leaving for the one entering
      // But, not before setting apropriate variable to zero
      access_vector(model->x, model->basic_index[k]) = 0;
      swap_i(&model->basic_index[k], &model->non_basic_index[reduced_cost_index]);
      access_vector(model->x, model->basic_index[k]) = min;

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
  destroy_vector(u);
  destroy_vector(p);
  return 0;
}

void destroy_model(SIMPLEX *model){
  int *index_r, *index_c;

  index_r = model->matrix->row_index;
  index_c = model->matrix->col_index;

  if(model->matrix != NULL) destroy_matrix(model->matrix);

  free(model->inv->col_index);
  if(model->inv != NULL)    destroy_matrix(model->inv);
  if(model->UL != NULL)     destroy_matrix(model->UL);

  if(model->y != NULL)    destroy_vector(model->y);
  if(model->x != NULL)    destroy_vector(model->x);
  if(model->b != NULL)    destroy_vector(model->b);
  if(model->cost != NULL) destroy_vector(model->cost);

  if(model->splitted_var != NULL) free(model->splitted_var);

  free(model->basic_index);
  free(model->non_basic_index);

  if(model != NULL) free(model);

  free(index_r);
  free(index_c);
}
