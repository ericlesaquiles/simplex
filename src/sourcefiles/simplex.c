#include "../headfiles/simplex.h"
#include "../headfiles/interface.h"
#include "../headfiles/matrix.h"
#include "dev.c"


int get_reduced_cost_index(SIMPLEX *model, VECTOR *p){
  double c;
  int i;

  for(i = 0; i < model->matrix->col_qtd - model->basis_size; i++){
    c = access_vector(model->cost, model->non_basic_index[i]) - mat_inner_product(model->matrix, p, model->non_basic_index[i]);
    if(TESTE) printf(" Custo reduzido %d: %lf\n\n", i, c);
    if(c < 0){
      return i;
    }
  }
  return -1;
}

// Resolve um problema de programacao linear pelo algoritmo SIMPLEX, a partir do modelo presente em model
// Retorna 1 quando encontra uma solucao otima, -1 quando encontra "otimizacao infinita"
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
    print_all(model);

    // Multiplies the inverse by the cost vector on the basic index and puts the result on p
    left_vet_mut(model->inv, model->cost, p, model->basic_index, NULL);

    printv("Vetor p: ", p);

    reduced_cost_index = get_reduced_cost_index(model, p);
    model->actual_cost = inner_product(model->cost, model->x);

    if(TESTE) printf(" Custo atual: %.2lf\n\n", model->actual_cost);
    if(TESTE) printf(" Entra na base: %d\n\n", reduced_cost_index);

    if(reduced_cost_index == -1){
        destroy_vector(p);
        destroy_vector(u);
        return 1;
    } else {
      mat_right_vet_mut(model->inv, model->matrix, u, model->non_basic_index[reduced_cost_index]);

      printv("Vetor u: ", u);

      for(i = 0; i < u->size; i++){
        if(u->vector[i] > 0) break;
      }
      // infinite optimization
      if( i == u->size ){
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
          if(TESTE) printf(" Variavel a sair: %d \n\n", model->basic_index[k]);
          if(TESTE) printf(" Theta minimo: %lf \n\n", j);
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
