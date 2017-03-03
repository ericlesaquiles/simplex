#include  "../headfiles/interface.h"
#include "../headfiles/matrix.h"
#include <ctype.h>

// Receives a simplex model and a string of symbols (=, < or > symbols) as parameters
// Turns the model into standard form
int phase_I(SIMPLEX *model, char *symbols){
  int i, j;

  // If the constraints are like Ax + s = b, with s != 0, s = b is an initial basic feasible solution
  // And that is the case when the quantity of added slack variables is equal to the one existing before
  if(model->matrix->row_qtd == model->new_cols){
    model->basic_index      =  malloc(sizeof(int) * model->new_cols);
    model->non_basic_index  =  malloc(sizeof(int) * model->matrix->col_qtd - model->new_cols);

    for(i = 0; i < model->matrix->col_qtd - model->new_cols; i++) model->non_basic_index[i] = i;
    for(j = 0; i < model->matrix->col_qtd; j++, i++) model->basic_index[j] = i;

    model->basis_size = model->new_cols;

    model->inv = make_matrix(model->basis_size, model->basis_size);
    model->inv->col_index = (int*)malloc(sizeof(int) * model->basis_size);
    for(i = 0; i < model->inv->col_qtd; i++) model->inv->col_index[i] = i;

    model->UL  = make_matrix(model->basis_size, model->basis_size);

  } else {
    // Let's then first create the artificial problem : min y'x s.t. Ax + y = b

    // Store the current cost vector into 'store'
    VECTOR *store = NULL;
    copy_vector(model->cost, &store);


    int dif_length = model->matrix->row_qtd;
    int flag = 1;

    // Then, update the current cost vector so that it has all current coefficients 0, and plus row_qtd - new_cols coefficients of 1
    resize_vector(model->cost, dif_length);
    model->cost->index = realloc(model->cost->index, sizeof(int) * (model->matrix->col_qtd + dif_length));

    for(i = 0; i < model->matrix->row_qtd; i++) access_vector(model->cost, i) = 0;
    for(i = model->matrix->col_qtd; i < model->matrix->col_qtd + dif_length; i++){
      model->cost->index[i] = i;
      access_vector(model->cost, i) = 1;
    }

    // Of course, it is also needed to update the A matrix, as well as the x and y vector
    model->matrix->col_index = model->cost->index;
    model->y->index = model->cost->index;
    model->x->index = model->cost->index;

    resize_matrix(model->matrix, 0, dif_length);
    resize_vector(model->x, dif_length);

    // Set the matrix
    for(i = 0;  i < dif_length; i++){
      for(j = 0; j < model->matrix->row_qtd; j++){
        if(access_vector(model->b, j) < 0) flag = -1;
        else flag = 1;
        access_matrix(j, i + model->matrix->col_qtd - dif_length, model->matrix) = i == j ? flag :  0;
      }
    }

    // Finish setting up the UL and inv matrices, as well as support vectors
    model->basic_index      =  malloc(sizeof(int) * model->matrix->row_qtd);
    model->non_basic_index  =  malloc(sizeof(int) * model->matrix->col_qtd - model->matrix->row_qtd);

    for(i = 0; i < model->matrix->col_qtd - model->matrix->row_qtd; i++){
      model->non_basic_index[i] = i;
    }
    for(j = 0; j < model->matrix->row_qtd; i++, j++){
      model->basic_index[j] = i;
    }

    model->basis_size = model->matrix->row_qtd;

    model->inv = make_matrix(model->basis_size, model->basis_size);
    model->inv->col_index = (int*)malloc(sizeof(int) * model->basis_size);
    for(i = 0; i < model->inv->col_qtd; i++) model->inv->col_index[i] = i;

    model->UL = make_matrix(model->basis_size, model->basis_size);

    // Well then, let's solve for it
    set_up_model_env(model);
    flag = simplex(model);

    if(flag == -1 || model->actual_cost > 0) exit(1); // problem is unfeasible
    else{
      // All problems lie on the basis. Let's check it, element by element
      for(i = 0; i < model->basis_size; i++){
        // Checks whether the ith basic index is artificial
        if(model->matrix->col_index[model->basic_index[i]] >= model->matrix->col_qtd - dif_length){
          // If it is, check the ith entry of inv(B) * Aj. If all of them are zero, the ith row is redundant
          // and can be eliminated
          for(j = 0; j < model->matrix->col_qtd - model->basis_size; j++){
            mat_left_vet_mut(model->inv, model->matrix, model->y, model->non_basic_index[j], NULL, NULL);

            if(access_vector(model->y, i)){
              // Column j is LI with the rest of the basis minus the artificial element there, so one can change basis
              flag = model->basic_index[i];
              model->basic_index[i] = model->non_basic_index[j];
              model->non_basic_index[j] = flag;

              flag = 1;
              break;
            }
          }
          // All the ith variables are at zero level, the the ith row is redundant and is to be dropped
          if(!flag) drop_row(model->matrix, i);
          i--;
        }
      }

      // REFACTOR
      model->matrix->col_qtd -= dif_length;
      // Drop the artificial variables and clean the mess
      for(i = 0; i < model->matrix->row_qtd; i++){
        // Remember to review this!
        // Push each non_basic variable that pointed to an artificial variable to the edge of the vector, so that it becomes unaccessible
        for(j = 0, flag = model->matrix->col_qtd - model->basis_size; j < model->matrix->col_qtd - model->basis_size; j++){
          // Checks whether the non_basic variables lie in the range
          if(model->non_basic_index[j] >= model->matrix->col_qtd){
            while(model->non_basic_index[flag] >= model->matrix->col_qtd) flag++;
            swap_i(&(model->non_basic_index[j]), &(model->non_basic_index[flag]));
            flag++;
          }
        }
      }
    }

    // Recover cost vector
    free(model->cost->vector);
    free(model->cost);
    model->cost = store;
    model->cost->index = model->matrix->col_index;

    free(symbols);
    return 0;
  }
  free(symbols);
  return 1;
}

// Receives a simplex model as parameter
// Creates the lup factorization, the inverse, the initial solution and its cost
void set_up_model_env(SIMPLEX *model){
  lup_factor(model->matrix, model->basic_index, model->UL);

  inverse(model->matrix, model->basic_index, model->UL, model->inv);
  left_vet_mut(model->inv, model->b, model->x, NULL, model->basic_index);

  model->actual_cost = inner_product(model->cost, model->x);
}
