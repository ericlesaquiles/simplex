#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "../headfiles/simplex.h"
#include "../headfiles/dual.h"
#include "../headfiles/interface.h"




/* Receives a file name and a char vector for symbols
 * as arguments and creates a simplex model for it by the following conventions:
 *
 *  1st line: int - quantity of coefficients (of variables) for model->cost function
 *  2nd line: int - quantity of equal signs
 *  3rd line: one int a, for the row size
 *  4th line: int(s) - indexes of the unbounded variables .. if the ith variable is unbounded, its index must appear
 *  5th line: coefficients for model->cost function
 *  6th line onward: b doubles followed by ">", "<" or "=" symbols followed by one double, for the value of b
 *
 *  obs.: It is assumed that the cost function will use all the variables, so it is used to define the
 *        quantity of columns on the constraint matrix
 *        (that is, if some variable does not contribute to the cost, its coefficient must be zero)
 *
 * */

/*
 * Remember:
 *    The columns index of the matrix is the same as the index of the cost vector and the onde of the x vector
 *    The rows index of the matrix is the same as the index of the b vector
 *
 * */

SIMPLEX *make_standard_form(const char *filename, char **symbols){
  int coef_size, qtd_row, qtd_col, i, j, eq;
  char temp;
  double tmp;

  SIMPLEX *model = NULL;

  // Open the file
  FILE *F = NULL;
  F = fopen(filename, "r");

  if(F == NULL){
    /*printf("\n Error reading from file.\n");*/
    exit(1);
  }

  model = malloc(sizeof(SIMPLEX));

  // Make first version of model->cost vector
  // by first taking the initial quantity of variables
  fscanf(F, "%d", &coef_size);

  fscanf(F, "%c", &temp);
  fscanf(F, "%d", &eq);

  fscanf(F, "%c", &temp);
  fscanf(F, "%d", &qtd_row);

  coef_size += qtd_row - eq;
  model->new_cols = qtd_row - eq;

  // then by taking the indexes of the positive variables
  model->splitted_var = calloc(coef_size - model->new_cols, sizeof(int));

  fscanf(F, "%c", &temp); // eats the '\n' of the first line.. this will be recurrent
  while(fscanf(F, "%c", &temp) && temp != '\n'){
    if(!isblank(temp)){
      model->splitted_var[atoi(&temp)] = 1;
      coef_size++;
    }
  }

  model->cost = make_vector(coef_size);
  model->cost->index = malloc(sizeof(int) * coef_size);

  // Initializes the cost index and reads the vector
  for(i = j = 0; i < coef_size - model->new_cols; i++){
    model->cost->index[i] = i;
    fscanf(F, "%lf", &tmp);
    access_vector(model->cost, i) = tmp;

    if(model->splitted_var[i-j]){
      i++;
      model->cost->index[i] = i;
      access_vector(model->cost, i) = -tmp;
      j++;
    }
  }
  for(i = coef_size - model->new_cols; i < coef_size; i++){
    model->cost->index[i] = i;
    access_vector(model->cost, i) = 0;
  }

  // Make first version of Matrix: model->matrix and of model->b
  qtd_col = coef_size;
  model->matrix = make_matrix(qtd_row, qtd_col);
  *symbols = malloc(sizeof(char) * qtd_row);

  model->b = make_vector(qtd_row);

  // Initializes index for rows of model matrix
  model->matrix->row_index = malloc(sizeof(int) * qtd_row);
  for(i = 0; i < qtd_row; i++) model->matrix->row_index[i] = i;
  model->matrix->col_index = model->cost->index;
  model->b->index = model->matrix->row_index;

  // Reads the matrix values from file
  for(i = temp = 0; i < qtd_row; i++){
    for(j = 0; j < qtd_col - model->new_cols; j++){
      fscanf(F, "%lf", &tmp);
      access_matrix(i, j, model->matrix) = tmp;
      if(model->splitted_var[j-temp]){
        access_matrix(i, ++j, model->matrix) = -tmp;
        tmp++;
      }
    }
    fgetc(F);
    fscanf(F, "%c", &(*symbols)[i]);
    fgetc(F);
    fscanf(F, "%lf", &access_vector(model->b, i));
  }

  // Initializes the slack part of the matrix
  for(i = coef_size - model->new_cols, temp = 1; i < coef_size; i++, temp = 1){
    temp = 1;
    for(j = 0; j < qtd_row; j++){
      if(temp){
        if((*symbols)[j] != '='){
          if((*symbols)[j] == '<'){
            access_matrix(j, i, model->matrix) = 1;
          }
          else if((*symbols)[j] == '>'){
            access_matrix(j, i, model->matrix) = -1;
          }
          temp = 0;
          (*symbols)[j] = '=';
        }
        else access_matrix(j, i, model->matrix) = 0;
      }
      else access_matrix(j, i, model->matrix) = 0;
    }
  }

  // Finish initialization of model
  model->x = make_vector(model->matrix->col_qtd);
  model->x->index = model->cost->index;
  model->y = make_vector(model->x->size);
  model->y->index = model->cost->index;

  fclose(F);

  return model;
}

typedef struct dic{
  char new_iteration[110];
  char lacking_file[50];
  char infinite_opt[50];
  char optimal_cost[50];
  char var_val     [50];
} DIC;

void init_dict(DIC *dict){
  strcpy(dict->new_iteration, " Voce gostaria de resolver o mesmo problema para um 'b' diferente? (s/n) ");
  strcpy(dict->lacking_file, " Please, send me the name of a file as a parameter.\n");
  strcpy(dict->infinite_opt,  " \n\n Infinite optimizatino.\n\n");
  strcpy(dict->optimal_cost,  " \n Custo otimo eh: %.4lf\n");
  strcpy(dict->var_val     ,  "\n Os valores das variaveis sao");
}

double read_d(){
  double r;
  scanf("%lf", &r);
  getc(stdin);
  return r;
}

void enter_vector(VECTOR *v){
   int i;

   printf(" Por favor, entre com o vetor a seguir (digite um numero seguido de um enter .. lembre-se de que o vetor b tem %d posicoes):\n", v->size);

   for(i = 0; i < v->size; i++){
    printf("    ");
    access_vector(v, i) = read_d();
   }
}


int main(int argc, char *argv[]){
  int r, flag, is_dual = 0;
  char *symbols, yn = 'n';
  SIMPLEX *model = NULL;
  DIC dict;
  init_dict(&dict);

  if(argc == 1){
    printf(" Please, send the name of a file as a parameter.\n");
  } else {
    model = make_standard_form(argv[1], &symbols);
    flag = phase_I(model, symbols);
    if(flag) set_up_model_env(model);

    do{
      r = is_dual ? dual(model) : simplex(model);

      if(r == -1){
        printf("\n\n Infinite optimization.\n\n");
        break;
      } else if(r == 1) {
        printf("\n Custo otimo eh: %.4lf\n", model->actual_cost);
        printf("\n Os valores das variaveis sao: \n");

        for(int i = 0; i < model->matrix->col_qtd - model->new_cols; i++) printf(" X[%d] = %lf ", i, model->x->vector[i]);
        printf("\n");

        // Checks whether user wants to run program with a different 'b' vector
        printf("%s", dict.new_iteration);
        scanf("%c", &yn);
        getc(stdin);
        if(yn == 's') enter_vector(model->b);
        is_dual = !is_dual;
        left_vet_mut(model->inv, model->b, model->x, NULL, model->basic_index);
      }
    } while(yn == 's');
    destroy_model(model);
  }

  return 0;
}
