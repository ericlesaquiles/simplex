#include <stdio.h>
#include "../headfiles/matrix.h"

int TESTE = 1;

void printv(const char* msg, VECTOR *v){
  if(TESTE){
    printf(" %s", msg);
    print_vector(v);
    printf("\n");
  }
}
void printm(const char* msg, MAT *mat){
  if(TESTE){
    printf(" %s", msg);
    print_matrix(mat);
    printf("\n");
  }
}

void print_all(SIMPLEX *model){
  if(TESTE){
    int i, j;
    printf(" Constraints matrix:\n");
    print_matrix(model->matrix);
    printf("\n");
    printf(" Constraints vector:\n");
    print_vector(model->b);
    printf("\n");

    printf(" Basic cost vector:\n");
    for(i = 0; i < model->basis_size; i++) printf(" %.2lf ", access_vector(model->cost, model->basic_index[i]));
    printf("\n\n");

    printf(" UL matrix:\n");
    print_matrix(model->UL);
    printf("\n");

    printf(" Basic matrix:\n");
    for(i = 0; i < model->matrix->row_qtd; i++){
     for(j = 0; j < model->basis_size; j++){
      printf(" %.2lf ", access_matrix(i, model->basic_index[j], model->matrix));
     }
     printf("\n");
    }
    printf("\n");

    printf(" Inverse basic matrix:\n");
    print_matrix(model->inv);
    printf("\n");

    printf(" Basis index:\n");
    for(i = 0; i < model->basis_size; i++){
      printf(" %d ", model->basic_index[i]);
    }
    printf("\n");
    printf("\n");

    printf(" Non basic indices:\n");
    for(i = 0; i < model->matrix->col_qtd - model->basis_size; i++){
      printf(" %d ", model->non_basic_index[i]);
    }
    printf("\n");
    printf("\n");

    printf(" x vector:\n");
    print_vector(model->x);
    printf("\n");
    printf(" b vector:\n");
    print_vector(model->b);
    printf("\n");
  }
}

void print_initial(SIMPLEX *model){
  if(TESTE){
    printf(" Constraints matrix:\n");
    print_matrix(model->matrix);
    printf("\n");
    printf(" Cost vector:\n");
    print_vector(model->cost);
    printf("\n");
    /*printf(" x vector:\n");*/
    /*print_vector(model->x);*/
    printf("\n");
    printf(" b vector:\n");
    print_vector(model->b);
    printf("\n");
  }
}

