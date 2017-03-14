#ifndef __matrix_h
#define __matrix_h

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct vector{
  double *vector;
  int size;
  int *index;
} VECTOR;

typedef struct matrix{
  double **matrix;

  int row_qtd;
  int col_qtd;

  int* row_index;
  int* col_index;
} MAT;


#define access_matrix(i, j, name) name->matrix[name->row_index[i]][name->col_index[j]]
#define access_vector(v,i) v->vector[v->index[i]]

//#define _swap(x,y) {typeof(x) tmp = x; x = y; y = tmp;}

// Returns an instance of a matrix structure with quantity of rows and columns as indicated in the variables
MAT *make_matrix(int row_qtd, int col_qtd);
VECTOR *make_vector(int size);

void initialize_vector(VECTOR *a, double (*value)(int));

void destroy_matrix(MAT *matrix);
void destroy_vector(VECTOR *v);

double ith_line_lth_column_mut(MAT *B, MAT *A, int i, int l);
void mat_mut(MAT *A, MAT *B, MAT *C, int *a, int *b, int *c, int size);
void mat_right_vet_mut(MAT * inv, MAT* mat, VECTOR *r, int i);
void lth_mult_row(MAT *B, MAT *A, VECTOR *u, int l);
void right_vet_mut(MAT *A, VECTOR *x, VECTOR *r);
void left_vet_mut(MAT *A, VECTOR *x, VECTOR *r, int *v1_index, int *v2_index);
void mat_left_vet_mut(MAT *A, MAT *C, VECTOR *r, int k, int *c_index, int *vr_index);
void inverse(MAT *mat, int *basic_index, MAT *UL, MAT *inv);

double mult_l_c(MAT *mat,  VECTOR *vet, int l, int offset);

void print_indexed_matrix(MAT *mat, int *index, int size);
void print_matrix(MAT *mat);

void print_vector(VECTOR *v);
void copy_vector(VECTOR *a, VECTOR **b);

void resize_vector(VECTOR *v, int additional);
void resize_matrix(MAT *mat, int addition_row, int addition_col);
void drop_row(MAT *mat, int row);

double inner_product(VECTOR *v, VECTOR *u);
double mat_inner_product(MAT* mat, VECTOR *v, int i);
void row_operation(MAT *mat, int row_1, int row_2, double multiplier);

VECTOR *get_column_vector(MAT *matrix, int c);
void compose_matrix(MAT *matrix, VECTOR *col, int i);
void inverse_comp(char *v, int i);

void swap_cols(MAT* mat, int col_1, int col_2);
void swap_rows(MAT* mat, int row_1, int row_2);
void swap_i(int *x, int *y);

void lup_factor(MAT *mat, int *a, MAT *UL);
void forward_elimination(MAT *UL, VECTOR *y, VECTOR *b, int *basic_index);
void backward_substitution(MAT *UL, VECTOR *x, VECTOR *y, int *basic_index);
void lup_solve(MAT *UL, VECTOR *x, VECTOR *y, int *basic_index, VECTOR *b);

double nil_value(int i);

#endif
