#include  "../headfiles/matrix.h"

void print_indexed_matrix(MAT *mat, int *index, int size){
  int i, j;

  for(i = 0; i < mat->row_qtd; i++){
    for(j = 0; j < size; j++){
      printf(" %.2lf ", access_matrix(i, index[j], mat));
    }
    printf("\n");
  }
}

void print_matrix(MAT *mat){
  int i, j;

  for(i = 0; i < mat->row_qtd; i++){
    for(j = 0; j < mat->col_qtd; j++){
      printf(" %.2lf ", access_matrix(i, j, mat));
    }
    printf("\n");
  }
}

void print_vector(VECTOR *v){
  int i;

  for(i = 0; i < v->size; i++){
    printf(" %.4lf ", access_vector(v, i));
  }
  printf("\n");
}

void swap_i(int *x, int *y){
  int tmp = *x;
  *x = *y;
  *y = tmp;
}

void swap_d(double *x, double *y){
  double tmp = *x;
  *x = *y;
  *y = tmp;
}

void swap_cols(MAT* mat, int col_1, int col_2){
  swap_i(&mat->col_index[col_1], &mat->col_index[col_2]);
}

void swap_rows(MAT* mat, int row_1, int row_2){
  swap_i(&mat->row_index[row_1], &mat->row_index[row_2]);
}

void swap_vector_elements(VECTOR *v, int row_1, int row_2){
  swap_i(&v->index[row_1], &v->index[row_2]);
}

void destroy_vector(VECTOR *v){
  if(v != NULL){
    if(v->vector != NULL) free(v->vector);
    free(v);
  }
}

void destroy_matrix(MAT *mat){
  int i;

  if(mat != NULL){
    if(mat->matrix != NULL){
      for(i = 0; i < mat->row_qtd; i++){
        if(mat->matrix[i] != NULL) free(mat->matrix[i]);
      }
      free(mat->matrix);
    }

    free(mat);
  }
}

VECTOR *make_vector(int size){
  VECTOR *v = NULL;
  v = malloc(sizeof(VECTOR));
  if(v){
    v->vector = (double *)malloc(sizeof(double) * size);
    v->size = size;
    v->index = NULL;
  }

  return v;
}

MAT *make_matrix(int row_qtd, int col_qtd){
  int i;

  MAT *mat = NULL;
  mat = (MAT *)malloc(sizeof(MAT));
  mat->col_qtd = col_qtd;
  mat->row_qtd = row_qtd;

  mat->matrix = calloc(mat->row_qtd, sizeof(double*));

  for(i = 0; i < mat->row_qtd; i++) mat->matrix[i] = malloc(mat->col_qtd * sizeof(double));

  mat->row_index = NULL;
  mat->col_index = NULL;

  return mat;
}

void resize_vector(VECTOR *v, int additional){
  v->size += additional;
  v->vector = realloc(v->vector, sizeof(double) * v->size);
}

void resize_matrix(MAT *mat, int addition_row, int addition_col){
  int i;

  if(addition_row > 0){
    mat->row_qtd += addition_row;
    mat->matrix = realloc(mat, sizeof(double *) * mat->row_qtd);
  }

  if(addition_col > 0){
    mat->col_qtd += addition_col;
    for(i = 0; i < mat->row_qtd; i++){
      mat->matrix[i] = realloc(mat->matrix[i], sizeof(double) * mat->col_qtd);
    }
  }
}

// Multiplies the vector x at the left of matrix A
// The result goes on r
void right_vet_mut(MAT *A, VECTOR *x, VECTOR *r){
  int i, j;

  for(i = 0; i < r->size; i++){
    access_vector(r, i) = 0;
    for(j = 0; j < r->size; j++){
      access_vector(r, i) += access_vector(x, j) * access_matrix(i, j, A);
    }
  }
}

// Yields, on u, the left multiplication of the ith column vector of matrix mat by matrix inv
void mat_right_vet_mut(MAT *inv, MAT* mat, VECTOR *r, int i){
  int j, k;

  for(j = 0; j < mat->row_qtd; j++){
    access_vector(r, j) = 0;
    for(k = 0; k < mat->row_qtd; k++){
      access_vector(r, j) += access_matrix(k, i, mat) * access_matrix(j, k, inv);
    }
  }
}




// Left multiplies vector v, indexed by the index v_index, by matrix A (does v * A)
// Puts result on r, indexed by v2_index
// If some of the indices vectors are NULL, start from 0
void left_vet_mut(MAT *A, VECTOR *v, VECTOR *r, int *v_index, int *vr_index){
  int i, j;
  for(i = 0; i < r->size; i++) r->vector[i] = 0;

  if(v_index && vr_index){
    for(i = 0; i < A->row_qtd; i++){
      for(j = 0; j < A->col_qtd; j++){
        access_vector(r, vr_index[j]) += access_matrix(i, j, A) * access_vector(v, v_index[i]);
      }
    }
  } else if(v_index){
    for(i = 0; i < A->row_qtd; i++){
      for(j = 0; j < A->col_qtd; j++){
        access_vector(r, j) += access_matrix(i, j, A) * access_vector(v, v_index[i]);
      }
    }
  } else if(vr_index){
    for(i = 0; i < A->row_qtd; i++){
      for(j = 0; j < A->col_qtd; j++){
        access_vector(r, vr_index[j]) += access_matrix(i, j, A) * access_vector(v, i);
      }
    }
  } else {
    for(i = 0; i < A->row_qtd; i++){
      for(j = 0; j < A->col_qtd; j++){
        access_vector(r, j) += access_matrix(i, j, A) * access_vector(v, i);
      }
    }
  }
}

/*void left_vet_mut(MAT *A, VECTOR *v, VECTOR *r, int *v_index, int *vr_index){*/
  /*int i, j;*/
  /*for(i = 0; i < A->row_qtd; i++){*/
    /*r->vector[i] = 0;*/
    /*for(j = 0; j < A->col_qtd; j++){*/
      /*if(v_index && vr_index) access_vector(r, vr_index[j]) += access_matrix(i, j, A) * access_vector(v, v_index[i]);*/
      /*else if(v_index)        access_vector(r, j) += access_matrix(i, j, A) *           access_vector(v, v_index[i]);*/
      /*else if(vr_index)       access_vector(r, vr_index[j]) += access_matrix(i, j, A) * access_vector(v, i);*/
      /*else                    access_vector(r, j) += access_matrix(i, j, A) *           access_vector(v, i);*/
    /*}*/
  /*}*/
/*}*/

// Left multiplies the ith column vector of matrix C, indexed by the index c_index, by matrix A
// Puts result on r, indexed by v2_index
// If some of the indices vectors are NULL, start from 0
void mat_left_vet_mut(MAT *A, MAT *C, VECTOR *r, int k, int *c_index, int *vr_index){
  int i, j;

  for(i = 0; i < A->row_qtd; i++){
    r->vector[i] = 0;
    for(j = 0; j < A->col_qtd; j++){
      if(c_index && vr_index) access_vector(r, vr_index[j]) += access_matrix(i, j, A) * access_matrix(c_index[i], k, C);
      else if(c_index)        access_vector(r, j) += access_matrix(i, j, A) * access_matrix(c_index[i], k, C);
      else if(vr_index)       access_vector(r, vr_index[j]) += access_matrix(i, j, A) * access_matrix(i, k, C);
      else                    access_vector(r, j) += access_matrix(i, j, A) * access_matrix(i, k, C);
    }
  }
}


// Makes u = lth row of B * A
void lth_mult_row(MAT *B, MAT *A, VECTOR *u, int l){
  int j, k, i;

  for(i = 0; i < B->col_qtd; i++){
    access_vector(u, i) = 0;
    for(j = 0; j < A->row_qtd; j++){
      access_vector(u, i) += access_matrix(j, i, A) * access_matrix(l, j, B);
    }
  }
}



void mat_mut(MAT *A, MAT *B, MAT *C, int *a, int *b, int *c, int size){
  int i, j, k;

  for(i = 0; i < size; i++){
    for(j = 0; j < size; j++){
      for(k = 0; k < size; k++){
        access_matrix(i,c[j], C) += access_matrix(k, b[j], B) * access_matrix(i, a[k], A);
        //C->matrix[i][j] += B->matrix[k][j] * A[i][k];
      }
    }
  }
}
double ith_line_lth_column_mut(MAT *B, MAT *A, int i, int l){
  int m;
  double r = 0;

  for(m = 0; m < A->col_qtd; m++){
    r += access_matrix(l, m, B) * access_matrix(m, i, A);
  }

  return r;
}

double inner_product(VECTOR *v, VECTOR *u){
  int i;
  double r = 0.0;

  for(i = 0; i < v->size; i++){
    r += access_vector(v, i) * access_vector(u, i);
  }

  return r;
}

// Returns the inner product of ith column of matrix m with vector v
double mat_inner_product(MAT* mat, VECTOR *v, int i){
  int j;
  double r;

  for(j = 0, r = 0; j < mat->row_qtd; j++) r += access_matrix(j, i, mat) * access_vector(v, j);

  return r;
}

void initialize_vector(VECTOR *a, double (*value)(int)){
  int i;

  for(i = 0; i < a->size; i++){
    access_vector(a, i) = value(i);
  }
}

// Multiplies row_1 by multiplier and sums it to row_2
void row_operation(MAT *mat, int row_1, int row_2, double multiplier){
  int i;

  for(i = 0; i < mat->col_qtd; i++){
    access_matrix(row_2, i, mat) += access_matrix(row_1, i, mat) * multiplier;
  }
}

double nil_value(int i){ return 0; }

// Return the cth column of matrix as a vector
VECTOR *get_column_vector(MAT *matrix,  int c){
  int i, j;
  VECTOR *col = NULL;
  col = make_vector(matrix->row_qtd);

  for(i = 0, j = c; i < matrix->row_qtd; i++){
    col->index = matrix->row_index;
    access_vector(col, i) = access_matrix(i, j, matrix);
  }
  return col;
}

void inverse(MAT *mat, int *basic_index, MAT *UL, MAT *inv){
  int i, j;

  VECTOR *ei = NULL;
  VECTOR *col = NULL;
  VECTOR *tmp = make_vector(UL->col_qtd);

  inv->row_index = mat->row_index;
  tmp->index = UL->col_index;

  /*Make ei, the ith component of I, the identity matrix*/
  ei = make_vector(mat->row_qtd);
  ei->index = mat->col_index;

  for(j = 0; j < ei->size; j++){
    access_vector(ei, j) = 0;
  }

  for(i = 0; i < mat->row_qtd; i++){
    col = get_column_vector(mat, basic_index[i]);

    access_vector(ei, i) = 1;
    if(i > 0) access_vector(ei, i-1) = 0;

    // Solves for BASIC_MATRIX * col = ei
    lup_solve(UL, col, tmp, NULL, ei);
    compose_matrix(inv, col, i);
  }

  destroy_vector(ei);
  destroy_vector(tmp);
}

// Turn the col vector into the ith column vector of matrix
// by copying it
void compose_matrix(MAT *matrix, VECTOR *col, int i){
  int j;

  for(j = 0; j < matrix->row_qtd; j++){
    access_matrix(j, i, matrix) = access_vector(col, j);
  }

  destroy_vector(col);
}

// Receives two vectors as arguments
// Copies vector a into vector b
// If b is non NULL, destroy it first
void copy_vector(VECTOR *a, VECTOR **b){
  int i;

  if(*b) destroy_vector(*b);

  (*b) = malloc(sizeof(VECTOR));
  (*b)->size = a->size;
  (*b)->vector = malloc(sizeof(double) * a->size);

  for(i = 0; i < a->size; i++){
    (*b)->vector[i] = a->vector[i];
  }
}

// Multiplies the Lth sub line of mat by sub vet, determined by the offset
double mult_l_c(MAT *mat, VECTOR *vet, int l,  int offset){
  int i;
  double r = 0.0;

  for(i = 0; i < offset; i++){
    r += access_matrix(l, i, mat) * access_vector(vet, i);
  }
  return r;
}


// Solves for Ux = y
void backward_substitution(MAT *UL, VECTOR *x, VECTOR *y, int *basic_index){
  int j, i;

  if(basic_index != NULL){
    for(i = UL->col_qtd - 1; i >= 0; i--){
      access_vector(x, basic_index[i]) = access_vector(y, basic_index[i]);

      for(j = UL->col_qtd - 1; j > i; j--){
        access_vector(x, basic_index[i]) -= access_matrix(i, j, UL)*access_vector(x, j);
      }

      access_vector(x, basic_index[i]) /= access_matrix(i, i, UL);
    }
  } else {
    for(i = UL->col_qtd - 1; i >= 0; i--){
      access_vector(x, i) = access_vector(y, i);

      for(j = UL->col_qtd - 1; j > i; j--){
        access_vector(x, i) -= access_matrix(i, j, UL)*access_vector(x, j);
      }

      access_vector(x, i) /= access_matrix(i, i, UL);
    }

  }
}

// Solves for Ly = b
// by doing:
//  yn = bn - inner_product of Ln(except for the nth element, that is, ln*yn) and y(n-1)
void forward_elimination(MAT *UL, VECTOR *y, VECTOR *b, int *basic_index){
  int i;
  if(basic_index != NULL){
    for(i = 0; i < UL->row_qtd ; i++){
      access_vector(y, basic_index[i]) = access_vector(b, i) - mult_l_c(UL, y, i, i-1);
    }
  } else {
    for(i = 0; i < UL->row_qtd ; i++){
      access_vector(y, i) = access_vector(b, i) - mult_l_c(UL, y, i, i-1);
    }
  }
}


// Receives as arguments one matrix mat, one matrix UL, and a int vector a
// Takes the UL factorization of mat and puts it into UL
// The vector 'a' is the vector for basis indices for 'mat's columns
void lup_factor(MAT *mat, int *a, MAT *UL){
  int pivot_index, pivot, multiplier, row, col;
  int i, j, size = UL->col_qtd;

  UL->row_index = mat->row_index;
  UL->col_index = mat->col_index;

  // Copy mat into UL
  for(i = 0; i < UL->row_qtd; i++)
    for(j = 0; j < UL->col_qtd; j++)
      access_matrix(i, j, UL) = access_matrix(i, a[j], mat);

  for(pivot_index = 0; pivot_index < size; pivot_index++){

    /*Find non-zero element for pivot position and switch rows if necessary*/
    if(access_matrix(pivot_index, pivot_index, UL) == 0){
      for(row = pivot_index + 1; row < size; row++){
        if(access_matrix(row, pivot_index, UL) != 0){

          //TEST
          swap_i(&UL->row_index[pivot_index], &UL->row_index[row]);
        }
      }
    }

    /* Make UL*/
    pivot = access_matrix(pivot_index, pivot_index, UL);

    /*Make U*/
    for(row = pivot_index + 1; row < size; row++){
      multiplier = -access_matrix(row, pivot_index, UL) / pivot;

      for(col = row + 1; col < size ; col++){
        access_matrix(row, col, UL) += multiplier * access_matrix(pivot_index, col, UL);
      }
    }

    /*Make L*/
    for(col = pivot_index + 1; col < size; col++){
      multiplier = -access_matrix(pivot_index, col, UL) / pivot;

      for(row = col + 1; row < size ; row++){
        access_matrix(row, col, UL) += multiplier * access_matrix(row, pivot_index, UL);
      }
    }
  }
}

// First solve for Ly = b
// then, for Ux = y
void lup_solve(MAT* UL, VECTOR *x, VECTOR *y, int *basic_index, VECTOR *b){
  forward_elimination(UL, y, b, basic_index);
  backward_substitution(UL, x, y, basic_index);
}

void drop_row(MAT *mat, int row){
  int j;
  for(j = row; j < mat->row_qtd + 1; j++) mat->row_index[j] = mat->row_index[j+1];
  free(mat->matrix[row]);
  mat->row_qtd--;
}

// Inverts the comparison operation on the ith position of vector v
// That is:
//  < becomes >
//  >  //     <
void invert_comp(char *v, int i){
  if(v){
    if(v[i] == '<') v[i] = '>';
    else if(v[i] == '>') v[i] = '<';
  }
}
