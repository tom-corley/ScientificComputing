#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mkl_lapacke.h>

/* TO RUN USE THE FOLLOWING COMMANDS:
  1. module load intel impi imkl
  2. gcc -Wall -Werror -std=c99 -o a5 -lmkl -liomp5 -lm  assign5.c
  3. ./a5
*/

/* Band matrix struct (copied from band_utility.c) */
struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
typedef struct band_mat band_mat;

/* Utility functions for implicit solver and printing matrices (mostly copied from band_utility.c)
  1. init_band_mat
  2. finalise_band_mat
  3. getp
  4. getv
  5. setv
  6. printmat 
  7. print_flattened_matrix
  8. print_banded_matrix
  9. solve_Ax_eq_b
*/ 

// 1. Sets up banded matrix for implicit solver
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }  
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};

// 2. Deallocates memory used for banded matrix
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

// 3. Returns a pointer to an entry in banded matrix, given row and column index
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

// 4. Returns entry in banded matrix, given row and column index
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}       

// 5. Sets a entry in banded matrix
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

// 6. Prints the matrix in banded format
int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}

// 7. Function to print a matrix that has been flattened into an array
void print_flattened_matrix(double* matrix, int N_x, int N_y) {
    for (int i = 0; i < N_x; i++) {
        for (int j = 0; j < N_y; j++) {
            printf("%8.4f ", matrix[i * N_y + j]); 
        }
        printf("\n"); 
    }
}

// 8. Function to print a banded matrix as a normal matrix
void print_banded_matrix(band_mat *bmat) {
    for (long i = 0; i < bmat->ncol; i++) {
        for (long j = 0; j < bmat->ncol; j++) {
            if (abs(i - j) > bmat->nbands_up && abs(i - j) > bmat->nbands_low) {
              // Index is outside bands so print 0
              printf("0.0 "); 
            } else {
              // Retrieve entry at index and print it
              printf("%g ", getv(bmat, i, j)); 
            }
        }
        // New line per row
        printf("\n");
    }
}

// 9. Solves the equation Ax = b for a matrix stored in band format and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) { 
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

/* Functions for memory management and file handling
  1. MemSwap
  2. read_inputs
  3. read_coefficients
  4. write_final_temp
*/

// 1. Function to swap pointers to arrays, effectively copying them to eachothers addresses
void MemSwap(double** a, double** b) {
  double * temp = *a;
  *a = *b;
  *b = temp;
}

// 2. Function to read input parameters from file
void read_inputs(const char *filename, int *N_x, int *N_y, double *t_f, int *I_min) {
    // Open File Connection
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file containing input parameters");
        exit(EXIT_FAILURE);
    }

    // Read in data
    fscanf(file, "%d %d %lf %d", N_x, N_y, t_f, I_min);

    // Close file connection
    fclose(file);
}

// 3. Function to read coefficients from file
void read_coefficients(const char *filename, double *T, double *Q11, double *Q22, double *R, double *S, int N_x, int N_y) {
    // Open File Connection
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening coefficients file");
        exit(EXIT_FAILURE);
    }

    // Read each element in row-major order into arrays
    for (int k = 0; k < (N_x*N_y); k++) {
      fscanf(file, "%lf %lf %lf %lf %lf", &T[k], &Q11[k], &Q22[k], &R[k], &S[k]);
    }

    // Close file connection
    fclose(file);
}

// 4. Function to write final temperature to output
void write_final_temp(const char *filename, int N_x, int N_y, double *T) {
    // Open File Connection
    FILE *output = fopen(filename, "w");
    if (!output) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }

    // Write final temperature to file
    for (int i = 0; i < N_x; i++) {
        for (int j = 0; j < N_y; j++) {
            fprintf(output, "%g\n", T[i*N_y + j]);
        }
    }

    // Close file connection
    fclose(output);
}

/* Functions for finite difference method logic and time evolution of PDE
  1. is_boundary
  2. construct_banded_matrix
  3. implicit_method
*/

// 1. Function to determine if a given index (i,j) is at the boundary of the grid
bool is_boundary(int i, int j, int N_x, int N_y) {
  return (i == 0 || i == N_x-1 || j == 0 || j == N_y-1);
}

// 2. Function to construct banded matrix for problem
void construct_banded_matrix(band_mat *bmat, double dx, double dy, double dt, int N_x, int N_y, double* Q_11, double* Q_22, double* S, double* R) {
  // Setting each of the (N_x * N_y) rows
  for (int i = 0; i < N_x; i++) {
    for (int j = 0; j < N_y; j++) {
      if (is_boundary(i, j, N_x, N_y)) {
        // Set band for T_i,j to 1, leave rest as 0 for T_new = T_old
        setv(bmat, i*N_y + j, i*N_y + j, 1.0);
      } else {
        // Calculate coefficients for each of the 5 bands T_i+1,j; T_i-1,j; T_i,j+1; T_i,j-1; T_i,j
        double a1 = (-dt/(dx*dx)) * (0.25*(Q_11[(i+1)*N_y + j] - Q_11[(i-1)*N_y + j]) + Q_11[i*N_y + j]);
        double a2 = (-dt/(dx*dx)) * (0.25*(-Q_11[(i+1)*N_y + j] + Q_11[(i-1)*N_y + j]) + Q_11[i*N_y + j]);
        double a3 = (-dt/(dy*dy)) * (0.25*(Q_22[i*N_y + (j+1)] - Q_22[i*N_y + (j-1)]) + Q_22[i*N_y + j]);
        double a4 = (-dt/(dy*dy)) * (0.25*(-Q_22[i*N_y + (j+1)] + Q_22[i*N_y + (j-1)]) + Q_22[i*N_y + j]);
        double a5 = 1.0 + dt*(2*(Q_11[i*N_y + j])/(dx*dx) + 2*(Q_22[i*N_y + j])/(dy*dy) + R[i*N_y + j]);

        // Set these coefficients using interface
        setv(bmat, i*N_y + j, (i+1)*N_y + j, a1);
        setv(bmat, i*N_y + j, (i-1)*N_y + j, a2);
        setv(bmat, i*N_y + j, i*N_y + (j+1), a3);
        setv(bmat, i*N_y + j, i*N_y + (j-1), a4);
        setv(bmat, i*N_y + j, i*N_y + j, a5);
      }
    }
  }
}

// 3.. Function to set and solve matrix equation representing method at each time step
void implicit_method(double* T_old, double* T_new, int N_x, int N_y, double dt, band_mat *bmat, double* S) {
  // Set up rhs of matrix equation
  for (int i = 0; i < N_x; i++) {
    for (int j = 0; j < N_y; j++) {
      if (is_boundary(i, j, N_x, N_y)) {
        continue;
      } else {
        T_old[i*N_y + j] += dt * S[i*N_y + j];
      }
    }
  }

  // Solve Ax = b, for A = bmat, x = T_new, b = T_old
  solve_Ax_eq_b(bmat, T_new, T_old);
}


/* Main Method */
int main() {
    // Initialise input and output files
    const char *input_file = "input.txt";
    const char *coefficents_file = "coefficients.txt";
    const char *output_file = "output.txt";

    // Initialise and read input parameters, then set space and time steps
    double t_f; int N_x, N_y, I_min;
    read_inputs(input_file, &N_x, &N_y, &t_f, &I_min);
    double dx = 1.0 / (N_x - 1);
    double dy = 1.0 / (N_y - 1);
    double dt = t_f / I_min;

    // Allocate memory for grid variables, stored in 1D in row-major order, validating pointers
    double *T = malloc((N_x * N_y) * sizeof(double));
    if (T == NULL) {
        printf("Memory allocation failed for T\n");
        exit(EXIT_FAILURE);
    }
    double *T_new = malloc((N_x * N_y) * sizeof(double));
    if (T_new == NULL) {
        printf("Memory allocation failed for T_new\n");
        exit(EXIT_FAILURE);
    }
    double *Q_11 = malloc((N_x * N_y) * sizeof(double));
    if (Q_11 == NULL) {
        printf("Memory allocation failed for Q_11\n");
        exit(EXIT_FAILURE);
    }
    double *Q_22 = malloc((N_x * N_y) * sizeof(double));
    if (Q_22 == NULL) {
        printf("Memory allocation failed for Q_22\n");
        exit(EXIT_FAILURE);
    }
    double *R = malloc((N_x * N_y) * sizeof(double));
    if (R == NULL) {
        printf("Memory allocation failed for R\n");
        exit(EXIT_FAILURE);
    }
    double *S = malloc((N_x * N_y) * sizeof(double));
    if (S == NULL) {
        printf("Memory allocation failed for S\n");
        exit(EXIT_FAILURE);
    }

    // Read coefficients and intial values
    read_coefficients(coefficents_file, T, Q_11, Q_22, R, S, N_x, N_y);

    // Construct banded matrix
    band_mat bmat;
    long ncols = (N_x*N_y);
    long nbands_low = N_y;
    long nbands_up  = N_y; 
    init_band_mat(&bmat, nbands_low, nbands_up, ncols);
    construct_banded_matrix(&bmat, dx, dy, dt, N_x, N_y, Q_11, Q_22, S, R);

    // Evolving equations on grid
    double t = 0;
    while (t < t_f) {
      // Do implicit method for this time step
      implicit_method(T, T_new, N_x, N_y, dt, &bmat, S);

      // Save new temps as current temps
      MemSwap(&T, &T_new);

      // Increment time
      t += dt;
    } 

    // Writing final temperatures
    write_final_temp(output_file, N_x, N_y, T);

    // Freeing memory
    finalise_band_mat(&bmat);
    free(T);
    free(T_new);
    free(Q_11);
    free(Q_22);
    free(R);
    free(S);

    // Program exit
    return 0;
}