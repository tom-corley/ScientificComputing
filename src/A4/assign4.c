#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_lapacke.h>

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

/* Utility functions for implicit solver (mostly copied from band_utility.c)
  1. init_band_mat
  2. finalise_band_mat
  3. getp
  4. getv
  5. setv
  6. printmat 
  7. solve_Ax_eq_b
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
  printf("\nMatrix Bands:\n");
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  printf("\n END MATRIX \n");
  return 0;
}

// 7. Solves the equation Ax = b for a matrix stored in band format and x and b real arrays                                          */
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

/* Utility functions for file handling and memory swap: */

// Function to read input parameters from file
void read_input_params(const char *filename, double *L, int *N, double *tf, double *tD, double *Ulb, double *Urb, double *Vlb, double *Vrb) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file containing input parameters");
        exit(EXIT_FAILURE);
    }
    fscanf(file, "%lf %d %lf %lf %lf %lf %lf %lf", L, N, tf, tD, Ulb, Urb, Vlb, Vrb);
    fclose(file);
}

// Function to read initial conditions from file
void read_initial_conditions(const char *filename, double *U, double *V, int N) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file containing initial conditions");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
        fscanf(file, "%lf %lf", &U[i], &V[i]);
    }
    fclose(file);
}

// Function to write output to file in required format
void write_current_to_file(FILE *file, double current_time, unsigned int num_grid_points, double dx, double *current_U, double *current_V) {
    // Write the current state of U and V for the given time step in desired format
    for (unsigned int i = 0; i < num_grid_points; i++) {
        double x = dx * i;
        fprintf(file, "%g %g %g %g\n", current_time, x, current_U[i], current_V[i]);
    }
}

// Function to swap pointers to arrays, effectively copying them to eachothers addresses
void MemSwap(double** a, double** b){
  double * temp = *a;
  *a = *b;
  *b = temp;
}

/* Functions for time splitting method */

// Implicit time step using LAPACKE to solve matrix equation adjusted for B.Cs
void implicit_linear(double *U, double *V, int N, double dx, double dt, double *U_new, double *V_new, band_mat *bmat, double *b_vec, double Ulb, double Urb, double Vlb, double Vrb) {
    // Set up rhs for U
    for (int i = 0; i < N; i++) {
        b_vec[i] = U[i];
    }

    // Boundary case adjustments for U
    b_vec[0]   +=  -2 * (dt/dx) * Ulb;
    b_vec[N-1] +=   2 * (dt/dx) * Urb; 

    // Solve for U_new, saves solved values into U_new
    solve_Ax_eq_b(bmat, U_new, b_vec);

    // Set up rhs for V
    for (int i = 0; i < N; i++) {
        b_vec[i] = V[i]; 
    }

    // Boundary case adjustments for V
    b_vec[0]   += -2 * (dt/dx) * Vlb;
    b_vec[N-1] +=  2 * (dt/dx) * Vrb;

    // Solve for V_new, saves solved values into V_new
    solve_Ax_eq_b(bmat, V_new, b_vec);
}

// Explicit time step with non-linear terms 
void explicit_nonlinear(double *U, double *V, int N, double dt) {
    for (int i = 0; i < N; i++) {
        // Calculate squared norm Z^2
        double Z_sq = U[i] * U[i] + V[i] * V[i];

        // Update U and V explicitly
        U[i] += dt * (2 * U[i] - 4 * U[i] * Z_sq);
        V[i] += dt * (2 * V[i] - 4 * V[i] * Z_sq);
    }
}

// Update matrix for implicit solver given a new dt
void update_matrix(band_mat *bmat, double dx, long double dt, int N) {
  // Setting -a, 1+2a, -a pattern
  double alpha = dt / (dx*dx);
  for(long i=0; i < N; i++) {
    if(i>0)       {setv(bmat,i,i-1, -alpha);}; 
    setv(               bmat,i,i,   1 + 2*alpha); 
    if(i< N - 1) {setv(bmat,i,i+1,-alpha);}; }

  // Ghost point vn bcs, set top and bottom row super and sub diagonal respectivley to -2a
  setv(bmat, 0, 1, -2*alpha); 
  setv(bmat, N-1, N-2, -2*alpha); 
}

// Main method, sets up and then runs time splitting method
int main() {
    // Initialise input and output files
    const char *input_params_file = "input_params.txt";
    const char *init_cond_file = "init_cond.txt";
    const char *output_file = "output.txt";

    // Initialise and read input parameters
    double L, tf, tD, Ulb, Urb, Vlb, Vrb; int N;
    read_input_params(input_params_file, &L, &N, &tf, &tD, &Ulb, &Urb, &Vlb, &Vrb);

    // Allocate memory and read initial values
    double *U = (double *)calloc(N, sizeof(double));
    double *V = (double *)calloc(N, sizeof(double));
    double *U_new = (double *)calloc(N, sizeof(double));
    double *V_new = (double *)calloc(N, sizeof(double));
    read_initial_conditions(init_cond_file, U, V, N);

    // Determine space step dx and time step dt
    double dx = L / (N - 1);
    double max_dt = 0.5 * dx * dx; // CFL Stability constraint
    long double dt = tD;
    while (dt > max_dt) {dt *= 0.5;} // Half dt until stable
    int steps_per_diag = ceil(tD / dt); // Ensure dt is a divisor of tD
    dt = tD / steps_per_diag;
    if (dt > max_dt) { // Check dt is stable
        fprintf(stderr, "Error satisfying stability constraint\n");
        exit(EXIT_FAILURE); }

    // Setting up implicit solver with LAPACKE
    band_mat bmat;
    long ncols = N;
    long nbands_low = 1;
    long nbands_up  = 1; // Matrix is tri-diagonal
    init_band_mat(&bmat, nbands_low, nbands_up, ncols);
    double *b_vec = malloc(sizeof(double)*ncols); // Allocate rhs vector b
    update_matrix(&bmat, dx, dt, N);

    // Initialise file stream
    FILE *output = fopen(output_file, "w");
    if (!output) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }

    // Main time loop
    long double t = 0.0; long double old_dt = 0; double next_output_time = 0.0;
    double Z_sq_max = 0;
    int diag_count = 0;
    while (t < tf + 1e-8) {
        // Write to file at next diagnostic time step and reconfigure dt and update matrix if necessary
        if (fabs(t - next_output_time) <= 1e-7 || t >= next_output_time) {
            // Write to file and set next file write time
            write_current_to_file(output, t, N, dx, U, V);
            diag_count++;
            next_output_time = diag_count * tD; 

            // Reconfigure time step if necessary 
            Z_sq_max = 0;
            old_dt = dt;
            for (int i = 0; i < N; i++) {
              // Find maximum magnitude squared of peturbed solution
              Z_sq_max = fmax(Z_sq_max, U[i] * U[i] + V[i] * V[i]);
            }
            if (Z_sq_max > 0) {
              // Halve until stable relative to magnitude of Z^2
              if (dt > (1/(2*Z_sq_max))) {
                while ((Z_sq_max > 0) && dt > (1.0/(2*Z_sq_max))) {
                  dt /= 2;
                }
              }
              // Double while keeping stable relative to CFL and magnitude of dt
              else if ((2 * dt < (1.0/(2*Z_sq_max))) && (2*dt < max_dt)) {
                while ((2 * dt < (1.0/(2*Z_sq_max))) && (2*dt < max_dt)) {
                  dt = 2*dt;
                }
              }
            }
            if (fabs(dt - old_dt) > 1e-10) {
              update_matrix(&bmat, dx, dt, N);
            }
          }

        // Solve linear part (implicit) (Saves U^lin, V^lin to U_new and V_new)
        implicit_linear(U, V, N, dx, dt, U_new, V_new, &bmat, b_vec, Ulb, Urb, Vlb, Vrb);
        
        // Solve nonlinear part (explicit) (Increments U_new and V_new according to fin. diff.)
        explicit_nonlinear(U_new, V_new, N, dt);

        // Swap pointers (Stores U_new in U etc)
        MemSwap(&U, &U_new);
        MemSwap(&V, &V_new);

        // Advance dt to next step
        t += dt;
    }

    // Cleanup
    finalise_band_mat(&bmat);
    free(U);
    free(V);
    free(U_new);
    free(V_new);
    free(b_vec);
    fclose(output);

    // Program exit
    return 0; 
}