/*******************************************************
* This program solves coupled equations
* 
* To compile, run: gcc -Wall -Werror -std=c99 -lm
*
* List of identified errors:
*  Line           Brief description of a fix
* Number (in original)
* Line Offset from original = 14
* -----------------------------------------
* Errors:
*  N/A            ..... changed filename to assign3.c
*  16             ..... added .h to <math>
*  26             ..... changed index from i+1 to i
*  35             ..... removed -1 from stopping condition of for loop
*  40,49,50,51,52 ..... changed index from num_grid_points to num_grid_points - 1
*  49,50,55,56    ..... Incorrect coefficients in calculate_next, u should have positive mv, v negative mu
*  55, 56         ..... Added constant m into interaction term
*  37             ..... changed dx to x in initial condition
*  62             ..... fixed array swap by adding a temporary pointer
*  78             ..... Changed dx to length/(num_grid_points-1) as for 100 points you step 99 times 
*  79             ..... Changed dt to 0.004 to satisfy stability condition, rearrange (K*dt)/(dx^2) <= 1/2
                  ..... to dt <= (dx^2)/(2K), substituting in gives dt <= 0.00403
*  84             ..... initialised current_time as 0
*  93             ..... added ! (not) to if statement checking heap allocation
*  103            ..... changed > to < in for loop over timesteps
*  110            ..... getting rid of redundant memory allocations which lead to leaks
********************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI  3.141593

void PrintCurrent(double current_time,unsigned int num_grid_points, double dx, double* current_U, double* current_V){
  //function to print current state of U and V for a given time step
    for(unsigned int i=0; i<num_grid_points; i++){
      double x = dx * i;
      printf("%g %g %g %g\n", current_time, x, current_U[i], current_V[i]);
    }
}

void Initialise(unsigned int num_grid_points,double dx, double* current_U, double* current_V, double length){ // CHECKED
  //function to initialise U(x,t) at t = 0
  for(unsigned int i=0; i<num_grid_points; i++){
    double x = dx*i; // evaluates to 0, dx, 2dx, ... (n-1)dx
    current_U[i] = 1.0+sin((2*PI*x)/length); // Should u(0) 
  }
  current_U[0] = current_U[num_grid_points-1]; // sin(0) = sin(2pi) = sin(2pi*(length/length)) is periodic
  
}


void CalculateNext(double K, double m,double dt,unsigned int num_grid_points, double dx, double* current_U, double* current_V,  double* next_U, double* next_V){
  //function to calculate next time step of U and V

    // Dealing with points on boundary 
    next_U[0] = current_U[0] + (((K*dt)/(dx*dx)) * (current_U[1] + current_U[num_grid_points-1] - (2*current_U[0]))) + (dt*m*current_V[0]); // +K(...) + mv
    next_V[0] = current_V[0] + (((K*dt)/(dx*dx)) * (current_V[1] + current_V[num_grid_points-1] - (2*current_V[0]))) - (dt*m*current_U[0]); // +K(...) - mu
    next_U[num_grid_points-1] = next_U[0];     
    next_V[num_grid_points-1] = next_V[0];

  for(unsigned int i=1; i<num_grid_points-1; i++){ // Looping over points that are not on the boundary
    next_U[i] = current_U[i] + (((K*dt)/(dx*dx)) * (current_U[i+1] + current_U[i-1] - (2*current_U[i]))) + (dt*m*current_V[i]); // Central difference scheme for second derivative, dt*V is interaction term
    next_V[i] = current_V[i] + (((K*dt)/(dx*dx)) * (current_V[i+1] + current_V[i-1] - (2*current_V[i]))) - (dt*m*current_U[i]);
  } 
  
}

void MemSwap(double** a, double** b){ // CHECKED
  //function to swap two double arrays
  double * temp = *a;
  *a = *b;
  *b = temp;
  
}


int main(){ // CHECKED

  // declaring constants K, m, domain length, number of grid points and final simulation time (do not change)
  const double K = 3.6; // DO NOT CHANGE
  const double m = 5; // SAME
  const double length = 16.873;
  const unsigned int num_grid_points = 100;
  const double final_time = 1.0;

  //calculating time and length step size and number of time steps
  double dx = length/(num_grid_points-1); //  -1 since 0 is first then add dx 99 times to get to 100th
  double dt = 0.004; // There is a lot of grid scale instability as kdt/dx^2 = 12.64 > 1/2 was originally 0.1

  unsigned int num_time_steps = (final_time/dt)+1; // step at 0 plus the number of time steps to get to t=1

  //initialisation of current time
  double current_time = 0;

  //allocating current and next step U and V arrays
  double* current_U = malloc((num_grid_points)*sizeof(double));
  double* current_V = calloc(num_grid_points, sizeof(double));
  double* next_U = malloc((num_grid_points)*sizeof(double));
  double* next_V = malloc((num_grid_points)*sizeof(double));

  //check to determine if memory allocation has been performed correctly
  if (!(current_U && current_V && next_U && next_V) ) {
    printf("Memory allocation failed\n");
    return 0;
  }
  
  Initialise(num_grid_points, dx, current_U, current_V, length);
  
  PrintCurrent(current_time, num_grid_points, dx, current_U, current_V);

  //loop over timesteps
  for(unsigned int j=1; j<num_time_steps; j++){
    
    current_time += dt;
    
    CalculateNext(K, m, dt, num_grid_points, dx, current_U, current_V, next_U, next_V);

    //making next step current step
    MemSwap(&current_U, &next_U);
    MemSwap(&current_V, &next_V);
    
    PrintCurrent(current_time, num_grid_points, dx, current_U, current_V);
    
  }

  //memory clean up
  free(current_U);
  free(current_V);
  free(next_U);
  free(next_V);

  return 0;

}
