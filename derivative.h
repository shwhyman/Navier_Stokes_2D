#ifndef DERIVATIVE_H_INCLUDED
#define DERIVATIVE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

extern int N_i, N_j;

#define INDEX(i,j) ((N_i+2)*(i)+(j))

extern int * proc;			                    //Process indexed by vertex
extern int * i_min, * i_max;		            //Min, Max vertex indices of processes
extern int * left_proc, * right_proc;	        //Processes to left and right

extern int my_rank;

double partial_U_Squared_wrt_x(double * U, double dx, int i, int j);

double partial_V_Squared_wrt_y(double * V, double dy, int i, int j);

double partial_UV_wrt_x(double * U, double * V, double dx, int i, int j);

double partial_UV_wrt_y(double * U, double * V, double dy, int i, int j);

double partial_squared_U_wrt_x(double * U, double dx, int i, int j);

double partial_squared_V_wrt_x(double * V, double dx, int i, int j);

double partial_squared_U_wrt_y(double * U, double dy, int i, int j);

double partial_squared_V_wrt_y(double * V, double dy, int i, int j);

#endif // DERIVATIVE_H_INCLUDED
