#ifndef ALGORITHM_H_INCLUDED
#define ALGORITHM_H_INCLUDED

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

int initialise_theta_function_field(double * U, double * V, double dx, double dy, int x_nodes, int y_nodes, double A, double a, double b);

int get_gaussian_nodes(double * U, double * V, double x_pos, int i_array, double y_pos, int j_array, double i_unit, double j_unit, int x_nodes, int y_nodes, double A, double a, double b);

void init_pressure(double * P);

double * allocate_NULL();

int set_F_and_G(double * U, double * V, double * F, double * G, double dx, double dy, double Re, double dt, double g_x, double g_y);

int set_RHS(double * R, double * F, double * G, double dx, double dy, double dt);

int updateField(double * U, double * V, double * F, double * G, double * P_new, double dx, double dy, double dt);

#endif // ALGORITHM_H_INCLUDED
