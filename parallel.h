#ifndef PARALLEL_H_INCLUDED
#define PARALLEL_H_INCLUDED

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

void make_domains(int num_procs);

void jacobi(int num_procs, double * R, double * P, double * P_new, double dx, double dy);

void SOR(int num_procs, double * R, double * P, double * P_new, double dx, double dy, double omega);

void set_F_G_and_R(int num_procs, double * F, double * G, double * R, double * U, double * V, double dx, double dy, double Re, double dt, double g_x, double g_y);

void timestamp();

#endif // PARALLEL_H_INCLUDED
