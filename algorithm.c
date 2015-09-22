#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "algorithm.h"
#include "derivative.h"

/***************************************************************************************************************************************/
int initialise_theta_function_field(double * U, double * V, double dx, double dy, int x_nodes, int y_nodes, double A, double a, double b)  //x, y nodes are the maximum node at which a gaussian will be placed
/*
  Purpose:
    Sum the contribution of a given number of 'Gaussian Fields' located at the centre of neighbouring unit cels
  Modified:
    August 2014
*/
/***************************************************************************************************************************************/
{
    int i, j;

    int i_cent = (int)(0.5*(N_i + 1));
    int j_cent = (int)(0.5*(N_j + 1));

    double i_unit = ((double)N_i)*dx;
    double j_unit = ((double)N_j)*dy;

    for(i = 1; i < N_i + 1; i++)
    {
        for(j = 1; j < N_j + 1; j++)
        {

            double x_pos = ((double)(i - i_cent))*dx;
            double y_pos = ((double)(j - j_cent))*dy;

            get_gaussian_nodes(U, V, x_pos, i, y_pos, j, i_unit, j_unit, x_nodes, y_nodes, A, a, b);

        }
    }

    //printf("i_cent = %d\nj_cent = %d\n\n", i_cent, j_cent);

    return 0;
}

/******************************************************************************************************************************************************************************************/
int get_gaussian_nodes(double * U, double * V, double x_pos, int i_array, double y_pos, int j_array, double i_unit, double j_unit, int x_nodes, int y_nodes, double A, double a, double b)
/*
  Purpose:
    Performs the summation described in the function above
  Modified
    August 2014
*/
/******************************************************************************************************************************************************************************************/
{
    int i, j;                                       //Now refers to gaussian node numbers

    for(i = -x_nodes; i < x_nodes + 1; i++)
    {
        for(j = -y_nodes; j < y_nodes + 1; j++)
        {
            double exp_to = (-1)*(((pow((x_pos - (((double)i)*i_unit)), 2.0))/(a)) + ((pow((y_pos - (((double)j)*j_unit)), 2.0))/(b)));
            double exp_term = (-2*A)*(exp(exp_to));

            U[INDEX(i_array, j_array)] += exp_term*((y_pos - (((double)j)*j_unit))/b);
            V[INDEX(i_array, j_array)] += (-1)*exp_term*((x_pos - (((double)i)*i_unit))/a);
        }
    }

    return 0;
}

/****************************************************************************************/
void init_pressure(double * P)
/*
  Purpose:
    Initialise the pressure fields
  Modified:
    August 2014
*/
/****************************************************************************************/
{
    int i;
    int ndof = ndof = (N_i + 2)*(N_j + 2);

    for(i = 0; i < ndof; i++)
    {
        P[i] = 0.0;
    }
    return;
}

/******************************************************************************/
double * allocate_NULL()
/*
  Purpose:
    ALLOCATE_ARRAYS creates and zeros out the arrays U and U_NEW.
  Modified:
    August 19 2014
*/
/******************************************************************************/
{
    int i;
    int ndof;

    double * array;

    ndof = (N_i + 2)*(N_j + 2);

    array = (double*)malloc(ndof*sizeof(double));

    for(i = 0; i < ndof; i++)
    {
        array[i] = 0.0;
    }

    return array;
}

/******************************************************************************/
void allocate_F_and_G(double * F, double * G)
/*
  Purpose:
    Allocate memory for F and G.
  Modified: August 19 2014
*/
/*******************************************************************************/
{
    int i;
    int ndof = (N_i + 2)*(N_j + 2);

    F = (double*)malloc(ndof*sizeof(double));
    G = (double*)malloc(ndof*sizeof(double));

    for(i = 0; i < ndof; i ++)
    {
        F[i] = 0.0;
        G[i] = 0.0;
    }

    return ;
}

/******************************************************************************************/
int set_F_and_G(double * U, double * V, double * F, double * G, double dx, double dy, double Re, double dt, double g_x, double g_y)
/*
  Purpose:
    Set the values of the F and G arrays
  Modified: August 19 2014
*/
/******************************************************************************************/
{
    int i, j;

    for(i = 1; i < N_i + 1; i++)
    {
        for(j = 1; j < N_j + 1; j++)
        {
            double d2U_dx = partial_squared_U_wrt_x(U, dx, i, j);
            double d2U_dy = partial_squared_U_wrt_y(U, dy, i, j);
            double dU2_dx = partial_U_Squared_wrt_x(U, dx, i, j);
            double dUV_dy = partial_UV_wrt_y(U, V, dy, i, j);

            F[INDEX(i, j)] = U[INDEX(i, j)] + (((d2U_dx + d2U_dy)/Re) - dU2_dx - dUV_dy + g_x)*dt;

            double d2V_dx = partial_squared_V_wrt_x(V, dx, i, j);
            double d2V_dy = partial_squared_V_wrt_y(V, dy, i, j);
            double dV2_dy = partial_V_Squared_wrt_y(V, dy, i, j);
            double dUV_dx = partial_UV_wrt_x(U, V, dx, i, j);

            G[INDEX(i, j)] = V[INDEX(i, j)] + (((d2V_dx + d2V_dy)/Re) - dV2_dy - dUV_dx + g_y)*dt;
        }
    }

    return 0;
}

/**************************************************************************/
int set_RHS(double * R, double * F, double * G, double dx, double dy, double dt)
/*
  Purpose:
    Set the R array values
  Modified: August 19 2014
*/
/**************************************************************************/
{
    int i, j;

    for(i = 1; i < N_i + 1; i++)
    {
        for(j = 1; j < N_j + 1; j++)
        {
            int i_minus = i - 1;
            int j_minus = j - 1;

            if(i == 1)
            {
                i_minus = N_i;
            }
            if(j == 1)
            {
                j_minus = N_j;
            }

            R[INDEX(i, j)] = (((F[INDEX(i, j)] - F[INDEX(i_minus, j)])/dx) +
                                ((G[INDEX(i, j)] - G[INDEX(i, j_minus)])/dy))/dt;
        }
    }

    return 0;
}

/***************************************************************************/
int updateField(double * U, double * V, double * F, double * G, double * P_new, double dx, double dy, double dt)
/*
  Purpose:
    Update the velocity field
  Modified: August 19 2014
*/
/***************************************************************************/
{
    int i, j;

    for(i = 1; i < N_i + 1; i++)
    {
        for(j = 1; j < N_j + 1; j++)
        {
            int i_plus = i + 1;
            int j_plus = j + 1;

            if(i == N_i)
            {
                i_plus = 1;
            }
            if(j == N_j)
            {
                j_plus = 1;
            }

            U[INDEX(i, j)] = F[INDEX(i, j)] - (dt/dx)*(P_new[INDEX(i_plus, j)] - P_new[INDEX(i, j)]);

            V[INDEX(i, j)] = G[INDEX(i, j)] - (dt/dy)*(P_new[INDEX(i, j_plus)] - P_new[INDEX(i, j)]);
        }
    }

    return 0;
}

