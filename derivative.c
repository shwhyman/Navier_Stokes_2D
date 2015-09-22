#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "derivative.h"

/****************************************************************/
/*
THE FOLLOWING FUNCTIONS ARE USED TO CALCULATE THE DERIVATIVES
REQUIRED FOR SOLVING THE NAVIER STOKES EQUATIONS.
THE PURPOSE OF EACH IS INCLUDED AS THE FUNCTION NAME.
*/
/****************************************************************/

double partial_U_Squared_wrt_x(double * U, double dx, int i, int j)      //  (d/dx)(U^2).  Ie, use for dU^2/dx at the point (i,j)
{
    int i_plus = i + 1;
    int i_minus = i - 1;

    if(i == N_i)
    {
        i_plus = 1;
    }
    if(i == 1)
    {
        i_minus = N_i;
    }

    double first = pow(((U[INDEX(i, j)] + U[INDEX(i_plus, j)])*0.5), 2.0);
    double second = pow(((U[INDEX(i_minus, j)] + U[INDEX(i, j)])*0.5), 2.0);

    double derivative = (first - second)/dx;

    return derivative;
}

double partial_V_Squared_wrt_y(double * V, double dy, int i, int j)      //  (d/dy)(V^2).  Ie, use for dV^2/dy at the point (i,j)
{
    int j_plus = j + 1;
    int j_minus = j - 1;

    if(j == N_j)
    {
        j_plus = 1;
    }
    if(j == 1)
    {
        j_minus = N_j;
    }

    double first = pow(((V[INDEX(i, j)] + V[INDEX(i, j_plus)])*0.5), 2.0);
    double second = pow(((V[INDEX(i, j_minus)] + V[INDEX(i, j)])*0.5), 2.0);

    double derivative = (first - second)/dy;

    return derivative;
}

double partial_UV_wrt_x(double * U, double * V, double dx, int i, int j)
{
    int i_plus = i + 1;
    int i_minus = i - 1;
    int j_plus = j + 1;

    if(i == N_i)
    {
        i_plus = 1;
    }
    if(i == 1)
    {
        i_minus = N_i;
    }
    if(j == N_j)
    {
        j_plus = 1;
    }

    double first = ((U[INDEX(i, j)] + U[INDEX(i, j_plus)])*0.5)*((V[INDEX(i, j)] + V[INDEX(i_plus, j)])*0.5);
    double second = ((U[INDEX(i_minus, j)] + U[INDEX(i_minus, j_plus)])*0.5)*((V[INDEX(i_minus, j)] + V[INDEX(i, j)])*0.5);

    double derivative = (first - second)/dx;

    return derivative;
}

double partial_UV_wrt_y(double * U, double * V, double dy, int i, int j)
{
    int i_plus = i + 1;
    int j_plus = j + 1;
    int j_minus = j - 1;

    if(i == N_i)
    {
        i_plus = 1;
    }
    if(j == N_j)
    {
        j_plus = 1;
    }
    if(j == 1)
    {
        j_minus = N_j;
    }

    double first = ((U[INDEX(i, j_plus)] + U[INDEX(i, j)])*0.5)*((V[INDEX(i, j)] + V[INDEX(i_plus, j)])*0.5);
    double second = ((U[INDEX(i, j_minus)] + U[INDEX(i, j)])*0.5)*((V[INDEX(i, j_minus)] + V[INDEX(i_plus, j_minus)])*0.5);

    double derivative = (first - second)/dy;

    return derivative;
}

double partial_squared_U_wrt_x(double * U, double dx, int i, int j)
{
    int i_plus = i + 1;
    int i_minus = i - 1;

    if(i == N_i)
    {
        i_plus = 1;
    }
    if(i == 1)
    {
        i_minus = N_i;
    }

    double derivative = (U[INDEX(i_plus, j)] - 2*U[INDEX(i, j)] + U[INDEX(i_minus, j)])/(pow(dx, 2.0));

    return derivative;
}

double partial_squared_V_wrt_x(double * V, double dx, int i, int j)
{
    int i_plus = i + 1;
    int i_minus = i - 1;

    if(i == N_i)
    {
        i_plus = 1;
    }
    if(i == 1)
    {
        i_minus = N_i;
    }

    double derivative = (V[INDEX(i_plus, j)] - 2*V[INDEX(i, j)] + V[INDEX(i_minus, j)])/(pow(dx, 2.0));

    return derivative;
}

double partial_squared_U_wrt_y(double * U, double dy, int i, int j)
{
    int j_plus = j + 1;
    int j_minus = j - 1;

    if(j == N_j)
    {
        j_plus = 1;
    }
    if(j == 1)
    {
        j_minus = N_j;
    }

    double derivative = (U[INDEX(i, j_plus)] - 2*U[INDEX(i, j)] + U[INDEX(i, j_minus)])/(pow(dy, 2.0));

    return derivative;
}

double partial_squared_V_wrt_y(double * V, double dy, int i, int j)
{
    int j_plus = j + 1;
    int j_minus = j - 1;

    if(j == N_j)
    {
        j_plus = 1;
    }
    if(j == 1)
    {
        j_minus = N_j;
    }

    double derivative = (V[INDEX(i, j_plus)] - 2*V[INDEX(i, j)] + V[INDEX(i, j_minus)])/(pow(dy, 2.0));

    return derivative;
}

