#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "parallel.h"
#include "derivative.h"

/*************************************************************************************/
void make_domains(int num_procs)
/*
  Purpose:
    MAKE_DOMAINS sets up the information defining the process domains.
  Modified:
    August 19 2014:  Accounts for periodic boundary conditions.
  Parameters:
    Input, int NUM_PROCS, the number of processes.
*/
/**************************************************************************************/
{
    double d;
    double eps;
    int i;
    int p;
    double x_max;
    double x_min;

    //---ALLOCATE ARRAYS FOR PROCESS INFORMATION---//

    proc = (int*)malloc((N_i + 2)*sizeof(int));
    i_min = (int*)malloc(num_procs*sizeof(int));
    i_max = (int*)malloc(num_procs*sizeof(int));
    left_proc = (int*)malloc(num_procs*sizeof(int));
    right_proc = (int*)malloc(num_procs*sizeof(int));

    //---DIVIDE THE RANGE [(-eps+1)..(N_i+eps)] EVENLY AMONG PROCESSES---//

    eps = 0.0001;
    d = (N_i - 1.0 + 2.0*eps)/(double)num_procs;

    for (p = 0; p < num_procs; p++)
    {
        //---PROCESS DOMAIN---//

        x_min = -eps + 1.0 + (double)(p*d);
        x_max = x_min + d;

        //---IDENTIFY VERTICES---//

        for(i = 1; i <= N_i; i++)
        {
            if(x_min <= i && i < x_max)
            {
                proc[i] = p;
            }
        }
    }

    for(p = 0; p < num_procs; p++)
    {
        for(i = 1; i <= N_i; i++)       //Find the smallest vertex index in the domain.
        {
            if(proc[i] == p)
            {
                break;
            }
        }
        i_min[p] = i;

        for(i = N_i; 1 <= i; i--)        //Find the largest vertex index in the domain.
        {
            if(proc[i] == p)
            {
                break;
            }
        }
        i_max[p] = i + 1;

        //printf("Process %d has i_min = %d, i_max = %d\n", p, i_min[p], i_max[p]);

        //---FIND PROCESSES TO THE LEFT AND RIGHT---//

        left_proc[p] = num_procs - 1;                   //Left and right boundaries are periodic.
        right_proc[p] = 0;

        if(proc[p] != -1)
        {
            if(1 < i_min[p] && i_min[p] <= N_i)
            {
                left_proc[p] = proc[i_min[p] - 1];
            }
            if( (0 < (i_max[p] - 1)) && ((i_max[p] - 1) < N_i))
            {
                right_proc[p] = proc[i_max[p]];
            }
        }

        //printf("Process %d has left = %d, right = %d\n", p, left_proc[p], right_proc[p]);
    }
    return;
}

/***********************************************************************************************/
void jacobi(int num_procs, double * R, double * P, double * P_new, double dx, double dy)
/*
  Purpose:
    JACOBI carries out the Jacobi iteration for the linear system.
  Modified:
    August 19 2014:  Accounts for periodic boundary conditions.
  Parameters:
    Input, int NUM_PROCS, the number of processes.
    Input, double F[(N_i+2)*(N_i+2)], the right hand side of the linear system.
  Notes:
    NOT CURRENTLY USED!
*/
/************************************************************************************************/
{
    int i, j;
    MPI_Request request[4];
    int requests;
    MPI_Status status[4];

    //---UPDATE TOP-BOTTOM GHOST LAYERS---//

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        P[INDEX(i,0)] = P[INDEX(i, N_j)];
        P[INDEX(i, N_j + 1)] = P[INDEX(i, 1)];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //---ADD SIDE GHOST LAYERS USING NON-BLOCKING SEND AND RECEIVE---//

    requests = 0;

    if(left_proc[my_rank] >= 0 && left_proc[my_rank] < num_procs)               //This doesn't change the top/bottom
    {
        MPI_Irecv(P + INDEX(i_min[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 0, MPI_COMM_WORLD,
                  request + requests++);

        MPI_Isend(P + INDEX(i_min[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 1, MPI_COMM_WORLD,
                  request + requests++);
    }

    if(right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs)                //Check validity
    {
        MPI_Irecv(P + INDEX(i_max[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 1, MPI_COMM_WORLD,
                  request + requests++);

        MPI_Isend(P + INDEX(i_max[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 0, MPI_COMM_WORLD,
                  request + requests++);
    }

    MPI_Waitall(requests, request, status);       //Wait for all non-blocking communications to complete before updating boundaries.
    MPI_Barrier(MPI_COMM_WORLD);

    //---JACOBI UPDATE FOR ALL NODES---//

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        for (j = 1; j < N_j + 1; j++)
        {
            P_new[INDEX(i,j)] =
                0.25*(P[INDEX(i-1,j)] + P[INDEX(i+1,j)] +
                      P[INDEX(i,j-1)] + P[INDEX(i,j+1)] +
                      dx*dy*R[INDEX(i,j)]);
        }
    }

    return;
}

/***********************************************************************************************/
void SOR(int num_procs, double * R, double * P, double * P_new, double dx, double dy, double omega)
/*
  Purpose:
    SOR carries out the SOR iteration for the linear system.
  Modified:
    August 28 2014:  Accounts for periodic boundary conditions.
  Parameters:
    Input, int NUM_PROCS, the number of processes.
    Input, double F[(N_i+2)*(N_i+2)], the right hand side of the linear system.
*/
/************************************************************************************************/
{
    int i, j;
    MPI_Request request[4];
    int requests;
    MPI_Status status[4];

    //---UPDATE TOP-BOTTOM GHOST LAYERS---//

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        P[INDEX(i,0)] = P[INDEX(i, N_j)];
        P[INDEX(i, N_j + 1)] = P[INDEX(i, 1)];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //---ADD SIDE GHOST LAYERS USING NON-BLOCKING SEND AND RECEIVE---//

    requests = 0;

    if(left_proc[my_rank] >= 0 && left_proc[my_rank] < num_procs)
    {
        MPI_Irecv(P + INDEX(i_min[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 0, MPI_COMM_WORLD,
                  request + requests++);

        MPI_Isend(P + INDEX(i_min[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 1, MPI_COMM_WORLD,
                  request + requests++);
    }

    if(right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs)                //Check validity
    {
        MPI_Irecv(P + INDEX(i_max[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 1, MPI_COMM_WORLD,
                  request + requests++);

        MPI_Isend(P + INDEX(i_max[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 0, MPI_COMM_WORLD,
                  request + requests++);
    }

    MPI_Waitall(requests, request, status);       //Wait for all non-blocking communications to complete before updating boundaries.
    MPI_Barrier(MPI_COMM_WORLD);

    //---SOR UPDATE FOR ALL NODES---//

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        for (j = 1; j < N_j + 1; j++)
        {
            P_new[INDEX(i, j)] = (P[INDEX(i,j)])*(1.0 - omega) + (omega/((2.0/(pow(dx, 2.0))) +
                                 (2.0/(pow(dy, 2.0)))))*(((P[INDEX(i + 1, j)] +
                                         P_new[INDEX(i - 1, j)])/(pow(dx, 2.0))) + ((P[INDEX(i, j + 1)] +
                                                 P_new[INDEX(i, j - 1)])/(pow(dy, 2.0))) - R[INDEX(i, j)]);
        }
    }

    return;
}

/*************************************************************************************************/
void set_F_G_and_R(int num_procs, double * F, double * G, double * R, double * U, double * V, double dx, double dy, double Re, double dt, double g_x, double g_y)
/*
  Purpose:
    Caculate the F, G and R fields in parallel
  Modifed:
    August 2014
*/
/*************************************************************************************************/
{
    int i, j;

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        for (j = 1; j < N_j + 1; j++)
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

    //---UPDATE TOP-BOTTOM GHOST LAYERS---//

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        F[INDEX(i,0)] = F[INDEX(i, N_j)];
        F[INDEX(i, N_j + 1)] = F[INDEX(i, 1)];

        G[INDEX(i,0)] = G[INDEX(i, N_j)];
        G[INDEX(i, N_j + 1)] = G[INDEX(i, 1)];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //---ADD SIDE GHOST LAYERS USING NON-BLOCKING SEND AND RECEIVE---//

    MPI_Request F_request[4];
    int F_requests = 0;
    MPI_Status F_status[4];

    MPI_Request G_request[4];
    int G_requests = 0;
    MPI_Status G_status[4];

    if(left_proc[my_rank] >= 0 && left_proc[my_rank] < num_procs)
    {
        MPI_Irecv(F + INDEX(i_min[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 0, MPI_COMM_WORLD,
                  F_request + F_requests++);

        MPI_Isend(F + INDEX(i_min[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 1, MPI_COMM_WORLD,
                  F_request + F_requests++);

        MPI_Irecv(G + INDEX(i_min[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 2, MPI_COMM_WORLD,
                  G_request + G_requests++);

        MPI_Isend(G + INDEX(i_min[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  left_proc[my_rank], 3, MPI_COMM_WORLD,
                  G_request + G_requests++);
    }

    if(right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs)                //Check validity
    {
        MPI_Irecv(F + INDEX(i_max[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 1, MPI_COMM_WORLD,
                  F_request + F_requests++);

        MPI_Isend(F + INDEX(i_max[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 0, MPI_COMM_WORLD,
                  F_request + F_requests++);

        MPI_Irecv(G + INDEX(i_max[my_rank], 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 3, MPI_COMM_WORLD,
                  G_request + G_requests++);

        MPI_Isend(G + INDEX(i_max[my_rank] - 1, 0), (N_j + 2), MPI_DOUBLE,
                  right_proc[my_rank], 2, MPI_COMM_WORLD,
                  G_request + G_requests++);
    }

    MPI_Waitall(F_requests, F_request, F_status);       //Wait for all non-blocking communications to complete before updating boundaries.
    MPI_Waitall(G_requests, G_request, G_status);

    for(i = i_min[my_rank]; i < i_max[my_rank]; i++)
    {
        for (j = 1; j < N_j + 1; j++)
        {
            R[INDEX(i, j)] = (((F[INDEX(i, j)] - F[INDEX(i - 1, j)])/dx) +
                                ((G[INDEX(i, j)] - G[INDEX(i, j - 1)])/dy))/dt;
        }
    }

    return;
}

/*************************************************************************************************/
void timestamp()
/*
  Purpose:
    TIMESTAMP prints the current YMDHMS date as a time stamp.
  Example:
    31 May 2001 09:45:54 AM
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 September 2003
  Author:
    John Burkardt
  Parameters:
    None
*/
/*************************************************************************************************/
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p\n", tm);

    printf("%s\n", time_buffer);

    return;
# undef TIME_SIZE
}

