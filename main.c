#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "parallel.h"
#include "algorithm.h"
#include "derivative.h"

/**********************************************************************************************/
    /*
    PARALEL NAVIER-STOKES SOLVER
    EMPLOYS THE SOR ALGORITHM TO SOLVE THE POISSON EQUATION


    MODIFIED: 29/08/14
    AUTHOR: SAM WHYMAN + see refs for credits.

    MODIFIED: 22/09/15
    Redundant pieces of code were removed.
    */
/**********************************************************************************************/

int N_i = 51, N_j = 51;			        //Number of internal grid points

#define INDEX(i,j) ((N_i+2)*(i)+(j))    //Macro to index into a 2-D (N_i+2)x(N_j+2) array

int my_rank;			                //Rank of this process

int * proc;			                    //Process indexed by vertex
int * i_min, * i_max;		            //Min, Max vertex indices of processes
int * left_proc, * right_proc;	        //Processes to left and right

int main(int argc, char * argv[])
{
    double * P, * P_new;                    //Pressure fields
    double * U, * V;                        //Velocity field
    double * F, * G;                        //Derivative fields
    double * R;                             //Source term

    int i, j;                               //Indices

    int i_cent = (int)(0.5*(N_i + 1));      //Centres
    int j_cent = (int)(0.5*(N_j + 1));

    double dx = 0.1;                        //Node separation
    double dy = 0.1;

    double change;                          //For convergence test
    double epsilon = 1.0E-3;
    double omega = 0.9;                     //SOR Relaxation parameter

    double Re = 100.00;                     //Reynolds number
    double t = 0.00, tMax = 1.0;
    double dt = 0.01;
    double g_x = 0.0, g_y = 0.0;            //Body forces

    double my_change;                       //For Poisson algorithm
    int my_n;
    int n;
    int num_procs;
    int sor_step = 0, step_max = 10000;
    double * swap;

    double wall_time;

    FILE * deviation_file;

    char flow_file_name[100];
    char mag_file_name[100];
    char dev_file_name[100];                  //For deviation data

    //---MPI INITIALISTION---//

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    make_domains(num_procs);                //Divides into strips

    MPI_Barrier(MPI_COMM_WORLD);

    strcpy(mag_file_name, "mag_field.txt");
    strcpy(flow_file_name, "flow_field.txt");
    strcpy(dev_file_name, "deviation.txt");

    //---PRINT INFORMATION---//

    if (my_rank == 0)
    {

        printf("\n");
        printf("2D NAVIER STOKES:\n");
        printf("  C version\n");
        printf("  2-D Navier Stokes solution with SOR Poisson equation\n");
        printf("  ===========================================\n");
        printf("  MPI version: 2-D domains, non-blocking send/receive\n");
        printf("  Number of processes         = %d\n", num_procs);
        printf("  Number of vertices = %d, %d\n", N_i, N_j);
        printf("  Desired fractional accuracy = %lf\n", epsilon);
        printf("  Maximum time = %.2lf s\n", tMax);
        printf("\n" );

        deviation_file = fopen(dev_file_name, "w");
    }

    //---ALLOCATE AND INITIALISE---//

    P = allocate_NULL();
    P_new = allocate_NULL();
    init_pressure(P);

    U = allocate_NULL();
    V = allocate_NULL();
    R = allocate_NULL();
    F = allocate_NULL();
    G = allocate_NULL();

    initialise_theta_function_field(U, V, dx, dy, 10, 10, 0.5, 6.0, 6.0);     //nodes = 10, A = 0.5, a = b = 6.0
    wall_time = MPI_Wtime();                                                  //Begin clock

    //---BEGIN ITERATION---//

    printf("START: ");
    timestamp();
    printf("================================\n");

    while(t < tMax)
    {
        if(my_rank == 0)
        {
            //---TEST DEVIATIONS---///


            if((sin(t*M_PI) > -0.001) && (sin(t*M_PI) < 0.001))
            {
                double mag_test = sqrt(pow(U[INDEX(i_cent, j_cent)], 2) + pow(V[INDEX(i_cent, j_cent)], 2));
                fprintf(deviation_file, "%lf %.15lf\n", t, mag_test);
                printf("Current Time: %lf\n", t);
            }


        }

        MPI_Bcast(P, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);               //Send arrays from main
        MPI_Bcast(P_new, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(U, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(V, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        sor_step = 0;                                                 //These must be initialised at each time step
        change = 0.0;
        n = 0;
        my_change = 0.0;
        my_n = 0;

        set_F_G_and_R(num_procs, F, G, R, U, V, dx, dy, Re, dt, g_x, g_y);                //Obtain F, G and R in parallel

        MPI_Barrier(MPI_COMM_WORLD);

        //---COPY F AND G TO ALL PROCESSES VIA SEND & RECEIVE---//

        MPI_Status F_status;
        MPI_Status G_status;
        MPI_Status R_status;

        if(my_rank != 0)
        {
            MPI_Send(F + INDEX(i_min[my_rank], 1), (i_max[my_rank] - i_min[my_rank])*(N_j + 2),
                     MPI_DOUBLE, 0, 43, MPI_COMM_WORLD);
            MPI_Send(G + INDEX(i_min[my_rank], 1), (i_max[my_rank] - i_min[my_rank])*(N_j + 2),
                     MPI_DOUBLE, 0, 44, MPI_COMM_WORLD);
            MPI_Send(R + INDEX(i_min[my_rank], 1), (i_max[my_rank] - i_min[my_rank])*(N_j + 2),
                     MPI_DOUBLE, 0, 45, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if((my_rank == 0) && (num_procs > 1))
        {
            int FGR_sender;
            for(FGR_sender = 1; FGR_sender < num_procs; FGR_sender++)
            {
                MPI_Recv(F + INDEX(i_min[FGR_sender], 1), (i_max[FGR_sender] - i_min[FGR_sender])*(N_j + 2),
                         MPI_DOUBLE, FGR_sender, 43, MPI_COMM_WORLD, &F_status);
                MPI_Recv(G + INDEX(i_min[FGR_sender], 1), (i_max[FGR_sender] - i_min[FGR_sender])*(N_j + 2),
                         MPI_DOUBLE, FGR_sender, 44, MPI_COMM_WORLD, &G_status);
                MPI_Recv(R + INDEX(i_min[FGR_sender], 1), (i_max[FGR_sender] - i_min[FGR_sender])*(N_j + 2),
                         MPI_DOUBLE, FGR_sender, 45, MPI_COMM_WORLD, &R_status);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(F, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);               //Send F, G and R arrays from main
        MPI_Bcast(G, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(R, (N_i + 2)*(N_j + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //---SOLVE POISSON EQUATION---//

        do
        {
            SOR(num_procs, R, P, P_new, dx, dy, omega);
            sor_step++;

            //---ESTIMATE ERROR---//

            change = 0.0;
            n = 0;

            my_change = 0.0;
            my_n = 0;

            for (i = i_min[my_rank]; i < i_max[my_rank]; i++)
            {
                for (j = 1; j < N_j + 1; j++)
                {
                    if (P_new[INDEX(i,j)] != 0.0)
                    {
                        my_change = my_change
                                    + fabs(1.0 - P[INDEX(i,j)]/P_new[INDEX(i,j)]);
                        my_n++;
                    }
                }
            }
            MPI_Allreduce(&my_change, &change, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&my_n, &n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (n != 0)
            {
                change = change/n;
            }

            swap = P;               //Swap P and P_new
            P = P_new;
            P_new = swap;

            if(sor_step > step_max)
            {
                printf("Did not converge within %d iterations!\n", sor_step);
            }
        }
        while(epsilon < change);

        MPI_Barrier(MPI_COMM_WORLD);

        //---COPY TO MAIN VIA SEND & RECEIVE---//

        MPI_Status P_status;

        if(my_rank != 0)
        {
            MPI_Send(P_new + INDEX(i_min[my_rank], 1), (i_max[my_rank] - i_min[my_rank])*(N_j + 2),
                     MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if((my_rank == 0) && (num_procs > 1))
        {
            int sender;
            for(sender = 1; sender < num_procs; sender++)
            {
                MPI_Recv(P_new + INDEX(i_min[sender], 1), (i_max[sender] - i_min[sender])*(N_j + 2),
                         MPI_DOUBLE, sender, 42, MPI_COMM_WORLD, &P_status);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if(my_rank == 0)
        {
            updateField(U, V, F, G, P_new, dx, dy, dt);
        }

        t += dt;

        if( ((int)(t*100)) % 50 == 0  ){
            printf("Current time: %.1f\n", t);
        }

    }

    //---GET WALLTIME AND PRINT TO FILE---//

    wall_time = MPI_Wtime() - wall_time;

    if(my_rank == 0)
    {
        printf ("\n");
        printf ("  Wall clock time = %f secs\n", wall_time);
        printf("================================\n");

        FILE * magPointer;
        FILE * flowPointer;
        magPointer = fopen(mag_file_name, "w");
        flowPointer = fopen(flow_file_name, "w");

        for(i = 1; i < N_i + 1; i++)
        {
            for(j = 1; j < N_j + 1; j++)
            {
                double press = P_new[INDEX(i, j)];
                double mag = sqrt(pow(U[INDEX(i, j)], 2) + pow(V[INDEX(i, j)], 2));
                fprintf(magPointer, "%lf %lf %.30lf\n", (double)(i - 1)*dx, (double)(j - 1)*dy, mag);
                fprintf(flowPointer, "%lf %lf %.10lf %.10lf\n", (double)(i - 1)*dx, (double)(j - 1)*dy, U[INDEX(i, j)], V[INDEX(i, j)]);
            }
        }
        fclose(magPointer);
        fclose(flowPointer);
    }

    //---TERMINATE AND FREE---//

    MPI_Finalize();

    free(U);
    free(V);
    free(F);
    free(G);
    free(R);
    free(P);
    free(P_new);

    if(my_rank == 0)
    {
        fclose(deviation_file);                         //Close file for deviation data
        printf ("\n");
        printf ("NAVIER_STOKES_MPI:\n");
        printf ("  Normal end of execution.\n");
        printf ("  ");
        timestamp();
    }

    return 0;
}
