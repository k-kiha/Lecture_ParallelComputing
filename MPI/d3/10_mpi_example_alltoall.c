// mpicc -Itools 10_mpi_example_alltoall.c tools/mpi_tools.c
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include"./tools/mpi_tools.h"

int main(int argc, char *argv[]) {

    int n=16;
    int n1=n;
    int n2=n;
    
    //--- MPI
    int         rank, nprocs;
    MPI_INFO    mpi_info;
    int         n1sub, n2sub;
    double      r1sub,r2sub;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    mpi_tool_cart_2D(rank, nprocs, &mpi_info); /* << new */
    mpi_tool_setsubdomain_2D(rank, n1, n2, &mpi_info);
    
    n1sub = mpi_info.local_n1;
    n2sub = mpi_info.local_n2;

    double *x;
    x = (double *)malloc((n1sub)*(n2sub) * sizeof(double));

    for (int jj = 0; jj < n2sub; jj++) {
        for (int ii = 0; ii < n1sub; ii++) {
            x[ii*(n2sub) + jj] = rank*0.1 
                                 + 100*((ii+1) + mpi_info.ista -1)
                                 +   1*((jj+1) + mpi_info.jsta -1);
        }
    }


    for (int i = nprocs-1; i >= 0; i--)
    {
        if (rank==i){
            printf(".  === Rank %d ===\n", rank);
            for(int jj=n2sub-1; jj>=0; jj--) {
                printf(".  .  Befor::  ");
                for(int ii=0; ii<n1sub; ii++) {
                    printf("%6.1f ", x[ii*(n2sub) + jj]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //--- Alltoall
    int n1subsub, n2subsub;
    n1subsub = n1sub / mpi_info.split_y.nprocs;
    n2subsub = n2sub / mpi_info.split_x.nprocs;    

    double *sendbuf;
    double *recvbuf;
    sendbuf = (double *)malloc(n1sub*n2sub * sizeof(double));
    recvbuf = (double *)malloc(n1subsub*n2 * sizeof(double));

    int count = 0;
    for (int rrr = 0; rrr < mpi_info.split_y.nprocs; rrr++)
    {
        for (int i = 0; i < n1subsub; i++)
        for (int j = 0; j < n2sub   ; j++)
        {
            sendbuf[count] = x[(i+rrr*n1subsub)*(n2sub) + j];
            count++;
        }
    }

    MPI_Alltoall(sendbuf, n1subsub*n2sub, MPI_DOUBLE,
                 recvbuf, n1subsub*n2sub, MPI_DOUBLE,
                 mpi_info.split_y.subcomm);


    double *x_yline;
    x_yline = (double *)malloc((n1subsub)*(n2) * sizeof(double));
    count = 0;
    for (int rrr = 0; rrr < mpi_info.split_y.nprocs; rrr++)
    {
        for (int i = 0; i < n1subsub; i++)
        for (int j = 0; j < n2sub   ; j++)
        {
            x_yline[(i)*(n2) + j+rrr*n2sub] = recvbuf[count];
            count++;
        }
    }

    for (int i = nprocs-1; i >= 0; i--)
    {
        if (rank==i){
            printf(".  === Rank %d ===\n", rank);
            for(int jj=n2-1; jj>=0; jj--) {
                printf(".  .  After::  ");
                for(int ii=0; ii<n1subsub; ii++) {
                    printf("%6.1f ", x_yline[ii*(n2) + jj]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(x);
    free(x_yline);

    mpi_tool_clean(&mpi_info);
    MPI_Finalize();

    return 0;
}

