// mpicc -Itools 12_mpi_ddt_subarray_alltoall.c tools/mpi_tools.c
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

    int n1subsub, n2subsub;
    n1subsub = n1sub / mpi_info.split_y.nprocs;
    n2subsub = n2sub / mpi_info.split_x.nprocs; 

    double *x;
    x = (double *)malloc((n1sub)*(n2sub) * sizeof(double));
    double *x_yline;
    x_yline = (double *)malloc((n1subsub)*(n2) * sizeof(double));

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
    MPI_Datatype mytype_send[mpi_info.split_y.nprocs];
    MPI_Datatype mytype_recv[mpi_info.split_y.nprocs];
    int sendcounts[mpi_info.split_y.nprocs];
    int recvcounts[mpi_info.split_y.nprocs];
    int senddisps[mpi_info.split_y.nprocs];
    int recvdisps[mpi_info.split_y.nprocs];

    for (int rrr = 0; rrr < mpi_info.split_y.nprocs; rrr++){
        int send_sizes[2]    = {    8,4}; // {    n1sub   ,n2sub};
        int send_subsizes[2] = {    2,4}; // {    n1subsub,n2sub};
        int send_starts[2]   = {rrr*2,0}; // {rrr*n1subsub,    0};

        MPI_Type_create_subarray(2, send_sizes, send_subsizes, send_starts, MPI_ORDER_C, MPI_DOUBLE, &mytype_send[rrr]);
        MPI_Type_commit(&mytype_send[rrr]);

        int recv_sizes[2]    = {2,   16}; // { n1subsub,    n2   };
        int recv_subsizes[2] = {2,    4}; // { n1subsub,    n2sub};
        int recv_starts[2]   = {0,rrr*4}; // {        0,rrr*n2sub};

        MPI_Type_create_subarray(2, recv_sizes, recv_subsizes, recv_starts, MPI_ORDER_C, MPI_DOUBLE, &mytype_recv[rrr]);
        MPI_Type_commit(&mytype_recv[rrr]);

        sendcounts[rrr] =1;
        recvcounts[rrr] =1;
        senddisps[rrr] =0;
        recvdisps[rrr] =0;

    }    

    MPI_Alltoallw(x      , sendcounts, senddisps, mytype_send,
                  x_yline, recvcounts, recvdisps, mytype_recv,
                  mpi_info.split_y.subcomm);

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

    for (int rrr = 0; rrr < mpi_info.split_y.nprocs; rrr++){
        MPI_Type_free(&mytype_send[rrr]);
        MPI_Type_free(&mytype_recv[rrr]);
    }
    mpi_tool_clean(&mpi_info);
    MPI_Finalize();

    return 0;
}

