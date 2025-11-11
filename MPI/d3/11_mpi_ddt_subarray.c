// mpicc -Itools 11_mpi_ddt_subarray.c tools/mpi_tools.c
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include"./tools/mpi_tools.h"

int main(int argc, char *argv[]) {

    int n=16;
    int n1=n;
    int n2=n;
    double *x;
    
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

    x = (double *)malloc((n1sub+2)*(n2sub+2) * sizeof(double));

    for (int jj = 1; jj < n2sub+1; jj++) {
        for (int ii = 1; ii < n1sub+1; ii++) {
            x[ii*(n2sub+2) + jj] = rank*0.1 
                                 + 100*(ii + mpi_info.ista -1)
                                 +   1*(jj + mpi_info.jsta -1);
        }
    }

    for (int i = nprocs-1; i >= 0; i--)
    {
        if (rank==i){
            printf("=== Rank %d ===\n", rank);
            for(int jj=n2sub+1; jj>=0; jj--) {
                printf("Befor::  ");
                for(int ii=0; ii<n1sub+2; ii++) {
                    printf("%6.1f ", x[ii*(n2sub+2) + jj]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    int ndims =2;
    int send_sizes[2]    = {8+2,4+2};
    int send_subsizes[2] = {8+2,  1};
    int send_starts[2]   = {  0,  4};
    int recv_sizes[2]    = {8+2,4+2};
    int recv_subsizes[2] = {8+2,  1};
    int recv_starts[2]   = {  0,  0};

    MPI_Datatype mytype_send;
    MPI_Datatype mytype_recv;

    MPI_Type_create_subarray(ndims,
                            send_sizes, 
                            send_subsizes, 
                            send_starts,
                            MPI_ORDER_C, MPI_DOUBLE, &mytype_send);

    MPI_Type_create_subarray(ndims,
                            recv_sizes, 
                            recv_subsizes, 
                            recv_starts,
                            MPI_ORDER_C, MPI_DOUBLE, &mytype_recv);

    MPI_Type_commit(&mytype_send);
    MPI_Type_commit(&mytype_recv);

    MPI_Sendrecv(x, 1, mytype_send, mpi_info.split_y.R, 11,
                 x, 1, mytype_recv, mpi_info.split_y.L, 11,
                 mpi_info.split_y.subcomm, MPI_STATUS_IGNORE);

    MPI_Type_free(&mytype_send);
    MPI_Type_free(&mytype_recv);


    for (int i = nprocs-1; i >= 0; i--)
    {
        if (rank==i){
            printf("=== Rank %d ===\n", rank);
            for(int jj=n2sub+1; jj>=0; jj--) {
                printf("After::  ");
                for(int ii=0; ii<n1sub+2; ii++) {
                    printf("%6.1f ", x[ii*(n2sub+2) + jj]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(x);

    mpi_tool_clean(&mpi_info);
    MPI_Finalize();

    return 0;
}