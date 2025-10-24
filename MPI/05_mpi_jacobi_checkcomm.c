// mpicc -Itools 05_mpi_jacobi_checkcomm.c tools/mpi_tools.c
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
    int rank, nprocs;
    MPI_INFO mpi_info;
    int n1sub, n2sub;
    double r1sub,r2sub;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    mpi_tool_setmap_2D(rank, nprocs, &mpi_info);
    mpi_tool_setsubdomain_2D(rank, n1, n2, &mpi_info);
    
    //
    double *x;
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


    mpi_tool_communicate_w2e(mpi_info, x, n1sub, n2sub);
    mpi_tool_communicate_e2w(mpi_info, x, n1sub, n2sub);
    mpi_tool_communicate_s2n(mpi_info, x, n1sub, n2sub);
    mpi_tool_communicate_n2s(mpi_info, x, n1sub, n2sub);

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

