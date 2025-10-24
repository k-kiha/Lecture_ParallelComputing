// mpicc -Itools 05_mpi_jacobi_check.c tools/mpi_tools.c
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
    
    for (int i = nprocs-1; i >= 0 ; i--)
    {
        if(rank==i) mpi_tool_monitor(mpi_info);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    mpi_tool_clean(&mpi_info);
    MPI_Finalize();

    return 0;
}
