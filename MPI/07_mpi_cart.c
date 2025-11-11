// mpicc -Itools 07_mpi_cart.c tools/mpi_tools.c
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

    mpi_tool_cart_2D(rank, nprocs, &mpi_info);
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



// void mpi_tool_cart_2D(int rank, int nprocs, MPI_INFO *mpi_info){
//     mpi_info->rank = rank;
//     mpi_info->nprocs = nprocs;

//     mpi_info->map_dims[0] = 2;
//     mpi_info->map_dims[1] = 4;
//     mpi_info->map_periods[0] = 0;
//     mpi_info->map_periods[1] = 1;

//     MPI_Cart_create(MPI_COMM_WORLD, 2, mpi_info->map_dims, mpi_info->map_periods, 0, &mpi_info->cart_comm); 

//     MPI_Cart_coords(mpi_info->cart_comm, rank, 2,               mpi_info->map_coords);
//     MPI_Cart_rank  (mpi_info->cart_comm, mpi_info->map_coords, &mpi_info->rank);

//     MPI_Cart_shift(mpi_info->cart_comm, 0, 1, &mpi_info->w, &mpi_info->e);
//     MPI_Cart_shift(mpi_info->cart_comm, 1, 1, &mpi_info->s, &mpi_info->n);

//     // Create sub-communicators for X and Y directions
//     MPI_Comm_split(mpi_info->cart_comm, mpi_info->map_coords[1], mpi_info->map_coords[0], &mpi_info->split_x.subcomm);
//     MPI_Comm_split(mpi_info->cart_comm, mpi_info->map_coords[0], mpi_info->map_coords[1], &mpi_info->split_y.subcomm);

//     MPI_Comm_rank(mpi_info->split_x.subcomm, &mpi_info->split_x.rank);
//     MPI_Comm_rank(mpi_info->split_y.subcomm, &mpi_info->split_y.rank);

//     MPI_Comm_size(mpi_info->split_x.subcomm, &mpi_info->split_x.nprocs);
//     MPI_Comm_size(mpi_info->split_y.subcomm, &mpi_info->split_y.nprocs);

//     MPI_Cart_create(mpi_info->split_x.subcomm, 1, &mpi_info->split_x.nprocs, &mpi_info->map_periods[0], 0, &mpi_info->split_x.subcomm); 
//     MPI_Cart_create(mpi_info->split_y.subcomm, 1, &mpi_info->split_y.nprocs, &mpi_info->map_periods[1], 0, &mpi_info->split_y.subcomm); 

//     MPI_Cart_shift(mpi_info->split_x.subcomm, 0, 1, &mpi_info->split_x.L, &mpi_info->split_x.R);
//     MPI_Cart_shift(mpi_info->split_y.subcomm, 0, 1, &mpi_info->split_y.L, &mpi_info->split_y.R);

// }