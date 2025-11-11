#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    MPI_Status status;
    MPI_Comm mpi_cart_comm;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    
    int dims[2]     = {2, 4};
    int periods[2]  = {0, 1};
    int cart_coords[2];
    int cart_rank;
    int cart_wrank,cart_erank;
    int cart_srank,cart_nrank;

    MPI_Comm mpi_subcomm_split_x;
    int split_x_nprocs,split_x_rank;
    int split_x_Lrank ,split_x_Rrank;

    MPI_Comm mpi_subcomm_split_y;
    int split_y_nprocs,split_y_rank;
    int split_y_Lrank ,split_y_Rrank;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &mpi_cart_comm);
    MPI_Cart_coords(mpi_cart_comm, rank, 2, cart_coords);
    MPI_Cart_rank  (mpi_cart_comm, cart_coords, &cart_rank);
    MPI_Cart_shift (mpi_cart_comm, 0, 1, &cart_wrank, &cart_erank);
    MPI_Cart_shift (mpi_cart_comm, 1, 1, &cart_srank, &cart_nrank);

    // Create sub-communicators for X and Y directions
    MPI_Comm_split(mpi_cart_comm, cart_coords[0], cart_coords[1], &mpi_subcomm_split_y);
    MPI_Comm_split(mpi_cart_comm, cart_coords[1], cart_coords[0], &mpi_subcomm_split_x);

    MPI_Comm_size(mpi_subcomm_split_x, &split_x_nprocs); //<<--2 
    MPI_Comm_size(mpi_subcomm_split_y, &split_y_nprocs); //<<--4

    MPI_Comm_rank(mpi_subcomm_split_x, &split_x_rank); //0~1
    MPI_Comm_rank(mpi_subcomm_split_y, &split_y_rank); //0~3

    MPI_Cart_create(mpi_subcomm_split_x, 1, &split_x_nprocs, &periods[0], 0, &mpi_subcomm_split_x);
    MPI_Cart_create(mpi_subcomm_split_y, 1, &split_y_nprocs, &periods[1], 0, &mpi_subcomm_split_y);
    MPI_Cart_shift(mpi_subcomm_split_x, 0, 1, &split_x_Lrank, &split_x_Rrank);
    MPI_Cart_shift(mpi_subcomm_split_y, 0, 1, &split_y_Lrank, &split_y_Rrank);

    for (int i = 0; i < nprocs; i++)
    {
        if(rank==i) {
            printf("=== Rank %d ", cart_rank);
            printf("   ::split (X,Y) nprocs: (%2d, %2d)", split_x_nprocs, split_y_nprocs);
            printf("   ::split (X,Y) rank: (%2d, %2d)", split_x_rank, split_y_rank);
            printf("---::X Left-Right : (%2d, %2d)", split_x_Lrank, split_x_Rrank);
            printf("---::Y Left-Right : (%2d, %2d)\n", split_y_Lrank, split_y_Rrank);
        };
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}