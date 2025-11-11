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

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &mpi_cart_comm);

    MPI_Cart_coords(mpi_cart_comm, rank, 2, cart_coords);
    MPI_Cart_rank  (mpi_cart_comm, cart_coords, &cart_rank);

    MPI_Cart_shift (mpi_cart_comm, 0, 1, &cart_wrank, &cart_erank);
    MPI_Cart_shift (mpi_cart_comm, 1, 1, &cart_srank, &cart_nrank);

    for (int i = 0; i < nprocs; i++)
    {
        if(rank==i) {
            printf("=== Rank %d ", cart_rank);
            printf("---::map_coords:(%2d, %2d)", cart_coords[0], cart_coords[1]);
            printf("---::W-E : (%2d, %2d)", cart_wrank, cart_erank);
            printf("---::S-N : (%2d, %2d)\n", cart_srank, cart_nrank);
        };
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}