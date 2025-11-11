#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    MPI_Status status;
    MPI_Comm mpi_cart_comm;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    
    int dims[3]     = {3, 2, 2};
    int periods[3]  = {0, 1, 0};
    int cart_coords[3];
    int cart_rank;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &mpi_cart_comm);
    MPI_Cart_coords(mpi_cart_comm, rank, 3, cart_coords);
    MPI_Cart_rank  (mpi_cart_comm, cart_coords, &cart_rank);

    MPI_Comm mpi_subcomm_split_x;
    int split_x_nprocs,split_x_rank;
    int split_x_Lrank ,split_x_Rrank;

    MPI_Comm mpi_subcomm_split_y;
    int split_y_nprocs,split_y_rank;
    int split_y_Lrank ,split_y_Rrank;

    MPI_Comm mpi_subcomm_split_z;
    int split_z_nprocs,split_z_rank;
    int split_z_Lrank ,split_z_Rrank;

    int remain_dims_x[3] = {1, 0, 0};
    int remain_dims_y[3] = {0, 1, 0};
    int remain_dims_z[3] = {0, 0, 1};
    MPI_Cart_sub(mpi_cart_comm, remain_dims_x, &mpi_subcomm_split_x);
    MPI_Cart_sub(mpi_cart_comm, remain_dims_y, &mpi_subcomm_split_y);
    MPI_Cart_sub(mpi_cart_comm, remain_dims_z, &mpi_subcomm_split_z);

    MPI_Comm_size(mpi_subcomm_split_x, &split_x_nprocs);
    MPI_Comm_size(mpi_subcomm_split_y, &split_y_nprocs);
    MPI_Comm_size(mpi_subcomm_split_z, &split_z_nprocs);
    MPI_Comm_rank(mpi_subcomm_split_x, &split_x_rank);
    MPI_Comm_rank(mpi_subcomm_split_y, &split_y_rank);
    MPI_Comm_rank(mpi_subcomm_split_z, &split_z_rank);
    
    MPI_Cart_shift(mpi_subcomm_split_x, 0, 1, &split_x_Lrank, &split_x_Rrank);
    MPI_Cart_shift(mpi_subcomm_split_y, 0, 1, &split_y_Lrank, &split_y_Rrank);
    MPI_Cart_shift(mpi_subcomm_split_z, 0, 1, &split_z_Lrank, &split_z_Rrank);

    for (int i = 0; i < nprocs; i++)
    {
        if(rank==i) {
            printf("=== Rank %2d ", cart_rank);
            printf("   ::split (X,Y,Z) nprocs: (%2d, %2d, %2d)", split_x_nprocs, split_y_nprocs, split_z_nprocs);
            printf("   ::split (X,Y,Z) rank: (%2d, %2d, %2d)", split_x_rank, split_y_rank, split_z_rank);
            printf("---::X Left-Right : (%2d, %2d)", split_x_Lrank, split_x_Rrank);
            printf("---::Y Left-Right : (%2d, %2d)", split_y_Lrank, split_y_Rrank);
            printf("---::Z Left-Right : (%2d, %2d)\n", split_z_Lrank, split_z_Rrank);
        };
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}