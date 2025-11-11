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

    MPI_Comm mpi_subcomm_split_xy, mpi_subcomm_split_yz;
    int split_xy_rank, split_yz_rank;
    int split_xy_coord[2], split_yz_coord[2];

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &mpi_cart_comm);
    MPI_Cart_coords(mpi_cart_comm, rank, 3, cart_coords);
    MPI_Cart_rank  (mpi_cart_comm, cart_coords, &cart_rank);

    MPI_Comm_split(mpi_cart_comm, cart_coords[2], rank, &mpi_subcomm_split_xy);
    MPI_Comm_split(mpi_cart_comm, cart_coords[0], rank, &mpi_subcomm_split_yz);

    MPI_Cart_create(mpi_subcomm_split_xy, 2, &dims[0], &periods[0], 0, &mpi_subcomm_split_xy);
    MPI_Cart_create(mpi_subcomm_split_yz, 2, &dims[1], &periods[1], 0, &mpi_subcomm_split_yz);

    MPI_Comm_rank(mpi_subcomm_split_xy, &split_xy_rank);
    MPI_Comm_rank(mpi_subcomm_split_yz, &split_yz_rank);

    MPI_Cart_coords(mpi_subcomm_split_xy, split_xy_rank, 2, split_xy_coord);
    MPI_Cart_coords(mpi_subcomm_split_yz, split_yz_rank, 2, split_yz_coord);


    MPI_Comm mpi_subcomm_split_x;
    int split_x_nprocs,split_x_rank;
    int split_x_Lrank ,split_x_Rrank;

    MPI_Comm mpi_subcomm_split_y;
    int split_y_nprocs,split_y_rank;
    int split_y_Lrank ,split_y_Rrank;

    MPI_Comm mpi_subcomm_split_z;
    int split_z_nprocs,split_z_rank;
    int split_z_Lrank ,split_z_Rrank;

    MPI_Comm_split(mpi_subcomm_split_xy, split_xy_coord[1], rank, &mpi_subcomm_split_x);
    MPI_Comm_split(mpi_subcomm_split_xy, split_xy_coord[0], rank, &mpi_subcomm_split_y);
    MPI_Comm_split(mpi_subcomm_split_yz, split_yz_coord[0], rank, &mpi_subcomm_split_z);

    MPI_Comm_size(mpi_subcomm_split_x, &split_x_nprocs);
    MPI_Comm_size(mpi_subcomm_split_y, &split_y_nprocs);
    MPI_Comm_size(mpi_subcomm_split_z, &split_z_nprocs);
    MPI_Comm_rank(mpi_subcomm_split_x, &split_x_rank);
    MPI_Comm_rank(mpi_subcomm_split_y, &split_y_rank);
    MPI_Comm_rank(mpi_subcomm_split_z, &split_z_rank);

    MPI_Cart_create(mpi_subcomm_split_x, 1, &split_x_nprocs, &periods[0], 0, &mpi_subcomm_split_x);
    MPI_Cart_create(mpi_subcomm_split_y, 1, &split_y_nprocs, &periods[1], 0, &mpi_subcomm_split_y);
    MPI_Cart_create(mpi_subcomm_split_z, 1, &split_z_nprocs, &periods[2], 0, &mpi_subcomm_split_z);
    
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