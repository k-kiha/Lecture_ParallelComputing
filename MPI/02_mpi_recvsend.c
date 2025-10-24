#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    int A,B;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    A = 100 + rank;
    B = -1;

    printf("Before: Rank %d has A = %d, B = %d\n", rank, A, B);

    if (rank == 0) MPI_Send(&A, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD);
    if (rank == 1) MPI_Recv(&B, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);

    printf("After: Rank %d has A = %d, B = %d\n", rank, A, B);

    MPI_Finalize();

    return 0;
}