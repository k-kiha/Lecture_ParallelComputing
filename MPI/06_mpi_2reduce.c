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

    for (int i = 0; i < nprocs; i++)
    {
        if (rank==i) printf("Before: Rank %d has A = %d, B = %d\n", rank, A, B);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /**---**/
    MPI_Reduce(&A, &B, 1, MPI_INT, MPI_SUM, 3, MPI_COMM_WORLD);
    /**---**/

    for (int i = 0; i < nprocs; i++)
    {
        if (rank==i) printf("After: Rank %d has A = %d, B = %d\n", rank, A, B);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}