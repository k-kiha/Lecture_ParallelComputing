#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    int A[10]={0,},B[10]={0,};
    MPI_Status status;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    for (int i = 0; i < 10; i++) {
        A[i] = 100*(i+1) + rank;
        B[i] = -1;
    }

    for (int rrr = 0; rrr < nprocs; rrr++)
    {
        if (rank==rrr)
        for (int i = 0; i < 5; i++) {
            printf("Before: Rank %d has A[%1d] = %4d, B[%1d] = %4d\n", rank, i, A[i], i, B[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    

    if (rank == 0) MPI_Send(&A[0], 5, MPI_INT, 1, 0, MPI_COMM_WORLD);
    if (rank == 1) MPI_Recv(&B[3], 5, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    for (int rrr = 0; rrr < nprocs; rrr++)
    {
        if (rank==rrr)
        for (int i = 0; i < 10; i++) {
            printf("After: Rank %d has A[%1d] = %4d, B[%1d] = %4d\n", rank, i, A[i], i, B[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}