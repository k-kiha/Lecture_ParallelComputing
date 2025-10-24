#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    int BUF,SUM;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    SUM = 0;
    printf("Before: Rank %d has SUM = %d\n", rank, SUM);

    if (rank == 0) {
        for (int i = 1; i < nprocs; i++)
        {
            MPI_Recv(&BUF, 1, MPI_INTEGER, i, 88, MPI_COMM_WORLD, &status);
            SUM += BUF;
        }
        
    }
    else
    {
        MPI_Send(&rank, 1, MPI_INTEGER, 0, 88, MPI_COMM_WORLD);
    }

    printf("After: Rank %d has SUM = %d\n", rank, SUM);

    MPI_Finalize();

    return 0;
}