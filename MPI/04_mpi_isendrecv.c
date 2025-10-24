#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    int A,B;
    MPI_Request req1, req2;
    MPI_Status status1, status2;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    A = 100 + rank;
    B = -1;

    printf("Before: Rank %d has A = %d, B = %d\n", rank, A, B);

    if (rank == 0){ 
        MPI_Irecv(&B, 1, MPI_INTEGER, 1, 88, MPI_COMM_WORLD, &req2);
        MPI_Isend(&A, 1, MPI_INTEGER, 1, 55, MPI_COMM_WORLD, &req1);
    }
    if (rank == 1){ 
        MPI_Irecv(&B, 1, MPI_INTEGER, 0, 55, MPI_COMM_WORLD, &req2);
        MPI_Isend(&A, 1, MPI_INTEGER, 0, 88, MPI_COMM_WORLD, &req1);
    }

    MPI_Wait(&req1, &status1);
    MPI_Wait(&req2, &status2);

    printf("After: Rank %d has A = %d, B = %d\n", rank, A, B);

    MPI_Finalize();

    return 0;
}