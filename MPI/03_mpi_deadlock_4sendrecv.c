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

    if (rank == 0){ 
        MPI_Sendrecv(&A, 1, MPI_INTEGER, 1, 55, &B, 1, MPI_INTEGER, 1, 88, MPI_COMM_WORLD, &status);
    }
    if (rank == 1){ 
        MPI_Sendrecv(&A, 1, MPI_INTEGER, 0, 88, &B, 1, MPI_INTEGER, 0, 55, MPI_COMM_WORLD, &status);
    }

    printf("After: Rank %d has A = %d, B = %d\n", rank, A, B);
    
    //---------------------------------------------
    int Aarray[512*512],Barray[512*512];

    if (rank == 0){ 
        MPI_Sendrecv(&Aarray, 512*512, MPI_INTEGER, 1, 55, &Barray, 512*512, MPI_INTEGER, 1, 88, MPI_COMM_WORLD, &status);
        printf("OK!! :: rank:%d\n", rank);
    }
    if (rank == 1){ 
        MPI_Sendrecv(&Aarray, 512*512, MPI_INTEGER, 0, 88, &Barray, 512*512, MPI_INTEGER, 0, 55, MPI_COMM_WORLD, &status);
        printf("OK!! :: rank:%d\n", rank);
    }

    MPI_Finalize();

    return 0;
}