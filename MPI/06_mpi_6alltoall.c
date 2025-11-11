#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    int A[10],B[10];
    MPI_Status status;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    for (int i = 0; i < 10; i++) {
        A[i] = 100*(i+1) + rank;
        B[i] = -1;
    }

    for (int i = 0; i < nprocs; i++)
    {
        if (rank==i) {
            printf("rank=%d A:", rank);
            for (int j = 0; j < 10; j++) {
                printf(" %5d", A[j]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /**---**/
    MPI_Alltoall(A, 2, MPI_INT, B, 2, MPI_INT, MPI_COMM_WORLD);
    /**---**/
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < nprocs; i++)
    {
        if (rank==i) {
            printf("rank=%d B:", rank);
            for (int j = 0; j < 10; j++) {
                printf(" %5d", B[j]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}