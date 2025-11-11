#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    int A[20],B[20];
    int counts_send[4],counts_recv[4];
    int displs_send[4],displs_recv[4];
    MPI_Status status;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    for (int i = 0; i < 20; i++) {
        A[i] = 100*(i+1) + rank;
        B[i] = -1;
    }

    for (int i = 0; i < nprocs; i++) {
        counts_send[i] = i+1;
        counts_recv[i] = rank+1;
    }

    displs_send[0] = 0;
    displs_recv[0] = 0;
    for (int i = 1; i < nprocs; i++) {
        displs_send[i] = displs_send[i-1] + counts_send[i-1];
        displs_recv[i] = displs_recv[i-1] + counts_recv[i-1];
    }
    /**---**/
    MPI_Alltoallv(
        A, counts_send, displs_send, MPI_INT,
        B, counts_recv, displs_recv, MPI_INT,
        MPI_COMM_WORLD
    );
    /**---**/

    for (int i = 0; i < nprocs; i++)
    {
        if (rank==i) {
            printf("rank=%d B:", rank);
            for (int j = 0; j < 20; j++) {
                printf(" %5d", B[j]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}