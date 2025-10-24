#include<stdio.h>
#include<unistd.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    double Aarray[512*512],Barray[512*512];
    MPI_Status status;

    double ta,tb;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    ta = MPI_Wtime();
    if (rank == 0){ 
        MPI_Send(&Aarray, 512*512, MPI_DOUBLE, 1, 55, MPI_COMM_WORLD);
        MPI_Recv(&Barray, 512*512, MPI_DOUBLE, 1, 88, MPI_COMM_WORLD, &status);
    }
    if (rank == 1){ 
        MPI_Recv(&Barray, 512*512, MPI_DOUBLE, 0, 55, MPI_COMM_WORLD, &status);
        MPI_Send(&Aarray, 512*512, MPI_DOUBLE, 0, 88, MPI_COMM_WORLD);
    }

    // Sleep for 1 second
    usleep(10000);

    tb = MPI_Wtime();
    
    printf("rank %d: Elapsed time: %f sec\n", rank, tb - ta);

    MPI_Finalize();

    return 0;
}