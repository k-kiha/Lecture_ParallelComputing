#include<stdio.h>
#include<unistd.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    int rank, nprocs;
    double Aarray[512*512],Barray[512*512];
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    double ta,tb;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    ta = MPI_Wtime();
    if (rank == 0){ 
        MPI_Irecv(&Barray, 512*512, MPI_DOUBLE, 1, 88, MPI_COMM_WORLD, &req2);
        MPI_Isend(&Aarray, 512*512, MPI_DOUBLE, 1, 55, MPI_COMM_WORLD, &req1);
    }
    if (rank == 1){ 
        MPI_Irecv(&Barray, 512*512, MPI_DOUBLE, 0, 55, MPI_COMM_WORLD, &req2);
        MPI_Isend(&Aarray, 512*512, MPI_DOUBLE, 0, 88, MPI_COMM_WORLD, &req1);
    }

    // Sleep for 1 second
    usleep(10000);

    MPI_Wait(&req1, &status1);
    MPI_Wait(&req2, &status2);
    tb = MPI_Wtime();
    
    printf("rank %d: Elapsed time: %f sec\n", rank, tb - ta);

    MPI_Finalize();

    return 0;
}