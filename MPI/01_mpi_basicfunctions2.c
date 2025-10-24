#include<stdio.h>
#include<unistd.h>
#include<mpi.h>

int main(int argc, char *argv[]) {
    int rank, nprocs;
    double t_a, t_b;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    t_a = MPI_Wtime();
    sleep(1); 
    t_b = MPI_Wtime();

    for (int i = 0; i < nprocs; i++)
    {
        if (i == rank) printf("Process %d elapsed time: %f seconds\n", rank, t_b - t_a);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    

    MPI_Finalize();

    return 0;
}