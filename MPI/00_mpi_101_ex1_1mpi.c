#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    double x,y;
    int rank, nprocs;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    x = (double)rank;
    y = x*x + x + 1;
    printf("i=%d, y=x^2+x+1=%5.2f\n", rank, y);

    MPI_Finalize();

    return 0;
}