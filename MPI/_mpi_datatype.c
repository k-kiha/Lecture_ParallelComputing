#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[]) {

    MPI_Datatype mytype;
    
    MPI_Init(&argc, &argv);

MPI_Datatype mytype_send;

int ndims = 2;
int array_of_sizes[2]    = {8+2,8+2};
int array_of_subsizes[2] = {8+2,  1};
int array_of_starts[2]   = {  0,  4};

MPI_Type_create_subarray(ndims
                        array_of_sizes, 
                        array_of_subsizes, 
                        array_of_starts,
                        MPI_ORDER_C, MPI_DOUBLE, &mytype_send);

    MPI_Type_commit(&mytype);
    
    MPI_Type_free(&mytype);


    MPI_Finalize();

    return 0;
}

