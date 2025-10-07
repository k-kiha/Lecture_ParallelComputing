#include <stdio.h>
#include <omp.h>

int main(void) {
    printf("AAA: Hello, serial!\n");

    #pragma omp parallel
    {
        printf("BBB0: Hello, OpenMP--Enviroment!\n");
    }

    #pragma omp parallel num_threads(4)
    {
        printf("BBB1: Hello, OpenMP--Directives!\n");
    }

    #pragma omp parallel
    {
        printf("BBB0: Hello, OpenMP--Enviroment!\n");
    }

    omp_set_num_threads(3);
    #pragma omp parallel
    {
        printf("BBB2: Hello, OpenMP--Runtime Library!\n");
    }

    #pragma omp parallel
    {
        printf("BBB0: Hello, OpenMP--Enviroment!\n");
    }

    printf("CCC: Hello, serial!\n");
    return 0;
}



// nvc -mp main.c -o main
// nvc++ -mp main.cpp -o main
// gcc -fopenmp main.c -o main
// icc -qopenmp main.c -o main
// clang -fopenmp main.c -o main
// export OMP_NUM_THREADS=4
// ./main