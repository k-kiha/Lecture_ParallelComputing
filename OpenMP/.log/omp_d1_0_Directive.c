#include <stdio.h>

int main(void) {
    printf("AAA: Hello, serial!\n");

    #pragma omp parallel num_threads(2)
    {
        printf("BBB: Hello, OpenMP!\n");
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