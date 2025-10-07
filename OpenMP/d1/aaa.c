// build: gcc -O2 -fopenmp demo.c -o demo
#include <stdio.h>
#include <omp.h>
#include <math.h>

int main(void) {
    int A;

    omp_set_num_threads(4);
    
    A=9;
    #pragma omp parallel for private(A)
    for (int i = 0; i < 8; ++i) {
        A += pow(10,i+1);
        printf("private(A) Thread %d: A=%12d @%2d\n", omp_get_thread_num(),A, i);
    }
    printf("After private(A).  : A=%12d\n\n", A);

    A=9;
    #pragma omp parallel for firstprivate(A)
    for (int i = 0; i < 8; ++i) {
        A += pow(10,i+1);
        printf("firstprivate(A) Thread %d: A=%12d @%2d\n", omp_get_thread_num(),A, i);
    }
    printf("After firstprivate(A).  : A=%12d\n\n", A);

    A=9;
    #pragma omp parallel for lastprivate(A)
    for (int i = 0; i < 8; ++i) {
        A += pow(10,i+1);
        printf("lastprivate(A) Thread %d: A=%12d @%2d\n", omp_get_thread_num(),A, i);
    }
    printf("After lastprivate(A).  : A=%12d\n\n", A);

    return 0;
}
