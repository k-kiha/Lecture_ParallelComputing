#include <stdio.h>
#include <omp.h>
#include <unistd.h>
#define N 1024

int main(void) {
    int a[N][N];
    double ta,tb;
    double ta_nw,tb_nw;
    omp_set_num_threads(4);
    
    ta=omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                a[i][j] = i + j;
            }
        }

        #pragma omp for
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                a[i][j] = i - j;
            }
        }
    }
    tb=omp_get_wtime();

    ta_nw=omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                a[i][j] = i + j;
            }
        }

        #pragma omp for
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                a[i][j] = i - j;
            }
        }
    }
    tb_nw=omp_get_wtime();

    printf("With barrier   : %f sec\n", tb - ta);
    printf("Without barrier: %f sec\n", tb_nw - ta_nw);

    return 0;
}


