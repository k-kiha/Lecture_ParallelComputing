#include <stdio.h>
#include <omp.h>
#define N 1024

int main(void) {
    int a[N][N];
    double ta_serial,tb_serial;
    double ta,tb;
    double ta_omp_dy,tb_omp_dy;
    omp_set_num_threads(2);

    ta_serial=omp_get_wtime();
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            a[i][j] = i + j;
        }
    }
    tb_serial=omp_get_wtime();

    ta=omp_get_wtime();
    #pragma omp parallel 
    {
        #pragma omp for
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                a[i][j] = i + j;
            }
        }
    }
    tb=omp_get_wtime();

    ta_omp_dy=omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1) 
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                a[i][j] = i + j;
            }
        }
    }
    tb_omp_dy=omp_get_wtime();

    printf("serial.    : %f sec\n", tb_serial - ta_serial);
    printf("omp.       : %f sec\n", tb - ta);
    printf("omp dynamic: %f sec\n", tb_omp_dy - ta_omp_dy);

    return 0;
}






