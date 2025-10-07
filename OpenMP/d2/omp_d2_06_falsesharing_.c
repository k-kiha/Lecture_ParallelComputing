#include <stdio.h>
#include <omp.h>

#define N 1024000

int main(void) {
    int A[N],B[N];
    int sum0_ser=0;
    int sum3_arr=0;
    int sum3_2arr=0;

    double t0_a,t0_b;
    double t1_a,t1_b;
    double t2_a,t2_b;
    
    omp_set_num_threads(2);

    for (int i = 0; i < N; i++) {
        A[i] = i;
        B[i] = i;
    }

    t0_a = omp_get_wtime();
    for (int i = 0; i < N; i++) {
        sum0_ser+= A[i]*B[i];
    }
    t0_b = omp_get_wtime();

    int tid;
    int tmp_sum[2] = {0,0};
    t1_a = omp_get_wtime();
    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < N; i++) {
            tmp_sum[tid] += A[i]*B[i];
        }
    }
    for (int i = 0; i < 2; i++) {
        sum3_arr += tmp_sum[i];
    }
    t1_b = omp_get_wtime();



    int tmp_sum2[128] = {0,};
    t2_a = omp_get_wtime();
    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < N; i++) {
            tmp_sum2[tid*32] += A[i]*B[i];
        }
    }
    for (int i = 0; i < 2; i++) {
        sum3_2arr += tmp_sum2[i*32];
    }
    t2_b = omp_get_wtime();
    


    printf("Serial   : %d in %f sec.\n", sum0_ser,(t0_b - t0_a));
    printf("Array1   : %d in %f sec. (%f sec.)\n", sum3_arr,(t1_b - t1_a),(t0_b - t0_a)/4.);
    printf("Array2   : %d in %f sec. (%f sec.)\n", sum3_2arr,(t2_b - t2_a),(t0_b - t0_a)/4.);



    return 0;
}
