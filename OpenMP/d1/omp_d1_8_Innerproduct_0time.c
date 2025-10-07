#include <stdio.h>
#include <omp.h>

#define N 1024000

int main(void) {
    int A[N],B[N];
    int sum0_ser=0;
    int sum1_cri=0;
    int sum2_ato=0;
    int sum3_arr=0;
    int sum3_2arr=0;
    int sum4_red=0;

    double t0_a,t0_b;
    double t1_a,t1_b;
    double t2_a,t2_b;
    double t3_a,t3_b;
    double t3_2a,t3_2b;
    double t4_a,t4_b;
    
    omp_set_num_threads(4);

    for (int i = 0; i < N; i++) {
        A[i] = i;
        B[i] = i;
    }

    t0_a = omp_get_wtime();
    for (int i = 0; i < N; i++) {
        sum0_ser+= A[i]*B[i];
    }
    t0_b = omp_get_wtime();
    
    t1_a = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        #pragma omp critical
        sum1_cri += A[i]*B[i];
    }
    t1_b = omp_get_wtime();

    t2_a = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        #pragma omp atomic
        sum2_ato += A[i]*B[i];
    }
    t2_b = omp_get_wtime();
        

    int tid;
    int tmp_sum[4] = {0,0,0,0};
    t3_a = omp_get_wtime();
    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < N; i++) {
            tmp_sum[tid] += A[i]*B[i];
        }
    }
    for (int i = 0; i < 4; i++) {
        sum3_arr += tmp_sum[i];
    }
    t3_b = omp_get_wtime();

    int tmp_sum2[128] = {0,};
    t3_2a = omp_get_wtime();
    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < N; i++) {
            tmp_sum2[tid*32] += A[i]*B[i];
        }
    }
    for (int i = 0; i < 4; i++) {
        sum3_2arr += tmp_sum2[i*32];
    }
    t3_2b = omp_get_wtime();

    t4_a = omp_get_wtime();
    #pragma omp parallel for reduction(+:sum4_red)
    for (int i = 0; i < N; i++) {
        sum4_red += A[i]*B[i];
    }
    t4_b = omp_get_wtime();
    
    printf("Serial   : %d in %f sec.\n", sum0_ser,(t0_b - t0_a));
    printf("Critical : %d in %f sec. (%f sec.)\n", sum1_cri,(t1_b - t1_a),(t0_b - t0_a)/4.);
    printf("Atomic   : %d in %f sec. (%f sec.)\n", sum2_ato,(t2_b - t2_a),(t0_b - t0_a)/4.);
    printf("Array1   : %d in %f sec. (%f sec.)\n", sum3_arr,(t3_b - t3_a),(t0_b - t0_a)/4.);
    printf("Array2   : %d in %f sec. (%f sec.)\n", sum3_2arr,(t3_2b - t3_2a),(t0_b - t0_a)/4.);
    printf("Reduction: %d in %f sec. (%f sec.)\n", sum4_red,(t4_b - t4_a),(t0_b - t0_a)/4.);

    return 0;
}
