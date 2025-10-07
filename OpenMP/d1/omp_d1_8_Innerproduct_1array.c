#include <stdio.h>
#include <omp.h>

#define N 20

int main(void) {
    int A[N],B[N];
    int sum=0;
    int tid;

    omp_set_num_threads(4);

    for (int i = 0; i < N; i++) {
        A[i] = 100 + i;
        B[i] = 200 + i;
    }

    int tmp_sum[4] = {0,0,0,0};
    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();

        #pragma omp for
        for (int i = 0; i < N; i++) {
            tmp_sum[tid] += A[i]*B[i];
        }
    }
    for (int i = 0; i < 4; i++) {
        sum += tmp_sum[i];
    }
    printf("Parallel: %d vs (459470)\n", sum);

    return 0;
}



