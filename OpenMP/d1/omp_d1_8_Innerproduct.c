#include <stdio.h>
#include <omp.h>

#define N 20

int main(void) {
    int A[N],B[N];
    int sum=0;

    omp_set_num_threads(4);

    for (int i = 0; i < N; i++) {
        A[i] = 100 + i;
        B[i] = 200 + i;
    }
    
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N; i++) {
            sum += A[i]*B[i];
        }
    }
    printf("Parallel: %d vs (459470)\n", sum);

    return 0;
}







