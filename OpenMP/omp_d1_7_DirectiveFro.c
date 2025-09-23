#include <stdio.h>
#include <omp.h>

#define N 20

int main(void) {
    int A[N];

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N; i++) {
            A[i] = 100 + i;
        }
    }

    for (int i = 0; i < 20; i++) {
        printf("%8d %0.4d\n", i,A[i]);
    }
    
    return 0;
}
