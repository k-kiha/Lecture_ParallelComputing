#include <stdio.h>
#include <omp.h>

#define N 20

int main(void) {
    int A[N];

    #pragma omp parallel
    {
        int thread_id,num_threads;        
        int chunk_size;
        int start,end;

        thread_id   = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        chunk_size  = N / num_threads;
        start       = thread_id * chunk_size;
        end         = (thread_id + 1) * chunk_size;
        
        for (int i = start; i < end; i++) {
            A[i] = 100 + i;
        }
    }

    for (int i = 0; i < 20; i++) {
        printf("%8d %0.4d\n", i,A[i]);
    }
    
    return 0;
}


