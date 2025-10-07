#include <stdio.h>
// <<< !!

#define N 20

int main(void) {
    int A[N];

    #pragma omp parallel
    {
        int thread_id,num_threads;        
        int chunk_size;
        int start,end;

        thread_id   = 0; // <<< !!
        num_threads = 0; // <<< !!
        chunk_size  = 0; // <<< !!
        start       = 0; // <<< !!
        end         = 0; // <<< !!
        
        for (int i = 0; i < 0; i++) { // <<< !!
            A[i] = 100 + i;
        }
    }

    for (int i = 0; i < 20; i++) {
        printf("%8d %0.4d\n", i,A[i]);
    }
    
    return 0;
}




