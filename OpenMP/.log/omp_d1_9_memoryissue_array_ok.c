#include <stdio.h>
#include <omp.h>

int main(void) {
    int tid;
    int a[10];

    omp_set_num_threads(4);

    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        a[tid] = tid;
    }

    for(int i=0; i<4; i++) 
    {
        printf("a[%d] = %d\n", i, a[i]);
    }

    return 0;
}
