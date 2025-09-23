#include <stdio.h>
#include <omp.h>

int g = 100;


int main(void) {
    omp_set_num_threads(4);

    #pragma omp parallel 
    {
        int tid = omp_get_thread_num();
        g = g*10 + tid;
        printf("R1: tid=%d, g=%d\n", tid, g);
    }
    printf("== g=%d\n", g);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        printf("R2: tid=%d, g=%d\n", tid, g);
    }
    printf("== g=%d\n", g);
    
    return 0;
}
