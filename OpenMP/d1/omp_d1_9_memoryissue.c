#include <stdio.h>
#include <omp.h>

int main(void) {
    int tid;

    omp_set_num_threads(4);

    #pragma omp parallel
    {
        tid = omp_get_thread_num();
        sleep(1); 
        printf("I am %d, tid = %d\n", omp_get_thread_num(), tid);
    }
    return 0;
}



