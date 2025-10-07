#include <stdio.h>
#include <omp.h>

int main(void) {
    omp_set_num_threads(4);
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        #pragma omp master
        {
            printf("AAA: %d\n", tid);
        }

        printf("BBB: %d\n", tid);
    }


    return 0;
}


