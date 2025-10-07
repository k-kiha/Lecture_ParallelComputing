#include <stdio.h>
#include <omp.h>

int main(void) {
    omp_set_num_threads(4);

    int a = 100;
    
    #pragma omp parallel private(a) 
    {
        int tid = omp_get_thread_num();

        #pragma omp single copyprivate(a)
        {
            a = 1000+tid;
            printf("AAA: %d, a=%d\n", tid, a);
        }

        printf("BBB: %d, a=%d\n", tid, a);
    }


    return 0;
}


