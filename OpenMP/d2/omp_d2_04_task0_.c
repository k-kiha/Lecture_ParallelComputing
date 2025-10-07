#include <stdio.h>
#include <omp.h>
#include <unistd.h>

int main(void) {
    double ta,tb;
    omp_set_num_threads(4);

    ta = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        {
            printf("Thread %d: creating tasks...\n", omp_get_thread_num());

            for (int i = 0; i < 16; i++) {
                #pragma omp task
                {
                    // sleep(1);
                    printf("  Task %d executed by thread %d\n", i, omp_get_thread_num());
                }
            }

            printf("Thread %d: all tasks created.\n", omp_get_thread_num());
        } 
    }

    tb = omp_get_wtime();
    printf("time=%f\n", tb-ta);

    return 0;
}


