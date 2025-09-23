#include <stdio.h>
#include <omp.h>

int main(void) {

    omp_set_num_threads(3);
    #pragma omp parallel
    {
        printf("AAA: my thread=%d, total=%d!\n", omp_get_thread_num(), omp_get_num_threads());
    }

    omp_set_num_threads(4);
    #pragma omp parallel
    {
        printf("BBB: my thread=%d, total=%d!\n", omp_get_thread_num(), omp_get_num_threads());
    }

    return 0;
}