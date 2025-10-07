#include <stdio.h>
#include <omp.h>

int main() {
    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel num_threads(4)
    {
        omp_set_lock(&lock);   // 한 번에 한 스레드만 통과
        printf("Thread %d in critical section\n", omp_get_thread_num());
        omp_unset_lock(&lock);
    }

    omp_destroy_lock(&lock);
    return 0;
}
