#include <stdio.h>
#include <omp.h>
#include <unistd.h>

int main(void) {
    int flag = 0;
    omp_set_num_threads(2);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        if (tid == 0) {
            sleep(1);  

            flag = 1;  

        }
            
        if (tid == 1) {            
            while (!flag) {
                
            }
            printf("Thread 1: detected flag = %d, proceeding!\n", flag);
        }
    }

    return 0;
}


