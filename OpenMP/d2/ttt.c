#include <stdio.h>
#include <omp.h>
#include <unistd.h>

int main(void) {
    double ta,tb;
    omp_set_num_threads(4);
    
    ta = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp for 
        for (int i = 0; i < 4; i++) {
            sleep(1);
            
            #pragma omp critical
            printf("ID=%d\n", omp_get_thread_num());
        }
    }
    tb = omp_get_wtime();
    printf("time=%f\n", tb-ta);

    return 0;
}




