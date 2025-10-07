#include <stdio.h>
#include <omp.h>

int main(void) {
    omp_set_num_threads(4);

    int a = 100;

    
    #pragma omp parallel
    {
        
        #pragma omp sections 
        {
            #pragma omp section 
            {
                printf("%d: I\n", omp_get_thread_num());
            }

            #pragma omp section
            {
                printf("%d: am\n", omp_get_thread_num());
            }

            #pragma omp section
            {
                printf("%d: Open\n", omp_get_thread_num());
            }

            #pragma omp section
            {
                printf("%d: MP\n", omp_get_thread_num());
            }
        }

    }

    return 0;
}




