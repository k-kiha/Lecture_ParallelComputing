#include <stdio.h>
#include <omp.h>
#include <unistd.h>

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
                sleep(1);
            }

            #pragma omp section
            {
                printf("%d: am\n", omp_get_thread_num());
                sleep(1);
            }

            #pragma omp section
            {
                printf("%d: Open\n", omp_get_thread_num());
                sleep(1);
            }

            #pragma omp section
            {
                printf("%d: MP\n", omp_get_thread_num());
                sleep(1);
            }


            #pragma omp section
            {
                printf("%d: xxx\n", omp_get_thread_num());
                sleep(1);
            }

            #pragma omp section
            {
                printf("%d: ooo\n", omp_get_thread_num());
                sleep(1);
            }
        }

    }


    return 0;
}


