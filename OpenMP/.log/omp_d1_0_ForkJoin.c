#include <stdio.h>

int main(void) {
    printf("AAA: Hello, serial!\n");

    #pragma omp parallel
    {
        printf("BBB: Hello, OpenMP!\n");
    }

    printf("CCC: Hello, serial!\n");
    return 0;
}