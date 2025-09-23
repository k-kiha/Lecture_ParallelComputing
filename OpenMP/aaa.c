#include <stdio.h>
#include <omp.h>

int main(void) {
    int a,b,c,d;
    
    a =10;
    b =5;
    c =2;
    d =4;
    
    
    #pragma omp parallel private(c,d) num_threads(10)
    {
        c = a+10;
        d = b+5;
        a = a+1;
    }

    printf("a = %d, b = %d, c = %d, d = %d\n", a, b, c, d);

    return 0;
}
