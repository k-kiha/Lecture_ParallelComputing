#include<stdio.h>

int main(int argc, char *argv[]) {

    double x,y;

    for(int i=0; i<4; i++) {
        x = (double)i;
        y = x*x + x + 1;
        printf("i=%d, y=x^2+x+1=%5.2f\n", i, y);
    }

    return 0;
}