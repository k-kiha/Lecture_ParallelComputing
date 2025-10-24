#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[]) {

    int n=128;
    int n1=n;
    int n2=n;
    double *xbuf;
    double *x;
    double *b;
    double h=1.0/(n+1);
    double r1,r2,buf;

    xbuf    = (double *)malloc((n1+2)*(n2+2) * sizeof(double));
    x       = (double *)malloc((n1+2)*(n2+2) * sizeof(double));
    b       = (double *)malloc((n1+2)*(n2+2) * sizeof(double));

    for(int i=0; i<n1+2; i++) {
        for(int j=0; j<n2+2; j++) {
            b[i*(n2+2) + j] = 0.0;
            x[i*(n2+2) + j] = 0.0;
            xbuf[i*(n2+2) + j] = x[i*(n2+2) + j];
        }
    }
    b[(n1/2)*(n2+2) + (n2/2)] = 1.0/(h*h);

    for(int iter=0; iter<1000000; iter++) {
        for(int i=1; i<=n1; i++) {
            for(int j=1; j<=n2; j++) {
                x[i*(n2+2) + j] = -h*h/4.0 * b[i*(n2+2) + j]
                                  +1.0/4.0 * ( xbuf[(i-1)*(n2+2) + (j  )]
                                             + xbuf[(i+1)*(n2+2) + (j  )]
                                             + xbuf[(i  )*(n2+2) + (j-1)]
                                             + xbuf[(i  )*(n2+2) + (j+1)] );
            }
        }


        for(int i=1; i<=n1; i++) {
            x[i*(n2+2) + 0   ] = x[i*(n2+2) + n2];
            x[i*(n2+2) + n2+1] = x[i*(n2+2) + 1 ];
        }

        r1 = 0.;
        r2 = 0.;
        for(int i=1; i<=n1; i++) {
            for(int j=1; j<=n2; j++) {
                buf = b[(i  )*(n2+2) + (j  )]
                      - (  x[(i-1)*(n2+2) + (j  )] + x[(i+1)*(n2+2) + (j  )] 
                         + x[(i  )*(n2+2) + (j-1)] + x[(i  )*(n2+2) + (j+1)] 
                         - x[(i  )*(n2+2) + (j  )]*4.0
                        )/(h*h);
                r1 = r1 + buf*buf;
                r2 = r2 + b[(i  )*(n2+2) + (j  )]*b[(i  )*(n2+2) + (j  )];
            }
        }

        if(iter%1000==0)printf("iter=%8d, r1/r2=%12.3e\n", iter, sqrt(r1/r2));
        if (sqrt(r1/r2)<1e-12)
        {
            printf("!!!CONVERGED!!!\n");
            printf("iter=%8d, r1/r2=%12.3e\n", iter, sqrt(r1/r2));
            break;
        }
        for(int i=0; i<n1+2; i++) {
            for(int j=0; j<n2+2; j++) {
                xbuf[i*(n2+2) + j] = x[i*(n2+2) + j];
            }
        }
        
    }


    free(xbuf);
    free(x);
    free(b);
    return 0;
}