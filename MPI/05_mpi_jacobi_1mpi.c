#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include"./tools/mpi_tools.h"

int main(int argc, char *argv[]) {

    int n=128;
    int n1=n;
    int n2=n;
    double *xbuf;
    double *x;
    double *b;
    double h=1.0/(n+1);
    double r1,r2,buf;
    
    //--- MPI
    int rank, nprocs;
    MPI_INFO mpi_info;
    int n1sub, n2sub;
    double r1sub,r2sub;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    mpi_tool_setmap_2D(rank, nprocs, &mpi_info);
    mpi_tool_setsubdomain_2D(rank, n1, n2, &mpi_info);
    
    n1sub = mpi_info.local_n1;
    n2sub = mpi_info.local_n2;

    xbuf = (double *)malloc((n1sub+2)*(n2sub+2) * sizeof(double));
    x = (double *)malloc((n1sub+2)*(n2sub+2) * sizeof(double));
    b = (double *)malloc((n1sub+2)*(n2sub+2) * sizeof(double));

    for(int i=0; i<n1sub+2; i++) {
        for(int j=0; j<n2sub+2; j++) {
            b[i*(n2sub+2) + j] = 0.0;
            x[i*(n2sub+2) + j] = 0.0;
            xbuf[i*(n2sub+2) + j] = x[i*(n2sub+2) + j];

            if (i+mpi_info.ista-1==n1/2 && j+mpi_info.jsta-1==n2/2) {
                b[i*(n2sub+2) + j] = 1.0/(h*h);
            }
        }
    }

    for(int iter=0; iter<1000000; iter++) {
        for(int i=1; i<=n1sub; i++) {
            for(int j=1; j<=n2sub; j++) {
                x[i*(n2sub+2) + j] = -h*h/4.0 * b[i*(n2sub+2) + j]
                                  +1.0/4.0 * ( xbuf[(i-1)*(n2sub+2) + (j  )]
                                             + xbuf[(i+1)*(n2sub+2) + (j  )]
                                             + xbuf[(i  )*(n2sub+2) + (j-1)]
                                             + xbuf[(i  )*(n2sub+2) + (j+1)] );
            }
        }

        mpi_tool_communicate_w2e(mpi_info, x, n1sub, n2sub);
        mpi_tool_communicate_e2w(mpi_info, x, n1sub, n2sub);
        mpi_tool_communicate_s2n(mpi_info, x, n1sub, n2sub);
        mpi_tool_communicate_n2s(mpi_info, x, n1sub, n2sub);

        r1 = 0.;
        r2 = 0.;
        r1sub = 0.;
        r2sub = 0.;
        for(int i=1; i<=n1sub; i++) {
            for(int j=1; j<=n2sub; j++) {
                buf = b[(i  )*(n2sub+2) + (j  )]
                      - (  x[(i-1)*(n2sub+2) + (j  )] + x[(i+1)*(n2sub+2) + (j  )]
                         + x[(i  )*(n2sub+2) + (j-1)] + x[(i  )*(n2sub+2) + (j+1)]
                         - x[(i  )*(n2sub+2) + (j  )]*4.0
                        )/(h*h);
                r1sub = r1sub + buf*buf;
                r2sub = r2sub + b[(i  )*(n2sub+2) + (j  )]*b[(i  )*(n2sub+2) + (j  )];
            }
        }
        mpi_tool_communicate_allreduce(mpi_info, r1sub, &r1);
        mpi_tool_communicate_allreduce(mpi_info, r2sub, &r2);

        if(rank==0 && iter%1000==0)printf("iter=%8d, r1/r2=%12.3e\n", iter, sqrt(r1/r2));
        if (sqrt(r1/r2)<1e-12)
        {
            if (rank==0) printf("!!!CONVERGED!!!\n");
            if (rank==0) printf("iter=%8d, r1/r2=%12.3e\n", iter, sqrt(r1/r2));

            break;
        }
        for(int i=0; i<n1sub+2; i++) {
            for(int j=0; j<n2sub+2; j++) {
                xbuf[i*(n2sub+2) + j] = x[i*(n2sub+2) + j];
            }
        }

    }  

    free(x);
    free(xbuf);
    free(b);

    mpi_tool_clean(&mpi_info);
    MPI_Finalize();

    return 0;
}


    // for (int i = 0; i < nprocs; i++)
    // {
    //     if (rank==i){
    //         printf("=== Rank %d ===\n", rank);
    //         for(int jj=0; jj<n2sub+2; jj++) {
    //             for(int ii=0; ii<n1sub+2; ii++) {
    //                 printf("%8.3f ", x[ii*(n2sub+2) + jj]);
    //             }
    //             printf("\n");
    //         }
    //         // printf(": %f\n", r1);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
