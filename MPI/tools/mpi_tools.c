#include<stdio.h>   
#include<stdlib.h>
#include<math.h>
#include"./mpi_tools.h"

void mpi_tool_setmap_2D(int rank, int nprocs, MPI_INFO *mpi_info){
    mpi_info->rank = rank;
    mpi_info->nprocs = nprocs;

    mpi_info->map_dims[0] = 2;
    mpi_info->map_dims[1] = 4;

    mpi_info->map_coords[0] = rank / mpi_info->map_dims[1];
    mpi_info->map_coords[1] = rank % mpi_info->map_dims[1];

    int x,y;
    int Wx,Wy;
    x = mpi_info->map_coords[0];
    y = mpi_info->map_coords[1];
    Wx = mpi_info->map_dims[0];
    Wy = mpi_info->map_dims[1];

    mpi_info->w = (x-1)*Wy + (y  );
    mpi_info->e = (x+1)*Wy + (y  );
    mpi_info->s = (x  )*Wy + (y-1);
    mpi_info->n = (x  )*Wy + (y+1);
    if(x-1 < 0  ) mpi_info->w = MPI_PROC_NULL;
    if(x+1 >= Wx) mpi_info->e = MPI_PROC_NULL;
    if(y-1 < 0  ) mpi_info->s = (x  )*Wy + Wy-1;
    if(y+1 >= Wy) mpi_info->n = (x  )*Wy + 0 ;
}

void mpi_tool_setsubdomain_2D(int rank, int n1, int n2, MPI_INFO *mpi_info){
    mpi_info->local_n1 = n1/mpi_info->map_dims[0];
    mpi_info->local_n2 = n2/mpi_info->map_dims[1];
    
    mpi_info->ista = mpi_info->map_coords[0]*mpi_info->local_n1 + 1;
    mpi_info->jsta = mpi_info->map_coords[1]*mpi_info->local_n2 + 1;

    mpi_info->iend = mpi_info->ista + mpi_info->local_n1 - 1;
    mpi_info->jend = mpi_info->jsta + mpi_info->local_n2 - 1;

    int maxnsub=((mpi_info->local_n1+2)>(mpi_info->local_n2+2)?(mpi_info->local_n1+2):(mpi_info->local_n2+2));
    mpi_info->buf_send = (double *)malloc(maxnsub * sizeof(double));
    mpi_info->buf_recv = (double *)malloc(maxnsub * sizeof(double));
}

void mpi_tool_monitor(MPI_INFO mpi_info){
    printf("=== Rank %d ===\n", mpi_info.rank);
    printf("::: map_coords   : (%4d, %4d)\n", mpi_info.map_coords[0], mpi_info.map_coords[1]);
    printf("::: (n1sub,n2sub): (%4d, %4d)\n", mpi_info.local_n1, mpi_info.local_n2);
    printf("::: ista - iend  : (%4d, %4d)\n", mpi_info.ista, mpi_info.iend);
    printf("::: jsta - jend  : (%4d, %4d)\n", mpi_info.jsta, mpi_info.jend);
    printf(":::      %2d     \n", mpi_info.n);
    printf(":::   %2d ~~ %2d  \n", mpi_info.w, mpi_info.e);
    printf(":::      %2d     \n", mpi_info.s);
}

void mpi_tool_communicate_w2e(MPI_INFO mpi_info, double *x, int n1sub, int n2sub){
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    for (int j = 0; j <= n2sub+1; j++)
    {
        mpi_info.buf_send[j] = x[n1sub*(n2sub+2) + j];
    }
    

    MPI_Isend(mpi_info.buf_send, n2sub+2, MPI_DOUBLE, mpi_info.e, 11, MPI_COMM_WORLD, &req1);
    MPI_Irecv(mpi_info.buf_recv, n2sub+2, MPI_DOUBLE, mpi_info.w, 11, MPI_COMM_WORLD, &req2);

    MPI_Wait(&req1, &status1);
    MPI_Wait(&req2, &status2);

    if(mpi_info.w>=0){
        for (int j = 0; j <= n2sub+1; j++)
        {
            x[0*(n2sub+2) + j] = mpi_info.buf_recv[j];
        }
    }

}

void mpi_tool_communicate_e2w(MPI_INFO mpi_info, double *x, int n1sub, int n2sub){
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    for (int j = 0; j <= n2sub+1; j++)
    {
        mpi_info.buf_send[j] = x[1*(n2sub+2) + j];
    }
    

    MPI_Isend(mpi_info.buf_send, n2sub+2, MPI_DOUBLE, mpi_info.w, 11, MPI_COMM_WORLD, &req1);
    MPI_Irecv(mpi_info.buf_recv, n2sub+2, MPI_DOUBLE, mpi_info.e, 11, MPI_COMM_WORLD, &req2);

    MPI_Wait(&req1, &status1);
    MPI_Wait(&req2, &status2);

    if(mpi_info.e>=0)
    {
        for (int j = 0; j <= n2sub+1; j++)
        {
            x[(n1sub+1)*(n2sub+2) + j] = mpi_info.buf_recv[j];
        }
    }

}

void mpi_tool_communicate_s2n(MPI_INFO mpi_info, double *x, int n1sub, int n2sub){
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    for (int i = 0; i <= n1sub+1; i++)
    {
        mpi_info.buf_send[i] = x[i*(n2sub+2) + n2sub];
    }
    

    MPI_Isend(mpi_info.buf_send, n1sub+2, MPI_DOUBLE, mpi_info.n, 11, MPI_COMM_WORLD, &req1);
    MPI_Irecv(mpi_info.buf_recv, n1sub+2, MPI_DOUBLE, mpi_info.s, 11, MPI_COMM_WORLD, &req2);

    MPI_Wait(&req1, &status1);
    MPI_Wait(&req2, &status2);

    if(mpi_info.s>=0){
        for (int i = 0; i <= n1sub+1; i++)
        {
            x[i*(n2sub+2) + 0] = mpi_info.buf_recv[i];
        }
    }

}

void mpi_tool_communicate_n2s(MPI_INFO mpi_info, double *x, int n1sub, int n2sub){
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    for (int i = 0; i <= n1sub+1; i++)
    {
        mpi_info.buf_send[i] = x[i*(n2sub+2) + 1];
    }
    

    MPI_Isend(mpi_info.buf_send, n1sub+2, MPI_DOUBLE, mpi_info.s, 11, MPI_COMM_WORLD, &req1);
    MPI_Irecv(mpi_info.buf_recv, n1sub+2, MPI_DOUBLE, mpi_info.n, 11, MPI_COMM_WORLD, &req2);

    MPI_Wait(&req1, &status1);
    MPI_Wait(&req2, &status2);

    if(mpi_info.n>=0){
        for (int i = 0; i <= n1sub+1; i++)
        {
            x[i*(n2sub+2) + n2sub+1] = mpi_info.buf_recv[i];
        }
    }
}

void mpi_tool_clean(MPI_INFO *mpi_info){
    free(mpi_info->buf_send);
    free(mpi_info->buf_recv);
}   

void mpi_tool_communicate_allreduce(MPI_INFO mpi_info, double local_value, double *global_value){
    double BUF;
    MPI_Status status;
    *global_value = local_value;

    // Gather all local_value to rank 0 and sum them up
    if (mpi_info.rank == 0) {
        for (int i = 1; i < mpi_info.nprocs; i++){
            MPI_Recv(&BUF, 1, MPI_DOUBLE, i, 88, MPI_COMM_WORLD, &status);
            *global_value += BUF;
        }
        
    }
    else{
        MPI_Send(global_value, 1, MPI_DOUBLE, 0, 88, MPI_COMM_WORLD);
    }  

    // Broadcast the global_value from rank 0 to all ranks
    if (mpi_info.rank == 0) {
        for (int i = 1; i < mpi_info.nprocs; i++){
            MPI_Send(global_value, 1, MPI_DOUBLE, i, 88, MPI_COMM_WORLD);
        }
        
    }
    else{
        MPI_Recv(global_value, 1, MPI_DOUBLE, 0, 88, MPI_COMM_WORLD, &status);
    }  

}




void mpi_tool_cart_2D(int rank, int nprocs, MPI_INFO *mpi_info){
    mpi_info->rank = rank;
    mpi_info->nprocs = nprocs;

    mpi_info->map_dims[0] = 2;
    mpi_info->map_dims[1] = 4;
    mpi_info->map_periods[0] = 0;
    mpi_info->map_periods[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, mpi_info->map_dims, mpi_info->map_periods, 0, &mpi_info->cart_comm); 

    MPI_Cart_coords(mpi_info->cart_comm, rank, 2,               mpi_info->map_coords);
    MPI_Cart_rank  (mpi_info->cart_comm, mpi_info->map_coords, &mpi_info->rank);

    MPI_Cart_shift(mpi_info->cart_comm, 0, 1, &mpi_info->w, &mpi_info->e);
    MPI_Cart_shift(mpi_info->cart_comm, 1, 1, &mpi_info->s, &mpi_info->n);
   
}

