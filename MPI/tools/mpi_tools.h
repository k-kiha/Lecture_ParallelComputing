#include<mpi.h>

typedef struct{
    int rank;
    int nprocs;

    int map_dims[2];
    int map_coords[2];
    int w,e,n,s;

    int local_n1, local_n2;
    int ista, iend;
    int jsta, jend;

    double *buf_send, *buf_recv;

    int map_periods[2];
    MPI_Comm cart_comm;
} MPI_INFO;

void mpi_tool_setmap_2D(int rank, int nprocs, MPI_INFO *mpi_info);
void mpi_tool_setsubdomain_2D(int rank, int n1, int n2, MPI_INFO *mpi_info);
void mpi_tool_monitor(MPI_INFO mpi_info);
void mpi_tool_clean(MPI_INFO *mpi_info);

void mpi_tool_communicate_w2e(MPI_INFO mpi_info, double *x, int n1sub, int n2sub);
void mpi_tool_communicate_e2w(MPI_INFO mpi_info, double *x, int n1sub, int n2sub);
void mpi_tool_communicate_s2n(MPI_INFO mpi_info, double *x, int n1sub, int n2sub);
void mpi_tool_communicate_n2s(MPI_INFO mpi_info, double *x, int n1sub, int n2sub);
void mpi_tool_communicate_allreduce(MPI_INFO mpi_info, double local_value, double *global_value);


void mpi_tool_cart_2D(int rank, int nprocs, MPI_INFO *mpi_info);