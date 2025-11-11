program main
    use mpi
    implicit none
    integer :: ierr, rank, size,i 
    real*8 :: A(1:9), B(1:9)

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Initialize arrays A and B
    A = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/)*100+rank
    B = 0.0

    ! if ( rank==0 ) then
    !     call MPI_Send(A(1), 4, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
    ! else if ( rank==1 ) then
    !     call MPI_Recv(B(1), 4, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    ! end if
    if ( rank==0 ) then
        call MPI_Send(A(1), 4, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
    else if ( rank==1 ) then
        call MPI_Recv(B(1), 2, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    do i = 0, size-1
        if ( rank==i ) then
            write(*,"(1A5,1I5,1A7,9F8.1)") "Rank ", rank, ": B = ", B
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)        
    end do

    call MPI_Finalize(ierr)
end program main