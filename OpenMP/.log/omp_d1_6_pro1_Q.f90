program main
    ! <<< !!
    implicit none

    integer, parameter :: N = 20
    integer :: A(0:N-1)

    integer :: i
    integer :: thread_id, num_threads
    integer :: chunk_size
    integer :: i_start, i_end

    !$omp parallel private(i,thread_id, num_threads, chunk_size, i_start, i_end)
        thread_id   = 0 ! <<< !!
        num_threads = 0 ! <<< !!
        chunk_size  = 0 ! <<< !!
        i_start     = 0 ! <<< !!
        i_end       = 0 ! <<< !!

        do i = 0, 0 ! <<< !!
            A(i) = 100 + i
        end do
    !$omp end parallel

    do i = 0, N-1
        write(*,'(I8,1X,I4.4)') i, A(i) 
    end do
end program main



