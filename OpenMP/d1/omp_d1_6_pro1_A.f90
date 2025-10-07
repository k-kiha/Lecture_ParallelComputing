program main
    use omp_lib
    implicit none

    integer, parameter :: N = 20
    integer :: A(0:N-1)

    integer :: i
    integer :: thread_id, num_threads
    integer :: chunk_size
    integer :: i_start, i_end

    !$omp parallel private(i,thread_id, num_threads, chunk_size, i_start, i_end)
        thread_id   = omp_get_thread_num()
        num_threads = omp_get_num_threads()
        chunk_size  = N / num_threads
        i_start     = thread_id * chunk_size
        i_end       = (thread_id + 1) * chunk_size - 1   

        do i = i_start, i_end
            A(i) = 100 + i
        end do
    !$omp end parallel

    do i = 0, N-1
        write(*,'(I8,1X,I4.4)') i, A(i) 
    end do
end program main


