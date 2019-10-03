program test_random_seed


    ! HOW TO: initialize the RNG such that each process has a different seed

    ! from https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    ! https://stackoverflow.com/questions/18754438/generating-random-numbers-in-a-fortran-module
    ! and my own modification to adapt it to MPI
    use mpi
    implicit none
    integer, allocatable :: seed(:)
    integer :: n, i, clock


    integer ierr
    integer world_rank
    integer world_size
    
    real(8) dx
    
    call MPI_INIT ( ierr )
    
    call MPI_COMM_SIZE( MPI_COMM_WORLD, world_size, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, world_rank, ierr )

    
    
    call random_seed(size = n)
!~     print *, n
!~     n = n + 100*world_rank
    !print *, n
    allocate(seed(n))
    call random_seed(get=seed)
    
    
    call system_clock(count=clock)
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /) + world_rank
    
    call random_seed(put=seed)
    !write (*,*) seed
    
    call random_number(dx)
    
    print *, dx
    
    deallocate(seed)
    
    call MPI_FINALIZE ( ierr )
  
end program test_random_seed