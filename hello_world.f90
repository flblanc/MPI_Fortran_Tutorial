PROGRAM hello_world
    implicit none
    
    ! compilation:
    ! mpif90 -o hello_world hello_world.f90 -lmpi
    
    include '/usr/include/mpi/mpif.h'
    integer ierr
    
    call MPI_INIT ( ierr )
    print *, "Hello World"
    call MPI_FINALIZE ( ierr )
    
    stop
end program hello_world