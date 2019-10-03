program parallel_MC
    
    use mpi
    use mc
    implicit none
    
!~     include '/usr/include/mpi/mpif.h'
    
    integer ierr
    
    integer world_size
    integer world_rank
    
    ! This program runs MC dynamics in a bistable quartic potential
    ! U(x) = (x-a)²(x+a)² for nsteps; after nsteps, all replicas send 
    ! their final position to the master process which takes the average
    ! then each replica starts a new MC dynamics from the average value;
    ! they do it for niter times
    
    
    
    integer :: nsteps=100
    integer :: niter=10
    integer :: is_accepted
    integer :: k,i,r,l
    real(8) :: beta=1./0.596, a=+5.0, E=0.0, E_new=0.0
    real(8) :: x, dx
    
    real(8), dimension(:), allocatable :: x_list
    
    
    
    character(len=1) :: replica_id
    character(len=50) :: file_name
    integer unit_id
    
    integer seed
    
    
    
    call MPI_INIT ( ierr )
    
    call MPI_COMM_SIZE( MPI_COMM_WORLD, world_size, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, world_rank, ierr )
    
    ! Prepare outputfile; one per replica
    
    ! HOW TO: cast a number into a string
    write(replica_id,'(i0)') world_rank
    
    ! HOW TO: concatenate strings
    file_name = 'replica' // trim(adjustl(replica_id)) // '.log'
    
    
    unit_id = world_rank
    open ( unit = unit_id, file=trim(file_name) )
    
    ! allocate memory for x_list
    allocate(x_list(world_size))
    
    !x = real(world_rank)
    
    x = -a*(1. - (real(world_rank)/real(world_size))) + ( real(world_rank)/real(world_size) )*a 
    
    
    
    
    
    
    do k=1,niter
    

        do i=1,nsteps
        

            ! Compute energy
            call Potential(x,a,E)
            ! Make MC move
            call random_number(dx)
            dx = 2.0*dx - 1.0

            
            ! Compute new energy
            call Potential(x+dx,a,E_new)
            
            call Metropolis_Criterion(E,E_new,beta,is_accepted)
            
            if ( is_accepted == 1 ) then
                !print *, "Accepted"
                x = x + dx
                
            end if
            
            write(unit_id,'(I5.3A4F8.6)') ((k-1)*nsteps+i),"    ", x
        
        end do
        
        ! wait until everyone is done
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
        ! MC dynamics complete
        ! Now all replicas send their final position to the master node
        
        if ( world_rank /= 0) then
        
            call MPI_SEND ( x, 1, MPI_REAL8, 0, world_rank, MPI_COMM_WORLD, ierr )
            call MPI_RECV ( x, 1, MPI_REAL8, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
        
        else
            
            x_list(1) = x
            
            ! gather the values
            do r=1,world_size-1
                
                call MPI_RECV( x_list(r+1), 1, MPI_REAL8, r, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
            
            end do 
            
            do l=2,world_size
                x = x + x_list(l)
            end do
            
            ! take the average
            x = x/real(world_size)
            
            ! send back the values
            do r=1,world_size-1
                
                call MPI_SEND ( x, 1, MPI_REAL8, r, r, MPI_COMM_WORLD, ierr )
            
            end do
        
        end if 
        
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        write(unit_id,*) '#Averaging done.'
    end do
    
    deallocate(x_list)
    close ( unit_id )
    call MPI_FINALIZE ( ierr )
    
    stop

end program parallel_MC
            
                
            
        
    
        