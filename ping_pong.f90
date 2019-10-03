! Florian Blanc August 1st 2017


program ping_pong

    implicit none
    include '/usr/include/mpi/mpif.h'
    
    integer ierr
    
    integer world_size
    integer world_rank
    
    integer partner_rank
    integer :: exchange_counter=0, max_exchange=10
    
    
    partner_rank = MODULO(world_rank+1,2)
    
    
    call MPI_INIT ( ierr )
    
    call MPI_COMM_SIZE( MPI_COMM_WORLD, world_size, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, world_rank, ierr )
    
    ! Terminate if there are not exactly two processes
    if (world_size /= 2) then
        print * , "Too many processes, exiting."
        call MPI_FINALIZE ( ierr )
        stop
    end if 
    
    
    do while ( exchange_counter < max_exchange )
    
        if (MODULO(exchange_counter,2) == 0) then
            
            exchange_counter = exchange_counter + 1
            
            if ( world_rank == 0 ) then ! sender
                
                call MPI_SEND ( exchange_counter, 1, MPI_INTEGER, partner_rank, 0, MPI_COMM_WORLD, ierr)
                
                WRITE (*,*) "Process ",world_rank," sent data to process ", partner_rank
                
            else if ( world_rank == 1 ) then
            
                call MPI_RECV ( exchange_counter, 1, MPI_INTEGER, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            
                WRITE (*,*) "Process ",world_rank," received data from process ", partner_rank
            
            end if 
            
        else 
        
            exchange_counter = exchange_counter + 1
            
            if ( world_rank == 1 ) then ! sender
                
                call MPI_SEND ( exchange_counter, 1, MPI_INTEGER, partner_rank, 0, MPI_COMM_WORLD, ierr)
                
                WRITE (*,*) "Process ",world_rank," sent data to process ", partner_rank
                
            else if ( world_rank == 0 ) then
            
                call MPI_RECV ( exchange_counter, 1, MPI_INTEGER, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            
                WRITE (*,*) "Process ",world_rank," received data from process ", partner_rank
            
            end if 
            
        end if 
        
    end do
    

    !print *, world_rank
    
    call MPI_FINALIZE ( ierr )
    
    stop
    
end program ping_pong
    
    