program string_method

    use mpi
    use utility_functions

    implicit none

    integer, parameter :: dp=kind(0.d0)
    integer ierr
    integer :: world_rank, world_size
    integer :: t,m,k,r

    ! Simulation parameters

    integer :: nsteps=100
    integer :: n_iter=15

    integer :: n_images, n_cv

    integer :: is_accepted

    character(len=1) :: replica_id
    character(len=50) :: file_name
    integer unit_id

    integer seed

    real(8) , allocatable, dimension(:) :: x,dx
    real(8), allocatable, dimension(:,:) :: string_array, center_array
    real(8) :: energy,energy_new, bias_energy, bias_energy_new, xc, yc

    real(8) :: beta=1.0/5.0, delta=0.1, smooth_param=0.1

    ! set up MPI world
    call MPI_INIT ( ierr )

    call MPI_COMM_SIZE( MPI_COMM_WORLD, world_size, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, world_rank, ierr )

    ! Initialize output file for this replica
    replica_id = world_rank
    write(replica_id,'(i0)') world_rank
    file_name = 'image' // trim(adjustl(replica_id)) // '.log'


    unit_id = world_rank
    open ( unit = unit_id, file=trim(file_name) )

    ! Initialize the image to the value read in the input file

    allocate(string_array(n_images, n_cv))
    allocate(center_array(n_images, n_cv))
    allocate(x(n_cv))
    allocate(dx(n_cv))
    call read_string("string.initial.dat", string_array, n_images, n_cv)
    center_array = string_array ! Initialize bias centers to the initial string

    ! Starting point for current image
    x = string_array(world_rank+1,:)

    do m=1,niter


        xc = center_array(world_rank+1,1)
        yc = center_array(world_rank+1,2)

        do t=1,nsteps

            call Potential(x(1),x(2), X0,Y0, AA2, aa,bb,cc,energy)
            call HarmonicRestraint(x(1),x(2), forceConstant_x,forceConstant_y, xc, yc, bias_energy)

            call mc_move1(dx,delta)

            call Potential(x(1)+dx(1),x(2)+dx(2), X0,Y0, AA2, aa,bb,cc,energy_new)
            call HarmonicRestraint(x(1)+dx(1),x(2)+dx(2), forceConstant_x,forceConstant_y, xc,yc, bias_energy_new)

            call Metropolis_Criterion(energy+bias_energy, energy_new+bias_energy_new, beta, is_accepted)

            if ( is_accepted == 1) then

                x(:) = x(:) + dx(:)

            end if

            write(unit_id,'(I5.3A4F8.6A4F8.6)') ((m-1)*nsteps+t),"    ", x(1),"    ",x(2)
        end do

        ! wait until everyone is done
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ! Now the master node needs to update the string

        ! Naive for loop over replicas; the binary tree will be implemented later

        if ( world_rank /= 0) then

            do k=1,n_cv

                call MPI_SEND(( x(k), 1, MPI_DOUBLE_PRECISION, 0, world_rank, MPI_COMM_WORLD, ierr ))

            end do

            do k=1,n_cv

                call MPI_RECV( x(k))
                
            end do


        else

            string_array(1,:) = x(:)

            do r=1,world_size-1
                do k=1,n_cv

                    call MPI_RECV( string_array(r+1,k), 1, MPI_DOUBLE_PRECISION, r, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

                end do

            end do

            ! String values are gathered and the master node can proceed to smoothing and reparametrization

            call smooth_string(string_array, smooth_param)
            call reparametrize_string(string_array)
