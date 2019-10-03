module utility_functions

    implicit none

    contains

        


! Mueller-Brown potential and energetics

        subroutine PotentialComponent(x,y,x0,y0,A2,a,b,c,E)

            implicit none
            real(8), intent(in) :: x,y
            real(8), intent(in) :: x0(:),y0(:),A2(:),a(:),b(:),c(:)
            real(8), intent(out) :: E(:)

            real(8), dimension(:), allocatable :: paraboloid
            
            
            allocate(paraboloid(size(x0)))
            paraboloid = a*(x-x0(:))**2 + b*(x-x0(:))*(y-y0(:)) + c*(y-y0(:))**2

            E(:) = A2(:)*exp(paraboloid(:))
            deallocate(paraboloid)

        end subroutine PotentialComponent

        subroutine Potential(x,y, X0,Y0, AA2, aa,bb,cc,E)

            real(8), intent(in) :: x,y
            real(8), dimension(4), intent(in) :: X0, Y0, AA2, aa, bb, cc

            real(8), intent(out) :: E
            real(8), dimension(4) :: ee

            call PotentialComponent(x,y,X0(:),Y0(:),AA2(:),aa(:),bb(:),cc(:),ee(:))

            E = sum(ee)

        end subroutine Potential

        subroutine HarmonicRestraint(x,y,forceConstant_x,forceConstant_y,xc,yc,V)

            implicit none
            real(8), intent(in) :: x,y,xc,yc,forceConstant
            real(8), intent(out) :: V

            V = 0.5*forceConstant_x*(x-xc)**2 + 0.5*forceConstant_y*(y-yc)**2

        end subroutine HarmonicRestraint




! String method utilities

        subroutine smooth_string(current_string, smooth_param)

            implicit none

            integer :: n_images, n_cv
            real(8), intent(in) :: smooth_param
            real(8), intent(in out) :: current_string(:,:)

            real(8), dimension(:,:), allocatable :: smoothed_string

            integer :: i, j

            n_images = shape(current_string)(1)
            n_cv = shape(current_string)(2)

            allocate(smoothed_string(n_images,n_cv))

            smoothed_string(1,:) = current_string(1,:)
            smoothed_string(n_images,:) = current_string(n_images,:)

            do j=2,n_images-1

                smoothed_string(j,:) = (1.-smooth_param)*current_string(j,:) & 
                    + 0.5*smooth_param*(current_string(j-1,:)+current_string(j+1,:))

            end do

            current_string = smoothed_string
            deallocate(smoothed_string)

        end subroutine smooth_string

        subroutine reparametrize_string(current_string)

            implicit none

            real(8), intent(in out) :: current_string(:,:)
            real(8), dimension(:,:), allocatable :: reparametrized_string
            real(8), dimension(:,:), allocatable :: g1, g2
            real(8) :: dx, scaling_factor


            integer :: n_images, n_cv
            integer, dimension(2) :: array_shape

            integer :: i,j

            array_shape = shape(current_string)
            n_images = array_shape(1)
            n_cv = array_shape(2)

            allocate(reparametrized_string(n_images,n_cv))

            allocate(g1(n_images))
            allocate(g2(n_images))

            reparametrized_string(1,:) = current_string(1,:)
            reparametrized_string(n_images,:) = current_string(n_images,:)

            ! g1 is linspace(0,1,M)

            !  A = (/ (I, I = 1, 10) /)

            g1 = (/ (double(I)/(double(n_images)-1.0)), I=0,n_images-1 /)

            ! g2 contains the unevenly spaced arc-lengths
            g2(1) = 0.0

            do i=2,n_images

                dx = sum((current_string(i,:)-current_string(i-1,:))**2)
                dx = sqrt(dx)

                g2(i) = g2(i-1) + dx

            end do

            g2 = g2 / g2(n_images)

            ! Once g2 is built, we interpolate the functions cv = S(g2)
            ! then compute cv_reparam = S(g1)

            do j=2,n_images-1

                if ( g1(i) < g2(j)) then

                    exit

                end if

            end do

            scaling_factor = (current_string(j,:) - current_string(j-1,:)) / (g2(j) - g2(j-1))

            reparametrized_string(i,:) = scaling_factor*(g1(i) - g2(j-1)) + current_string(j-1,:)

            current_string = reparametrized_string
            deallocate(reparametrized_string)

        end subroutine reparametrize_string

        subroutine read_string(string_file_name, string_array, n_images, n_cv)

            implicit none

            real(8), allocatable, dimension(:,:), intent(out) :: string_array
            integer, intent(in) :: n_images, n_cv
            character(len=*), intent(in) :: string_file_name

            integer :: unit_id=1000, i,j
            integer, dimension(n_images) :: ii

            allocate(string_array(n_images,n_cv))
            open(unit=unit_id, file=string_file_name)

            do i=1,n_images

                read(unit_id,*) ii(:), string_array(i,:)

            end do

            close(unit=unit_id)

        end subroutine read_string




! Monte Carlo dynamics

        subroutine Metropolis_Criterion(E_old, E_new, beta, is_accepted)
            ! Evaluate Metropolis criterion
            ! Return: is_accepted (boolean)

            implicit none

            real(8), intent(in) :: E_old ! unit is kcal/mol
            real(8), intent(in) :: E_new
            real(8) :: dE=0.0, boltzmann_factor=1.0


            real(8), intent(in) :: beta ! unit is 1/(kcal/mol)

            real(8) :: draw

            integer, intent(out) :: is_accepted

            ! Compute energy variation
            dE = E_new - E_old

            if ( dE <= 0.0 ) then
                is_accepted = 1
            else

                ! call random_seed
                call random_number(draw)

                ! Compute Boltzmann factor
                boltzmann_factor = EXP(-beta*dE)
                if ( boltzmann_factor >= draw ) then
                    ! Move accepted
                    is_accepted = 1
                else
                    is_accepted = 0
                endif
            endif

            ! if (is_accepted == 0) then
            !     print *, "Move Rejected"
            ! else
            !     print *, "Move Accepted"
            ! endif
        end subroutine Metropolis_Criterion



    ! MC Move

    subroutine mc_move1 ( dx, delta )


        real(8), intent(in) :: delta
        real(8), intent(out) :: dx(:)

        real(8), dimension(:), allocatable :: draw


        integer :: i,n

        n= size(dx)


        allocate(draw(n))

        do i=1,n

            call random_number(dx(i))

            dx(i) = delta*(2.*dx(i) - 1.0)


        end do

    end subroutine mc_move1
    
    
end module utility_functions
