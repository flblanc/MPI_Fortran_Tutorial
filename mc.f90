module mc
    implicit none
    
    
    contains 
    
    subroutine Potential(x,a,E)
    
        implicit none
        
        real(8), intent(in) :: x, a
        real(8), intent(out) :: E
        
        E = (x-a)*(x-a)*(x+a)*(x+a)
        
    end subroutine Potential
    
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


    subroutine Initialize_Random_Seed()
    
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n , un , istat, dt(8), pid 
        
        call random_seed(size=n)
        allocate(seed(n))
        
        call system_clock(t)
        pid = getpid()
        
        t = ieor(t,int(pid, kind(t)))
        
        call random_seed(put=seed)
        
            
        

end module mc