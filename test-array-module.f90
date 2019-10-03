module test_array_module

    implicit none
    contains

    subroutine TEST(in_array, out_array)

        implicit none

        real(8), intent(in) :: in_array(:)
        real(8), intent(out), dimension(:), allocatable :: out_array

        integer n

        n=size(in_array)

        allocate(out_array(n))

        out_array = 2*in_array

        

    end subroutine TEST

end module
