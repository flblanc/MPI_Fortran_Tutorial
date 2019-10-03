PROGRAM test_array

    use test_array_module

    implicit none
    real(8), dimension(4):: a1
    real(8), dimension(:), allocatable :: a2

    a1 = (/1.,2.,4.,8./)

    call TEST(a1,a2)

    write(*,*) a2

    deallocate(a2)

end program test_array
