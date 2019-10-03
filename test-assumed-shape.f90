program test

      implicit none

      real(8) :: my_array(4,8)
      integer, dimension(:), allocatable :: my_shape
      integer :: a,b

      my_shape = shape(my_array)

      a=my_shape(1)
      b=my_shape(2)


      write(*,*) my_shape
      end program test
