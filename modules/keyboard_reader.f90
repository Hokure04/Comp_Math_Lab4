module KeyboardReader
    implicit none

contains
    subroutine keyboard_input(x, y, n)
        real, allocatable, intent(out) :: x(:), y(:)
        integer, intent(out) :: n
        integer :: i
        print *, 'Enter array size:'
        read(*,*) n

        allocate(x(n))
        allocate(y(n))

        print *, 'Please input x values'
        do i = 1,n
            read(*,*) x(i)
        end do

        print *, 'Please input y values'
        do i=1,n
            read(*,*) y(i)
        end do
    end subroutine
end module KeyboardReader