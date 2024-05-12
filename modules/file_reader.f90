module FileReader
    implicit none
    
contains
    subroutine file_input(file, x, y, n)
        character(len=*), intent(in) :: file
        real, allocatable, intent(out) :: x(:), y(:)
        integer, intent(out) :: n
        integer :: i
        real :: a,b
        open(unit=10, file = file, status='old', action='read')
        n = 0
        do
            read(10, *, end=100) a, b
            n = n+1
        end do
    100 continue
        allocate(x(n), y(n))
        rewind(10)
        do i=1,n
            read(10,*) x(i), y(i)
        end do

        close(10)
    end subroutine file_input
end module FileReader