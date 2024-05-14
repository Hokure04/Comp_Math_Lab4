module FileReader
    implicit none
    
contains
    subroutine file_input(x, y, n)
        ! character(len=*), intent(in) :: file
        real, allocatable, intent(out) :: x(:), y(:)
        integer, intent(out) :: n
        character(len=100) :: file_name
        integer :: i, file_found
        real :: a,b

        print *, "Enter path to file: "
        read(*, '(A)') file_name
        open(unit=10, file = file_name, status='old', action='read', iostat=file_found)

        if (file_found /= 0) then
            print *, "File not found"
            return
        end if

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