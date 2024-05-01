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

program main
    use FileReader
    use KeyboardReader
    implicit none
    real, allocatable :: x(:), y(:)
    integer :: n
    integer :: func_number
    integer :: i
    do
        print *, '1. Read from file'
        print *, '2. Input values by myself'
        print *, 'Please input number what yo want to do:'
        read(*,*) func_number
        if (func_number > 2 .or. func_number < 1) then
            print *, 'Error: Uncorrect input, try again'
        else if ( func_number == 1 ) then
            call file_input('test.txt', x, y, n)
            if (n < 8 .or. n > 12) then
                print *, 'Points count incorrect shold to be from 8 till 12'
            else
                print *, 'All data correct'
                print *, "x:", x
                print *, "y:", y
                print *, "n:", n 
            end if
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
