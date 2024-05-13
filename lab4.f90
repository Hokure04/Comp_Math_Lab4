


program main
    use FileReader
    use KeyboardReader
    use LinearAproximation
    use QuadraticApproximation
    use CubicApproximation
    use ExponentialApproximation
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
            print *, 'ERROR: Uncorrect input, try again'
        else if ( func_number == 1 ) then
            call file_input('tests/present.txt', x, y, n)
            if (n < 8 .or. n > 12) then
                print *, 'Points count incorrect shold to be from 8 till 12'
            else
                print *, 'All data correct'
                print *, "x:", x
                print *, "y:", y
                print *, "n:", n 
                print *, ""
            end if
            call linear_aproximation(x,y,n)
            call quadratic_aproximation(x,y,n)
            call cubic_approximation(x, y, n)
            call exp_approximation(x, y, n)
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
