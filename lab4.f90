program main
    use FileReader
    use KeyboardReader
    use LinearApproximation
    use QuadraticApproximation
    use CubicApproximation
    use ExponentialApproximation
    use PowerApproximation
    use LogarithmicApproximation
    implicit none
    real, allocatable :: x(:), y(:)
    integer :: n
    integer :: func_number
    integer :: i
    logical :: valid_input

    do
        valid_input = .false.
        do while (.not. valid_input)
            print *, '1. Read from file'
            print *, '2. Input values by myself'
            print *, 'Please input number what you want to do:'
            read(*,*, iostat=i) func_number
            if (i /= 0) then
                print *, 'ERROR: Incorrect input, try again'
                cycle
            end if
            if (func_number > 2 .or. func_number < 1) then
                print *, 'ERROR: Incorrect input, try again'
                cycle
            end if

            if (func_number == 1) then
                call file_input(x, y, n)
            else if (func_number == 2) then
                call keyboard_input(x, y, n)
            end if

            if (n < 7 .or. n > 12) then
                print *, 'ERROR: Points count incorrect, should be from 7 to 12'
                deallocate(x, y)
                cycle
            end if

            valid_input = .true.
            print *, 'All data correct'
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
            print *, ""
            call linear_approximation(x, y, n)
            call quadratic_approximation(x, y, n)
            call cubic_approximation(x, y, n)
            call exp_approximation(x, y, n)
            call power_approximation(x, y, n)
            call logarithmic_approximation(x, y, n)
        end do
    end do

end program
