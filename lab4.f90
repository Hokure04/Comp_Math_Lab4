
module ExponentialApproximation
    use Parameters
    implicit none
    
contains

    subroutine solve_linear_system(matrix, rhs, a, b)
        real, intent(in) :: matrix(2,2)
        real, intent(in) :: rhs(2)
        real, intent(out) :: a, b
        real :: det

        det = matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1)

        if (abs(det) < 1e-10) then
            print *, "Error: Matrix is singular or near-singular. Cannot solve the system."
            return
        end if

        b = (matrix(2,2) * rhs(1) - matrix(1,2) * rhs(2)) / det
        a = (matrix(1,1) * rhs(2) - matrix(2,1) * rhs(1)) / det

    end subroutine solve_linear_system

    subroutine exp_approximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx = 0, sx_2 = 0, sxy = 0, sy = 0, m, a, b
        real, allocatable :: P_x(:), e_i(:)
        integer :: i, non_positive_y = 0 
        real :: matrix(2,2)
        real :: rhs(2)

        print *, "Exponential approximation"

        do i=1,n
            if (y(i) <= 0) then
                print *, y(i)
                non_positive_y = 1
            endif
        end do
        

        if (non_positive_y == 1) then
            print *, "ERROR: Exponential approxiamtion isn't possible, file mustn't have y less or equal zero"
        else
            m = n
            sx = sum(x)
            sx_2 = sum(x**2)
            sy = sum(log(y))
            sxy = sum(x*log(y))

            matrix = reshape([sx_2, sx, & 
                            sx, m], [2, 2])

            rhs = [sxy, sy]

            call solve_linear_system(matrix, rhs, a, b)
            
            a = exp(a)
            allocate(P_x(n), e_i(n))
            do i=1,n
                P_x(i) = a*exp(b*x(i))
                e_i(i) = P_x(i) - y(i)
            end do
            print *, 'X values:', x
            print *, 'Y values:', Y
            print *, 'P_x:', P_x
            print *, 'e_i:', e_i
            call  calc_determination(y, P_x, n)
            call calc_deviation(y, P_x, n)
        endif
        print *, ""
        
    end subroutine
end module ExponentialApproximation


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
