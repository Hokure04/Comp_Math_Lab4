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

module CorrelationCoefficient
    implicit none

contains
    subroutine calc_correlation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: r, mean_x, mean_y, sum_x, sum_y, numerator
        integer :: i
        sum_x = 0
        sum_y = 0
        do i=1,n
            sum_x = sum_x + x(i)
            sum_y = sum_y + y(i)
        end do
        mean_x = sum_x/n
        mean_y = sum_y/n
        sum_x = 0
        sum_y = 0
        do i=1,n
            sum_x = sum_x + (x(i)-mean_x)**2
            sum_y = sum_y + (y(i)-mean_y)**2
            numerator = numerator + (x(i) - mean_x)*(y(i) - mean_y)
        end do
        r = numerator/sqrt(sum_x*sum_y)
        print *, 'Correlation coefficient: ', r
    end subroutine
end module CorrelationCoefficient

module DeterminationCoefficient
    implicit none
    
contains

    subroutine calc_determination(y, fi, n)
        real, allocatable, intent(in) :: y(:), fi(:)
        integer, intent(in) :: n
        real :: mean_fi, numerator, denominator, determination
        integer :: i
        mean_fi = 0
        do i =1,n
            mean_fi = mean_fi + fi(i)
        end do
        mean_fi = mean_fi/n
        
        numerator = 0
        denominator = 0
        do i=1,n
            numerator = numerator + (y(i) - fi(i))**2
            denominator = denominator + (y(i)-mean_fi)**2
        end do
        determination = 1 - numerator/denominator
        print *,'Determination coefficient: ', determination
    end subroutine
end module DeterminationCoefficient

module StandatrDeviation
    implicit none
    
contains

    subroutine calc_deviation(y, fi, n)
        real, allocatable, intent(in) :: y(:), fi(:)
        integer, intent(in) :: n
        real :: sum, deviation
        integer :: i
        do i=1,n
            sum = sum + (fi(i) - y(i))**2
        end do
        deviation = sqrt(sum/n)
        print *, 'Standart deviation: ', deviation
    end subroutine
end module StandatrDeviation

module LinearAproximation
    use CorrelationCoefficient
    use DeterminationCoefficient
    use StandatrDeviation
    implicit none
    
contains
    subroutine linear_approximation(x,y,n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx, sxx, sy, sxy, delta, delta1, delta2, a, b
        real, allocatable :: P_x(:), e_i(:)
        integer :: i
        sx = 0.0
        sxx = 0.0
        sy = 0.0
        sxy = 0.0
        do i =1,n
            sx = sx+x(i)
            sxx = sxx+x(i)**2
            sy = sy+y(i)
            sxy = sxy + y(i)*x(i)
        end do
        print *, 'SX:', sx, 'SXX:', sxx, 'SY:', sy, 'SXY', sxy
        delta = sxx*n - sx*sx
        delta1 = sxy*n -sx*sy
        delta2 = sxx*sy - sx*sxy
        print *, 'delta:', delta, 'delta1:', delta1, 'delta2:', delta2
        if(delta /= 0) then
            a = delta1/delta
            b = delta2/delta
            print *, 'Value a:', a
            print *, 'Value b:', b
            allocate(P_x(n))
            allocate(e_i(n))
            do i=1,n
                P_x(i) = a*x(i)+b
                e_i(i) = P_x(i) - y(i)
            end do
            print *, 'P_x:',P_x
            print *,'e_i:', e_i
            call calc_correlation(x,y,n)
            call calc_determination(y, P_x, n)
            call calc_deviation(y,P_x,n)
        else
            print *, 'ERROR: Delta parameter is 0, but we cannnot divide by zero'
        end if
    end subroutine
end module LinearAproximation
module Matrix
    implicit none
    
contains

function determinant(matrix) result(det)
    real, intent(in) :: matrix(:,:)
    real :: det
    det = matrix(1,1) * (matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
        - matrix(1,2) * (matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
        + matrix(1,3) * (matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
end function determinant

subroutine inverse_matrix(matrix, inverse)
    real, intent(in) :: matrix(:,:)
    real, intent(out) :: inverse(:,:)
    real :: det
    ! integer :: i, j

    det = determinant(matrix)
    if (det == 0.0) then
        print *, 'ERROR: The matrix is singular, inverse cannot be computed.'
        return
    endif

    inverse(1,1) = (matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) / det
    inverse(1,2) = -(matrix(1,2)*matrix(3,3) - matrix(1,3)*matrix(3,2)) / det
    inverse(1,3) = (matrix(1,2)*matrix(2,3) - matrix(1,3)*matrix(2,2)) / det
    inverse(2,1) = -(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) / det
    inverse(2,2) = (matrix(1,1)*matrix(3,3) - matrix(1,3)*matrix(3,1)) / det
    inverse(2,3) = -(matrix(1,1)*matrix(2,3) - matrix(1,3)*matrix(2,1)) / det
    inverse(3,1) = (matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)) / det
    inverse(3,2) = -(matrix(1,1)*matrix(3,2) - matrix(1,2)*matrix(3,1)) / det
    inverse(3,3) = (matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)) / det
end subroutine inverse_matrix

subroutine solve_linear_system(matrix, rhs, solution)
    real, intent(in) :: matrix(:,:), rhs(:)
    real, intent(out) :: solution(:)
    real :: det
    det = determinant(matrix)
    
    if (det == 0.0) then
        print *, 'ERROR: The matrix of the system is degenerate, the solution is impossible'
        return
    endif
    
    solution(1) = (rhs(1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
        - matrix(1,2)*(rhs(2)*matrix(3,3) - matrix(2,3)*rhs(3)) &
        + matrix(1,3)*(rhs(2)*matrix(3,2) - matrix(2,2)*rhs(3))) / det
    
    solution(2) = (matrix(1,1)*(rhs(2)*matrix(3,3) - matrix(2,3)*rhs(3)) &
        - rhs(1)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
        + matrix(1,3)*(matrix(2,1)*rhs(3) - rhs(2)*matrix(3,1))) / det
    
    solution(3) = (matrix(1,1)*(matrix(2,2)*rhs(3) - rhs(2)*matrix(3,2)) &
        - matrix(1,2)*(matrix(2,1)*rhs(3) - rhs(2)*matrix(3,1)) &
        + rhs(1)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)) ) / det
end subroutine solve_linear_system
end module Matrix

module QuadraticApproximation
    use DeterminationCoefficient
    use StandatrDeviation
    use Matrix
    implicit none
    
contains
    
    subroutine quadratic_approximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx = 0, sx_2 = 0, sx_3 = 0, sx_4 = 0, sy = 0, sxy = 0, sx_2y = 0
        real :: matrix(3,3)
        real :: rhs(3)
        real, allocatable :: a(:)
        integer :: i

        do i=1,n
            sx = sx + x(i)
            sx_2 = sx_2 + x(i)**2
            sx_3 = sx_3 + x(i)**3
            sx_4 = sx_4 + x(i)**4
            sy = sy + y(i)
            sxy = sxy + x(i)*y(i)
            sx_2y = sx_2y + x(i)**2*y(i)
        end do
        print *, 'SX=', sx, 'SX_2=', sx_2, 'SX_3=', sx_3, 'SX_4=', sx_4
        print *, 'SY=', sy, 'SXY=', sxy, 'SX_2Y=', sx_2y

        matrix(1,1) = n
        matrix(1,2) = sx
        matrix(1,3) = sx_2
        matrix(2,1) = sx
        matrix(2,2) = sx_2
        matrix(2, 3) = sx_3
        matrix(3, 1) = sx_2
        matrix(3, 2) = sx_3
        matrix(3, 3) = sx_4

        rhs(1) = sy
        rhs(2) = sxy
        rhs(3) = sx_2y

        allocate(a(3))
        call solve_linear_system(matrix, rhs, a)
        print *, a
    end subroutine
end module

program main
    use FileReader
    use KeyboardReader
    use LinearAproximation
    use QuadraticApproximation
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
            call file_input('present2.txt', x, y, n)
            if (n < 8 .or. n > 12) then
                print *, 'Points count incorrect shold to be from 8 till 12'
            else
                print *, 'All data correct'
                print *, "x:", x
                print *, "y:", y
                print *, "n:", n 
            end if
            call linear_approximation(x,y,n)
            call quadratic_approximation(x,y,n)
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
