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




module Matrix3x3
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
end module Matrix3x3


module Matrix4x4
    implicit none

contains

    subroutine solve_linear_system(matrix, rhs, solution)
        real, intent(in) :: matrix(4,4), rhs(4)
        real, intent(out) :: solution(4)
        real :: det

        ! Находим определитель матрицы
        det = determinant(matrix)

        ! Проверка на невырожденность матрицы
        if (det == 0.0) then
            print *, 'ERROR: The matrix of the system is degenerate, the solution is impossible'
            return
        endif

        ! Вычисление коэффициентов с помощью метода Крамера
        solution(1) = determinant(get_minor(matrix, 1, 1)) / det
        solution(2) = determinant(get_minor(matrix, 2, 1)) / det
        solution(3) = determinant(get_minor(matrix, 3, 1)) / det
        solution(4) = determinant(get_minor(matrix, 4, 1)) / det

    end subroutine solve_linear_system

    ! Функция для вычисления определителя матрицы
    function determinant(matrix) result(det)
        real, intent(in) :: matrix(4,4)
        real :: det
        real :: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p
        
        a = matrix(1,1)
        b = matrix(1,2)
        c = matrix(1,3)
        d = matrix(1,4)
        e = matrix(2,1)
        f = matrix(2,2)
        g = matrix(2,3)
        h = matrix(2,4)
        i = matrix(3,1)
        j = matrix(3,2)
        k = matrix(3,3)
        l = matrix(3,4)
        m = matrix(4,1)
        n = matrix(4,2)
        o = matrix(4,3)
        p = matrix(4,4)
        
        det = a * (f * (k * p - l * o) - g * (j * p - l * n) + h * (j * o - k * n)) &
            - b * (e * (k * p - l * o) - g * (i * p - l * m) + h * (i * o - k * m)) &
            + c * (e * (j * p - l * n) - f * (i * p - l * m) + h * (i * n - j * m)) &
            - d * (e * (j * o - k * n) - f * (i * o - k * m) + g * (i * n - j * m))
            
    end function determinant

    ! Вспомогательная функция для получения минора матрицы
    pure function get_minor(matrix, row, col) result(minor)
        real, intent(in) :: matrix(4,4)
        integer, intent(in) :: row, col
        real :: minor(3,3)
        integer :: i, j, m, n

        m = 1
        do i = 1, 4
            if (i == row) then
                cycle
            endif
            n = 1
            do j = 1, 4
                if (j == col) then
                    cycle
                endif
                minor(m,n) = matrix(i,j)
                n = n + 1
            end do
            m = m + 1
        end do

    end function get_minor

end module Matrix4x4


module QuadraticApproximation
    use DeterminationCoefficient
    use StandatrDeviation
    use Matrix3x3
    implicit none
    
contains
    
    subroutine quadratic_approximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx = 0, sx_2 = 0, sx_3 = 0, sx_4 = 0, sy = 0, sxy = 0, sx_2y = 0
        real :: matrix(3,3)
        real :: rhs(3)
        real, allocatable :: a(:), P_x(:), e_i(:)
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
        print *, a(1)

        allocate(P_x(n), e_i(n))
        do i =1,n
            P_x(i) = a(1)+a(2)*x(i)+a(3)*x(i)**2
            e_i(i) = P_x(i) - y(i)
        end do
        print *, 'X values:', x
        print *, 'Y values: ', y
        print *, 'P_x values:', P_x
        print *, 'e_i values:', e_i
        call calc_determination(x, y, n)
        call calc_deviation(x, y, n)
    end subroutine
end module

module CubicApproximation
    use DeterminationCoefficient
    use StandatrDeviation
    use Matrix4x4
    implicit none
    
contains
    
subroutine cubic_approximation(x, y, n)
    real, allocatable, intent(in) :: x(:), y(:)
    integer, intent(in) :: n
    real :: sx = 0, sx_2 = 0, sx_3 = 0, sx_4 = 0, sy = 0, sxy = 0, sx_2y = 0, sx_4y = 0, m
    real :: matrix(4,4)
    real :: rhs(4)
    real, allocatable :: a(:), P_x(:), e_i(:)
    integer :: i

    m = n 

    do i=1,n
        sx = sx + x(i)
        sx_2 = sx_2 + x(i)**2
        sx_3 = sx_3 + x(i)**3
        sx_4 = sx_4 + x(i)**4
        sy = sy + y(i)
        sxy = sxy + x(i)*y(i)
        sx_2y = sx_2y + x(i)**2*y(i)
        sx_4y = sx_4y + x(i)**4*y(i)
    end do

    matrix = reshape([m, sx, sx_2, sx_3, &
                      sx, sx_2, sx_3, sx_4, &
                      sx_2, sx_3, sx_4, sx_4*sx_4, &
                      sx_3, sx_4, sx_4*sx_4, sx_4**2], [4, 4])

    rhs = [sy, sxy, sx_2y, sx_4y]

    allocate(a(4))
    call solve_linear_system(matrix, rhs, a)

    allocate(P_x(n), e_i(n))

    do i=1,n
        P_x(i) = a(1) + a(2)*x(i) + a(3)*x(i)**2 + a(4)*x(i)**3
        e_i(i) = P_x(i) - y(i)
    end do

    print *, 'X values: ', x
    print *, 'Y values: ', y
    print *, 'P_x values: ', P_x
    print *, 'e_i values: ', e_i

    call calc_determination(x, y, n)
    call calc_deviation(x, y, n)

end subroutine cubic_approximation
end module

module ExponentialApproximation
    implicit none
    
contains

    subroutine exp_approximation(x, y, n)
        
    end subroutine
end module ExponentialApproximation

program main
    use FileReader
    use KeyboardReader
    use LinearAproximation
    use QuadraticApproximation
    use CubicApproximation
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
            ! call linear_approximation(x,y,n)
            call quadratic_approximation(x,y,n)
            call cubic_approximation(x, y, n)
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
