
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

module CubicApproximation
    use Parameters
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

! module ExponentialApproximation
!     implicit none
    
! contains

!     subroutine exp_approximation(x, y, n)
        
!     end subroutine
! end module ExponentialApproximation


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
            call file_input('tests/present2.txt', x, y, n)
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
            call cubic_approximation(x, y, n)
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
