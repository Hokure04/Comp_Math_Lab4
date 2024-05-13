module Matrix4x4
    implicit none
    
contains

    subroutine solve_linear_system(matrix, rhs, solution)
        real, intent(in) :: matrix(4,4)
        real, intent(in) :: rhs(4)
        real, intent(out) :: solution(4)
        real :: local_matrix(4,4)  ! Declare local matrix variable
        real :: local_rhs(4)       ! Declare local rhs variable
        real :: pivot, factor
        integer :: i, j, k

        local_matrix = matrix  ! Copy input matrix to local variable
        local_rhs = rhs        ! Copy input rhs to local variable

        do k = 1, 4
            pivot = local_matrix(k,k)
            do j = k+1, 4
                factor = local_matrix(j,k) / pivot
                local_matrix(j,k:4) = local_matrix(j,k:4) - factor * local_matrix(k,k:4)
                local_rhs(j) = local_rhs(j) - factor * local_rhs(k)
            end do
        end do

        solution(4) = local_rhs(4) / local_matrix(4,4)
        do i = 3, 1, -1
            solution(i) = (local_rhs(i) - dot_product(local_matrix(i,i+1:4), solution(i+1:4))) / local_matrix(i,i)
        end do
    end subroutine solve_linear_system

end module Matrix4x4


! module Matrix4x4
!     implicit none

! contains

!     subroutine solve_linear_system(matrix, rhs, solution)
!         real, intent(in) :: matrix(4,4), rhs(4)
!         real, intent(out) :: solution(4)
!         real :: det

!         ! Находим определитель матрицы
!         det = determinant(matrix)

!         ! Проверка на невырожденность матрицы
!         if (abs(det) < 1e-10) then
!             print *, 'ERROR: The matrix of the system is degenerate, the solution is impossible'
!             return
!         endif

!         ! Вычисление коэффициентов с помощью метода Крамера
!         solution(1) = determinant(get_minor(matrix, 1, 1)) / det
!         solution(2) = -determinant(get_minor(matrix, 1, 2)) / det
!         solution(3) = determinant(get_minor(matrix, 1, 3)) / det
!         solution(4) = -determinant(get_minor(matrix, 1, 4)) / det

!     end subroutine solve_linear_system

!     ! Функция для вычисления определителя матрицы
!     function determinant(matrix) result(det)
!         real, intent(in) :: matrix(4,4)
!         real :: det
!         real :: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p
        
!         a = matrix(1,1)
!         b = matrix(2,1)
!         c = matrix(3,1)
!         d = matrix(4,1)
!         e = matrix(1,2)
!         f = matrix(2,2)
!         g = matrix(3,2)
!         h = matrix(4,2)
!         i = matrix(1,3)
!         j = matrix(2,3)
!         k = matrix(3,3)
!         l = matrix(4,3)
!         m = matrix(1,4)
!         n = matrix(2,4)
!         o = matrix(3,4)
!         p = matrix(4,4)
        
!         det = a * (f * (k * p - l * o) - g * (j * p - l * n) + h * (j * o - k * n)) &
!             - e * (b * (k * p - l * o) - c * (j * p - l * n) + d * (j * o - k * n)) &
!             + i * (b * (g * p - h * o) - c * (f * p - h * m) + d * (f * o - g * m)) &
!             - m * (b * (g * l - h * k) - c * (f * l - h * i) + d * (f * k - g * i))
            
!     end function determinant

!     ! Вспомогательная функция для получения минора матрицы
!     pure function get_minor(matrix, row, col) result(minor)
!         real, intent(in) :: matrix(4,4)
!         integer, intent(in) :: row, col
!         real :: minor(3,3)
!         integer :: i, j, m, n

!         m = 1
!         do i = 1, 4
!             if (i == row) then
!                 cycle
!             endif
!             n = 1
!             do j = 1, 4
!                 if (j == col) then
!                     cycle
!                 endif
!                 minor(m,n) = matrix(i,j)
!                 n = n + 1
!             end do
!             m = m + 1
!         end do

!     end function get_minor

! end module Matrix4x4

module CubicApproximation
    use Parameters
    use Matrix4x4
    implicit none
    
contains
    
    subroutine cubic_approximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx = 0, sx_2 = 0, sx_3 = 0, sx_4 = 0, sx_5 = 0, sx_6 = 0, sy = 0, sxy = 0, sx_2y = 0, sx_3y = 0, m
        real :: matrix(4,4)
        real :: rhs(4)
        real, allocatable :: a(:), P_x(:), e_i(:)
        integer :: i

        print *, "Cubic Approximation"

        m = n 
        sx = sum(x)
        sx_2 = sum(x**2)
        sx_3 = sum(x**3)
        sx_4 = sum(x**4)
        sx_5 = sum(x**5)
        sx_6 = sum(x**6)
        sy = sum(y)
        sxy = sum(x*y)
        sx_2y = sum(x**2 * y)
        sx_3y = sum(x**3 * y)

        matrix = reshape([m, sx, sx_2, sx_3, &
                          sx, sx_2, sx_3, sx_4, &
                          sx_2, sx_3, sx_4, sx_5, &
                          sx_3, sx_4, sx_5, sx_6], [4, 4])

        rhs = [sy, sxy, sx_2y, sx_3y]

        allocate(a(4))
        call solve_linear_system(matrix, rhs, a)

        print *, "Coefficients a: ", a

        allocate(P_x(n), e_i(n))

        do i=1,n
            P_x(i) = a(1) + a(2)*x(i) + a(3)*x(i)**2 + a(4)*x(i)**3
            e_i(i) = P_x(i) - y(i)
        end do

        print *, 'X values: ', x
        print *, 'Y values: ', y
        print *, 'P_x values: ', P_x
        print *, 'e_i values: ', e_i

        call calc_determination(y, P_x, n)
        call calc_deviation(y, P_x, n)

        print *, ""

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
            call file_input('tests/present.txt', x, y, n)
            if (n < 8 .or. n > 12) then
                print *, 'Points count incorrect shold to be from 8 till 12'
            else
                print *, 'All data correct'
                print *, "x:", x
                print *, "y:", y
                print *, "n:", n 
            end if
            call linear_aproximation(x,y,n)
            call quadratic_aproximation(x,y,n)
            call cubic_approximation(x, y, n)
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
