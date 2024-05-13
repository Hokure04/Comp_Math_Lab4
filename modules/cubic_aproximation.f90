module Matrix4x4
    implicit none
    
contains

    subroutine solve_linear_system(matrix, rhs, solution)
        real, intent(in) :: matrix(4,4)
        real, intent(in) :: rhs(4)
        real, intent(out) :: solution(4)
        real :: local_matrix(4,4)  ! объявление локальной матрицы
        real :: local_rhs(4)       ! объявление локального вектора решений
        real :: pivot, factor
        integer :: i, j, k

        local_matrix = matrix  ! копируем матрицу
        local_rhs = rhs        ! копируем вектор решений

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