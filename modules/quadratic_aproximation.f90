
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


module QuadraticApproximation
    use Parameters
    use Matrix3x3
    implicit none
    
contains
    
    subroutine quadratic_aproximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: m, sx = 0, sx_2 = 0, sx_3 = 0, sx_4 = 0, sy = 0, sxy = 0, sx_2y = 0
        real :: matrix(3,3)
        real :: rhs(3)
        real, allocatable :: a(:), P_x(:), e_i(:)
        integer :: i

        print *, "Quadratic aproximation"
    
        m = n
        sx = sum(x)
        sx_2 = sum(x**2)
        sx_3 = sum(x**3)
        sx_4 = sum(x**4)
        sy = sum(y)
        sxy = sum(x*y)
        sx_2y = sum(x**2 * y)
        
        print *, 'SX=', sx, 'SX_2=', sx_2, 'SX_3=', sx_3, 'SX_4=', sx_4
        print *, 'SY=', sy, 'SXY=', sxy, 'SX_2Y=', sx_2y

        matrix = reshape([m, sx, sx_2, &
                          sx, sx_2, sx_3, &
                          sx_2, sx_3, sx_4], [3, 3])

        rhs = [sy, sxy, sx_2y]

        allocate(a(3))
        call solve_linear_system(matrix, rhs, a)
        print *, a

        allocate(P_x(n), e_i(n))
        do i =1,n
            P_x(i) = a(1)+a(2)*x(i)+a(3)*x(i)**2
            e_i(i) = P_x(i) - y(i)
        end do
        print *, 'X values:', x
        print *, 'Y values: ', y
        print *, 'P_x values:', P_x
        print *, 'e_i values:', e_i
        call calc_determination(y, P_x, n)
        call calc_deviation(y, P_x, n)

        print *, ""
    end subroutine
end module

