module PowerApproximation
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

    subroutine power_approximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx = 0, sx_2 = 0, sxy = 0, sy = 0, m, a, b
        real, allocatable :: P_x(:), e_i(:)
        integer :: i, non_positive = 0
        real :: matrix(2,2)
        real :: rhs(2)

        print *, "Power approximation"

        do i=1,n
            if(y(i) <= 0 .OR. x(i) <= 0) then
                print *, y(i), x(i)
                non_positive = 1
            endif
        end do

        if(non_positive == 1) then
            print *, "ERROR: Power approxiamtion isn't possible, file mustn't have y,x less or equal zero"
        else
            m = n
            sx = sum(log(x))
            sx_2 = sum(log(x)**2)
            sy = sum(log(y))
            sxy = sum(log(x)*log(y))
            matrix = reshape([sx_2, sx, &
                              sx, m], [2, 2])

            rhs = [sxy, sy]
            call solve_linear_system(matrix, rhs, a, b)

            a = exp(a)
            allocate(P_x(n), e_i(n))
            do i=1,n
                P_x(i) = a*x(i)**(b)
                e_i(i) = P_x(i) - y(i)
            end do
            print *, 'X values:', x
            print *, 'Y values:', Y
            print *, 'P_x:', P_x
            print *, 'e_i:', e_i
            call  calc_determination(y, P_x, n)
            call calc_deviation(y, P_x, n)
        end if
    end subroutine

end module PowerApproximation