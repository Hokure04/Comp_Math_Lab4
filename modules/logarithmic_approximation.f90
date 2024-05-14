module LogarithmicApproximation
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

    subroutine logarithmic_approximation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx = 0, sx_2 = 0, sy = 0, sxy = 0, m, a, b
        real, allocatable :: P_x(:), e_i(:)
        integer :: i, non_positive
        real :: matrix(2,2)
        real :: rhs(2)

        print *, "Logarithmic approximation"

        do i=1,n
            if(x(i) <= 0) then
                print *, x(i)
                non_positive = 1
            endif
        end do

        if(non_positive == 1) then
            print *, "ERROR: Logarithmic approxiamtion isn't possible, file mustn't have x less or equal zero"
        else
            m = n
            sx = sum(log(x))
            sx_2 = sum(log(x)**2)
            sy = sum(y)
            sxy = sum(log(x) * y)

            matrix = reshape([sx_2, sx, &
                                sx, m], [2, 2])
                        
            rhs = [sxy, sy]

            call solve_linear_system(matrix, rhs, a, b)
            
            allocate(P_x(n), e_i(n))
            do i=1,n
                P_x(i) = b*log(x(i)) + a
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

end module LogarithmicApproximation