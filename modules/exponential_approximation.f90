
module ExponentialApproximation
    use Parameters
    use SolveSystem
    implicit none
    
contains

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
            call calc_deviation_measure(e_i)
            call  calc_determination(y, P_x, n)
            call calc_deviation(y, P_x, n)
        endif
        print *, ""
        
    end subroutine
end module ExponentialApproximation