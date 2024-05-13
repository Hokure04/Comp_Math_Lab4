module LinearAproximation
    use Parameters
    implicit none
    
contains
    subroutine linear_approximation(x,y,n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: sx, sxx, sy, sxy, delta, delta1, delta2, a, b
        real, allocatable :: P_x(:), e_i(:)
        integer :: i
        sx = sum(x)
        sxx = sum(x**2)
        sy = sum(y)
        sxy = sum(x*y)
        
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