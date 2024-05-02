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

program main
    use FileReader
    use KeyboardReader
    use LinearAproximation
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
            call linear_approximation(x,y,n)
        else if (func_number == 2) then
            call keyboard_input(x, y, n)
            print *, "x values:", x
            print *, "y values:", y
            print *, "n:", n 
        end if
    end do

end program
