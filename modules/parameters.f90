module Parameters
    implicit none

contains
    subroutine calc_correlation(x, y, n)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in) :: n
        real :: r, mean_x, mean_y, sum_x, sum_y, numerator
        
        mean_x = sum(x)/n
        mean_y = sum(y)/n
        sum_x = sum((x-mean_x)**2)
        sum_y = sum((y-mean_y)**2)
        numerator = sum((x-mean_x)*(y-mean_y))
        r = numerator/sqrt(sum_x*sum_y)
        print *, 'Correlation coefficient: ', r
    end subroutine

    subroutine calc_determination(y, fi, n)
        real, allocatable, intent(in) :: y(:), fi(:)
        integer, intent(in) :: n
        real :: mean_fi, numerator, denominator, determination

        mean_fi = sum(fi)/n
        
        numerator = sum((y-fi)**2)
        denominator = sum((y-mean_fi)**2)
        determination = 1 - numerator/denominator
        print *,'Determination coefficient: ', determination
    end subroutine

    subroutine calc_deviation(y, fi, n)
        real, allocatable, intent(in) :: y(:), fi(:)
        integer, intent(in) :: n
        real :: sum, deviation
        
        deviation = sqrt(sum((fi-y)**2)/n)
        print *, 'Standart deviation: ', deviation
    end subroutine


end module Parameters