module SolveSystem
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

end module SolveSystem