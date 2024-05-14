module DrawGraphics
    implicit none
    
contains 

    subroutine draw(file)
        character(len=100) :: file
        if (trim(file) == 'tests\present.txt') then
            call execute_command_line('C:\Users\Asus\Desktop\Выч Мат\лабораторная №4\graphics\graph1.jpg')
        end if
    end subroutine
end module DrawGraphics