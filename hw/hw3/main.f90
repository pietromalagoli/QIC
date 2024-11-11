program main 
    
    ! Using matrix_timing module
    use matrix_timing

    implicit none

    ! Variable declaration
    integer :: n

    ! Initial input: size of the square matrix 
    read(*,*) n
    
    ! Call the timing subroutine from the matrix_timing module
    call timing(n)

end program main