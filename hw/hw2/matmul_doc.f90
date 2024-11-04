!>===============================================================
!>-----------SUBROUTINE FOR MATRIX MULTIPLICATION----------------
!>---------------------------------------------------------------
!>      This subroutine implements matrix multiplication
!>      in three different fashions: row-major order,
!>      column-major order and using the intrinsic 
!>      function matmul implemented in Fortran.
!>  -------------------------------------------------------------
!>
!>  matrix_mult(m1, m2, m3, mode)
!>
!>
!>          Inputs | m1, m2, m3 (non-negative integer) Dimensions
!>                  of the matrices
!>                 | mode (string) Mode of matrix multiplication
!>
!>         Outputs | C (real) Result of the matrix multiplication
!>
!>===============================================================
subroutine matrix_mult(m1, m2, m3, mode)
    use debugger !< debugger module

    implicit none
 
    integer, intent(in) :: m1, m2, m3 !< dimensions of the matrices
    real, dimension(m1,m2) :: A !< matrix A
    real, dimension(m2, m3) :: B !< matrix B
    real, dimension(m1, m3) :: C !< matrix C (result)
    character(len=*), intent(in) :: mode !< mode of matrix multiplication (row, column, intrinsic)
    integer :: i, j, k !< loop indices
    real :: start_time, end_time !< time variables for measuring the time of the matrix multiplication

    !> Check for the dimensions of the matrices

    if (m1 <= 0 .or. m2 <= 0 .or. m3 <= 0) then !< check if the dimensions are valid
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Invalid dimensions')
        return
    end if

    if (size(A, 1) /= m1 .or. size(A, 2) /= m2) then !< check if the dimensions of A are valid
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Invalid dimensions of A')
        return
    end if

    if (size(B, 1) /= m2 .or. size(B, 2) /= m3) then !< check if the dimensions of B are valid
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Invalid dimensions of B')
        return
    end if

    if (mode /= "row" .and. mode /= "column" .and. mode /= "intrinsic") then !< check if the mode is valid
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Invalid mode')
        return
    end if


    !> Initialize A and B
    !> One could also use a random number generator and seed to initialize the matrices.
    do i = 1, m1
        do j = 1, m2
            A(i,j) = 1.0
        end do
    end do

    do i = 1, m2
        do j = 1, m3
            B(i,j) = 1.0
        end do
    end do

    !> Initialize C to zero
    do i = 1, m1
        do j = 1, m3
            C(i,j) = 0.0
        end do
    end do
    
    !> Check the mode of matrix multiplication

    if (mode == "row") then !< row-major order mode
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Row-major order')

            !> Measure the time of the method (in seconds)
            call cpu_time(start_time)
        
            !> this loop access the elements in a row-major order, but Fortran stores 
            !> the elements in a column-major order (differently from C++, for example). This means that this order of 
            !> execution is not cache-friendly, and the performance will be worse than the second method.
            do i = 1, m1
                do j = 1, m3
                    do k = 1, m2
                        C(i,j) = C(i,j) + A(i,k) * B(k,j)
                    end do
                end do
            end do
        
            call cpu_time(end_time)

        !> Print the result and the time of the first method
        print *, "Row-major order method result:", C(1,1), C(m1,m3)
    
        print *, "Time of the row-major order method:", real(end_time - start_time)

    
    else if (mode == "column") then !< column-major order mode
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Column-major order')

            !> Measure the time of the second method (in seconds)
            call cpu_time(start_time)
            
            !> this loop access the elements in a column-major order, which is the same order in which Fortran
            !> stores the elements. This means that this order of execution is cache-friendly, and the performance will be better.
            do i = 1, m1
                do k = 1, m2
                    do j = 1, m3
                        C(i,j) = C(i,j) + A(i,k) * B(k,j)
                    end do
                end do
            end do
        
            call cpu_time(end_time)
        
        !> Print the result and the time of the second method
        print *, "Column-major order method result:", C(1,1), C(m1,m3)
    
        print *, "Time of the column-major order method:", real(end_time - start_time)

    else if (mode == "intrinsic") then
        call checkpoint(debug = .TRUE., verb = 2, msg = 'Intrinsic function')

        !> Implent the matrix multiplication using the intrinsic function matmul
        call cpu_time(start_time)

        C = matmul(A, B) !< matrix multiplication using the intrinsic function matmul implement in Fortran 
    
        call cpu_time(end_time)
    
        !> Print the result and the time of the second method
        print *, "Intrinsic function result:", C(1,1), C(m1,m3)
    
        print *, "Time using the intrinsic function :", real(end_time - start_time)

    end if

    !> Check the result of the matrix multiplication through the intrinsic function matmul

    if (all(C == matmul(A, B))) then !< check if the result is correct (Note: you need to use the 
                                     !> all() function to compare arrays)
                                     !> Note: this is a viable check for the correcteness of the result
                                     !> for small matrices, but for large matrices, one should use a tolerance.

        call checkpoint(debug = .TRUE., verb = 2, msg = 'The result is correct')
    else
        call checkpoint(debug = .TRUE., verb = 2, msg = 'The result is incorrect')
    end if

end subroutine matrix_mult


!>===============================================================
!>---------------------MAIN PROGRAM------------------------------
!>---------------------------------------------------------------
!>      This program tests the checkpoint function and the
!>      matrix multiplication subroutine.
!>===============================================================
program main
    use debugger

    implicit none

    !> test messages
    call checkpoint(debug = .TRUE., verb = 2, msg = 'Debug mode is on')

    !> this is not executed because DEBUG is false
    call checkpoint(debug = .FALSE., verb = 2, msg = 'Debug mode is off')

    !> test matrix multiplication
    call matrix_mult(100, 100, 100, "row")
    call matrix_mult(100, 100, 100, "column")
    call matrix_mult(100, 100, 100, "intrinsic")

end program

