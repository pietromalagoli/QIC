module c_mat
!>===============================================================
!>---------------------------------------------------------------
!>  This module defines a new type, designed to handle
!>  matrices. It contains the necessaries variables,
!>  the methods for retriving the trace and the adjoint
!>  and to initiaòize the matrix. Lastly, it writes the
!>  matrix on a file in a readable format.
!>---------------------------------------------------------------   
!> 
!> FUNCTIONS AND SUBROUTINES:   
!>  - cmat_trace: returns the trace of a square matrix
!>      Inputs: cmat: a matrix
!>      Outputs: trace: the trace of the matrix
!>
!>  - c_mat_adj: returns the adjoint of a matrix
!>      Inputs: cmat: a matrix
!>      Outputs: adj: the adjoint of the matrix
!>
!>  - c_mat_init: initializes a matrix
!>      Inputs: cmat: a matrix
!>              size: the size of the matrix
!>
!>  - c_mat_write: writes a matrix to a file
!>      Inputs: cmat: a matrix
!>              filename: the name of the file
!>
!>===============================================================

type c_matrix
    !> store the matrix dimension in a double integer
    integer, dimension(2) :: size

    !> store the matrix elements in 8 byte complex numbers
    complex*8, dimension(:,:), allocatable :: elem !< the allocatable attribute is used to declare arrays
                                                    !< whose size can be determined during runtime rather than at compile time
end type

    interface operator(.Tr.)
        module procedure cmat_trace
    end interface

    interface operator(.Adj.)
        module procedure cmat_adj
    end interface

contains

function cmat_trace(cmat) result(trace)
        type(c_matrix), intent(in) :: cmat
        complex*8 :: trace
        integer :: i

        trace = (0d0,0d0) !< initialize the trace to zero
        do i = 1, cmat%size(1)
            trace = trace + cmat%elem(i,i)
        end do

    end function 

function cmat_adj(cmat) result(adj)
        type(c_matrix), intent(in) :: cmat
        type(c_matrix) :: adj

        adj%size(1) = cmat%size(2);     adj%size(2) = cmat%size(1);
        allocate( adj%elem(adj%size(1),adj%size(2)) )
        adj%elem = conjg(transpose(cmat%elem))

    end function 

subroutine c_mat_init(cmat, size) 
        type(c_matrix), intent(inout) :: cmat
        integer, dimension(2), intent(in) :: size

        cmat%size = size
        allocate(cmat%elem(cmat%size(1),cmat%size(2)))

    end subroutine c_mat_init

subroutine c_mat_write(cmat, filename)
        type(c_matrix), intent(in) :: cmat
        character(len=*), intent(in) :: filename
        integer :: i, j
        integer :: unit

        !> Open the file for writing
        open(newunit=unit, file=filename, status='replace', action='write')

        ! Write the matrix dimensions
        write(unit, '(A, 2I10)') 'Matrix dimensions: ', cmat%size

        !> Write the matrix elements
        !> Note: particular attention must be paid to the fact that this 
        !> custom matrix type is designed to store 8-byte complex numbers.
        write(unit, '(A)') 'Matrix elements:'
        do i = 1, cmat%size(1)
            do j = 1, cmat%size(2)
                write(unit, '(F10.4, " + i*", F10.4)', advance='no') real(cmat%elem(i, j)), aimag(cmat%elem(i, j))

                  ! Add a comma between elements on the same row, for readibility
                if (j < cmat%size(2)) then
                    write(unit, '(A)', advance='no') ', '
                end if
            end do

        ! New line after each row, for readibility
            write(unit, *) '' 
        end do

        !> Close the file
        close(unit)
    end subroutine c_mat_write

end module c_mat