module matrix_timing

        implicit none

    ! Write the matri multiplication loop in two different orders 
    contains
        subroutine timing(n)

            integer :: n
            integer :: i, j, k
            real :: start_time, end_time, elapsed1, elapsed2, elapsed3
            real, allocatable :: A(:,:), B(:,:), C(:,:)     ! Declare the arrays A, B, and C as allocatable, otherwise the compiler will 
                                                            ! not know the size of the arrays at compile time since n is given as input                                    
                                                            ! This is necessary to have the value of n given as input from command line

            ! Allocate arrays after reading n
            allocate(A(n,n), B(n,n), C(n,n))
            
            ! Initialize A and B
            do i = 1, n
                do j = 1, n
                    A(i,j) = 1.0
                    B(i,j) = 2.0
                end do
            end do

            ! Initialize C to zero
            do i = 1, n
                do j = 1, n
                    C(i,j) = 0.0
                end do
            end do

            ! Measure the time of the first method (in seconds)
            call cpu_time(start_time)

            ! First method; this loop access the elements in a row-major order, but Fortran stores 
            ! the elements in a column-major order (differently from C++, for example). This means that this order of 
            ! execution is not cache-friendly, and the performance will be worse than the second method.
            do i = 1, n
                do j = 1, n
                    do k = 1, n
                        C(i,j) = C(i,j) + A(i,k) * B(k,j)
                    end do
                end do
            end do

            call cpu_time(end_time)

            elapsed1 = real(end_time - start_time)
            
            ! Initialize C to zero, so not to accumulate the results of the previous multiplication.
            do i = 1, n
                do j = 1, n
                    C(i,j) = 0.0
                end do
            end do

            ! Measure the time of the second method (in seconds)
            call cpu_time(start_time)
            
            ! Second method; this loop access the elements in a column-major order, which is the same order in which Fortran
            ! stores the elements. This means that this order of execution is cache-friendly, and the performance will be better.
            do i = 1, n
                do k = 1, n
                    do j = 1, n
                        C(i,j) = C(i,j) + A(i,k) * B(k,j)
                    end do
                end do
            end do

            call cpu_time(end_time)

            elapsed2 = real(end_time - start_time)

            ! Implent the matrix multiplication using the intrinsic function matmul

            ! Initialize C to zero, so not to accumulate the results of the previous multiplication.
            do i = 1, n
                do j = 1, n
                    C(i,j) = 0.0
                end do
            end do

            call cpu_time(start_time)

            C = matmul(A, B)

            call cpu_time(end_time)

            elapsed3 = real(end_time - start_time)

            ! Output the elapsed times in a structured format.
            print *, "RC", ",", n, ",", elapsed1
            print *, "CR", ",", n, ",", elapsed2
            print *, "I", ",", n, ",", elapsed3

            ! Final print to show that the process was successful
            print *, "Matrix multiplication completed successfully!"

        end subroutine timing

end module matrix_timing