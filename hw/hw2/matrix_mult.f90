program matrix_mult

    implicit none

    ! Write the matri multiplication loop in two different orders 

    integer, parameter :: n = 1000
    real, dimension(n,n) :: A, B, C
    integer :: i, j, k
    real :: start_time, end_time

    print *, "Matrix multiplication for n =", n
    
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

    ! Print the result and the time of the first method
    print *, "First method result:", C(1,1), C(n,n)

    print *, "Time of the first method:", real(end_time - start_time)

    
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

    ! Print the result and the time of the second method
    print *, "Second method result:", C(1,1), C(n,n)

    print *, "Time of the second method:", real(end_time - start_time)

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

    ! Print the result and the time of the second method
    print *, "Intrinsic function result:", C(1,1), C(n,n)

    print *, "Time using the intrinsic function :", real(end_time - start_time)


end program matrix_mult