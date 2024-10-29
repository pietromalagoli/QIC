program sum_integer_types

    ! You need to add this line so to not have Fortran set the type of the variables automatically (it looks at the first letter of the name)
    implicit none
    
    ! Declare the variable via the protocol type :: name_var
    integer*2 :: small_res
    integer*4 :: large_res
    
    ! Do the sum
    small_res = 2000000 + 1     ! Should overflow
    large_res = 2000000 + 1     ! Should not overflow
    
    ! Display the results
    print *, "Sum with integer*2 (overflow expected):", small_res
    print *, "Sum with integer*4:", large_res
    
end program sum_integer_types