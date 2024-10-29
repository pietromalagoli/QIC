program sum_real

    ! You need to add this line so to not have Fortran set the type of the variables automatically (it looks at the first letter of the name)
    implicit none
    
    ! Declare variables with single and double precision
    real*4 :: single_res, single_pi, single_sqrt2
    real*8 :: double_res, double_pi, double_sqrt2

    ! Assign values of π and √2 to single and double precision variables
    single_pi = 3.1415927e+32  ! Approximate π in single precision with exponent notation
    single_sqrt2 = 1.4142136e+21  ! Approximate √2 in single precision with exponent notation
    
    double_pi = 3.141592653589793d+32  ! More precise π in double precision
    double_sqrt2 = 1.414213562373095d+21  ! More precise √2 in double precision

    ! Sum
    single_res = single_pi + single_sqrt2
    double_res = double_pi + double_sqrt2

    ! Display the results
    print *, "Sum with single precision:", single_res
    print *, "Sum with double precision:", double_res
    
end program sum_real