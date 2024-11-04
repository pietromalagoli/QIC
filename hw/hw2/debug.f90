module debugger

    !>============================================================= 
    !>  
    !>    This module implements a checkpoint function with 
    !>    levels of verbosity. 
    !>  ----------------------------------------------------------------- 
    !>  
    !>  SUBROUTINES: 
    !>     
    !>  checkpoint(debug, verb, msg) 
    !> 
    !> 
    !>           Inputs  | debug  (logical)  If true, the checkpoints  
    !>                            are printed in output         
    !>                   |    
    !>                   |  msg    (string) (optional)
    !>                   |  verb   (integer) (optional); possible values are 1,2 (default: 2)
    implicit none
    
    
contains
    subroutine checkpoint(debug, verb, msg)
        logical, intent(in) :: debug
        integer, intent(in) :: verb

        !> optional arguments
        !> character(len=*) is used to tell to the compiler that the length of the variable msg is not fixed
        !> intent(in) is used to tell that the variable msg is an input argument
        !> optional is used to tell that the variable msg is optional
        character(len=*), intent(in), optional :: msg 
        if (debug) then
            select case (verb)
                case (1)
                    print *, "Checkpoint"
                case (2)
                    print *, "Checkpoint"
                    if (present(msg)) then
                        print *, msg
                    end if
                case default
                    print *, "Checkpoint"
                    if (present(msg)) then
                        print *, msg
                    end if
            end select
        end if
    end subroutine checkpoint
    
end module debugger