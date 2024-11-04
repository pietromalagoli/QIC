program main
    use debugger

    implicit none

    ! test messages
    call checkpoint(debug = .TRUE., verb = 1, msg = 'Debug mode is on')
    ! this is not executed because DEBUG is false
    call checkpoint(debug = .FALSE., verb = 2, msg = 'your message')

   end program