!>time module
module md_time
    implicit none
    !>time when the program starts
    integer(kind = 4), dimension(8) :: n_time_start
    !>time when the program ends
    integer(kind = 4), dimension(8) :: n_time_end
    !>time when the program runs
    integer(kind = 4), dimension(8) :: n_time_run

    contains
    !>calculate the time that the programme runs
    subroutine sb_time(time_start, time_end, time_run)
        !subroutine arguments
        !>time program starts
        integer(kind = 4), dimension(8), intent(in)  :: time_start
        !>time program ends
        integer(kind = 4), dimension(8), intent(in)  :: time_end
        !>time program runs
        integer(kind = 4), dimension(8), intent(out) :: time_run
    
        !local variables
        integer(kind = 4) :: n_day

        !days of the month
        select case (time_start(2))
            case (1)
                n_day = 31
            case (2)
                n_day  = 28
                if (MOD(time_start(1), 4) == 0) then
                    n_day = 29
                    if ((MOD(time_start(1), 100) == 0) .and. (MOD(time_start(1), 400) /= 0)) then
                        n_day = 28
                    end if
                end if
            case (3)
                n_day = 31
            case (4)
                n_day = 30
            case (5)
                n_day = 31
            case (6)
                n_day = 30
            case (7)
                n_day = 31
            case (8)
                n_day = 31
            case (9)
                n_day = 30
            case (10)
                n_day = 31
            case (11)
                n_day = 30
            case (12)
                n_day = 31
            case default
                write(*, *) "Incorrect time."
                return
        end select
    
        time_run = time_end - time_start
        !>milliseconds program runs
        if (time_run(8) < 0) then
            time_run(8) = time_run(8)+1000
            time_run(7) = time_run(7)-1
        end if
        !>seconds program runs
        if (time_run(7) < 0) then
            time_run(7) = time_run(7)+60
            time_run(6) = time_run(6)-1
        end if
        !>minutes program runs
        if (time_run(6) < 0) then
            time_run(6) = time_run(6)+60
            time_run(5) = time_run(5)-1
        end if
        !>hours program runs
        if (time_run(5) < 0) then
            time_run(5) = time_run(5)+24
            time_run(3) = time_run(3)-1
        end if
        !>days program runs
        if (time_run(3) < 0) then
            time_run(3) = time_run(3)+n_day
            time_run(2) = time_run(2)-1
        end if
        !>monthes program runs
        if (time_run(2) < 0) then
            time_run(2) = time_run(2)+12
            time_run(1) = time_run(1)-1
        end if

        return
    end subroutine sb_time
end module md_time
