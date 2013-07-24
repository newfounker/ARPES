program main
    use md_time
    use md_system_data
    use md_basis
    !use md_hamiltonian
    
    implicit none
    call date_and_time(values = n_time_start)
    write(*, "(A24, I4, A1, I2.2, A1, I2.2, 1X, A3, I3.3, 1X, 3(I2.2, A1), I3.3)") &
        "The programme starts at ", n_time_start(1), "/", n_time_start(2), "/", n_time_start(3), "UTC", n_time_start(4), &
        n_time_start(5), ":", n_time_start(6), ":", n_time_start(7), ".", n_time_start(8)

    !>switch for the test subroutines
    !logical :: b_Test   

    call sb_gen_basis()
    !call sb_Gen_Ham()
    !call Lanczos()
    
    !>Do some tests
    !!to verify the subroutines above
!    b_test = .false.
!    if(b_test) then
!        call LWtest()
!        call testAii()
!        call testSumRule()
!    end if

!    call InitSpec()
!    call Spectrum()
    write(*,*) 'n_basis_dn'
    write(*,"(B20)") n_basis_dn
    write(*,*) 'n_basis_up'
    write(*,"(B20)") n_basis_up
    call date_and_time(values = n_time_end)
    call sb_time(n_time_start, n_time_end, n_time_run)
    write(*, "(A22, I4, A1, I2.2, A1, I2.2, 1X, A3, I3.3, 1X, 3(I2.2, A1), I3.3)") &
        "The programme ends at ", n_time_end(1), "/", n_time_end(2), "/", n_time_end(3), "UTC", n_time_end(4), &
        n_time_end(5), ":", n_time_end(6), ":", n_time_end(7), ".", n_time_end(8)
    write(*, "(A18, I4.2, A1, I2.2, A1, I2.2, 1X, 3(I2.2, A1), I3.3)") &
        "The programme runs ", n_time_run(1), "/", n_time_run(2), "/", n_time_run(3), &
        n_time_run(5), ":", n_time_run(6), ":", n_time_run(7), ".", n_time_run(8)
    stop
end program main
