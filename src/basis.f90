!>choose bases: cd1...cdm...cdn|0>, 1<m<n, cdm stands for create operator on site m
!!site index increases from left to right. 
!!Set all spin down operators on the left and all spin up operators on the right
!<
module md_basis
    implicit none
    public  ::  sb_gen_basis
    private ::  sb_init_basis
    private ::  sb_gen_basis_spin
    private ::  sb_cal_num_basis

    contains
    !>genetate the bases of model which is occupied by n_num_elec_up/n_num_elec_dn spin up/down electrons
    !>the full bases {|up>|dn>} can be noted with b = b_up * 2^n_num_site + b_down 
    !!while {b_up} is stored as n_basis_up(:) and {b_down} is stored as n_basis_dn(:)
    subroutine sb_gen_basis()
        use md_system_data
        write(*,*)"no interface"
        call sb_init_basis(n_num_site, n_num_elec, n_num_elec_up, n_num_elec_dn, n_num_basis, &
                           n_num_basis_up, n_num_basis_dn, n_basis_up, n_basis_dn)
        call sb_gen_basis_spin(n_num_site, n_num_elec_up, n_num_basis_up, n_basis_up)
        call sb_gen_basis_spin(n_num_site, n_num_elec_dn, n_num_basis_dn, n_basis_dn)
        call sb_gen_full_basis(n_num_site, n_num_basis_up, n_num_basis_dn, n_basis_up, n_basis_dn)

        open(10, file = "basis_down.dat")
        write(10, "(B20)") n_basis_dn(:)
        close(10)
    
        open(10, file = "basis_up.dat")
        write(10, "(B20)") n_basis_up(:)
        close(10)
        write(*, *) "Bases generated"
    
        return
    end subroutine sb_gen_basis
    
    !>initialize the bases of a n_num_site sites Hubbard model 
    subroutine sb_init_basis(num_site, num_elec, num_elec_up, num_elec_dn, num_basis, &
                            num_basis_up, num_basis_dn, basis_up, basis_dn)
        !subroutine arguments
        !>number of sites
        integer(kind = 4), intent(in)   ::  num_site 
        !>number of electrons
        integer(kind = 4), intent(in)   ::  num_elec
        !>number of spin up electrons
        integer(kind = 4), intent(in)   ::  num_elec_up 
        !>number of spin down electrons
        integer(kind = 4), intent(in)   ::  num_elec_dn
        !>number of bases
        integer(kind = 4), intent(out)  ::  num_basis
        !>number of spin up bases
        integer(kind = 4), intent(out)  ::  num_basis_up
        !>number of spin down bases
        integer(kind = 4), intent(out)  ::  num_basis_dn
        !>spin up bases
        integer(kind = 4), dimension(:), allocatable, intent(out) :: basis_up
        !>spin down bases
        integer(kind = 4), dimension(:), allocatable, intent(out) :: basis_dn

        call sb_cal_num_basis(num_site, num_elec, num_elec_up, num_elec_dn, num_basis, &
                              num_basis_up, num_basis_dn)
        
        allocate(basis_up(num_basis_up))
        allocate(basis_dn(num_basis_dn))
        
        basis_up = 0
        basis_dn = 0
    
        return
    end subroutine sb_init_basis

    !>calculate the number of bases
    subroutine sb_cal_num_basis(num_site, num_elec, num_elec_up, num_elec_dn, num_basis, num_basis_up, num_basis_dn)
        !subroutine arguments
        !>number of sites
        integer(kind = 4), intent(in)   ::  num_site 
        !>number of electrons
        integer(kind = 4), intent(in)   ::  num_elec
        !>number of spin up electrons
        integer(kind = 4), intent(in)   ::  num_elec_up 
        !>number of spin down electrons
        integer(kind = 4), intent(in)   ::  num_elec_dn
        !>number of bases
        integer(kind = 4), intent(out)  ::  num_basis
        !>number of spin up bases
        integer(kind = 4), intent(out)  ::  num_basis_up
        !>number of spin down bases
        integer(kind = 4), intent(out)  ::  num_basis_dn
        
        !local variables
        integer(kind = 4)   ::  i, numerator, denominator
        
        !calculate the number of bases of up down electron
        numerator = 1
        denominator = 1
        do i = num_site, num_site-num_elec_up+1, -1
            numerator = numerator*i
        end do
        do i = num_elec_up, 1, -1
            denominator = denominator*i
        end do
        num_basis_up = numerator/denominator
        write(*, *) "Dimension of spin up Hilbert subspace is", num_basis_up
        
        !calculate the number of bases of spin down electron
        numerator = 1
        denominator = 1
        do i = num_site, num_site-num_elec_dn+1, -1
            numerator = numerator*i
        end do
        do i = num_elec_dn, 1, -1
            denominator = denominator*i
        end do
        num_basis_dn = numerator/denominator
        write(*, *) "Dimension of spin down Hilbert subspace is", num_basis_dn
    
        num_basis = num_basis_up*num_basis_dn
        write(*, *) "Dimension of Hilbert space is", num_basis
    
        return
    end subroutine sb_cal_num_basis
    
    !>genetate the basis of num_site sites lattice 
    !!which is occupied by num_elec_up/n_num_elec_dn spin up/down electrons
    recursive subroutine sb_gen_basis_spin(num_site, num_elec_spin, num_basis_spin, basis_spin)
        !subroutine arguments
        !>number of sites
        integer(kind = 4), intent(in)   ::  num_site 
        !>number of electrons of specific spin
        integer(kind = 4), intent(in)   ::  num_elec_spin
        !>bases of specific spin
        integer(kind = 4), intent(out)  ::  basis_spin(2**num_site)
        !>number of bases of specific spin
        integer(kind = 4), intent(inout)  ::  num_basis_spin
        
        !local variables
        !>number of bases when one site is removed and one electron of specific spin is removed
        integer(kind = 4)   ::  n_num_basis_spin_remv
        !>number of bases when one site is removed but electrons are kept
        integer(kind = 4)   ::  n_num_basis_spin_keep
        !>bases where one site and one electron removed
        integer(kind = 4)   ::  n_basis_spin_remv(2**(num_site-1))
        !>bases where one site removed but electrons kept
        integer(kind = 4)   ::  n_basis_spin_keep(2**(num_site-1))
        
        integer(kind = 4)   ::  i
        
        !IBSET: creation operator
        !IBCLR: annihilation operator
        if(num_site == 1) then
            !site is occupied by an electron
            if(num_elec_spin == 1) then
                basis_spin(1) = IBSET(basis_spin(1), 0)
                num_basis_spin = 1
            !site is not occupied by an electron
            else if(num_elec_spin == 0) then
                basis_spin(1) = IBCLR(basis_spin(1), 0)
                num_basis_spin = 1
            else
                write(*, *) "The number of electrons is more than number of sites or negtive"
                stop
            end	if
        else if(num_elec_spin == num_site) then
            basis_spin(1) = 2**num_site-1
            num_basis_spin = 1
        else if(num_elec_spin == 0) then
            basis_spin(1) = 0
            num_basis_spin = 1
        else
            !last site is occupied
            call sb_gen_basis_spin(num_site-1, num_elec_spin-1, n_num_basis_spin_remv, n_basis_spin_remv)        
            !last site is vacant
            call sb_gen_basis_spin(num_site-1, num_elec_spin, n_num_basis_spin_keep, n_basis_spin_keep)
            
            num_basis_spin = n_num_basis_spin_keep+n_num_basis_spin_remv
            do i = 1, n_num_basis_spin_remv
                basis_spin(i) = n_basis_spin_remv(i)
                basis_spin(i) = IBSET(basis_spin(i), num_site-1)
            end do
            do i = n_num_basis_spin_remv+1, num_basis_spin
                basis_spin(i) = n_basis_spin_keep(i-n_num_basis_spin_remv)
                basis_spin(i) = IBCLR(basis_spin(i), num_site-1)
            end do
        end if
        do i=1, num_basis_spin
            basis_spin(i)=IAND(basis_spin(i), 2**num_site-1)
        end do

        return
    end subroutine sb_gen_basis_spin
    
    !>generate the full bases of the n_num_site Hubbard model
    subroutine sb_gen_full_basis(num_site, num_basis_up, num_basis_dn, basis_up, basis_dn)
        implicit none
        !>number of sites
        integer(kind = 4), intent(in)   ::  num_site 
        !>number of spin up bases
        integer(kind = 4), intent(in)  ::  num_basis_up
        !>number of spin down bases
        integer(kind = 4), intent(in)  ::  num_basis_dn
        !>spin up bases
        integer(kind = 4), dimension(num_basis_up), intent(in) :: basis_up
        !>spin down bases
        integer(kind = 4), dimension(num_basis_dn), intent(in) :: basis_dn
        
        !local variables
        integer(kind = 4) :: n_full_basis
        integer(kind = 4) :: i, j, k
        k = 2**num_site
        open(10, file = "basis.dat")
        do i = 1, num_basis_up
            do j = 1, num_basis_dn
                n_full_basis = basis_up(i)*k+basis_dn(j)
                write(10, "(B20)") n_full_basis
            end do
        end do
        close(10)
        return
    end subroutine sb_gen_full_basis
end module md_basis
