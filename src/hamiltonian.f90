!>generate Hamiltonian of the Hubbard model
!!the bases |up>|dn> are generated in the following sequence
!!|1>|1>,|1>|2>,...|1>|num_basis_dn>;|2>|1>,|2>|2>...|2>|num_basis_dn>;...;
!!|num_basis_up>|1>,|num_basis_up>|2>...|num_basis_up>|num_basis_dn>
!!H' = K_up X 1_dn + 1_up X K_dn + V_int

!>The algorithm to generate the Hamiltonian matrices
module md_ham
    use md_system_data
    implicit none
    public  ::  sb_gen_ham
    private ::  sb_init_ham
    private ::  sb_gen_hu
    private ::  sb_gen_hu_1d
    private ::  sb_gen_hu_2d
    private ::  sb_gen_ht
    private ::  sb_gen_ht_1d
    private ::  sb_gen_ht_2d

    contains
    !>Initialize the Hamiltonian matrices
    subroutine sb_init_ham(hu, ht_up, ht_dn, num_site, num_basis, num_basis_up, num_basis_dn)
        !subroutine arguments
        !>values of the diagonal elements of Hu of the 1d Hamiltonian
        integer(kind = 4), dimension(:), allocatable, intent(out)    ::  hu
        !>matrix elements of Ht of the 1d Hamiltonian in spin up subspace with hopping between nearest-neighbouring sites
        real(kind = dp), dimension(:, :), allocatable, intent(out)   ::  ht_up
        !>matrix elements of Ht of the 1d Hamiltonian in spin down subspace with hopping between nearest-neighbouring sites
        real(kind = dp), dimension(:, :), allocatable, intent(out)   ::  ht_dn
        !>number of sites
        integer(kind = 4), intent(in)   ::  num_site
        !>number of bases
        integer(kind = 4), intent(in)   ::  num_basis
        !>number of spin up bases
        integer(kind = 4), intent(in)   ::  num_basis_up
        !>number of spin down bases
        integer(kind = 4), intent(in)   ::  num_basis_dn

        allocate(hu(num_basis))
        allocate(ht_up(num_basis_up, num_site))
        allocate(ht_dn(num_basis_dn, num_site))
        !allocate(ht_pbc_dn(num_basis_dn))
        !allocate(ht_pbc_up(num_basis_up))

        hu(:) = 0
        ht_dn(:, :) = 0
        ht_up(:, :) = 0
        !ht_pbc_dn(:) = 0
        !ht_pbc_up(:) = 0

        return
    end subroutine sb_init_ham

    !!>Initialize the 1d Hamiltonian matrices
    !subroutine sb_init_ham_1d(hu, ht_up_1d, ht_dn_1d, num_site, num_basis, num_basis_up, num_basis_dn)
        !!subroutine arguments
        !!>values of the diagonal elements of Hu of the 1d Hamiltonian
        !integer(kind = 4), dimension(:), allocatable, intent(out)    ::  hu
        !!>matrix elements of Ht of the 1d Hamiltonian in spin up subspace with hopping between nearest-neighbouring sites
        !real(kind = dp), dimension(:, :), allocatable, intent(out)   ::  ht_up_1d
        !!>matrix elements of Ht of the 1d Hamiltonian in spin down subspace with hopping between nearest-neighbouring sites
        !real(kind = dp), dimension(:, :), allocatable, intent(out)   ::  ht_dn_1d
        !!>number of sites
        !integer(kind = 4), intent(in)   ::  num_site
        !!>number of bases
        !integer(kind = 4), intent(in)   ::  num_basis
        !!>number of spin up bases
        !integer(kind = 4), intent(in)   ::  num_basis_up
        !!>number of spin down bases
        !integer(kind = 4), intent(in)   ::  num_basis_dn

        !allocate(hu(num_basis))
        !allocate(ht_up_1d(num_basis_up, num_site))
        !allocate(ht_dn_1d(num_basis_dn, num_site))
        !!allocate(ht_pbc_dn(num_basis_dn))
        !!allocate(ht_pbc_up(num_basis_up))

        !hu(:) = 0
        !ht_dn_1d(:, :) = 0
        !ht_up_1d(:, :) = 0
        !!ht_pbc_dn(:) = 0
        !!ht_pbc_up(:) = 0

        !return
    !end subroutine sb_init_ham_1d

    !!>Initialize the 2d Hamiltonian matrices
    !subroutine sb_init_ham_2d(hu, ht_up_2d, ht_dn_2d, num_site, num_basis, num_basis_up, num_basis_dn)
        !!subroutine arguments
        !!>values of the diagonal elements of Hu of the 2d Hamiltonian
        !integer(kind = 4), dimension(:), allocatable, intent(out)    ::  hu
        !!>matrix elements of Ht of the 2d Hamiltonian in spin up subspace with hopping between nearest-neighbouring sites
        !real(kind = dp), dimension(:, :), allocatable, intent(out)   ::  ht_up_2d
        !!>matrix elements of Ht of the 2d Hamiltonian in spin down subspace with hopping between nearest-neighbouring sites
        !real(kind = dp), dimension(:, :), allocatable, intent(out)   ::  ht_dn_2d
        !!>number of sites
        !integer(kind = 4), intent(in)   ::  num_site
        !!>number of bases
        !integer(kind = 4), intent(in)   ::  num_basis
        !!>number of spin up bases
        !integer(kind = 4), intent(in)   ::  num_basis_up
        !!>number of spin down bases
        !integer(kind = 4), intent(in)   ::  num_basis_dn

        !allocate(hu(num_basis))
        !allocate(ht_up_2d(num_basis_up, num_site))
        !allocate(ht_dn_2d(num_basis_dn, num_site))
        !!allocate(ht_pbc_dn(num_basis_dn))
        !!allocate(ht_pbc_up(num_basis_up))

        !hu(:) = 0
        !ht_dn_2d(:, :) = 0
        !ht_up_2d(:, :) = 0
        !!ht_pbc_dn(:) = 0
        !!ht_pbc_up(:) = 0

        !return
    !end subroutine sb_init_ham_2d

    !>generate Ht part of the Hamiltonian of 1D Hubbard Model
    subroutine sb_gen_ht_1d()
        implicit none
        ! code
    end subroutine sb_gen_ht_1d

    !>generate Ht part of the Hamiltonian of 2D Hubbard Model
    subroutine sb_gen_ht_2d()
        implicit none
        ! code
    end subroutine sb_gen_ht_2d

    !>generate Hu part of the Hamiltonian of Hubbard Model
    !Hu = U * sum_i(n_i_up,n_i_dn)
    subroutine sb_gen_hu(num_site, basis_up, basis_dn, num_basis, num_basis_up, num_basis_dn, hu)
        !subroutine arguments
        !>number of sites
        integer(kind = 4), intent(in)   :: num_site
        !>number of bases
        integer(kind = 4), intent(in)   :: num_basis
        !>number of spin up bases
        integer(kind = 4), intent(in)   :: num_basis_up
        !>number of spin down bases
        integer(kind = 4), intent(in)   :: num_basis_dn
        !>spin up bases
        integer(kind = 4), intent(in)   :: basis_up(num_basis_up)
        !>spin down bases
        integer(kind = 4), intent(in)   :: basis_dn(num_basis_dn)
        !>Hamiltonian of electron's repuslsion
        integer(kind = 4), intent(out)  :: hu(num_basis)

        !local variables
        !>site occupied by a spin up electron
        logical ::  b_elec_up
        !>site occupied by a spin down electron
        logical ::  b_elec_dn
        !>index of spin up basis
        integer(kind = 4)   ::  n_no_up
        !>index of spin down basis
        integer(kind = 4)   ::  n_no_dn

        integer(kind = 4)   ::  i, j

        !traverse the rows of Ht
        do i = 1, num_basis
            n_no_dn = mod(i, num_basis_dn)
            n_no_up = i/num_basis_dn+1
            if(n_no_dn == 0) then
                n_no_dn = num_basis_dn
                n_no_up = n_no_up-1
            end if
            !traverse the sites of present basis
            do j = 1, num_site
                b_elec_up = BTEST(basis_up(n_no_up), j-1)
                b_elec_dn = BTEST(basis_dn(n_no_dn), j-1)
                !one site is occupied by a spin up electron and a spin down electron
                if(b_elec_up .and. b_elec_dn) then
                    hu(i) = hu(i)+1
                end if
            end do
        end do

        return
    end subroutine sb_gen_hu

    !>generate Ht part of the Hamiltonian 1D Hubbard Model
    subroutine sb_gen_ht(num_site, basis_up, basis_dn, num_basis_up, num_basis_dn, ht_up, ht_dn, &
                         ht_pbc_up, ht_pbc_dn, num_elec_up, num_elec_dn)
        !subroutine arguments
        !>number of sites
        integer(kind = 4), intent(in)   ::  num_site
        !>number of spin up bases
        integer(kind = 4), intent(in)   ::  num_basis_up
        !>number of spin up bases
        integer(kind = 4), intent(in)   ::  num_basis_dn
        !>number of spin up electrons
        integer(kind = 4), intent(in)   ::  num_elec_up
        !>number of spin down electrons
        integer(kind = 4), intent(in)   ::  num_elec_dn
        !>spin up bases
        integer(kind = 4), intent(in)   ::  basis_up(num_basis_up)
        !>spin down bases
        integer(kind = 4), intent(in)   ::  basis_dn(num_basis_dn)
        !>Ht in spin up subspace
        integer(kind = 4), intent(out)  ::  ht_up(num_basis_up, num_site)
        !>Ht in spin down subspace
        integer(kind = 4), intent(out)  ::  ht_dn(num_basis_dn, num_site)
        !>index of none-zero  element of Ht matrix of spin up subspace when periodic boundary condition is used
        !!if ht_pbc_up(j)=i, then the none-zero element of matrix Ht in spin up subspace is (i,j)
        integer(kind = 4), intent(out)  ::  ht_pbc_up(num_basis_up)
        !>index of none-zero  element of Ht matrix of spin down subspace when periodic boundary condition is used
        !!if ht_pbc_dn(j)=i, then the none-zero element of matrix Ht in spin down subspace is (i,j)
        integer(kind = 4), intent(out)  ::  ht_pbc_dn(num_basis_dn)

        !local variables
        integer(kind = 4)   ::  i, j, k, l, m, n, n_counter       !,nHtLength

        logical ::  b_sign, b_pbc_ht, b_temp
        !b_site_i: the bit on site i
        logical ::  b_site_i
        !b_site_j: the bit on site j=i+1
        logical ::  b_site_j
        !store the basis  which spin up electrons hop to the adjacent site
        integer(kind = 4)   ::  n_up_tmp
        !store the basis  which spin Dn electrons hop to the adjacent site
        integer(kind = 4)   ::  n_dn_tmp

        !spin up electron
        open(10, file = "ht_up.dat")
        open(20, file = "ht_pbc_up.dat")
        do k = 1, num_basis_up
            n_counter = 0
            b_pbc_ht = .false.
            if(num_site == 2) then            !only two sites
                n_up_tmp = basis_up(k)
                b_site_i = BTEST(n_up_tmp, 0)
                b_site_j = BTEST(n_up_tmp, 1)
                if(b_site_i) then
                    n_up_tmp = IBSET(n_up_tmp, 1)
                else
                    n_up_tmp = IBCLR(n_up_tmp, 1)
                end if
                if(b_site_j) then
                    n_up_tmp = IBSET(n_up_tmp, 0)
                else
                    n_up_tmp = IBCLR(n_up_tmp, 0)
                end if
                do l=1,num_basis_up
                    if(n_up_tmp == basis_up(l)) then
                        n_counter = n_counter + 1
                        ht_up(k, n_counter) = l
                        exit
                    end if
                end do
            !more than two sites
            else if(num_site > 2) then
                do i = 1, num_site
                    n_up_tmp = basis_up(k)
                    !boundary condition
                    if(i == num_site) then
                        !periodic boundary condition
                        if(b_pbc) then
                            j = 1
                            if(mod(num_elec_up, 2) == 0) then
                                b_pbc_ht = .true.
                            end if
                        !open boundary condition
                        else
                            exit
                        end if
                    else
                        !spin up electrons hop between adjacent sites in one dimension lattice
                        j = i+1
                    end if

                    b_site_i = BTEST(n_up_tmp, i-1)
                    b_site_j = BTEST(n_up_tmp, j-1)
                    !there's only one electron on the adjacent two sites
                    if(b_site_i /= b_site_j) then
                        !one electron hops between site i and j
                        if(b_site_i) then
                            n_up_tmp = IBSET(n_up_tmp, j-1)
                        else
                            n_up_tmp = IBCLR(n_up_tmp, j-1)
                        end if
                        if(b_site_j) then
                            n_up_tmp = IBSET(n_up_tmp, i-1)
                        else
                            n_up_tmp = IBCLR(n_up_tmp, i-1)
                        end if

                        !traverse the bases
                        do l= 1, num_basis_up
                            if(n_up_tmp == basis_up(l)) then
                                n_counter = n_counter + 1
                                !k is the column index and l is the row index of none-zero element of Ht matrix in spin up subspace
                                ht_up(k, n_counter) = l
                                if(b_pbc_ht) then
                                    ht_pbc_up(k) = l
                                end if
                                exit
                            end if
                        end do
                    end if
                end do
            else
                write(*, *) "the number of sites is wrong"
                stop
            end if

            do m = 1, num_site - 1
                write(10, "(I4)", advance = "NO") ht_up(k, m)
            end do
            write(10, "(I4)", advance = "YES") ht_up(k, num_site)
            write(20, "(I4)") ht_pbc_up(k)
        end do
        close(10)
        close(20)

        !spin down electron
        open(10, file = "ht_dn.dat")
        open(20, file = "ht_pbc_dn.dat")
        do k=1, num_basis_dn
            n_counter = 0
            b_pbc_ht = .false.
            if(num_site == 2) then            !only two sites
                n_dn_tmp = basis_dn(k)
                b_site_i = BTEST(n_dn_tmp, 0)
                b_site_j = BTEST(n_dn_tmp, 1)
                if(b_site_i) then
                    n_dn_tmp = IBSET(n_dn_tmp, 1)
                else
                    n_dn_tmp = IBCLR(n_dn_tmp, 1)
                end if
                if(b_site_j) then
                    n_dn_tmp = IBSET(n_dn_tmp, 0)
                else
                    n_dn_tmp = IBCLR(n_dn_tmp, 0)
                end if
                do l = 1, num_basis_dn
                    if(n_dn_tmp == basis_dn(l)) then
                        n_counter = n_counter + 1
                        ht_dn(k, n_counter) = l
                        exit
                    end if
                end do
            else if(num_site > 2) then        !more than two sites
                do i = 1, num_site
                    n_dn_tmp = basis_dn(k)
                    !boundary condition
                    if(i == num_site) then
                        !periodic boundary condition
                        if(b_pbc) then
                            j = 1
                            if(mod(num_elec_dn, 2) == 0) then
                                b_pbc_ht = .true.
                            end if
                        !open boundary condition
                        else
                            exit
                        end if
                    else
                        !spin down electrons hop between adjacent sites in one dimension lattice
                        j=i+1
                    end if

                    b_site_i = BTEST(n_dn_tmp, i - 1)
                    b_site_j = BTEST(n_dn_tmp, j - 1)
                    !there's only one electron on the adjacent two sites
                    if(b_site_i /= b_site_j) then
                        !one electron hops between site i and j
                        if(b_site_i) then
                            n_dn_tmp = IBSET(n_dn_tmp, j - 1)
                        else
                            n_dn_tmp = IBCLR(n_dn_tmp, j - 1)
                        end if    
                        if(b_site_j) then
                            n_dn_tmp = IBSET(n_dn_tmp, i - 1)
                        else
                            n_dn_tmp = IBCLR(n_dn_tmp, i - 1)
                        end if

                        !traverse the bases
                        do l = 1, num_basis_dn
                            if(n_dn_tmp == basis_dn(l)) then
                                n_counter = n_counter + 1
                                !k is the column index and l is the row index of none-zero element of Ht matrix in spin down subspace
                                ht_dn(k, n_counter) = l
                                if(b_pbc_ht) then
                                    ht_pbc_dn(k) = l
                                end if
                                exit
                            end if
                        end do
                    end if
                end do
            else
                write(*, *) "the number of sites is wrong"
                stop
            end if

            do m = 1, num_site - 1
                write(10, "(I4)", advance = "NO") ht_dn(k, m)
            end do
            write(10, "(I4)", advance = "YES") ht_dn(k, num_site)
            write(20, "(I4)") ht_pbc_dn(k)
        end do
        close(10)
        close(20)

        return
    end subroutine sb_gen_ht

    !!>print the Hu matrix
    !subroutine sb_print_hu()
        !implicit none
        !! code
    !end subroutine sb_print_hu

    !!>print the Hu matrix
    !subroutine sb_print_hu()
        !implicit none
        !! code
    !end subroutine sb_print_hu

    !!>print the Ht_up matrix
    !subroutine sb_print_ht_up()
        !implicit none
        !! code
    !end subroutine sb_print_ht_up

    !!>print the Ht_dn matrix
    !subroutine sb_print_ht_dn()
        !implicit none
        !! code
    !end subroutine sb_print_ht_dn

    !!>print the H matrix
    !subroutine sb_print_h()
        !implicit none
        !! code
    !end subroutine sb_print_h

end module md_ham

!>Practically generate the Hamiltonian matrices
module md_gen_ham
    use md_system_data
    use md_ham
    implicit none

    public  ::  sb_gen_ham
    private ::  sb_gen_ham_1d
    private ::  sb_gen_ham_2d

    contains
        !>generate the Hamiltonian matrices
        subroutine sb_gen_ham()
            if (n_dim == 1) then
                call sb_gen_ham_1d()
            end if

            if (n_dim == 2) then
                call sb_gen_ham_2d()
            end if

            return
        end subroutine sb_gen_ham

        !>generate the 1d Hamiltonian matrices
        subroutine sb_gen_ham_1d()
            call sb_init_ham_1d()
            call sb_gen_hu(n_num_site, n_basis_up, n_basis_dn, n_num_basis, n_num_basis_up, n_num_basis_dn, n_hu)
            call sb_gen_ht_1d(n_num_site, n_basis_up, n_basis_dn, n_num_basis_up, n_num_basis_dn, n_ht_up, n_ht_dn, &
                       n_ht_pbc_up, n_ht_pbc_dn, n_num_elec_up, n_num_elec_dn)

            open(20, file = "hu_1d.dat")
            write(20, "(I4)") n_hu(:)
            close(20)

            write(*, *) "Hamiltonian Generated"
            return
        end subroutine sb_gen_ham_1d

        !>generate the 2d Hamiltonian matrices
        subroutine sb_gen_ham_2d()
            call sb_init_ham_2d()
            call sb_gen_hu(n_num_site, n_basis_up, n_basis_dn, n_num_basis, n_num_basis_up, n_num_basis_dn, n_hu)
            call sb_gen_ht_2d(n_num_site, n_basis_up, n_basis_dn, n_num_basis_up, n_num_basis_dn, n_ht_up, n_ht_dn, &
                       n_ht_pbc_up, n_ht_pbc_dn, n_num_elec_up, n_num_elec_dn)

            open(20, file = "hu_2d.dat")
            write(20, "(I4)") n_hu(:)
            close(20)

            write(*, *) "Hamiltonian Generated"
            return
        end subroutine sb_gen_ham_2d
end module md_gen_ham
