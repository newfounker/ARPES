!>set constants
module md_constant
    implicit none
    !>data type 
    integer(kind = 4), parameter :: dp = selected_real_kind(15, 99)
    !>set precision of PI
    real(kind = dp), parameter   :: PI = 3.141592653589793238462643383d0
end module md_constant

!>public variables of system
module md_system_data
    use md_constant, only: dp	
    implicit none
    !>number of sites
    integer(kind = 4), parameter  ::  n_num_site = 4
    !>number of electrons
    integer(kind = 4), parameter  ::  n_num_elec = 4
    !>number of spin up electrons
    integer(kind = 4), parameter  ::  n_num_elec_up = 1
    !>number of spin down electrons                            
    integer(kind = 4), parameter  ::  n_num_elec_dn = 2

    !>coefficient T, the hopping energy
    real(kind = dp), parameter  ::  r_t = 1.0d0
    !>coefficient U, the on-site Coulomb repulsion
    real(kind = dp), parameter  ::  r_u = 5.0d0   !known bug with ground state at u<0.1d0 and with Aii sum rule at u>10
    !>lattice constant
    real(kind = dp), parameter  ::  r_a = 1.0d0                       
    
    !>number of bases
    integer(kind = 4)   ::  n_num_basis    
    !>number of bases of spin up electrons
    integer(kind = 4)   ::  n_num_basis_up
    !>number of bases of spin down electrons
    integer(kind = 4)   ::  n_num_basis_dn
    !>bases of spin up electrons
    integer(kind = 4), dimension(:), allocatable  ::  n_basis_up
    !>bases of spin down electrons
    integer(kind = 4), dimension(:), allocatable  ::  n_basis_dn
    
    !>extreme eigen values
    real(kind = dp), dimension(:), allocatable  ::  r_eigenvalue
    !>extreme eigen vectors
    real(kind = dp), dimension(:, :), allocatable, target  ::  r_eigenvector

    !real(dp),pointer    ::  rpgroundstate(:)                

    !>Hu of the Hamilton which stores the diagonal elements of the Hamilton
    integer(kind = 4), dimension(:), allocatable    ::  n_hu
    !>Ht of the Hamilton in spin up subspace
    integer(kind = 4), dimension(:, :), allocatable ::  n_ht_up
    !>Ht of the Hamilton in spin down subspace
    integer(kind = 4), dimension(:, :), allocatable ::  n_ht_dn
    !>index of none-zero  element of ht matrix of spin up subspace when periodic boundary condition is on
    integer(kind = 4), dimension(:), allocatable    ::  n_ht_pbc_up
    !>index of none-zero  element of ht matrix of spin down subspace when periodic boundary condition is on
    integer(kind = 4), dimension(:), allocatable    ::  n_ht_pbc_dn
    
    !>switch for test code
    logical, parameter  ::  b_test = .false.
    !>switch for whether use the periodic boundary condition
    logical, parameter  ::  b_pbc  = .false.
    
    !>ground energy
    real(kind = dp)   ::  r_ground_energy
    !>ground state	
    real(kind = dp), dimension(:), allocatable   ::  r_ground_state
    
end module md_system_data
