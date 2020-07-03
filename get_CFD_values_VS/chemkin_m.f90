module chemkin_params
    ! unit nuber for input
    integer, parameter :: unit_tplink = 20
    integer, parameter :: unit_cklink = 21
    integer, parameter :: unit_stdout = 6
    
    ! length of work array
    integer len_logi_ckwk
    integer len_int_ckwk
    integer len_real_ckwk
    integer len_char_ckwk

    ! chemkin work array
    integer,      allocatable :: int_ckwk(:)
    real(8),      allocatable :: real_ckwk(:)
    character(6), allocatable :: char_ckwk(:)*16
    
    ! transport work array
    integer,      allocatable :: int_tpwk(:)
    real(8),      allocatable :: real_tpwk(:)

    ! index of chemistry
    integer mm ! number of elements
    integer kk ! number of species
    integer ii ! number of reactions
    integer nift ! number of fitting coefficients
    
    ! array pointers
    integer nsys, neq, lsdas, lidas, lrdas

    contains

    subroutine initialize_chemkin_workarray()

        ! open input unit
        open(unit_cklink, form='unformatted', file='./link/cklink')
        open(unit_tplink, form='unformatted', file='./link/tplink')


        !   ------- initialize chemkin work array ---------
        
        call cklen(unit_cklink, unit_stdout, len_int_ckwk, len_real_ckwk, len_char_ckwk)
        
        allocate(int_ckwk(len_int_ckwk))
        allocate(real_ckwk(len_real_ckwk))
        allocate(char_ckwk(len_char_ckwk))
        
        call ckinit(len_int_ckwk, len_real_ckwk, len_char_ckwk, unit_cklink, &
                    unit_stdout, int_ckwk, real_ckwk, char_ckwk)

        !   ------- get chemkin index ---------
        
        call ckindx(int_ckwk, real_ckwk, mm, kk, ii, nfit)
        nsys = kk + 1
        neq = nsys
        lsdas = 1
        lidas = 20 + neq
        lrdas = 40 + neq*9 + nsys*nsys

        !   ------- initialize transport work array ---------

        call mclen(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk)
        
        allocate(int_tpwk(len_int_tpwk))
        allocate(real_tpwk(len_real_tpwk))

        call mcinit(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk, &
                      int_tpwk, real_tpwk)

    end subroutine initialize_chemkin_workarray

    subroutine get_tranport_data(t_cfd, p_cfd, y_cfd, D_mix, Lambda_mix, c_p)

        ! input values
        real(8), intent(in) :: t_cfd     ! K
        real(8), intent(in) :: p_cfd     ! Pa
        real(8), intent(in) :: y_cfd(kk) ! Mass fractions

        ! output transport data
        ! mixture diffusion coefficient [CM**2/S]
        real(8), intent(out) :: D_mix(kk) 
        ! mixture thermal conductivity [ERG/CM*K*S]
        real(8), intent(out) :: Lambda_mix
        !  mean specific heat at constant pressure [ergs/(gm*K)]
        real(8), intent(out) :: c_p

        ! variables for calculations
        real(8) p_calc ! dyne/cm**2
        real(8) :: x_calc(kk) ! Mole fractions
        p_calc = p_cfd*10.0d0       ! Pa to dyne/cm**2
        call ckytx(y_cfd, int_ckwk, real_ckwk, x_calc)

        call mcadif(p_calc, t_cfd, x_calc, real_tpwk, D_mix) 
        call mcacon (t_cfd, x_calc, real_tpwk, Lambda_mix)
        call ckcpbs(t_cfd, y_cfd, int_ckwk, real_ckwk, c_p)

    end subroutine get_tranport_data

    subroutine get_next_TY(p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

        real(8), intent(inout) :: t_cfd
        real(8), intent(in)    :: p_cfd
        real(8), intent(inout) :: y_cfd(kk)
        real(8), intent(in)    :: delta_t_cfd
        real(8), intent(in)    :: tols_cfd(4)

        integer :: icase = 1
        integer  lin, lout, lsave, lign, lrest
        logical :: lsens = .false.

        DATA LIN/5/, LOUT/6/, LINKCK/25/, LSAVE/7/, LIGN/9/, LREST/10/

        !call begin(nsys, neq, icase, ii, kk, len_int_ckwk, len_real_ckwk,      &
        !    len_char_ckwk, unit_cklink, lin, lout, lsave, lign, lrest, lsens,   &
        !    lidas, lrdas, lsdas, int_ckwk(nidas), real_ckwk(nrdas),              &
        !    real_ckwk(nsdas), real_ckwk(nrpar), int_ckwk(nipar), real_ckwk(nz),  &
        !    real_ckwk(nzp), real_ckwk(nrtol), real_ckwk(natol), real_ckwk(nxmol),&
        !    char_ckwk(nksym), char_ckwk(ipcck),                                  &
        !    p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

    end subroutine get_next_TY

end module chemkin_params