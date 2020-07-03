module senkin_params
    implicit none
    
    integer, parameter :: unit_skout = 20
    
    real(8) :: tols(4) = (/1.d-8, 1.d-20, 1.d-5, 1.d-5/)
    
    real(8), parameter :: pressure_Pa = 1.01325d5*40 ! Pa
    real(8), parameter :: pressure = pressure_Pa*10d0 ! Dyne/cm**2

    real(8), allocatable :: z(:)
    real(8), allocatable :: zp(:)
    real(8), allocatable :: atol(:)
    real(8), allocatable :: rtol(:)
    
    real(8), allocatable :: dwork(:)
    integer, allocatable :: idwork(:)
    real(8), allocatable :: sdwork(:)
    
    integer :: isen(5) = 0
    integer :: info(15) = 0
            
    real(8) :: time = 0.0d0
    real(8), parameter :: end_time = 1.0d0
    real(8), parameter :: dtime = 1.0d-2

    contains
    
        subroutine allocate_arrays()
            use chemkin_params, only: nsys, lrdas, lidas, lsdas
            implicit none
            
            allocate(z(nsys))
            allocate(zp(nsys))
            
            allocate(atol(nsys))
            allocate(rtol(nsys))
            
            allocate(dwork(lrdas))
            allocate(idwork(lidas))
            allocate(sdwork(lsdas))
            
        end subroutine allocate_arrays
        
        subroutine init_variables()
            implicit none
            
            open(unit_skout, form='formatted', file='skout.csv')
            
            z(1) = 1000d0
            z(2:) = (/0.00d+00,0.00d+00,0.00d+00,2.20d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,5.51d-02,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,7.24d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00/)
            
            rtol = tols(1)
            atol = tols(2)
            
        end subroutine init_variables

        subroutine timeloop()
        
            implicit none
            
            call allocate_arrays()
            call init_variables()
            call init_dasac()

            do while (time < end_time)
                                
                call timestep_dasac()
                call output_result()
                
            end do

        end subroutine timeloop
        
        subroutine init_dasac()
            use chemkin_params, only: nsys, int_ckwk, real_ckwk
            implicit none
            
            integer :: dummy = 1 ! dummy argument
            
            external dres
            
            info(1) = 0
            info(2) = 1
            info(3) = 1
            
            call ddinit(nsys, isen(5), time, real_ckwk, int_ckwk, &
                        zp, z, dummy, isen(1), isen(4), rconp_re, dres)
        
        end subroutine init_dasac

        subroutine timestep_dasac()
            use chemkin_params, only: nsys, lrdas, lidas, lsdas, int_ckwk, real_ckwk
            implicit none
            
            integer :: idid
            integer :: jac    ! dummy argument
            real(8) :: dfdyp  ! dummy argument
            
            external dres
            
            call ddasac(rconp_re, nsys, time, z, zp, end_time, info, isen,   &
                        rtol, atol, idid, sdwork, lsdas, dwork, lrdas,    &
                        idwork, lidas, real_ckwk, int_ckwk, &
                        jac, dres, dfdyp)
            
            if (idid <= 0) then
                print *, 'idid = ', idid
                stop
            end if

        end subroutine timestep_dasac
        
        subroutine rconp_re(time, z, zp, delta, ires, real_ckwk, int_ckwk)
            use chemkin_params, only: nsys, kk
            implicit none
            
            real(8), intent(in) :: time
            real(8), intent(in) :: z(nsys)
            real(8), intent(in) :: zp(nsys)
            real(8), intent(out) :: delta(nsys)
            integer, intent(in) :: ires
            real(8), intent(in) :: real_ckwk(*)
            integer, intent(in) :: int_ckwk(*)
            
            real(8) :: wt(kk)
            real(8) :: wdot(kk)
            real(8) :: rho
            real(8) :: cpb
            real(8) :: hms(kk)
            real(8) :: volume
            real(8) :: sum
            integer :: k
            
            call ckwt(int_ckwk, real_ckwk, wt)
            call ckrhoy(pressure, z(1), z(2), int_ckwk, real_ckwk, rho)
            call ckcpbs(z(1), z(2), int_ckwk, real_ckwk, cpb)
            call ckwyp(pressure, z(1), z(2), int_ckwk, real_ckwk, wdot)
            call ckhms(z(1), int_ckwk, real_ckwk, hms)
            
            if (rho == 0.0d0) then
                print *, 'Stop, zero density in RCONP'
                stop
            end if
            volume = 1.0d0/rho
            
            ! energy eq.
             delta(1) = zp(1) + volume*sum(hms*wdot*wt)/cpb
             
            ! species conservation eq.
             delta(2:) = zp(2:) - wdot*wt*volume

        end subroutine rconp_re
        
        subroutine output_result()
            implicit none
            
            write (unit_skout, '(100((E8.3)","))') time, z       
        
        end subroutine output_result


end module senkin_params