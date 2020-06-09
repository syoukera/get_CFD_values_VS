!     test code for using CHEMKIN library

program test_main
      use chemkin_params, only: initialize_chemkin_workarray, &
                                get_tranport_data, get_next_TY
      use output, only: make_output

      !   ------- user input section ---------

      integer,parameter :: num_spec = 53
      
      ! phisycal values from CFD calculation
      real(8) :: p_cfd = 1.01325d5*40   ! Pa
      real(8) :: t_cfd = 1000d0         ! K
      real(8) y_cfd(num_spec)           ! Mass fractions
      real(8) :: delta_t_cfd = 1.0d0    ! s
      real(8) :: tols_cfd(4)            ! Tolerances
      
      ! output transport data
      ! mixture diffusion coefficient [CM**2/S]
      real(8) :: D_mix(num_spec) 
      ! mixture thermal conductivity [ERG/CM*K*S]
      real(8) :: Lambda_mix
      !  mean specific heat at constant pressure [ergs/(gm*K)]
      real(8) c_p
      
      ! Assurme Y has a same secuence as species in ckout
      data y_cfd /0.00d+00,0.00d+00,0.00d+00,2.20d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,5.51d-02,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,7.24d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00/

      ! tolerance values for ODE solver
      data tols_cfd /1.d-8, 1.d-20, 1.d-5, 1.d-5/

      ! flag to generate output
      make_output = .false.

      !   ------- initialize section ---------

      call initialize_chemkin_workarray()
      
      !   ------- transport section ---------

      call get_tranport_data(t_cfd, p_cfd, y_cfd, D_mix, Lambda_mix, c_p)

      write(6, *) 'mixture diffusion coefficient [CM**2/S]'
      write(6, *) D_mix
      write(6, *) 'mixture thermal conductivity [ERG/CM*K*S]'
      write(6, *) Lambda_mix
      write(6, *) 'mean specific heat at constant pressure [ergs/(gm*K)]'
      write(6, *) c_p
      
      !   ------- chemistry section ---------

      call get_next_TY(p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

      write(6, *) 'temperature [K]'
      write(6, *) t_cfd
      write(6, *) 'mole fractions [-]'
      write(6, *) y_cfd

end program test_main