program monte

  implicit none

  integer, parameter :: dp = selected_real_kind (15,307)
  real (dp), parameter :: boltz_const = 1
  real (dp), parameter  :: pi=3.1415926535897932384626433

  integer, parameter :: unit1 = 10, unit_pressures = 20, unit_potential_en = 30

  !no dimensions
  integer, parameter :: d = 3

  integer  :: i, j, k,    N, N_moves, Equil_moves, accept_count

  real (dp) :: start_time, stop_time

  real (dp) :: sigmalj, epsilonlj, volume,  L, density, r_cut, temperature, beta, max_move

  real (dp) :: lj_pot, lj_pot_diff, lj_pot_corr, lj_pot_corr_new, lj_pot_mean, lj_pot_mean_prev, lj_pot_pre_stdev, lj_pot_stdev

  real (dp) :: factor, P, P_mean, P_stdev, P_corr
  real (dp) :: L_new, pressure_factor, pressure_factor_mean, pressure_factor_mean_prev, pre_stdev


  !coord array
  real (dp), allocatable, dimension (:,:) :: r


  !path to save files to
  character*(*), parameter :: directory = 'output'
  character*(*), parameter :: path_out = './' // directory // '/'

  character*(*), parameter :: makedirectory = 'mkdir ' // directory

  ! filename/path of the input/control files

  character*(*), parameter :: directory_in = 'input'
  character*(*), parameter :: path_in = './' // directory_in // '/'

  character(len=100) :: param_file, input_coord_file
  input_coord_file = path_in//'input.txt'
  param_file = path_in//'params.txt'



  !
  !
  !       Program
  !
  !


  call read_param_file(param_file, sigmalj, epsilonlj, r_cut, max_move, N_moves, &
  temperature, Equil_moves, density)

  call init_random_seed()


  pre_stdev = 0.0_dp
  P_corr = 0.0_dp
  pressure_factor_mean_prev = 0.0_dp
  pressure_factor_mean = 0.0_dp
  P_stdev = 0.0_dp


  lj_pot_corr = 0.0_dp
  lj_pot_mean_prev = 0.0_dp
  lj_pot_pre_stdev = 0.0_dp
  lj_pot_stdev = 0.0_dp
  lj_pot_mean = 0.0_dp

  factor = 1.000005_dp

  beta = 1.0_dp / (temperature * boltz_const )


  !compile with fortran2008 standard for support for the execute_command_line
  call execute_command_line (makedirectory)

  call cpu_time (start_time)

  call read_no_particles (input_coord_file, d,N)

  allocate ( r(d,N) )

  L = (  (real (N, dp))/ (density) ) **(1.0_dp/3.0_dp)

  L_new = L * factor



  call read_in_coords (input_coord_file, r)



  !start running main code once all parameters read-in:
  call pbc ( r, L )

  call init_pot ( r, lj_pot , sigmalj, epsilonlj, d, N, L, r_cut )

  call pressure_and_lj_pot_correction (pi, sigmalj, epsilonlj, r_cut, density, P_corr, lj_pot_corr, factor, lj_pot_corr_new )




  call write_out_params (path_out, sigmalj, epsilonlj, L, r_cut, max_move, N_moves, Equil_moves, &
  temperature, N, density, r, lj_pot, lj_pot_corr)


  !set up files for data write-out
  open (unit1,file=path_out//'traj.xyz')

  open (unit_pressures,file=path_out//'pressures.txt')
  write (unit_pressures, '(1x,a,10x,a)') 'loop', 'Pressure'

  open (unit_potential_en,file=path_out//'potential.txt')
  write (unit_potential_en, '(1x,a,10x,a,10x,a,10x,a,10x,a)') 'loop', 'uncorrected potential', 'corrected potential', &
  'running mean of corrected', 'standard deviation'
  write (unit_potential_en, *) 'Before equil', lj_pot, lj_pot+lj_pot_corr



  !Equilibration
  if (Equil_moves < (2*N) ) then
    Equil_moves =  0
    write (*,*) 'Ignoring equilibration step'
  end if

  do i = 1, Equil_moves

      call equil_new_position_and_potential (r, sigmalj, epsilonlj, d, N, max_move, L, &
      beta, P, accept_count, lj_pot, r_cut )

  end do

  i = 0
  write (unit_potential_en,*) i, lj_pot, lj_pot+lj_pot_corr

  write(*,*) pressure_factor, pressure_factor_mean


  call pressure_calc ( beta, r, lj_pot, lj_pot_corr, lj_pot_corr_new, L, N, d, r_cut, sigmalj, epsilonlj, &
  pressure_factor, factor)

  write(*,*) pressure_factor, pressure_factor_mean



  !Production run
  accept_count = 0

  do i = 1, N_moves

      call new_position_and_potential (r, sigmalj, epsilonlj, d, N, max_move, L, &
      beta, P, accept_count, lj_pot, r_cut, lj_pot_corr, lj_pot_corr_new, pressure_factor, factor )



      !write out positions to trajectory file, every N, (no. atoms) loops
      if (  modulo(i, N) == 0  ) then

        call write_out_s_to_file (i, r, N, unit1 )

      end if


      !stop errors accumulating in lj_pot during typical move cycle
      if (  modulo(i, (N*10)) == 0  ) then

        call init_pot ( r, lj_pot , sigmalj, epsilonlj, d, N, L, r_cut )

      end if

      call equil_lj_pot_stats (i, lj_pot, lj_pot_mean, lj_pot_mean_prev, lj_pot_pre_stdev, lj_pot_stdev)


      write (unit_potential_en,*) i, lj_pot, lj_pot+lj_pot_corr, lj_pot_mean+lj_pot_corr, lj_pot_stdev


      pressure_factor_mean = pressure_factor_mean + (  (pressure_factor - pressure_factor_mean) / ( real (i, dp) )   )

      !Welford algorithm
      ! pressure_factor_mean = pressure_factor_mean_prev + (  (pressure_factor - pressure_factor_mean_prev) / ( real (i, dp) )   )
      ! pre_stdev = pre_stdev + (  (pressure_factor - pressure_factor_mean_prev)*(pressure_factor - pressure_factor_mean)  )
      ! pressure_factor_mean_prev = pressure_factor_mean


      P = P_corr +  (    (log (pressure_factor_mean) ) / (beta * ( (L_new**3) - (L**3) ) )   )

      ! write(*,*) pressure_factor, pressure_factor_mean, P


      write (unit_pressures,*) i, P



  end do

  write (*,'(/,a)') 'Results:'

  write (*,*) 'accepted moves', accept_count
  write (*,*) 'total moves', N_moves
  write (*,*) 'ratio', (real(accept_count))/(real(N_moves))
  write (*,*) 'average pressure', P
  write (*,*) 'pressure correction', P_corr
  write (*,*) 'lj correction', lj_pot_corr



  close (unit1)
  close (unit_pressures)
  close (unit_potential_en)

  open (unit1,file=path_out//'final_coords.xyz')

    write (unit1,*) N
    write (unit1,*) i
    do j = 1, N
        write (unit1,*) 'Ar', r (:,j)
    end do
  close (unit1)

  open (unit1,file=path_out//'final_coords.txt')
    do j = 1, N
        write (unit1,*)  r (:,j)
    end do
  close (unit1)

  deallocate ( r )

  call cpu_time(stop_time)
  write (*,*) 'Time taken:', stop_time - start_time, "seconds"



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !           Subroutines
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains


    subroutine read_param_file(param_file, sigmalj, epsilonlj, r_cut, max_move, N_moves, &
      temperature, Equil_moves, density)
      ! new variable in needs: change int/real param no., a string & assignment at end of subroutine

      implicit none

      character(len=100), intent(in) :: param_file

      ! Control file variables
      real (dp), intent (out) :: sigmalj, epsilonlj, r_cut, max_move, temperature, density
      integer, intent (out) :: N_moves, Equil_moves


      ! Input related variables
      character(len=100) :: buffer, label
      integer :: pos, i, param_no, no_int_params, no_real_params
      integer, parameter :: unit5 = 15
      integer :: ios = 0
      integer :: iosfile = 0
      integer :: line = 0

      ! Timer related variables
      real(dp) :: start_time, stop_time

      integer, allocatable, dimension (:) :: chk
      character(len=100), allocatable, dimension (:) :: real_strings
      character(len=100), allocatable, dimension (:) :: int_strings
      real (dp), allocatable, dimension (:) :: reals
      integer, allocatable, dimension (:) :: ints

      ! Parameter type/no. variables
      no_int_params = 2
      no_real_params = 6

      param_no = no_int_params + no_real_params

      allocate( chk(param_no) )
      allocate( real_strings(no_real_params), reals(no_real_params) )
      allocate( int_strings(no_int_params), ints(no_int_params) )

      chk (:) = 0


      int_strings(1) = 'N_moves'
      int_strings(2) = 'Equil_moves'

      real_strings(1) = 'sigmalj'
      real_strings(2) = 'epsilonlj'
      real_strings(3) = 'density'
      real_strings(4) = 'r_cut'
      real_strings(5) = 'max_move'
      real_strings(6) = 'temperature'


      open(unit5, file=param_file, status='old', iostat=iosfile)

      if (iosfile /= 0) then

        write (*,'(a,/,a)') 'parameter file not found','program terminating...'
        STOP

      end if


      ! start timer in case an error occurs or parameter file reading takes a very long time
      call cpu_time (start_time)
      call cpu_time (stop_time)



      !read first line into "buffer"
      read(unit5, '(A)', iostat=ios) buffer


      !do until reach end of file or find all variables (or an error occurs)
      do while (  (ios == 0)   .and.  ( sum(chk) < param_no ) .and. ( (stop_time-start_time) < 5_dp )   )

            !count line no.
            line = line + 1

            !remove blank space at the start of a line
            buffer = adjustl (buffer)

            !pos = place where whitepace starts (1 for blank lines)
            pos = index(buffer, ' ')

            !label excludes blank space and includes data's 'label'
            label = buffer(1:pos)

            !pos+1 to end -  holds the variable string and blank space
            buffer = buffer(pos+1:)

            !take no action if a line doesn't have the variable listed in front of everything and then a space
            !e.g.:
            !N    7000    - correct
            !  N   7000   - correct
            !a   7000     - incorrect

            !check line for integers
            do i = 1, no_int_params

              !if "label" matches a variable name, read value in "buffer" in
              if ( label == int_strings(i) ) then

                !get a value for the integer; if already have a value for the integer, don't read in another
                if (chk(i) == 0) then
                  read(buffer, *, iostat=ios)  ints (i)
                  if (ios == 0) then
                    chk(i) = 1
                  end if
                !otherwise set pos =1 , so that the param counter below treats the line as blank
                else
                  pos =1
                end if

              !change pos value to 1 for any non-variable lines
              else
                pos =1
              end if
            end do

            !repeat for real values
            do i = 1, no_real_params

              if ( label == real_strings(i) ) then

                if (chk(i+no_int_params) == 0) then
                  read(buffer, *, iostat=ios)  reals (i)
                  if (ios == 0) then
                    chk(i+no_int_params) = 1
                  end if
                else
                  pos =1
                end if
              else
                pos =1
              end if
            end do

            !read next line into buffer
            read(unit5, '(A)', iostat=ios) buffer
            !write (*,*) line, pos, param_check

            call cpu_time (stop_time)

      end do


      close (unit5)

      !timeout stop
      if  ( (stop_time-start_time) >= 5_dp  ) then
        write (*,'(a,/,a)') 'parameter file has taken over 5 seconds to read in', 'program terminating...'
        STOP
      end if


      !warning message if can't find/read-in all parameters
      if ( param_no /= sum(chk) ) then
        write (*,'(a,/,a,/,a,I0)') 'WARNING, could not read-in all parameters', 'Running calculation with unknown number', &
        'read parameter file down to line number  ', line
      else
        write (*,'(a,/)') 'all parameters read-in  -  see output/input_params for their values'
      end if

      !specific warning messages
      do i = 1, no_int_params
        if ( chk(i) /= 1) then
          write (*,'(a,3x,a)') 'Failed to read-in:', int_strings(i)
        end if
      end do

      do i = 1, no_real_params
        if ( chk(i+no_int_params) /= 1) then
          write (*,'(a,3x,a)') 'Failed to read-in a value for:', real_strings(i)
        end if
      end do

      ! assign the values read into the arrays to their respective variables
      !integer variables
      N_moves = ints(1)
      Equil_moves = ints(2)

      !real variables
      sigmalj = reals (1)
      epsilonlj = reals (2)
      density = reals (3)
      r_cut = reals (4)
      max_move = reals (5)
      temperature = reals (6)

    end subroutine read_param_file



    subroutine init_random_seed()

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ ((i - 1)**2, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)

    end subroutine init_random_seed



  subroutine read_no_particles (input_coord_file,d,N)

    implicit none

    character(len=100), intent(in) :: input_coord_file

    integer, intent (in) :: d
    integer, intent (out) :: N

    integer  ::  io, unit

    real (dp) :: r (d,1)


    !read-in coord file
    unit = 20

    open (unit, file=input_coord_file, status='old')
        !n = no. of lines - start at -1 beacuse of while loop structure below
        N = -1
        io = 0
        do while ( io == 0 )
            N = N + 1
            read(unit,*,iostat=io) r
        end do

    close (unit)

  end subroutine read_no_particles



  subroutine read_in_coords (input_coord_file, r)

    implicit none

    character(len=100), intent(in) :: input_coord_file

    real (dp), intent (inout) :: r (:,:)

    integer :: i,j,    d, N, unit

    d = size (r,1)
    N = size (r,2)
    unit = 20

    open (unit, file=input_coord_file, status='old')

      read (unit,*)   ( (   r(i,j),   i=1,d   ),   j=1,N  )

    close (unit)

  end subroutine read_in_coords



  subroutine init_pot ( r, lj_pot , sigmalj, epsilonlj, d, N, L, r_cut )

    real (dp), intent (in) :: r (:,:), sigmalj, epsilonlj, L, r_cut
    integer, intent (in) :: d, N
    real (dp), intent (out) :: lj_pot

    integer :: i, j, k

    real (dp), dimension (d) :: del_dim
    real (dp) :: r_mod_sq, r_cut_sq, sigma_6_div_r_mod_6, L_half


    lj_pot = 0.0_dp

    L_half = (L/2.0_dp)

    r_cut_sq = r_cut**2


    do i = 1, N-1
      do j = i+1, N

        del_dim(:) = r (:,i) - r (:,j)

            !pbc conditions
            do k = 1, d
                if (del_dim(k) > L_half ) then
                  del_dim(k) = del_dim(k) - L
                end if
                if (del_dim(k) < - L_half ) then
                  del_dim(k) = del_dim(k) + L
                end if
            end do

            r_mod_sq = sum (del_dim**2)


            !cutoff condition
            if  (r_mod_sq  < r_cut_sq)  then

                  sigma_6_div_r_mod_6 = (sigmalj**6)/(r_mod_sq**3)

                  lj_pot = lj_pot + (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - (1.0_dp)) )

            end if

        end do
    end do

    !multiply by a correction factor
    lj_pot = lj_pot * (epsilonlj * 4.0_dp )

  end subroutine init_pot



  subroutine pbc (r, L)

    real (dp), intent (inout) :: r(:, :)
    real (dp), intent (in) :: L
    real (dp) :: L_half


    integer :: n , m

    m = size(r, 1)      !m = dimension no
    n = size(r, 2)      !n = No particles

    L_half = (L/2.0_dp)

    !check all particles coords/dimensions - move back inside of box
    do i = 1, m
      do j = 1, n
              r(i,j) = r(i,j) - L * (nint( r(i,j) / L ))
      end do
    end do

  end subroutine pbc


    subroutine equil_new_position_and_potential (r, sigmalj, epsilonlj, d, N, max_move, L, &
      beta, P, accept_count, lj_pot, r_cut )

      real (dp), intent (in) :: sigmalj, epsilonlj, max_move, L, r_cut, beta
      integer, intent (in) :: d, N
      real (dp), intent (inout) :: r(:,:)
      real (dp), intent (inout) :: lj_pot, P
      integer, intent (inout) :: accept_count

      integer :: i, k,  particle_to_move
      real (dp), dimension (d) :: del_dim, del_dim_trial, random_no_arr, trial_positions

      real (dp) :: r_mod_sq, r_cut_sq, sigma_6_div_r_mod_6, L_half, r_mod_sq_trial, lj_pot_trial, lj_pot_old, lj_pot_diff


      call random_number ( random_no_arr(1) )

      particle_to_move = ( int (random_no_arr(1) * N ) ) + 1

      do i = 1, d
          call random_number ( random_no_arr(i) )
      end do

      trial_positions(:) = r(:,particle_to_move) + (  (2.0_dp) * max_move * (random_no_arr(:) - (0.5_dp) )    )

      L_half = (L/2.0_dp)

      lj_pot_old = 0
      lj_pot_trial = 0
      r_cut_sq = r_cut**2


      !check trial_positions - move back inside of box
      do i = 1, d

          trial_positions(i) = trial_positions(i) - L * (nint( trial_positions(i) / L ))

      end do


      do i = 1, N

        if ( i /= particle_to_move ) then


            del_dim (:) = r (:,i) - r (:, particle_to_move)

            del_dim_trial (:) = r (:,i) - trial_positions(:)

            !pbc conditions
            do k = 1, d
                if (del_dim(k) > L_half ) then
                  del_dim(k) = del_dim(k) - L
                end if
                if (del_dim(k) < - L_half ) then
                  del_dim(k) = del_dim(k) + L
                end if

                if (del_dim_trial(k) > L_half ) then
                  del_dim_trial(k) = del_dim_trial(k) - L
                end if
                if (del_dim_trial(k) < - L_half ) then
                  del_dim_trial(k) = del_dim_trial(k) + L
                end if
            end do

            r_mod_sq = sum (del_dim**2)
            r_mod_sq_trial = sum (del_dim_trial**2)


            !cutoff condition
            if  (r_mod_sq_trial  < r_cut_sq)  then

                  sigma_6_div_r_mod_6 = (sigmalj**6)/(r_mod_sq_trial**3)

                  lj_pot_trial = lj_pot_trial + (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - (1.0_dp)) )

            end if

            if  (r_mod_sq  < r_cut_sq)  then

                  sigma_6_div_r_mod_6 = (sigmalj**6)/(r_mod_sq**3)

                  lj_pot_old = lj_pot_old + (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - (1.0_dp)) )

            end if

        end if


      end do


      !difference in energy of old and trial configurations
      lj_pot_diff = lj_pot_trial - lj_pot_old

      !multiply by a factor to get
      lj_pot_diff = lj_pot_diff * (epsilonlj * 4.0_dp )


      !Metropolis test
      call random_number ( random_no_arr(1) )


      if (  random_no_arr(1)  <  exp( (-beta) * lj_pot_diff)      ) then

          r (:, particle_to_move) = trial_positions (:)

          lj_pot = lj_pot + lj_pot_diff

          accept_count = accept_count + 1

      end if

    end subroutine equil_new_position_and_potential


  subroutine new_position_and_potential (r, sigmalj, epsilonlj, d, N, max_move, L, &
    beta, P, accept_count, lj_pot, r_cut, lj_pot_corr, lj_pot_corr_new, pressure_factor, factor )


    real (dp), intent (in) :: sigmalj, epsilonlj, max_move, L, r_cut, beta, lj_pot_corr, lj_pot_corr_new, factor, P
    integer, intent (in) :: d, N
    real (dp), intent (inout) :: r(:,:)
    real (dp), intent (inout) :: lj_pot, pressure_factor
    integer, intent (inout) :: accept_count

    integer :: i, k,  particle_to_move
    real (dp), dimension (d) :: del_dim, del_dim_trial, random_no_arr, trial_positions

    real (dp) :: r_mod_sq, r_cut_sq, sigma_6_div_r_mod_6, L_half, r_mod_sq_trial, lj_pot_trial, lj_pot_old, lj_pot_diff


    call random_number ( random_no_arr(1) )

    particle_to_move = ( int (random_no_arr(1) * N ) ) + 1

    do i = 1, d
        call random_number ( random_no_arr(i) )
    end do

    trial_positions(:) = r(:,particle_to_move) + (  (2.0_dp) * max_move * (random_no_arr(:) - (0.5_dp) )    )

    L_half = (L/2.0_dp)

    lj_pot_old = 0
    lj_pot_trial = 0
    r_cut_sq = r_cut**2


    !check trial_positions - move back inside of box
    do i = 1, d

        trial_positions(i) = trial_positions(i) - L * (nint( trial_positions(i) / L ))

    end do


    do i = 1, N

      if ( i /= particle_to_move ) then


          del_dim (:) = r (:,i) - r (:, particle_to_move)

          del_dim_trial (:) = r (:,i) - trial_positions(:)

          !pbc conditions
          do k = 1, d
              if (del_dim(k) > L_half ) then
                del_dim(k) = del_dim(k) - L
              end if
              if (del_dim(k) < - L_half ) then
                del_dim(k) = del_dim(k) + L
              end if

              if (del_dim_trial(k) > L_half ) then
                del_dim_trial(k) = del_dim_trial(k) - L
              end if
              if (del_dim_trial(k) < - L_half ) then
                del_dim_trial(k) = del_dim_trial(k) + L
              end if
          end do

          r_mod_sq = sum (del_dim**2)
          r_mod_sq_trial = sum (del_dim_trial**2)


          !cutoff condition
          if  (r_mod_sq_trial  < r_cut_sq)  then

                sigma_6_div_r_mod_6 = (sigmalj**6)/(r_mod_sq_trial**3)

                lj_pot_trial = lj_pot_trial + (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - (1.0_dp)) )

          end if

          if  (r_mod_sq  < r_cut_sq)  then

                sigma_6_div_r_mod_6 = (sigmalj**6)/(r_mod_sq**3)

                lj_pot_old = lj_pot_old + (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - (1.0_dp)) )

          end if

      end if


    end do


    !difference in energy of old and trial configurations
    lj_pot_diff = lj_pot_trial - lj_pot_old

    !multiply by a factor to get
    lj_pot_diff = lj_pot_diff * (epsilonlj * 4.0_dp )


    !Metropolis test
    call random_number ( random_no_arr(1) )

    !  - if lt 0 faster? - requires nested if statements

    if (  random_no_arr(1)  <  exp( (-beta) * lj_pot_diff)      ) then

        r (:, particle_to_move) = trial_positions (:)

        lj_pot = lj_pot + lj_pot_diff

        accept_count = accept_count + 1

        call pressure_calc ( beta, r, lj_pot, lj_pot_corr, lj_pot_corr_new, L, N, d, r_cut, sigmalj, epsilonlj, &
        pressure_factor, factor)

    end if

  end subroutine new_position_and_potential



  subroutine write_out_params (path_out, sigmalj, epsilonlj, L, r_cut, max_move, N_moves, Equil_moves, &
    temperature, N, density, r, lj_pot, lj_pot_corr)

    character*(*) , intent (in)  :: path_out
    integer, intent (in) :: N_moves, N, Equil_moves
    real (dp), intent (in) ::  sigmalj, epsilonlj, L, r_cut, max_move, temperature, r(:,:), lj_pot, &
    lj_pot_corr, density

    integer :: i, unit

    unit = 10

    open (unit,file=path_out//'input_coords.xyz')

      write (unit, *) N
      write (unit,*) ''

      do i = 1, N
          write (unit,*)   'Ar',   r (:,i)
      end do

    close (unit)

    unit = 30
    open (unit,file=path_out//'input_params.txt')

      write (unit, 100 ) 'sigma', sigmalj
      write (unit, 100 ) 'epsilon', epsilonlj
      write (unit, *) 'box length', L
      write (unit, *) 'cutoff distance', r_cut
      write (unit, 100) 'maximum move distance', max_move
      write (unit, 200) 'Number of moves', N_moves
      write (unit, 200) 'Number of equilabration moves', Equil_moves
      write (unit, *) 'Temperature', temperature
      write (unit, 200 ) 'Number of atoms', N
      write (unit, *) 'Density', density
      write (unit, *) 'potential energy correction', lj_pot_corr
      write (unit, *) 'initial potential', lj_pot+lj_pot_corr


      100 format (1x,a,5x,F7.5)
      200 format (1x,a,5x,I0)

    close (unit)

  end subroutine write_out_params



  subroutine pressure_and_lj_pot_correction (pi, sigmalj, epsilonlj, r_cut, density, P_corr, lj_pot_corr, factor, lj_pot_corr_new )

    implicit none

    real(dp), intent (in) :: pi, sigmalj, epsilonlj, r_cut, factor, density
    ! integer, intent (in) :: N

    real(dp), intent (out) ::P_corr, lj_pot_corr, lj_pot_corr_new

    real (dp) :: r_cut_new, sigma_div_rcut_3

    sigma_div_rcut_3 = (sigmalj/r_cut)**3


    lj_pot_corr = (sigmalj**3)*pi*epsilonlj*density*(8.0_dp/3.0_dp) * ( ( (sigma_div_rcut_3**3) / (3.0_dp) ) - (sigma_div_rcut_3) )


    ! P corr for truncated energy systems
     P_corr = ( (sigma_div_rcut_3**3) - sigma_div_rcut_3 ) * (pi*epsilonlj*(sigmalj**3))  * (density**2) * (8.0_dp/3.0_dp)


    r_cut_new = r_cut * factor
    sigma_div_rcut_3 = (sigmalj/r_cut_new)**3

    ! lj_pot_corr_new = ( ( sigma_div_rcut_3*(sigmalj**3)*pi*epsilonlj* N_factor ) / (L**3) ) &
    !  * ( ( (sigma_div_rcut_3**2) / (3.0_dp) ) - (1.0_dp) )
    lj_pot_corr_new = (sigmalj**3)*pi*epsilonlj*density*(8.0_dp/3.0_dp) &
    * ( ( (sigma_div_rcut_3**3) / (3.0_dp) ) - (sigma_div_rcut_3) )


  end subroutine pressure_and_lj_pot_correction



  subroutine write_out_s_to_file (i, r, N, unit1 )

    integer, intent (in) :: i, N, unit1
    real (dp), intent (in) :: r(:,:)

    integer :: j

    write (unit1,*) N
    write (unit1,*) i
    do j = 1, N
        write (unit1,*) 'Ar', r (:,j)
    end do

  end subroutine write_out_s_to_file



  subroutine equil_lj_pot_stats (i, lj_pot, lj_pot_mean, lj_pot_mean_prev, lj_pot_pre_stdev, lj_pot_stdev)

    implicit none

    integer, intent (in) :: i
    real (dp), intent (in) :: lj_pot
    real (dp), intent (inout) :: lj_pot_mean, lj_pot_mean_prev, lj_pot_pre_stdev, lj_pot_stdev

    !Welford algorithm
    lj_pot_mean = lj_pot_mean_prev + (  (lj_pot - lj_pot_mean_prev) / ( real (i, dp) )   )
    lj_pot_pre_stdev = lj_pot_pre_stdev + (  (lj_pot - lj_pot_mean_prev)*(lj_pot - lj_pot_mean)  )
    lj_pot_mean_prev = lj_pot_mean

    if  (   i > 1   ) then
        lj_pot_stdev = sqrt( lj_pot_pre_stdev / (real((i-1), dp) ) )
    end if

  end subroutine equil_lj_pot_stats



  subroutine pressure_calc ( beta, r, lj_pot, lj_pot_corr, lj_pot_corr_new, L, N, d, r_cut, sigmalj, epsilonlj, &
    pressure_factor, factor)

    integer, intent (in) :: N, d
    real (dp), intent (in) :: beta, lj_pot, lj_pot_corr, lj_pot_corr_new, L, r_cut, sigmalj, epsilonlj, factor
    real (dp), intent (in) :: r (:,:)
    real (dp), intent (out) :: pressure_factor

    integer :: i, j
    real (dp) :: lj_pot_new, L_new, r_cut_new
    real (dp) :: r_new (d,N)


    L_new = L * factor
    r_cut_new = r_cut * factor
    r_new =  r * factor

    write (*,*) 'new', r_new(1,1), lj_pot_new , sigmalj, epsilonlj, d, N, L_new, r_cut_new
    write (*,*) 'old', r(1,1), lj_pot , sigmalj, epsilonlj, d, N, L, r_cut

    call init_pot ( r_new, lj_pot_new , sigmalj, epsilonlj, d, N, L_new, r_cut_new )

    ! write (*,*) ((L_new/L)**(3*N))
    write (*,*) (exp (beta * ((lj_pot)-(lj_pot_new)) ))
    write (*,*) ((lj_pot)-(lj_pot_new))
    write (*,*) lj_pot, lj_pot_new



    pressure_factor = ((L_new/L)**(3*N))   * (exp (beta * ((lj_pot)-(lj_pot_new)) ))


    ! pressure_factor = ((L_new/L)**(3*N))   * (exp (beta * ((lj_pot+lj_pot_corr)-(lj_pot_new+lj_pot_corr_new)) ))

  end subroutine pressure_calc


end program monte
