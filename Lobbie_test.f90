!  Lobbie_test.f90 
!
!  FUNCTIONS:
!  Lobbie_test - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Lobbie_test
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Lobbie_test
    use lobbie2_mod
    implicit none

    ! Variables

    real(8), parameter :: pi = 4.d0*atan(1.d0)
    real(8), parameter :: k_gauss = 0.017202098949957226d0
    real(8), parameter :: days_in_year = 365.25636556 !  365.2421904020d0
    real(8), parameter :: au_km = 149597870.691d0
    real(8), parameter :: rad = 180.0_dp/pi
    
	! For Lobbie
	real(8), dimension(:) :: posvel(6), elements(6), el(11)
    real(8), dimension(:), allocatable :: pos, vel, etol, pos_ini, vel_ini, t_collision
    real(qp),dimension(:), allocatable :: pos32, vel32
    integer(4), dimension(:), allocatable :: indices, err_ar
	real(8) :: ts !, tf	! time start and finish
	real(8) :: step
	integer(4) :: nxy,nst,ncf !,ns,ni
    integer(4) :: num_particles
	
	integer(4) :: i, j, k, num_steps, index, ios, stat
	
    real(8) :: r_vector, t_curr, step_tt, period,period_ini, major_axis, velocity, r_perihelion, alpha
    real(dp), dimension(:), allocatable :: semi_major_axes, posvel_pl, radius
    real(8) :: t_prog_start, t_prog_final, elapsed_time
    integer(8) :: start, finish, rate, count_max
	
    ! constants that we read from a file
    real(8) :: step_t, tf, accuracy_step, max_time
    integer(4) :: ns, ni    ! number of nodes of integrator; number of itirations in the integrator
    integer(4) :: current_file, checkpoint_output
    character(20) :: elements_char, line, elem_output, short_output
    character(80) :: file_check_name
    logical :: checkpoint_exists
    integer(4) :: n_lines, start_k
    
    
    
    ! read the parameters of the integrator from a file
    open(14, file = 'parameters.dat')
    read(14,*) ns               ! number of nodes of integrator
    read(14,*) ni               ! number of itirations in the integrator
    read(14,*) tf               ! final time
    read(14,*) step_t           ! step in time for output
    read(14,*) accuracy_step    ! accuracy on each step (1e-8, 1e-10...)
    read(14,*) elements_char    ! either cartesian state vectors or orbital elements (a,e,i,Om,g,M). Default - cartesians
    read(14,*) elem_output      ! if we need to compute orbital elements and print them in a file
    read(14,*) short_output     ! if we need write only the particles when they are NEO (significanlty save space on a drive)
    read(14,*) max_time         ! maximal time for the run of the program in minutes
    read(14,*) checkpoint_output ! once in how many times should we print checkpoints
    close(14)
    ! multiply by 2pi because GM = 1
    tf = tf * 2.d0*pi
    step_t = step_t * 2.d0*pi
    max_time = max_time * 60.0d0    ! transform to seconds
    if(checkpoint_output <=0) checkpoint_output = 1
    
    if_relativity_include = .false.; if_sun_radiation_press_include=.false.; if_yarkovskii_include = .false.
    !if_relativity_include = .true.; if_sun_radiation_press_include=.true.; if_yarkovskii_include = .true.
    
    !--------------------------------------------------------------------------------------------------------------------
    ! read data associated with planets
    !--------------------------------------------------------------------------------------------------------------------
    
    ! file with coordinates and masses of planets (ecliptical in au, au/days, mass in mass of the Sun)
    open(178, file = 'planets.in')
    read(178,*) N_planets
    allocate(mass_pl(N_planets), semi_major_axes(N_planets), posvel_pl(6*N_planets), radius_pl(0:N_planets), sphere_of_influence(N_planets))

    ! read the planets masses and initial coordinates and velocities
    do i = 1, N_planets
        read(178,*) mass_pl(i)
        read(178,*) posvel_pl( (i-1)*6+1:i*6 )
    enddo
    close(178)
    
    ! file with radiuses of planets and the Sun (Sun is 0-th planet)
    open(178, file = "planets_radiuses.in")
    do i = 0, N_planets
        read(178,*) radius_pl(i)
    enddo
    radius_pl = radius_pl / au_km
    close(178)
    
    ! file with schere of influences of planets
    open(178, file = "sphere_of_influences.in")
    do i = 1, N_planets
        read(178,*) sphere_of_influence(i)
    enddo
    sphere_of_influence = sphere_of_influence / au_km
    close(178)
    
    !--------------------------------------------------------------------------------------------------------------------
    ! read coordinates and velocities of test particles
    !--------------------------------------------------------------------------------------------------------------------
    
    write(*,*) 'Input the filename of asteroid orbits in [filename].txt (without .txt):'
    read(*,*) file_char !current_file
    
    !call get_command_argument(1, file_char)
    
    open(14,file = 'asteroids/' // trim(file_char) // '.txt', status='old', action='read')
    
    ! compute the number of particles
    N_particles = 0
    do
        read(14, '(A)', iostat=ios) line
        if (ios /= 0) exit  ! Exit loop on end-of-file or error
        N_particles = N_particles + 1
    end do
    rewind(14)
    
    Num_particles = N_particles + N_planets
    
    step = 5.d-2		! initial step
	nxy = 3 * Num_particles
	ts = 0.d0
    
    ! allocate the data
    allocate(pos(nxy), vel(nxy), etol(2*nxy)) !, posvel(2*nxy))
    allocate(pos32(nxy), vel32(nxy),indices(nxy/3),err_ar(nxy/3),t_collision(nxy/3))
    
    ! put as accuracy on each step the same for each of the particle and planet
    etol = accuracy_step
    
    
    !-------------------------------------------------------------------------------------------------------------------------
    ! CHECKPOINTS
    !-------------------------------------------------------------------------------------------------------------------------
    
    file_check_name = 'results/' // trim(file_char) // '_checkpoints.txt'
    ! check if it already exists
    inquire(file=file_check_name, exist=checkpoint_exists)
    
    ! if _checkpoint already exists then we read the data from it
    if (checkpoint_exists) then

        ! First count how many lines are in the file
        open(unit=167, file=file_check_name, status='old', iostat=ios)
        ! count number of lines in the file
        n_lines = 0
        do
            read(167, '(A)', iostat=stat) line
            if (stat /= 0) exit
            n_lines = n_lines + 1
        end do
        rewind(167)
    
        ! if the file is empty or has less lines that Num_particles we close it
        ! and read the data from original files
        if(n_lines < Num_particles) then
            checkpoint_exists = .false.
            close(167)
        endif
    endif
    
    
    if(checkpoint_exists) then

        close(14)
        
        ! Read again to store all lines
        do i = 1, n_lines - Num_particles
            read(167, *) 
        end do
        
        ! read the data from _checkpoint file
        do i = 1, Num_particles
            read(167,*) indices(i), t_collision(i), pos32(i*3-2:i*3), vel32(i*3-2:i*3), err_ar(i), step
        enddo
    
        pos = pos32; vel = vel32
        
        ! count nxy
        nxy = count(err_ar == 0) * 3
        
    ! _checkpoint DOESN'T exist
    else
    
    
        ! coordinates of planets in system where GM = 1
        do i = 1, N_planets        
            pos(i*3-2:i*3) = posvel_pl( (i-1)*6+1:(i-1)*6+3 )
            vel(i*3-2:i*3) = posvel_pl( (i-1)*6+4:i*6 ) / k_gauss
        enddo 
        
        ! if we read orbital elements for particles we transform them into ecliptical state vectors 
        if(elements_char(1:4) == 'elem' ) then
            do i = 1, N_particles
                ! read elements: a(au), e, i(deg), Om(deg), g(deg), M(deg)
                read(14,*) elements(:)
                ! transpose to radians
                elements(3:6) = elements(3:6) / rad
                ! compute state vectors
                call ElementsPosVel_new(elements(1)*(1._dp - elements(2)), elements(2), elements(5), elements(4), elements(3), elements(6), elements(1), 0.0_dp, 0.0_dp, posvel)
                j = i + N_planets
                pos(j*3-2:j*3) = posvel(1:3); vel(j*3-2:j*3) = posvel(4:6) / k_gauss
            enddo
        else
            ! read coordinates and velocities of particles
            do i = N_planets+1, Num_particles
                read(14,*) pos(i*3-2:i*3), vel(i*3-2:i*3)
                vel(i*3-2:i*3) = vel(i*3-2:i*3) / k_gauss
            enddo
        endif
        close(14)

        err_ar(:) = 0
        t_collision(:) = 0.d0
        indices(:) = [(i,i=1,nxy/3)]
    
        ! create file for checkpoints
        open(unit=167, file=file_check_name, iostat=ios)
        
        ! copy coordinates and velocities to real(16) data
        pos32 = pos; vel32 = vel
        
    ! endif _checkpoint exists
    endif
    
    !--------------------------------------------------------------------------------------------------------------------
    ! open file for output the results
    !--------------------------------------------------------------------------------------------------------------------
    
    ! files for writing the results
    open(17, file = 'results/' // trim(file_char) // '.txt', position='append', action='write')


    CALL SYSTEM_CLOCK(start, rate, count_max)


    ! compute the starting point if _checkpoint exists
    if(.not. checkpoint_exists) then
        start_k = 0
    else
        start_k = nint(maxval(t_collision(:), mask = err_ar(:) == 0) / step_t)
    endif
    
 
    do k = start_k, nint(tf/step_t)-1
        t_curr = step_t*k
        
        call lobbie2_2(pos32(:nxy),vel32(:nxy),0.0_dp,step_t,step,etol,nxy,ns,ni,nst,ncf,fun32,indices(:nxy),err_ar(:nxy),t_collision(:nxy),t_curr)
        pos(:nxy) = pos32(:nxy); vel(:nxy) = vel32(:nxy)
        
        ! change the number of particles that are integrated
        nxy = count(err_ar == 0) * 3
        
        ! writing the results in a file
        do j = 1, Num_particles
            
            ! find where j-th particle is now (because we rearrange the arrays and put collisional particles at the end
            i = findloc(indices, j, dim = 1)
            
            ! if the particle didn't collide we plot the current time, otherwise the collisional time
            if(err_ar(i) == 0) then
                t_collision(i) = (t_curr+step_t)
            else
                if(abs(t_collision(i)) <= abs(step_t)) t_collision(i) = t_collision(i) + t_curr
            endif
            
            ! if we only plot NEOs then we first compute the orbital elements
            if(short_output(1:1) == 'y') then
                
                ! if it is a planet we cycle
                if(i <= N_planets) cycle
                
                ! if the particle didn't collide or escape
                if(err_ar(i) == 0) then
                    
                    ! get orbital elements
                    call PosVelElements_new(0.d0, [pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss], 0.d0, 0.d0, el(1), el(2), el(3), el(4), el(5), el(6), el(7), el(8), el(9), el(10), el(11)) 
                    
                    ! if q < 1.3
                    if( el(10)*(1.d0 - el(3)) <= 1.3d0) then
                    ! write the results                             N of particle  time           coordinates     velocities              err         a       e       i       Om      g   M 
                    write(17,'(i0,",",7(g0.16,","), *(g0.8,:,","))') j,t_collision(i)/(2.d0*pi), pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss, err_ar(i), el(10), el(3), [el(7), el(6), el(5), el(8)]*rad
                    endif
                    
                else
                    ! if the particle collided in this run
                    if( (t_collision(i) - t_curr)/step_t <= 1.d0 .and. (t_collision(i) - t_curr)/step_t >= 0.d0 ) then
                        ! get orbital elements
                        call PosVelElements_new(0.d0, [pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss], 0.d0, 0.d0, el(1), el(2), el(3), el(4), el(5), el(6), el(7), el(8), el(9), el(10), el(11)) 
                        
                        ! write the results                             N of particle  time           coordinates     velocities              err         a       e       i       Om      g   M 
                        write(17,'(i0,",",7(g0.16,","), *(g0.8,:,","))') j,t_collision(i)/(2.d0*pi), pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss, err_ar(i), el(10), el(3), [el(7), el(6), el(5), el(8)]*rad
                    endif                    
                endif                
                
            ! if we print everything
            else
            
                if( elem_output(:2) == 'no') then
                    write(17,'(i0,",",*(g0.16,:,","))') j,t_collision(i)/(2.d0*pi), pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss, err_ar(i)
                else
                
                    ! get orbital elements
                    !call PosVelElements_new(ETIME, [pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss], 0.d0, 0.d0, PERIT, Q,     E,     P,     W,    OM,     I,     M,     N,     A,      Q2)
                    call PosVelElements_new(0.d0, [pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss], 0.d0, 0.d0, el(1), el(2), el(3), el(4), el(5), el(6), el(7), el(8), el(9), el(10), el(11)) 
                    
                    ! write the results                  N of particle  time        coordinates     velocities              err         a       e       i       Om      g   M 
                    write(17,'(i0,","g0.16,",",6(g0.5,","),*(g0.16,:,","))') j,t_collision(i)/(2.d0*pi), pos(i*3-2:i*3), vel(i*3-2:i*3)*k_gauss, err_ar(i), el(10), el(3), [el(7), el(6), el(5), el(8)]*rad
                endif
            endif

        enddo
        
        ! if all the particles collided or left we stop integration
        if(nxy == 3*N_planets) then
            ! write the checkpoint of data to a file
            do i = 1, Num_particles
                write(167,'(i0,",", g0.16,",",*(g0.8,:,","))') indices(i),t_collision(i), pos32(i*3-2:i*3), vel32(i*3-2:i*3), err_ar(i), step
            enddo
            write(*,*) 'current time/final time', real(step_t*(k+1)/(2.d0*pi)), real(tf/(2.d0*pi)), nxy/3
            
            print *, 'integrations finished since all particles left'
            
            close(167)
            close(17)
            
            exit
        endif
        
        ! write the checkpoint of data to a file
        if(mod(k,checkpoint_output) == 0) then
            do i = 1, Num_particles
                write(167,'(i0,",", g0.16,",",*(g0.38,:,","))') indices(i),t_collision(i), pos32(i*3-2:i*3), vel32(i*3-2:i*3), err_ar(i), step
            enddo
        endif
        
        
        write(*,*) 'current time/final time', real(step_t*(k+1)/(2.d0*pi)), real(tf/(2.d0*pi)), nxy/3
        
        ! compute the time we are here
        CALL SYSTEM_CLOCK(COUNT=finish)
        elapsed_time = real(modulo(finish - start, count_max),16) / real(rate,16)
        
        write(*,*) 'computation time', elapsed_time
        
        if(elapsed_time > max_time .and. max_time > 0.d0) then
            ! write the checkpoint of data to a file
            do i = 1, Num_particles
                write(167,'(i0,",", g0.16,",",*(g0.38,:,","))') indices(i),t_collision(i), pos32(i*3-2:i*3), vel32(i*3-2:i*3), err_ar(i), step
            enddo
            ! stop the program
            stop
        endif
        
        call flush(6)
        call flush(167)
        call flush(17)

    enddo
    

    ! if we are here then we finished the integrations!
    ! create a file that indicates that the program finished
    open(756, file = 'done/' // trim(file_char) // '.flag')
    write(756,*) 'done'
    close(756)

    
    
    !CALL CPU_TIME(t_prog_final)
    CALL SYSTEM_CLOCK(COUNT=finish)
    elapsed_time = real(modulo(finish - start, count_max),16) / real(rate,16)
    
    write(17,*) 'computation time', elapsed_time
    print *, 'computation time', elapsed_time
    
    ! Body of Lobbie_test
    print *, 'well done'
    
    end program Lobbie_test

