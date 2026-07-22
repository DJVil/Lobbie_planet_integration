module lobbie2_mod
integer(4) :: N_planets, N_particles
real(8), dimension(:), allocatable :: mass_pl, radius_pl, sphere_of_influence
integer(4), parameter :: qp = 16, dp = 8
!real(dp), dimension(8) :: sphere_of_influence = [0.117_dp,0.616_dp,0.929_dp,0.578_dp,48.2_dp,54.5_dp,51.9_dp,86.2_dp]/149.597870691_dp
character(20) :: file_char
logical :: if_relativity_include, if_sun_radiation_press_include, if_yarkovskii_include
    
real(dp), dimension(56), private :: spacing
data spacing/                        &
 0.27639320225002103035908263312687, &
 0.17267316464601142810085377187657, &
 0.11747233803526765357449851302033, &
 0.35738424175967745184292450297956, &
 0.08488805186071653506398389301626, &
 0.26557560326464289309811405904562, &
 0.06412992574519669233127711938966, &
 0.20414990928342884892774463430102, &
 0.39535039104876056561567136982732, &
 0.05012100229426992134382737779083, &
 0.16140686024463112327705728645432, &
 0.31844126808691092064462396564567, &
 0.04023304591677059308553366958883, &
 0.13061306744724746249844691257008, &
 0.26103752509477775216941245363437, &
 0.41736052116680648768689011702091, &
 0.03299928479597043283386293195030, &
 0.10775826316842779068879109194577, &
 0.21738233650189749676451801526112, &
 0.35212093220653030428404424222047, &
 0.02755036388855888829620993084839, &
 0.09036033917799666082567920914154, &
 0.18356192348406966116879757277817, &
 0.30023452951732553386782510421652, &
 0.43172353357253622256796907213015, &
 0.02334507667891804405154726762227, &
 0.07682621767406384156703719645062, &
 0.15690576545912128696362048021682, &
 0.25854508945433189912653138318153, &
 0.37535653494688000371566314981289, &
 0.02003247736636954932244991899228, &
 0.06609947308482637449988989854586, &
 0.13556570045433692970766379973955, &
 0.22468029853567647234168864707046, &
 0.32863799332864357747804829817916, &
 0.44183406555814806617061164513192, &
 0.01737703674808071360207430396519, &
 0.05745897788851185058729918425888, &
 0.11824015502409239964794076201185, &
 0.19687339726507714443823503068163, &
 0.28968097264316375953905153063071, &
 0.39232302231810288088716027686354, &
 0.01521597686489103352387863081627, &
 0.05039973345326395350268586924007, &
 0.10399585406909246803445586451842, &
 0.17380564855875345526605839017970, &
 0.25697028905643119410905460707656, &
 0.35008476554961839595082327263885, &
 0.44933686323902527607848349747704, &
 0.01343391168429084292151024906313, &
 0.04456000204221320218809874680113, &
 0.09215187438911484644662472338123, &
 0.15448550968615764730254032131377, &
 0.22930730033494923043813329624797, &
 0.31391278321726147904638265963237, &
 0.40524401324084130584786849262344/
    
    contains

subroutine lobbie2_2(x32,y32,ts,tf,step,etol_ini,nxy,ns,ni,nst,ncf,fun,indices,err_ar,t_collision,t_global)

!     The integrator Lobbie II numerically solves ordinary              *
!     differential equations of second order x''=f(t,x,x').             *
!     The integrator is a software implimentation of a collocation      *
!     Lobatto method based on Newton's polynomial interpolants          *
!     for the function f. The theoretical background of the method is   *
!     presented in Avdyushev V.A. New Collocation Integrator for        *
!     Solving Dynamic Problems // Russ. Phys. J. Vol. 63.               *
!     P. 1977–1988 (2021).                                              *
!                                                                       *
!     Input parameters :                                                *
!                                                                       *
!     x,y    initial values for x,x' at ts                              *
!     ts     initial value for t                                        *
!     tf     final value for t                                          *
!     step   start stepsize                                             *
!     etol   local tolerances                                           *
!     nxy    dimension of x & y                                         *
!     ns     number of nodes (2 <= ns <= 17)                            *
!     ni     number of iterations                                       *
!     fun    name of subroutine for f                                   *
!                                                                       *
!     Notes :                                                           *
!                                                                       *
!     1. The order of the integration scheme is p = 2 ns - 2.           *
!     2. If ni < 1, iterative process runs until convergence            *
!        (up to 25 iterations).                                         *
!                                                                       *
!     Output parameters :                                               *
!                                                                       *
!     x,y    final values for x,x' at tf                                *
!     nst    number of steps                                            *
!     ncf    number of calls of fun                                     *
!                                                                       *
!     Subroutine for f :                                                *
!                                                                       *
!     fun(t,x,y,f)                                                      *
!                                                                       *
!     Input parameters :                                                *
!                                                                       *
!     x,y    values for x,x' at t                                       *
!                                                                       *
!     Output parameters :                                               *
!                                                                       *
!     f      values for f (dim f = nxy)                                 *
!                                                                       *
!************************************************************************
      
      implicit real(8)(a-h,o-z)
      logical ls,vs, r_l_rmn, r_g_rmx
      parameter (ni0=25,sgm=3.16_dp,eps=1.0E-14_dp)
      real(dp), parameter :: pi = 4._dp * atan(1._dp)
      integer*4, dimension(nxy) :: ietol
      real(dp), dimension(nxy)    :: x,y,f,etol,u,v,us, f_ns, y_ns
      real(qp), dimension(nxy)    :: x32, y32
      real(qp) :: t32
      real(dp), dimension(nxy/3,3) :: xc
      integer(4), dimension(nxy/3) :: indices, err_ar
      real(dp),dimension(nxy/3),intent(out) :: t_collision
      real(dp), intent(in) :: t_global
      real(dp), dimension(ns)     :: c
      real(dp), dimension(ns+1)   :: g
      real(dp), dimension(ns,nxy) :: a,da,ax, tmp_a(ns,3), a_previous, ax_previous
      real(dp), dimension(nxy,ns) :: contributions
      real(dp), dimension(ns,ns)  :: b,d,e
      real(dp), dimension(nxy)    :: tmp_array
      integer(4) :: npc             ! number of variables to be computed 3*(N_planets + particles that didn't collide yet)
      real(dp), dimension(nxy/3) :: dist_sun
      real(dp), dimension(nxy/3,N_planets) :: dist_pl
      real(dp) :: ecal, r
      real(dp), dimension(nxy) :: etol_ini

      real(qp), dimension(nxy)    :: x32_previous, y32_previous
      real(dp) :: t_previous, t_inside
      logical, dimension(nxy/3,N_planets) :: if_close, if_close_saved
      real(dp), dimension(nxy/3,14) :: close_enc_inf        ! close encounter information (time, distance, state vect of particle, state vect of planet)
      logical :: encounter_mode, was_inside, inside
      !real(dp), dimension(56)     :: spacing
      !common/spacing/spacing

      external fun
      
      npc = nxy
      
      etol = etol_ini
      !etol_ini = etol
      if_close = .false.; if_close_saved = .false.
      encounter_mode = .false.; was_inside = .false.; inside = .false.
      close_enc_inf(:,2) = 100000_dp
      x32_previous = x32; y32_previous = y32; t_previous = ts

      x = x32; y = y32
      
      if(ts==tf) return

      if(ns<2.or.ns>17) stop 'Inadmissible Order!'
      
      ! Spacing

      nr=ns/2-1; l=((ns-1)/2-1)*(ns/2-1)
      c(1)=0.0_dp; c(ns)=1.0_dp; if(ns>2) c(nr+2)=0.5_dp
      c(2:nr+1)=spacing(l+1:l+nr)
      c(ns-1:ns-nr:-1)=1.0_dp-spacing(l+1:l+nr)
      r_l_rmn=.false.; r_g_rmx=.false.

      ! Node Differences
      
      d=0.0_dp
      do i=1,ns
      do j=1,ns
        d(i,j)=c(i)-c(j)
      end do
      end do

      ! Constants

      do i=1,ns

      prd=1.0_dp
      do k=1,ns+1
        prd=prd*c(i)/dfloat(k) 
        g(k)=prd
      end do
      b(i,1)=g(1)
      e(i,1)=g(2)

      do j=2,ns
        do k=1,ns-j+2
            g(k)=d(i,j-1)*g(k)-dfloat(k)*g(k+1)
        end do
        b(i,j)=g(1)
        e(i,j)=g(2)
      end do

      end do

      ! Initialization
      
      t=ts; t32=ts; h=abs(step)*(tf-ts)/abs(tf-ts)
      ncf=1; ncf0=ni0*(ns-1); nst=0
      a=0.0_dp; ax=0.0_dp; da=0.0_dp; us=0.0_dp

      ! Variable Step-Size Mode

      netol=0
      do i=1,nxy
      if(etol(i)>0.0_dp) then
        netol=netol+1; ietol(netol)=i
      end if
      end do

      vs=netol/=0; deg=1.0_dp/dfloat(ns)
      r=1.0_dp; rmx=sgm**deg; rmn=1.0_dp/rmx

      call fun(npc,real(t,dp),real(x,dp),real(y,dp),a(1,:))

      ! Start Step-Size

      if(h==0.0_dp) then
        if(.not.vs) stop 'Zero Step!'
        eta=sqrt(eps); ecal=0.0_dp
        do while(ecal==0.0_dp)
            eta=eta*10.0_dp
            u=x+eta*y+eta**2*a(1,:)/2.0_dp
            v=y+eta*a(1,:)
            call fun(npc,real(t+eta,dp),real(u,dp),real(v,dp),f); ncf=ncf+1
            do i=1,netol
                j=ietol(i)
                ecal=ecal+((f(j)-a(1,j))/etol(j))**2
            end do
        end do
        h=sqrt(2.0_dp*eta/sqrt(ecal))
        h=h*(tf-ts)/abs(tf-ts)
      end if

    1 h=r*h; h2=h**2; ls=.false.

      ! First & Last Step

      if((tf-t)/h<1.0_dp) then
        h=tf-t; h2=h**2; ls=.true.
      end if

      ! Integration
      
      do

      ! Iterations

    2 do it=1,ni

        do i=2,ns
            
            ! Collocation Solution
            
            f=0.0_dp
            !do j=1,ns
            do j=ns,1,-1 
                f(:npc)=f(:npc)+e(i,j)*a(j,:npc)
            end do
            u(:npc)=x(:npc)+h*c(i)*y(:npc)+h2*f(:npc)       ! important that the order of summation is this
            !u = u + x
            
            f=0.0_dp
            !do j=1,ns
            do j=ns,1,-1 
                f(:npc)=f(:npc)+b(i,j)*a(j,:npc)
            end do
            v(:npc)=y(:npc)+h*f(:npc)
            
            ! save y on i = ns
            if(i .eq. ns) y_ns(:npc) = v(:npc)

            ! Function

            call fun(npc,real(t+h*c(i),dp),real(u(:npc),dp),real(v(:npc),dp),f(:npc)); ncf=ncf+1

            ! save f on i = ns
            if(i .eq. ns) f_ns(:npc) = f(:npc)
            
            ! Divided Difference

            do j=1,i-1
                f(:npc)=(f(:npc)-a(j,:npc))/d(i,j)
            end do
            a(i,:npc)=f(:npc)

        end do

        ! Convergence Check

        if(all(abs(u(:npc)-us(:npc))<=abs(eps*u(:npc)))) exit   
        us(:npc)=u(:npc)

      end do

      ! Forced Convergence
      
      if(nst==0.and.ncf<ncf0.and.it>ni) goto 2
      
      ! compute y_ns as the 3d vector velocity lenght
      do i = 1, npc/3
          y_ns((i-1)*3+1:i*3) = sqrt(sum( (y_ns((i-1)*3+1:i*3))**2 ))
      enddo
      
      ! Scale Factor

      if(vs) then

        ecal=0.0_dp
        do i=1,netol
            j=ietol(i)
            !ecal=ecal+(a(ns,j)/etol(j))**2
            ! this is the dimensionless step for velocity
            ecal = max(ecal, abs(h) * abs( &
                a(ns,j) / merge(y_ns(j),1.d0,y_ns(j)/=0.d0) * &
                merge(1.d0,0.d0,y_ns(j)/=0.d0) ) * deg / etol(j))
        end do
        if(ecal/=0.0_dp) then
            !ecal=abs(h)*sqrt(ecal)*deg
            r=(1.0_dp/ecal)**deg
        else
            r=rmx
        end if

        ! First Step-Size Adjustment

        if(nst==0) then
            if(rmx<r.and.ls  ) then
                goto 3
            endif
            
            ! if r was > rmx before and was < rmn, then we force to proceed
            if(rmx<r .and. r_g_rmx .and. r_l_rmn) goto 3
            
            if(rmx<r) then
                r_g_rmx = .true.
                goto 1
            endif
            if(r<rmn) then
                r_l_rmn = .true.
                goto 1
            endif
            !if(rmx<r.or.r<rmn) goto 1
        end if

        if(r>rmx) r=rmx ! Damping

      end if

      ! Solution

3     f(:npc)=0.0_dp
      
      x32_previous = x32; y32_previous = y32; t_previous = t; 
      
      !do j=1,ns
      do j=ns,1,-1 
        f(:npc)=f(:npc)+e(ns,j)*a(j,:npc)
      end do
      x32(:npc) = x32(:npc) + h*y(:npc)+h2*f(:npc)
      x(:npc)=x32(:npc)
      
      f(:npc)=0.0_dp
      !do j=1,ns
      do j=ns,1,-1  
        f(:npc)=f(:npc)+b(ns,j)*a(j,:npc)
      end do
      y32(:npc)=y32(:npc)+h*f(:npc)
      y(:npc) = y32(:npc)
      
      !t32 = t32 + h; t = t32
      t=t+h; 
      nst=nst+1

	  !========================================================================================================
	  ! Here we check that the particle did not collide with a planet or the Sun
	  ! nor it left the Solar system
	  ! If a particle did so we place it in the end of the array and keep the track of particles indices
	  
      do i = 1, 3
          xc(:,i) = x(i::3)
      enddo
      
      ! Compute the distances to the Sun
      dist_sun(N_planets+1:npc/3) = sqrt( xc(N_planets+1:npc/3,1)**2 + xc(N_planets+1:npc/3,2)**2 + xc(N_planets+1:npc/3,3)**2  )
      
      ! compute the distances to planets    
      do j = 1, N_planets
          dist_pl(N_planets+1:npc/3,j) = sqrt((xc(j,1) - xc(N_planets+1:npc/3,1))**2 + (xc(j,2) - xc(N_planets+1:npc/3,2))**2 + (xc(j,3) - xc(N_planets+1:npc/3,3))**2)
      end do
      
      
      ! Check if we close to a planet we need to make the integration step smaller
      
      ! if something is inside 0.1au sphere to any planet?
      inside = any(dist_pl(N_planets+1:npc/3,:) <= 0.03_dp)
      dmin_dist_to_planets = minval(dist_pl(N_planets+1:npc/3,:))
      
      if (.not. encounter_mode .and. inside) then
              etol = 1.d-14
              !r = min(r, 2._dp*pi/365.25_dp / 100._dp / h)
              r = r / 3.0_dp
              t_inside = t
              x32 = x32_previous; y32 = y32_previous; t = t_previous; 
              x = x32; y = y32
              a = a_previous; ax = ax_previous
              encounter_mode = .true.
              nst = nst - 1
              ls = .false. ! we need to insure this is not the last step (so r = r / 3 is important)
              go to 5
      endif
      
      
      if (encounter_mode) then
          if (inside) was_inside = .true.
          etol = 1.d-14
      
          ! if we were inside but not now we restore the step
          ! or the inside time was 50 days ago and we are not inside
          if (was_inside .and. .not. inside .or. (abs(t_inside-t) > 50._dp/365.25*pi*2 .and. .not. inside) ) then
              ! We have exited the sphere after truly being inside
              encounter_mode = .false.
              was_inside = .false.
              etol = etol_ini
          endif
      endif
      
      
      ! Check for close encounters (distance < 50 * radius_planet)
      do i = N_planets+1,npc/3
      do j = 1, N_planets
          if( dist_pl(i,j) <= sphere_of_influence(j) ) then !50.d0 * radius_pl(j) ) then
              !r = min(r,  2._dp*pi/365.25_dp / 100._dp / h)
              ! if the new distance is smaller we save it
              if( dist_pl(i,j) < close_enc_inf(i,2) .or. (.not. if_close(i,j)) ) then
                  close_enc_inf(i,:) =[ t_global+t,dist_pl(i,j), x(3*i-2:3*i), y(3*i-2:3*i),x(3*j-2:3*j), y(3*j-2:3*j) ]
              else
                  ! if we didn't save the close encounter info before
                  if( .not. if_close_saved(i,j) ) then
                      open(143,file='results/' // trim(file_char) // '_encounter.txt', position='append', action='write')
                      write(143,'(2(i0,","),*(g0.16,:,","))') indices(i),j,close_enc_inf(i,1)/(2._dp*pi),close_enc_inf(i,2), close_enc_inf(i,3:)
                      close(143)
                      if_close_saved(i,j) = .true.
                  endif
              endif
              if_close(i,j) = .true.
          else
              if_close(i,j) = .false.
              if_close_saved(i,j) = .false.
          endif
      enddo
      enddo

      
      ! check that particles didn't collide with planets of left the Solar system
      i = N_planets+1
      do
          if(i > npc/3) exit
          ! if any collision or leaving the solar system
          if( any(dist_pl(i,:) <= radius_pl(1:N_planets)) .or. dist_sun(i) <= radius_pl(0) .or. dist_sun(i) >= 100.d0 ) then
   
              if(dist_sun(i) <= radius_pl(0)) err_ar(i) = -100
              if(dist_sun(i) >= 100.d0) then
                  ! if the energy says it is unbounding orbit
                  if( sum(y(3*i-2:3*i)**2)/2._dp - 1._dp /sqrt(sum(x(3*i-2:3*i)**2)) >= 0.d0 .or. dist_sun(i) >= 200.d0) then
                    err_ar(i) = -99
                  else
                      i = i + 1
                      cycle
                  endif
              endif
              
              do j = 1, N_planets
                  if(dist_pl(i,j) <= radius_pl(j)) then
                      err_ar(i) = j
                      open(143,file='results/' // trim(file_char) // '_collision.txt', position='append', action='write')
                      write(143,'(2(i0,","),*(g0.16,:,","))') indices(i),j,(t_global+t)/(2._dp*pi),x(3*i-2:3*i), y(3*i-2:3*i),x(3*j-2:3*j), y(3*j-2:3*j)
                      close(143)
                  endif
              enddo
              
              ! save the time when it happened
              t_collision(i) = t
              
              ! put this particle to the end
              x32(3*i-2:npc) = [x32(3*i+1:npc), x32(3*i-2:3*i)]
              y32(3*i-2:npc) = [y32(3*i+1:npc), y32(3*i-2:3*i)]
              x = x32; y = y32
              dist_sun(i:npc/3) = [dist_sun(i+1:npc/3), dist_sun(i)]
              dist_pl(i:npc/3-1,:) = dist_pl(i+1:npc/3,:);
              err_ar(i:npc/3) = [err_ar(i+1:npc/3), err_ar(i)]
              indices(i:npc/3) = [indices(i+1:npc/3), indices(i)]
              t_collision(i:npc/3) = [t_collision(i+1:npc/3), t_collision(i)]
              
              close_enc_inf(i:npc/3-1,:) = close_enc_inf(i+1:npc/3,:)
              if_close(i:npc/3-1,:) = if_close(i+1:npc/3,:)
              if_close_saved(i:npc/3-1,:) = if_close_saved(i+1:npc/3,:)
              
              tmp_a(:,:3) = a(:,3*i-2:3*i); a(:,3*i-2:npc-3) = a(:,3*i+1:npc); a(:,npc-2:npc) = tmp_a(:,:3)
              tmp_a(:,:3) = ax(:,3*i-2:3*i); ax(:,3*i-2:npc-3) = ax(:,3*i+1:npc); ax(:,npc-2:npc) = tmp_a(:,:3)
              npc = npc - 3
              netol = netol - 3
              i = i - 1
          endif
          i = i + 1
          if(i > npc/3) exit
      enddo
      
    !==============================================================================================================
      
5      if(ls) exit ! Exit

      h=r*h; 
      !if(h < 1.d-8) h = 1.d-8
      h2=h**2; us(:npc)=0.0_dp

      ! Last Step

      if((tf-t)/h<1.0_dp) then
        r=r*(tf-t)/h; step=abs(h) 
        h=tf-t; h2=h**2; ls=.true.
      end if

      ! Extrapolation
      
      a_previous = a; ax_previous = ax ! Save a and ax
      
      do i=1,ns
      
      f(:npc)=0.0_dp
      do j=ns,1,-1
        f(:npc)=a(j,:npc)+(1.0_dp+r*c(i)-c(j))*f(:npc)
      end do

      do j=1,i-1
        f(:npc)=(f(:npc)-ax(j,:npc))/d(i,j)
      end do

      ! Correction

      if(nst/=1) da(i,:npc)=a(i,:npc)-ax(i,:npc)
      
      ax(i,:npc)=f(:npc)
      
      end do

      ! Predictor

      a(:,:npc)=ax(:,:npc)+da(:,:npc)

      end do

      end subroutine lobbie2_2
!************************************************************************



!************************************************************************

subroutine fun32_old(n,t,x,y,f)

!*************************************************************************
!*                                                                       *
!*     This subroutine calculates the function f of the differential     *
!*     equation x''=f(x) of the N-body problem.    *
!*                                                                       *
!*************************************************************************
implicit none 
!implicit real*8(a-h,o-z)
real(8), parameter :: k_gauss =  0.017202098949957226d0 !0.01720209895d0
real(8), parameter :: GM = K_gauss**2
integer(4), intent(in) :: n
real(dp), intent(in) :: t
real(dp),dimension(n), intent(in) :: x,y    ! 3*(N_planets+N_particles)
real(dp), dimension(n), intent(out) :: f     ! 3*(N_planets+N_particles)
integer(4) :: i, j, n_part
real(dp), dimension(:,:) :: s2(n/3,N_planets)  ! (N_planets+N_particles,N_planets)
real(dp), dimension(:) :: s1(n/3), fot(3), beta(n/3)   ! (N_planets+N_particles)
real(dp) :: rho = 2000, radius = 50
real(dp), dimension(n/3, 3) :: angul_momentum, transvers_direction
real(dp) :: A2 = - 1.e-13 / GM

n_part = size(x,1)/3 - N_planets

f = 0.0_dp

! first create the array of **3/2 distances, but we use 1/distance**(3/2)
! so that we need to multiply by it and not divide. 
! Therefore just make it zero if it is from planet i to itself

! 1 / distances from the Sun of each object
forall(i=1:N_planets+n_part) s1(i) = 1.0_dp / sqrt(sum(X(3*i-2:3*i)**2))**3

s2 = 0.0_dp
! pairwise distances between planets
do i=2,N_planets
    forall(j=1:i-1) s2(i,j) =  1.0_dp / ( sqrt(sum((x(3*i-2:3*i)-x(3*j-2:3*j))**2))**3 )
enddo
! transpose the matrix
do j = 1, N_planets-1
    s2(j, j+1:) = s2(j+1:, j)
end do
! 1 / pairwise distance between planets and particles
do j=1, N_planets
    forall(i=N_planets+1:N_planets+n_part) s2(i,j) = 1.0_dp / (sqrt(sum((x(3*i-2:3*i)-x(3*j-2:3*j))**2))**3)
enddo


 !loop for coordinates x,y and z [ vectorized ]
do i = 1, 3
! loop for each planet
do j = 1, N_planets
    f(i::3) = f(i::3) + mass_pl(j) * ( (x(3*j-3+i)-x(i::3))*s2(:,j) - x(3*j-3+i)*s1(j))
enddo    
enddo


!! add relativity
!if(if_relativity_include) then
!    do i = 1, N_planets+n_part
!        call RELATIVITY(x(3*i-2:3*i),y(3*i-2:3*i)*k_gauss,fot)
!        f(3*i-2:3*i) = f(3*i-2:3*i) + fot(:) / GM
!    enddo
!endif


!! add Yarkovsii effect (only A2)
!if(if_yarkovskii_include) then
!    do i = N_planets+1,N_planets+n_part
!        ! angul_momentum = r vector mult by V
!        angul_momentum(i,:) = vect_product(x(3*i-2:3*i), y(3*i-2:3*i))
!        
!        ! make it a unit vector
!        angul_momentum(i,:) = angul_momentum(i,:) / sqrt( sum( (angul_momentum(i,:))**2 ) )
!        
!        ! transverse direction = angul_momentum vector mult by r (we don't make it a unit vector)
!        transvers_direction(i,:) = vect_product( angul_momentum(i,:), x(3*i-2:3*i) )
!        
!        ! add the Yarkovsii acceleration ( s1 is 1/r^3, that is why transvers_direction is proportional to r so in the end it is 1/r^2)
!        f(3*i-2:3*i) = f(3*i-2:3*i) + A2 * transvers_direction(i,:) * s1(i)
!    enddo
!endif

beta(:) = 0.0_dp
!! add solar radiation pressure
!if(if_sun_radiation_press_include) then
!    beta(N_planets+1:) = 5.7e-4_dp / (rho*radius)
!else
!    beta(:) = 0.0_dp
!endif


! add the gravitation from the Sun vectorized
f(1::3) = f(1::3) - (1.0_dp - beta(:)) * x(1::3) * s1(:)
f(2::3) = f(2::3) - (1.0_dp - beta(:)) * x(2::3) * s1(:)
f(3::3) = f(3::3) - (1.0_dp - beta(:)) * x(3::3) * s1(:)


contains

function vect_product(a, b) result(c)
  implicit none
  real(8), intent(in) :: a(3), b(3)
  real(8) :: c(3)

  c(:) = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function vect_product

end subroutine fun32_old



subroutine fun32(n,t,x,y,f)

!*************************************************************************
!*                                                                       *
!*     This subroutine calculates the function f of the differential     *
!*     equation x''=f(x) of the N-body problem.    *
!*                                                                       *
!*************************************************************************
implicit none 
!implicit real*8(a-h,o-z)
real(8), parameter :: k_gauss =  0.017202098949957226d0 !0.01720209895d0
real(8), parameter :: GM = K_gauss**2
integer(4), intent(in) :: n
real(dp), intent(in) :: t
real(dp),dimension(n), intent(in) :: x,y    ! 3*(N_planets+N_particles)
real(dp), dimension(n), intent(out) :: f     ! 3*(N_planets+N_particles)
real(dp), dimension(n/3,3) :: xc, fc
real(dp), dimension(n/3,N_planets) :: dx, dy, dz
integer(4) :: i, j, n_part
real(dp), dimension(:,:) :: s2(n/3,N_planets)  ! (N_planets+N_particles,N_planets)
real(dp), dimension(:) :: s1(n/3), fot(3), beta(n/3)   ! (N_planets+N_particles)
real(dp) :: rho = 2000, radius = 50
real(dp), dimension(n/3, 3) :: angul_momentum, transvers_direction
real(dp) :: A2 = - 1.e-13 / GM

n_part = size(x,1)/3 - N_planets

f = 0.0_dp

do i = 1, 3
    xc(:,i) = x(i::3)
enddo

fc = 0.0_dp

! first create the array of **3 distances, but we use 1/distance**(3)
! so that we need to multiply by it and not divide. 
! Therefore just make it zero if it is from planet i to itself

! 1 / distances^3 from the Sun of each object
s1(:) = 1.0_dp / sqrt( xc(:,1)**2 + xc(:,2)**2 + xc(:,3)**2  )**3

!--------------------------------------------------
! Pairwise planet-body differences (vectorized in i)
!--------------------------------------------------
do j = 1, N_planets

    dx(:,j) = xc(j,1) - xc(:,1)
    dy(:,j) = xc(j,2) - xc(:,2)
    dz(:,j) = xc(j,3) - xc(:,3)

    s2(:,j) = 1.0_dp / sqrt(dx(:,j)**2 + dy(:,j)**2 + dz(:,j)**2)**3

    ! avoid self singularity 
    s2(j,j) = 0.0_dp
end do


! loop for each planet
do j = 1, N_planets  
    fc(:,1) = fc(:,1) + mass_pl(j) * ( dx(:,j)*s2(:,j) - xc(j,1)*s1(j))
    fc(:,2) = fc(:,2) + mass_pl(j) * ( dy(:,j)*s2(:,j) - xc(j,2)*s1(j))
    fc(:,3) = fc(:,3) + mass_pl(j) * ( dz(:,j)*s2(:,j) - xc(j,3)*s1(j))
enddo    


!! add relativity
!if(if_relativity_include) then
!    do i = 1, N_planets+n_part
!        call RELATIVITY(x(3*i-2:3*i),y(3*i-2:3*i)*k_gauss,fot)
!        !f(3*i-2:3*i) = f(3*i-2:3*i) + fot(:) / GM
!        fc(i,:) = fc(i,:) + fot(:) / GM
!    enddo
!endif


!! add Yarkovsii effect (only A2)
!if(if_yarkovskii_include) then
!    do i = N_planets+1,N_planets+n_part
!        ! angul_momentum = r vector mult by V
!        angul_momentum(i,:) = vect_product(x(3*i-2:3*i), y(3*i-2:3*i))
!        
!        ! make it a unit vector
!        angul_momentum(i,:) = angul_momentum(i,:) / sqrt( sum( (angul_momentum(i,:))**2 ) )
!        
!        ! transverse direction = angul_momentum vector mult by r (we don't make it a unit vector)
!        transvers_direction(i,:) = vect_product( angul_momentum(i,:), x(3*i-2:3*i) )
!        
!        ! add the Yarkovsii acceleration ( s1 is 1/r^3, that is why transvers_direction is proportional to r so in the end it is 1/r^2)
!        !f(3*i-2:3*i) = f(3*i-2:3*i) + A2 * transvers_direction(i,:) * s1(i)
!        fc(i,:) = f(i,:) + A2 * transvers_direction(i,:) * s1(i)
!    enddo
!endif

beta(:) = 0.0_dp
!! add solar radiation pressure
!if(if_sun_radiation_press_include) then
!    beta(N_planets+1:) = 5.7e-4_dp / (rho*radius)
!else
!    beta(:) = 0.0_dp
!endif


! add the gravitation from the Sun vectorized
fc(:,1) = fc(:,1) - (1.0_dp - beta(:)) * xc(:,1) * s1(:)
fc(:,2) = fc(:,2) - (1.0_dp - beta(:)) * xc(:,2) * s1(:)
fc(:,3) = fc(:,3) - (1.0_dp - beta(:)) * xc(:,3) * s1(:)


do i = 1, 3
    f(i::3) = fc(:,i)
enddo

contains

function vect_product(a, b) result(c)
  implicit none
  real(8), intent(in) :: a(3), b(3)
  real(8) :: c(3)

  c(:) = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function vect_product

end subroutine fun32

      
!-----------------------------------------------------------------------
!   Ýòà ïîäïðîãðàììà ó÷èòûâàåò ðåëÿòèâèñòñêèå ýôôåêòû íà îñíîâå òåîðèè
!   îòíîñèòåëüíîñòè.
!
!   Âõîäíûå ïàðàìåòðû:
!
!   1) X - ìàññèâ êîîðäèíàò.
!   2) Y - ìàññèâ ñêîðîñòåé
!   3) FOT - ìàññèâ ïîïðàâîê çà ðåëÿòèâèñòñêèå ýôôåêòû 
!-----------------------------------------------------------------------

SUBROUTINE RELATIVITY(X,V,FOT)

IMPLICIT REAL*8 (A - H, O - Z)

DIMENSION X(*), V(*), FOT(*)

AMU = 9.8709893D-9
GAUSS_SQUARED = 0.000295912208286D0   !Ãàóññ **2
!XV = X(1) * V(1) + X(2) * V(2) + X(3) * V(3)
!R = DSQRT(X(1)**2 + X(2)**2 + X(3)**2)
!V2 = V(1)**2 + V(2)**2 + V(3)**2

XV = sum(X(1:3)*V(1:3))
R = dsqrt(sum(X(1:3)**2))
V2 = sum(V(1:3)**2)

FOT(1:3) = AMU * (4.0D0 * GAUSS_SQUARED * X(1:3) / R - V2 * X(1:3) + 4.0D0 * XV * V(1:3)) / R**3


RETURN
END SUBROUTINE RELATIVITY



function kahan_sum(arr) result(sum_result)
    implicit none
    real(dp), intent(in) :: arr(:)   ! Input array of any size
    real(dp) :: sum_result           ! Final result of the summation
    real(dp) :: sum, compensation    ! Variables for Kahan summation
    real(dp) :: temp, y              ! Temporary variables for calculation
    integer(4) :: i, n

    n = size(arr)
    if (n == 0) then
        sum_result = 0.0_dp
        return
    end if

    sum = 0.0_dp
    compensation = 0.0_dp

    ! Loop over the array using Kahan summation
    do i = 1, n
        y = arr(i) - compensation
        temp = sum + y
        compensation = (temp - sum) - y
        sum = temp
    end do

    sum_result = sum
end function kahan_sum

! vectorized version of kahan summation
! arr is a matrix m*n. The result is a vector of dimension m
! NOTE that we add by lines and the result is the dimension of the number of lines!
subroutine kahan_sum_vector(arr, sum)
    implicit none
    real(dp), dimension(:,:), intent(in) :: arr  ! Input: array of vectors (size [m, n])
    real(dp), dimension(:), intent(out) :: sum ! Output: result vector (size m)
    integer(4) :: i !, n, m
    real(dp), dimension(size(arr, 1)) :: compensation, y, t
    
    !! Get the dimensions of the input array
    !n = size(arr, 2)  ! Number of elements in the first dimension
    !m = size(arr, 1)  ! Length of each vector
    
    ! Initialize variables
    sum = 0.0_dp
    compensation = 0.0_dp
    
    ! Perform Kahan summation for each vector component
    do i = 1, size(arr, 2)
        y = arr(:,i) - compensation    ! Apply the compensation
        t = sum + y                    ! Temporarily add the compensated value
        compensation = (t - sum) - y   ! Calculate the new compensation
        sum = t                        ! Update the sum
    end do
    
    !! Return the result
    !result = sum

end subroutine kahan_sum_vector


end module lobbie2_mod