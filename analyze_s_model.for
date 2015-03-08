      module thermo_data
c
c	Fundamental constants are from Novotny, Introduction to
c    Stellar Atmospheres and Interiors, 1973, Appendix II.
c						  2/16/83  WD Pesnell
c
      real*8, parameter :: zero = 0.d0
      real*8, parameter :: one = 1.d0
      real*8, parameter :: two = 2.d0
      real*8, parameter :: thre = 3.d0
      real*8, parameter :: four = 4.d0
      real*8, parameter :: ten = 10.d0
      real*8, parameter :: ahf = 0.5d0
      real*8, parameter :: qrt = 0.25d0
      real*8, parameter :: grav = 6.689067384d-8
      real*8, parameter :: ac3 = 7.5595D-5
      real*8, parameter :: sigma = 5.66961D-5
      real*8, parameter :: pi = 3.1415926535897932D0
      real*8, parameter :: pi2 = 2.d0*pi !6.2831853072D0
      real*8, parameter :: pi4 = 4.d0*pi !12.5663706144D0
      real*8, parameter :: pi8 = 8.d0*pi !25.1327412288D0
      real*8, parameter :: pi43 = 4.d0*pi/3.d0 !4.1887902048D0
c
      real*8, parameter :: L_solar = 3.846d33    ! ergs/s
      real*8, parameter :: M_solar = 1.98855d33   ! grams
      real*8, parameter :: R_solar = 6.955d10    ! cm
c
      end

      program analyze_s_model
      use thermo_data
      implicit real*8(a-h,o-z)
      real*8 nu_gt
      character*32 plot_filename
c
c	   Polytrope integrator. Integrates the second order equation
c	governing polytropes using a Runge-Kutta method and calculates
c       the linear, adiabatic eigenvalues of the resulting model. The
c       adiabatic index is set to 5/3 in this version.
c
      parameter ( nmax=4096 )
      integer g_max(1)
      common/phypar/ r(nmax),T_K(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax), c_s(nmax)
c      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/corevl/ rz0,p0,rho0
      common/stnmod/ nmode
      dimension c2(nmax),sl2(nmax)
c
      open (unit=1,file='model_s.sum',status='unknown')
c      open (unit=5,file='sys$input',status='unknown')
c      open (unit=6,file='sys$output',status='unknown')
      open (unit=11,file='model_s.out',status='unknown')
c
c          read in the input variables from terminal.
c
c          iomod := output flag. Iout greater than unity produces
c                 voluminous lines of output.
c
c          npts := number of points in model. Npts must be less than
c               nmax-2, which is the dimension of the matrices
c               used.
c
c          fnpol := polytropic index of the initial model. Default
c                  value for fnpol .lt. 0 Is 3.
c
c          r0 := initial guess at the radius of the model. This
c               parameter is found as an eigenvalue in the solution
c               and a good guess speeds the convergence. Default
c               value is 5 (used if the inputted value is .le. 0).
c
c      write(*,*) grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
   1  continue
      write(6,1000)
c      read(5,*,end=910) iomod,npts,fnpol,r0
c      if( r0 .le. zero ) r0 = 5.d0
c      if( fnpol .lt. zero ) fnpol = thre
c      if( npts+2 .gt. nmax ) npts = nmax - 2
c
c          12/28/2002: Added a unique filename for the plot file
c
      fnpol = 3.d0
c      write(plot_filename, 1010) int( 1000.*fnpol ), npts
      totmas = 1.991e33
      write(plot_filename, 1010) int( 1000.*totmas/M_solar ), npts
      close(12)
      open (unit=12,file=plot_filename,status='unknown')
c
      idonpl = 1
      call read_S_model(npts)
c
c       Core conditions are in the first point
c
      pc = p(1)
      rz0 = ahf*xi(2)*R_solar
      rho0 = rho(1)
      Tc = T_K(1)
      p0 = pc
      npuls = npts-1
c
c       thermodynamic quantities are already centered
c
      rzone(1:npuls) = xi(2:npts)*R_solar
      c2 = c_s**2
      do 10 i=1,npuls
         r(i) = (rzone(i-1)+rzone(i))/2.d0
         v(i) = one/rho(i)
         P_rad = 7.56471d-15*T_K(i)**4/3.d0
c	 write(*,*) T_K(i), P_rad/P(i)
  10  continue
      r(npts) = r(npts-1) + (rzone(npts-1)-rzone(npts-2))
c      write(*,*) r(npts), r(npts-1), rzone(npts-1), rzone(npts-2)
c
      asymfr = zero
      do i=1,npuls
         asymfr = asymfr + (r(i+1)-r(i))/c_s(i)
      enddo
      asymfr = one/asymfr
c      bvfac = fnpol - (fnpol+one)/g1(1)
c      om02 = pi*g*rhobar
      dm1(1) = 4.d0*pi/3.d0*rho(1)*r(1)**3
      tmas = dm1(1)
      do 30 i=2, npts
         dm1(i) = (r(i+1)**3-r(i)**3)*pi43*rho(i)
         tmas = tmas + dm1(i)
  30  continue
      dm2(1) = ahf*(dm1(1))
      do 40 i=1,npts
         dm2(i) = ahf*(dm1(i)+dm1(i-1))
  40  continue
      dm2(npuls+1) = ahf*dm1(npuls)
      rm(1) = dm1(1)
      gor(1) = grav*rm(1)/r(1)**2
      do 15 i=2, npts+1
         rm(i) = rm(i-1) + dm1(i)
         gor(i) = grav*rm(i)/r(i)**2
c         if( i .eq. npuls + 1 ) goto 15
c         if( dabs(bvfac) .gt. 1.d-12 ) then
c            rai = bvfac*dthdr(i+1)/(T_K(i)*xi(i+1))
c            ra(i) = rai
c            bv(i) = dlog10((dabs(rai)*gor(i)/r(i))/om02)
c         endif
  15  continue
      rmch = (rm(npuls+1)-tmas)/rm(npuls+1)
      rhobar = tmas/pi43/R_solar**3
c
      g_max = maxloc( gor(1:npuls) )
c      write(*,*) g_max(1)
c      igmax = 0
c      do 16 i=2,npuls
c         if( gor(i+1) .lt. gor(i) ) goto 17
c          igmax = i+1
c  16  continue
c  17  continue
      igmax = g_max(1)
      write(6,1700) igmax,xi(igmax), rm(igmax)/rm(npuls+1)
      omega_ac = c_s(npts)/(P(npts)/rho(npts)/gor(npts))/2.
      write(6,1099) asymfr*0.5e6, omega_ac/pi2*1.e6
c
c      call pltdmp(bv,nmax,npuls,'bv  ')
      call pltdmp(gor,nmax,npuls,'g   ')
      call pltdmp(c2,nmax,npuls,'c2  ')
      do 20 i=1,npuls
         sl2(i) = 6.d0*c2(i)/rzone(i)**2
  20  continue
      call pltdmp(sl2,nmax,npuls,'sl2 ')
      write(11,2000) npts
      write(11,2003) fnpol,r0,alfa,rlam,rhobar,rm(npuls+1),tmas,rmch
      write(11,2006) rz0,p0,rho0
      write(11,1700) igmax,xi(igmax),rm(igmax)/rm(npuls+1)
      write(11,1099) asymfr*0.5e6, omega_ac/pi2*1.e6
      write(11,2004)
      write(11,2005) (i,xi(i),T_K(i),c_s(i),r(i),p(i),rho(i),
     $   gor(i),rm(i),dm1(i),rzone(i),i=1,npts)
c
      qch = rlam/rhobar
      write(1,2007) npuls,fnpol,r0,alfa,rlam,rhobar,rm(npuls+1),
     $      tmas,rmch
      write(1,2008) rz0,p0,rho0,qch
      write(1,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
      write(1,1099) asymfr*0.5e6, omega_ac/pi2*1.e6
c
c          pulsation analysis
c
c          enter <crtl/z> to stop program.
c
c          iopuls := output flag. Iout equal to unity produces
c                   convergence information on file 11, iopuls
c		    equal to 2 prints conv. info., eigenvectors,
c                   and matrices.
c
c	   lval := l-value of the pulsation, l < 0 stops program.
c
c          nomega := number of points in scan calculation (l>0),
c                   lowest order mode in the radial case.
c
c          ihmax := highest order mode calculated in the l=0 case,
c                  must be greater than ihmin.
c
c          omlow := shortest omsq value used in scan calculation (l>0),
c                   not used in radial case.
c
c          omhigh := largest omsq value used in scan calculation (l>0),
c                  not used in radial case.
c	   (both omlow and omhigh are entered in units of pi*g*rhobar.)
c
 100  continue
         write(6,1200,advance='no')
         read(5,*,end=910) iopuls,lval,nomega,ihmax !,omlow,omhigh
         if( lval .lt. 0 ) then
            goto 910
         else if( lval .eq. 0 ) then
            nmode = 0
	    omlow = 0.
	    omhigh = 0.
          call lnanon(npuls,iopuls,lval,nomega,ihmax,omlow,omhigh,fnpol)
         else
	    call cpu_time( start_time )
            do lval=2, 800
c               omlow = (pi2*5.d-6*(float(lval)/2.))**2
c               omhigh = (pi2*2.d-5*(1000. + float(lval)/2.))**2
c
c       Asympototic fits from model runs, 09-Oct-2014
c
               omlow = (pi2*(7.2d-5*float(lval)**0.4805d0))**2
	       nu_gt = 8.79d-3 + 1.7d-5*(100. + float(lval)/2.+0.25d0)
               omhigh = (pi2*nu_gt)**2
            nmode = 0
          call lnanon(npuls,iopuls,lval,nomega,ihmax,omlow,omhigh,fnpol)
            enddo
	    call cpu_time( end_time )
      write(6,*) end_time - start_time, ' seconds of CPU time'
         endif
      goto 100
 900  continue
      goto 1
 910  continue
      stop
c
 1000 format(1x,'Analyzing Model S')
 1010 format('Model_S_',i4.4,'_'i5.5,'.plt')
 1200 format(1x,'enter pulsation parameters:',/,
     $  1x,' io, l, nom, and ihmax: ')
c     $  1x,' io, l, nom, ihmax, omlow, and omhigh: ')
 1700 format(1x,'max. g in zone ',i4,' radius = ',1pe10.3,
     $       ' mass = ',e10.3)
 1099 format(1x,'asymptotic frequency for p modes is ',1pe11.4,
     $       ' muHz', /
     $       1x,'Acoustic cutoff frequency is ', 1pe11.4,' muHz')
 2000 format('S Model of the current Sun (Christensen-Dalsgaard)',//,
     $      1x,i5,' zones',/)
 2003 format(7h fnpol=,f6.2,26H guess for surface radius=,1pe10.3,
     $   7H alpha=,e10.3,/,2X,8h lambda=,e10.3,8H rhobar=,e10.3,
     $   8H totmas=,e10.3,11H int. Mass=,e10.3,7H diff.=,E10.3)
 2004 format(3x,1hi,5x,2hxi,7x,'T (K)',5x,' c_s ',4x,6hradius,6x,1hp,
     $   8x,3hrho,8x,1hg,9x,2hrm,7x,3hdm1,7x,5hrzone)
 2005 format(1x,i4,1p,10e10.3)
 2006 format(1x,5h rz0=,1pe10.3,4H p0=,e10.3,6H rho0=,e10.3)
 2007 format(1h1,/,1x,18houtput from polyrk,i5,6h zones,/,
     $   7h fnpol=,f6.2,26H guess for surface radius=,1pe10.3,
     $   7H alpha=,e10.3,/,2X,8h lambda=,e10.3,8H rhobar=,e10.3,/,
     $   8H totmas=,e10.3,11H int. Mass=,e10.3,7H diff.=,E10.3)
 2008 format(1x,5h rz0=,1pe10.3,4H p0=,e10.3,6H rho0=,e10.3,/
     $   1X,29h central/average densities = ,e10.3)
      end
      block data
      implicit real*8(a-h,o-z)
c
c	Fundamental constants are from Novotny, Introduction to
c    Stellar Atmospheres and Interiors, 1973, Appendix II.
c						  2/16/83  WD Pesnell
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
c      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726D-8,7.5595D-5,
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.689067384e-8,7.5595D-5,
     $  5.66961D-5,3.1415926536D0,6.2831853072D0,12.5663706144D0,
     $   25.1327412288D0,4.1887902048D0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,10.0D0,
     $      0.5D0,0.25D0      /
      end
      subroutine read_s_model( npts )
      implicit real*8(a-h,o-z)
      character*80 ititl, a
c
c          read the S model.
c
      parameter ( nmax=4096 )
      common/phypar/ r(nmax),T_K(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax), c_s(nmax)
c      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      lunit = 50
      open(lunit,
     $  file='../../solar/evolution/cptrho.l5bi.d.15c_Solar_Model.txt',
     $   action='read')
      read( lunit, *) ititl
      do i=0, 2
         read( lunit, *) a
c         write(*,*) a
      enddo
c
      npts = 2482
c
c       Columns are r/R_s, c_s, rho, P, G1, T(i)
c
      do i=npts, 1, -1
         read(lunit,*) xi(i), c_s(i), rho(i), P(i), G1(i), T_K(i)
c         write(*,*) xi(i), c_s(i), rho(i), P(i), G1(i), T_K(i)
      enddo
      close(lunit)
c
      idonpl = 1
      call pltint(nmax,npts, 0., xi(npts), idonpl)
      write(*,2000) npts
 2000 format(1x,'Read the Model S with ', i6, ' shells')
c
      return
      end
      subroutine pltint(nmaxin,n,fnpol,r0,idonpl)
      implicit real*8(a-h,o-z)
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=4096 )
      common/phypar/ r(nmax),T_K(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax), c_s(nmax)
c      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      if( idonpl .gt. 0 ) then
         idonpl =-1
         write(12,1000) r0
      endif
      call pltdmp(xi(2),nmaxin,n,'xi  ')
      do  20 i=2,n
         r(i) = xi(i)/xi(n)
  20  continue
      call pltdmp(r(2),nmaxin,n-1,'x   ')
      return
c 1000 format(1x,18hpolytrope with n =,f5.2,20H^h2}x^lh.5%0 ^Lxhx%=,f7.4)
 1000 format(1x,'Solar Model with outer radius= ',f7.4)
      end

      subroutine lnanon(nin,iout,lin,nomega,ihmax,omlow,omhigh,fnpoln)
      implicit real*8(a-h,o-z)
c
      type puls_sum
         integer*2 :: mass_flag
         integer*4 :: ell, k, n_p, n_g
         real*4 :: omega_ad, omsq_n, nu_ad, period, acc_ad_w, acc_ad_e
	 real*4 :: rbar
         complex*8 :: del_om_qa
         complex*8 :: comega
         real*4 :: kappa, acc_na_k
      end type
      type(puls_sum)::nonradial
c
      parameter ( nmax=4096 )
      parameter ( nmax3=3*nmax )
      common/phypar/ r(nmax),T_K(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax), c_s(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk37/  drdm(nmax,2)
      common/const/  zero,one,two,thre,four,ten,ahf,qrt
      common/corevl/ rz0,p0,rho0
      common/scrtch/ ag1(nmax,3),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),ap4(nmax)
      common/linear/ wtherm(nmax),wgrav(nmax),wcross(nmax),wdiag(nmax),
     $      dr(nmax),dh(nmax),gam(nmax),strke(nmax),z(nmax),
     $      dp(nmax),adrho(nmax),weight(nmax),xo(nmax),
     $      yo(nmax),stwait(nmax),epstwt(nmax),spac(nmax,5)
      dimension drz(nmax),rkl2(nmax),drint(nmax)
      dimension x(nmax3),stomeg(200), stomsq(200)
      dimension st_k_p(200), st_k_g(200)
      data accur/1.d-11/
c
c          zero the eigenvector matrix to avoid indefinite operands.
c
      do 5 i=1,nmax3
         x(i) = zero
   5  continue
c
c          define the differences of the two radii (zone
c       and interface centered) for ease of future work. The
c       array gor is the local gravity and rho is the zone density.
c
      n = nin
      np1 = n + 1
      lval = lin
      rl = dfloat(lval)
      rl1 = rl*(rl+one)
      fnpol = fnpoln
      rhom = rm(np1)/(pi43*r(np1)**3)
      do 10 i=1,n
         rkl2(i) = rl1/rzone(i)**2
         if( i.gt.1 ) drz(i) = rzone(i)-rzone(i-1)
         drint(i) = r(i+1)-r(i)
         drdm(i,1) = forpi*rzone(i)/(dm1(i)*v(i)) -
     $    ahf*(rl+one)/rzone(i)**2
         drdm(i,2) =-(forpi*rzone(i)/(dm1(i)*v(i)) +
     $    ahf*(rl+one)/rzone(i)**2)
         z(i) = rzone(i)/r(np1)
  10  continue
c
      call pltdmp(rzone,nmax,n,'rz  ')
      drz(1) = rzone(1) - rz0
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      drz(np1) = p(n)*v(n)/gor(np1)
      pi4g = forpi*g
      if( lval .le. 0 ) goto 50
c
c          Initialize the matrices.
c
c          Inner boundary conditions.
c
      afac = forpi*r(1)**3/dm2(1)
      drdm01 =-ahf*(rl+one)/rz0**2
      drdm02 =-ahf*(rl+one)/rz0**2
      ag11 = drdm01*g1(1)*p0*(ahf*rl/rho0 - afac)
      ag1(1,1) = zero
      ag1(1,2) =-gor(1)/r(1) +
     $     drdm(1,1)*g1(1)*p(1)*(ahf*rl*v(1) + afac) +
     $     drdm02*g1(1)*p0*(ahf*rl/rho0 - afac )
      ag1(1,3) = drdm(1,2)*g1(1)*p(1)*(ahf*rl*v(1)+afac)
      ag1(1,3) = ag1(1,3) + ag11
      ag3(1,1) = zero
      ag3(1,2) = rl
      ag41 = rl1*(g1(1)*p0*(ahf*rl/rho0-afac))/rz0**2
      ag4(1,1) = zero
      ag4(1,2) = rkl2(1)*(g1(1)*p(1)*(ahf*rl*v(1)+afac)) +
     $     rl1*gor(1)/r(1)
      ag4(1,2) = ag4(1,2) + ag41
c
      ap1(1,1) =-pi4g*rho(1)*rzone(1)**2*drdm(1,1)
      ap1(1,2) =-pi4g*rho(1)*(rzone(1)**2*drdm(1,2)-
     $     ahf*dlog(rho(2)/rho(1))/dlog(rzone(2)/rzone(1)) )
      ap311 = r(1)**2/(drint(1)*drz(1))
      ap313 = r(2)**2/(drint(1)*drz(2))
      ap3(1,1) = zero
      ap3(1,2) =-ap313 - ap311*(drz(1)/drz(2))**2
      ap3(1,3) = ap313 + ap311*(drz(1)/drz(2))**2
      ap4(1)   =-pi4g*rl1*rho(1)
c
      do 20 i=1,n
c
c          horizontal momentum equation.
c
         gzone = ahf*(gor(i)+gor(i+1))/rzone(i)
         ah1(i,1) = drdm(i,1)*g1(i)*p(i)*v(i) + ahf*gzone
         ah1(i,2) = drdm(i,2)*g1(i)*p(i)*v(i) + ahf*gzone
         ah3(i) = one
         ah4(i) = rkl2(i)*g1(i)*p(i)*v(i)
         if( i .eq. n ) goto 20
c
c          radial momentum equation.
c
         afac = forpi*r(i+1)**3/dm2(i+1)
         ag1(i+1,1) = drdm(i,1)*g1(i)*p(i)*(ahf*rl*v(i) - afac)
         ag1(i+1,2) =-four*gor(i+1)/r(i+1)+pi4g*ahf*(rho(i+1)+rho(i))+
     $     drdm(i+1,1)*g1(i+1)*p(i+1)*(ahf*rl*v(i+1) + afac) +
     $     drdm(i,2)*g1(i)*p(i)*(ahf*rl*v(i) - afac )
         ag1(i+1,3) = drdm(i+1,2)*g1(i+1)*p(i+1)*(ahf*rl*v(i+1)+afac)
         ag3(i+1,1) = ahf*rl - r(i+1)/drz(i+1)
         ag3(i+1,2) = ahf*rl + r(i+1)/drz(i+1)
         ag4(i+1,1) = rkl2(i)*(g1(i)*p(i)*(ahf*rl*v(i)-afac)) +
     $     ahf*rl1*gor(i+1)/r(i+1)
         ag4(i+1,2) = rkl2(i+1)*(g1(i+1)*p(i+1)*(ahf*rl*v(i+1)+afac)) +
     $     ahf*rl1*gor(i+1)/r(i+1)
         if( i .eq. 1 ) goto 20
c
c          poisson equation.
c
         ddr11 =-drz(i+1)/(drz(i)*(drz(i)+drz(i+1)))
         ddr12 = (drz(i+1)-drz(i))/(drz(i)*drz(i+1))
         ddr13 = drz(i)/(drz(i+1)*(drz(i)+drz(i+1)))
         ap1(i,1) =-pi4g*rho(i)*(rzone(i)**2*drdm(i,1)-
     $     ahf*dlog(rho(i)/rho(i-1))/dlog(rzone(i)/rzone(i-1)) )
         ap1(i,2) =-pi4g*rho(i)*(rzone(i)**2*drdm(i,2)-
     $     ahf*dlog(rho(i+1)/rho(i))/dlog(rzone(i+1)/rzone(i)) )
         ap3i1 = r(i)**2/(drint(i)*drz(i))
         ap3i3 = r(i+1)**2/(drint(i)*drz(i+1))
         ap3(i,1) = ap3i1 +           ddr11*two*rl*rzone(i)
         ap3(i,2) =-(ap3i1 + ap3i3) + ddr12*two*rl*rzone(i)
         ap3(i,3) = ap3i3 +           ddr13*two*rl*rzone(i)
         ap4(i) =-pi4g*rl1*rho(i)
  20  continue
c
c          outer boundary conditions.
c
      dlnr = (rl+ahf)*drz(np1)/r(np1)
      f = (drz(np1)/r(np1))/(one+dlnr)
      vl =-g1(n)*gor(np1)*r(np1)
      vl =-g1(n)*p(n)*forpi*r(np1)**3/dm2(np1)
      ag1(np1,1) = vl*drdm(n,1)
c
c       12/20/99: Changing FDR to 1 had increased DRHO near the surface by
c                 a small amount. Eigenvalue stayed the same.
c
      fdr = (one + ahf*rl*drz(np1)/r(np1))/(one + dlnr)
      ag1(np1,2) =-four*gor(np1)/r(np1) + vl*drdm(n,2) +
     $     pi4g*rho(n)*ahf*fdr
      ag1(np1,3) = zero
      ag3(np1,1) =-(rl + one)/(one+dlnr)
      ag3(np1,2) = zero
      ag4(np1,1) = rkl2(n)*vl + rl1*gor(np1)/r(np1)
      ag4(np1,2) = zero
c
c          set up the outer boundary condition for the
c       poisson equation.
c
      ddr11 =-drz(np1)/(drz(n)*(drz(n)+drz(np1)))
      ddr12 = (drz(np1)-drz(n))/(drz(n)*drz(np1))
      ddr13 = drz(n)/(drz(np1)*(drz(n)+drz(np1)))
      ap3n1 = r(n)**2/(drint(n)*drz(n))
      ap3n3 = r(np1)**2/(drint(n)*drz(np1))
c
      ap1(n,1) =-pi4g*rho(n)*(rzone(n)**2*drdm(n,1)-
     $   ahf*dlog(rho(n)/rho(n-1))/dlog(rzone(n)/rzone(n-1)) )
      dlrhdr =-(rho(n)*gor(n+1)*r(n+1))/p(n)*fnpol/(one+fnpol)
      ap1(n,2) =-pi4g*rho(n)*( (ddr13*two*rl*rzone(n) + ap3n3)*f +
     $    rzone(n)**2*drdm(n,2) - dlrhdr )
c
      ap3(n,1) = ap3n1 + ddr11*two*rl*rzone(n)
      fp3 = (rl+ahf)*(drz(np1)/r(np1)) - 1
      ap3(n,2) =-(ap3n1+ap3n3) + ddr12*two*rl*rzone(n) -
     $     (ddr13*two*rl*rzone(n) + ap3n3)*fp3/(one+dlnr)
      ap3(n,3) = zero
c
c       12/20/99: Dividing AP4(N) by 2 decreased DRHO near the surface by
c                 a small amount. Eigenvalue increased as well.
c
      ap4(n) =-pi4g*rl1*rho(n)
      ah4(np1) = zero
c
c        write the matrices to file 6.
c
      if( iout .ge. 2 ) then
         write(11,2200)
         do 27 i=1,np1
            write(11,2201) i,ag1(i,1),ag1(i,2),ag1(i,3),drdm(i,1),
     $         drdm(i,2),ag3(i,1),ag3(i,2),ag4(i,1),ag4(i,2),ah4(i)
  27     continue
         write(11,2300)
         do 28 i=1,n
            write(11,2201) i,ah1(i,1),ah1(i,2),ah3(i),ap1(i,1),
     $         ap1(i,2),ap3(i,1),ap3(i,2),ap3(i,3),ap4(i)
  28     continue
      endif
c
c	   Convergence to a normal mode is governed by the momentum
c       equation for the outermost zone. As in the radial case, this
c       function (divided by any x value) is used in a secant method
c       for controlling the iterations. As the equations are not
c       nonlinear (as in non-non) in the frequency, it is hoped the
c       calculation is, in some sense, more stable.
c
c						10/4/83 WD Pesnell
c
      write(11,3001)
      xnorm = r(np1)**2
      iroot = 0
      write(*,*) omlow, omhigh, nomega
      dom = (omhigh - omlow)/dfloat(nomega - 1)
      do 30 i=1,nomega
         omsq = omlow + dfloat(i-1)*dom
c         omsq = (omlow + dfloat(i-1)*dom)**2
         fomeg = fomsq(omsq,x,xnorm,n)
         if( i .eq. 1 ) then
            fomi = fomeg
         else
            if( fomeg*fomi .lt. zero ) then
               iroot = iroot + 1
               stomeg(iroot) = omsq - dom*fomeg/(fomeg-fomi)
            endif
            if( iout .gt. 0 ) write(11,4400) i,omsq,fomeg
            fomi = fomeg
         endif
  30  continue
      if( iroot .eq. 0 ) goto 900
      write(11,3000) iroot,lval,(stomeg(i),i=1,iroot)
      write(1,3000) iroot,lval,(stomeg(i),i=1,iroot)
      write(6,3000) iroot,lval,(stomeg(i),i=1,iroot)
c
c          attempt to converge on each of the frequencies found in
c       discriminant search. A secant method to find the zeroes of
c       error function is used.
c
      do 100 nroot=1,iroot
         omsq = stomeg(nroot)
         do 45 icount=1,30
            err = fomsq(omsq,x,xnorm,n)
            if( icount .eq. 1 ) then
               omsq1 = omsq
               erra  = err
               omsq  = omsq*(one + 1.d-7)
            else
               afac = dabs((omsq-omsq1)/omsq)
               if(iout.gt.0) Write(11,4400) icount,omsq,afac,err,erra
               if( afac .le. Accur ) goto 47
               omsq2 = (erra*omsq-err*omsq1)/(erra-err)
               omsq1 = omsq
               erra  = err
               omsq  = omsq2
            endif
  45     continue
         write(6,4500) nroot,omsq
         write(1,4500) nroot,omsq
         write(11,4500) nroot,omsq
         goto 100
c
c          converged to a root, calculate the eigenvectors.
c
  47     continue
         x(3*n+1) = xnorm
         do 40 i=1,n
c            write(11,*) i, x(3*i-2), (r(i)/r(np1))**lval
            dr(i) = (x(3*i-2)/r(i)**2)*(r(i)/r(np1))**lval
            dh(i) = (x(3*i)/rzone(i)**2)*(rzone(i)/r(np1))**lval
            adrho(i) = (drdm(i,1)*x(3*i-2)+drdm(i,2)*x(3*i+1) +
     $               rkl2(i)*x(3*i) )*(rzone(i)/r(np1))**lval
            dp(i) = g1(i)*adrho(i)
            gam(i) = x(3*i-1)*(rzone(i)/r(np1))**lval/(gor(i)*r(i))
  40     continue
         eta0 = dr(1)/rl
         gam0 = x(2)*(rz0/r(np1))**lval
         drho0 = adrho(1)*(rz0/rzone(1))**lval
         dp0 = drho0*g1(1)
         dr(n+1) = one
         drtst = (ag1(n+1,1)*x(3*n-2)+ag3(n+1,1)*x(3*n-1)+
     $      ag4(n+1,1)*x(3*n))/(omsq-ag1(n+1,2))
c         write(6,*) drtst,xnorm,dabs((drtst-xnorm)/xnorm)
         dlnr = ahf*(rl+one)*drz(np1)/r(np1)
         gam(n+1) = gam(n)*(one-dlnr)/(one+dlnr)
         dh(n+1) = zero
         qdr0 = dr(2)*(r(1)/r(2))**(lval-2)
c
c          weight function for the adiabatic oscillation.
c
         wgrav(1) = dr(1)*dm2(1)*rl*gam(1)
         wcross(1) = two*rl1*gor(1)*r(1)*dr(1)*dh(1)*dm2(1)
         wdiag(1) = (pi4g*rho(1)-gor(1)*four/r(1) )*
     $               (dr(1)*r(1))**2*dm2(1)
         stwait(1) =  wgrav(1) + wcross(1) + wdiag(1)
c
         rsum = (dr(1)*r(1))**2*dm2(1)
         hsum = zero
         rbar = r(1)*(dr(1)*r(1))**2*dm2(1)
         hbar = zero
         crsum = zero 
         do 48 i=1,n
            wtherm(i) = g1(i)*p(i)*v(i)*dm1(i)*adrho(i)**2
            rsum = rsum + (dr(i+1)*r(i+1))**2*dm2(i+1)
            hsum = hsum + (dh(i)*rzone(i))**2*dm1(i)
            strke(i) = (dr(i+1)*r(i+1))**2*dm2(i+1) +
     $       (dh(i)*rzone(i))**2*dm1(i)
            rbar = rbar + r(i+1)*(dr(i+1)*r(i+1))**2*dm2(i+1)
            hbar = hbar + rzone(i)*(dh(i)*rzone(i))**2*dm1(i)
            crsum = crsum + two*dr(i+1)*dm2(i+1)*r(i+1)**2*
     $                      ahf*(dh(i)+dh(i+1))
            wgrav(i+1) = rl1*dh(i)*dm1(i)*gam(i) +
     $       (dr(i+1)*r(i+1))*dm2(i+1)*(gam(i+1)-gam(i))/drz(i+1)
            wcross(i+1) = gor(i+1)*r(i+1)*dr(i+1)*dm2(i+1)*
     $         rl1*(dh(i)+dh(i+1))
            wdiag(i+1) = ( pi4g*ahf*(rho(i)+rho(i+1)) -
     $         four*gor(i+1)/r(i+1) )*(dr(i+1)*r(i+1))**2*dm2(i+1)
           weight(i) = wtherm(i) + wgrav(i+1) + wcross(i+1) + wdiag(i+1)
            stwait(i+1) = stwait(i) + weight(i)
  48     continue
c
c          normalize the weight function.
c
         rke = rsum + rl1*hsum
         if( rke .le. zero ) rke = one
         rbar = (rbar + rl1*hbar)/(rke*r(np1))
         do 46 i=1,n
            stwait(i) = stwait(i)/(omsq*rke)
            weight(i) = weight(i)/(omsq*rke)
            wtherm(i) = wtherm(i)/(omsq*rke)
            wdiag(i) = wdiag(i)/(omsq*rke)
            wgrav(i) = wgrav(i)/(omsq*rke)
            wcross(i) = wcross(i)/(omsq*rke)
            strke(i) = strke(i)/rke
  46     continue
         stwait(np1) = stwait(np1)/rke
         qch = dabs((omsq-stwait(np1))/omsq)
c
c	   Evaluate the rotational spitting coefficients.
c
         crsum = crsum/rke
         hrsum = hsum/rke
         ckl = crsum + hrsum
c
c          Normalize the eigenvalue for direct comparison with Robe,
c       Ann. d'Astrophysique 31 (475) 1968.
c
         qomsq = omsq/(pi*g*rhom)
         index = knode(dr,dh,nmax,n,np,ng)
         stomsq(nroot) = omsq
         st_k_p(nroot) = float(np)
         st_k_g(nroot) = float(ng)
c
c       Assemble the information into a variable and write to a file.
c    The non-adiabatic values are set so the file is always the same.
c
         nonradial%mass_flag = nint(10.d0*totmas/1.991d33)
         nonradial%ell       = lval
         nonradial%k         = index
         nonradial%n_p       = np
         nonradial%n_g       = ng
         nonradial%omega_ad  = sqrt(omsq)
         nonradial%omsq_n    = omsq
         nonradial%nu_ad     = sqrt(omsq)/twopi
         nonradial%period    = twopi/sqrt(omsq)
         nonradial%acc_ad_w  = qch
	 nonradial%rbar      = rbar
         nonradial%del_om_qa = 0.d0
         nonradial%comega   = 0.d0
         nonradial%kappa    = 0.d0
         nonradial%acc_na_k = 0.d0
c
         write(6,5503) qomsq,lval,index,np,ng
         write(6,9999) dr(1),qdr0
         write(6,5507) omsq,stwait(np1),qch,rbar
         write(1,5504) nroot,lval,omsq,stwait(np1),qch
         if( lval .eq. 2 ) call robeig(index,qomsq,fnpol,lval)
         write(1,5503) qomsq,lval,index,np,ng
         write(1,5506) ckl,crsum,hrsum,rbar
         write(1,9999) dr(1),qdr0
         write(11,5500) nroot,lval,omsq,stwait(np1),qch
         write(11,5503) qomsq,lval,index,np,ng
         write(11,5506) ckl,crsum,hrsum,rbar
c
         stwait(np1) = stwait(np1)/omsq
         wdiag(np1) = wdiag(np1)/(omsq*rke)
         if( iout .gt. 0 ) then
            write(11,5502)
            write(11,5501) index,dr(1),eta0,drho0,gam0,dp0
            write(11,5501) (i,dr(i+1),dh(i),adrho(i),gam(i),dp(i),
     $                      weight(i),stwait(i),i=1,n)
         endif
         call orthog(dr,dh,n,omsq,index,nmax,1,lval)
         call pltdmp(dr,nmax,np1,'dr/r')
         call pltdmp(dh,nmax,n,'dh/h')
         call pltdmp(adrho,nmax,n,'drho')
         call pltdmp(dp,nmax,n,'dp  ')
         call pltdmp(gam,nmax,n,'gam ')
         call pltdmp(weight,nmax,n,'wait')
         call pltdmp(strke,nmax,n,'rke ')
         call pltdmp(stwait,nmax,np1,'wint')
         call pltdmp(wtherm,nmax,n,'wthr')
         call pltdmp(wgrav,nmax,np1,'wgrv')
         call pltdmp(wdiag,nmax,np1,'wdia')
         call pltdmp(wcross,nmax,np1,'wcrs')
c
         call summary( nonradial )
c
 100  continue
      call pltdmp(stomsq,200,iroot,'omsq')
      call pltdmp(st_k_p,200,iroot,'k_p ')
      call pltdmp(st_k_g,200,iroot,'k_g ')
      call orthog(dr,dh,n,omsq,index,nmax,-1,lval)
      return
  50  continue
      ihmin = nomega
c
c      Radial, adiabatic pulsation equation solver.
c      Basic references -- Castor, Ap. J. 166 (109) 1971.
c			   Pesnell, PASP, 99 (975), 1987.
c
      drdm0 =-thre/(r(1)*dsqrt(dm2(1)))
      afacp = forpi*r(1)**2/dsqrt(dm2(1))
      afac = forpi*r(2)**2/dsqrt(dm2(2))
      drdm(1,1) = afacp/(v(1)*dm1(1))
      drdm(1,2) =-afac/(v(1)*dm1(1))
      ag1(1,1) = zero
      ag1(1,2) =-four*gor(1)/r(1) + afacp*g1(1)*(p(1)*drdm(1,1)
     $      - drdm0*p0 )
      ag1(1,3) = afacp*g1(1)*p(1)*drdm(1,2)
      do 60 i=2,n
         afacp = afac
         afac = forpi*r(i+1)**2/dsqrt(dm2(i+1))
         drdm(i,1) = afacp/(v(i)*dm1(i))
         drdm(i,2) =-afac/(v(i)*dm1(i))
c
c          set up the mechanical matrices.
c
         ag1(i,1) =-afacp*g1(i-1)*p(i-1)*drdm(i-1,1)
         ag1(i,2) =-four*gor(i)/r(i) + afacp*(g1(i)*p(i)*drdm(i,1)-
     $               g1(i-1)*p(i-1)*drdm(i-1,2) )
         ag1(i,3) = afacp*g1(i)*p(i)*drdm(i,2)
  60  continue
      g1p = g1(n)*p(n)
      afacp = forpi*r(np1)**2/dsqrt(dm2(np1))
      ag1(np1,1) =-afacp*g1p*drdm(n,1)
      ag1(np1,2) =-four*gor(np1)/r(np1)-afacp*g1p*drdm(n,2)
      ag1(np1,3) = zero
c
c          write the matrices to file tape11, if iout is less
c       than 2, this write is not done.
c
      if( iout .ge. 1 ) then
         write(11,7000)
         do 75 i=1,np1
            write(11,7001) i,ag1(i,1),ag1(i,2),ag1(i,3)
  75     continue
      endif
c
c          calculate the acoustic tavel time from surface to the
c       innermost zone. This is stored in transt.
c
      transt = zero
      do 77 i=1,n
         transt = transt+(r(i+1)-r(i))/dsqrt(p(i)*v(i)*g1(i))
  77  continue
      termq = dsqrt(rhom/1.41d0)
      tcon = (pi/transt)**2
      x0np = r(np1)*dsqrt(dm2(np1))
      omsqp = zero
      omsqc = dfloat(ihmin)*tcon
      iv = 0
c
c          start of adiabatic loop
c
      do 79 ih=ihmin,ihmax
      write(11,7900) ih
  78  continue
      omsqc = omsqc*two
      do 80 iterad=1,20
c
      if(iout.ge.1) Write(11,7800) iterad,iv,ih,omsq,omsqc,omsqp
      omsq = ahf*(omsqc+omsqp)
c
c          insure that initial omsqc gives iv .ge. Iv-1
c
      if( iterad .eq. 1) omsq = omsqc
      omsqs = omsq
c
c          iteration on the outer boundary condition
c
      do 81 icount=1,30
c
         do 82 i=1,np1
            xo(i) = zero
            yo(i) = zero
  82     continue
         yo(n) =-ag1(n,3)*x0np
         call trisol(ag1,1,n,omsq,yo,xo,nmax)
         if( icount .eq. 1 ) then
c
c          set up initial derivative
c
            erra = ag1(np1,1)*xo(n)+(ag1(np1,2)-omsq)*x0np
            erra = erra/xo(2)
            omsql = omsq
            omsq  = omsq*(one+1.d-7)
         else
            err = ag1(np1,1)*xo(n)+(ag1(np1,2)-omsq)*x0np
            afacp = dabs((omsq-omsql)/omsql)
            err = err/xo(2)
            if(iout.ge.1) Write(11,4400) icount,omsq,omsql,afacp
            if( afacp .lt. 1.d-10 ) goto 86
            omsq2 = (erra*omsq-err*omsql)/(erra-err)
            omsql = omsq
            domsq = omsq2-omsq
            domsq = dsign(dmin1(dabs(domsq),qrt*dabs(omsq)),domsq)
            omsq  = omsq + domsq
            erra  = err
         endif
  81  continue
c
c          no convergence in the adiabatic eigenvalue.
c
      write(6,8100) omsq,iterad,iv,ih
      omsqp = zero
      omsqc = tcon*dfloat(ih+1)
      goto 79
c
c	   Converged to omsq value, check if positive and with
c      with proper number of nodes in the displacement eigen-
c      vector. If either is not true, try again.
c
  86  continue
      if( omsq .gt. zero ) goto 265
c
c	   Negative omsq value, stop working on this mode and
c       move on.
c
      write(6,8600) omsq,iv,ih
      omsqp = zero
      omsqc = tcon*dfloat(ih+1)
      goto 79
c
 265  continue
      do 90 i=1,n
         yo(i) = zero
  90  continue
      xo(np1) = x0np
      yo(n)  =-ag1(n,3)*x0np
      call trisol(ag1,1,n,omsq,yo,xo,nmax)
      iv = nodes(xo,np1,iout)
      if( iv .eq. Ih-1 ) goto 103
c
c      if initial omsgc gives iv .lt. Ih-1, increase omsgc and try again
c
      if( omsqs .eq. Omsqc .And. Iv .lt. Ih-1 ) goto 78
      if( iv .lt. Ih-1 ) omsqp = omsqs
      if( iv .gt. Ih-1 ) omsqc = omsqs
  80  continue
c
c	   Not converged to a frequency with right number of nodes,
c       try for next mode.
c
      qomsq = omsq/(pi43*g*rhom)
      write(11,8000) iv,ih,qomsq
      omsqp = zero
      omsqc = tcon*dfloat(ih)
      goto 79
c
 103  continue
      xo(n+2) = zero
      stwait(1) = zero
      rke = zero
      do 200 i=1,n
         dr(i) = xo(i)/(r(i)*dsqrt(dm2(i)))
         adrho(i) = drdm(i,1)*xo(i)+drdm(i,2)*xo(i+1)
         dp(i) = g1(i)*adrho(i)
c
c          set up the weight function calculation.
c
         weight(i) =-four*gor(i+1)*xo(i+1)**2/r(i+1) +
     $      p(i)*v(i)*g1(i)*adrho(i)**2*dm1(i)
         stwait(i+1) = stwait(i) + weight(i)
         strke(i) = strke(i) + xo(i+1)*xo(i+1)
         rke = rke + strke(i)
 200  continue
      dr(np1) = one
c
c	    Given the radial eigenvectors, calculate the Epstein
c	weight functions. See Epstein Ap. J. 112(6)1950.
c
      epstwt(1) = (thre*g1(1)-four)*g*rm(2)*dm2(2)/r(2)*dr(1)**2
      chwait = epstwt(1)
      do 104 i=2,n
       epstwt(i)=(thre*g1(i)-four)*g*rm(i+1)*dm2(i+1)/r(i+1)*dr(i+1)**2+
     $     p(i)*v(i)*g1(i)*((dr(i+1)-dr(i))/dlog(r(i+1)/r(i)))**2*dm1(i)
         chwait = chwait + epstwt(i)
 104  continue
      chwait = chwait/rke
      ech = dabs((omsq-chwait)/omsq)
      qomsq = omsq/(pi43*g*rhom)
      zomsq = omsq/tcon
      stwait(np1) = stwait(np1)/rke
      qch = dabs((stwait(np1)-omsq)/omsq)
      write(6,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      write(1,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      write(1,1002) dr(1),zomsq
      write(11,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      iv = nodes(dr,np1,1)
c
c          normalize the weight function per zone to omsq.
c
      strke(1:n) = strke(1:n)/rke
      do 105 i=1,n
         weight(i) = weight(i)/(omsq*rke)
         stwait(i) = stwait(i)/(omsq*rke)
         epstwt(i) = epstwt(i)/(omsq*rke)
 105  continue
         stwait(np1) = stwait(np1)/omsq
c
c       Assemble the information into a variable and write to a file.
c    The non-adiabatic values are set so the file is always the same.
c
         nonradial%mass_flag = nint(10.d0*totmas/1.991d33)
         nonradial%ell       = lval
         nonradial%k         = iv+1
         nonradial%n_p       = iv+1
         nonradial%n_g       = 0
         nonradial%omega_ad  = sqrt(omsq)
         nonradial%omsq_n    = omsq
         nonradial%nu_ad     = sqrt(omsq)/twopi
         nonradial%period    = twopi/sqrt(omsq)
         nonradial%acc_ad_w  = qch
         nonradial%del_om_qa = 0.d0
         nonradial%comega   = 0.d0
         nonradial%kappa    = 0.d0
         nonradial%acc_na_k = 0.d0
c
      if(iout.ge.2) write(11,1000) omsq,qomsq,(i,dr(i+1),adrho(i),
     $    dp(i),weight(i),stwait(i),i=1,n)
c
          call orthog(dr,dh,n,omsq,iv,nmax,1,0)
          call pltdmp(dr,nmax,np1,'dr/r')
          call pltdmp(adrho,nmax,n,'drho')
          call pltdmp(dp,nmax,n,'dp  ')
          call pltdmp(weight,nmax,n,'wait')
          call pltdmp(strke,nmax,n,'rke ')
          call pltdmp(epstwt,nmax,n,'epwt')
          call pltdmp(stwait,nmax,np1,'wint')
c
c	   Finished with this mode, reset the search interval and
c       start on the next one.
c
         omsqp = omsq
         omsqc = omsq
         stomsq(ih) = omsq
c
         call summary( nonradial )
  79  continue
      call pltdmp(stomsq, 200, ihmax-ihmin+1,'omsq')
      call orthog(dr,dh,n,omsq,iv,nmax,-1,0)
      return
c
c          error recovery.
c
 900  continue
      write(1,9000) lval,omlow,omhigh
      write(11,9000) lval,omlow,omhigh
      return
c
 1000 format(1h1,16h radial solution,/,8h  omsq =,1pe14.7,5X,
     $  6hqomsq=,e11.3,//,3X,1hi,5x,4hdr/r,10x,4hdrho,11x,
     $  2hdp,11x,6hweight/(1x,i5,1p,5e14.6))
 1001 format(/,' linear radial mode',/,
     $ 8h  omsq =,1pe13.6,2X,7hqomsq =,e11.3,4H iv=,i3,4h ih=,i3,/,
     $ 2x,'eigenvalue from weight function =',e12.5,' error = ',e10.3,/,
     $ 2X,'eigenvalue from epstein weight function = ',e12.5,
     $ 7H error=,e10.3)
 1002 format(1x,9h dr(1) = ,1pe11.4,26H omsq*(sound t. Time)**2 =,e11.3)
c
 2200 format(1h1,//,21h ag1,drdm,ag3,ag4,ah4)
 2201 format(1x,i5,1p,10e12.4)
 2300 format(1h1,//,20h ah1,ah3,ap1,ap3,ap4)
c
 3001 format(1h1,1x,26hbegin discriminant search.)
 3000 format(1x,i4,24h roots for analysis: l =,i4,/,(3x,1p,5e15.6))
 4400 format(1x,i5,1p,4e15.6)
 4500 format(2x,26hno convergence root number,i5,11h frequency=,1pe15.6)
c
 5500 format(1h1,//,5x,4hroot,i5,3h l=,i4,10x,8h omsq = ,1pe15.8,/,
     $   38H eigenfrequency from weight function =,e12.4,
     $   7H error=,e12.4)
 5501 format(1x,i5,1p,7e15.6)
 5502 format(//4x,1hi,8x,4hdr/r,10x,4hdh/h,10x,8hdrho/rho,8x,5hgamma,
     $   11x,4hdp/p,10x,6hweight)
 5503 format(1x,20h Norm. eigenvalue = ,1pe12.5, ',  ell: ', i4,/,
     $   1X,23h radial quantum number=,i5,25h number of nodes, p-type=,
     $   i5,8h g-type=,i5)
 5504 format(//,5x,4hroot,i5,3h l=,i4,10x,8h omsq = ,1pe15.8,/,
     $   39H eigenfrequency from weight function = ,e12.3,
     $   7H error=,e12.3)
 5506 format(1x,31h rotational splitting, total = ,1pe11.4,5H a = ,
     $   e11.4,5H b = ,e11.4,/,1X,28h effective radius of mode = ,e11.4)
 5507 format(1x,7homsq = ,1pe15.8,13H int. omsq = ,e12.4,8H error =,
     $   e11.4,/,1X,28h effective radius of mode = ,e11.4)
c
 7000 format(1h1,/,5x,26h l=0 pulsation matrix ag1.,/)
 7001 format(1x,i4,1x,1p,3e10.3)
 7800 format(8h iterad=,i3,4h iv=,i3,4h ih=,i3,6h omsq=,1pe12.4,
     $   7H omsqc=,e12.4,7H omsqp=,e12.4)
 7900 format(1h1,//,20x,28h***** begin search for mode ,i2,6h *****)
c
 8000 format(1x,'Cant find no. of nodes = ih-1 ... Nodes=',i3,4h ih=,i3,
     $ 6h omsq=,1pe13.6)
 8600 format(14h adiabatic f2=,1pe17.10,18H is .lt. 0 For iv=,
     $  i3,4h ih=,i3)
 8100 format(14h adiabatic f2=,1pe15.6,20H not converged after,i4,
     $   16h iterations, ih=,i4,17h number of nodes=,i4)
c
 9000 format(43h from lnanon...No adiabatic roots found, l=,i4,/,
     $   17h frequency range:,1p,2e15.6)
 9999 format(1x,7hdr(1) =,1pe11.3,15H scaled dr(1) =,e11.3)
      end
      function knode(x,y,nmax,n,nps,ngs)
      implicit real*8(a-h,o-z)
c
c          returns the number of zero-crossings in the counter-
c       clockwise direction as np and the clockwise crossings
c       in ng.
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax),y(nmax)
      np = 0
      ng = 0
      knode = 0
      if( n .eq. 1 ) Return
      nm1 = n - 1
      do 10 i=2,nm1
         if( x(i)*x(i+1) .lt. zero ) then
            if( y(i)*y(i+1) .lt. zero ) write(11,1000)
            if( x(i+1) .lt. x(i) ) then
               if( y(i+1) .ge. zero ) np = np + 1
               if( y(i+1) .lt. zero ) ng = ng + 1
            else
               if( y(i) .ge. zero ) ng = ng + 1
               if( y(i) .lt. zero ) np = np + 1
            endif
         endif
  10  continue
      if( x(n-1)*x(n) .lt. zero ) write(11,1001)
      knode = np - ng
      nps = np
      ngs = ng
      return
 1000 format(1x,50hfrom knode...Quadrant jump default rotations used.)
 1001 format(1x,51hfrom knode...Node at outer boundary is not counted.)
      end
      function fomsq(omsq,x,xnorm,n)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=4096 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),ap4(nmax)
      common/linear/ atrix(nmax3,7)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax3)
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega**2.
c
      do 200 i=1,n
c
c      ix is index for radial component, idh for horizontal and
c    igam is for the poisson equation.
c
         ix = 3*i - 2
         idh = ix + 2
         igam = ix + 1
c
c          radial component of motion.
c
         x(ix) = zero
         atrix(ix,1) = ag1(i,1)
         atrix(ix,2) = ag3(i,1)
         atrix(ix,3) = ag4(i,1)
         atrix(ix,4) = ag1(i,2) - omsq
         atrix(ix,5) = ag3(i,2)
         atrix(ix,6) = ag4(i,2)
         atrix(ix,7) = ag1(i,3)
c
c          horizontal component of motion.
c
         x(idh) = zero
         atrix(idh,1) = zero
         atrix(idh,2) = ah1(i,1)
         atrix(idh,3) = ah3(i)
         atrix(idh,4) = ah4(i) - omsq
         atrix(idh,5) = ah1(i,2)
         atrix(idh,6) = zero
         atrix(idh,7) = zero
c
c          poisson equation.
c
         x(igam) = zero
         atrix(igam,1) = ap3(i,1)
         atrix(igam,2) = zero
         atrix(igam,3) = ap1(i,1)
         atrix(igam,4) = ap3(i,2)
         atrix(igam,5) = ap4(i)
         atrix(igam,6) = ap1(i,2)
         atrix(igam,7) = ap3(i,3)
 200  continue
      x(3*n) =-ah1(n,2)*xnorm
      x(3*n-1) =-ap1(n,2)*xnorm
      x(3*n-2) =-ag1(n,3)*xnorm
      call rbmles(atrix,nmax3,1,3*n,7,x)
      err = ag1(n+1,1)*x(3*n-2) + (ag1(n+1,2)-omsq)*xnorm +
     $      ag3(n+1,1)*x(3*n-1) + ag4(n+1,1)*x(3*n)
      fomsq = err/x(4)
      return
      end
      subroutine rbmles(a,nmax,imin,imax,m,y)
      implicit real*8(a-h,o-z)
c
c	   Linear equation solver. (R)eal (b)and (m)atrix (l)inear
c       (e)quation (s)olver. Solves a*x = y, destroying a and
c       returns the value of x in y.
c
      dimension a(nmax,m),y(nmax)
      id = imax-imin+1
      idm = id-1
      mmid = (m+1)/2
      mmm = mmid-1
      do 20 ii=1,idm
         I = imax+1-ii
         den = 1.0D0/a(i,mmid)
         y(i) = y(i)*den
         kmax = min0(i-imin,mmm)
         do 10 j=1,mmm
            a(i,j) = a(i,j)*den
            do 10 kk=1,kmax
             jk = j+kk
             k = mmid+kk
             ik = i-kk
             a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  10     continue
         do 20 kk=1,kmax
            ik = i-kk
            k = mmid+kk
            y(ik) = y(ik)-y(i)*a(ik,k)
  20  continue
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
         kmax = min0(i-imin,mmm)
         do 30 kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
30    continue
      return
      end
      subroutine trisol(a,imin,imax,x,y,z,nmaxin)
      implicit real*8(a-h,o-z)
c
      dimension a(nmaxin,3), y(nmaxin), z(nmaxin)
      parameter( nmax=4096 )
      dimension e(nmax),f(nmax)
      id = imax+1-imin
      e(imax) =-a(imax,1)/(a(imax,2)-x)
      f(imax) = y(imax)/(a(imax,2)-x)
      do 10 ii=2,id
         i = imax+1-ii
         den = (a(i,2)-x)+a(i,3)*e(i+1)
         if( ii .ne. id ) then
            e(i) =-a(i,1)/den
         endif
         f(i) = (y(i)-a(i,3)*f(i+1))/den
  10  continue
      imp = imin+1
      z(imin) = f(imin)
      do 20 i=imp,imax
         z(i) = e(i)*z(i-1)+f(i)
  20  continue
      return
      end
      function nodes(w,nz,lout)
      implicit real*8(a-h,o-z)
c
c        Finds no. of nodes(sign changes) in nz elements of array w
c
      dimension w(nz),inode(200)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      nodes = 0
      nodcnt = 0
      nz1 = nz-1
      do 10 i=3,nz1
         if( w(i-1)*w(i) .lt. zero ) then
c
c           zero crossing found. Is it general enough to be a node?
c
            if( w(i-2)*w(i-1) .ge. zero ) then
               if( w(i)*w(i+1) .ge. zero ) then
                  nodcnt = nodcnt + 1
                  inode(nodcnt) = i-1
               endif
            endif
         endif
  10  continue
      inode(nodcnt+1) = nz
      nodes = nodcnt
      if( lout .gt. 0 ) then
         write(11,1000) nodcnt
         if( nodcnt .gt. 0 ) then
            write(11,1001) (inode(i),i=1,nodcnt)
         endif
      endif
 1000 format(/,1x,i2,12h nodes found)
 1001 format(' nodes at ',25i5)
      return
      end
      subroutine orthog(dr,dh,n,omsq,index,nmax1,imode,lval)
      implicit real*8(a-h,o-z)
c
c          check for orthogonality of the wave functions. Part
c      one stores the vertical and horizontal wavefunctions,
c      part two does the integrations.
c
      parameter ( nmax=4096 )
      common/phypar/ r(nmax),T_K(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax), c_s(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/eigstr/ sdr(nmax,4),sdh(nmax,4),stome(4),iomeg(4)
      common/stnmod/ nmode
      dimension dr(nmax1),dh(nmax1),rke(4),cross(4,4)
c
      if( lval .le. 0 ) then
c
c          radial modes.
c
         if( nmode .gt. 4 ) Return
         if( imode .le. 0 ) then
c
c          part two: the integrals.
c
            if( nmode .le. 1 ) Return
            do 160 j1=1,nmode
               do 160 j2=1,nmode
                  cross(j2,j1) = cross(j2,j1) +
     $                              sdr(1,j2)*sdr(1,j1)*dm2(1)*r(1)**2
 160        continue
            do 170 j1=1,nmode
               do 170 j2=1,nmode
                  do 170 i=1,n+1
                     cross(j2,j1) = cross(j2,j1) +
     $                      sdr(i+1,j2)*sdr(i+1,j1)*dm2(i+1)*r(i+1)**2
 170        continue
            do 175 j1=1,nmode
               rke(j1) = dsqrt(cross(j1,j1))
 175        continue
            do 180 j1=1,nmode
               do 180 j2=1,nmode
                  cross(j2,j1) = cross(j2,j1)/(rke(j1)*rke(j2))
 180        continue
            write(1,1800)
            do 185 j1=1,nmode
               do 185 j2=1,nmode
                  write(1,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
 185        continue
         else
c
c          Store the eigenvectors
c
            nmode = nmode + 1
            do 110 i=1,n
               sdr(i,nmode) = dr(i)
 110        continue
            sdr(n+1,nmode) = dr(n+1)
            stome(nmode) = omsq
            iomeg(nmode) = index
         endif
      else
         if( imode .le. 0 ) then
c
c          part two: the integrals.
c
            if( nmode .le. 1 ) Return
            rl1 = dfloat(lval)*dfloat(lval+1)
            do 60 j1=1,nmode
               rke(j1) = zero
                  do 60 j2=1,nmode
                     cross(j2,j1) = sdr(1,j2)*sdr(1,j1)*dm2(1)*r(1)**2
  60        continue
            do 70 j1=1,nmode
               do 70 j2=1,nmode
                  do 70 i=1,n
                     cross(j2,j1) = cross(j2,j1) +
     $                   rl1*sdh(i,j2)*sdh(i,j1)*dm1(i)*rzone(i)**2 +
     $                   sdr(i+1,j2)*sdr(i+1,j1)*dm2(i+1)*r(i+1)**2
  70        continue
            do 75 j1=1,nmode
               rke(j1) = dsqrt(cross(j1,j1))
  75        continue
            do 80 j1=1,nmode
               do 80 j2=1,nmode
                  cross(j2,j1) = cross(j2,j1)/(rke(j1)*rke(j2))
  80        continue
            write(1,8000) lval
            do 85 j1=1,nmode
               do 85 j2=1,nmode
                  write(1,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
  85        continue
         else
            if( nmode .ge. 4 ) Return
            nmode = nmode + 1
            do 10 i=1,n
               sdr(i,nmode) = dr(i)
               sdh(i,nmode) = dh(i)
  10        continue
            sdr(n+1,nmode) = dr(n+1)
            stome(nmode) = omsq
            iomeg(nmode) = index
         endif
      endif
      return
c
 1800 format(1x,34h cross integrals for radial modes:)
 8000 format(1x,4hl = ,i3,27h cross integrals for modes:)
 8001 format( 2(3x,1h<,i3,1h,,i3,1h>,1p,2e12.5) )
      end
      subroutine robeig(index,omsq,fnpol,lin)
      implicit real*8(a-h,o-z)
c
c          comparison of omsq with the eigenvalues given by robe.
c
      dimension omsq1(21),omsq2(21),omsq3(21),!omsq33(21),
     $         omsq35(21),omsq4(21) !,omsq15(21)
      data omsq1/-1.583D-2,-1.898D-2,-2.319D-2,-2.899D-2,-3.731D-2,
     $      -4.989D-2,-7.029D-2,-0.1070D0,-.1844D0,-.4039D0,1.997D0,
     $      12.41D0,31.32D0,57.11D0,89.44D0,128.2D0,173.2D0,224.5D0,
     $      282.2D0,346.d0,416.d0/
      data omsq2/4.061D-2,4.840D-2,5.868D-2,7.265D-2,9.238D-2,.1215D0,
     $      .1672D0,0.2452D0,0.3957D0,0.7510D0,4.151D0,15.41D0,32.10D0,
     $    54.17D0,81.51D0,114.0D0,151.6D0,194.3D0,242.d0,295.d0,353.d0/
      data omsq3/ 0.429D0,0.510D0,0.616D0,0.7588D0,0.9584D0,1.248D0,
     $      1.694D0,2.430D0,3.771D0,6.553D0,10.90D0,20.35D0,35.63D0,
     $  55.29D0,79.23D0,107.4D0,139.8D0,176.5D0,217.d0,262.d0,312.d0/
      data omsq35/1.40D0,1.66D0,2.d0,2.464D0,3.099D0,4.023D0,5.417D0,
     $      7.655D0,11.39D0,16.13D0,21.55D0,27.91D0,39.68D0,57.88D0,
     $      80.39D0,106.8D0,137.2D0,171.d0,210.d0,252.d0,297.d0/
      data omsq4/ 6.2D0,7.34D0,8.82D0,10.78D0,13.44D0,17.02D0,20.48D0,
     $      23.99D0,30.67D0,36.79D0,45.77D0,56.18D0,67.75D0,83.82D0,
     $      102.2D0,116.6D0,140.3D0,171.d0,206.d0,245.d0,288.d0/
c
      lval = lin
      npol = int(10.d0*fnpol)
      if( npol .eq. 0 ) then
c
c          pekeris model
c
         rl = dfloat(lval)
         rl1 = rl*(rl+1.d0)
         k = iabs(index)
         if( index .eq. 0 .and. omsq .gt. 0. ) then
c
c          The f mode
c
            zomsq = 2.d0*rl*(rl-1.d0)/(2.d0*rl+1.d0)
         elseif( index .le. 0 .and. omsq .lt. 0. ) then
c
c          The g mode branch has an altered node count.
c
            k = k + 1
            dn =-2.d0 + dfloat(k)*(rl + 0.5d0 + dfloat(k))*(5.d0/3.d0)
            zomsq = dn - dsqrt(dn**2 + rl1)
         else
c
c          The p mode branch
c
            dn =-2.d0 + dfloat(k)*(rl + 0.5d0 + dfloat(k))*(5.d0/3.d0)
            zomsq = dn + dsqrt(dn**2 + rl1)
         endif
         zomsq = zomsq*(4.d0/3.d0)
         qch = dabs((omsq-zomsq)/omsq)
         write(6,1001) index,omsq,zomsq,qch
         write(1,1001) index,omsq,zomsq,qch
      elseif( lval .eq. 2 ) then
         if( npol .eq. 10 ) then
c
c          n = 1 polytrope, l=2.
c
            call robcmp(index,omsq,omsq1)
         elseif( npol .eq. 20 ) then
c
c          n = 2 polytrope, l=2.
c
            call robcmp(index,omsq,omsq2)
         elseif( npol .eq. 30 ) then
c
c          n = 3 polytrope, l=2.
c
            call robcmp(index,omsq,omsq3)
         elseif( npol .eq. 35 ) then
c
c          n = 3.5 Polytrope, l=2.
c
            call robcmp(index,omsq,omsq35)
         elseif( npol .eq. 40 ) then
c
c          n = 4 polytrope, l=2.
c
            call robcmp(index,omsq,omsq4)
         else
            write(6,1000) fnpol
         endif
      else
         write(6,1010) fnpol,lval
      endif
      return
c
 1000 format(1x,'No comparison n =',f6.2)
 1010 format(1x,'No comparison, n =',f6.2,' l =',i3)
 1001 format(1x,4h k =,i3,8h omsq = ,1pe10.3,
     $   2X,12hreal omsq = ,e10.3,2X,7h error=,e12.3)
      end
      subroutine robcmp(index,omsq,omsqar)
      implicit real*8(a-h,o-z)
      dimension omsqar(21)
      if( iabs(index) .gt. 10 ) Return
      k = index + 11
      zomsq = omsqar(k)
      qch = dabs((omsq-zomsq)/omsq)
      write(6,1001) index,omsq,zomsq,qch
      write(1,1001) index,omsq,zomsq,qch
      return
c
 1001 format(1x,4h k =,i3,8h omsq = ,1pe10.3,
     $   2X,12hreal omsq = ,e10.3,2X,7h error=,e12.3)
      end
