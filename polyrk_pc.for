      program newpol
      implicit real*8(a-h,o-z)
c
c          Polytrope integrator. Integrates the second order equation
c       governing polytropes using a Runge-Kutta method and calculates
c       the linear, adiabatic eigenvalues of the resulting model. The
c       adiabatic index is set to 5/3 in this version.
c
      parameter ( nmax=4096 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/corevl/ rz0,p0,rho0
      common/stnmod/ nmode
      dimension c2(nmax),sl2(nmax)
c
      open (unit=1,file='sumout.pol',status='unknown')
c      open (unit=5,file='sys$input',status='unknown')
c      open (unit=6,file='sys$output',status='unknown')
cVAX      open (unit=11,file='polout.rk',status='new')
cVAX      open (unit=12,file='polplt.rk',status='new')
      open (unit=11,file='polout.rk',status='unknown')
      open (unit=12,file='polplt.rk',status='unknown')
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
      idonpl = 1
   1  continue
      write(6,1000)
      read(5,*,end=910) iomod,npts,fnpol,r0
      if( r0 .le. zero ) r0 = 5.d0
      if( fnpol .lt. zero ) fnpol = thre
      if( npts+2 .gt. nmax ) npts = nmax - 2
c
      call polint(iomod,npts,fnpol,r0,idonpl)
      rhobar = 1.d0/pi43
      rlam = rhobar*(xi(npts)**3/(thre*dthdr(npts)))
      alfa = 1.d0/xi(npts)
      pc = forpi*g*(alfa*rlam)**2/(one+fnpol)
      npuls = npts-2
c
c          Define quantities for NONRK
c
      vn(1) = zero
      un(1) = thre
      do 2 i=2,npts-1
         vn(i) = (one+fnpol)*xi(i)*dthdr(i)/theta(i)
         un(i) = xi(i)**3*theta(i)**fnpol/dthdr(i)
   2  continue
c
c      call tau_dump( npts, xi, theta, rlam, fnpol )
      do 5 i=1,npts
         theta(i) = ahf*(theta(i)+theta(i+1))
   5  continue
c
      asymfr = zero
      do 10 i=1,npuls
         rho(i) = rlam*theta(i+1)**fnpol
         v(i) = one/rho(i)
         p(i) = theta(i+1)*rho(i)*pc/rlam
         g1(i) = 5.d0/3.d0
         c2(i) = p(i)*v(i)*g1(i)
         asymfr = asymfr + (r(i+1)-r(i))/dsqrt(c2(i))
         c2(i) = dlog10(c2(i))
  10  continue
      asymfr = one/asymfr
      bvfac = fnpol - (fnpol+one)/g1(1)
      om02 = pi*g*rhobar
      do 15 i=1,npuls+1
         r(i) = alfa*xi(i+1)
         rm(i) = forpi*alfa**3*rlam*dthdr(i+1)
         gor(i) = g*forpi*alfa*rlam*dthdr(i+1)/xi(i+1)**2
         if( i .eq. npuls + 1 ) goto 15
         if( dabs(bvfac) .gt. 1.d-12 ) then
            rai = bvfac*dthdr(i+1)/(theta(i)*xi(i+1))
            ra(i) = rai
            bv(i) = dlog10((dabs(rai)*gor(i)/r(i))/om02)
         endif
  15  continue
      igmax = 0
      do 16 i=1,npuls
         if( gor(i+1) .lt. gor(i) ) goto 17
          igmax = i+1
  16  continue
  17  continue
      write(6,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
      write(6,1099) asymfr
c      call cjhdmp( npuls, fnpol)
      call pltdmp(bv,nmax,npuls,'bv  ')
      call pltdmp(gor,nmax,npuls+1,'g   ')
      call pltdmp(c2,nmax,npuls,'c2  ')
      rz0 = ahf*r(1)
      rho0 = rlam
      p0 = pc
      do 20 i=1,npuls
         rzone(i) = ahf*(r(i)+r(i+1))
         sl2(i) = dlog10(6.d0/rzone(i)**2) + c2(i)
  20  continue
      call pltdmp(sl2,nmax,npuls,'sl2 ')
      r(npuls+1) = alfa*r0
      tmas = pi43*r(1)**3*rlam
      do 30 i=1,npuls
         dm1(i) = (r(i+1)**3-r(i)**3)*pi43*rho(i)
         tmas = tmas + dm1(i)
  30  continue
      dm2(1) = ahf*(rm(1)+dm1(1))
      do 40 i=2,npuls
         dm2(i) = ahf*(dm1(i)+dm1(i-1))
  40  continue
      dm2(npuls+1) = dm1(npuls)/2.
      call pltdmp(dm2, nmax, npuls, 'dm2 ')
      rmch = (rm(npuls+1)-tmas)/rm(npuls+1)
      write(11,2000) npuls
      write(11,2003) fnpol,r0,alfa,rlam,rhobar,rm(npuls+1),tmas,rmch
      write(11,2006) rz0,p0,rho0
      write(11,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
      write(11,1099) asymfr
      write(11,2004)
      write(11,2005) (i,xi(i),theta(i),dthdr(i),r(i),p(i),rho(i),
     $   gor(i),rm(i),dm1(i),rzone(i),i=1,npts)
c
      qch = rlam/rhobar
      write(1,2007) npuls,fnpol,r0,alfa,rlam,rhobar,rm(npuls+1),
     $      tmas,rmch
      write(1,2008) rz0,p0,rho0,qch
      write(1,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
      write(1,1099) asymfr
c
c          pulsation analysis
c
c          enter <crtl/z> to stop program.
c
c          iopuls := output flag. Iout equal to unity produces
c                   convergence information on file 11, iopuls
c                   equal to 2 prints conv. info., eigenvectors,
c                   and matrices.
c
c          lval := l-value of the pulsation, l < 0 stops program.
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
c          (both omlow and omhigh are entered in units of pi*g*rhobar.)
c
 100  continue
         write(6,1001)
         read(5,*,end=900) iopuls,lval,nomega,ihmax,omlow,omhigh
         if( lval .lt. 0 ) goto 910
         om02 = pi*g*rhobar
         omlow = omlow*om02
         omhigh = omhigh*om02
         nmode = 0
         call lnanon(npuls,iopuls,lval,nomega,ihmax,
     $                                omlow,omhigh,fnpol)
      goto 100
 900  continue
      goto 1
 910  continue
      stop
c
 1000 format(1x,'Enter iout, no. Zones, fnpol, and r0. ',$)
 1001 format(1x,27henter pulsation parameters:,/,
     $  1x,38h io, l, nom, ihmax, omlow, and omhigh.)
 1700 format(1x,15hmax. g in zone ,i4,9h radius =,1pe10.3,
     $      8H mass = ,e10.3)
 1099 format(1x,36hasymptotic frequency for p modes is ,1pe11.4)
 2000 format(41h1 polytropic model in hydrostatic balance,//,
     $      1x,i5,6h zones,/)
 2003 format(7h fnpol=,f6.2,26H guess for surface radius=,1pe10.3,
     $   7H alpha=,e10.3,/,2X,8h lambda=,e10.3,8H rhobar=,e10.3,
     $   8H totmas=,e10.3,11H int. Mass=,e10.3,7H diff.=,E10.3)
 2004 format(3x,1hi,5x,2hxi,7x,5htheta,5x,5hdthdr,4x,6hradius,6x,1hp,
     $   8x,3hrho,8x,1hg,9x,2hrm,7x,3hdm1,7x,5hrzone)
 2005 format(1x,i4,1p10e10.3)
 2006 format(1x,5h rz0=,1pe10.3,4H p0=,e10.3,6H rho0=,e10.3)
 2007 format(1h1,/,1x,18houtput from polyrk,i5,6h zones,/,
     $   7h fnpol=,f6.2,26H guess for surface radius=,1pe10.3,
     $   7H alpha=,e10.3,/,2X,8h lambda=,e10.3,8H rhobar=,e10.3,/,
     $   8H totmas=,e10.3,11H int. Mass=,e10.3,7H diff.=,E10.3)
 2008 format(1x,5h rz0=,1pe10.3,4H p0=,e10.3,6H rho0=,e10.3,/
     $   1X,29h central/average densities = ,e10.3)
      end
      subroutine polint(iout,nin,fnpol,r0,idonpl)
      implicit real*8(a-h,o-z)
      logical lout
c
      parameter ( nmax=4096 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(2),f(2),y(2)
      data accur/1.d-6/
c
      n = nin
      lout = .false.
      if( iout .ge. 1 ) lout = .True.
C
      itry =-1
      itermx = 20
      do 40 iter=1,itermx
c
         dr = r0/dfloat(n+1)
c
         np = 0
         xi(1) = zero
         theta(1) = one
         dthdr(1) = zero
         xi(2) = dr
         theta(2) = one - dr**2*(one-fnpol*dr**2/20.d0)/6.d0
         dthdr(2) = dr**3*(one-fnpol*dr**2/10.d0)/3.d0
         x(1) = dr
         y(1) = theta(2)
         y(2) =-dthdr(2)/dr**2
         do 30 ic=2,n+1
            do 35 m=1,5
               kst = irunge(2,y,f,x(1),dr,m)
               if( kst .eq. 0 ) goto 31
               if( y(1) .lt. zero ) goto 31
                f(1) = y(2)
                f(2) =-two*y(2)/x(1) - y(1)**fnpol
  35        continue
  31       continue
            np = ic + 1
            xi(np) = x(1)
            theta(np) = y(1)
            dthdr(np) =-x(1)*x(1)*y(2)
            if( y(1) .lt. zero ) goto 45
  30     continue
         write(11,4500) fnpol,r0
         if( itry .le. 0 ) goto 49
c
c          Converged to a solution, put the surface at r(n).
c
  45     continue
         if( np .le. 1 ) goto 70
         itry = 1
         dr0 = theta(np)*xi(np)*xi(np)/dthdr(np)
         r0 = xi(np) + dr0
         if( dabs(dr0/r0) .lt. Accur ) goto 50
         write(6,4700) np,fnpol,r0,dr0
         if( lout ) write(11,5501) (i,xi(i),theta(i),dthdr(i),i=1,np)
         goto 40
  49    continue
         r0 = r0*1.25E0
         write(6,4700) np,fnpol,r0,r0
  40  continue
c
c          not converged after itermx tries, write error message and
c       stop the code.
c
      write(11,4000) itermx,fnpol,r0
      write(6,4000) itermx,fnpol,r0
      stop
  50  continue
c
c          converged to a solution, set number of zones and outer
c       boundary value of theta.
c
      theta(np) = zero
      write(11,5500) fnpol,xi(np)
      if( lout ) write(11,5501) (i,xi(i),theta(i),dthdr(i),i=1,np)
      n = np
      nin = np
      call pltint(nmax,n,fnpol,xi(np),idonpl)
      call pltdmp(theta,nmax,n,'thet')
      call pltdmp(dthdr,nmax,n,'dtdr')
      return
  70  continue
c
c          Unknown problem in the integration. Run away, run away...
c
      write(6,7000) iter,np,x(1),y(1),y(2),f(1),f(2)
      stop
c
 4000 format(1x,20hno convergence after,i4,18h tries for fnpol =,f7.3,
     $      9H radius =,1pe11.4)
 4500 format(1x,22hno convergence fnpol =,f7.3,9H radius =,1pe11.4)
 4700 format(2x,3hnp=,i4,6h fnpol,f7.3,8H radius=,1pe11.3,5H dr0=,e11.3)
 5500 format(/,5x,5hfnpol,f7.3,8H radius=,1pe15.6)
 5501 format(1x,i5,1p3e15.6)
 7000 format(1x,2i3,1p5e11.3)
c
      end
      function irunge(n,y,f,x,h,m)
      implicit real*8(a-h,o-z)
      dimension phi(20),savey(20),y(n),f(n)
c
      if( m .eq. 1 ) then
c
c          step 1.
c
         irunge = 1
      elseif( m .eq. 2 ) then
c
c          step 2.
c
         do 25 j=1,n
            savey(j) = y(j)
            phi(j) = f(j)
            y(j) = y(j) + h*f(j)/2.d0
  25     continue
         x = x + h/2.d0
         irunge = 1
      elseif( m .eq. 3 ) then
c
c         step 3.
c
         do 35 j=1,n
            phi(j) = phi(j) + 2.d0*f(j)
            y(j) = savey(j) + f(j)*h/2.d0
  35     continue
         irunge = 1
      elseif( m .eq. 4 ) then
c
c          step 4.
c
         do 45 j=1,n
            phi(j) = phi(j) + 2.d0*f(j)
            y(j) = savey(j) + f(j)*h
  45     continue
         x = x + h/2.d0
         irunge = 1
      elseif( m .eq. 5 ) then
c
c         final pass.
c
         do 55 j=1,n
            y(j) = savey(j) + (phi(j) + f(j))*h/6.d0
  55     continue
         irunge = 0
      endif
      return
      end
      block data
      implicit real*8(a-h,o-z)
c
c       Fundamental constants are from Novotny, Introduction to
c    Stellar Atmospheres and Interiors, 1973, Appendix II.
c                                                 2/16/83  WD Pesnell
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
c      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726D-8,7.5595D-5,
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/1.d0,7.5595D-5,
     $  5.66961D-5,3.1415926536D0,6.2831853072D0,12.5663706144D0,
     $   25.1327412288D0,4.1887902048D0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,10.0D0,
     $      0.5D0,0.25D0      /
      end
      subroutine pltint(nmaxin,n,fnpol,r0,idonpl)
      implicit real*8(a-h,o-z)
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=4096 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      if( idonpl .gt. 0 ) then
         idonpl =-1
         write(12,1000) fnpol,r0
      endif
      call pltdmp(xi(2),nmaxin,n,'xi  ')
      do  20 i=2,n
         r(i) = xi(i)/xi(n)
  20  continue
      call pltdmp(r(2),nmaxin,n-1,'x   ')
      return
c 1000 format(1x,18hpolytrope with n =,f5.2,20H^h2}x^lh.5%0 ^Lxhx%=,f7.4)
 1000 format(1x,'Polytrope with n =',f5.2,' outer radius= ',f7.4)
      end
c      subroutine tau_dump( npts, xi, theta, rlam, fnpol )
c      implicit real*8(a-h,o-z)
c      character*80 filename
c      dimension xi(npts), theta(npts)
cc
c      write(filename, 1000) nint(10.*fnpol)
c      open (31, file=filename, status='unknown' )
c      write(31,1105) npts
c      do 5 i=1,npts
c         write(31,1100) xi(i)/xi(npts), rlam*theta(i)**fnpol
c   5  continue
c      close(31)
c      return
cc
c 1000 format('poly_',i2.2,'.txt')
c 1105 format(i6)
c 1100 format(1x,1p2g12.5)
c      end

      subroutine lnanon(nin,iout,lin,nomega,ihmax,omlow,omhigh,fnpoln)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=4096 )
      parameter ( nmax3=3*nmax )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk37/  drdm(nmax,2)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/corevl/ rz0,p0,rho0
      common/scrtch/ ag1(nmax,3),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),ap4(nmax)
      common/linear/ wtherm(nmax),wgrav(nmax),wcross(nmax),wdiag(nmax),
     $      dr(nmax),dh(nmax),gam(nmax),strke(nmax),z(nmax),
     $      dp(nmax),adrho(nmax),weight(nmax),xo(nmax),
     $      yo(nmax),stwait(nmax),epstwt(nmax),spac(nmax,5)
      dimension drz(nmax),rkl2(nmax),drint(nmax)
      dimension x(nmax3),stomeg(200), stomsq(200,4)
      data accur/1.d-12/
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
      drz(1) = drz(2)
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      drz(np1) = p(n)*v(n)/gor(np1)
      drz(np1) = drz(n)
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
c
c          Second form
c
c      drdm01 =-(rl+one)/rz0**2*(23.d0/3.d0)
c      drdm02 = (rl+one)/rz0**2*(5.d0/3.d0)
c
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
c         if( i .eq. 1 ) then
c            ah1(i,1) = drdm01*g1(i)*p(i)*v(i) + ahf*gzone
c            ah1(i,2) = drdm02*g1(i)*p(i)*v(i) + ahf*gzone
c         else
            ah1(i,1) = drdm(i,1)*g1(i)*p(i)*v(i) + ahf*gzone
            ah1(i,2) = drdm(i,2)*g1(i)*p(i)*v(i) + ahf*gzone
c         endif
         ah3(i) = one
         ah4(i) = rkl2(i)*g1(i)*p(i)*v(i)
         if( i .eq. n ) goto 20
c
c          radial momentum equation.
c
         afac = forpi*r(i+1)**3/dm2(i+1)
         ag1(i+1,1) = drdm(i,1)*g1(i)*p(i)*(ahf*rl*v(i) - afac)
         ag1(i+1,2) =-for*gor(i+1)/r(i+1)+pi4g*ahf*(rho(i+1)+rho(i))+
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
c      vl =-g1(n)*gor(np1)*r(np1)
c      vl =-g1(n)*p(n)*forpi*r(np1)**3/dm2(np1)
c      ag1(np1,1) = vl*drdm(n,1)
c      fdr = (one + ahf*rl*drz(np1)/r(np1))/(one + dlnr)
c      ag1(np1,1) = drdm(i,1)*g1(i)*p(i)*(ahf*rl*v(i) - afac)
c      ag1(np1,2) =-for*gor(np1)/r(np1) + vl*drdm(n,2) -
c     $    pi4g*rho(n)*ahf*fdr
c
c          02/04/2000
c
      afac = forpi*r(np1)**3/dm2(np1)
c      vl = 2.d0*g1(n)*p(n)*(ahf*rl*v(n) - afac)
      vl = -two*g1(n)*p(n)*afac
      ag1(np1,1) = drdm(n,1)*vl
c      fdr = -ahf*(one + rl)*(drz(np1)/r(np1))/(one + dlnr)
      fdr = (one + ahf*rl*(drz(np1)/r(np1)))/(one + dlnr)
      ag1(np1,2) =-for*gor(np1)/r(np1) + drdm(n,2)*vl +
     $    pi4g*rho(n)*ahf*fdr
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
      dlrhdr =-rho(n)*gor(n+1)/p(n)*fnpol/(one+fnpol)
      ap1(n,2) =-pi4g*rho(n)*( (ddr13*two*rl*rzone(n) + ap3n3)*f +
     $    rzone(n)**2*drdm(n,2) - ahf*dlrhdr )
c
      ap3(n,1) = ap3n1 + ddr11*two*rl*rzone(n)
      fp3 = (rl+ahf)*(drz(np1)/r(np1)) - 1
      ap3(n,2) =-(ap3n1+ap3n3) + ddr12*two*rl*rzone(n) -
     $    (ddr13*two*rl*rzone(n) + ap3n3)*fp3/(one+dlnr)
      ap3(n,3) = zero
c
      ap4(n) =-pi4g*rl1*rho(n)
      ah4(np1) = zero
c
c        write the matrices to file 6.
c
      if( iout .ge. 2 ) then
         write(11,2200) lval
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
c          Convergence to a normal mode is governed by the momentum
c       equation for the outermost zone. As in the radial case, this
c       function (divided by any x value) is used in a secant method
c       for controlling the iterations. As the equations are not
c       nonlinear (as in non-non) in the frequency, it is hoped the
c       calculation is, in some sense, more stable.
c
c                                               10/4/83 WD Pesnell
c
      write(11,3001)
      xnorm = r(np1)**2
      iroot = 0
      dom = (omhigh - omlow)/dfloat(nomega - 1)
      do 30 i=1,nomega
         omsq = omlow + dfloat(i-1)*dom
         fomeg = fomsq(omsq,x,xnorm,n)
         if( i .eq. 1 ) then
            fomi = fomeg
         else
            if( fomeg*fomi .lt. zero ) then
               iroot = iroot + 1
               stomeg(iroot) = omsq - dom*fomeg/(fomeg-fomi)
            endif
            fomi = fomeg
         endif
         if( iout .gt. 0 ) write(11,4400) i,omsq,fomeg
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
      uc = 0.1d0
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
               if( afac .le. accur ) goto 47
               omsq2 = (erra*omsq-err*omsq1)/(erra-err)
               omsq1 = omsq
               domsq = omsq2-omsq
               domsq = dsign(dmin1(dabs(domsq),uc*dabs(omsq)),domsq)
               omsq  = omsq + domsq
               erra  = err
            endif
            if(iout.gt.0) write(11,4400) icount,omsq,afac,err,erra
  45     continue
         write(6,4500) nroot,omsq
         write(1,4500) nroot,omsq
         write(11,4500) nroot,omsq
         goto 100
c
c          converged to a root, calculate the eigenvectors.
c
  47     continue
         drtst = (ag1(n+1,1)*x(3*n-2)+ag3(n+1,1)*x(3*n-1)+
     $      ag4(n+1,1)*x(3*n))/(omsq-ag1(n+1,2))
         write(6,*) drtst,xnorm,dabs((drtst-xnorm)/xnorm)
c
         x(3*n+1) = drtst
         do 40 i=1,n
            dr(i) = (r(i)/r(np1))**lval*(x(3*i-2)/r(i)**2)
            dh(i) = (x(3*i)/rzone(i)**2)*(rzone(i)/r(np1))**lval
            adrho(i) = (drdm(i,1)*x(3*i-2)+drdm(i,2)*x(3*i+1) +
     $               rkl2(i)*x(3*i) )*(rzone(i)/r(np1))**lval
            dp(i) = g1(i)*adrho(i)
            gam(i) = x(3*i-1)*(rzone(i)/r(np1))**lval
  40     continue
         eta0 = dr(1)/rl
         gam0 = x(2)*(rz0/r(np1))**lval
         drho0 = adrho(1)*(rz0/rzone(1))**lval
         dp0 = drho0*g1(1)
c         dr(n+1) = one
         dr(n+1) = x(3*n+1)/r(np1)**2
c
c          Variables now refer to physical mesh and need different factors
c
         dlnr = (rl+one)*ahf*drz(np1)/r(np1)
         fdr = (drz(np1)/r(np1))/(one + dlnr)
         gam(n+1) = gam(n)*(one-dlnr)/(one+dlnr)
     $             - pi4g*(rho(n)/2.d0)*dr(np1)*fdr
         dh(n+1) = dh(n)
         qdr0 = dr(2)*(r(1)/r(2))**(lval-2)
c
c          weight function for the adiabatic oscillation.
c
c         wgrav(1) = dr(1)*dm2(1)*rl*gam(1)
         wgrav(1) = zero
         wcross(1) = two*rl1*gor(1)*r(1)*dr(1)*dh(1)*dm2(1)
         wdiag(1) = (pi4g*rho(1)-gor(1)*for/r(1) )*
     $               (dr(1)*r(1))**2*dm2(1)
         stwait(1) =  wgrav(1) + wcross(1) + wdiag(1)
c
         rsum = (dr(1)*r(1))**2*dm2(1)
         hsum = zero
         strke(1) = (dr(1)*r(1))**2*dm2(1)
         rbar = r(1)*(dr(1)*r(1))**2*dm2(1)
         hbar = zero
         crsum = zero 
         istart = 0
         do 48 i=1,n
            if( i .lt. istart ) then
               wtherm(i) = zero
               strke(i+1) = zero
               wgrav(i) = zero
               wcross(i+1) = zero
               wdiag(i+1) = zero
            else
            rsum = rsum + (dr(i+1)*r(i+1))**2*dm2(i+1)
            hsum = hsum + (dh(i)*rzone(i))**2*dm1(i)
            strke(i+1) = (dr(i+1)*r(i+1))**2*dm2(i+1) +
     $       (dh(i)*rzone(i))**2*dm1(i)
            rbar = rbar + r(i+1)*(dr(i+1)*r(i+1))**2*dm2(i+1)
            hbar = hbar + rzone(i)*(dh(i)*rzone(i))**2*dm1(i)
            crsum = crsum + two*dr(i+1)*dm2(i+1)*r(i+1)**2*
     $                      ahf*(dh(i)+dh(i+1))
c
c          Calculate the various weight function terms
c
            wtherm(i) = g1(i)*p(i)*v(i)*dm1(i)*adrho(i)**2
c            wgrav(i+1) = rl1*dh(i)*dm1(i)*gam(i) +
c     $       (dr(i+1)*r(i+1))*dm2(i+1)*(gam(i+1)-gam(i))/drz(i+1)
            drho_Eul = adrho(i)+fnpol/(one+fnpol)*
     $               (rho(i)/P(i))*ahf*(gor(i)*dr(i) + gor(i+1)*dr(i+1))
            wgrav(i) = drho_Eul*gam(i)*dm1(i) 
            if( i .eq. n ) then
               wcross(i+1) = gor(i+1)*r(i+1)*dr(i+1)*dm2(i+1)*
     $            2.d0*rl1*dh(i)
               wdiag(i+1) = ( pi4g*ahf*(rho(i)+rho(i+1))
     $         - for*gor(i+1)/r(i+1) )*(dr(i+1)*r(i+1))**2*dm2(i+1)
            else
               wcross(i+1) = gor(i+1)*r(i+1)*dr(i+1)*dm2(i+1)*
     $            rl1*(dh(i)+dh(i+1))
               wdiag(i+1) = ( pi4g*ahf*(rho(i)+rho(i+1)) -
     $            for*gor(i+1)/r(i+1) )*(dr(i+1)*r(i+1))**2*dm2(i+1)
            endif
           endif
            weight(i) = wtherm(i) + wgrav(i) + wcross(i+1) + wdiag(i+1)
            stwait(i+1) = stwait(i) + weight(i)
  48     continue
c
c          normalize the weight function.
c
         rke = rsum + rl1*hsum
         if( rke .le. zero ) rke = one
         rbar = (rbar + rl1*hbar)/(rke*r(np1))
         omsq_Weight = stwait(np1)/rke
         qch = dabs((omsq-omsq_Weight)/omsq)
c
c          Evaluate the rotational spitting coefficients.
c
         crsum = crsum/rke
         hrsum = hsum/rke
         ckl = crsum + hrsum
c
c          Normalize the eigenvalue for direct compaison with Robe,
c       Ann. d'Astrophysique 31 (475) 1968.
c
         stomsq(nroot,1) = omsq
         stomsq(nroot,2) = omsq_Weight
c
         qomsq = omsq/(pi*g*rhom)
         index = knode(dr,dh,nmax,n,np,ng)
         write(6,5503) qomsq,index,np,ng
         write(6,9999) lval, x(1), x(3)
         write(6,5507) omsq,omsq_Weight,qch,rbar
         write(1,5504) nroot,lval,omsq,omsq_Weight,qch
         call robeig(index,qomsq,fnpol,lval, omsq_comp)
         stomsq(nroot,4) = omsq_comp
         write(1,5503) qomsq,index,np,ng
         write(1,5506) ckl,crsum,hrsum,rbar
         write(1,9999) lval, x(1), x(3)
         write(11,5500) nroot,lval,omsq,omsq_Weight,qch
         write(11,5503) qomsq,index,np,ng
         write(11,5506) ckl,crsum,hrsum,rbar
c
         do i=1,np1
            stwait(i) = stwait(i)/abs(omsq*rke)
            weight(i) = weight(i)/abs(omsq*rke)
            wtherm(i) = wtherm(i)/abs(omsq*rke)
            wdiag(i) = wdiag(i)/abs(omsq*rke)
            wgrav(i) = wgrav(i)/abs(omsq*rke)
            wcross(i) = wcross(i)/abs(omsq*rke)
            strke(i) = strke(i)/rke
         enddo
c
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
         call pltdmp(strke,nmax,np1,'rke ')
         call pltdmp(stwait,nmax,np1,'wint')
         call pltdmp(wtherm,nmax,n,'wthr')
         call pltdmp(wgrav,nmax,np1,'wgrv')
         call pltdmp(wdiag,nmax,np1,'wdia')
         call pltdmp(wcross,nmax,np1,'wcrs')
         if( fnpol .eq. 0. ) then
            k = iabs(index)
            call pekeris( k, lval, r, wcross, n)
          endif
c         call wndnum(dr,dh,nmax,n)
c
 100  continue
      call pltdmp( stomsq(1,1), 200, iroot, 'omsq' )
      call pltdmp( stomsq(1,2), 200, iroot, 'wpes' )
      call pltdmp( stomsq(1,3), 200, iroot, 'wsch' )
      do i=1, iroot
         write(6,5600) (stomsq(i,j), j=1,4)
      enddo
 5600 format(1x,4(' &  ',f11.5))
      call orthog(dr,dh,n,omsq,index,nmax,-1,lval)
      return
  50  continue
      ihmin = nomega
c
c      Radial, adiabatic pulsation equation solver.
c      Basic references -- Castor, Ap. J. 166 (109) 1971.
c                          Pesnell, PASP, 99 (975), 1987.
c
      drdm0 =-thre/(r(1)*dsqrt(dm2(1)))
      afacp = forpi*r(1)**2/dsqrt(dm2(1))
      afac = forpi*r(2)**2/dsqrt(dm2(2))
      drdm(1,1) = afacp/(v(1)*dm1(1))
      drdm(1,2) =-afac/(v(1)*dm1(1))
      ag1(1,1) = zero
      ag1(1,2) =-for*gor(1)/r(1) + afacp*g1(1)*(p(1)*drdm(1,1)
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
         ag1(i,2) =-for*gor(i)/r(i) + afacp*(g1(i)*p(i)*drdm(i,1)-
     $               g1(i-1)*p(i-1)*drdm(i-1,2) )
         ag1(i,3) = afacp*g1(i)*p(i)*drdm(i,2)
  60  continue
      g1p = g1(n)*p(n)
      afacp = 2.d0*forpi*r(np1)**2/dsqrt(dm2(np1))
      ag1(np1,1) =-afacp*g1p*drdm(n,1)
      ag1(np1,2) =-for*gor(np1)/r(np1)-afacp*g1p*drdm(n,2)
      ag1(np1,3) = zero
c
c          write the matrices to file tape11, if iout is less
c       than 2, this write is not done.
c
      if( iout .ge. 1 ) then
         write(11,7000)
         do 75 i=1,np1
            write(11,7001) i,ag1(i,1),ag1(i,2),ag1(i,3),
     $                       drdm(i,1),drdm(i,2)
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
      if(iout.ge.1) write(11,7800) iterad,iv,ih,omsq,omsqc,omsqp
      omsq = ahf*(omsqc+omsqp)
c
c          insure that initial omsqc gives iv .ge. iv-1
c
      if( iterad .eq. 1) Omsq = omsqc
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
            if(iout.ge.1) write(11,4400) icount,omsq,omsql,afacp
            if( afacp .lt. 1.d-12 ) goto 86
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
c          Converged to omsq value, check if positive and with
c      with proper number of nodes in the displacement eigen-
c      vector. If either is not true, try again.
c
  86  continue
      if( omsq .gt. zero ) goto 265
c
c          Negative omsq value, stop working on this mode and
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
      if( iv .eq. ih-1 ) goto 103
c
c      if initial omsgc gives iv .lt. ih-1, increase omsgc and try again
c
      if( omsqs .eq. omsqc .and. iv .lt. ih-1 ) goto 78
      if( iv .lt. ih-1 ) omsqp = omsqs
      if( iv .gt. ih-1 ) omsqc = omsqs
  80  continue
c
c          Not converged to a frequency with right number of nodes,
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
      strke(1) = zero
      stwait(1) = zero
      rke = zero
      do 200 i=1,n
         dr(i) = xo(i)/(r(i)*dsqrt(dm2(i)))
         adrho(i) = drdm(i,1)*xo(i)+drdm(i,2)*xo(i+1)
         dp(i) = g1(i)*adrho(i)
c
c          set up the weight function calculation.
c
         weight(i) =-for*gor(i+1)*xo(i+1)**2/r(i+1) +
     $      p(i)*v(i)*g1(i)*adrho(i)**2*dm1(i)
         stwait(i+1) = stwait(i) + weight(i)
         strke(i+1) = xo(i+1)*xo(i+1)
         rke = rke + xo(i+1)*xo(i+1)
 200  continue
      dr(np1) = one
c
c           Given the radial eigenvectors, calculate the Epstein
c       weight functions. See Epstein Ap. J. 112(6)1950.
c
      epstwt(1) = (thre*g1(1)-for)*g*rm(2)*dm2(2)/r(2)*dr(1)**2
      chwait = epstwt(1)
      do 104 i=2,n
        epstwt(i)=(thre*g1(i)-for)*g*rm(i+1)*dm2(i+1)/r(i+1)*dr(i+1)**2+
     $     p(i)*v(i)*g1(i)*((dr(i+1)-dr(i))/dlog(r(i+1)/r(i)))**2*dm1(i)
         chwait = chwait + epstwt(i)
 104  continue
      chwait = chwait/rke
      ech = dabs((omsq-chwait)/omsq)
      qomsq = omsq/(pi43*g*rhom)
      zomsq = omsq/tcon
      stwait(np1) = stwait(np1)/rke
c
      stomsq(ih,1) = omsq
      stomsq(ih,2) = stwait(np1)
      stomsq(ih,3) = chwait
c
      qch = dabs((stwait(np1)-omsq)/omsq)
      write(6,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      write(1,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      write(1,1002) dr(1),zomsq
      write(11,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      iv = nodes(dr,np1,1)
      call robeig(ih, qomsq,fnpol,lval, omsq_comp)
c
c          normalize the weight function per zone to omsq.
c
      strke(1:np1) = strke(1:np1)/rke
      do 105 i=1,n
         weight(i) = weight(i)/abs(omsq*rke)
         stwait(i) = stwait(i)/abs(omsq*rke)
         epstwt(i) = epstwt(i)/abs(omsq*rke)
 105  continue
         stwait(np1) = stwait(np1)/omsq
c
      if(iout.ge.2) write(11,1000) omsq,qomsq,(i,dr(i+1),adrho(i),
     $    dp(i),weight(i),stwait(i),i=1,n)
c
          call orthog(dr,dh,n,omsq,iv,nmax,1,0)
          call pltdmp(dr,nmax,np1,'dr/r')
          call pltdmp(adrho,nmax,n,'drho')
          call pltdmp(dp,nmax,n,'dp  ')
          call pltdmp(weight,nmax,n,'wait')
          call pltdmp(epstwt,nmax,n,'epwt')
          call pltdmp(strke, nmax, np1, 'rke ')
          call pltdmp(stwait,nmax,np1,'wint')
c
c          Finished with this mode, reset the search interval and
c       start on the next one.
c
         omsqp = omsq
         omsqc = omsq
  79  continue
      iroot = ihmax-ihmin+1
      call pltdmp( stomsq(1,1), 200, iroot, 'omsq' )
      call pltdmp( stomsq(1,2), 200, iroot, 'wpes' )
      call pltdmp( stomsq(1,3), 200, iroot, 'wsch' )
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
     $  2hdp,11x,6hweight/(1x,i4,1p5e14.6))
 1001 format(/,19h linear radial mode,/,
     $ 8h  omsq =,1pe13.6,2X,7hqomsq =,e11.3,4H iv=,i3,4h ih=,i3,/,
     $ 2x,33heigenvalue from weight function =,e12.5,8H error =,e10.3,/,
     $ 1X,42h eigenvalue from epstein weight function =,e12.5,
     $ 7H error=,e10.3)
 1002 format(1x,9h dr(1) = ,1pe11.4,26H omsq*(sound t. Time)**2 =,e11.3)
c
 2200 format(1h1,1x,'Matrices for l=', i4//,1x,'ag1,drdm,ag3,ag4,ah4')
 2201 format(1x,i5,1p10e12.4)
 2300 format(1h1,//,20h ah1,ah3,ap1,ap3,ap4)
c
 3001 format(1h1,1x,26hbegin discriminant search.)
 3000 format(1x,i4,24h roots for analysis: l =,i4,/,(3x,1p5e15.6))
 4400 format(1x,i5,1p4e15.6)
 4500 format(2x,'No convergence root number ',i5,' frequency=',1pe15.6)
c
 5500 format(1h1,//,5x,'root',i5,' l =',i4,10x,' omsq = ',1pe15.8,/,
     $   ' eigenfrequency from weight function =',e12.4,
     $   ' error =',e12.4)
 5501 format(1x,i5,1p7e15.6)
 5502 format(//4x,1hi,8x,4hdr/r,10x,4hdh/h,10x,8hdrho/rho,8x,5hgamma,
     $   11x,4hdp/p,10x,6hweight)
 5503 format(/1x,'Norm. eigenvalue = ',1pe12.5,/,
     $   1x,' radial quantum number =',i5,' number of nodes, p-type =',
     $   i5,' g-type =',i5)
 5504 format(//,5x,'root',i5,' l=',i4,10x,' omsq = ',1pe15.8,/,
     $   ' eigenfrequency from weight function = ',e12.3,
     $   ' error =',e12.3)
 5506 format(1x,31h rotational splitting, total = ,1pe11.4,5H a = ,
     $   e11.4,5H b = ,e11.4,/,1X,28h effective radius of mode = ,e11.4)
 5507 format(1x,7homsq = ,1pe15.8,13H int. omsq = ,e12.4,8H error =,
     $   e11.4,/,1X,28h effective radius of mode = ,e11.4)
c
 7000 format(1h1,/,5x,'l=0 pulsation matrix ag1 and drdm:',/)
 7001 format(1x,i4,1x,1p3e11.3,2x,3e11.3)
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
     $   17h frequency range:,1p2e15.6)
 9999 format(1x,'l = ',i3,' u(1) = ',1pe11.3,' y(1) = ',g10.3)
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
c      dimension al(nmax3,7), indx(nmax3)
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
         igam = ix + 1
         idh = ix + 2
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
 200  continue
      x(3*n-2) =-ag1(n,3)*xnorm
      x(3*n-1) =-ap1(n,2)*xnorm
      x(3*n)   =-ah1(n,2)*xnorm
      call rbmles(atrix,nmax3,1,3*n,7,x)
c
c      m1 = 3
c      m2 = 3
c      mp = 7
c      mpl = 7
c      call bandec(atrix,3*n,m1,m2,nmax3,mp,al,mpl,indx,d)
c	call banbks(atrix,3*n,m1,m2,nmax3,mp,al,mpl,indx, x)
      err = ag1(n+1,1)*x(3*n-2) + (ag1(n+1,2)-omsq)*xnorm +
     $      ag3(n+1,1)*x(3*n-1) + ag4(n+1,1)*x(3*n)
      fomsq = err/x(4)
      return
      end
      subroutine rbmles(a,nmax,imin,imax,m,y)
      implicit real*8(a-h,o-z)
c
c          Linear equation solver. (R)eal (b)and (m)atrix (l)inear
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
c         do 30 kk=1,kmax
         do kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
         enddo
30    continue
      return
      end
      subroutine trisol(a,imin,imax,x,y,z,nmaxin)
      implicit real*8(a-h,o-z)
c
      parameter( nmax=4096 )
      dimension a(nmaxin,3), y(nmaxin), z(nmaxin)
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
c           zero crossing found. is it general enough to be a node?
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
 1001 format(7h nodes=,30i4)
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
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $      gor(nmax),xi(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $      un(nmax), drhdr(nmax), ra(nmax)
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
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
            if( nmode .le. 1 ) return
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
 1800 format(1x,'Cross integrals for radial modes:')
 8000 format(1x,'l = ',i3,' cross integrals for modes:')
 8001 format( 2(3x,'<',i3,',',i3,'>',1p2e12.5) )
      end
      subroutine robeig(index,omsq,fnpol,lin, omsq_comp)
      implicit real*8(a-h,o-z)
c
c          comparison of omsq with the eigenvalues given by robe.
c
      dimension omsq1(21),omsq2(21),omsq3(21), !omsq33(21),
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
      dimension n_HRW(7)
      data n_HRW/100, 150, 200, 300, 325, 350, 400/
      dimension om2_HRW(3,7)
      data om2_HRW/ 1.892, 12.09, 27.08, ! n=1.
     $              2.706, 12.54, 26.58, ! n=1.5
     $              4.001, 13.34, 26.58, ! n=2.
     $              9.255, 16.98, 28.48, ! n=3.
     $              11.03, 18.89, 29.87, ! n=3.25
     $              12.64, 21.21, 32.08, ! n=3.50
     $              15.15, 24.94, 37.07/ ! n=4.
c
      lval = lin
      npol = ifix(sngl(10.e0*fnpol))
      if( npol .eq. 0 ) then
         if( lval .eq. 0 ) then
c
c          pekeris model
c
            k = iabs(index)
c
c          The p mode branch
c
            dn =-2.d0 + dfloat(k)*(0.5d0 + dfloat(k))*(5.d0/3.d0)
            zomsq = dn + dsqrt(dn**2 + rl1)
         else
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
               dn =-2.d0 + dfloat(k)*(rl + 0.5d0 +dfloat(k))*(5.d0/3.d0)
               zomsq = dn - dsqrt(dn**2 + rl1)
            else
c
c          The p mode branch
c
               dn =-2.d0 + dfloat(k)*(rl + 0.5d0 +dfloat(k))*(5.d0/3.d0)
               zomsq = dn + dsqrt(dn**2 + rl1)
            endif
            omsq_comp = zomsq
            zomsq = zomsq*(4.d0/3.d0)
            qch = dabs((omsq-zomsq)/omsq)
            write(6,1001) index,omsq,zomsq,qch
            write(1,1001) index,omsq,zomsq,qch
         endif
      elseif( lval .eq. 0 ) then
         npol = ifix(sngl(100.e0*fnpol))
         k = index
         if( k .le. 3 ) then
            do i=1, 7
               if( n_HRW(i) .eq. npol ) then
                  zomsq = om2_HRW(k,i)
                  goto 10
               endif
            enddo
         endif
         return
  10     continue
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
      if( omsq .le. 0. ) index = index-1
      if( iabs(index) .gt. 10 ) return
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
c      subroutine wndnum(x,y,nmax,nptsin)
c      implicit real*8(a-h,o-z)
c      logical iwflag
c      save iwflag
c
c      common/winstr/ iwflag
c      dimension x(nmax),y(nmax)
c      data iwflag/.true./
c
c      if( iwflag ) then
c         open (unit=30,file='polyrk_eig',status='new')
c         iwflag = .false.
c      endif
c      if( iwflag .le. 0 ) Iwflag = 1
c      write(30,1000) nptsin
c      write(30,1001) (x(i+1),i=1,nptsin)
c      write(30,1001) (y(i),i=1,nptsin)
c      return
c 1000 format(i5)
c 1001 format(1p6e12.5)
c      end
