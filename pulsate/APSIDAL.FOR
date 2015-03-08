      program apsidal
      implicit real*8(a-h,o-z)
      character*80 ititl
      integer*2 i, j, k, imin, imax
c
c          Calculate the apsidal constant for the model.
c
c                                                 05/05/95 WD Pesnell
c
      parameter ( nmax=512 )
      common/phypar/ rp(nmax),tp(nmax),vp(nmax),cv(nmax),dkdr(nmax),
     $   dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),ru(nmax),U(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,pi4,pi8,pi43
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/aps/    np1
c
c          read in the input variables from file 10.
c
      open (unit=1,file='sumout',status='unknown',access='append')
      open (unit=11,file='nonout',status='new')
c
c          read in the input variables from file 10.
c
      open (unit=10,file='lnrad',status='old',readonly)
      read(10,1001,end=100) ititl
      read(10,*) npts,rlumgv,totmas,teff,rphoto,corlum
      read(10,*) irad,noburn,onemq0,onemq1,rlums,rlumc
      np1 = npts+1
      read(10,1000) (rp(i),i=1,np1)
      read(10,1000) (tp(i),i=1,npts)
      read(10,1000) (vp(i),i=1,npts)
      read(10,1000) (cv(i),i=1,npts)
      read(10,1000) (dkdr(i),i=1,npts)
      read(10,1000) (dkdt(i),i=1,npts)
      read(10,1000) (dm1(i),i=1,npts)
      read(10,1000) (akap(i),i=1,npts)
      read(10,1000) (dm2(i),i=1,np1)
      read(10,1000) (rm(i),i=1,np1)
      read(10,1000) (p(i),i=1,npts)
      read(10,1000) (g1(i),i=1,npts)
      read(10,1000) (g3m1(i),i=1,npts)
      read(10,1000) (frft(i),i=1,np1)
      read(10,1000) (sorce(i),i=1,npts)
      read(10,1000) (dtsorc(i),i=1,npts)
      read(10,1000) (dvsorc(i),i=1,npts)
      read(10,1000) (bv(i),i=1,np1)
      read(10,1000,end=900) pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,
     $   g10,g3m10,cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      close (unit=10)
      totmas = rm(npts+1)
c
c          output the initial model data
c
c      write(11,2000) ititl,npts
c      write(11,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
c      write(11,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
c      write(11,2004)
c      write(11,2005) i,rp(1),tc,1.e0/rhoc,pc,g10,g3m10,opac0,
c     $      dkdr0,dkdt0,cv0,frft(1),rm(1)
c      write(11,2005) (i,rp(i+1),tp(i),vp(i),p(i),g1(i),g3m1(i),akap(i),
c     $      dkdr(i),dkdt(i),cv(i),frft(i+1),rm(i+1),i=1,npts)
c
      call pltint(npts,ititl)
      write(1,2000) ititl,npts
      write(1,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(1,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
c
      write(6,1003)
      ru(1) = 0.d0
      U(1) = 3.d0
      do i=1, npts
         ru(i+1) = rp(i+1)/rp(npts+1)
         rho_i = 2.d0/(vp(i) + vp(i+1))
         U(i+1) = pi4*rp(i+1)**3*rho_i/rm(i+1)
      enddo
      call apsidal_k( npts+1, rk_2 )
      write(6,3000) rk_2, dlog10(rk_2)
      stop
c
c          if end of file, try again.
c
 100  continue
      write(1,1002)
      stop
c
 900  continue
      stop
c
 1000 format(1p,4e20.13)
 1001 format(a)
 1002 format(1x,'End-of-File on title, try again.')
 1003 format(1x,27hEnter pulsation parameters.,/,
     $   1x,' io, l, nom, ihmax, permin, and permx: ',$)
c
 2000 format(50h1 initial model in hydrostatic and thermal balance,//,
     $      1x,a,1x,i5,6h zones,/)
 2001 format(4x,18h total luminosity=,1pe11.4,12h total mass=,e11.4,
     $ 17h effective temp.=,e11.4,/,6x,21h photospheric radius=,e11.4,
     $ 13h core radius=,e11.4,11h core lum.=,e11.4)
 2002 format(6h irad=,i2,8h noburn=,i2,8h onemq0=,1pe10.3,
     $ 8h onemq1=,e10.3,7h rlums=,e10.3,7h rlumc=,e10.3)
 2003 format(6h iout=,i2,13h l-value used,i4,
     $ 25h number of points in scan,i4,15h starting freq.,1pe10.3,
     $ 13h ending freq.,e10.3)
 2004 format(3x,1hi,3x,6hradius,5x,4htemp,5x,6hsp.vol,3x,8hpressure,
     $ 5x,2hg1,7x,4hg3m1,5x,7hopacity,3x,6hdlkdlr,4x,6hdlkdlt,6x,
     $ 2hcv,7x,4hfrft,4x,8hint.mass)
 2005 format(1x,i4,1p,12e10.3)
 2100 format(/,30h model has no nuclear burning.)
 2101 format(/,29h following zones are ignited.)
 2102 format(3(1x,i4,1p,3e10.3))
 3000 format(1x,'Apsidal constant is: ',1pe12.4,3x,0pf8.4)
      end
      block data
      implicit real*8(a-h,o-z)
c
c       Fundamental constants are from Novotny, Introduction to
c    Stellar Atmospheres and Interiors, 1973, Appendix II.
c                                                 2/16/83  WDP
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/thermo/ r,a,bk,avagd,ad3
      data r,a,bk,avagd,ad3 / 8.31434d7,7.56471d-15,8.6170837d-5,
     $     6.02217d23,2.52157d-15 /
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726d-8,7.5595d-5,
     $  5.66961d-5,3.1415926536d0,6.2831853072d0,12.5663706144d0,
     $   25.1327412288d0,4.1887902048d0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,10.0d0,
     $      0.5d0,0.25d0      /
      end
      subroutine pltint(n,ititl)
      implicit real*8(a-h,o-z)
      character*80 ititl
      integer*2 i
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=512 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	 dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),omq(nmax),tlog(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      dimension plog(nmax),xrad(nmax)
c
      open (unit=12,file='apsidal.plt',status='new')
      write(12,1000) ititl
      np1 = n + 1
      call pltdmp(r,nmax,n+1,'r   ')
      sfact = one - 5.d-14
      do 10 i=1,n
         omq(i) =-dlog10(one - sfact*rm(i)/totmas)
         tlog(i) = dlog10(t(i))
         plog(i) = dlog10(p(i))
         xrad(i) = r(i)/r(np1)
  10  continue
      xrad(np1) = one
      omq(n+1) =-dlog10(one-sfact)
      call pltdmp(tlog,nmax,n,'t   ')
      call pltdmp(plog,nmax,n,'p   ')
      call pltdmp(omq,nmax,n+1,'1-q ')
      call pltdmp(xrad,nmax,n+1,'x   ')
      call pltdmp(g1,nmax,n,'g1   ')
      call pltdmp(g3m1,nmax,n,'g3m1')
      do 20 i=1,n
         xrad(i) = dfloat(i)
  20  continue
      call pltdmp(xrad,nmax,n,'zone')
      return
 1000 format(a)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
      character*4 ititl
      integer*2 i
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p,6e12.4) )
      end

      subroutine apsidal_k( npts, retreg )
      implicit real*8(a-h,o-z)
      integer*2 i, j
      logical convrg
c
c          Integration of the external function F at the using the
c       adaptive Runge-Kutta routine in Numerical Recipes. This
c       integral is from 0 to +\infty, and is stopped by convergence 
c       of the value.
c
c                                                 01/27/92 W. D. Pesnell
c
      parameter ( pi=3.141592653589793d0, zero=0.d0 )
      parameter ( one=1.d0, pi2=pi*2.d0 )
      parameter ( tiny=1.d-20, eps=1.d-6 )
      parameter ( nvarmx=6)
      dimension ystart(nvarmx), yscal(nvarmx), y(nvarmx)
      dimension reg(nvarmx), dydx(nvarmx)
      parameter (nmax=512)
      common/blk4/   y2(nmax),g1(nmax),g3m1(nmax),ru(nmax),U(nmax)
      dimension fi(nmax)
c
      nint = 1000
      nvar = 2
      do 5 i=1, nint
         fi(i) = zero
   5  continue
      do 10 i=1, nvarmx
         y(i) = zero
         reg(i) = zero
  10  continue
c
c           Set up the spline coefficients.
c
      call spline(ru, U, npts, yp1, ypn, y2)
      y(1) = 1.d0
      y(2) = 0.d0
      x = 1.d-2*ru(2)
      H = x/10.d0
      convrg = .false.
      do 90 i=2,nint
c
c          Find the F function for current value of X
c
         call Fc( x, y, dydx )
         do 80 j=1,nvar
            yscal(j) = abs(y(j)) + abs(H*dydx(j)) + tiny
  80     continue
         call rkqc( y, dydx, nvar, x, H, eps, yscal, Hdid, Hnext, Fc)
         fi(i) = y(1)
         H = Hnext
         Nverg = 0
         do 85 j=1, nvar
            write(6,8500) x,H,y(j),dydx(j)
            reg(j) = y(j)
  85     continue
         if( convrg ) goto 20
         if( x + H .gt. ru(npts) ) then
            convrg = .true.
            H = ru(npts) - x
         endif
  90  continue
c
  20  continue
      retreg = (3.d0*y(1) - ru(npts)*y(2))/
     $         (4.d0*y(1) + 2.d0*ru(npts)*y(2))
      return
 8500 format(1x,'X: ',1pe11.3,' H: ',e11.3,' Y: ',e11.3,' dydx: ',e11.3)
c
      end
      subroutine Fc( Xin, Y, dcde )
      implicit real*8(a-h,o-z)
c
c          Definition of apsidal motion constant. See Schwarzschild,
c       Structure and Evolution of Stars, p. 153.
c
      parameter (nmax=512)
      common/blk4/   y2(nmax),g1(nmax),g3m1(nmax),ru(nmax),U(nmax)
      common/aps/    npts
      dimension Y(1), dcde(1)
c
      call intspl(Xin, U_x, dydx, ru, U, y2, npts, istart)
      write(6,2000) Xin, U_x
      dcde(1) = y(2)
      dcde(2) = 2.d0*(3.d0 - U_x)*y(1)/Xin**2 - 2.d0*U_x*y(2)/Xin
c
      return
 2000 format(1x,1pe12.4,3x,4e12.4)
      end

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      implicit real*8(a-h,o-z)
      integer*2 i, nmax
c
c          Set up the coefficients for a cubic spline. Values
c       of yp1 and ypn are the first derivatives at the first 
c       and last mesh points respectively. 
c
      parameter (nmax=1024)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
c
      IF (YP1 .GT. .99d30) THEN
	Y2(1) = 0.d0
	U(1) = 0.d0
      ELSE
	Y2(1) =-0.5d0
	U(1) = (3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
	SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
	P = SIG*Y2(I-1)+2.d0
	Y2(I) = (SIG-1.d0)/P
	U(I) = (6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF( YPN .GT. .99d30 ) THEN
	QN = 0.d0
	UN = 0.d0
      ELSE
	QN = 0.5d0
	UN = (3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      end
      subroutine intspl(xin, yout, dydx, xa, ya, y2a, n, istart)
      implicit real*8(a-h,o-z)
      integer*2 klow, khigh, k
c
c          Uses the arrays from spline to interpolate a value for
c       for y at x.
c
      dimension xa(n), ya(n), y2a(n)
      save klow, khigh
      x = xin
      if( istart .gt. 0 ) then
         klow = 1
         khigh = n 
      else
         klow = 1
         khigh = n 
      endif
      if( x .gt. xa(n) ) then
         y = ya(n)
      elseif( x .lt. xa(1) ) then
         y = ya(1)
      else
  10     continue
            if( khigh-klow .gt. 1 ) then
               k = (khigh+klow)/2
               if( xa(k) .gt. x ) then
                  khigh = k
               else
                  klow = k
               endif
               goto 10
            endif
c
         xfound = xa(khigh) - xa(klow)
         a = (xa(khigh)-x)/xfound
         b = (x-xa(klow))/xfound
         y = a*ya(klow) + b*ya(khigh) +
     $	       ( a*(a*a - 1.d0)*y2a(klow) + b*(b*b-1.d0)*y2a(khigh) )*
     $             (xfound**2)/6.d0
         dydx = (ya(khigh) - ya(klow) )/xfound +
     $	       ( (3.d0*a*a-1.d0)*y2a(klow)+(3.d0*b*b-1.d0)*y2a(khigh) )*
     $             (xfound)/6.
      endif
c      write(6,1000) x, klow, khigh, y
c      write(6,1010) x, klow, khigh, y,
c     $		       ya(klow), ya(khigh), y2a(klow), y2a(khigh),
c     $		       a, b
      yout = y
      return
c 1000 format(1x,'x=',1pe10.3,' indices are:',2i4,' Yout=',e10.3)
c 1010 format(1x,'x=',1pe10.3,' indices are:',2i4,' Yout=',e10.3,/,
c     $	 3x,'y is bracketed by:',2e11.3,/,
c     $	 3x,'y2 is bracketed by:',2e11.3, 3x, 2e11.3 )
      end
