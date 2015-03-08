      program polywd
      implicit real*8(a-h,o-z)
      character*11 ig1ttl(2)
c
c	   White dwarf polytrope integrator. the various constants
c       in the totally degenerate eos have been taken from the
c       routine degen, not chandrasekhar's book. the value of the
c       mean molecular weight/electron is assumed to be unity. the
c       second part performs an adiabatic pulsation analysis on the
c	model. For l=0 the castor code is used and for l>0 the Lagrangian
c	code developed at JILA dominates.
c
c                                                    wd pesnell 5/19/86
c
      parameter ( aelec=6.003d22 , rhod=9.8104d5 )
      parameter ( onthd=0.333333333333333d0, qmsol=1.991d33 )
      parameter ( nmax=600 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),eta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/corevl/ rz0,p0,rho0,g10
      common/stnmod/ nmode
      dimension glog(nmax),sl2(nmax)
      data ig1ttl(1)/'consistent.'/,ig1ttl(2)/'5/3.       '/
c
c          files in this program.
c  unit   routine    name                purpose
c    1      main    sumout.dat          summary output file
c    5      main    sys$input           terminal input file
c    6      main    sys$output          terminal output file
c   11      main    polout.dat          main output file
c
c   12     pltint   polplt.dat          has vectors of model
c                                          quantities for pltcmd program
c
c   41     cjhdmp   cjh_in.dat          has vectors of model
c                                          quantities for c. j. hansen's
c                                          pulsation code.
c
cVAX	  open (unit=1,file='sumout',status='unknown',access='append')
      open (unit=1,file='sumout',status='unknown')
cVAX	  open (unit=5,file='sys$input',status='old')
cVAX	  open (unit=6,file='sys$output',status='old')
      open (unit=11,file='polout',status='unknown')
c
c          iout:= output flag. iout greater than unity produces
c                 voluminous lines of output.
c
c          ir:= number of points in model. ir must be less than
c               nmax-1, which is the dimension of the matrices
c               used.
c
c          y02:= parameter describing the central condensation
c                of the model. the values of y02 range from 0.01
c                to 0.95 or so. default value is 0.01 (used if the
c                inputted value is .le. 0).
c
c          r0:= initial guess at the radius of the model. this
c               parameter is found as an eigenvalue in the solution
c               and a good guess speeds the convergence. default
c               value is 5 (used if the inputted value is .le. 0).
c
c          rmuein:= mean molecular weight per electron. if .le. 0
c                   the default is unity. note that harper and rose,
c                   (aP. j. 162 (963) 1970) use rmue = 2. the omsq
c                   (omega squared) variable scales as rmue, thus the
c                   period scales as 1/sqrt(rmue).
c
c          ig1:= greater than zero for g1=5/3, .le. 0 g1 given by the
c                usual completely degenerate formula.
c
c          nmax is the maximum dimension of many of the matrices.
c
c
c          read in the input variables from file 5.
c
   1  continue
      write(6,1000)
      read(5,*,end=900) iout,ir,y02,r0,rmuein,ig1
      if( r0 .le. zero ) r0 = 5.d0
      if( y02 .le. zero ) y02 = 0.01d0
      if( ir .ge. nmax-1 ) ir = nmax-1
c
      call polint(iout,ir,y02,r0)
      q = (one-y02)**1.5d0
      rmue = rmuein
      if( rmuein .le. zero ) rmue = one
      rhomue = rhod*rmue
      rlam = rhomue*q/y02**1.5
      rhobar = rlam*thre*dthdr(ir)/(q*eta(ir)**3)
      alfa = 7.71381d8*dsqrt(y02)/rmue
      npuls = ir-2
c
c          calculate the dimensionless stellar structure functions.
c
      vfac = g*forpi*rlam*alfa**2/aelec/q
      vn(1) = zero
      un(1) = thre
      drhdr(1) = zero
      do  2 i=2,ir-1
         rhoi = rlam*(theta(i)**2 - y02)**1.5d0/q
         x3 = rhoi/rhomue
         x = x3**onthd
         x2 = x*x
         sqx2p1 = dsqrt(one+x2)
         fx = x*(two*x2-thre)*sqx2p1+thre*dlog(x+sqx2p1)
         vn(i) = vfac*dthdr(i)*rhoi/(eta(i)*fx)
         un(i) = eta(i)**3*(theta(i)**2-y02)**1.5d0/dthdr(i)
         drhdr(i) =-thre*theta(i)*dthdr(i)/(eta(i)*(theta(i)**2-y02))
   2  continue
      call pltdmp(vn,nmax,ir-1,'vn  ')
      call pltdmp(un,nmax,ir-1,'un  ')
      do 5 i=1,ir
         theta(i) = ahf*(theta(i)+theta(i+1))
   5  continue
      fivths = 5.d0/thre
      asymp = zero
      r(npuls+1) = alfa*eta(npuls+2)
      do 10 i=1,npuls
         rho(i) = rlam*(theta(i+1)**2 - y02)**1.5d0/q
         v(i) = one/rho(i)
         x3 = rho(i)/rhomue
         x = x3**onthd
         x2 = x*x
         sqx2p1 = dsqrt(one+x2)
         fx = x*(two*x2-thre)*sqx2p1+thre*dlog(x+sqx2p1)
         p(i) = aelec*fx
         g1(i) = 8.d0*x2*x3/(thre*sqx2p1*fx)
         r(i) = alfa*eta(i+1)
         sl2(i) = dlog10(6.d0*g1(i)*p(i)*v(i)/r(i)**2)
         rm(i) = forpi*alfa**3*rlam*dthdr(i+1)/q
         gor(i) = g*forpi*alfa*rlam*dthdr(i+1)/(q*eta(i+1)**2)
         glog(i) = dlog10(gor(i))
         ra(i) = zero
	 if( ig1 .gt. 0 ) then
	    g1(i) = fivths
	    sl2(i) = dlog10(6.d0*g1(i)*p(i)*v(i)/r(i)**2)
	    ra(i) = dthdr(i+1)*(thre/x3-8.d0*x2/(g1(i)*sqx2p1*fx))
	    bv(i) = ra(i)*gor(i)/r(i)
	    dlnr = (rzone(i) - rzone(i-1))/r(i)
	    if( i .gt. 1 )
     $	       asymp = asymp + dlnr/dsqrt(dabs(bv(i)))
	    bv(i) = dlog10(dabs(bv(i)))
	 endif
  10  continue
      write(6,*) asymp
      rm(npuls+1) = forpi*alfa**3*rlam*dthdr(npuls+2)/q
      gor(npuls+1)=g*forpi*alfa*rlam*dthdr(npuls+2)/(q*eta(npuls+2)**2)
      glog(npuls+1) = dlog10(gor(npuls+1))
      if( ig1 .gt. 0 ) bv(npuls+1) = dlog10(dabs(gor(npuls+1)/r(npuls+1)
     $     *dthdr(npuls+2)*(thre/x3-8.d0*x2/(g1(npuls)*sqx2p1*fx)) ) )
      igmax = 0
      do 16 i=1,npuls
         if( gor(i+1) .lt. gor(i) ) goto 17
          igmax = i
  16  continue
  17  continue
      write(6,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
      call cjhdmp(npuls,y02)
      if( ig1 .gt. 0 ) call pltdmp(bv,nmax,npuls,'bv  ')
      call pltdmp(sl2,nmax,npuls,'sl2 ')
      call pltdmp(glog,nmax,npuls,'glog')
      rz0 = ahf*r(1)
      rho0 = rlam
      x3 = rho0/rhomue
      x = x3**onthd
      x2 = x*x
      sqx2p1 = dsqrt(one+x2)
      fx = x*(two*x2-thre)*sqx2p1+thre*dlog(x+sqx2p1)
      p0 = aelec*fx
      g10 = 8.d0*x2*x3/(thre*fx*sqx2p1)
      if( ig1 .gt. 0 ) g10 = fivths
      do 20 i=1,npuls
         rzone(i) = ahf*(r(i)+r(i+1))
  20  continue
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
      dm2(npuls+1) = ahf*dm1(npuls)
      rmch = (rm(npuls+1)-tmas)/rm(npuls+1)
      tmas = tmas/qmsol
      qmas = rm(npuls+1)/qmsol
      write(6,*) qmas
      write(11,2000) npuls
      write(11,2003) y02,r0,alfa,rlam,rhobar,rmue,qmas,tmas,rmch
      write(11,2006) rz0,p0,rho0
      write(11,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
      write(11,2004)
      write(11,2005) (i,eta(i),theta(i),dthdr(i),r(i),p(i),rho(i),
     $   g1(i),gor(i),rm(i),rzone(i),i=1,ir)
c
c          initialize the summary file.
c
      qch = rlam/rhobar
      write(1,2007) npuls,y02,r0,alfa,rlam,rhobar,rmue,qmas,tmas,rmch
      write(1,2008) rz0,p0,rho0,qch
      if( ig1 .gt. 0 ) write(1,2009) ig1ttl(2)
      if( ig1 .le. 0 ) write(1,2009) ig1ttl(1)
      write(1,1700) igmax,r(igmax)/r(npuls+1),rm(igmax)/rm(npuls+1)
c
c          start the pulsation calculation.
c
c          type <crtl/z> to end program.
c
c          iopuls:= output flag. iout equal to unity produces
c                   convergence information on file 11, iopuls
c                   equal to 2 prints conv. info., eigenvectors,
c                   and matrices.
c
c          lval:= l-value of the pulsation.
c
c          nomega:= number of points in scan calculation (l>0),
c                   lowest order mode in the radial case.
c
c          ihmax:= highest order mode calculated in the l=0 case,
c                  must be greater than ihmin.
c
c          perlow:= shortest period used in scan calculation (l>0),
c                   not used in radial case.
c
c          perhi:= longest period used in scan calculation (l>0),
c                  not used in radial case.
c          (it is not necessary to order perlow and perhi.)
c
 100  continue
	 write(6,1001)
	 read(5,*,end=950) iopuls,lval,nomega,ihmax,perlow,perhi
	 if( lval .gt. 0 ) then
	    omlow = (twopi/perhi)**2
	    omhigh = (twopi/perlow)**2
	 endif
	 nmode = 0
	 call lnanon(npuls,iopuls,lval,nomega,ihmax,omlow,omhigh)
      goto 100
 950  continue
      close (unit=12)
      goto 1
 900  continue
      stop
c
 1000 format(1x,46hEnter iout, no. zones, yo2, r0, rmue, and ig1.)
 1001 format(1x,27hEnter pulsation parameters:,/,
     $	1x,38h io, l, nom, ihmax, perlow, and perhi.)
 1700 format(1x,15hmax. g in zone ,i4,9h radius =,1pe10.3,
     $       7h mass =,e10.3)
 2000 format(41h1 polytropic model in hydrostatic balance,//,
     $	    1x,i5,6h zones,/)
 2003 format(5h yo2=,1pe12.5,26h guess for surface radius=,e10.3,
     $   7h alpha=,e10.3,/,2x,8h lambda=,e10.3,8h rhobar=,e10.3,
     $   7h rmue =,e10.3,/,8h totmas=,e10.3,11h int. mass=,e10.3,
     $   15h (solar masses),7h diff.=,e10.3)
 2004 format(3x,1hi,3x,6hradius,5x,5htheta,5x,5hdthdr,7x,1hr,9x,1hp,
     $   8x,3hrho,7x,2hg1,9x,1hg,8x,2hmr,7x,5hrzone)
 2005 format(1x,i4,1p10e10.3)
 2006 format(1x,5h rz0=,1pe10.3,4h p0=,e10.3,6h rho0=,e10.3)
 2007 format(1h1,/,1x,18houtput from polywd,i5,11h zone model,/,
     $	 1x,4hy02=,1pe12.5,26h guess for surface radius=,e10.3,
     $   7h alpha=,e10.3,/,2x,8h lambda=,e10.3,8h rhobar=,e10.3,
     $   7h rmue =,e10.3,/,8h totmas=,e10.3,11h int. mass=,e10.3,
     $   15h (solar masses),7h diff.=,e10.3)
 2008 format(1x,5h rz0=,1pe10.3,4h p0=,e10.3,6h rho0=,e10.3,/,
     $	 1x,28h central/average densities =,e11.4)
 2009 format(1x,41hthe value of the adiabatic index (g1) is ,a11)
      end
      subroutine polint(iout,nin,y02,r0)
      implicit real*8(a-h,o-z)
      logical lout
c
c          integrate the polytropic equation with index 1.5 and
c       the additional y0 term. the results represent white dwarf
c       stars with varying degrees of central condensation and
c       radius. for equations see chandrasekhar "an intro. to
c       stellar structure", (1967 dover) ch. xii, p. 417. the
c       solution method is a fourth-order runge-kutta integration
c       (carnahan, luther, and wilkes) with the outer radius of
c       model treated as an eigenvalue.
c          written december 16-21, 1985 at jila, boulder, co.
c
c                                            12/21/85 wd pesnell
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),eta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension f(2),y(2)
      data accur/1.d-6/
c
      n = nin
      lout = .false.
      if( iout .ge. 1 ) lout = .true.
      itry =-1
      q = dsqrt(one - y02)
      sqy02 = dsqrt(y02)
      do 40 iter=1,20
         dr = r0/dfloat(n+1)
         eta(1) = zero
         theta(1) = one
         dthdr(1) = zero
         eta(2) = dr
         theta(2) = one - dr**2*q**3*(one-q*thre*dr**2/20.d0)/6.d0
         dthdr(2) = dr**3*q**3*(1.d0/thre-q*dr**2/10.d0)
	 x = dr
         y(1) = theta(2)
         y(2) =-dthdr(2)/dr**2
         do 30 ic=2,n+1
            do 35 m=1,5
	       kst = irunge(2,y,f,x,dr,m)
               if( kst .eq. 0 ) goto 31
               if( y(1) .lt. sqy02 ) goto 31
               f(1) = y(2)
	       f(2) =-two*y(2)/x - (y(1)**2 - y02)**1.5
  35        continue
  31        continue
            np = ic + 1
	    eta(np) = x
            theta(np) = y(1)
	    dthdr(np) =-y(2)*x*x
            if( y(1) .lt. sqy02 ) goto 45
  30     continue
         write(11,4500) y02,r0
         if( itry .le. 0 ) goto 49
c
c          converged to a solution, put the surface
c       at eta(n).
c
  45     continue
         if( np .le. 2 ) stop
         itry = 1
         dr0 =-(sqy02-theta(np))*eta(np)*eta(np)/dthdr(np)
         r0 = eta(np) + dr0
         if( dabs(dr0/r0) .lt. accur ) goto 50
         write(6,4700) y02,r0,dr0
         if( lout ) write(11,5501) (i,eta(i),theta(i),dthdr(i),i=1,np)
         goto 40
  49     continue
         r0 = r0*1.05e0
         write(6,4700) y02,r0,r0
  40  continue
      write(11,4500) y02,r0
      return
  50  continue
      theta(np) = sqy02
      write(11,5500) y02,eta(np)
      n = np
      nin = np
      call pltint(np,y02,r0)
      call pltdmp(theta,nmax,np,'thet')
      call pltdmp(dthdr,nmax,np,'dtdr')
      return
  60  continue
      write(11,4500) y02,r0
      write(6,4500) y02,r0
      stop
c
 4400 format(1x,i5,1p4e15.6)
 4500 format(20h no convergence y02=,f7.2,8h radius=,1p2e15.6)
 4700 format(5x,3hy02,f7.2,8h radius=,1pe11.3,5h dr0=,e11.3)
 5500 format(/,5x,3hy02,f7.2,8h radius=,1pe15.8)
 5501 format(1x,i5,1p3e15.6)
c
      end
      function irunge(n,y,f,x,h,m)
      implicit real*8(a-h,o-z)
      dimension phi(10),savey(10),y(n),f(n)
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
  25	 continue
	 x = x + h/2.d0
	 irunge = 1
      elseif( m .eq. 3 ) then
c
c         step 3.
c
	 do 35 j=1,n
	    phi(j) = phi(j) + 2.d0*f(j)
	    y(j) = savey(j) + f(j)*h/2.d0
  35	 continue
	 irunge = 1
      elseif( m .eq. 4 ) then
c
c          step 4.
c
	 do 45 j=1,n
	    phi(j) = phi(j) + 2.d0*f(j)
	    y(j) = savey(j) + f(j)*h
  45	 continue
	 x = x + h/2.d0
	 irunge = 1
      elseif( m .eq. 5 ) then
c
c         final pass.
c
	 do 55 j=1,n
	    y(j) = savey(j) + (phi(j) + f(j))*h/6.d0
  55	 continue
	 m = 0
	 irunge = 0
      endif
      return
      end
      block data
      implicit real*8(a-h,o-z)
c
c       fundamental constants are from novotny, introduction to
c    stellar atmospheres and interiors, 1973, appendix ii.
c                                                 2/16/83  wdp
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
c
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726d-8,7.5595d-5,
     $  5.66961d-5,3.1415926536d0,6.2831853072d0,12.5663706144d0,
     $   25.1327412288d0,4.1887902048d0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,10.0d0,
     $      0.5d0,0.25d0      /
      end
      subroutine pltint(n,y02,r0)
      implicit real*8(a-h,o-z)
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),eta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
c
      open (unit=12,file='polplt',status='unknown')
c
      write(12,1000) y02,r0
      call pltdmp(eta,nmax,n,'eta ')
      do 10 i=1,n
         rzone(i) = eta(i)/eta(n)
  10  continue
      call pltdmp(rzone,nmax,n,'x   ')
      return
c 1000 format(1x,32hwhite dwarf with Y^lh.5\0^lxhx\=,f5.3,
c     $      17h {C^lh.5\0^lxhx\=,f6.4)
c
 1000 format(1x,'White dwarf: y02 = ',f5.3,
     $		' Outer radius = ',f6.4)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
      character*4 ititl
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p6e12.4) )
      end
      subroutine lnanon(nin,iout,lin,nomega,ihmax,omlow,omhigh)
      implicit real*8(a-h,o-z)
c
      parameter ( accur=1.d-12 )
      parameter ( nmax=600 )
      parameter ( nmax3=3*nmax )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),ceta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/blk37/  drdm(nmax,2)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/corevl/ rz0,p0,rho0,g10
      common/scrtch/ ag1(nmax,3),ag3(nmax,2),ag4(nmax,2),
     $		     ah1(nmax,2),ah3(nmax),  ah4(nmax),
     $		     ap1(nmax,2),ap3(nmax,3),ap4(nmax)
      common/linear/ wtherm(nmax),wgrav(nmax),wdiag(nmax),
     $	    dr(nmax),dh(nmax),gam(nmax),y(nmax),
     $	    z(nmax),dp(nmax),adrho(nmax),weight(nmax),
     $	    xo(nmax),yo(nmax),stwait(nmax),
     $	    spac(nmax,7)
      dimension drz(nmax),rkl2(nmax),drint(nmax)
      dimension x(nmax3),stomeg(200)
c      nmax = (nmax)
c
c          zero the eigenvector matrix to avoid indefinite operands.
c
      do 5 i=1,nmax3
         x(i) = zero
   5  continue
c
c          define the zone-centered radii and the differences
c       of the two radii (zone and interface centered) for ease
c       of future work. the array gor is the local gravity
c       divided by the radius squared and rho is the zone density.
c
      n = nin
      np1 = n + 1
      lval = lin
      rl = dfloat(lval)
      rl1 = rl*(rl+one)
      rhom = rm(np1)/(pi43*r(np1)**3)
      do 10 i=1,n
         rkl2(i) = rl1/rzone(i)**2
         if(i.gt.1) drz(i) = rzone(i)-rzone(i-1)
         drint(i) = r(i+1)-r(i)
         drdm(i,1) = forpi*rzone(i)/(dm1(i)*v(i)) -
     $       ahf*(rl+one)/rzone(i)**2
         drdm(i,2) =-(forpi*rzone(i)/(dm1(i)*v(i)) +
     $       ahf*(rl+one)/rzone(i)**2)
         z(i) = rzone(i)/r(np1)
  10  continue
      call pltdmp(z,nmax,n,'rz  ')
c
      drz(1) = rzone(1) - rz0
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      drz(np1) = p(n)*v(n)/gor(np1)
      pi4g = forpi*g
      if( lval .le. 0 ) goto 50
c
c          initialize the matrices.
c
c          inner boundary conditions.
c
      afac = forpi*r(1)**3/dm2(1)
      drdm01 =-ahf*(rl+one)/rz0**2
      drdm02 =-ahf*(rl+one)/rz0**2
      ag11 = drdm01*g1(1)*p0*(ahf*rl/rho0 - afac)
      ag1(1,1) = zero
      ag1(1,2) =-gor(1)/r(1) +
     $       drdm(1,1)*g1(1)*p(1)*(ahf*rl*v(1) + afac) +
     $          drdm02*g10*p0*(ahf*rl/rho0 - afac)
      ag1(1,3) = drdm(1,2)*g1(1)*p(1)*(ahf*rl*v(1) + afac) + ag11
      ag3(1,1) = zero
      ag3(1,2) = rl
      ag4(1,1) = zero
      ag41 = rl1*(g10*p0*(ahf*rl/rho0 - afac))/rz0**2
      ag4(1,2) = rkl2(1)*g1(1)*p(1)*(ahf*rl*v(1) + afac) +
     $      rl1*gor(1)/r(1) + ag41
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
         ag1(i+1,2) =-for*gor(i+1)/r(i+1)+pi4g*ahf*(rho(i+1)+rho(i))+
     $     drdm(i+1,1)*g1(i+1)*p(i+1)*(ahf*rl*v(i+1) + afac) +
     $     drdm(i,2)*g1(i)*p(i)*(ahf*rl*v(i) - afac )
         ag1(i+1,3) = drdm(i+1,2)*g1(i+1)*p(i+1)*(ahf*rl*v(i+1)+afac)
         ag3(i+1,1) = ahf*rl - r(i+1)/drz(i+1)
         ag3(i+1,2) = ahf*rl + r(i+1)/drz(i+1)
         ag4(i+1,1) = rkl2(i)*(g1(i)*p(i)*(ahf*rl*v(i) - afac)) +
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
      fdr = one + ahf*rl*(drz(np1)/r(np1))
      ag1(np1,2) =-for*gor(np1)/r(np1) + vl*drdm(n,2) -
     $	  pi4g*rho(n)*ahf*fdr*f
      ag1(np1,3) = zero
      ag3(np1,1) =-(rl+one)/(one+dlnr)
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
      dlrhdr =-(rho(n)*gor(n+1)*r(n+1))/p(n)*0.6d0
      ap1(n,2) =-pi4g*rho(n)*( (ddr13*two*rl*rzone(n) + ap3n3)*f +
     $    rzone(n)**2*drdm(n,2) - dlrhdr )
c
      ap3(n,1) = ap3n1 + ddr11*two*rl*rzone(n)
      fp3 = (rl+ahf)*(drz(np1)/r(np1)) - one
      ap3(n,2) =-(ap3n1+ap3n3) + ddr12*two*rl*rzone(n) -
     $	  (ddr13*two*rl*rzone(n) + ap3n3)*fp3/(one+dlnr)
      ap3(n,3) = zero
c
      ap4(n) =-pi4g*rl1*rho(n)
      ah4(np1) = zero
c
c        write the matrices to file 6.
c
      if( iout .gt. 1 ) then
	 write(11,2200)
	 do 27 i=1,np1
	    write(11,2201) i,ag1(i,1),ag1(i,2),ag1(i,3),drdm(i,1),
     $	       drdm(i,2),ag3(i,1),ag3(i,2),ag4(i,1),ag4(i,2),ah4(i)
  27	 continue
	 write(11,2300)
	 do 28 i=1,n
	    write(11,2201) i,ah1(i,1),ah1(i,2),ah3(i),ap1(i,1),
     $	       ap1(i,2),ap3(i,1),ap3(i,2),ap3(i,3),ap4(i)
  28	 continue
      endif
c
c          convergence to a normal mode is governed by the momentum
c       equation for the outermost zone. as in the radial case, this
c       function (divided by any x value) is used in a secant method
c       for controlling the iterations. as the equations are not
c       nonlinear (as in non-non) in the frequency, it is hoped the
c       calculation is, in some sense, more stable.
c
c                                               10/4/83 wd pesnell
c
      write(11,3001)
      xnorm = r(np1)**2
      iroot = 0
      dom = (omhigh - omlow)/dfloat(nomega - 1)
      do 30 i=1,nomega
         omsq = omlow + dfloat(i-1)*dom
	 fomeg = fomsq(omsq,x,xnorm,n)
	 if( iout .gt. 0 ) write(11,4400) i,omsq,fomeg
	 if( i .gt. 1 ) then
	    if( fomeg*fomi .lt. zero ) then
	       iroot = iroot + 1
	       stomeg(iroot) = omsq - dom*fomeg/(fomeg-fomi)
	    endif
	 endif
         fomi = fomeg
  30  continue
      if( iroot .eq. 0 ) goto 900
      write(11,3000) iroot,lval,(stomeg(i),i=1,iroot)
      write(6,3000) iroot,lval,(stomeg(i),i=1,iroot)
      write(1,3000) iroot,lval,(stomeg(i),i=1,iroot)
c
c          attempt to converge on each of the frequencies found in
c       discriminant search. a secant method to find the zeroes of
c       error function is used.
c
      omsq1 = zero
      erra  = zero
      do 100 nroot=1,iroot
         write(11,7900) nroot
         acctst = accur
         omsq = stomeg(nroot)
         do 45 icount=1,30
	    err = fomsq(omsq,x,xnorm,n)
	    if( icount .eq. 1 ) then
	       omsq1 = omsq
	       erra  = err
	       omsq  = omsq*(one + 1.d-7)
	    else
	       afac = dabs((omsq-omsq1)/omsq)
	       if(iout.gt.0) write(11,4400) icount,omsq,afac,err,erra
	       if( icount .ge. 20 ) acctst = acctst*two
	       if( afac .le. acctst ) goto 47
	       omsq2 = (erra*omsq-err*omsq1)/(erra-err)
	       domsq = omsq2 - omsq
	       domsq = dsign(dmin1(dabs(domsq),qrt*dabs(omsq)),domsq)
	       omsq1 = omsq
	       erra  = err
	       omsq  = omsq + domsq
	    endif
  45	 continue
c
c          no convergence for this root, skip to next one.
c
         write(1,4500) nroot,omsq,accur,acctst
         write(6,4500) nroot,omsq,accur,acctst
         write(11,4500) nroot,omsq,accur,acctst
         goto 100
c
c          converged to a root, calculate the eigenvectors.
c
  47     continue
         x(3*n+1) = xnorm
c
         do 40 i=1,n
            dr(i) = (r(i)/r(np1))**lval*(x(3*i-2)/r(i)**2)
            dh(i) = (x(3*i)/rzone(i)**2)*(rzone(i)/r(np1))**lval
            adrho(i) = (drdm(i,1)*x(3*i-2)+drdm(i,2)*x(3*i+1) +
     $             rkl2(i)*x(3*i) )*(rzone(i)/r(np1))**lval
            dp(i) = g1(i)*adrho(i)
            gam(i) = x(3*i-1)*(rzone(i)/r(np1))**lval
  40     continue
         dh0 = dr(1)/rl
         gam0 = x(2)*(rz0/r(np1))**lval
         drho0 = adrho(1)*(rz0/rzone(1))**lval
         dp0 = drho0*g1(1)
         dr(n+1) = one
         dlnr = ahf*(rl+one)*drz(np1)/r(np1)
         gam(n+1) = gam(n)*(one-dlnr)/(one+dlnr)
         dh(n+1) = zero
         qdr0 = dr(2)*(r(1)/r(2))**(lval-2)
c
c          weight function for the adiabatic oscillation.
c
         wgrav(1) = dr(1)*r(1)*dm2(1)*rl*gam(1)
         wdiag(1) = drhdr(1)*gor(1)*r(1)*dr(1)**2*dm2(1)
         stwait(1) = wgrav(1) + wdiag(1)
         rsum = (dr(1)*r(1))**2*dm2(1)
         hsum = zero
         rbar = r(1)*(dr(1)*r(1))**2*dm2(1)
         hbar = zero
         crsum = zero
         do 48 i=1,n
            wtherm(i) = g1(i)*p(i)*v(i)*dm1(i)*adrho(i)**2
            rsum = rsum + (dr(i+1)*r(i+1))**2*dm2(i+1)
            hsum = hsum + (dh(i)*rzone(i))**2*dm1(i)
            rbar = rbar + r(i+1)*(dr(i+1)*r(i+1))**2*dm2(i+1)
            hbar = hbar + rzone(i)*(dh(i)*rzone(i))**2*dm1(i)
            crsum = two*dr(i+1)*dm2(i+1)*r(i+1)**2*
     $         rl1*ahf*(dh(i)+dh(i+1))
            drhoel = adrho(i) - ahf*
     $                (drhdr(i)*dr(i) + drhdr(i+1)*dr(i+1))
            wgrav(i+1) = dm1(i)*gam(i)*drhoel
            wdiag(i+1) = drhdr(i+1)*gor(i+1)*r(i+1)*
     $         dr(i+1)**2*dm2(i+1)
            wdiag(i+1) = two*gor(i+1)*r(i+1)*dr(i+1)*dm2(i+1)*
     $         rl1*ahf*(dh(i)+dh(i+1)) +
     $      ( pi4g*ahf*(rho(i)+rho(i+1)) -
     $         for*gor(i+1)/r(i+1) )*(dr(i+1)*r(i+1))**2*dm2(i+1)
            weight(i) = wtherm(i) + wgrav(i+1) + wdiag(i+1)
            stwait(i+1) = stwait(i) + weight(i)
  48     continue
c
c          normalize the weight function.
c
         rke = rsum + rl1*hsum
         if( rke .le. zero ) rke = one
c
c          rbar is the effective radius of the pulsation mode.
c
         rbar = (rbar + rl1*hbar)/(rke*r(n+1))
         do 46 i=1,n
            stwait(i) = stwait(i)/(omsq*rke)
            weight(i) = weight(i)/(omsq*rke)
            wtherm(i) = wtherm(i)/(omsq*rke)
            wgrav(i) = wgrav(i)/(omsq*rke)
            wdiag(i) = wdiag(i)/(omsq*rke)
  46     continue
         stwait(np1) = (stwait(np1)/rke)
         qch = dabs((omsq-stwait(np1))/omsq)
c
c          evaluate the rotational splitting coefficient.
c
         crsum = crsum/rke
         hrsum = hsum/rke
         ckl = crsum + hrsum
c
c          find the analytic eigenvalue and compare to the
c       one calculated in the mesh.
c
         period = twopi/dsqrt(dabs(omsq))
         qomsq = omsq/(pi*g*rhom)
         index = knode(dr,dh,nmax,n,np,ng)
         write(6,5504) nroot,lval,omsq,stwait(np1),qch,period
         write(6,5505) index,np,ng,qomsq,dr(1),qdr0,rbar
         write(1,5504) nroot,lval,omsq,stwait(np1),qch,period
         write(1,5503) index,np,ng,qomsq,dr(1),qdr0,rbar
         write(1,5506) ckl,crsum,hrsum
         write(11,5500) nroot,lval,omsq,stwait(np1),qch,period
         write(11,5503) index,np,ng,qomsq,dr(1),qdr0,rbar
         write(11,5506) ckl,crsum,hrsum
         stwait(np1) = stwait(np1)/omsq
         wgrav(np1) = wgrav(np1)/(omsq*rke)
         wdiag(np1) = wdiag(np1)/(omsq*rke)
	 if( iout .gt. 0 ) then
	    write(11,5502)
	    write(11,5501) index,dr(1),dh0,drho0,gam0,dp0
	    write(11,5501) (i,dr(i+1),dh(i),adrho(i),gam(i),dp(i),
     $			    weight(i),stwait(i),i=1,n)
	 endif
	 iv = nodes(dr,np1,1)
         iv = nodes(dh,n,1)
	 call orthog(dr,dh,n,omsq,index,1)
	 call pltdmp(dr,nmax,np1,'dr/r')
	 call pltdmp(dh,nmax,n,'dh/h')
	 call pltdmp(adrho,nmax,n,'drho')
	 call pltdmp(dp,nmax,n,'dp  ')
	 call pltdmp(gam,nmax,n,'gam ')
	 call pltdmp(weight,nmax,n,'wait')
	 call pltdmp(stwait,nmax,np1,'wint')
	 call pltdmp(wtherm,nmax,n,'wthr')
	 call pltdmp(wdiag,nmax,np1,'wdia')
	 call pltdmp(wgrav,nmax,np1,'wgrv')
	 call eigenf(omsq,n,lval,rke)
c
 100  continue
      call orthog(dr,dh,n,omsq,index,-lval)
      return
  50  continue
      ihmin = nomega
      if( ihmin.le.0 .or. ihmax.le.0 ) return
c
c          pulsation equation solver
c       basic reference -- castor, ap.j. 166 (109) 1971.
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
      afacp = forpi*r(np1)**2/dsqrt(dm2(np1))
      ag1(np1,1) =-afacp*g1p*drdm(n,1)
      ag1(np1,2) =-for*gor(np1)/r(np1)-afacp*g1p*drdm(n,2)
      ag1(np1,3) = zero
c
c          write the matrices to file tape11, if io is less
c       than 3, this write is not done.
c
      if( iout .gt. 0 ) then
	 write(11,7000)
	 do 75 i=1,np1
	    write(11,7001) i,ag1(i,1),ag1(i,2),ag1(i,3)
  75	 continue
      endif
c
c          calculate the acoustic tavel time from surface to the
c       innermost zone. this is stored in transt.
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
	 call trisol(ag1,1,n,omsq,yo,xo,nmax,y,z)
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
	    if(iout.ge.1) write(11,1400) icount,omsq,omsql,afacp
	    if( afacp .lt. 1.d-12 ) goto 86
	    if( icount.gt.20 .and. afacp.lt.1.d-10 ) goto 86
	     omsq2 = (erra*omsq-err*omsql)/(erra-err)
	     omsql = omsq
	     domsq = omsq2-omsq
	     domsq = dsign(dmin1(dabs(domsq),1.1d0*dabs(omsq)),domsq)
	     omsq  = omsq + domsq
	     erra  = err
	 endif
  81  continue
c
c          no convergence in the adiabatic eigenvalue.
c
      write(11,8100) omsq,iterad,iv,ih
      omsqp = zero
      omsqc = tcon*dfloat(ih+1)
      goto 79
c
c          converged to omsq value, check if positive and with
c      with proper number of nodes in the displacement eigen-
c      vector. if either is not true, try again.
c
  86  continue
      if( omsq .le. zero ) then
c
c          negative omsq value, stop working on this mode and
c       move on.
c
	 write(11,8600) omsq,iv,ih
	 omsqp = zero
	 omsqc = tcon*dfloat(ih+1)
	 goto 79
      endif
c
      do 90 i=1,n
         yo(i) = zero
  90  continue
      xo(np1) = x0np
      yo(n)  =-ag1(n,3)*x0np
      call trisol(ag1,1,n,omsq,yo,xo,nmax,y,z)
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
c          not converged to a frequency with right number of nodes,
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
         weight(i) =-for*gor(i+1)*xo(i+1)**2/r(i+1) +
     $      p(i)*v(i)*g1(i)*adrho(i)**2*dm1(i)
         stwait(i+1) = stwait(i) + weight(i)
         rke = rke + xo(i+1)*xo(i+1)
 200  continue
      dr(np1) = one
      qomsq = omsq/(pi43*g*rhom)
      stwait(np1) = stwait(np1)/rke
      qch = dabs((stwait(np1)-omsq)/omsq)
      period = twopi/dsqrt(dabs(omsq))
      write(1,1001) omsq,qomsq,iv,ih,stwait(np1),qch,period
      write(1,1002) dr(1)
      write(6,1001) omsq,qomsq,iv,ih,stwait(np1),qch,period
      write(11,1001) omsq,qomsq,iv,ih,stwait(np1),qch,period
      iv = nodes(dr,np1,1)
c
c          normalize the weight function per zone to omsq.
c
      do 105 i=1,n
         weight(i) = weight(i)/(omsq*rke)
         stwait(i) = stwait(i)/(omsq*rke)
 105  continue
c
      if(iout.ge.2) write(11,1000) omsq,qomsq,(i,dr(i+1),adrho(i),
     $    dp(i),weight(i),stwait(i),i=1,n)
      call orthog(xo,dh,-n,omsq,iv,1)
      call pltdmp(dr,nmax,np1,'dr/r')
      call pltdmp(adrho,nmax,n,'drho')
      call pltdmp(dp,nmax,n,'dp  ')
      call pltdmp(weight,nmax,n,'wait')
      omsqp = omsq
      omsqc = omsq
  79  continue
      call orthog(dr,dh,-n,omsq,index,-1)
      return
c
c          error recovery.
c
 900  continue
      write(6,9000) lval,omlow,omhigh
      write(11,9000) lval,omlow,omhigh
      return
c
 2200 format(1h1,//,21h ag1,drdm,ag3,ag4,ah4)
 2201 format(1x,i5,1p10e12.4)
 2300 format(1h1,//,20h ah1,ah3,ap1,ap3,ap4)
c
 3001 format(1h1,1x,26hbegin discriminant search.)
 3000 format(1x,i4,23h roots for analysis: l=,i4,/,(3x,1p5e15.6))
 4400 format(1x,i5,1p4e15.6)
 4500 format(1x,28hno convergence root number ,i4,11h frequency=,1pe15.6,
     $   /,3x,13haccuracy was ,e10.3,8h is now ,e10.3)
 5500 format(//,5x,4hroot,i5,3h l=,i4,7h omsq =,1pe15.6,
     $   /,2x,37h eigenfrequency from weight function=,e12.3,
     $	 7h error=,e12.3,/,1x,10h period = ,e11.4,8h (secs) )
 5501 format(1x,i5,1p7e15.6)
 5502 format(4x,1hi,8x,4hdr/r,10x,4hdh/h,10x,8hdrho/rho,8x,5hgamma,
     $   11x,4hdp/p,10x,6hweight)
 5503 format(23h radial quantum number=,i5,25h number of nodes, p-type=,
     $	 i5,8h g-type=,i5,/,1x,24hnormalized eigenvalue = ,1pe12.4,
     $	 /,1x,8hdr(1) = ,e10.3,16h scaled dr(1) = ,e10.3,/,
     $	 1x,28h effective radius of mode = ,e11.4)
 5504 format(//,5x,4hroot,i5,3h l=,i4,7h omsq =,1pe15.6,
     $   /,2x,38h eigenfrequency from weight function =,e12.4,
     $	 7h error=,e12.3,/,1x,10h period = ,e11.4,8h (secs) )
 5505 format(1x,4h k =,i5,9h p-type =,i5,9h g-type =,i5,
     $	 1x,9h qomsq = ,1pe12.4,
     $	 /,1x,8hdr(1) = ,e10.3,16h scaled dr(1) = ,e10.3,/,
     $	 1x,28h effective radius of mode = ,e11.4)
 5506 format(1x,31h rotational splitting, total = ,1pe11.4,
     $   5h a = ,e11.4,5h b = e11.4)
c
 1400 format(1x,i3,1p7e11.4)
c
 7000 format(1h1,/,5x,26h l=0 pulsation matrix ag1.,/)
 7001 format(1x,i4,1x,1p3e10.3)
 7800 format(8h iterad=,i3,4h iv=,i3,4h ih=,i3,6h omsq=,1pe12.4,
     $   7h omsqc=,e12.4,7h omsqp=,e12.4)
 7900 format(1h1,//,20x,28h***** begin search for mode ,i2,6h *****)
 8000 format(40h cant find no. of nodes.eq.ih-1...nodes=,i3,4h ih=,i3,
     $ 6h omsq=,1pe13.6)
 8600 format(1x,13hadiabatic f2=,1pe17.10,18h is .lt. 0 for iv=,
     $  i3,4h ih=,i3)
 8100 format(1x,13hadiabatic f2=,1pe15.6,20h not converged after,i4,
     $   16h iterations, ih=,i4,17h number of nodes=,i4)
 1000 format(1h1,19h adiabatic solution,/,8h  omsq =,1pe17.9,5x,
     $  6hqomsq=,e11.3,//,3x,1hi,5x,4hdr/r,10x,4hdrho,11x,
     $	2hdp,11x,6hweight/(1x,i3,1p5e14.6))
 1001 format(/,31h linear, adiabatic, radial mode,/,
     $ 8h  omsq =,1pe16.9,2x,7hqomsq =,e12.4,4h iv=,i3,4h ih=,i3,/,
     $ 37h eigenfrequency from weight function=,e12.4,7h error=,e12.4,
     $	 /,1x,10h period = ,e11.4,8h (secs) )
 1002 format(1x,9h dr(1) = ,1pe10.3)
c
 9000 format(43h from lnanon...no adiabatic roots found, l=,i4,/,
     $   17h frequency range:,1p2e15.6)
 9999 format(1x,1p3e12.3)
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
      if( n .eq. 1 ) return
      nm1 = n - 1
      do 10 i=2,nm1
         if( x(i)*x(i+1) .ge. zero ) goto 10
         if( y(i)*y(i+1) .lt. zero ) write(11,1000)
	 if( x(i+1) .gt. x(i) ) then
	    if( y(i) .ge. zero ) ng = ng + 1
	    if( y(i) .lt. zero ) np = np + 1
	 else
	    if( y(i+1) .ge. zero ) np = np + 1
	    if( y(i+1) .lt. zero ) ng = ng + 1
	 endif
  10  continue
      if( x(n-1)*x(n) .lt. zero ) write(11,1001)
      knode = np - ng
      nps = np
      ngs = ng
      return
 1000 format(1x,50hfrom knode...quadrant jump default rotations used.)
 1001 format(1x,51hfrom knode...node at outer boundary is not counted.)
      end
      function fomsq(omsq,x,xnorm,n)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=600 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag3(nmax,2),ag4(nmax,2),
     $		     ah1(nmax,2),ah3(nmax),  ah4(nmax),
     $		     ap1(nmax,2),ap3(nmax,3),ap4(nmax)
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
      call rbmles(atrix,3*nmax,1,3*n,7,x)
      err = ag1(n+1,1)*x(3*n-2) + (ag1(n+1,2)-omsq)*xnorm +
     $      ag3(n+1,1)*x(3*n-1) + ag4(n+1,1)*x(3*n)
      fomsq = err/x(1)
      return
      end
      subroutine rbmles(a,nmax,imin,imax,m,y)
      implicit real*8(a-h,o-z)
c
c          linear equation solver. (r)eal (b)and (m)atrix (l)inear
c       (e)quation (s)olver. solves a*x = y, destroying a and
c       returns the value of x in y.
c
      dimension a(nmax,m),y(nmax)
      id = imax-imin+1
      idm = id-1
      mmid = (m+1)/2
      mmm = mmid-1
      do 20 ii=1,idm
         i = imax+1-ii
         den = 1.0d0/a(i,mmid)
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
      subroutine trisol(a,imin,imax,x,y,z,nmax,e,f)
      implicit real*8(a-h,o-z)
c
      dimension a(nmax,3), y(nmax), z(nmax)
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
c          finds no. of nodes(sign changes) in nz elements of array w
c
      dimension w(nz),inode(100)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      nodes = 0
      if( nz .le. 4 ) return
      nz1 = nz-1
      nodest = 0
      do 10 i=3,nz1
	 if( w(i-1)*w(i) .lt. zero ) then
c
c       zero crossing found. is it general enough to be a node?
c
	    if( w(i-2)*w(i-1) .gt. zero ) then
	       if( w(i)*w(i+1) .gt. zero ) then
		  nodest = nodest + 1
		  inode(nodest) = i-1
	       endif
	    endif
	 endif
  10  continue
      inode(nodest+1) = nz
      nodes = nodest
      if( lout .eq. 0 ) return
      write(11,1000) nodest
      if( nodes .le. 0 ) return
      write(11,1001) (inode(i),i=1,nodest)
 1000 format(/,1x,i3,12h nodes found)
 1001 format(1x,7hnodes =,30i4)
      return
      end
      subroutine orthog(dr,dh,n,omsq,index,imode)
      implicit real*8(a-h,o-z)
c
c          check for orthogonality of the wave functions. part
c      one stores the vertical and horizontal wavefunctions,
c      part two does the integrations.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),qeta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/eigstr/ sdr(nmax,4),sdh(nmax,4),stome(4),iomeg(4)
      dimension dr(nmax),dh(nmax),rke(4),cross(4,4)
      common/stnmod/ nmode
c
      if( n .lt. 0 ) goto 100
      if( imode .le. 0 ) goto 50
      if( nmode .ge. 4 ) return
      if( nmode .eq. 0 ) goto 6
      do 5 i=1,nmode
         if( dabs((omsq-stome(i))/stome(i)) .gt. 1.d-6 ) goto 5
          return
   5  continue
   6  continue
      nmode = nmode + 1
      stome(nmode) = omsq
      iomeg(nmode) = index
      do 10 i=1,n
         sdr(i,nmode) = dr(i)
         sdh(i,nmode) = dh(i)
  10  continue
      sdr(n+1,nmode) = dr(n+1)
      return
  50  continue
c
c          part two: the integrals.
c
      if( nmode .le. 1 ) return
      lval = iabs(imode)
      rl1 = dfloat(lval)*dfloat(lval+1)
      do 60 j1=1,nmode
         rke(j1) = zero
         do 60 j2=1,nmode
            cross(j2,j1) = sdr(1,j2)*sdr(1,j1)*dm2(1)*r(1)**2
  60  continue
      do 70 j1=1,nmode
         do 70 j2=1,nmode
            do 70 i=1,n
               cross(j2,j1) = cross(j2,j1) +
     $                rl1*sdh(i,j2)*sdh(i,j1)*dm1(i)*rzone(i)**2 +
     $                sdr(i+1,j2)*sdr(i+1,j1)*dm2(i+1)*r(i+1)**2
  70  continue
      do 75 j1=1,nmode
         rke(j1) = dsqrt(cross(j1,j1))
  75  continue
      do 80 j1=1,nmode
         do 80 j2=1,nmode
            cross(j2,j1) = cross(j2,j1)/(rke(j1)*rke(j2))
  80  continue
      write(1,8000) lval
      write(11,8000) lval
      do 85 j1=1,nmode
         do 85 j2=1,nmode
            write(11,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
            write(1,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
  85  continue
      return
c
c          special case of l=0.
c
 100  continue
      nuse =-n
      if( imode .le. 0 ) goto 150
      if( nmode .ge. 4 ) return
      nmode = nmode + 1
      do 110 i=1,nuse
         sdr(i,nmode) = dr(i)
 110  continue
      sdr(n+1,nmode) = dr(n+1)
      stome(nmode) = omsq
      iomeg(nmode) = index
      return
 150  continue
c
c          part two: the integrals.
c
      if( nmode .le. 1 ) return
      do 160 j1=1,nmode
         rke(j1) = zero
         do 160 j2=1,nmode
            cross(j2,j1) = sdr(1,j2)*sdr(1,j1)
 160  continue
      do 170 j1=1,nmode
         do 170 j2=1,nmode
            do 170 i=1,nuse
               cross(j2,j1) = cross(j2,j1) + sdr(i+1,j2)*sdr(i+1,j1)
 170  continue
      do 175 j1=1,nmode
         rke(j1) = dsqrt(cross(j1,j1))
 175  continue
      do 180 j1=1,nmode
         do 180 j2=1,nmode
            cross(j2,j1) = cross(j2,j1)/(rke(j1)*rke(j2))
 180  continue
      write(1,1800)
      write(11,1800)
      do 185 j1=1,nmode
         do 185 j2=1,nmode
            write(11,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
            write(1,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
 185  continue
      return
c
 1800 format(1x,34h radial cross integrals for modes:)
 8000 format(1x,4hl = ,i3,27h cross integrals for modes:)
 8001 format( 2(3x,1h<,i3,1h,,i3,1h>,1p2e12.5) )
      end
      subroutine cjhdmp(nin,y02)
      implicit real*8(a-h,o-z)
      character*1 iyorn
c
c          write out the input file for carl hansen's runge-kutta
c       pulsation analysis.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),eta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/corevl/ rz0,p0,rho0,g10
      dimension xi(nmax),opv(nmax)
      data lnafil/41/
c
      write(6,1100)
      read(5,1101,end=900) iyorn
      if( iyorn .ne. 'y' .and. iyorn .ne. 'Y' ) goto 900
c
      open (unit=lnafil,file='cjhinp',status='unknown')
c
c      write information to file 41.
c
      n = nin
      np1 = n + 1
      p(np1) = zero
      vi = zero
      opv(1) = 1.d0
      do 20 i=2,n
         opv(i) = 1.d0/(one+vn(i))
         xi(i) = dlog(r(i)/p(i))
  20  continue
  10  continue
      xi(1) = xi(2) - 0.2d0
      qmass = rm(n)/1.991d33
      write(lnafil,1001) y02,n-1,qmass,n-1
      write(lnafil,1000) (xi(i),r(i),gor(i),rho(i),i=2,n)
      write(lnafil,1002) n-1
      do 30 i=2,n
         write(lnafil,1000) xi(i),gor(i)/r(i),vn(i)/g1(i),ra(i)
         write(lnafil,1000) un(i),opv(i)
  30  continue
      close(unit=lnafil)
      return
 900  continue
      return
c
 1000 format(1p4e20.12)
 1001 format(1pe12.5,i5/1pe12.4/i4)
 1002 format(i5)
 1100 format(1x,'Do you want the eulerian code dump? (y/n): ',$)
 1101 format(a1)
      end
      subroutine eigenf(omsq,nin,lin,rke)
      implicit real*8(a-h,o-z)
c
c          find the normalized eigenvectors and the weight functions
c       as a function of position.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),theta(nmax),dthdr(nmax),v(nmax),rm(nmax),
     $	    gor(nmax),ceta(nmax),dm1(nmax),dm2(nmax),bv(nmax),vn(nmax),
     $	    un(nmax),drhdr(nmax),ra(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),rho(nmax),rzone(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/corevl/ rz0,p0,rho0,g10
      common/linear/ spac(nmax,3),dr(nmax),dh(nmax),gam(nmax),
     $	    qnyr(nmax),gyr(nmax),dp(nmax),adrho(nmax),weight(nmax),
     $	    cyr(nmax),xo(nmax),yo(nmax),qwait(nmax),y1(nmax),
     $	    y2(nmax),y3(nmax),y4(nmax),stwait(nmax),spac1(nmax)
      dimension vq(6)
c
      nsurf = nin
      np1 = nsurf + 1
      lval = lin
      rl = dfloat(lval)
      rl1 = rl*(rl+one)
      do 10 i=1,nsurf
         xo(i) = dlog(rzone(i)/p(i))
         y1(i) = ahf*(dr(i) + dr(i+1))
         y2(i) = omsq*dh(i)/(ahf*(gor(i)/r(i)+gor(i+1)/r(i+1)))
         y3(i) = gam(i)/(ahf*(gor(i)*r(i)+gor(i+1)*r(i+1)))
         y4(i) = ((gam(i+1)-gam(i))/(rzone(i+1)-rzone(i)))/gor(i)
  10  continue
c
c          weight functions
c
      stwait(1) = zero
      do 20 i=1,nsurf
         vq(1) = ahf*(gor(i)/r(i) + gor(i+1)/r(i+1))
         vq(2) = ahf*(vn(i) + vn(i+1))/g1(i)
         vq(3) = ra(i)
         vq(4) = ahf*(un(i)+un(i+1))
         vq(5) = one/(one + vq(2)*g1(i))
         vq(6) = rzone(i)
	 if( i .eq. 1 ) then
	    dx = xo(2)-xo(1)
	 elseif( i .eq. nsurf ) then
	    dx = xo(nsurf) - xo(nsurf-1)
	 else
	    dx = ahf*(xo(i+1)-xo(i-1))
	 endif
	 temp = forpi*dx*rho(i)*vq(1)*vq(5)*vq(6)**5
         cyr(i) = temp*vq(2)*(y2(i) - y3(i))**2
         qnyr(i) = temp*ra(i)*y1(i)**2
         gyr(i) =-temp*(y4(i) + (rl+one)*y3(i))**2/vq(4)
         weight(i) =  cyr(i) + qnyr(i) + gyr(i)
         stwait(i+1) = stwait(i) + weight(i)
  20  continue
      do 30 i=1,nsurf
	  weight(i) = weight(i)/(omsq*rke)
          stwait(i) = stwait(i)/(omsq*rke)
          cyr(i) = cyr(i)/(omsq*rke)
          qnyr(i) = qnyr(i)/(omsq*rke)
          gyr(i) = gyr(i)/(omsq*rke)
  30  continue
      stwait(np1) = stwait(np1)/rke
      qch = dabs((omsq-stwait(np1))/omsq)
      qrke = rke*omsq/2.0d0
      write(6,2000) omsq,stwait(np1),qch,qrke
      write(1,2000) omsq,stwait(np1),qch,qrke
      write(11,2000) omsq,stwait(np1),qch,qrke
      call pltdmp(weight,nmax,nsurf,'twat')
      stwait(np1) = stwait(np1)/omsq
      call pltdmp(stwait,nmax,np1,'twnt')
      call pltdmp(cyr,nmax,nsurf,'cyr ')
      call pltdmp(qnyr,nmax,nsurf,'nyr ')
      call pltdmp(gyr,nmax,nsurf,'gyr ')
      return
c
 2000 format(2x,32heigenvalue from matrix solution=,1pe12.3,
     $     /,2x,32heigenvalue from weight function=,e12.3,
     $       7h error=,e12.3,/,
     $	 1x,21hkinetic energy amp. =,e12.4)
      end
