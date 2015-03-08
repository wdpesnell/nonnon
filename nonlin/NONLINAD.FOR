c$large
c$notruncate
c$storage:2
      program nonlin
      implicit real*8(a-h,o-z)
      character*80 title
      logical iprt
c
c          The Lagrangian radial adiabatic nonlinear analysis
c       for stellar pulsations. The notation is lifted from the model
c       code and radial stability analysis. Begun July 1987.
c
c                                                7/23/87 WD Pesnell
c
      parameter ( nmax=200, nmaxtd=5000 )
      common/phypar/ rn(nmax),tn(nmax),vn(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn(nmax),un(nmax),g1(nmax),g3m1(nmax),
     $               cs(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),p(nmax),up(nmax)
      common/scrtch/ ag1(nmax,3),xo(nmax),yo(nmax),w(nmax),glamp(nmax),
     $       spac1(nmax,7)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      dimension r2(nmaxtd), rmax(nmaxtd), umax(nmaxtd), sttime(nmaxtd)
c
c          files in this program.
c  unit   routine    name                purpose
c    5      main    sys$input           terminal input file
c    6      main    sys$output          terminal output file
c    9      main    sumout              summary output file
c
c   10      main    lnrad               contains the model quantities
c                (logical to nondmp.Dat)   for use by the stability
c                                          analyses.
C   11      main    casout              main output file
c
c   12     pltint   casplt              has vectors of eigenvectors
c                                          quantities for pltcmd program
c
      open (unit=11,file='casout',status='new')
c
c          read in the input variables from file 10.
c
      open (unit=10,file='lnrad',status='old',readonly)
      read(10,1001) title
      read(10,*) npts,rlumgv,totmas,teff,rphoto,corlum
      read(10,*) irad,noburn,onemq0,onemq1,rlums,rlumc
      np1 = npts+1
      read(10,1000) (rn(i),i=1,np1)
      read(10,1000) (tn(i),i=1,npts)
      read(10,1000) (vn(i),i=1,npts)
      read(10,1000) (dum,i=1,npts)
      read(10,1000) (dum,i=1,npts)
      read(10,1000) (dum,i=1,npts)
      read(10,1000) (dm1(i),i=1,npts)
      read(10,1000) (dum,i=1,npts)
      read(10,1000) (dm2(i),i=1,np1)
      read(10,1000) (rm(i),i=1,np1)
      read(10,1000) (pn(i),i=1,npts)
      read(10,1000) (g1(i),i=1,npts)
      read(10,1000) (g3m1(i),i=1,npts)
      close (unit=10)
c
c          iout:= 0 for little output (no eigenvectors)
c                 1 for eigenvectors
c                 2 for eigenvectors, matrices and conver. Info.
c
c          ihmin:= number of nodes+1 in lowest eigenvector.
c
c          ihmax:= number of nodes+1 in highest eigenvector.
c
      write(11,2000) title,npts
      write(11,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(11,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
      write(11,2003) iout,ihmin,ihmax
      write(11,2004)
      write(11,2005) (i,rn(i+1),tn(i),vn(i),pn(i),g1(i),
     $                  rm(i+1),i=1,npts)
c
      call pltint(nmax,npts)
      iout = 3
      ihmax = 1
      ihmin = 1
      call linrad( iout, ihmin, ihmax, omega, period)
      write(6,9000)
      read(5,*,end=900) ntime,nprt,amp,dtfac
      if( ntime .le. 0 ) goto 900
      if( ntime .gt. nmaxtd ) ntime = nmaxtd
      if( amp .eq. 0.0 ) amp = 10.D9/(omega*rn(np1))
      up(1) = zero
      rp(1) = zero
      do 90 i=1,npts
          vp(i) = vn(i)
          tp(i) = tn(i)
          p(i) = pn(i)
c         g1(i) = 5.D0/3.D0
c         g3m1(i) = 2.D0/3.D0
         cs(i) = sqrt(g1(i)*pn(i)*vn(i))
         un(i+1) =-omega*rn(i+1)*w(i+1)*amp
  90  continue
      dtmax = period/1.D2
      eltime = zero
      Rinit = rn(npts)
      do 100 itime=1,ntime
         iprt = .false.
         If( mod(itime,nprt) .eq. 0 ) iprt = .true.
c
c          Transfer the last time step to storage.
c
         do 110 i=1,npts
            rp(i+1) = rn(i+1)
            up(i+1) = un(i+1)
c            vp(i) = vn(i)
c            tp(i) = tn(i)
c            p(i) = pn(i)
 110     continue
         dtmin = qtime(rn,cs,nmax,npts,iprt)*dtfac
         dt = dmin1(dtmax,dtmin)
         call hydro(iout, dt, itime, ierrhy)
         if( ierrhy .lt. 0 ) goto 910
         eltime = eltime + dt
         phase = eltime/period
         sttime(itime) = phase
         rmax(itime) = rn(npts)/Rinit - 1.d0
         r2(itime) = rn(2)/Rinit - 1.d0
         umax(itime) = un(npts)*1.d-5
         if( iprt ) then
            write(11,2006) (i,rn(i+1),un(i+1),tn(i),vn(i),
     $       pn(i),i=npts-5,npts)
            write(6,2006) itime,phase,rn(1),rn(2),rn(npts),rn(np1)
         endif
 100  continue
 910  continue
      call pltdmp(sttime,nmaxtd,itime-1,'time')
      call pltdmp(rmax, nmaxtd, itime-1,'rt  ')
      call pltdmp(r2, nmaxtd, itime-1,'r2  ')
      call pltdmp(umax, nmaxtd, itime-1,'ut  ')
 900  continue
      stop
c
 1000 format(1p,4e20.13)
 1001 format(a)
 1010 format(1p,4e20.13)
 1011 format(1x,'Enter iout, ihmin, and ihmax. ',$)
c
 2000 format(1h1,/,1x,19houtput from nonlin.,/,
     $      50H  initial model in hydrostatic and thermal balance,//,
     $      2x,a,1x,i5,6h zones,/)
 2001 format(4x,18h total luminosity=,1pe11.4,12H total mass=,e11.4,
     $ 17H effective temp.=,E11.4,/,6X,21h photospheric radius=,e11.4,
     $ 13H core radius=,e11.4,11H core lum.=,E11.4)
 2002 format(6h irad=,i2,8h noburn=,i2,8h onemq0=,1pe10.3,
     $ 8H onemq1=,e10.3,7H rlums=,e10.3,7H rlumc=,e10.3)
 2003 format(10x,6h iout=,i2,7h ihmin=,i3,7h ihmax=,i3)
 2004 format(3x,1hi,3x,6hradius,5x,4htemp,5x,6hsp.Vol,3x,8hpressure,
     $ 5x,2hg1,7x,4hg3m1,5x,7hopacity,3x,6hdlkdlr,4x,6hdlkdlt,6x,
     $ 2hcv,7x,4hfrft,4x,8hint.Mass)
 2005 format(1x,i4,1p,6e10.3)
 2006 Format(1x,i4,1p,5e15.8)
C
 9000 format(1x,'Enter ntime, nprt, amplitude, and dtfac. ',$)
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
      data r,a,bk,avagd,ad3 / 8.31434D7,7.56471D-15,8.6170837D-5,
     $     6.02217D23,2.52157D-15 /
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726D-8,7.5595D-5,
     $  5.66961D-5,3.1415926536D0,6.2831853072D0,12.5663706144D0,
     $   25.1327412288D0,4.1887902048D0 /
      data zero,one,two,thre,for,ten,ahf,qrt /
     $      0.d0,1.d0,2.d0,3.d0,4.d0,10.d0,
     $      0.5D0,0.25D0      /
      end
      subroutine linrad(io,ihmin,ihmax,omeout,perout)
      implicit real*8(a-h,o-z)
c
c     *****************************************************************
c      Pulsation equation solver
c      Basic references -- Castor, Ap.J. 166, (109) 1971.
c                          Pesnell, 1983, thesis, Univ. of Florida.
c                          Pesnell, 1987, PASP, submitted.
c     *****************************************************************
c
c          This version of the LNA radial stability analysis is suitable
c       for envelope models generated with the "model" building code.
c       The inner boundary is invalid for a complete model and models
c       of low central condensation. However, the eigenvalue is somewhat
c       insensitive to these details and is reasonable no matter what.
c
      parameter ( accur=1.d-10 )
      parameter ( nmax=200 )
      common/phypar/ r(nmax),tn(nmax),v(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),p(nmax),un(nmax),g1(nmax),g3m1(nmax),
     $               cs(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),press(nmax),up(nmax)
      common/scrtch/ ag1(nmax,3),xo(nmax),yo(nmax),w(nmax),glamp(nmax),
     $       dro(nmax),dp(nmax),dt(nmax),dr(nmax,2),cvt(nmax),dl(nmax)
      common/blk1/   irad,nflag,n
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      np = n+1
      nm1 = n-1
c
      afac = forpi*r(1)**2/sqrt(dm2(1))
      do 60 i=1,np
         afacp = afac
         if( i .eq. np ) goto 40
         afac = forpi*r(i+1)**2/sqrt(dm2(i+1))
c
         dr(i,1) = afacp/v(i)/dm1(i)
         dr(i,2) =-afac/v(i)/dm1(i)
c
c          set up the mechanical matrices.
c
         ag1(i,1) =-afacp*g1(i-1)*p(i-1)*dr(i-1,1)
         ag1(i,2) =-for*g*rm(i)/r(i)**3 + afacp*(g1(i)*p(i)*dr(i,1)-
     $               g1(i-1)*p(i-1)*dr(i-1,2) )
         ag1(i,3) = afacp*g1(i)*p(i)*dr(i,2)
  40    continue
         g1p = g1(n)*p(n)
         ag1(np,1) =-afacp*g1p*dr(n,1)
         ag1(np,2) =-for*g*rm(np)/r(np)**3 - afacp*g1p*dr(n,2)
         ag1(np,3) = zero
  60  continue
c
c          write the matrices to file tape11, if io is less
c       than 3, this write is not done.
c
      if( io .lt. 3 ) goto 76
      write(11,7000)
      do 75 i=1,np
         write(11,7001) i,ag1(i,1),ag1(i,2),ag1(i,3)
  75  continue
  76  continue
c
c          calculate the acoustic tavel time from surface to the
c       innermost zone. This is stored in transt.
C
      transt = zero
      do 77 i=2,np
         transt = transt+(r(i)-r(i-1))/sqrt(p(i-1)*v(i-1)*g1(i-1))
  77  continue
      rhom = totmas/(pi43*r(np)**3)
      rhom_sol = 1.41d0
      termq = sqrt(rhom/rhom_sol)
      tcon = pi43*g*rhom*dfloat(n)**2
      x0np = r(np)*sqrt(dm2(np))
      omsqp = zero
      omsqc = 5.D-4*dfloat(ihmin)*tcon
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
      if(io.ge.1) write(11,7800) iterad,iv,ih,omsq,omsqc,omsqp
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
         do 82 i=1,np
            xo(i) = zero
            yo(i) = zero
  82     continue
         yo(n) =-ag1(n,3)*x0np
         call trisol(ag1,2,n,omsq,yo,xo,nmax)
         if( icount .eq. 1 ) then
c
c          set up initial derivative
c
            qerra = ag1(np,1)*xo(n)+(ag1(np,2)-omsq)*x0np
            qerra = qerra/xo(2)
            omsql = omsq
            omsq  = omsq*(one+1.d-7)
         else
            qerr = ag1(np,1)*xo(n)+(ag1(np,2)-omsq)*x0np
            afacp = dabs((omsq-omsql)/omsql)
            qerr = qerr/xo(2)
            if(io.ge.1) write(11,1400) icount,omsq,omsql,afacp
            if( afacp .lt. 1.d-12 ) goto 86
            omsq2 = (qerra*omsq-qerr*omsql)/(qerra-qerr)
            omsql = omsq
            domsq = omsq2-omsq
            domsq = dsign(dmin1(dabs(domsq),qrt*dabs(omsq)),domsq)
            omsq  = omsq + domsq
            qerra  = qerr
         endif
  81  continue
c
c          no convergence in the adiabatic eigenvalue.
c
      write(6,8100) omsq,iterad,iv,ih
      omsqp = zero
      omsqc = 5.D-6*tcon*dfloat(ih+1)
      goto 79
c
c          converged to omsq value, check if positive and with
c      with proper number of nodes in the displacement eigen-
c      vector. If either is not true, try again.
c
  86  continue
      if( omsq .le. zero ) then
c
c          negative omsq value, stop working on this mode and
c       move on.
c
         write(6,8600) omsq,iv,ih
         omsqp = zero
         omsqc = 5.D-6*tcon*dfloat(ih+1)
         goto 79
      endif
c
c          Positive OMSQ value, find eigenvalues
c
      do 90 i=2,nm1
         yo(i) = zero
  90  continue
      xo(np) = x0np
      yo(n)  =-ag1(n,3)*x0np
      call trisol (ag1,2,n,omsq,yo,xo,nmax)
      iv = nodes(xo,np,io)
      if( iv .eq. ih-1 ) goto 103
c
c      if initial omsgc gives iv .lt. ih-1, increase omsgc and try again
c
      if( omsqs .eq. omsqc .and. iv .lt. ih-1 ) goto 78
      if( iv .lt. ih-1 ) omsqp = omsqs
      if( iv .gt. ih-1 ) omsqc = omsqs
  80  continue
c
c          not converged to a period with right number of nodes,
c       try for next mode.
c
      omega = sqrt(omsq)
      period = twopi/(omega*86400.D0)
      write(11,8000) iv,ih,omega,period
      omsqp = zero
      omsqc = 5.D-6*tcon*dfloat(ih)
      goto 79
c
 103  continue
      omega = sqrt(omsq)
      xo(1) = zero
      xo(n+2) = zero
      stwait = zero
      rke = zero
      do 100 i=1,n
         w(i) = xo(i)/(r(i)*sqrt(dm2(i)))
         dro(i) = dr(i,1)*xo(i)+dr(i,2)*xo(i+1)
         dp(i) = g1(i)*dro(i)
c
c          Set up the weight function calculation.
c
         glamp(i) =-for*(g*rm(i+1)/r(i+1))*(xo(i+1)/r(i+1))**2 +
     $      p(i)*v(i)*g1(i)*dro(i)**2*dm1(i)
         stwait = stwait + glamp(i)
         rke = rke + xo(i+1)*xo(i+1)
 100  continue
      w(1) = zero
      w(np) = one
      period = twopi/omega
c      call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax)
      period = period/86400.D0
      qvalue = period*termq
      stwait = stwait/rke
      qch = dabs((stwait-omsq)/omsq)
      write(6,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch
      write(11,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch
      iv = nodes(w,np,1)
c
c          normalize the weight function per zone to omsq.
C
      do 105 i=1,n
         glamp(i) = glamp(i)/(omsq*rke)
 105  continue
c
      if(io.ge.2) write(11,1000) omega,period,(i,w(i+1),dro(i),
     $ dp(i),cvt(i),glamp(i),i=1,n)
c
c          dump adiabatic eigenvectors to casplt.
c
      call pltdmp(w,nmax,np,'dr/r')
      call pltdmp(dro,nmax,n,'drho')
      call pltdmp(glamp,nmax,n,'wait')
c
      omsqp = omsq
      omsqc = omsq
  79  continue
      omeout = omega
      perout = period*86400.D0
      return
c
 1400 format(1x,i4,1p,7e11.4)
C
 7000 format(1h1,/1x,26hpulsation matrices: ag1(3))
 7001 format(1x,i4,1x,1p,3e10.3)
 7800 Format(8h iterad=,i3,4h iv=,i3,4h ih=,i3,6h omsq=,1pe12.4,
     $   7H omsqc=,e12.4,7H omsqp=,e12.4)
 8000 Format(40h cant find no. Of nodes.eq.Ih-1...Nodes=,i3,4h ih=,i3,
     $ 7h omega=,1pe13.6,10H period = ,e13.6)
 8600 Format(14h adiabatic f2=,1pe17.10,18H is .lt. 0 For iv=,
     $  i3,4h ih=,i3)
 8100 format(14h adiabatic f2=,1pe15.6,20H not converged after,i4,
     $   16h iterations, ih=,i4,17h number of nodes=,i4)
 1000 format(1h1,19h adiabatic solution,/,9h  omega =,1pe17.9,5X,
     $ 8hperiod =,0pf12.6,//,3X,1hi,7x,2hdr,12x,2hdl,11x,4hdrho,11x,
     $  2hdt,12x,2hdp,11x,4hcvtm,10x,6hweight/(1x,i4,1p,5e14.6))
 1001 Format(/,22h linear adiabatic mode,4h iv=,i3,4h ih=,i3,/,
     $ 1x,8h omega =,1pe16.8,9H sec(-1),,1x,8hperiod =,e12.5,
     $ 6H days,,1x,9hq-value =,e10.3,/,
     $ 1X,38h eigenfrequency from matrix solution =,e12.4,/,
     $ 1X,38h eigenfrequency from weight function =,e12.4,
     $ 7H error=,e12.4)
 7900 Format(1h1,//,20x,28h***** begin search for mode ,i2,6h *****)
      end
      subroutine cvtm(n,p,rl,cv,t,dm,x,nmax)
      implicit real*8(a-h,o-z)
      dimension cv(nmax),t(nmax),dm(nmax),x(nmax)
c
c       compute x=sum(cv*t*dm/rl*p) from outside in
c
      sum = 0.d0
      do 10 i=1,n
         j = n+1-i
         sum = sum + cv(j)*t(j)*dm(j)/(rl*p)
         x(j) = sum
  10  continue
      return
      end
      subroutine trisol(a,imin,imax,x,y,z,nmaxin)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=200 )
      dimension a(nmax,3), y(nmax), z(nmax)
      common/blk37/  useno(nmax), e(nmax), f(nmax)
      id = imax+1-imin
      e(imax) =-a(imax,1)/(a(imax,2)-x)
      f(imax) = y(imax)/(a(imax,2)-x)
      do 10 ii=2,id
         I = imax+1-ii
         den = (a(i,2)-x)+a(i,3)*e(i+1)
         if( ii .eq. Id ) goto 10
         e(i) =-a(i,1)/den
  10     f(i) = (y(i)-a(i,3)*f(i+1))/den
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
c finds no. Of nodes(sign changes) in nz elements of array w
c
      dimension w(nz),inode(30)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      nodes = 0
      nz1 = nz-1
      do 10 i=3,nz1
         if( w(i-1)*w(i) .ge. zero ) goto 10
c zero crossing found. Is it general enough to be a node#
         if( w(i-2)*w(i-1) .lt. zero ) goto 10
         if( w(i)*w(i+1) .lt. zero ) goto 10
          nodes = nodes + 1
          inode(nodes) = i-1
  10  continue
      inode(nodes+1) = nz
      if( lout .eq. 0 ) Return
      write(11,1000) nodes
      if( nodes .le. 0 ) Return
      write(11,1001) (inode(i),i=1,nodes)
 1000 format(/,1x,i2,12h nodes found)
 1001 format(7h nodes=,30i4)
      return
      end
      subroutine pltint(nmaxin,n)
      implicit real*8(a-h,o-z)
c
c          initialize the plot file and write the title line.
C
      parameter ( nmax=200 )
      common/phypar/ rn(nmax),tn(nmax),vn(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn(nmax),un(nmax),g1(nmax),g3m1(nmax),
     $               cs(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/blk37/  x(nmax),rholog(nmax),onemq(nmax)
c
      open (unit=12,file='casplt',status='new')
c
      rms = totmas/1.991D33
      rl = dlog10(rlumgv/3.82D33)
      write(12,1000) rms,rl,teff
      qm = ahf*dm1(n)
      do 10 i=1,n
         x(i+1) = rn(i+1)/rn(n+1)
         rholog(i) =-dlog10(vn(i))
         onemq(n+2-i) =-dlog10(qm/totmas)
         qm = qm + dm1(n+1-i)
  10  continue
      onemq(1) =-dlog10(qm/totmas)
      x(1) = rn(1)/rn(n+1)
      call pltdmp(x,nmax,n+1,'x   ')
      call pltdmp(tn,nmax,n,'t   ')
      call pltdmp(pn,nmax,n,'p   ')
      call pltdmp(rholog,nmax,n,'rho ')
      call pltdmp(onemq,nmax,n+1,'1-q ')
      call pltdmp(g1,nmax,n,'g1  ')
      return
 1000 format(1x,f7.2,2HMs,3h L=,f7.3,4HLsun,6h Teff=,f8.1)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
      character*4 ititl
c
c          Write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p,6e12.4) )
      end
      subroutine hydro(iout, dtin, itime, ierrhy)
      implicit real*8(a-h,o-z)
c
      parameter ( onth=1.d0/3.d0 )
      parameter ( nmax=200 )
      common/phypar/ rn(nmax),tn(nmax),vn(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn(nmax),un(nmax),g1(nmax),g3m1(nmax),
     $               cs(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),p(nmax),up(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/blk8/   g,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
c
      dt = dtin
      do 100 itry=1,20
         do 10 i=2,npts
c            un(i) = up(i) - dt*(G*rm(i)/rp(i)**2 + pi4*rp(i)**2/dm2(i)*
c     $          (p(i)-p(i-1)))
            rn(i) = rp(i) + dt*up(i) - ahf*dt*dt*
     $            (G*rm(i)/rp(i)**2 + pi4*rp(i)**2*(p(i)-p(i-1)))/dm2(i)
  10     continue
         pnpts = p(npts) - 2.52157D-15*tp(npts)**4
c         un(npts+1) = up(npts+1) - dt*( g*rm(npts+1)/rp(npts+1)**2 -
c     $                pi4*rp(npts+1)**2*pnpts/dm2(npts+1) )
         rn(npts+1) = rp(npts+1) + dt*up(npts+1) - ahf*dt*dt*
     $                (G*rm(npts+1)/rp(npts+1)**2 -
     $                 pi4*rp(npts+1)**2*pnpts/dm2(npts+1) )
         do 20 i=2,npts
c            rn(i) = rp(i) + dt*un(i)
            un(i) = two*(rn(i)-rp(i))/dt - up(i) 
            if( rn(i) .le. rn(i-1) ) then
               write(6,2000) i, rn(i-1), rn(i), vn(i), un(i)
               is = i
               dt = dt/2.d0
               goto 100
            endif
  20     continue
         un(npts+1) = two*(rn(npts+1)-rp(npts+1))/dt - up(npts+1) 
c         rnp1 = rp(npts+1) + dt*un(npts+1)
         if( rn(npts+1) .le. rn(npts) ) then
            write(6,2000) npts+1, rn(npts), rn(npts+1), vn(npts),
     $                    un(npts+1)
            rn(npts+1) = (vn(npts)*dm1(npts)/pi43 + rn(npts)**3)**onth
         endif
         do 30 i=1,npts
            vn(i) = pi43*(rn(i+1)**3-rn(i)**3)/dm1(i)
            pn(i) = p(i)*(vp(i)/vn(i))**g1(i)
            tn(i) = tp(i)*(vp(i)/vn(i))**g3m1(i)
            cs(i) = sqrt(g1(i)*vn(i)*pn(i))
  30     continue
         dtin = dt
         ierrhy = 1
         return
 100  continue
      write(6,1000) itime,dt,is
      ierrhy =-1
      return
 1000 format(1x,'No go in time step ',i3,8h dtime =,1pe11.3,
     $          ' zone: ',i4)
 2000 format(1x,i4,1p6e12.5)
c 5000 format(1x,i4,1p,2e12.4)
      end
      function qtime(r,cs,nmax,n,iprtin)
      implicit real*8(a-h,o-z)
      logical iprt,iprtin
c
c          Returns the minimum time step allowed using an
c       acoustic propagation velocity across each zone.
c
      dimension r(nmax),cs(nmax)
      iprt = iprtin
      reg = (r(2)-r(1))/cs(1)
      is = 1
      do 10 i=2,n-1
         if( (r(i+1)-r(i))/cs(i) .lt. reg ) then
            reg = (r(i+1)-r(i))/cs(i)
            is = i
         endif
  10  continue
      qtime = reg
c      write(6,1000) is,reg
      if( iprt ) write(11,1000) is,reg,r(is),r(is-1),cs(is)
      return
 1000 format(21h minimum time in zone,i4,11h time step=,1pe10.3,
     $   3X,1p,3e11.4)
      end
