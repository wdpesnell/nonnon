From:	ESSDP2::CCVAX2::090245       "David Kratzer, ICN Consulting Office" 31-AUG-1990 13:14:42.75
To:	ESSDP2::DE614::PESNELL
CC:	
Subj:	NONAD.FOR;1

*byp3
*rcft i=%me, x=%mex
*destroy b%me m%me
      program nonad( tty, tape5=tty, tape6=tty)
      implicit real(a-h,o-z)
      character*80 ititl
      character*1 iyorn
      logical bigout
c
c          the lagrangian non-radial adiabatic eigenvalue analysis
c       for stellar pulsations. the notation is lifted from the model
c       code and radial stability analysis. written summer 1986.
c          this version of the code treats the "realistic" stellar
c       models with the r**l transformation in the eigenvectors.
c
c                                                12/20/83 wd pesnell
c
c          Adapted to do nonadiabatic pulsations, with central 
c       ball included, July 9-20, 1989, Los Alamos Nat'l. Lab.
c
c                                                 7/30/89 WD Pesnell
c
      parameter ( nmax=600 )
      common/phypar/ rp(nmax),tp(nmax),vp(nmax),cv(nmax),dkdr(nmax),
     $   dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/stnmod/ nmode
c
c          read in the input variables from file 10.
c
c      open (unit=1,file='sumout',status='unknown',acess='append')
c      open (unit=5,file='sys$input',status='unknown')
c      open (unit=6,file='sys$output',status='unknown')
      open (unit=10,file='lnrad',status='old')
cPC      open (unit=10,file=' ',status='old')
c
      open (unit=1,file='sumout.adn',status='unknown')
      open (unit=11,file='adnout',status='unknown')
c
      call filerep
c
      read(10,1001,end=100) ititl
c
c          read in the input variables from file 10.
c
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
c          Does user want to write the initial model?
c
      bigout = .false.
      write(6,1100)
      read(5,1101,end=110) iyorn
      if( iyorn .eq. 'y' .or. iyorn .eq. 'Y' ) bigout = .true.
 110  continue
c
c          output the initial model data
c
      write(11,2000) ititl,npts
      write(11,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(11,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
      write(11,2004)
      write(11,2005) npts,rp(1),tc,1.e0/rhoc,pc,g10,g3m10,opac0,
     $      dkdr0,dkdt0,cv0,frft(1),rm(1)
      if( bigout ) then
         write(11,2005) (i,rp(i+1),tp(i),vp(i),p(i),g1(i),g3m1(i),
     $      akap(i),dkdr(i),dkdt(i),cv(i),frft(i+1),rm(i+1),i=1,npts)
c
c          check for nuclear burning and print out those
c       zones participating.
c
         if( tp(1) .ge. 3.d6 ) then
            do 10 i=1,npts
               imax = np1-i
               if( sorce(imax) .ne. zero ) goto 15
  10        continue
            write(11,2100)
            goto 20
  15        continue
            do 16 i=1,imax
               imin = i
               if( sorce(i) .ne. zero ) goto 17
  16        continue
            write(11,2100)
            goto 20
  17        continue
            write(11,2101)
            write(11,2102) (i,sorce(i),dtsorc(i),dvsorc(i),i=imin,imax)
  20        continue
         endif
      endif
      call pltint(npts,ititl)
      write(1,2000) ititl,npts
      write(1,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(1,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
c
c          The pulsation control loop.
c
c          iout:= 0 for little output (no eigenvectors)
c                 1 for convergence information
c                 2 for eigenvectors, matrices and conver. info.
c
c          lin:= laditudinal quantum number.
c
c          nomega:= number of points in frequency scan. (lin>0)
c                   number of nodes -1 in first mode (lin=0)
c
c          ihmax:= not used (lin>0), enter a 0.
c                   number of nodes - 1 in last mode (lin=0)
c
c          permin:= starting period for scan (sec).
c
c          permax:= ending period for scan (sec).
c
c          (it is not necessary for permin to be less than permax.)
c
c          Example input line: 0 0 1 3 0 0
c
c          This would calculate the first three radial modes (fundamental
c       and first and second overtone, with little output.
c
c          Example input line: 3 2 10 0 4000 4500
c
c          This would search 10 periods between 4000 and 4500 seconds
c       for l=2 modes, producing a great deal of output.
c
c          In every case, a plot file is produced with the eigenvectors,
c       weight functions, and global quantities such as pressure.
c
 200  continue
         write(6,1003)
         read(5,*,end=900) iout,lin,nomega,ihmax,permin,permax
         if( lin .lt. 0 ) then
            goto 900
         elseif( lin .eq. 0 ) then
            write(6,3000)
            call cmplna(iout,nomega,ihmax,npts)
         else
            omlow = (twopi/permax)**2
            omhigh = (twopi/permin)**2
            write(11,2003) iout,lin,nomega,omlow,omhigh
            lval = lin
            nmode = 0
            call lnanon(npts,iout,lval,nomega,omlow,omhigh)
         endif
      goto 200
c
c          if end of file, try again.
c
 100  continue
      write(1,1002)
      stop
 900  continue
      stop
c
 1000 format(1p,4e20.13)
 1001 format(a80)
 1002 format(1x,'End-of-File on title, try again.')
 1003 format(1x,27henter pulsation parameters:,/,
     $   1x,38h io, l, nom, ihmax, permin, and permx.)
 1100 format(1x,'Write initial model data? (y/n) '  )
 1101 format(a)
c
 2000 format(50h1 initial model in hydrostatic and thermal balance,//,
     $      a80,1x,i5,6h zones,/)
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
 3000 format(1x,'No radial analysis available at this time.')
      end
      block data
      implicit real(a-h,o-z)
c
c       fundamental constants are from novotny, introduction to
c    stellar atmospheres and interiors, 1973, appendix ii.
c                                                 2/16/83  wdp
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/thermo/ r,a,bk,avagd,ad3
      data r,a,bk,avagd,ad3 / 8.31434e7,7.56471e-15,8.6170837e-5,
     $     6.02217d23,2.52157e-15 /
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726e-8,7.5595e-5,
     $  5.66961e-5,3.1415926536e0,6.2831853072e0,12.5663706144e0,
     $   25.1327412288e0,4.1887902048e0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.e0,1.e0,2.e0,3.e0,4.e0,10.e0,
     $      0.5e0,0.25e0      /
      end
      subroutine pltint(n,ititl)
      implicit real(a-h,o-z)
      character*80 ititl
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $   dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),omq(nmax),tlog(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      dimension plog(nmax),xrad(nmax)
c
      open (unit=12,file='nonplt',status='unknown')
      write(12,1000) ititl
      np1 = n + 1
      call pltdmp(r,nmax,n+1,'r   ')
      sfact = one - 5.d-14
      do 10 i=1,n
         omq(i) =-alog10(one - sfact*rm(i)/totmas)
         tlog(i) = alog10(t(i))
         plog(i) = alog10(p(i))
         xrad(i) = r(i)/r(np1)
  10  continue
      xrad(np1) = one
      omq(n+1) =-alog10(one-sfact)
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
      do 30 i=1,n
         xrad(i) = 0.
         if( bv(i) .ne. 0. ) xrad(i) = alog10(abs(bv(i)))
  30  continue
      call pltdmp(xrad,nmax,n,'bv  ')
      return
 1000 format(a80)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real(a-h,o-z)
      character*4 ititl
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p,6e12.4) )
      end
      subroutine lnanon(nin,iout,lin,nomgin,omlow,omhigh)
      implicit real(a-h,o-z)
      logical iout0, iout1
c
c          Stability analysis of linear, non-adiabatic nonradial
c      stellar pulsations. First version has the non-adiabatic
c      effects included in the momentum equation but the entropy
c      change will be zero everywhere for now. The notation is
c      based on the lagrangian variations of Castor Ap.J. 162(166)1971,
c      and the actual linearization is based on the technique of
c      Lynden-Bell and Ostriker.
c
c                                                9/15/83 WD Pesnell
c
      parameter ( accur=1.e-12 )
      parameter ( onth=0.333333333333333e0 )
      parameter ( nmax=600 )
      parameter ( nmax3=nmax*3 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $        dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
     $               bv(nmax)
      common/blkdrs/ drz(nmax),drint(nmax),gor(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk37/  drdm(nmax,2),dl1(nmax,2),dl2(nmax,2)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),            ap3(nmax,3),ap4(nmax)
      common/linear/ dr(nmax),dh(nmax),gam(nmax),y(nmax),z(nmax),
     $               dp(nmax),adrho(nmax),weight(nmax),cvt(nmax),
     $               wtherm(nmax),wgrav(nmax),wdiag(nmax),wcross(nmax),
     $               rkepl(nmax),stwait(nmax),wint(nmax),
     $               spac(nmax,11)
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      dimension rkl2(nmax)
      dimension x(nmax3),stomeg(50),stomsq(50)
c
      n = nin
      np1 = n + 1
      nomega = nomgin
      lval = lin
      rl = float(lval)
      rl1 = rl*(rl+one)
      iout0 = .false.
      if( iout .gt. 0 ) iout0 = .true.
      iout1 = .false.
      if( iout .gt. 1 ) iout1 = .true.
      rhobar = totmas/(r(np1)**3*pi43)
c
c          The routine matset puts the matrices into order.
c
      call matset(n,iout,lval,rkl2,z)
      rz0 = ahf**onth*r(1)
      rkl20 = rl1/rz0**2
      pi4g = forpi*g
      drdm01 =-ahf*(rl+one)/rz0**2
      drdm02 =-ahf*(rl+one)/rz0**2
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
      dom = (omhigh - omlow)/float(nomega - 1)
      do 30 i=1,nomega
         omsq = omlow + float(i-1)*dom
         fomeg = fomsq(omsq,x,xnorm,n)
         if( iout0 ) write(11,4400) i,omsq,fomeg
         if( i .gt. 1 ) then
            if( fomeg*fomi .lt. zero ) then
               iroot = iroot + 1
               stomeg(iroot) = omsq - dom*fomeg/(fomeg-fomi)
            endif
         endif
         fomi = fomeg
  30  continue
      if( iroot .eq. 0 ) then
c
c          no roots in this interval, return to main program
c
         write(1,9000) lval,omlow,omhigh
         write(6,9000) lval,omlow,omhigh
         write(11,9000) lval,omlow,omhigh
         return
      else
         write(1,3000) lval,(stomeg(i),i=1,iroot)
         write(6,3000) lval,(stomeg(i),i=1,iroot)
         write(11,3000) lval,(stomeg(i),i=1,iroot)
      endif
c
c          attempt to converge on each of the frequencies found in
c       discriminant search. a secant method to find the zeroes of
c       error function is used.
c
      nomsq = 0
      do 100 nroot=1,iroot
         omsq = stomeg(nroot)
         do 45 icount=1,30
            err = fomsq(omsq,x,xnorm,n)
            if( icount .eq. 1 ) then
               omsq1 = omsq
               erra  = err
               omsq  = omsq*(one + 1.e-7)
            else
               afac = abs((omsq-omsq1)/omsq)
               if( iout0 ) write(11,4400) icount,omsq,afac,err,erra
               if( afac .le. accur ) goto 47
               omsq2 = (erra*omsq-err*omsq1)/(erra-err)
               omsq1 = omsq
               erra  = err
               omsq  = omsq2
            endif
  45     continue
         write(11,4500) nroot,omsq
         goto 100
c
c          converged to a root, calculate the eigenvectors.
c
  47     continue
         x(3*(n+1)) = xnorm
c         call invint(omsq,x,finerr,xnorm,n)
c
c          do not repeat a root.
c
         nomsq = nomsq + 1
         if( nomsq .gt. 1 ) then
            do 491 ido=1,nomsq-1
               if( abs((omsq-stomsq(ido))/omsq) .lt. accur ) then
                  nomsq = nomsq - 1
                  goto 100
               endif
 491        continue
         endif
         stomsq(nomsq) = omsq
         do 40 i=1,np1
            dr(i) = (r(i)/r(np1))**lval*(x(3*i)/r(i)**2)
            if( i .eq. 1 ) then
               dh(i) = (x(3*i-1)/rz0**2)*(rz0/r(np1))**lval
               adrho(i) = ((drdm01+drdm02)*x(3*i) +
     $             rl1*x(3*i-1)/rz0**2 )*(rz0/r(np1))**lval
               dp(i) = g10*adrho(i)
               gam(i) = x(3*i-2)*(rz0/r(np1))**lval
            else
               dh(i)=(x(3*i-1)/rzone(i-1)**2)*(rzone(i-1)/r(np1))**lval
               adrho(i) = (drdm(i-1,1)*x(3*i-3)+drdm(i-1,2)*x(3*i) +
     $           rl1*x(3*i-1)/rzone(i-1)**2)*(rzone(i-1)/r(np1))**lval
               dp(i) = g1(i-1)*adrho(i)
               gam(i) = x(3*i-2)*(rzone(i-1)/r(np1))**lval
            endif
  40     continue
         dlnr = ahf*(rl+one)*drz(np1)/r(np1)
         gam(n+2) = gam(n)*(one-dlnr)/(one+dlnr)
         dh(n+2) = zero
         qdr0 = dr(2)*(r(1)/r(2))**(lval-2)
         idr = nodes(dr, n+1, 1)
         idh = nodes(dh, n+1, 1)
c
c          weight function for the adiabatic oscillation.
c
         wtherm(1) = (g10*pc/rhoc)*cormas*adrho(1)**2
         wgrav(1) = dr(1)*dm2(1)*rl*gam(1)
         wcross(1) = two*rl1*gor(1)*r(1)*dr(1)*dh(1)*dm2(1)
         wdiag(1) = (pi4g*rho(1)-gor(1)*for/r(1) )*
     $                   (dr(1)*r(1))**2*dm2(1)
         weight(1) = wtherm(1) + wgrav(1) + wcross(1) + wdiag(1)
         rgwait = weight(1)
         wint(1) = weight(1)
         rsum = (dr(1)*r(1))**2*dm2(1)
         hsum = (dh(1)*rz0)**2*cormas
         rkepl(1) = (dr(1)*r(1))**2*dm2(1) +
     $                   (dh(1)*rz0)**2*cormas
         rbar = r(1)*(dr(1)*r(1))**2*dm2(1)
         hbar = rz0*(dh(1)*rz0)**2*cormas
         crsum = zero
         do 48 i=1,n
c
c          weight function for the adiabatic oscillation.
c
            wtherm(i+1) = g1(i)*p(i)*v(i)*dm1(i)*adrho(i+1)**2
            wgrav(i+1) = rl1*dh(i+1)*dm1(i)*gam(i+1) +
     $      (dr(i+1)*r(i+1))*dm2(i+1)*(gam(i+2)-gam(i+1))/drz(i+1)
            wcross(i+1) = two*gor(i+1)*r(i+1)*dr(i+1)*dm2(i+1)*
     $            rl1*ahf*(dh(i+1)+dh(i+2))
            wdiag(i+1) = ( pi4g*ahf*(rho(i)+rho(i-1)) -
     $            for*gor(i+1)/r(i+1) )*(dr(i+1)*r(i+1))**2*dm2(i+1)
            rsum = rsum + (dr(i+1)*r(i+1))**2*dm2(i+1)
            hsum = hsum + (dh(i+1)*rzone(i))**2*dm1(i)
            rkepl(i+1) = (dr(i+1)*r(i+1))**2*dm2(i+1) +
     $                   (dh(i+1)*rzone(i))**2*dm1(i)
            rbar = rbar + r(i+1)*(dr(i+1)*r(i+1))**2*dm2(i+1)
            hbar = hbar + rzone(i)*(dh(i+1)*rzone(i))**2*dm1(i)
            weight(i+1) = wtherm(i+1)+wgrav(i+1)+wcross(i+1)+wdiag(i+1)
            rgwait = rgwait + weight(i+1)
            wint(i+1) = wint(i) + weight(i+1)
  48     continue
c
c          normalize the weight function.
c
         rke = rsum + rl1*hsum
         if( rke .le. zero ) rke = one
         rbar = (rbar + rl1*hbar)/(rke*r(n+1))
         crsum = zero
         dh(n+2) = zero
         do 46 i=1,np1
            crsum = crsum + dr(i)*r(i)**2*dm2(i)*(dh(i)+dh(i+1))
            rkepl(i) = rkepl(i)/rke
            wint(i) = wint(i)/(rke*omsq)
            weight(i) = weight(i)/(rke*omsq)
            wtherm(i) = wtherm(i)/(rke*omsq)
            wgrav(i) = wgrav(i)/(rke*omsq)
            wdiag(i) = wdiag(i)/(rke*omsq)
            wcross(i) = wcross(i)/(rke*omsq)
  46     continue
         rgwait = (rgwait/rke)
         qch = abs((omsq-rgwait)/omsq)
         if( rgwait*omsq .le. 0. ) then
            qch = abs(abs(omsq)-abs(rgwait))/abs(omsq)
         endif
c
c          evaluate the rotational spitting coefficients.
c
         crsum = rl1*crsum/rke
         hrsum = hsum/rke
         ckl = crsum + hrsum
         index = knode(dr,dh,nmax,np1,np,ng)
         omega = sqrt(omsq)
         period = twopi/omega
         call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax,itrans)
         pdays = period/86400.e0
         write(1,5505) nroot,omega,period,pdays,omsq,rgwait,qch
         write(1,5503) index,np,ng,itrans,rbar
         write(1,5506) ckl,crsum,hrsum
         write(6,5505) nroot,omega,period,pdays,omsq,rgwait,qch
         write(6,5503) index,np,ng,itrans,rbar
         write(11,5500) nroot,omega,period,pdays,omsq,rgwait,qch
         write(11,5503) index,np,ng,itrans,rbar
         write(11,5506) ckl,crsum,hrsum
c
         if( iout1 ) then
            write(11,5502)
            write(11,5501) (i,dr(i),dh(i),adrho(i),gam(i),
     $             dp(i),weight(i),cvt(i),i=1,np1)
         endif
         call pltdmp(dr,nmax,np1,'dr/r')
         call pltdmp(dh,nmax,np1,'dh/h')
         call pltdmp(adrho,nmax,np1,'drho')
         call pltdmp(rkepl,nmax,np1,'ke  ')
         call pltdmp(weight,nmax,np1,'wait')
         call pltdmp(wint  ,nmax,np1,'wint')
         call pltdmp(wtherm,nmax,np1,'wthr')
         call pltdmp(wgrav,nmax,np1,'wgrv')
         call pltdmp(wdiag,nmax,np1,'wdia')
         call pltdmp(wcross,nmax,np1,'wcrs')
c
         call orthog(dr,dh,n,omsq,index,1)
         call eigenf(omsq,n,lval,rke,gor)
c
c          non-adiabatic loop removed 10/30/89.
c
 100  continue
      call orthog(dr,dh,n,omsq,index,-1)
      return
c
 2200 format(1h1,//,20h ag1,ag2,ag3,ag4,ah4)
 2201 format(1x,i5,1p,10e12.4)
 2300 format(1h1,//,24h ah1,ah2,ah3,ap1,ap3,ap4)
 2600 format(1h1,//,10x,36h the pulsation matrices, ak1,ak2,ak4,/)
c
 3000 format(1x,23h roots for analysis: l=,i3,/,(3x,1p,5e15.6) )
 3001 format(1h1,1x,26hBegin discriminant search.)
c
 4400 format(1x,i5,1p,4e15.6)
 4500 format(1x,27hno convergence root number ,i4,
     $       11h frequency=,1p,2e15.6)
 5500 format(1h1,//,5x,4hroot,i5,/,8h omega =,1pe13.6,9h period =,e13.6,
     $   7h (secs),e15.6,7h (days),
     $   /,2x,33heigenvalue from matrix solution =,e12.4,
     $   /,2x,33heigenvalue from weight function =,e12.4,
     $   8h error =,e10.3)
 5501 format(1x,i5,1p,7e15.6)
 5502 format(/4x,1hi,8x,4hdr/r,10x,4hdh/h,10x,8hdrho/rho,8x,5hgamma,
     $   11x,4hdp/p,10x,6hweight,10x,4hcvtm/)
 5503 format(1x,24h radial quantum number =,i4,
     $   26h number of nodes, p-type =,i4,9h g-type =,i4,/,
     $   1x,18h transition zone =,i4,/,
     $   1x,28h effective radius of mode = ,1pe11.4)
 5505 format(/,5x,4hroot,i5,/,8h omega =,1pe13.6,9h period =,e13.6,
     $   7h (secs),e15.6,7h (days),
     $   /,2x,33heigenvalue from matrix solution =,e12.4,
     $   /,2x,33heigenvalue from weight function =,e12.4,
     $   8h error =,e10.3)
 5506 format(1x,31h rotational splitting, total = ,1pe11.4,
     $   5h a = ,e11.4,5h b = ,e11.4)
c
 9000 format(1x,44h from lnanon...no adiabatic roots found, l =,i3,/,
     $   1x,17h frequency range:,1p,2e15.6)
      end
      function knode(x,y,nmax,nin,nps,ngs)
      implicit real(a-h,o-z)
c
c          returns the number of zero-crossings in the counter-
c       clockwise direction as np and the clockwise crossings
c       in ng.
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax),y(nmax)
      n = nin
      np = 0
      ng = 0
      knode = 0
      if( n .eq. 1 ) return
      nm1 = n - 1
      do 10 i=1,nm1
         if( x(i)*x(i+1) .lt. zero ) then
            if( y(i)*y(i+1) .lt. zero ) write(11,1000)
            if( x(i+1) .gt. x(i) ) then
               if( y(i) .ge. zero ) ng = ng + 1
               if( y(i) .lt. zero ) np = np + 1
            else
               if( y(i+1) .ge. zero ) np = np + 1
               if( y(i+1) .lt. zero ) ng = ng + 1
            endif
         endif
  10  continue
      if( x(n-1)*x(n) .lt. zero ) write(11,1001)
      knode = np - ng
      nps = np
      ngs = ng
      return
 1000 format(1x,50hFrom knode...quadrant jump default rotations used.)
 1001 format(1x,51hFrom knode...node at outer boundary is not counted.)
      end
      subroutine matset(nin,iout,lin,rkl2,z)
      implicit real(a-h,o-z)
c
c          Stability analysis of linear, non-adiabatic nonradial
c      stellar pulsations. First version has the non-adiabatic
c      effects included in the momentum equation but the entropy
c      change will be zero everywhere for now. The notation is
c      based on the lagrangian variations of Castor Ap.J. 162(166)1971,
c      and the actual linearization is based on the technique of
c      Lynden-Bell and Ostriker.
c
c                                                9/15/83 WD Pesnell
c
c          Modified to full nonadiabatic calculations summer 1989.
c       Uses the r**l transformation on all variables.
c
c                                                8/23/89 WD Pesnell
c
      parameter( onth=0.333333333333333e0 )
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $        dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
     $               bv(nmax)
      common/blkdrs/ drz(nmax),drint(nmax),gor(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk37/  drdm(nmax,2),dl1(nmax,2),dl2(nmax,2)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),           ap4(nmax)
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/thermo/ rther,a,bk,avagd,ad3
      common ak1(nmax,4),ak2(nmax,3),ak4(nmax,3)
      dimension rkl2(nmax),z(nmax)
c
      n = nin
      np1 = n + 1
      n4 = nin
      lval = lin
      rl = float(lval)
      rl1 = rl*(rl+one)
      rhobar = totmas/(r(np1)**3*pi43)
c
c          define useful quantities.
c
      rz0 = (ahf**onth)*r(1)
      rkl20 = rl1/rz0**2
      dkdr0 = dkdr0 + g3m10*dkdt0
      do 10 i=1,n
         rho(i) = one/v(i)
         dkdr(i) = dkdr(i) + g3m1(i)*dkdt(i)
         rzone(i) = (r(i)**3 + ahf*dm1(i)*v(i)/pi43)**onth
         rkl2(i) = rl1/rzone(i)**2
         if( i .eq. 1 ) then
            drz(i) = rzone(i)-rz0
         else
            drz(i) = rzone(i) - rzone(i-1)
         endif
         drint(i) = r(i+1)-r(i)
         gor(i) = g*rm(i)/r(i)**2
         drdm(i,1) = forpi*rzone(i)/(dm1(i)*v(i)) -
     $                  ahf*(rl+one)/rzone(i)**2
         drdm(i,2) =-(forpi*rzone(i)/(dm1(i)*v(i)) +
     $                  ahf*(rl+one)/rzone(i)**2)
         z(i) = rzone(i)/r(np1)
  10  continue
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      call pltdmp(z,nmax,n4,'rz  ')
      gor(np1) = g*rm(np1)/r(np1)**2
      drz(np1) = p(n)*v(n)/gor(np1)
      pi4g = forpi*g
      rlumx = corlum
c
c          initialize the matrices.
c
c          all matrices are prefixed with the letter a, the second
c       letter denotes the equation ( g=radial component of the
c       momentum equation, h=horizontal component of the momentum
c       equation, p=poisson equation and k= thermal equation)
c       in which each matrix is used. the number tells which
c       variable is to be acted on the the matrix: 1= radial
c       component of the motion eigenvector (mechanical
c       eigenvector), 2= entropy variations (to within a factor
c       of the temperature, actuallly the thermal eigenvector),
c       3= gravitation variations multiplied by the square root
c       of the interface radius differences), 4= horizontal
c       component of the mechanical eigenvector. it should be
c       noted that there is no ap2 matrix nor is ak3 present.
c
c                                              11/28/83 wd pesnell
c
      do 20 i=1,np1
      if( i .eq. 1 ) then
c
c          inner boundary conditions.
c
         afac = forpi*r(1)**3/dm2(1)
         vl = g1(1)*p(1)*(ahf*rl*v(1) + afac)
         drdm01 =-ahf*(rl+one)/rz0**2
c         drdm01 = zero
         drdm02 =-ahf*(rl+one)/rz0**2
c         drdm02 =-(rl+one)/rz0**2
         ag11 = drdm01*g10*pc*(ahf*rl/rhoc - afac)
         ag1(1,1) = zero
         ag1(1,2) =-gor(1)/r(1) + drdm(1,1)*vl +
     $             drdm02*g10*pc*(ahf*rl/rhoc - afac)
         ag1(1,3) = drdm(1,2)*vl + ag11
         ag2(1,1) = g3m10*(ahf*rl - rhoc*afac)
         ag2(1,2) = g3m1(1)*(ahf*rl + rhoc*afac)
         ag3(i,1) = ahf*rl - r(i)/drz(i)
         ag3(i,2) = ahf*rl + r(i)/drz(i)
         ag4(1,1) = zero
         ag41 = g10*pc*(ahf*rl/rhoc - afac)*rl1/rz0**2
         ag4(i,1) = ag41 + ahf*rl1*gor(1)/r(1)
c         ag4(1,2) = rkl2(1)*vl + ag41 + rl1*gor(1)/r(1)
         ag4(i,2) = rkl2(i)*vl + ahf*rl1*gor(i)/r(i)
c
c          horizontal momentum equation.
c
c         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
c     $         ahf*gor(i-1)/r(i-1)
         ah1(i,1) = zero
         ah1(i,2) = (drdm01+drdm02)*g10*pc/rhoc + gor(i)/r(i)
         ah2(i) = g3m10
         ah3(i) = one
         ah4(i) = rkl20*g10*pc/rhoc
c         ap1(1,1) =-pi4g*rho(1)*rzone(1)**2*drdm(1,1)
         ap1(1,1) = zero
c         ap1(1,2) =-pi4g*rho(1)*(rzone(1)**2*drdm(1,2)-
c     $     ahf*alog(rho(2)/rho(1))/alog(rzone(2)/rzone(1)) )
         ap1(1,2) =-pi4g*rhoc*(rz0**2*drdm02-
     $     alog(rho(1)/rhoc)/alog(rzone(1)/rz0) )
         ap311 = r(1)**2/(drint(1)*drz(1))
         ap313 = r(2)**2/(drint(1)*drz(2))
         ap3(1,1) = zero
         ap3(1,2) =-ap313 - ap311*(drz(1)/drz(2))**2
         ap3(1,3) = ap313 + ap311*(drz(1)/drz(2))**2
         ap4(1)   =-pi4g*rl1*rhoc
c
c          the thermal equation matrices.
c
         dtl1 = tc**4
         dtl2 = t(i)**4
         t4oki = dtl2/gkap(i)
         t4oki1 = dtl1/opac0
         diff = t4oki-t4oki1
         dl1(i,1) = zero
         dl1(i,2) = (for*g3m10-dkdr0)*(t4oki1/diff)
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rhoc)))
         dl2(i,1) = zero
         dl2(i,2) = ((for-dkdt0)*(t4oki1/diff))/(cv0*tc)
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rhoc)))
c
c          keep track of the luminosity as a function of mass
c
            glom = zero
            rlumx1 = sorc0*cormas
            if( noburn .ne. 0 ) then
               ipass = i
               rlumx1 = cool(ipass,cormas/totmas)
            endif
            glom1 = frft(i)*(rlumx1/cormas)
     $         *(one+ahf*rl*cormas/(rhoc*forpi*rz0**3))
            dldm = rlumx1/cormas
            aconlm = (one-frft(i))*rlumx1
c            rlumx = rlumx1
c          convert the nuclear energy derivatives.
c
         dedr = (-dedv0+g3m10*dedt0)*sorc0
         deds = dedt0*sorc0/(cv0*tc)
           ak1(i,1) = zero
           ak1(i,2) = zero
           ak1(i,3) = - glom1*for/r(i)**2 -
     $        glom1*(dl1(i,1)*drdm02 + dl1(i,2)*drdm(i,1) )
     $        + dedr*drdm02
           ak1(i,4) =-glom1*dl1(i,2)*drdm(i,2)
           ak4(i,1) = zero
           ak4(i,2) =-rkl20*(glom1*dl1(i,1)-dedr)
           ak4(i,3) =-rkl2(i)*glom1*dl1(i,2)
           ak2(i,1) = zero
           ak2(i,2) =-glom1*dl2(i,1) + deds
           ak2(i,3) =-glom1*dl2(i,2)
c
c          add horizontal heat flow terms.
c
         rlom = zero
         rlom1 = frft(i)*rlumx1
         horiz = rkl20*for*ac3*t4oki1/rhoc**2
c         ak1(i,2) = ak1(i,2) - horiz*g3m10*drdm01 -
c     $       rkl20*ahf*rlom/(rhoc*forpi*rz0**3)
         ak1(i,3) = ak1(i,3) - horiz*g3m10*drdm02 -
     $       rkl20*ahf*rlom1/(rhoc*forpi*r(i)**3)
     $       - horiz*g3m10*drdm01 -
     $       rkl20*ahf*rlom/(rhoc*forpi*rz0**3)
         ak2(i,2) = ak2(i,2) - horiz/(cv0*tc)
         ak4(i,2) = ak4(i,2) + rkl20*(dldm-horiz*g3m10 -
     $      aconlm/(rhoc*forpi*rz0**3) )
      elseif( i .eq. np1 ) then
c
c          outer boundary conditions
c
         tau = gkap(n)*dm1(n)/(for*forpi*r(np1)**2)
         tau = 1.5e0*tau/(one+1.5e0*tau)
         fac = (r(np1)-rphoto)/r(np1)
         ef1 = (one+ahf*fac)/(one+fac)
         ef3 = two/(one+fac)
         dl1(np1,1) = ef1*(for*g3m1(n)-tau*dkdr(n))
c     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*rho(i-1)))
         dl1(np1,2) = dl1(np1,1)
         dl1(np1,2) = ef1*(for*g3m1(n)-tau*dkdr(n))
c     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*rho(i-1)))
         dl2(np1,1) = ef1*(for-tau*dkdt(n))/(t(n)*cv(n))
c     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*rho(i-1)))
         dl2(np1,2) = zero
         glom = rlumgv/dm1(n)
     $      *(one-ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
         glom1 = rlumgv/dm1(n)
     $      *(one+ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
         ak1(np1,1) = glom*dl1(n,1)*drdm(n-1,1)
         ak1(np1,2) =-glom1*dl1(np1,1)*drdm(n,1)+glom*for/r(n)**2 +
     $     glom*(dl1(n,1)*drdm(n-1,2)-dl1(n,2)*drdm(n,1) )
         ak1(np1,3) = glom*dl1(n,2)*drdm(n,2) -
     $      glom1*(dl1(np1,1)*drdm(n,2) + (ef3+tau*ef1)/r(np1)**2 )
         ak1(np1,4) = zero
         ak2(np1,1) = glom*dl2(n,1)
         ak2(np1,2) = glom*dl2(n,2) - glom1*dl2(np1,1)
         ak2(np1,3) = zero
         ak4(np1,1) = glom*rkl2(n-1)*dl1(n,1)
         ak4(np1,2) =-rkl2(n)*(glom1*dl1(np1,1) - glom*dl1(n,2))
         ak4(np1,3) = zero
c
c          Add horizontal heat flow terms. Note the absence of the
c       convective luminosity term. The outer boundary zone is 
c       assumed to be radiative, with no incoming luminosity.
c
         rlumx = rlumgv
         horiz = rkl2(n)*for*ac3*t(n)**4*v(n)**2/gkap(n)
         ak1(np1,2) = ak1(np1,2) - horiz*g3m1(n)*drdm(n,1) -
     $          rkl2(n)*ahf*frft(n)*rlumx*v(n)/(forpi*r(n)**3)
         ak1(np1,3) = ak1(np1,3) - horiz*g3m1(n)*drdm(n,2) -
     $          rkl2(n)*ahf*frft(n+1)*rlumx*v(n)/(forpi*r(n+1)**3)
         ak2(np1,2) = ak2(np1,2) - horiz/(cv(n)*t(n))
         ak4(np1,2) = ak4(np1,2) - rkl2(n)*horiz*g3m1(n)
c
c          outer boundary conditions.
c
         afac = forpi*r(np1)**3/dm2(np1)
         g1p = g1(n)*p(n)-for*ad3*t(n)**4*g3m1(n)
         g3m1ro = g3m1(n)/v(n)-for*ad3*t(n)**3/cv(n)
         ag1(np1,1) = drdm(n,1)*g1p*(ahf*rl*v(i-1) - afac)
         ag1(np1,2) =-for*gor(np1)/r(np1) + pi4g*rho(n) +
     $       drdm(n,2)*g1p*(ahf*rl*v(i-1) - afac)
         ag1(np1,3) = zero
         ag2(np1,1) = g3m1ro*(ahf*rl*v(i-1) - afac)
         ag2(np1,2) = zero
c         ag3(np1,1) =-(rl+one)/(drz(np1)*drint(n) )
         dlnr = (rl+ahf)*drz(np1)/r(np1)
         ag3(np1,1) =-(rl+one)/(one + dlnr )
         ag3(np1,2) = zero
         ag4(np1,1) = rl1*gor(np1)/r(np1) +
     $                rkl2(n)*g1p*(ahf*rl*v(i-1) - afac)
         ag4(np1,2) = zero
c
c          The Poisson equation boundary condition is very messy.
c
         ddr11 =-drz(np1)/(drz(n)*(drz(n)+drz(np1)))
         ddr12 = (drz(np1)-drz(n))/(drz(n)*drz(np1))
         ddr13 = drz(n)/(drz(np1)*(drz(n)+drz(np1)))
         ap3n1 = r(n)**2/(drint(n)*drz(n))
         ap3n3 = r(np1)**2/(drint(n)*drz(np1))
c
         ap1(np1,1) =-pi4g*rho(n)*(rzone(n)**2*drdm(n,1)-
     $      ahf*alog(rho(n)/rho(n-1))/alog(rzone(n)/rzone(n-1)) )
         f = (drz(np1)/r(np1))/(one+dlnr)
         ap1(np1,2) =-pi4g*rho(n)*( (ddr13*two*rl*rzone(n) + ap3n3)*f +
     $      ahf*rho(n)*gor(n+1)*r(n+1)/p(n) + rzone(n)**2*drdm(n,2) )
c
         ap3(np1,1) = ap3n1 + ddr11*two*rl*rzone(n)
         ap3(np1,2) =-(ap3n1 + ap3n3) + ddr12*two*rl*rzone(n) -
     $         (ddr13*two*rl*rzone(n) + ap3n3)*(rl+one)/(one+dlnr)
         ap3(np1,3) = zero
         ap4(np1) =-pi4g*rl1*rho(n)
c
c          horizontal momentum equation.
c
         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i-1)/r(i-1)
         ah1(i,2) = drdm(i-1,2)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i)/r(i)
         ah2(i) = g3m1(i-1)
         ah3(i) = one
         ah4(i) = rkl2(i-1)*g1(i-1)*p(i-1)*v(i-1)
      else
c
c          horizontal momentum equation.
c
         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i-1)/r(i-1)
         ah1(i,2) = drdm(i-1,2)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i)/r(i)
         ah2(i) = g3m1(i-1)
         ah3(i) = one
         ah4(i) = rkl2(i-1)*g1(i-1)*p(i-1)*v(i-1)
c
c          radial momentum equation.
c
         afac = forpi*r(i)**3/dm2(i)
         ag1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*(ahf*rl*v(i-1)-afac)
         ag1(i,2)=-for*gor(i)/r(i)+pi4g*ahf*(rho(i-1)+rho(i))+
     $     drdm(i,1)*g1(i)*p(i)*(ahf*rl*v(i) + afac) +
     $     drdm(i-1,2)*g1(i-1)*p(i-1)*(ahf*rl*v(i-1) - afac )
         ag1(i,3) = drdm(i,2)*g1(i)*p(i)*(ahf*rl*v(i)+afac)
         ag2(i,1) = g3m1(i-1)*(ahf*rl - rho(i-1)*afac)
         ag2(i,2) = g3m1(i)*(ahf*rl + rho(i)*afac)
         ag3(i,1) = ahf*rl - r(i)/drz(i)
         ag3(i,2) = ahf*rl + r(i)/drz(i)
         ag4(i,1) = rkl2(i-1)*(g1(i-1)*p(i-1)*(ahf*rl*v(i-1) - afac)) +
     $     ahf*rl1*gor(i)/r(i)
        ag4(i,2) = rkl2(i)*(g1(i)*p(i)*(ahf*rl*v(i)+afac)) +
     $     ahf*rl1*gor(i)/r(i)
c
c          poisson equation.
c
         ddr11 =-drz(i)/(drz(i-1)*(drz(i-1)+drz(i)))
         ddr12 = (drz(i)-drz(i-1))/(drz(i-1)*drz(i))
         ddr13 = drz(i-1)/(drz(i)*(drz(i-1)+drz(i)))
         if( i .eq. 2 ) then
            ap1(i,1) =-pi4g*rho(i-1)*(rzone(i-1)**2*drdm(i-1,1)-
     $     ahf*alog(rho(i-1)/rhoc)/alog(rzone(i-1)/rz0) )
         else
            ap1(i,1) =-pi4g*rho(i-1)*(rzone(i-1)**2*drdm(i-1,1)-
     $     ahf*alog(rho(i-1)/rho(i-2))/alog(rzone(i-1)/rzone(i-2)) )
         endif
         ap1(i,2) =-pi4g*rho(i-1)*(rzone(i-1)**2*drdm(i-1,2)-
     $     ahf*alog(rho(i)/rho(i-1))/alog(rzone(i)/rzone(i-1)) )
         ap3i1 = r(i-1)**2/( drint(i-1)*drz(i-1))
         ap3i3 = r(i)**2/(drint(i)*drz(i))
         ap3(i,1) = ap3i1 +           ddr11*two*rl*rzone(i-1)
         ap3(i,2) =-(ap3i1 + ap3i3) + ddr12*two*rl*rzone(i-1)
         ap3(i,3) = ap3i3 +           ddr13*two*rl*rzone(i-1)
         ap4(i) =-pi4g*rl1*rho(i-1)
c
c          the thermal equation matrices.
c
         dtl1 = t(i-1)**4
         dtl2 = t(i)**4
         t4oki = dtl2/gkap(i)
         t4oki1 = dtl1/gkap(i-1)
         diff = t4oki-t4oki1
         wiowi1 = alog(dtl2/dtl1)
         wowsq = wiowi1**2
         gkogk1 = alog(gkap(i)/gkap(i-1))
         denom = one-gkogk1/wiowi1
         dl1(i,1) = ((-for*g3m1(i-1)+dkdr(i-1))*(t4oki1/diff)
     $    + (-dkdr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom)
     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
         dl1(i,2) = ((for*g3m1(i)-dkdr(i))*(t4oki/diff)
     $    + (dkdr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom)
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
         dl2(i,1) = ((-for+dkdt(i-1))*(t4oki1/diff) -
     $    (dkdt(i-1)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i-1)*t(i-1))
     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
         dl2(i,2) = ((for-dkdt(i))*(t4oki/diff) +
     $    (dkdt(i)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i)*t(i))
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
c
c          keep track of the luminosity as a function of mass
c
            rlumx = rlumx1
            glom = frft(i-1)*rlumx/dm1(i-1)
     $         *(one-ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
            rlumx1 = rlumx + sorce(i-1)*dm1(i-1)
            if( noburn .ne. 0 ) then
               ipass = i
               rlumx1 = cool(ipass,rm(ipass)/totmas)
            endif
            glom1 = frft(i)*rlumx1/dm1(i-1)
     $         *(one+ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
            dldm = (rlumx1-rlumx)/dm1(i-1)
            aconlm = ahf*((one-frft(i-1))*rlumx+(one-frft(i))*rlumx1)
c
c          convert the nuclear energy derivatives.
c
         dedr = (-dvsorc(i-1)+g3m1(i-1)*dtsorc(i-1))*sorce(i-1)
         deds = dtsorc(i-1)*sorce(i-1)/(cv(i-1)*t(i-1))
         if( i .eq. 2 ) then
            ak1(i,1) = zero
         else
            ak1(i,1) = glom*dl1(i-1,1)*drdm(i-2,1)
         endif
c         if( i .eq. 3 ) ak1(i,1) = zero
         ak1(i,2) =-glom1*dl1(i,1)*drdm(i-1,1)+glom*for/r(i-1)**2 +
     $     glom*(dl1(i-1,1)*drdm(i-2,2)+dl1(i-1,2)*drdm(i-1,1) )
     $     +  dedr*drdm(i-1,1)
         ak1(i,3) = glom*dl1(i-1,2)*drdm(i-1,2) - glom1*for/r(i)**2 -
     $     glom1*(dl1(i,1)*drdm(i-1,2)+dl1(i,2)*drdm(i,1) )
     $     + dedr*drdm(i-1,2)
         ak1(i,4) =-glom1*dl1(i,2)*drdm(i,2)
         if( i .eq. 2 ) then
            ak4(i,1) = zero
         else
            ak4(i,1) = glom*rkl2(i-2)*dl1(i-1,1)
         endif
         ak4(i,2) =-rkl2(i-1)*(glom1*dl1(i,1) -  dedr
     $     - glom*dl1(i-1,2) )
         ak4(i,3) =-glom1*rkl2(i)*dl1(i,2)
         ak2(i,1) = glom*dl2(i-1,1)
         ak2(i,2) = glom*dl2(i-1,2)-glom1*dl2(i,1) + deds
         ak2(i,3) =-glom1*dl2(i,2)
c
c          add horizontal heat flow terms.
c
         rlom = frft(i-1)*rlumx
         rlom1 = frft(i)*rlumx1
         horiz = rkl2(i-1)*for*ac3*t4oki1*v(i-1)**2
         ak1(i,2) = ak1(i,2) - horiz*g3m1(i-1)*drdm(i-1,1) -
     $       rkl2(i-1)*ahf*rlom*v(i-1)/(forpi*r(i-1)**3)
         ak1(i,3) = ak1(i,3) - horiz*g3m1(i-1)*drdm(i-1,2) -
     $       rkl2(i-1)*ahf*rlom1*v(i-1)/(forpi*r(i)**3)
         ak2(i,2) = ak2(i,2) - horiz/(cv(i-1)*t(i-1))
         ak4(i,2) = ak4(i,2) + rkl2(i-1)*(dldm - horiz*g3m1(i-1) -
     $      aconlm*(v(i-1)/(forpi*rzone(i-1)**3)) )
      endif
  20  continue
c
c        write the matrices to file 6.
c
      if( iout .ge. 2 ) then
         write(11,2200)
         do 27 i=1,np1
            write(11,2201) i,ag1(i,1),ag1(i,2),ag1(i,3),ag2(i,1),
     $         ag2(i,2),ag3(i,1),ag3(i,2),ag4(i,1),ag4(i,2),ah4(i)
  27     continue
         write(11,2300)
         do 28 i=1,np1
            write(11,2201) i,ah1(i,1),ah1(i,2),ah2(i),ah3(i),ap1(i,1),
     $         ap1(i,2),ap3(i,1),ap3(i,2),ap3(i,3),ap4(i)
  28     continue
         write(11,2600)
         do 26 i=1,np1
            write(11,2201) i,ak1(i,1),ak1(i,2),ak1(i,3),ak1(i,4),
     $          ak2(i,1),ak2(i,2),ak2(i,3),ak4(i,1),ak4(i,2),ak4(i,3)
  26     continue
      endif
      return
c
 2200 format(1h1,//,20h ag1,ag2,ag3,ag4,ah4)
 2201 format(1x,i5,1p,10e12.4)
 2300 format(1h1,//,24h ah1,ah2,ah3,ap1,ap3,ap4)
 2600 format(1h1,//,10x,36h the pulsation matrices, ak1,ak2,ak4,/)
      end
      function nodes(w,nz,lout)
      implicit real(a-h,o-z)
c
c finds no. of nodes(sign changes) in nz elements of array w
c
      dimension w(nz),inode(300)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      nodcnt = 0
      nz1 = nz-1
      do 10 i=3,nz1
         if( w(i-1)*w(i) .ge. zero ) goto 10
c zero crossing found. is it general enough to be a node#
         if( w(i-2)*w(i-1) .lt. zero ) goto 10
         if( w(i)*w(i+1) .lt. zero ) goto 10
          nodcnt = nodcnt + 1
          inode(nodcnt) = i-1
  10  continue
      inode(nodcnt+1) = nz
      nodes = nodcnt
      if( lout .eq. 0 ) return
      write(11,1000) nodes
      if( nodcnt .le. 0 ) return
      write(11,1001) (inode(i),i=1,nodcnt)
 1000 format(/,1x,i2,12h nodes found)
 1001 format(7h nodes=,(30i4))
      return
      end
      function cool(qc)
      implicit real(a-h,o-z)
c
c      returns the value of luminosity out of a zone as
c   a given function of qc, q1, and q0.
c
c          noburn:= 0, normal nuclear profile.
c                   1, linear ramp.
c                   2, quadratic ramp.
c                   3, square root ramp.
c
c          needs in common block coolum:
c                 onemq0:= 1.-q0(inner mass fraction)
c                 onemq1:= 1.-q1(outer mass fraction)
c                 rlumc:= core luminosity
c                 rlums:= surface luminosity
c
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      if( noburn .gt. 0 ) then
         if( qc .gt. 1.e0-onemq1 ) then
            cool = rlums
         elseif( qc .lt. 1.e0-onemq0 ) then
            cool = rlumc
         else
            f = (qc+onemq0-1.e0)/(onemq0-onemq1)
            if( noburn .eq. 1 ) cool = rlumc+(rlums-rlumc)*f
            if( noburn .eq. 2 ) cool = rlumc+(rlums-rlumc)*f*f
            if( noburn .eq. 3 ) cool = rlumc+(rlums-rlumc)*sqrt(f)
         endif
      endif
      return
      end
      subroutine orthog(dr,dh,n,omsq,index,imode)
      implicit real(a-h,o-z)
c
c          check for orthogonality of the wave functions. part
c      one stores the vertical and horizontal wavefunctions,
c      part two does the integrations.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $   dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/eigstr/ sdr(nmax,4),sdh(nmax,4),stome(4),iomeg(4)
      dimension dr(nmax),dh(nmax),rke(4),cross(4,4)
      common/stnmod/ nmode
c
      if( n .lt. 0 ) goto 100
      if( imode .le. 0 ) goto 50
      if( nmode .ge. 4 ) return
      if( nmode .gt. 0 ) then
         do 5 i=1,nmode
            if( abs((omsq-stome(i))/stome(i)) .lt. 1.e-6 ) return
   5     continue
      endif
      nmode = nmode + 1
      do 10 i=1,n+1
         sdr(i,nmode) = dr(i)
         sdh(i,nmode) = dh(i)
  10  continue
      stome(nmode) = omsq
      iomeg(nmode) = index
      return
  50  continue
c
c          part two: the integrals.
c
      if( nmode .le. 1 ) return
      lval = iabs(imode)
      rl1 = float(lval)*float(lval+1)
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
         rke(j1) = sqrt(cross(j1,j1))
  75  continue
      do 80 j1=1,nmode
         do 80 j2=1,nmode
            cross(j2,j1) = cross(j2,j1)/(rke(j1)*rke(j2))
  80  continue
      write(1,8000)
      do 85 j1=1,nmode
         do 85 j2=1,nmode
            write(1,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
  85  continue
      return
c
c          special case for l=0.
c
 100  continue
      nuse = iabs(n)
      if( imode .le. 0 ) goto 150
      if( nmode .ge. 4 ) return
      nmode = nmode + 1
      do 110 i=1,nuse+1
         sdr(i,nmode) = dr(i)
 110  continue
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
         rke(j1) = sqrt(cross(j1,j1))
 175  continue
      do 180 j1=1,nmode
         do 180 j2=1,nmode
            cross(j2,j1) = cross(j2,j1)/(rke(j1)*rke(j2))
 180  continue
      write(1,1800)
      do 185 j1=1,nmode
         do 185 j2=1,nmode
            write(1,8001) iomeg(j1),iomeg(j2),cross(j1,j2)
 185  continue
      return
c
 1800 format(1x,34h radial cross integrals for modes:)
 8000 format(1x,26hcross integrals for modes:)
 8001 format( 2(3x,1h<,i3,1h,,i3,1h>,1p,2e12.5) )
      end
      subroutine eigenf(omsq,nin,lin,rke,gor)
      implicit real(a-h,o-z)
c
c          find the normalized eigenvectors and the weight functions
c       as a function of position.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $   dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/linear/ dr(nmax),dh(nmax),gam(nmax),qnyr(nmax),gyr(nmax),
     $               dp(nmax),adrho(nmax),weight(nmax),cyr(nmax),
     $               xo(nmax),yo(nmax),qwait(nmax),ra(nmax),vn(nmax),
     $               un(nmax),y1(nmax),y2(nmax),y3(nmax),y4(nmax),
     $               spac1(nmax,8)
      dimension vq(6),gor(nmax)
c
      n = nin
      lval = lin
      rl = float(lval)
      rl1 = rl*(rl + one)
      do 10 i=1,n
         vn(i) = ahf*(r(i+1)*gor(i+1) + r(i)*gor(i))/(p(i)*v(i))
         un(i) = forpi*r(i+1)**3/rm(i+1)/v(i)
         ra(i) =-bv(i+1)*r(i+1)/gor(i+1)
         xo(i) = alog(rzone(i)/p(i))
         y1(i) = dr(i+1)
         y2(i) = omsq*dh(i)/(ahf*(gor(i)/r(i)+gor(i+1)/r(i+1)))
         y3(i) = gam(i)/(ahf*(gor(i)*r(i)+gor(i+1)*r(i+1)))
         y4(i) = ((gam(i+1)-gam(i))/(rzone(i+1)-rzone(i)))/gor(i)
  10  continue
c
c          weight functions
c
      stwait = zero
      do 20 i=1,n
         vq(1) = ahf*(gor(i)/r(i) + gor(i+1)/r(i+1))
         vq(2) = ahf*(vn(i) + vn(i+1))/g1(i)
         vq(3) =-ra(i)
         vq(4) = ahf*(un(i)+un(i+1))
         vq(5) = one/(one + vq(2)*g1(i))
         vq(6) = rzone(i)
         if( i .eq. 1 ) then
            dx = xo(2) - xo(1)
         elseif( i .eq. n ) then
            dx = xo(n) - xo(n-1)
         else
            dx = ahf*(xo(i+1)-xo(i-1))
         endif
         temp = forpi*dx*rho(i)*vq(1)*vq(5)*vq(6)**5
         cyr(i) = temp*vq(2)*(y2(i) - y3(i))**2
         qnyr(i) = temp*vq(3)*y1(i)**2
         gyr(i) =-temp*(y4(i) + (rl+one)*y3(i))**2/vq(4)
         gyr(i) =-temp*(y4(i)**2 + rl1*y3(i)**2)/vq(4)
         weight(i) =  cyr(i) + qnyr(i) + gyr(i)
         stwait = stwait + weight(i)
  20  continue
      do 30 i=1,n
          weight(i) = weight(i)/(omsq*rke)
          cyr(i) = cyr(i)/(omsq*rke)
          qnyr(i) = qnyr(i)/(omsq*rke)
          gyr(i) = gyr(i)/(omsq*rke)
  30  continue
      stwait = stwait/rke
      qch = abs((omsq-stwait)/omsq)
      qrke = rke*omsq/2.0e0
      write(6,2000) omsq,stwait,qch,qrke
      write(1,2000) omsq,stwait,qch,qrke
      write(11,2000) omsq,stwait,qch,qrke
      call pltdmp(weight,nmax,n,'epst')
      call pltdmp(cyr,nmax,n,'cyr ')
      call pltdmp(qnyr,nmax,n,'nyr ')
      call pltdmp(gyr,nmax,n,'gyr ')
      return
c
 2000 format(2x,32heigenvalue from matrix solution=,1pe12.3,/,
     $       2x,32heigenvalue from weight function=,e12.3,
     $       7h error=,e12.3,/,
     $       1x,21hkinetic energy amp. =,e12.4)
      end
      subroutine cvtm(n,p,rl,cv,t,dm,x,nmax,itrans)
      implicit real(a-h,o-z)
c
c          compute x=sum(cv*t*dm/(rl*p)) from outside in.
c
      dimension cv(nmax),t(nmax),dm(nmax),x(nmax)
      x(n) = cv(n)*t(n)*dm(n)/(rl*p)
      ist = n
      ido = n-1
      do 10 i=1,ido
         j = n - i
         x(j) = x(j+1) + cv(j)*t(j)*dm(j)/(rl*p)
         if( x(j) .gt. 1.e0 ) goto 10
         ist = j
  10  continue
      itrans = ist
      return
      end
      function fomsq(omsq,x,xnorm,n)
      implicit real(a-h,o-z)
c
      parameter ( nmax=600 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),            ap4(nmax)
      common/linear/ atrix(nmax3,7), spac(nmax3,2)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax3)
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega**2.
c
      np1 = n + 1
      do 200 i=1,np1
c
c      ix is index for radial component, idh for horizontal and
c    igam is for the poisson equation.
c
         ix = 3*i
         idh = ix - 1
         igam = ix - 2
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
      x(3*np1-1) =-ah1(np1,2)*xnorm
      x(3*np1-2) =-ap1(np1,2)*xnorm
      x(3*np1-3) =-ag1(n,3)*xnorm
      do 210 i=1,7
         atrix(3*np1,i) = zero
 210  continue
      call rbmles(atrix,3*nmax,1,3*np1-1,7,x)
      err = ag1(np1,1)*x(3*np1-3) + (ag1(np1,2)-omsq)*xnorm +
     $      ag3(np1,1)*x(3*np1-2) + ag4(np1,1)*x(3*np1-1)
      fomsq = err/x(3)
      return
      end
      subroutine rbmles(a,nmax,imin4,imax4,m4,y)
      implicit real (a-h,o-z)
c
c          linear equation solver. (r)eal (b)and (m)atrix (l)inear
c       (e)quation (s)olver. solves a*x = y, destroying a and
c       returns the value of x in y.
c
      dimension a(nmax,m4),y(nmax)
      imin = imin4
      imax = imax4
      m = m4
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
            do 15 kk=1,kmax
               jk = j+kk
               k = mmid+kk
               ik = i-kk
               a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  15        continue
  10     continue
         do 25 kk=1,kmax
            ik = i-kk
            k = mmid+kk
            y(ik) = y(ik)-y(i)*a(ik,k)
  25     continue
  20  continue
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
         kmax = min0(i-imin,mmm)
         do 35 kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
  35     continue
  30  continue
      return
      end
      subroutine trisol(a,imin4,imax4,x,y,z,nmax,e,f)
      implicit real (a-h,o-z)
c
c          Solves the matrix equation (A-xI)y = z, returning
c       y given A, x, and z. Arrays e and f are scratch matrices
c       for storing the recursion relation.
c
      dimension a(nmax,3), y(nmax), z(nmax)
      dimension e(nmax),f(nmax)
      imin = imin4
      imax = imax4
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
      subroutine cbmles(a,nmax,imin,imax,m,y)
      complex*16 a,y,den
c
c          linear equation solver. (c)omplex (b)and (m)atrix (l)inear
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
c         den = (1.d0,0.d0)
c         if( dreal(a(i,mmid)).ne.0.d0 .or. dimag(a(i,mmid)).ne.0.d0)
c     $      den = 1.d0/a(i,mmid)
         den = 1.d0/a(i,mmid)
         y(i) = y(i)*den
         kmax = min0(i-imin,mmm)
         do 10 j=1,mmm
            a(i,j) = a(i,j)*den
            do 15 kk=1,kmax
               jk = j+kk
               k = mmid+kk
               ik = i-kk
               a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  15        continue
  10     continue
         do 25 kk=1,kmax
            ik = i-kk
            k = mmid+kk
            y(ik) = y(ik)-y(i)*a(ik,k)
  25     continue
  20  continue
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
         kmax = min0(i-imin,mmm)
         do 35 kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
  35     continue
  30  continue
      return
      end
      subroutine trcsol(a,nmax,imin,imax,x,y,z,e,f)
      implicit real (a-h,o-z)
      complex*16 x,z,e,f,den
      dimension a(nmax,3),y(nmax),z(nmax)
      dimension e(nmax),f(nmax)
c
c          tridiagonal matrix solver for complex arrays.
c
      ido = imax+1 - imin
      e(imax) =-a(imax,1)/(a(imax,2)-x)
      f(imax) = y(imax)/(a(imax,2)-x)
      do 10 ii=2,ido
         i = imax+1 - ii
         den = (a(i,2)-x) + a(i,3)*e(i+1)
         if( ii .ne. ido ) then
            e(i) =-a(i,1)/den
         endif
         f(i) = (y(i)-a(i,3)*f(i+1))/den
  10  continue
      z(imin) = f(imin)
      ido = imin+1
      do 20 i=ido,imax
         z(i) = e(i)*z(i-1) + f(i)
  20  continue
      return
      end
      subroutine cmplna(io,ihmin,ihmax,nin)
      implicit real(a-h,o-z)
c
c     *****************************************************************
c      pulsation equation solver
c      basic references -- castor, ap.j. 166 (109) 1971.
c                          pesnell, 1983, thesis, univ. of florida.
c                          pesnell, pasp, 99 (975), 1987.
c     *****************************************************************
c
c          this version of the lna radial stability analysis is suitable
c       for envelope models generated with the "model" building code.
c       the inner boundary is valid for a complete model and models
c       of low central condensation. If the values of the central ball
c       are not present, or inconsistent, the results of the stability
c       analysis are questionable.
c
      complex x0,y0,y1,z,xs,dror,dlum,omsq1,comega,drho,dpress,
     $  dtemp,amat,oma,erra,err,omb,czero,czi,dlth,dlmech
c
c        arrays are rearranged from the original lna code
c
      parameter ( nmax=600, nmax2=2*nmax)
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $       dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
     $               drz(nmax),drint(nmax),gor(nmax)
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ak1(nmax,4),ak2(nmax,3),
     $       bl2(nmax,2),bl1(nmax,3),dr(nmax,2)
      common xs(nmax2),y0(nmax),y1(nmax),x0(nmax),dror(nmax),
     $   dlth(nmax),dlmech(nmax),dlum(nmax),drho(nmax),xamp(nmax),
     $   xpha(nmax),rhoamp(nmax),rhopha(nmax),pramp(nmax),prpha(nmax),
     $   temamp(nmax),tempha(nmax),glamp(nmax),glpha(nmax),
     $   w(nmax),cvt(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk37/  yo(nmax),xo(nmax),work(nmax),dro(nmax),dp(nmax),
     $       dt(nmax),dl(nmax)
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      dimension amat(nmax2,7)
      dimension zzero(nmax,36),szero(nmax,19)
      equivalence (zzero(1,1),xs(1)),(szero(1,1),ag1(1,1))
      data czero/(0.0e0,0.0e0)/,accur/1.e-10/
c
c   zero common blocks scrtch and (blank).
c
      do 5 j=1,36
         do 5 i=1,nmax
            zzero(i,j) = zero
   5  continue
      do 6 j=1,19
         do 6 i=1,nmax
            szero(i,j) = zero
   6  continue
c
      n = nin
      np = n+1
      nm1 = n-1
c
      rlumx = corlum
      afac = forpi*r(1)**2/sqrt(dm2(1))
      do 60 i=1,np
         afacp = afac
         if( i .eq. np ) goto 40
         afac = forpi*r(i+1)**2/sqrt(dm2(i+1))
c
c          change dln(kappa)/dln(rho)<t> to dln(kappa)/dln(rho)<s>
c
         dkdr(i) = dkdr(i)+g3m1(i)*dkdt(i)
         dr(i,1) = afacp/v(i)/dm1(i)
         dr(i,2) =-afac/v(i)/dm1(i)
         if( i .eq. 1 ) then
            dtl1 = tc**4
            dtl2 = t(i)**4
            t4oki = dtl2
            t4oki1 = dtl1
            diff = t4oki-t4oki1
            dkdr0 = dkdr0 + g3m10*dkdt0
            bl1(i,1) = zero
            bl1(i,2) = ((for*g3m10-dkdr0)*(t4oki1/diff)*thre
     $             + for)/(r(i)*sqrt(dm2(i))) +
     $     dr(i,1)*(for*g3m1(i)-dkdr(i))*(t4oki/diff)
            bl1(i,3) = dr(i,2)*(for*g3m1(i)-dkdr(i))*(t4oki/diff)
            bl2(i,1) = (-for+dkdt(i))*(t4oki1/diff)/(cv0*tc)
            bl2(i,2) = (for-dkdt(i))*(t4oki/diff)/(cv(i)*t(i))
            drdm0 =-thre/(r(1)*sqrt(dm2(1)))
            qafacp = forpi*r(1)**2/sqrt(dm2(1))
            qafac = forpi*r(2)**2/sqrt(dm2(2))
            ag1(1,1) = zero
            ag1(1,2) =-for*g*rm(1)/r(1)**3 + qafacp*(g1(1)*p(1)*dr(1,1)
     $               - drdm0*pc*g10 )
            ag1(1,3) = qafacp*g1(1)*p(1)*dr(1,2)
            ag2(i,1) =-qafacp*g3m10*rhoc
            ag2(i,2) = qafacp*g3m1(1)/v(1)
           goto 50
        endif
c
c          set up the mechanical matrices.
c
         g1p = g1(i)*p(i)
         g3m1ov = g3m1(i)/v(i)
         if( i .eq. n ) then
            g1p = g1p-1.00855d-14*t(i)**4*g3m1(i)
            g3m1ov = g3m1ov-1.00855d-14*t(i)**3/cv(i)
         endif
         ag1(i,1) =-afacp*g1(i-1)*p(i-1)*dr(i-1,1)
         ag1(i,2) =-for*g*rm(i)/r(i)**3 + afacp*(g1p*dr(i,1)-
     $               g1(i-1)*p(i-1)*dr(i-1,2) )
         ag1(i,3) = afacp*g1p*dr(i,2)
         ag2(i,1) =-afacp*g3m1(i-1)/v(i-1)
         ag2(i,2) = afacp*g3m1ov
c
c          initialize the thermal matrices, if irad=1 use the
c       stellingwerf interpolation formula.
c
         dtl1 = t(i-1)**4
         dtl2 = t(i)**4
         t4oki = dtl2/gkap(i)
         t4oki1 = dtl1/gkap(i-1)
         diff = t4oki-t4oki1
         wiowi1 = alog(dtl2/dtl1)
         wowsq = wiowi1**2
         gkogk1 = alog(gkap(i)/gkap(i-1))
         denom = one - gkogk1/wiowi1
       bl1(i,1) = dr(i-1,1)*((-for*g3m1(i-1)+dkdr(i-1))*(t4oki1/diff)
     $  + (-dkdr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom)
c        if( i .eq. 2 ) bl1(2,1) = zero
         bl1(i,2) = for/(r(i)*sqrt(dm2(i))) +
     $     dr(i,1)*((for*g3m1(i)-dkdr(i))*(t4oki/diff) +
     $    (dkdr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom ) +
     $  dr(i-1,2)*((-for*g3m1(i-1)+dkdr(i-1))*(t4oki1/diff) +
     $    (-dkdr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom )
         bl1(i,3) = dr(i,2)*((for*g3m1(i)-dkdr(i))*(t4oki/diff)
     $     + (dkdr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom)
         bl2(i,1) = ((-for+dkdt(i-1))*(t4oki1/diff)
     $  - (dkdt(i-1)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i-1)*t(i-1))
         bl2(i,2) = ((for-dkdt(i))*(t4oki/diff)
     $  + (dkdt(i)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i)*t(i))
         goto 50
  40    continue
         tau = gkap(n)*dm1(n)/(for*forpi*r(np)*r(np))
         tau = 1.5d0*tau/(one+1.5d0*tau)
         bbb = (r(np)-rphoto)/r(np)
         ef1 = (one+ahf*bbb)/(one+bbb)
         ef3 = two/(one+bbb)
         bl1(np,1) = ef1*dr(n,1)*( for*g3m1(n) - tau*(one+dkdr(n)) )
         bl1(np,2) = (ef3-tau*ef1*two)/(r(np)*sqrt(dm2(np))) +
     $      dr(n,2)*ef1*(for*g3m1(n) - tau*(one+dkdr(n)))
         bl1(np,3) = zero
         bl2(np,1) = ef1*(for - tau*dkdt(n))/(cv(n)*t(n))
         bl2(np,2) = zero
         g1p = g1(n)*p(n)
         g3m1ov = g3m1(n)/v(n)
         g1p = g1p-1.00855d-14*t(n)**4*g3m1(n)
         g3m1ov = g3m1ov-1.00855d-14*t(n)**3/cv(n)
         ag1(np,1) =-afacp*g1p*dr(n,1)
         ag1(np,2) =-for*g*rm(np)/r(np)**3-afacp*g1p*dr(n,2)
         ag1(np,3) = zero
         ag2(np,1) =-afacp*g3m1ov
         ag2(np,2) = zero
  50    continue
         if( i .eq. 1 ) then
            glom1 = frft(i)*sorc0
            dedr =-(dedv0/rhoc-dedt0*tc*g3m10)
            deds = dedt0/cv0
            ak1(i,1) = zero
            ak1(i,2) = zero
            ak1(i,3) =-glom1*bl1(i,2)+dedr*drdm0
            ak1(i,4) =-glom1*bl1(i,3)
            ak2(i,1) = zero
            ak2(i,2) =-glom1*bl2(i,1) + deds
            ak2(i,3) =-glom1*bl2(i,2)
         elseif( i .eq. np ) then
            glom = frft(i-1)*rlumx/dm1(i-1)
            glom1 = rlumx/dm1(i-1)
            ak1(i,1) = glom*bl1(i-1,1)
            ak1(i,2) = glom*bl1(i-1,2)-glom1*bl1(i,1)
            ak1(i,3) = glom*bl1(i-1,3)-glom1*bl1(i,2)
            ak1(i,4) = zero
            ak2(i,1) = glom*bl2(i-1,1)
            ak2(i,2) = glom*bl2(i-1,2)-glom1*bl2(i,1)
            ak2(i,3) = zero
         else
            glom = frft(i-1)*rlumx/dm1(i-1)
            rlumx = rlumx + sorce(i-1)*dm1(i-1)
            glom1 = frft(i)*rlumx/dm1(i-1)
            dedr =-(dvsorc(i-1)*v(i-1)-dtsorc(i-1)*t(i-1)*g3m1(i-1))
            deds = dtsorc(i-1)/cv(i-1)
            ak1(i,1) = glom*bl1(i-1,1)
            ak1(i,2) = glom*bl1(i-1,2)-glom1*bl1(i,1)+dedr*dr(i-1,1)
            ak1(i,3) = glom*bl1(i-1,3)-glom1*bl1(i,2)+dedr*dr(i-1,2)
            ak1(i,4) =-glom1*bl1(i,3)
            ak2(i,1) = glom*bl2(i-1,1)
            ak2(i,2) = glom*bl2(i-1,2)-glom1*bl2(i,1) + deds
            ak2(i,3) =-glom1*bl2(i,2)
         endif
  60  continue
c
c          write the matrices to file tape11, if io is less
c       than 3, this write is not done.
c
      if( io .ge. 3 ) then
         write(11,7000)
         do 75 i=1,np
          write(11,7001) i,ag1(i,1),ag1(i,2),ag1(i,3),ag2(i,1),ag2(i,2),
     $    ak1(i,1),ak1(i,2),ak1(i,3),ak1(i,4),ak2(i,1),ak2(i,2),ak2(i,3)
  75     continue
      endif
c
c          calculate the acoustic tavel time from surface to the
c       innermost zone. this is stored in transt.
c
      transt = zero
      do 77 i=2,np
         transt = transt+(r(i)-r(i-1))/sqrt(p(i-1)*v(i-1)*g1(i-1))
  77  continue
      rhom = totmas/(pi43*r(np)**3)
      termq = sqrt(rhom/1.41e0)
      tcon = pi43*g*rhom*float(n)**2
      x0np = r(np)*sqrt(dm2(np))
      omsqp = zero
      omsqc = 5.e-4*float(ihmin)*tcon
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
         call trisol(ag1,1,n,omsq,yo,xo,nmax,dt,dl)
        if( icount .eq. 1 ) then
c
c          set up initial derivative
c
            qerra = ag1(np,1)*xo(n)+(ag1(np,2)-omsq)*x0np
            qerra = qerra/xo(2)
            omsql = omsq
            omsq  = omsq*(one+1.e-7)
         else
            qerr = ag1(np,1)*xo(n)+(ag1(np,2)-omsq)*x0np
            afacp = abs((omsq-omsql)/omsql)
            qerr = qerr/xo(2)
           if(io.ge.1) write(11,1400) icount,omsq,omsql,afacp
           if( afacp .lt. 1.e-11 ) goto 86
            omsq2 = (qerra*omsq-qerr*omsql)/(qerra-qerr)
            omsql = omsq
            domsq = omsq2-omsq
            domsq = sign(amin1(abs(domsq),qrt*abs(omsq)),domsq)
            omsq  = omsq + domsq
            qerra  = qerr
         endif
  81  continue
c
c          no convergence in the adiabatic eigenvalue.
c
      write(6,8100) omsq,iterad,iv,ih
      omsqp = zero
      omsqc = 5.e-6*tcon*float(ih+1)
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
         write(6,8600) omsq,iv,ih
         omsqp = zero
         omsqc = 5.e-6*tcon*float(ih+1)
         goto 79
      endif
c
      do 90 i=1,n
         yo(i) = zero
  90  continue
      xo(np) = x0np
      yo(n)  =-ag1(n,3)*x0np
      call trisol(ag1,1,n,omsq,yo,xo,nmax,dt,dl)
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
      period = twopi/(omega*86400.e0)
      write(11,8000) iv,ih,omega,period
      omsqp = zero
      omsqc = 5.e-6*tcon*float(ih)
      goto 79
c
 103  continue
      omega = sqrt(omsq)
      w(1) = xo(1)/(r(1)*sqrt(dm2(1)))
      dro(1) = xo(1)*drdm0
      dp(1) = g10*dro(1)
      dt(1) = g3m10*dro(1)
      glamp(1) =-for*(g*rm(1)/r(1))*(xo(1)/r(1))**2 +
     $      pc*g10*dro(1)**2*cormas/rhoc
      stwait = glamp(1)
      rke = xo(1)*xo(1)
      dl(1) = bl1(1,2)*xo(1) + bl1(1,3)*xo(2)
      do 100 i=2,np
         w(i) = xo(i)/(r(i)*sqrt(dm2(i)))
         dro(i) = dr(i-1,1)*xo(i-1)+dr(i-1,2)*xo(i)
         dp(i) = g1(i-1)*dro(i)
         dt(i) = g3m1(i-1)*dro(i)
c
c          set up the weight function calculation.
c
         glamp(i) =-for*(g*rm(i)/r(i))*(xo(i)/r(i))**2 +
     $      p(i-1)*v(i-1)*g1(i-1)*dro(i)**2*dm1(i-1)
         stwait = stwait + glamp(i)
         rke = rke + xo(i)*xo(i)
         if( i .eq. np ) goto 100
         dl(i) = bl1(i,1)*xo(i-1)+bl1(i,2)*xo(i)+bl1(i,3)*xo(i+1)
 100  continue
      dl(np) = bl1(np,1)*xo(n)+bl1(np,2)*x0np
      w(np) = one
      period = twopi/omega
      call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax)
      period = period/86400.e0
      qvalue = period*termq
      stwait = stwait/rke
      qch = abs((stwait-omsq)/omsq)
      write(6,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch
      write(1,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch
      write(11,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch
      iv = nodes(w,np,1)
c
c          normalize the weight function per zone to omsq.
c
      do 105 i=1,np
         glamp(i) = glamp(i)/(omsq*rke)
 105  continue
c
      if(io.ge.2) write(11,1000) omega,period,(i,w(i),dl(i),dro(i),
     $ dt(i),dp(i),cvt(i),glamp(i),i=1,np)
c
c          dump adiabatic eigenvectors to casplt.
c
      call pltdmp(w,nmax,np,4hdr/r)
      call pltdmp(dro,nmax,np,4hdrho)
      call pltdmp(dl,nmax,np,4hdl/l)
      call pltdmp(glamp,nmax,np,4hwait)
c
c          find first guess to the imaginary part of omega by the use
c       of the quasi-adiabatic approximation. i.e.:
c
c               del(omega**2) = xo*ag2*(i*omega-ak2)**-1*ak1*xo.
c
c       note that the minus signs in yo(i) are compensated for by the
c       routine 'trcsol' solving the resolvent as ak2 - i*omega.
c
      yo(1) =-ak1(1,2)*xo(1)-ak1(1,3)*xo(2)-ak1(1,4)*xo(3)
      y1(1) = czero
      xo(n+2) = zero
      do 110 i=2,n
         yo(i) =-ak1(i,1)*xo(i-1)-ak1(i,2)*xo(i)-ak1(i,3)*xo(i+1)-
     $         ak1(i,4)*xo(i+2)
         y1(i) = czero
 110  continue
      yo(np) =-ak1(np,1)*xo(nm1)-ak1(np,2)*xo(n)-ak1(np,3)*x0np
      call trcsol(ak2,nmax,1,np,(0.e0,1.e0)*omega,yo,y1,dlmech,dlth)
      czi = czero
      do 120 i=1,n
         czi = czi + xo(i)*(ag2(i,1)*y1(i)+ag2(i,2)*y1(i+1))
 120  continue
      czi = czi + x0np*ag2(np,1)*y1(np)
      omsq1 = omsq + czi/rke
      if( cabs(omsq1) .gt. 1.5d0*omsq ) omsq1 = omsq*(1.e0,1.e-5)
      comega = csqrt(omsq1)
      growth =-forpi*aimag(comega)/real(comega)
      write(1,1250) rke,comega,growth
      write(11,1250) rke,comega,growth
c
c          nonadiabatic iteration on mode iv starts here
c
      do 130 icount=1,20
      do 140 i=1,np
         iy = 2*i-1
         ix = 2*i
         xs(ix) = czero
         xs(iy) = czero
c
c          amat reset every iteration
c
         amat(iy,1) = ak1(i,1)
         amat(iy,2) = ak2(i,1)
         amat(iy,3) = ak1(i,2)
         amat(iy,4) = ak2(i,2)-(0.0e0,1.0e0)*comega
         amat(iy,5) = ak1(i,3)
         amat(iy,6) = ak2(i,3)
         amat(iy,7) = ak1(i,4)
         amat(ix,1) = czero
         amat(ix,2) = ag1(i,1)
         amat(ix,3) = ag2(i,1)
         amat(ix,4) = ag1(i,2)-omsq1
         amat(ix,5) = ag2(i,2)
         amat(ix,6) = ag1(i,3)
         amat(ix,7) = czero
 140  continue
      write(11,1400) icount,omsq,comega,omsq1
      xs(2*np-3) =-ak1(np-1,4)*x0np
      xs(2*np-2) =-ag1(n,3)*x0np
      xs(2*np-1) =-ak1(np,3)*x0np
      call cbmles(amat,2*nmax,1,2*np-1,7,xs)
      rke = zero
      do 150 i=1,n
         y1(i) = xs(2*i-1)
         x0(i) = xs(2*i)
         rke = rke + conjg(x0(i))*x0(i)
 150  continue
      y1(np) = xs(2*np-1)
      x0(np) = x0np
      rke = sqrt(abs(rke + x0np+x0np))
      if( icount .eq. 1 ) then
         oma = comega
         erra = ag1(np,1)*x0(n)+(ag1(np,2)-omsq1)*x0np+ag2(np,1)*y1(np)
         erra = erra/sign(rke,real(x0(2)))
         omsq1 = omsq1*(one + 1.e-7)
         comega = csqrt(omsq1)
      else
         err = ag1(np,1)*x0(n)+(ag1(np,2)-omsq1)*x0np+ag2(np,1)*y1(np)
         afacp = cabs((comega-oma)/oma)
         err = err/sign(rke,real(x0(2)))
         if( io .ge. 1 ) write(11,1851) icount,afacp,comega,err
         if( afacp .lt. accur ) goto 190
         omb = (erra*comega-err*oma)/(erra-err)
         oma = comega
         comega = omb
         omsq1 = comega*comega
         erra = err
      endif
 130  continue
      write(6,1300) iv,afacp,accur
      write(1,1300) iv,afacp,accur
      write(11,1300) iv,afacp,accur
c
c          iteration loop ends here, even if omega is not converged,
c       print the eigenvectors.
c
c          converged on mode iv
c
  190 continue
      rke = x0(1)*conjg(x0(1))
      drho(1) = drdm0*x0(1)
      call amppha ( drho(1),rhoamp(1),rhopha(1))
      dpress = drho(1)*g10+y1(1)*g3m10*rhoc/pc
      call amppha ( dpress,pramp(1),prpha(1))
      work(1) = aimag(conjg(y1(1))*drho(1))*cormas*g3m10
      dtemp = drho(1)*g3m10+y1(1)/(cv0*tc)
      call amppha ( dtemp,temamp(1),tempha(1))
      do 200 i=2,np
         dror(i) = x0(i)/(r(i)*sqrt(dm2(i)))
         call amppha ( dror(i),xamp(i),xpha(i))
         drho(i) = dr(i-1,1)*x0(i-1)+dr(i-1,2)*x0(i)
         call amppha ( drho(i),rhoamp(i),rhopha(i))
         dpress = drho(i)*g1(i-1)+y1(i)*g3m1(i-1)/v(i-1)/p(i-1)
         call amppha ( dpress,pramp(i),prpha(i))
         rke = rke + x0(i)*conjg(x0(i))
         work(i) = aimag(conjg(y1(i))*drho(i))*dm1(i-1)*g3m1(i-1)
         dtemp = drho(i)*g3m1(i-1)+y1(i)/(cv(i-1)*t(i-1))
         call amppha ( dtemp,temamp(i),tempha(i))
c
c          dlth and dlmech are the thermal and mechanical
c       contributions to the luminosity variations.
c
         dlmech(i) = bl1(i,1)*x0(i-1)+bl1(i,2)*x0(i)+bl1(i,3)*x0(i+1)
         dlth(i) = bl2(i,1)*y1(i-1) + bl2(i,2)*y1(i)
         dlum(i) = dlth(i) + dlmech(i)
         call amppha ( dlum(i),glamp(i),glpha(i))
         y1(i-1) = y1(i-1)/(cv(i-1)*t(i-1))
 200  continue
      dror(1) = x0(1)/(r(1)*sqrt(dm2(1)))
      call amppha(dror(1),xamp(1),xpha(1))
      dlum(1) = bl1(1,2)*x0(1)+bl1(1,3)*x0(2)
     $        + bl2(1,1)*y1(1)+bl2(1,2)*y1(2)
      call amppha(dlum(1),glamp(1),glpha(1))
      call cptdmp(dror,nmax,np,4hdr/r)
      call cptdmp(drho,nmax,n,4hdrho)
      call cptdmp(dlum,nmax,np,4hdl/l)
c
      omega = real(comega)
      rke = rke*omega*omega*ahf
      omimag = aimag(comega)
      growth =-forpi*aimag(comega)/omega
c
c          normalized work integral and integrated work function.
c
      work(np+1) = zero
      yo(np+1) = zero
      do 205 ii=1,np
         i = np - ii + 1
         yo(i) = yo(i+1) + work(i)*pi/rke
         work(i) = work(i)*pi/rke/abs(growth)
 205  continue
      qch = abs((growth-yo(1))/growth)
c
      do 215 i=1,np
         yo(i) = yo(i)/abs(growth)
 215  continue
      call pltdmp(work,nmax,np,4hwork)
      call pltdmp(yo,nmax,np,4hwint)
c
c          define the incarnations of the period.
c
      period = twopi/omega
      phase = transt*two/period
      pdays = period/86400.e0
c
      s = yo(1)*abs(growth)
      write(6,3000) comega,growth,s,qch,period,pdays,phase
      write(1,3000) comega,growth,s,qch,period,pdays,phase
      write(11,3000) comega,growth,s,qch,period,pdays,phase
c
      if(io.gt.1) write(11,3100) omega,period,growth,(i,xamp(i+1),
     $  xpha(i+1),glamp(i+1),glpha(i+1),rhoamp(i),rhopha(i),temamp(i),
     $  tempha(i),pramp(i),prpha(i),work(i),yo(i),i=1,n)
      afac = omega**2/(pi43*g*rhom)
      afacp = pdays*termq
      write(11,3200) afac,afacp
      if(io.gt.1) write(11,3300) (i,dror(i+1),dlth(i+1),dlmech(i+1),
     $  dlum(i+1),drho(i),y1(i),i=1,n)
      call workmx(yo,np,nmax)
      write(11,3600) rke
      call zie(work,n,nmax)
      omsqp = omsq
      omsqc = omsq
  79  continue
      return
c
 1250 format(25h quasi-adiabatic results:,/,5x,3h j=,1pe11.3,
     $   7h omega=,2e12.4,14h growth rate =,e12.4)
 1400 format(1x,i4,1p,7e11.4)
 1851 format(1x,10hlna iter =,i3,7h afacp=,1pe10.3,8h comega=,
     $   2e10.3,5h err=,2e10.3)
 3000 format(/,1x,25hlinear non-adiabatic mode,/,8h omega =,1pe16.9,
     $   1x,9h(sec**-1),3x,15h omega (imag) =,e16.9,/,
     $   9x,24h predicted growth rate =,e14.6,/,
     $   1x,32hgrowth rate from work integral =,e14.6,7h error=,e11.3,
     $   //,3x,9hperiods =,e12.4,7h (secs),e13.4,7h (days),/,
     $   1x,23h acoustic travel time =,e13.6,8h periods)
 3100 format(1h1,27h linear, non-adiabatic mode,//,9h omega = ,
     $ 1pe16.8,9h sec(-1),,2x,9hperiod = ,e15.7,6h days,,3x,
     $ 34hfraction energy gain per period = ,e14.6,//,
     $ 3x,1hi,10x,2hdx,20x,2hdl,18x,3hdro,18x,2hdt,19x,2hdp,17x,4hwork,
     $ /,4x,5(5x,3hamp,6x,3hpha,4x),3x,4hzone,4x,8hintegral,/,
     $ (1x,i3,6(1x,2e10.3)) )
 3200 format(/,25h dimensionless frequency=,1pe12.4,
     $ /,17h q-value in days=,e12.5)
 3300 format(1h1,28h cartesian form of dx and dl,/,
     $ 3x,1hi,9x,2hdx,17x,6hdlther,15x,6hdlmech,17x,4hdl/l,17x,
     $ 4hdrho,16x,5hds/cv,/,4x,6(3x,4hreal,6x,5himag.,3x),/,
     $  (1x,i3,6(1x,1p,2e10.3)) )
 3600 format(/,15h kinetic energy,1pe16.8,7h (ergs))
 1300 format(/,43h can only converge on non-adiabatic period ,i2,3h to,
     $ 1pe10.3,7h accur=,e10.3,20h after 20 iterations)
c
 7000 format(1h1,/1x,43hpulsation matrices: ag1(3), ag2(2), ak1(4),,
     $              12h and ak2(3).)
 7001 format(1x,i4,1x,1p,3e10.3,1x,2e10.3,1x,4e10.3,1x,3e10.3)
 7800 format(8h iterad=,i3,4h iv=,i3,4h ih=,i3,6h omsq=,1pe12.4,
     $   7h omsqc=,e12.4,7h omsqp=,e12.4)
 8000 format(40h cant find no. of nodes.eq.ih-1...nodes=,i3,4h ih=,i3,
     $ 7h omega=,1pe13.6,10h period = ,e13.6)
 8600 format(14h adiabatic f2=,1pe17.10,18h is .lt. 0 for iv=,
     $  i3,4h ih=,i3)
 8100 format(14h adiabatic f2=,1pe15.6,20h not converged after,i4,
     $   16h iterations, ih=,i4,17h number of nodes=,i4)
 1000 format(1h1,19h adiabatic solution,/,9h  omega =,1pe17.9,5x,
     $ 8hperiod =,0pf12.6,//,3x,1hi,7x,2hdr,12x,2hdl,11x,4hdrho,11x,
     $  2hdt,12x,2hdp,11x,4hcvtm,10x,6hweight/(1x,i4,1p,7e14.6))
 1001 format(/,22h linear adiabatic mode,4h iv=,i3,4h ih=,i3,/,
     $ 1x,8h omega =,1pe16.8,9h sec(-1),,1x,8hperiod =,e12.5,
     $ 6h days,,1x,9hq-value =,e10.3,/,
     $ 1x,38h eigenfrequency from matrix solution =,e12.4,/,
     $ 1x,38h eigenfrequency from weight function =,e12.4,
     $ 7h error=,e12.4)
 7900 format(1h1,//,20x,28h***** begin search for mode ,i2,6h *****)
      end
      subroutine amppha(ca,amp,phase)
      implicit real(a-h,o-z)
      complex ca
c
c          returns the magnitude of ca in amp and the phase
c       of ca in phase.
c
      amp = cabs(ca)
      phase = atan2(aimag(ca),real(ca))
      return
      end
      subroutine zie(t,n,nzmax)
      implicit real(a-h,o-z)
c
c          find the maxima and minima in the work per zone function.
c
      dimension t(nzmax),imax(2001),imin(2001)
c
      nmax = 0
      nmin = 0
      ido = n - 2
      do 30 i=3,ido
         max = 0
         if( t(i) .gt. t(i-2) ) max = max + 1
         if( t(i) .gt. t(i-1) ) max = max + 1
         if( t(i) .gt. t(i+1) ) max = max + 1
         if( t(i) .gt. t(i+2) ) max = max + 1
         if( max .eq. 0 ) goto 20
         if( max .lt. 4 ) goto 30
c
c          record a maximum
c
         nmax = nmax + 1
         if( nmax .gt. 100 ) goto 30
         imax(nmax) = i
         goto 30
c
c          record a minimum
c
  20  continue
         nmin = nmin + 1
         if( nmin .gt. 100 ) goto 30
         imin(nmin) = i
  30  continue
      write(11,3000) nmax
      if( nmax .le. 0 ) goto 40
      write(11,3001) (imax(i),t(imax(i)),i=1,nmax)
  40  continue
      write(11,4000) nmin
      if( nmin .le. 0 ) return
      write(11,3001) (imin(i),t(imin(i)),i=1,nmin)
      return
c
 3000 format(/,1x,i4,40h maxima found in work per zone function ,
     $  42hwhere zone(i) .gt. (i-2),(i-1),(i+1),(i+2),
     $  30h...listed below by zone, value)
 4000 format(/,1x,i4,40h minima found in work per zone function ,
     $  42hwhere zone(i) .le. (i-2),(i-1),(i+1),(i+2),
     $  30h...listed below by zone, value)
 3001 format(8(1x,i4,1p,e10.3))
      end
      subroutine workmx(w,nz,nzmax)
      implicit real(a-h,o-z)
c
c          find all maxima in normalized work function.
c
      dimension w(nzmax),imax(2001)
      nmax = 0
      ido = nz - 1
      do 10 i=2,ido
         if( w(i).le.w(i-1) .or. w(i).le.w(i+1) ) goto 10
c
c          record maximum
c
         nmax = nmax + 1
         imax(nmax) = i
  10  continue
      write(11,1000) nmax
      if( nmax .le. 0 ) return
      write(11,1001) (imax(i),w(imax(i)),i=1,nmax)
      return
c
 1000 format(/,1x,i4,41h maxima found in normalized work function,
     $  30h...listed below by zone, value)
 1001 format(8(2x,i4,1pe10.3))
      end
      subroutine cptdmp(vec,nmax,n,ititl)
      implicit real(a-h,o-z)
      complex vec
c
c         write the vector vec to the plot file (file12) with appended
c       title ititl. In the version the real and imaginary parts are
c       dumped separately.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(real(vec(i)),i=1,n)
      write(12,1001) n,ititl,(aimag(vec(i)),i=1,n)
      return
c
 1000 format(i4,10x,a4,1x,4hreal,/,(1p,6e12.4) )
 1001 format(i4,10x,a4,1x,4himag,/,(1p,6e12.4) )
      end
