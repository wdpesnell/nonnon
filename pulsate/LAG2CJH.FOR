*byp3
*rcft i=%me,x=xcon,on=o,cs=170140
      program lag2cjh(tty,tape5=tty,tape6=tty,lnrad,tape10=lnrad,
     $   conout,tape11=conout,conplt,tape12=conplt,
     $   sumout,tape1=sumout,rkdmp,tape71=rkdmp)
      implicit real(a-h,o-z)
c
c
c                                                12/20/83 wd pesnell
c
      common/phypar/ rp(2001),tp(2001),vp(2001),cv(2001),dkdr(2001),
     $   dkdt(2001),dm1(2001),akap(2001),dm2(2001),rm(2001),bv(2001)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(2001),g1(2001),g3m1(2001),rho(2001),rzone(2001)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/lumins/ frft(2001),sorce(2001),dtsorc(2001),dvsorc(2001)
      common/stnmod/ nmode
      dimension ititl(20)
c
c          read in the input variables from file 10.
c
c     open (unit=1,file='sumout',status='unknown',acess='append')
c      open (unit=1,file='sumout',status='unknown')
c     open (unit=5,file='sys$input',status='unknown')
c     open (unit=6,file='sys$output',status='unknown')
c      open (unit=10,file='lnrad',status='old')
c      open (unit=11,file='nonout',status='new')
c
      call filerep
c
      nmax = (2001)
      read(10,1001,end=100) (ititl(i),i=1,20)
c
c          read in the input variables from file 10.
c
      read(10,*) np1,rlumgv,totmas,teff,rphoto,corlum
      read(10,*) irad,noburn,onemq0,onemq1,rlums,rlumc
      npts = np1 - 1
      read(10,1000) rdum,(rp(i),i=1,np1)
      read(10,1000) tc,(tp(i),i=1,npts)
      read(10,1000) vp1,(vp(i),i=1,npts)
      rhoc = one/vp1
      read(10,1000) cv0,(cv(i),i=1,npts)
      read(10,1000) dkdr0,(dkdr(i),i=1,npts)
      read(10,1000) dkdt0,(dkdt(i),i=1,npts)
      read(10,1000) cormas,(dm1(i),i=1,npts)
      read(10,1000) opac0,(akap(i),i=1,npts)
      read(10,1000) ddum,(dm2(i),i=1,np1)
      read(10,1000) rdum,(rm(i),i=1,np1)
      read(10,1000) pc,(p(i),i=1,npts)
      read(10,1000) g10,(g1(i),i=1,npts)
      read(10,1000) g3m10,(g3m1(i),i=1,npts)
      read(10,1000) fdum,(frft(i),i=1,np1)
      read(10,1000) sorc0,(sorce(i),i=1,npts)
      read(10,1000) dedt0,(dtsorc(i),i=1,npts)
      read(10,1000) dedv0,(dvsorc(i),i=1,npts)
c      read(10,1000,end=900) pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,
c     $   g10,g3m10,cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      close (unit=10)
      totmas = rm(npts+1)
c
c          output the initial model data
c
      write(11,2000) (ititl(i),i=1,20),npts
      write(11,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(11,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
      write(11,2004)
      write(11,2005) i,rp(1),tc,1.e0/rhoc,pc,g10,g3m10,opac0,
     $      dkdr0,dkdt0,cv0,frft(1),rm(1)
      write(11,2005) (i,rp(i+1),tp(i),vp(i),p(i),g1(i),g3m1(i),akap(i),
     $      dkdr(i),dkdt(i),cv(i),frft(i+1),rm(i+1),i=1,npts)
c
c          check for nuclear burning and print out those
c       zones participating.
c
      if( tp(1) .lt. 3.d6 ) goto 20
      do 10 i=1,npts
         imax = np1-i
         if( sorce(imax) .ne. zero ) goto 15
  10  continue
      write(11,2100)
      goto 20
  15  continue
      do 16 i=1,imax
         imin = i
         if( sorce(i) .ne. zero ) goto 17
  16  continue
      write(11,2100)
      goto 20
  17  continue
      write(11,2102) (i,sorce(i),dtsorc(i),dvsorc(i),i=imin,imax)
      write(11,2101)
  20  continue
      call pltint(nmax,npts,ititl)
c
c
c          write the file for the Runga-Kuuta analysis.
c
      call cjhdmp(bv,npts)
      call pltdmp(bv,nmax,npts+1,4hbv  )
      stop
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
 1001 format(20a4)
 1002 format(1x,'End-of-File on title, try again.')
 1003 format(1x,27henter pulsation parameters:,/,
     $   1x,38h io, l, nom, ihmax, permin, and permx.)
c
 2000 format(50h1 initial model in hydrostatic and thermal balance,//,
     $      20a4,1x,i5,6h zones,/)
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
     $     6.02217e23,2.52157e-15 /
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726e-8,7.5595e-5,
     $  5.66961e-5,3.1415926536e0,6.2831853072e0,12.5663706144e0,
     $   25.1327412288e0,4.1887902048e0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0e0,1.0e0,2.0e0,3.0e0,4.0e0,10.0e0,
     $      0.5e0,0.25e0      /
      end
      subroutine pltint(nmax,n,ititl)
      implicit real(a-h,o-z)
c
c          initialize the plot file and write the title line.
c
      common/phypar/ r(2001),t(2001),v(2001),cv(2001),dkdr(2001),
     $   dkdt(2001),dm1(2001),akap(2001),dm2(2001),rm(2001),bv(2001)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(2001),g1(2001),g3m1(2001),omq(2001),tlog(2001)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      dimension plog(2001),xrad(2001),zoneb(6003)
      dimension ititl(20)
c
c      open (unit=12,file='nonplt',status='new')
      write(12,1000) ititl
      np1 = n + 1
      call pltdmp(r,nmax,n+1,4hr   )
      sfact = one - 5.e-14
      do 10 i=1,n
         omq(i) =-alog10(one - sfact*rm(i)/totmas)
         tlog(i) = alog10(t(i))
         plog(i) = alog10(p(i))
         xrad(i) = r(i)/r(np1)
  10  continue
      xrad(np1) = one
      omq(n+1) =-alog10(one-sfact)
      call pltdmp(tlog,nmax,n,4ht   )
      call pltdmp(plog,nmax,n,4hp   )
      call pltdmp(omq,nmax,n+1,4h1-q )
      call pltdmp(xrad,nmax,n+1,4hx   )
      call pltdmp(g1,nmax,n,4hg1   )
      call pltdmp(g3m1,nmax,n,4hg3m1)
      do 20 i=1,n
         xrad(i) = float(i)
  20  continue
      call pltdmp(xrad,nmax,n,4hzone)
      do 30 i=1,3*(n+1)
         zoneb(i) = float(i)
  30  continue
c      call pltdmp(zoneb,nmax,3*(n+1),4hbzon)
      return
 1000 format(20a4)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real(a-h,o-z)
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p,6e12.4) )
      end
      function knode(x,y,nmax,n,nps,ngs)
      implicit real(a-h,o-z)
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
         if( x(i+1) .lt. x(i) ) goto 15
          if( y(i) .ge. zero ) ng = ng + 1
          if( y(i) .lt. zero ) np = np + 1
         goto 10
  15     if( y(i+1) .ge. zero ) np = np + 1
         if( y(i+1) .lt. zero ) ng = ng + 1
  10  continue
      if( x(n-1)*x(n) .lt. zero ) write(11,1001)
      knode = np - ng
      nps = np
      ngs = ng
      return
 1000 format(1x,50hfrom knode...quadrant jump default rotations used.)
 1001 format(1x,51hfrom knode...node at outer boundary is not counted.)
      end
      subroutine cjhdmp(bv,ir)
      implicit real(a-h,o-z)
c
c          write out the input file for carl hansen's runge-kutta
c       pulsation analysis.
c
      common/phypar/ r(2001),tp(2001),vp(2001),cv(2001),dkdr(2001),
     $      dkdt(2001),dm1(2001),akap(2001),dm2(2001),rm(2001)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(2001),g1(2001),g3m1(2001),entpry(2001)
      common/blk8/   g,ac3,sigma,pi,pi2,forpi,pi8,pi43
      common/coretc/ pc,rhoc,tc,pc0,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/scrtch/ chr(2001),chit(2001),cp(2001),del(2001),
     $      delad(2001),q(2001),cvel(2001),fconv(2001),spac(2001,5)
      dimension xi(2001),rnorm(2001),qg(2001),vog1(2001),ra(2001),
     $      u(2001),opv(2001),rho(2001),bv(2001)
      data lnafil/71/
c
c      open (unit=lnafil,file='cjhinp',status='new')
c
c      write information to file 11.
c
      qlsol = alog10(rlumgv/3.82e33)
      qmass = totmas/1.991e33
      irp1 = ir+1
      write(lnafil,1001) qlsol,ir-1,qmass,ir-1
      p(irp1) = zero
      do 10 i=1,irp1
         rnorm(i) = r(i)/r(irp1)
         qg(i) = g*rm(i)/r(i)**2
  10  continue
      xi(1) = alog(r(1)/(ahf*(pc0+p(1))))
      u(1) = forpi*r(1)**3*rhoc/rm(1)
      vog1(1) = qg(1)*r(1)/(vp(1)*p(1)*g1(1))
      vi = qg(1)*r(1)/(vp(1)*p(1))
      opv(1) = 1.e0/(one+vi)
      rho(1) = ahf*(rhoc + 1.e0/vp(1))
      dlnrho =-alog(rhoc*vp(1))
      xin = alog(alog(p(1))/alog(pc))
      denom = one - alog(g1(1)/g10)/xin
      dlpdlr = alog(p(1))/g1(1) - alog(pc0)/g10
      ra(1) =-forpi*r(1)**3*rho(1)*( dlnrho - dlpdlr/denom )/dm2(1)
      bv(1) = alog10(qg(1)/r(1)*abs(ra(1)))
      do 20 i=2,ir
         vog1(i) = qg(i)*r(i)*ahf*(1.e0/(vp(i-1)*p(i-1)*g1(i-1)) +
     $                                1.e0/(vp(i)*p(i)*g1(i)) )
         rho(i) = ahf*(1.e0/vp(i-1) + 1.e0/vp(i))
         dlnrho =-alog(vp(i)/vp(i-1))
         xin = alog(alog(p(i))/alog(p(i-1)))
         denom = one - alog(g1(i)/g1(i-1))/xin
         dlpdlr = alog(p(i))/g1(i) - alog(p(i-1))/g1(i-1)
         ra(i) =-forpi*r(i)**3*rho(i)*( dlnrho - dlpdlr/denom )/dm2(i)
         bv(i) = alog10(qg(i)/r(i)*abs(ra(i)))
         vi = qg(i)*r(i)*ahf*(1.e0/(vp(i-1)*p(i-1)+1.e0/(vp(i)*p(i))) )
         opv(i) = 1.e0/(one+vi)
         xi(i) = alog(r(i)/(ahf*(p(i-1)+p(i))))
         u(i) = forpi*r(i)**3*rho(i)/rm(i)
  20  continue
      rho(irp1) = ahf/vp(ir)
      xi(irp1) = alog(r(irp1)/(ahf*p(ir)))
      u(irp1) = forpi*rho(irp1)*r(irp1)**3/rm(irp1)
      vog1(irp1) = qg(irp1)*r(irp1)/(vp(ir)*p(ir)*g1(ir))
      vi = qg(irp1)*r(irp1)/(vp(ir)*p(ir))
      opv(irp1) = 1.e0/(one+ahf*vi)
      chiti = g3m1(ir)*cv(ir)*rho(ir)*tp(ir)/p(ir)
      chri = g1(ir) - g3m1(ir)*chiti
      ra(irp1) = vi*(g3m1(ir)/g1(ir))*chiti/chri
      bv(irp1) = alog10(qg(irp1)/r(irp1)*abs(ra(irp1)))
      write(lnafil,1000) (xi(i),rnorm(i),qg(i),rho(i),i=2,ir)
      write(lnafil,1002) ir-1
      do 30 i=2,ir
         write(lnafil,1000) xi(i),qg(i)/r(i),vog1(i),ra(i)
         write(lnafil,1000) u(i),opv(i)
  30  continue
      close(unit=lnafil)
      return
 1000 format(1p,4e20.12)
 1001 format(1pe12.5,i5/1pe12.4/i4)
 1002 format(i5)
      end
