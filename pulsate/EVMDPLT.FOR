From:	ESSDP2::CCVAX2::090245       "David Kratzer, ICN Consulting Office" 31-AUG-1990 13:12:02.01
To:	ESSDP2::DE614::PESNELL
CC:	
Subj:	EVMDPLTS.FOR;1

      program evmnplt
      implicit real*8(a-h,o-z)
      real*8 n2sc,n2le,n2
      character*40 infile
      integer*2 lfl
c
      dimension r(700),rm(700),rlum(700),t(700),rho(700),p(700),
     $  akap(700),cv(700),chr(700),cht(700),epsr(700),epst(700),
     $  dkdt(700),del(700),delad(700),rmu(700),
     $  eps(700),dkdr(700)
      dimension xx(19,700),tthl(700),g1(700),sl(700,3)
      dimension nip(700)
      dimension xi(700),u(700),v(700),voga1(700),ra(700),
     $  c1(700),c4(700),b1(700),b2(700),b3(700),b4(700),
     $  dpdrho(700),derdel(700),eta(700),xhe(700),derro(700),
     $  chtor(700),alfa(700),opv(700),gor(700),qg(700),
     $   n2sc(700),n2le(700),ratn2(700),n2(700)
c
      open (unit=5,file='sys$input',status='old')
      open (unit=6,file='sys$output',status='old')
      open (unit=7,file='sumout',status='new')
c
      istat = lib$get_input(infile,'enter input file:',lfl)
      if( .not. istat ) call lib$signal(%val(istat))
c
      open (unit=20,file=infile(1:lfl),defaultfile='evol_mod',
     $      status='old')
c
      nmax = (700)
      read(20,*) qmass,dum2,dum3,rfit,qlsol,qrsol,teff,qtau,qdm1np
      read(20,1000) np
 1000 format(i5)
      do 50 i=1,18
         in=i
         if(i.ge.8) in=i+1
         if(i.eq.17) in=8
         if(i.eq.18) in=i
         read(20,1001) (xx(in,n),n=1,np)
   50 continue
      do 51 n=1,np
         rho(n) = xx(5,n)
         r(n) = xx(1,n)
  51  continue
      call dert(rho,r,derro,np)
      do 52 i=1,np
         xx(19,i)=derro(i)
  52  continue
 1001 format(4e22.15)
c
c***** take first shell in convection zone as surface point
c     np=447
c
      write(6,5200)
      read(5,*) ist
      if( ist .eq. -1 ) goto 120
      if( ist .eq. 0 ) goto 64
c**** select shell spacing so delta xi approx. 0.20 ****
c
      do 60 i=2,np
         xi(i) = dlog(xx(1,i)/xx(6,i))
  60  continue
      xi(1) = xi(2)-0.18d0
      write(6,6000)
      read(5,*) dxmin
      if( dxmin .le. 0.0d0 ) dxmin = 0.2d0
c
c          ip is the zone number in the new model.
c          n is the pointer to the zone in the original model.
c          n1 is the running pointer to the next zone, when
c       x(n1) - x(n) is greater than dxmin, ip is increased and
c       the position of the next point stored in nip(ip).
c
      n = 1
      nip(1) = 1
      n1 = 1
      ip = 1
   61 continue
      n1 = n1+1
      dxi = xi(n1)-xi(n)
      if( n1 .eq. np ) goto 62
      if( dxi .lt. dxmin ) goto 61
  62  continue
      ip = ip+1
      dx2 = dxi - dxmin
      n1m = n1-1
      dx1 = dxmin - (xi(n1m)-xi(n))
      if( n1m .eq. n ) dx1 = 1.0d0
      n = n1
      if( dx2 .gt. dx1 ) n = n1m
      n1 = n
      nip(ip) = n
      if( n1m+1 .lt. np ) goto 61
      np = ip
      do 63 i=1,19
         do 63 n=1,np
            nn = nip(n)
            xx(i,n) = xx(i,nn)
  63  continue
c
c******************************************************
c
   64 continue
c
      iso=0
      if(iso.eq.0) goto 68
c**** select original shell spacing in evol. model ****
c
      do 67 j=1,19
      i=0
      do 65 in=1,190,10
      i=i+1
   65 xx(j,i)=xx(j,in)
      do 66 in=191,536,3
      i=i+1
   66 xx(j,i)=xx(j,in)
      do 69 in=537,584
      i=i+1
   69 xx(j,i)=xx(j,in)
   67 continue
      np=i
c
c******************************************************
c
   68 continue
c
 120  continue
      xx(2,np-1) = xx(2,np-2)
      xx(2,np) = xx(2,np-1)
c
      open (unit=18,file='mon_mod',status='new',form='unformatted')
c
      write(18) np,qmass,qlsol,teff
      write(6,1200) np,qmass,qlsol,teff
      write(7,1200) np,qmass,qlsol,teff
      do 100 i=1,19
         write(18) (xx(i,n),n=1,np)
 100  continue
      rewind 18
      write(6,2005) teff,np
      write(6,2006)
      read(5,*) ist
      if( ist .le. 0 ) stop
c
      open (unit=30,file='evol_plt',status='new')
      read(18) np
      write(30,1201) np,qlsol,teff
      read(18) (    r(n),n=1,np)
      call pltdmp(r,nmax,np,4hr   ,2)
      read(18) (   rm(n),n=1,np)
      call pltdmp(rm,nmax,np,4hrm  ,4)
      read(18) (   rlum(n),n=1,np)
      call pltdmp(rlum,nmax,np,4hrl  ,2)
      read(18) (    t(n),n=1,np)
      call pltdmp(t,nmax,np,4ht   ,3)
      read(18) (  rho(n),n=1,np)
      call pltdmp(rho,nmax,np,4hrho ,3)
      read(18) (    p(n),n=1,np)
      call pltdmp(p,nmax,np,4hp   ,3)
      read(18) (  eps(n),n=1,np)
      call pltdmp(eps,nmax,np,4heps ,1)
      read(18) ( akap(n),n=1,np)
      call pltdmp(akap,nmax,np,4hopac,3)
      read(18) (   cv(n),n=1,np)
      call pltdmp(cv,nmax,np,4hcv  ,1)
      read(18) (  chr(n),n=1,np)
      call pltdmp(chr,nmax,np,4hchr ,1)
      read(18) (  cht(n),n=1,np)
      call pltdmp(cht,nmax,np,4hchit,1)
      read(18) ( epsr(n),n=1,np)
      call pltdmp(epsr,nmax,np,4hepsr,1)
      read(18) ( epst(n),n=1,np)
      call pltdmp(epst,nmax,np,4hepst,1)
      read(18) ( dkdr(n),n=1,np)
      call pltdmp(dkdr,nmax,np,4hdkdr,1)
      read(18) ( dkdt(n),n=1,np)
      call pltdmp(dkdt,nmax,np,4hdkdt,1)
      read(18) (  del(n),n=1,np)
      call pltdmp(del,nmax,np,4hdel ,1)
      read(18) (delad(n),n=1,np)
      call pltdmp(delad,nmax,np,4hdlad,1)
      read(18) (xhe(n),n=1,np)
      call pltdmp(xhe,nmax,np,4hxhe ,1)
      read(18) (derro(n),n=1,np)
      call pltdmp(derro,nmax,np,4hdero,1)
      close (unit=18)
c
      open (unit=11,file='modelk_out',status='new')
c
      write(11,2000) teff
      write(11,2001) np
      write(11,2002)
      do 10 i=1,np
         write(11,2003) i,r(i),rm(i),rlum(i),t(i),rho(i),p(i),eps(i),
     $      akap(i),cv(i)
  10  continue
      write(11,2004)
      do 20 i=1,np
         write(11,2003) i,chr(i),cht(i),epsr(i),epst(i),dkdr(i),dkdt(i),
     $      del(i),delad(i),xhe(i)
  20  continue
      g = 6.6720d-8
      pi = 3.141593d0
      pi4 = 4.0d0*pi
      a = 7.56566d-15
      c = 2.997925d10
c
c          compute logarithmic derivatives
c
      call dert(p,rho,dpdrho,np)
      call pltdmp(dpdrho,nmax,np,4hdpdr,1)
      call dert(delad,r,derdel,np)
      call dert(rlum,r,eta,np)
c
      np = np-1
      write(6,4005) np
      read(5,*) imin
      if( imin .le. 0 ) imin = nmax + 1
      totr = r(np)
      totm = rm(np)
      gorsr = dsqrt(g*totm/totr)/totr
c
c          calculate the input quantities for the eulerian pulsation
c       codes. the important quantity is ra, the convective stablity
c       coefficient. for now, when i .lt. imin ra is purely radiative,
c       when i .ge. imin ra is given by the mesh derivatives.
c
      do 40 i=2,np
         xi(i) = dlog(r(i)/p(i))
         chtor(i) = cht(i)/chr(i)
         v(i) = g*rho(i)*rm(i)/p(i)/r(i)
         opv(i) = 1.0d0/(1.0d0 + v(i))
         qg(i) = g*rm(i)/r(i)**2
         gor(i) = qg(i)/r(i)
         g1(i) = chr(i)+p(i)*cht(i)**2/(rho(i)*t(i)*cv(i))
         g3m10 = p(i)*cht(i)/(rho(i)*t(i)*cv(i))
         delad0 = g3m10/g1(i)
         qch = dabs( (delad(i)-delad0)/delad(i) )
         write(7,5311) i,delad0,delad(i),qch
         voga1(i) = v(i)/g1(i)
         u(i) = pi4*rho(i)*r(i)**3/rm(i)
         ra(i) = v(i)*chtor(i)*(del(i)-delad(i))
         if( (i.lt.imin) .or. (ra(i).gt.0.0d0) ) goto 53
         ra(i) = (derro(i)+voga1(i))
   53   continue
c
c          calculate the two versions of ra, schwarzschild and ledoux,
c       and write thier ratio to the plot file.
c
         n2sc(i) = v(i)*chtor(i)*(del(i)-delad(i))
         n2le(i) = (derro(i)+voga1(i))
         ratn2(i) = dlog10(dabs(n2le(i)/ra(i)))
         n2sc(i) = dlog10(dabs(gor(i)*n2sc(i)))
         n2(i) = dlog10(dabs(gor(i)*ra(i)))
         n2le(i) = dlog10(dabs(gor(i)*n2le(i)))
c
         cp = cv(i)*g1(i)/chr(i)
         c4(i) = gorsr*cp/(a*c*akap(i)*t(i)**3)
         alfa(i) = delad(i)/del(i)
         bfac = pi4*r(i)**3*rho(i)/rlum(i)
         bfacep = bfac*eps(i)
         if( ra(i) .lt. 0.0d0 ) eta(i) = bfacep
         b1(i) = (delad(i)*(dkdt(i)-4.0d0)+dkdr(i)/g1(i))*v(i) +
     $           (v(i)+derdel(i))*alfa(i)
         b2(i) = (epsr(i)*voga1(i)+epst(i)*v(i)*delad(i))*bfacep
     $         - (eta(i)-bfacep)*voga1(i)
         b3(i) = (epst(i)-epsr(i)*chtor(i))*bfacep +
     $           (eta(i)-bfacep)*chtor(i)
         b4(i) = bfac*gorsr*t(i)*cp
         r(i) = r(i)/totr
         rm(i) = rm(i)/totm
         c1(i) = r(i)**3/rm(i)
   40 continue
      call pltdmp(n2(2),nmax,np-1,4hn2  ,1)
      call pltdmp(n2sc(2),nmax,np-1,4hn2sc,1)
      call pltdmp(n2le(2),nmax,np-1,4hn2le,1)
      call pltdmp(ratn2(2),nmax,np-1,4hrtn2,1)
      call pltdmp(g1,nmax,np,4hg1  ,1)
      call bvcomp(rho,p,g1,chtor,v,del,delad,xhe,np,nmax)
c
c       **   compute thermal timescale  **
c
      npm1 = np-1
      tthl(np) = 0.0d0
      tthl(1) = 0.0d0
      tth = 0.0d0
      alstar = rlum(np)
      write(6,4200) np,alstar
      do 41 i=2,npm1
         n=np-i+1
         dxi=xi(n+1)-xi(n)
         dtn = pi4*cv(n)*t(n)*rho(n)*(r(n)*totr)**3
         dtd = alstar*(1.0d0+v(n))
         dtth = dtn/dtd
         tth = tth+(dtth*dxi)
         tthl(n) = dlog10(tth)
   41 continue
      write(7,4200) (i,rlum(i),i=2,npm1)
c**********************************
      write(11,3000)
      do 80 n=2,np
         write(11,3001) n,r(n),rm(n),xi(n),u(n),v(n),voga1(n),ra(n),
     $      chtor(n),c1(n)
  80  continue
      write(11,3002)
      do 90 n=2,np
         write(11,3001) n,c4(n),b1(n),b2(n),b3(n),b4(n),alfa(n),tthl(n),
     $      derdel(n),eta(n)
  90  continue
      write(11,3003)
      write(11,3004)
c
      open (unit=19,file='puls_mod',status='new',form='unformatted')
c
      npm1 = np-1
      write (19) npm1,totr,totm
      write (19) (    r(n),n=2,np)
      write (19) (   xi(n),n=2,np)
      write (19) (    u(n),n=2,np)
      write (19) (    v(n),n=2,np)
      write (19) (voga1(n),n=2,np)
      write (19) (   ra(n),n=2,np)
      write (19) (   c1(n),n=2,np)
      write (19) ( rho(n),n=2,np)
      write (19) ( delad(n),n=2,np)
      write (19) ( tthl(n),n=2,np)
      write (19) (chtor(n),n=2,np)
      write (19) (  eta(n),n=2,np)
      write (19) (   c4(n),n=2,np)
      write (19) (   b1(n),n=2,np)
      write (19) (   b2(n),n=2,np)
      write (19) (   b3(n),n=2,np)
      write (19) (   b4(n),n=2,np)
      write (19) (  del(n),n=2,np)
      write (19) ( dkdr(n),n=2,np)
      write (19) ( dkdt(n),n=2,np)
      write (19) ( alfa(n),n=2,np)
      rewind 19
c
      nsurf = np - 1
      write(7,4200) (i,tthl(i),i=1,nsurf)
c
c           write out file for carl's code.
c
      open (unit=40,file='cjh_in',status='new')
      write(40,4000) qlsol,imin
      write(40,4001) qmass
      write(40,4002) nsurf
      write(40,4003) (xi(i),r(i),qg(i),rho(i),i=2,np)
      write(40,4004) nsurf
      do 110 i=2,np
         ra(i) =-ra(i)
         write(40,4003) xi(i),gor(i),voga1(i),ra(i)
         write(40,4003) u(i),opv(i)
 110  continue
      nte = nsurf
      call pltdmp(xi(2),nmax,nte,4hxi  ,1)
      call pltdmp(u(2),nmax,nte,4hu   ,1)
      call pltdmp(voga1(2),nmax,nte,4hv   ,1)
      call pltdmp(gor(2),nmax,nte,4hg   ,1)
      call pltdmp(ra(2),nmax,nte,4hra  ,1)
c
      nte = nte - 1
      do 70 i=2,nte
         r2 = (r(i)*totr)**2
         ra(i) = dlog10(dabs(gor(i)*ra(i)))
         sl(i,1) = dlog10(2.0d0*p(i)*g1(i)/rho(i)/r2)
         sl(i,2) = dlog10(6.0d0*p(i)*g1(i)/rho(i)/r2)
         sl(i,3) = dlog10(12.0d0*p(i)*g1(i)/rho(i)/r2)
  70  continue
      nte = nte - 1
      call pltdmp(ra(2),nmax,nte,4hbv  ,1)
      call pltdmp(sl(2,1),nmax,nte,4hsl1 ,1)
      call pltdmp(sl(2,2),nmax,nte,4hsl2 ,1)
      call pltdmp(sl(2,3),nmax,nte,4hsl3 ,1)
      stop
 1200 format(x,3hnp=,i4,8h m/msun=,f7.2,8h l/lsun=,1pe10.3,
     $      6h teff=,e10.3)
 1201 format(x,3hnp=,i4,8h l/lsun=,1pe10.3,6h teff=,e10.3)
c
 2000 format(1h1,50x,36hmodel data for mon 60510t ml1 teff =,f8.1)
 2001 format(64x,i3,7h points)
 2002 format(//,1hn,7x,1hr,13x,2hrm,12x,2hlr,12x,1ht,13x,3hrho,11x,
     $  1hp,13x,3heps,11x,3hkap,11x,2hcv/)
 2003 format(i4,1p9e14.6)
 2004 format(2h- ,1hn,7x,3hchr,11x,3hcht,11x,4hepsr,10x,4hepst,10x,
     $  4hdkdr,10x,4hdkdt,10x,3hdel,11x,5hdelad,9x,3hxhe/)
 2005 format(x,36hmodel data for mon 60510t ml1 teff =,f8.1,/,
     $   x,i4,7h points)
 2006 format(x,41hdo you want the model dumps? (1=yes,0=no) )
c
 3000 format(1h1,59x,19hcomputed quantities/,45x,
     $   47h(defined in saio & cox, ap. j. 236, 558 (1980)),//,
     $   52x,33hcenter and surface points deleted,//,
     $   3h n ,6x,2hfr,12x,3hfrm,11x,7hln(r/p),7x,1hu,13x,
     $   1hv,13x,5hvoga1,9x,2hra,12x,5hchtor,9x,2hc1/)
 3001 format(i4,1p9e14.6)
 3002 format(1h-,3h n ,6x,2hc4,12x,2hb1,12x,2hb2,12x,2hb3,12x,2hb4,12x,
     $  4halfa,10x,6h tthl ,8x,8hderdelad,6x,3heta/)
 3003 format(129h0* a discontinuity in delmu may occur about the point w
     $here the interior and envelope are joined.  this is a defect in th
     $e model,/44h   but is included in the calculation of ra.)
 3004 format('0  the jump in delmu starting around point 145 is due to i
     $onization, and is not included in the calculation of ra.')
c
 4000 format(1pe12.5,i5)
 4001 format(1pe12.4)
 4002 format(i4)
 4003 format(1p4e20.12)
 4004 format(i5)
 4005 format(x,9hmodel has,i4,20h points. enter imin.)
 4200 format( 4(2x,i4,1pe10.3) )
c
 5200 format(x,34henter control flag -1, full model.)
 5311 format(x,3hi =,i4,10h delad0 = ,1pe10.3,9h delad = ,e10.3,
     $         8h diff = ,e10.3)
 6000 format(x,34henter dxmin ( .le. 0 dxmin = 0.2.) )
c
      end
      subroutine dert(y,h,der,np)
      implicit real*8(a-h,o-z)
c
      dimension y(700),h(700),der(700)
c
c          program allows y or h = 0 for n = 1 only
c
      if (y(1) .eq. 0.) goto 5
      if (h(1) .eq. 0.) goto 5
      ym = dlog(y(2)/y(1))
      if( dabs(ym) .lt. 1.d-6) ym = 0.0d0
      hm = dlog(h(2)/h(1))
    5 npm2=np-2
      do 10 n=1,npm2
         yp = dlog(y(n+2)/y(n+1))
         if( dabs(yp) .lt. 1.d-6 ) yp = 0.0d0
         hp = dlog(h(n+2)/h(n+1))
         if( y(n) .eq. 0.0d0) goto 9
         if( h(n) .eq. 0.0d0) goto 9
         htot = dlog(h(n+2)/h(n))
         der(n+1) = (hm*yp/hp+hp*ym/hm)/htot
         ym=yp
         hm=hp
         goto 10
    9    der(2)=yp/hp
         ym=yp
         hm=hp
   10 continue
      der(np) = ym/hm
      return
      end
      subroutine pltdmp(vec,nmax,n,ititl,ityp)
      implicit real*8(a-h,o-z)
c
      dimension vec(nmax),temp(700)
      goto (10,20,30,40),ityp
  10  continue
c
c         write vector as inputted.
c
      write(30,1000) n,ititl,(vec(i),i=1,n)
      return
  20  continue
c
c         write vector normaized to outermost point.
c
      do 25 i=1,n
         temp(i) = vec(i)/vec(n)
  25  continue
      write(30,1000) n,ititl,(temp(i),i=1,n)
      return
  30  continue
c
c         write logarithm of vector.
c
      do 35 i=1,n
         temp(i) = dlog10(vec(i))
  35  continue
      write(30,1000) n,ititl,(temp(i),i=1,n)
      return
  40  continue
c
c         write logarithm of vector.
c
      sfact = 1.d0 - 1.d-15
      do 45 i=1,n
         temp(i) =-dlog10( 1.d0 - sfact*(vec(i)/vec(n)) )
  45  continue
      write(30,1000) n,ititl,(temp(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p6e12.4) )
      end
      subroutine bvcomp(rho,p,g1,chtor,v,del,delad,xhe,np,nmax)
      implicit real*8(a-h,o-z)
      dimension rho(nmax),p(nmax),g1(nmax),chtor(nmax),v(nmax),
     $      del(nmax),delad(nmax),xhe(nmax)
c
c          compare the various incarnations of the brunt-vaisala
c       frequencies.
c
      dimension plog(700),bvrad(700),bvnum(700),bvcor(700),rhol(700)
      dimension qmul(700)
c
      n = np
      do 10 i=2,n
         j = i
         if( (xhe(i).gt.xhe(i-1)) .and. (xhe(i).gt.xhe(i+1)) ) goto 15
  10  continue
  15  continue
      nmax = j
      do 20 i=2,nmax
         plog(i) = dlog10(p(i))
         rhol(i) = dlog10(rho(i))
         qmul(i) = dlog(6.0d0/(0.5d0 + xhe(i)))
         bvrad(i) = v(i)*chtor(i)*(del(i)-delad(i))
  20  continue
c
c          set up for numeric derivatives.
c
      do 30 i=3,nmax-2
         dlpi = plog(i) - plog(i-1)
         dlpi1 = plog(i+1)-plog(i)
         ddlp1 =-dlpi1/(dlpi*(dlpi+dlpi1))
         ddlp2 = (dlpi1-dlpi)/(dlpi*dlpi1)
         ddlp3 = dlpi/(dlpi1*(dlpi+dlpi1))
         dlrdlp = ddlp1*rhol(i-1) + ddlp2*rhol(i) + ddlp3*rhol(i+1)
         dlmdlp = ddlp1*qmul(i-1) + ddlp2*qmul(i) + ddlp3*qmul(i+1)
         bvnum(i) = (-dlrdlp + 1.0d0/g1(i))*v(i)
         bvcor(i) = v(i)*( chtor(i)*(del(i)-delad(i)) -
     $      dlmdlp*(xhe(i)+0.5d0)/3.d0 )
         write(7,3000) i,bvrad(i),bvnum(i),bvcor(i)
  30  continue
      call pltdmp(plog(3),nmax,nmax-2,4hplog,1)
      call pltdmp(bvrad(3),nmax,nmax-2,4hbvrd,1)
      call pltdmp(bvnum(3),nmax,nmax-2,4hbvnm,1)
      call pltdmp(bvcor(3),nmax,nmax-2,4hbvcr,1)
      call pltdmp(xhe(3),nmax,nmax-2,4hxheb,1)
      call pltdmp(v(3),nmax,nmax-2,4hvbv ,1)
      call pltdmp(qmul(3),nmax,nmax-2,4hmu  ,1)
      return
 3000 format(x,i4,1p3e12.4)
      end
