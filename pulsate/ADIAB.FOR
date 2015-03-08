From:	ESSDP2::CCVAX2::090245       "David Kratzer, ICN Consulting Office" 31-AUG-1990 13:05:17.34
To:	ESSDP2::DE614::PESNELL
CC:	
Subj:	ADIAB.FOR

      program adiab
      implicit real*8(a-h,o-z)
      character*80 idstar
c
c          the program that started it all. the eulerian adiabatic
c       pulsation code that appears in papers from the early 1970's.
c       this particular version uses the "gnr1" method for solving
c       the system of nonlinear, first order differential equations
c       as listed in osaki and hansen (ap. j. 185[277]1973) as well
c       as brad carroll's dissertation (this version is alleged to be
c       that very one) and d.e. winget's white dwarf work.
c
      common u(700),v(700),voga1(700),ra(700),c1(700)
      common/blrdg/ result(700)
      dimension xi(700),yizr(4,700),zkzr(1),ap(2),fd(4,6,700),ph(5,6),
     $  bg(5,4),b(5,6),d(5,6)
      dimension fr(700),rho(700),xint1(700),xint2(700),h(700),delad(700)
     $  ,dtot(700),stwait(700),stroot(20)
      real*8 mr,jzr
      integer s,smlp2,ap
      external co1,co2
c
c          files in this program.
c  unit   routine    name                purpose
c    5      main    sys$input           terminal input file
c    6      main    sys$output          terminal output file
c
c   11      main    adiab_out           vast amounts of output
c
c   12     pltint   adaib_plt           has vectors of pulsational
c                                          quantities for pltcmd program.
c
c   13      main    adiab_out           vast amounts of output
c
c   19      main    puls_mod            input quantities from the program
c                                          modelk.
c
      open (unit=5,file='sys$input',status='unknown')
      open (unit=6,file='sys$output',status='unknown')
      open (unit=11,file='adiab_out',status='new')
      open (unit=19,file='puls_mod',status='old',form='unformatted')
c
      nmax = (700)
      g = 6.6720d-08
      pi = 3.141592653d0
      pi2 = 2.d0*pi
      write(6,1112)
      read(5,1111,end=900) idstar
 1111 format(a80)
 1112 format(x,11henter title)
      read (19) n2,totr,totm
      n1=n2-1
      read (19) (   fr(n),n=1,n2)
      read (19) (   xi(n),n=1,n2)
      read (19) (    u(n),n=1,n2)
      read (19) (    v(n),n=1,n2)
      read (19) (voga1(n),n=1,n2)
      read (19) (   ra(n),n=1,n2)
      read (19) (   c1(n),n=1,n2)
      read (19) (  rho(j),j=1,n2)
      read (19) (delad(j),j=1,n2)
      close (unit=19)
      call pltint(nmax,fr,n2,idstar)
c
c          nint is the number of points in the eigenvalue search
c        routine, should be an even number.
c        nint .le. 0 suppresses the eigenvalue search routine
c        s is the central r-dependence of y1 used in gnr1
c        verg is the convergence criterion used in gnr1
c
   11 continue
      write(6,1000)
      read(5,*,end=900) nint,l,s,verg
      ap(1) = l
      rl1 = dfloat(l)*dfloat(l+1)
      ap(2) = s
      smlp2 = s-l+2
      uc = 1.0d0
      if( nint .le. 0 ) goto 49
c
c          permin is the shortest pulsation period (hours)
c          permax is the longest pulsation period (hours)
c
      write(6,1100)
      read(5,*,end=900) permin,permax
      if( permin .le. 0.0d0 ) stop
      omlow = pi2**2*totr**3/g/totm/(3600.d0*permax)**2
      omhigh = pi2**2*totr**3/g/totm/(3600.d0*permin)**2
      write(11,2000)
      iroot = 0
      disc1 = 0.0d0
      sint = (omhigh - omlow)/dfloat(nint-1)
      do 20 i=1,nint
         do 201 j=1,n2
            yizr(1,j) = fr(j)
            yizr(2,j) = fr(j)
            yizr(3,j) = 0.0d0
            yizr(4,j) = 0.0d0
  201    continue
         zkzr(1) = omlow + dfloat(i-1)*sint
         call gnr1(4,0,n1,4,1,5,6,700,2,xi,yizr,zkzr,ap,uc,fd,ph,
     $             bg,b,d,co1)
         disc = c1(1)*zkzr(1)*yizr(1,1)-l*yizr(2,1)
        if( disc*disc1 .gt. 0.0d0 ) goto 202
        if( i .eq. 1 ) goto 202
         iroot = iroot + 1
         stroot(iroot) = zkzr(1) - sint*disc/(disc-disc1)
 202    continue
         disc1 = disc
         smsig = sqrt(zkzr(1)*g*totm/totr**3)
         per = pi2/3600.d0/smsig
         write(11,2001) i,per,zkzr(1),disc
  20  continue
      if( iroot .eq. 0 ) goto 11
      write(6,2010) iroot,(stroot(i),i=1,iroot)
      goto 45
   49 continue
c
c          given a guess for a root, attempt to converge a solution
c        to the full adiabatic system.
c        per is the pulsation period (hours)
c
      write(6,1101)
      read(5,*,end=900) per
      if( per .le. 0.0d0 ) stop
      nroot = 1
      stroot(1) = pi2**2*totr**3/g/totm/(3600.d0*per)**2
  45  continue
      do 100 nroot=1,iroot
         zkzr(1) = stroot(nroot)
         do 40 n=1,n2
            yizr(1,n) = fr(n)
            yizr(2,n) = fr(n)
            yizr(3,n) = 0.0d0
            yizr(4,n) = 0.0d0
   40    continue
         call gnr1(4,0,n1,4,1,5,6,700,2,xi,yizr,zkzr,ap,uc,fd,ph,
     $             bg,b,d,co1)
         disc = c1(1)*zkzr(1)*yizr(1,1)-dfloat(l)*yizr(2,1)
         smsig = dsqrt(zkzr(1)*g*totm/totr**3)
         per = pi2/3600.d0/smsig
         write(11,2001) i,per,zkzr(1),disc
      do 50 n=1,2
         call gnr1(4,0,n1,4,1,5,6,700,2,xi,yizr,zkzr,ap,uc,fd,ph,
     $             bg,b,d,co1)
         disc = c1(1)*zkzr(1)*yizr(1,1)-dfloat(l)*yizr(2,1)
         write(11,5000) n,zkzr(1),disc
   50 continue
      do 60 n=1,10
         call gnr1(4,1,n1,4,1,5,6,700,2,xi,yizr,zkzr,ap,uc,fd,ph,
     $             bg,b,d,co2)
         disc = c1(1)*zkzr(1)*yizr(1,1)-dfloat(l)*yizr(2,1)
         q0 = dabs(d(5,1))
         q1 = dabs(b(1,1))
         q2 = dabs(b(2,1))
         q3 = dabs(b(3,1))
         q4 = dabs(b(4,1))
         write(11,6000) q0,q1,q2,q3,q4
         write(11,6001) (i,fr(i),yizr(1,i),yizr(2,i),yizr(3,i),
     $      yizr(4,i),i=1,10)
         if( q0.gt.verg .or. q1.gt.verg .or. q2.gt.verg .or.
     $      q3.gt.verg .or. q4.gt.verg) goto 60
         goto 70
  60  continue
      write(6,6010)
      write(11,6010)
      goto 100
c
c          converged to a root, print out eigenvectors.
c
  70  continue
      omsq = zkzr(1)*g*totm/totr**3
      smsig = dsqrt(omsq)
      per = pi2/3600.d0/smsig
      sunm = 1.989d+33
      sunr = 6.96d+10
      q = per*dsqrt(totm/sunm*(sunr/totr)**3)/24.d0
c
      if( smlp2 .eq. 0 ) goto 751
      do 75 j=1,n2
         do 75 i=1,4
         yizr(i,j) = yizr(i,j)/fr(j)**smlp2
  75  continue
 751  continue
      gmor3 = g*totm/totr**3
      pi4g = 4.d0*pi*g
      do 80 j=1,n2
         rr2 = rho(j)*(totr*fr(j))**3/(v(j)+1.d0)
         smg = gmor3*fr(j)*totr/c1(j)
         aa = fr(j)*totr*yizr(1,j)
         bb = fr(j)*totr*yizr(2,j)/c1(j)/zkzr(1)
         alf = voga1(j)*(yizr(1,j)+yizr(3,j)-yizr(2,j))
         dadr = (voga1(j)-2.)*yizr(1,j)+(rl1/c1(j)/zkzr(1)-voga1(j))
     $     *yizr(2,j)+voga1(j)*yizr(3,j)
         aa2 = aa**2
         bb2 = bb**2
         alf2 = alf**2
         xint1(j) = rr2*(aa2+bb2*rl1)
      xint2(j)=rr2*(-smg*(2.*(aa*alf-aa*dadr)-aa*yizr(4,j)-bb*yizr(3,j)
     $  *rl1-fr(j)*totr*alf2/voga1(j))+pi4g*rho(j)*aa2)
         dtot(j) = delad(j)*v(j)*(yizr(2,j)-yizr(3,j)-yizr(1,j))
         if (j .eq. 1) goto 80
         h(j) = xi(j)-xi(j-1)
  80  continue
      ititl = 4hdr/r
      write(12,8100) n2,ititl,(yizr(1,j),j=1,n2)
      ititl = 4hy2
      write(12,8100) n2,ititl,(yizr(2,j),j=1,n2)
      ititl = 4hy3
      write(12,8100) n2,ititl,(yizr(3,j),j=1,n2)
      ititl = 4hy4
      write(12,8100) n2,ititl,(yizr(4,j),j=1,n2)
      call sumit(xint2,h,n2)
      do 85 i=1,n2
         stwait(i) = result(i)
  85  continue
      twait = result(n2)
      call sumit(xint1,h,n2)
      rke = result(n2)
      do 86 i=1,n2
         stwait(i) = stwait(i)/(rke*omsq)
         result(i) = result(i)/rke
  86  continue
      call pltdmp(stwait,nmax,n2,4hwait)
      call pltdmp(result,nmax,n2,4hrke )
      qomsq = twait/rke
      qch = dabs((omsq-qomsq)/omsq)
      omint = dsqrt(dabs(qomsq))
      if( qomsq .lt. 0.d0 ) omint =-omint
      perint = pi2/3600.d0/omint
      persec = per*3600.d0
      write(11,8000) idstar
      write(11,8010) l,s
      write(11,8015) smsig,per,persec
      write(11,8020) q
      write(11,8025) perint,aa
      write(11,8030) rke,bb
      write(11,8035)
      write(11,8040) (i,fr(i),yizr(1,i),yizr(2,i),yizr(3,i),yizr(4,i),
     $  dtot(i),i=1,n2)
c
 100  continue
      goto 11
c
c          end-of-file during input -- give up.
c
 900  continue
      stop
 1000 format(x,52henter pts. in scan (.le. 0 no scan), l, s, and verg.)
 1100 format(x,51henter permin and permax in hours (per .le. 0 stop).)
 1101 format(x,43henter period in hours (per .le. 0 to stop).)
c
 2000 format(1h1,60x,17heigenvalue search)
 2001 format(x,i3,19h  period (hours) = ,1pe13.6,5x,10hzkzr(1) = ,
     $       e13.6,5x,7hdisc = ,e13.6)
 2010 format(x,i4,20h roots for analysis:,/, (3x,1p5e15.6) )
c
 5000 format(/,i2,12h  zkzr(1) = ,1pe13.6,2x,7hdisc = ,e13.6)
c
 6000 format(5x,5hq0 = ,1pe12.6,2x,5hq1 = ,e12.6,2x,5hq2 = ,e12.6,
     $      2x,5hq3 = ,e12.6,2x,5hq4 = ,e12.6)
 6001 format(x,i3,1p5e18.10)
 6010 format(/,x,13hnot converged)
c
 8000 format(1h1,59x,18hadiabatic solution,//,23x,a80)
 8010 format(/,65x,4hl = ,i2,/,66x,4hs = ,i2)
 8015 format(/,28x,32hpulsation frequency (rad/sec) = ,
     $  1pe13.6,5x,17hperiod (hours) = ,e13.6,/,e12.4)
 8020 format(/,47x,43h(zero rotation solution written on tape11),
     $   /,x,11hq (days) = ,1pe10.3,
     $   47x,35hsurface values for velocity field -)
 8025 format(/,x,28hintegrated period (hours) = ,1pe13.6,
     $     27x,23hzero rotation a (cm) = ,e13.6)
 8030 format(/,x,14hj (surface) = ,1pe13.6,
     $     41x,23hzero rotation b (cm) = ,e13.6)
 8035 format(1h-,3h n ,6x,2hfr,26x,4hy1zr,10x,4hy2zr,10x,4hy3zr,10x,
     $  4hy4zr,24x,4hdt/t,/)
 8040 format(i4,1pe14.6,14x,4e14.6,14x,e14.6)
 8100 format(i4,10x,a4,/,(1p6e12.4) )
      end
      subroutine sumit(f,h,n)
      implicit real*8(a-h,o-z)
c
      dimension f(n),h(n)
      common/blrdg/  result(700)
      nm = n-1
      result(1) = 0.0d0
      sum=h(2)*(f(2)+f(1))/2.d0
      fm=2.*((f(3)-f(2))/h(3)-(f(2)-f(1))/h(2))/(h(2)+h(3))
      fp=2.*((f(4)-f(3))/h(4)-(f(3)-f(2))/h(3))/(h(3)+h(4))
      endcor=(fm+(h(2)+h(3))*(fm-fp)/(h(2)+2.*h(3)+h(4)))*h(2)**3
      result(2)=sum-endcor/12.d0
      sum=h(3)*(f(3)+f(2))/2.d0
      dm=h(3)+h(2)
      dp=h(4)+h(3)
      sumc=(dm*fp+dp*fm)*h(3)**3/(dm+dp)
      result(3)=result(2)+sum-sumc/12.d0
      do 1 j=4,nm
         fm = fp
         dm = dp
         fp = 2.d0*((f(j+1)-f(j))/h(j+1)-(f(j)-f(j-1))/h(j))/
     $                                  (h(j)+h(j+1))
         dp = h(j+1)+h(j)
         sum = h(j)*(f(j)+f(j-1))/2.d0
         sumc = (dm*fp+dp*fm)*h(j)**3/(dm+dp)
         result(j) = result(j-1)+sum-sumc/12.d0
    1 continue
      sum = h(n)*(f(n)+f(n-1))/2.d0
      endcor=(fp+(h(n)+h(n-1))*(fp-fm)/(h(n)+2.*h(n-1)+h(n-2)))*h(n)**3
      result(n) = result(nm)+sum-endcor/12.d0
      return
      end
      subroutine gnr1(i1,k1,n1,i1m,k1m,ik1m,ik2m,n2m,np,xi,yi,zk,ap,
     $  uc,fd,ph,bg,b,d,coeffs)
      implicit real*8(a-h,o-z)
c
      dimension xi(n2m),yi(i1m,n2m),zk(k1m),ap(np),fd(i1m,ik2m,n2m),
     $  ph(ik1m,ik2m),bg(ik1m,i1m),b(ik1m,ik2m),d(ik1m,ik2m)
      integer ap
      n2=n1+1
      ik1=i1+k1
      ik2=ik1+1
      do 21 l=1,ik1
         do 21 j=1,i1
            bg(l,j)=0.d0
   21 continue
      do 25 m=1,ik2
         do 22 l=1,ik1
            ph(l,m)=0.d0
   22    continue
         do 25 i=1,i1
            do 25 n=1,n2
               fd(i,m,n)=0.d0
   25 continue
      call coeffs(n2,i1m,k1m,ik1m,ik2m,n2m,np,xi,yi,zk,ap,fd,ph,bg)
      do 38 n=1,n1
      dx=-0.5d0*(xi(n+1)-xi(n))
      do 32 i=1,i1
      do 30 j=1,i1
      b(i,j)=-dx*fd(i,j,n)
   30 continue
      b(i,i)=b(i,i)+1.d0
      do 31 m=1,ik2
      d(i,m)=dx*(fd(i,m,n+1)+fd(i,m,n))
   31 continue
      d(i,ik2)=yi(i,n+1)-yi(i,n)+d(i,ik2)
   32 continue
      call leq(b,d,i1,ik2,ik1m,ik1m,err)
      if (err .lt. 0.) write(6,1000) n
 1000 format(x,27hfrom gnr1... err .lt. 0, n=,i3)
      do 38 i=1,i1
         do 38 m=1,ik2
            fd(i,m,n)=d(i,m)
   38 continue
      do 51 l=1,ik1
      do 41 j=1,i1
      if (ph(l,j) .ne. 0.) goto 42
   41 continue
      goto 45
   42 do 44 n=1,n1
      do 43 j=1,i1
      b(j,1)=ph(l,j)
   43 continue
      do 44 m=1,ik2
      do 44 j=1,i1
      ph(l,m)=ph(l,m)+fd(j,m,n)*b(j,1)
   44 continue
   45 do 46 j=1,i1
      ph(l,j)=ph(l,j)+bg(l,j)
   46 continue
      d(l,1)=-ph(l,ik2)
   51 continue
      call leq(ph,d,ik1,1,ik1m,ik1m,err)
      if (err .lt. 0.) write(6,1000) n
      do 55 i=1,i1
         ph(i,1) = 0.0d0
         bg(i,1) = 0.0d0
   55 continue
      b(1,ik2) = 1.0d0
      do 63 nn=1,n2
      n=n2-nn+1
      do 60 i=1,i1
      fd(i,ik2,n)=d(i,1)
      yi(i,n)=yi(i,n)+d(i,1)*uc
      ph(i,1)=ph(i,1) + dabs(d(i,1))
      bg(i,1)=bg(i,1) + dabs(yi(i,n))
   60 continue
      if (n .le. 1) goto 65
      do 62 m=1,ik1
      b(1,m)=d(m,1)
   62 continue
      do 63 i=1,i1
      do 63 m=1,ik2
      d(i,1)=d(i,1)+fd(i,m,n-1)*b(1,m)
   63 continue
   65 if (k1) 68,68,66
   66 do 67 k=1,k1
      zk(k)=zk(k)+d(i1+k,1)*uc
   67 continue
   68 do 70 i=1,i1
      b(i,1)=ph(i,1)/bg(i,1)
   70 continue
      return
      end
      subroutine leq(a,b,nn,mm,na,nb,err)
      implicit real*8(a-h,o-z)
c
      dimension a(na,nn),b(nb,mm)
      n=nn
      m=mm
      n1=n-1
      do 19 i=1,n
      r = dabs(a(i,1))
      do 16 j=2,n
      q = dabs(r)
      r = dmax1(q,dabs(a(i,j)))
   16 continue
      do 17 j=1,n
      a(i,j)=a(i,j)/r
   17 continue
      do 18 j=1,m
      b(i,j)=b(i,j)/r
   18 continue
   19 continue
      do 9 j=1,n1
      imax=j
      r = dabs(a(j,j))
      l=j+1
      do 2 i=l,n
      if (i .gt. n) goto 2
      if (r-dabs(a(i,j))) 1,2,2
    1 imax=i
      r = dabs(a(i,j))
    2 continue
      if (imax-j) 6,6,3
    3 do 4 k=j,n
         if( k .gt. n ) goto 4
         r=a(j,k)
         a(j,k)=a(imax,k)
         a(imax,k)=r
    4 continue
      do 5 k=1,m
         r=b(j,k)
         b(j,k)=b(imax,k)
         b(imax,k)=r
    5 continue
    6 do 9 i=l,n
      if (i .gt. n) goto 9
      r=-a(i,j)/a(j,j)
      do 7 k=l,n
      if (k .gt. n) goto 7
      a(i,k)=a(i,k)+r*a(j,k)
    7 continue
      do 8 k=1,m
      b(i,k)=b(i,k)+r*b(j,k)
    8 continue
    9 continue
      err=1.
      do 11 i=1,n
      if (a(i,i)) 11,10,11
   10 err=-1.0d0
      goto 15
   11 continue
      do 14 k=1,m
         b(n,k)=b(n,k)/a(n,n)
         do 13 l=1,n1
            i=n-l
            r = 0.0d0
            imax=i+1
            do 12 j=imax,n
               if (j .gt. n) goto 12
               jj=i+n+1-j
               r=r+a(i,jj)*b(jj,k)
   12      continue
           b(i,k)=(b(i,k)-r)/a(i,i)
   13    continue
   14 continue
   15 return
      end
      subroutine co1(n2,i1m,k1m,ik1m,ik2m,n2m,np,xi,yi,zk,ap,fd,ph,bg)
      implicit real*8(a-h,o-z)
c
      common u(700),v(700),voga1(700),ra(700),c1(700)
      dimension xi(n2m),yi(i1m,n2m),zk(k1m),ap(np),fd(i1m,ik2m,n2m),
     $  ph(ik1m,ik2m),bg(ik1m,i1m)
      integer s,smlp2,ap
c
      l = ap(1)
      rl1 = dfloat(l)*dfloat(l+1)
      s = ap(2)
      smlp2 = s-l+2
      do 100 n=1,n2
         vp1 = v(n) + 1.0d0
         c1zk = c1(n)*zk(1)
         fd(1,1,n) = (voga1(n)-3.0d0+smlp2)/vp1
         fd(1,2,n) = (rl1/c1zk-voga1(n))/vp1
         fd(1,3,n) = voga1(n)/vp1
       fd(1,5,n)=fd(1,1,n)*yi(1,n)+fd(1,2,n)*yi(2,n)+fd(1,3,n)*yi(3,n)
         fd(2,1,n) = (c1zk+ra(n))/vp1
         fd(2,2,n) = (1.0d0-u(n)-ra(n)+smlp2)/vp1
         fd(2,3,n) = ra(n)/vp1
       fd(2,5,n)=fd(2,1,n)*yi(1,n)+fd(2,2,n)*yi(2,n)+fd(2,3,n)*yi(3,n)
         fd(3,3,n) = (1.0d0-u(n)+smlp2)/vp1
         fd(3,4,n) = 1.0d0/vp1
         fd(3,5,n) = fd(3,3,n)*yi(3,n)+fd(3,4,n)*yi(4,n)
         fd(4,1,n) =-ra(n)*u(n)/vp1
         fd(4,2,n) = u(n)*voga1(n)/vp1
         fd(4,3,n) = (rl1-u(n)*voga1(n))/vp1
         fd(4,4,n) = (-u(n)+smlp2)/vp1
         fd(4,5,n)=fd(4,1,n)*yi(1,n)+fd(4,2,n)*yi(2,n) +
     $      fd(4,3,n)*yi(3,n) + fd(4,4,n)*yi(4,n)
  100 continue
      ph(1,3) = dfloat(l)
      ph(1,4) =-1.0d0
      ph(1,5) = ph(1,3)*yi(3,1)+ph(1,4)*yi(4,1)
      bg(2,1) = v(n2)-c1zk-4.0d0
      bg(2,2) = rl1/c1zk-v(n2)
      bg(2,3) = v(n2) - dfloat(l+1)
      ph(2,5) = bg(2,1)*yi(1,n2)+bg(2,2)*yi(2,n2)+bg(2,3)*yi(3,n2)
      bg(3,1) = u(n2)
      bg(3,3) = dfloat(l+1)
      bg(3,4) = 1.0d0
      ph(3,5) = bg(3,1)*yi(1,n2)+bg(3,3)*yi(3,n2)+bg(3,4)*yi(4,n2)
      bg(4,1) = 1.0d0
      ph(4,5) = bg(4,1)*yi(1,n2)-1.0d0
      return
      end
      subroutine co2 (n2,i1m,k1m,ik1m,ik2m,n2m,np,xi,yi,zk,ap,fd,ph,bg)
      implicit real*8(a-h,o-z)
c
      common u(700),v(700),voga1(700),ra(700),c1(700)
      dimension xi(n2m),yi(i1m,n2m),zk(k1m),ap(np),fd(i1m,ik2m,n2m),
     $  ph(ik1m,ik2m),bg(ik1m,i1m)
      integer s,smlp2,ap
      l = ap(1)
      rl1 = dfloat(l)*dfloat(l+1)
      s = ap(2)
      smlp2 = s-l+2
      do 100 n=1,n2
         vp1 = v(n)+1.0d0
         c1zk = c1(n)*zk(1)
         fd(1,1,n) = (voga1(n)-3.0d0+smlp2)/vp1
         fd(1,2,n) = (rl1/c1zk-voga1(n))/vp1
         fd(1,3,n) = voga1(n)/vp1
         fd(1,5,n) =-rl1/c1zk/zk(1)*yi(2,n)/vp1
       fd(1,6,n) = fd(1,1,n)*yi(1,n)+fd(1,2,n)*yi(2,n)+fd(1,3,n)*yi(3,n)
         fd(2,1,n) = (c1zk+ra(n))/vp1
         fd(2,2,n) = (1.0d0-u(n)-ra(n)+smlp2)/vp1
         fd(2,3,n) = ra(n)/vp1
         fd(2,5,n) = c1(n)*yi(1,n)/vp1
      fd(2,6,n) = fd(2,1,n)*yi(1,n)+fd(2,2,n)*yi(2,n)+fd(2,3,n)*yi(3,n)
         fd(3,3,n) = (1.0d0-u(n)+smlp2)/vp1
         fd(3,4,n) = 1.0d0/vp1
      fd(3,6,n) = fd(3,3,n)*yi(3,n)+fd(3,4,n)*yi(4,n)
         fd(4,1,n) =-ra(n)*u(n)/vp1
         fd(4,2,n) = u(n)*voga1(n)/vp1
         fd(4,3,n) = (rl1-u(n)*voga1(n))/vp1
         fd(4,4,n) = (-u(n)+smlp2)/vp1
      fd(4,6,n) = fd(4,1,n)*yi(1,n)+fd(4,2,n)*yi(2,n) +
     $      fd(4,3,n)*yi(3,n) + fd(4,4,n)*yi(4,n)
  100 continue
      ph(1,3) = dfloat(l)
      ph(1,4) =-1.0d0
      ph(1,6) = ph(1,3)*yi(3,1)+ph(1,4)*yi(4,1)
      bg(2,1) = v(n2)-c1zk-4.0d0
      bg(2,2) = rl1/c1zk-v(n2)
      bg(2,3) = v(n2) - dfloat(l+1)
      ph(2,5) =-c1(n2)*yi(1,n2)-rl1/c1zk/zk(1)*yi(2,n2)
      ph(2,6) = bg(2,1)*yi(1,n2)+bg(2,2)*yi(2,n2)+bg(2,3)*yi(3,n2)
      bg(3,1) = u(n2)
      bg(3,3) = dfloat(l+1)
      bg(3,4) = 1.0d0
      ph(3,6) = bg(3,1)*yi(1,n2)+bg(3,3)*yi(3,n2)+bg(3,4)*yi(4,n2)
      bg(4,1) = 1.0d0
      ph(4,6) = bg(4,1)*yi(1,n2)-1.0d0
      ph(5,1) = c1(1)*zk(1)
      ph(5,2) =-dfloat(l)
      ph(5,5) = c1(1)*yi(1,1)
      ph(5,6) = ph(5,1)*yi(1,1)+ph(5,2)*yi(2,1)
      return
      end
      subroutine pltint(nmax,fr,n,idstar)
      implicit real*8(a-h,o-z)
      character*80 idstar
c
c          initialize the plot file and write the title line.
c
      open (unit=12,file='adiab_plt',status='new')
      write(12,1000) idstar
      call pltdmp(fr,nmax,n,4hx   )
      return
 1000 format(x,80a)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p6e12.4) )
      end
