      subroutine cjhdmp(nin,y02)
      implicit real*8(a-h,o-z)
      character*1 iyorn
c
c          write out the input file for carl hansen's runge-kutta
c       pulsation analysis.
c
      parameter ( nmax=4096 )
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
