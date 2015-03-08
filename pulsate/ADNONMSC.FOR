      function nodes(w,nz,lout)
      implicit real*8(a-h,o-z)
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
      subroutine orthog(dr,dh,n,omsq,index,imode)
      implicit real*8(a-h,o-z)
c
c          check for orthogonality of the wave functions. part
c      one stores the vertical and horizontal wavefunctions,
c      part two does the integrations.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	 dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
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
	    if( dabs((omsq-stome(i))/stome(i)) .lt. 1.d-6 ) return
   5	 continue
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
      implicit real*8(a-h,o-z)
c
c          find the normalized eigenvectors and the weight functions
c       as a function of position.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	 dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/linear/ dr(nmax),dh(nmax),gam(nmax),qnyr(nmax),gyr(nmax),
     $		     dp(nmax),adrho(nmax),weight(nmax),cyr(nmax),
     $		     xo(nmax),yo(nmax),qwait(nmax),ra(nmax),vn(nmax),
     $		     un(nmax),y1(nmax),y2(nmax),y3(nmax),y4(nmax),
     $		     spac1(nmax,8)
      dimension vq(6),gor(nmax)
c
      n = nin
      lval = lin
      rl = float(lval)
      rl1 = rl*(rl+one)
      do 10 i=1,n
         vn(i) = ahf*(r(i+1)*gor(i+1) + r(i)*gor(i))/(p(i)*v(i))
         un(i) = forpi*r(i+1)**3/rm(i+1)/v(i)
	 ra(i) =-bv(i+1)*r(i+1)/gor(i+1)
         xo(i) = dlog(rzone(i)/p(i))
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
      qch = dabs((omsq-stwait)/omsq)
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
      implicit real*8(a-h,o-z)
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
         if( x(j) .gt. 1.d0 ) goto 10
         ist = j
  10  continue
      itrans = ist
      return
      end
      function fomsq(omsq,x,xnorm,n)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=600 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $		     ap1(nmax,2),ap3(nmax,3),		 ap4(nmax)
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
