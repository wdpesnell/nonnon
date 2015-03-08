      subroutine invint(omsq,x,finerr,xnorm,nmaxin,n)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=512 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),            ap4(nmax)
      common/linear/ atrix(nmax3,7), errvec(nmax3), xold(nmax3),
     $               spac(nmax,77)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax3)
      parameter ( accur=1.d-6 )
c
c          Inverse interation to reduce the global error in the
c       eigensystem. A description can be found in Press et al.
c       Numerical Recipes, section 11.7.
c
       do 100 iverg=1,10
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega**2.
c
      err = 0.d0
      xoldn = 0.d0
      do 200 i=1,n+1
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
         atrix(ix,1) = ag1(i,1)
         atrix(ix,2) = ag3(i,1)
         atrix(ix,3) = ag4(i,1)
         atrix(ix,4) = ag1(i,2) - omsq
         atrix(ix,5) = ag3(i,2)
         atrix(ix,6) = ag4(i,2)
         atrix(ix,7) = ag1(i,3)
         errvec(ix) = atrix(ix,1)*x(ix-3) + atrix(ix,4)*x(ix) +
     $      atrix(ix,7)*x(ix+3) + atrix(ix,2)*x(ix-2) +
     $      atrix(ix,5)*x(ix+1) + atrix(ix,3)*x(ix-1) +
     $      atrix(ix,6)*x(ix+2)
         err = err + dabs(errvec(ix))
         xoldn = xoldn + dabs(x(ix))
c
c          horizontal component of motion.
c
         atrix(idh,1) = zero
         atrix(idh,2) = ah1(i,1)
         atrix(idh,3) = ah3(i)
         atrix(idh,4) = ah4(i) - omsq
         atrix(idh,5) = ah1(i,2)
         atrix(idh,6) = zero
         atrix(idh,7) = zero
         errvec(idh) = atrix(idh,3)*x(idh-1) + atrix(idh,4)*x(idh) +
     $         atrix(idh,2)*x(idh-2) + atrix(idh,5)*x(idh+1)
         err = err + dabs(errvec(idh))
         xoldn = xoldn + dabs(x(idh))
c
c          poisson equation.
c
         atrix(igam,1) = ap3(i,1)
         atrix(igam,2) = zero
         atrix(igam,3) = ap1(i,1)
         atrix(igam,4) = ap3(i,2)
         atrix(igam,5) = ap4(i)
         atrix(igam,6) = ap1(i,2)
         atrix(igam,7) = ap3(i,3)
         errvec(igam) = atrix(igam,1)*x(igam-3) + atrix(igam,4)*x(igam)
     $       + atrix(igam,7)*x(igam+3) + atrix(igam,3)*x(igam-1) +
     $         atrix(igam,6)*x(igam+2) + atrix(igam,5)*x(igam+1)
         err = err + dabs(errvec(igam))
         xoldn = xoldn + dabs(x(igam))
 200  continue
      err = err/xoldn
      write(6,2000) err
      write(11,2000) err
      do 25 i=1,3*(n+1)
         xold(i) = x(i)
  25  continue
      call rbmles(atrix,nmax3,1,3*(n+1),7,x)
c
c          evaluate the error in the new eigenvector.
c
      err = 0.d0
      xoldn = 0.d0
      xnewn = 0.d0
      dlamt = 0.d0
      dlamb = 0.d0
      do 30 i=1,3*(n+1)
         xnewn = xnewn + dabs(x(i))
         xoldn = xoldn + dabs(xold(i))
         dlamt = dlamt + x(i)*x(i)
         dlamb = dlamb + x(i)*xold(i)
  30  continue
      if( x(3*(n+1))*xold(3*(n+1)) .lt. zero ) xoldn =-xoldn
      do 35 i=1,3*(n+1)
         err = err + dabs((x(i)/xnewn)-(xold(i)/xoldn))
  35  continue
      write(6,3000) err
      if( err .lt. accur ) then
         xn = xnorm/x(3*(n+1))
         do 40 i=1,3*(n+1)
            x(i) = x(i)*xn
  40     continue
         write(6,4000) omsq,domsq
         finerr = err
         return
      endif
         domsq = dlamt/dlamb*(xoldn/xnewn)
 100  continue
         xn = xnorm/x(3*(n+1))
         do 50 i=1,3*(n+1)
            x(i) = x(i)*xn
  50     continue
         write(6,4000) omsq,domsq
      finerr =-1.d0
      return
c
 2000 format(1x,'INVINT ... residual error is',1pe12.4)
 3000 format(1x,'INVINT... error =',1pe12.4)
 4000 format(1x,'INVINT... omsq and domsq are:',1p2e11.4)
      end
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
      subroutine orthog(dr,dh,n,omsq,index,nmaxin,imode)
      implicit real*8(a-h,o-z)
c
c          check for orthogonality of the wave functions. part
c      one stores the vertical and horizontal wavefunctions,
c      part two does the integrations.
c
      parameter ( nmax=512 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	 dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/eigstr/ sdr(nmax,4),sdh(nmax,4),stome(4),iomeg(4)
      dimension dr(nmaxin),dh(nmaxin),rke(4),cross(4,4)
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
      do 10 i=1,n
         sdr(i,nmode) = dr(i)
         sdh(i,nmode) = dh(i)
  10  continue
      sdr(n+1,nmode) = dr(n+1)
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
      subroutine eigenf(nmaxin,omsq,nin,lin,rke,gor, echout)
      implicit real*8(a-h,o-z)
      integer*2 i
c
c          Find the normalized eigenvectors and the weight functions
c       as a function of position.
c
      parameter ( nmax=512 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	 dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/linear/ spac(nmax,22),ra(nmax),vn(nmax),un(nmax),
     $	    dr(nmax),dh(nmax),gam(nmax),
     $	    qnyr(nmax),gyr(nmax),dp(nmax),adrho(nmax),weight(nmax),
     $	    cyr(nmax),xo(nmax),yo(nmax),qwait(nmax),
     $	    y1(nmax),y2(nmax),y3(nmax),y4(nmax),spac1(nmax,63)
      dimension vq(6),gor(nmax)
c
      n = nin
      n4 = nin
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
c          Schwank weight functions
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
c         gyr(i) =-temp*(y4(i) + (rl+one)*y3(i))**2/vq(4)
c4Horsemen         gyr(i) =-temp*(y4(i)**2 + rl1*y3(i)**2)/vq(4)
cKawaler         gyr(i) =-temp*(y4(i) + rl1*y3(i))**2/vq(4)
         gyr(i) =-temp*(y4(i) + rl*y3(i))**2/vq(4)
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
      qrke = rke*omsq/two
      write(6,2000) omsq,stwait,qch,qrke
      write(1,2000) omsq,stwait,qch,qrke
      write(11,2000) omsq,stwait,qch,qrke
      call pltdmp(weight,nmax,n4,'epst')
      call pltdmp(cyr,nmax,n4,'cyr ')
      call pltdmp(qnyr,nmax,n4,'nyr ')
      call pltdmp(gyr,nmax,n4,'gyr ')
c
      echout = qch
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
      function fomsq(omsq,x,xnorm,nmaxin,n)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=512 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),            ap4(nmax)
      common/linear/ atrix(nmax3,7),spac(nmax,83)
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
      function fomega(omega,x,xnorm,nmaxin,n)
      implicit real*8(a-h,o-z)
      complex*16 fomega,atrix,x,omsq,omega,err,czero,ci
c
      parameter ( czero=(0.d0,0.d0), ci=(0.d0,1.d0) )
      parameter ( nmax=512 )
      parameter ( nmax4=4*nmax )
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),            ap3(nmax,3),ap4(nmax)
      common/linear/ atrix(nmax4,13)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common ak1(nmax,4),ak2(nmax,3),ak4(nmax,3)
      dimension x(nmax4)
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega.
c
      omsq = omega*omega
      np1 = n + 1
      do 200 i=1,np1
c
c      ix is index for radial component, idh for horizontal,
c    igam is for the poisson equation and ient is the entropy index.
c
         ix = 4*i
         idh = ix - 1
         igam = ix - 3
         ient = ix - 2
c
c          radial component of motion.
c
         x(ix) = czero
         atrix(ix,1) = czero
         atrix(ix,2) = czero
         atrix(ix,3) = ag1(i,1)
         atrix(ix,4) = ag3(i,1)
         atrix(ix,5) = ag2(i,1)
         atrix(ix,6) = ag4(i,1)
         atrix(ix,7) = ag1(i,2) - omsq
         atrix(ix,8) = ag3(i,2)
         atrix(ix,9) = ag2(i,2)
         atrix(ix,10) = ag4(i,2)
         atrix(ix,11) = ag1(i,3)
         atrix(ix,12) = czero
         atrix(ix,13) = czero
c
c          horizontal component of motion.
c
         x(idh) = czero
         atrix(idh,1) = czero
         atrix(idh,2) = czero
         atrix(idh,3) = czero
         atrix(idh,4) = ah1(i,1)
         atrix(idh,5) = ah3(i)
         atrix(idh,6) = ah2(i)
         atrix(idh,7) = ah4(i) - omsq
         atrix(idh,8) = ah1(i,2)
         atrix(idh,9) = czero
         atrix(idh,10) = czero
         atrix(idh,11) = czero
         atrix(idh,12) = czero
         atrix(idh,13) = czero
c
c          Poisson equation.
c
         x(igam) = czero
         atrix(igam,1) = czero
         atrix(igam,2) = czero
         atrix(igam,3) = ap3(i,1)
         atrix(igam,4) = czero
         atrix(igam,5) = czero
         atrix(igam,6) = ap1(i,1)
         atrix(igam,7) = ap3(i,2)
         atrix(igam,8) = czero
         atrix(igam,9) = ap4(i)
         atrix(igam,10) = ap1(i,2)
         atrix(igam,11) = ap3(i,3)
         atrix(igam,12) = czero
         atrix(igam,13) = czero
c
c          thermal equations.
c
         x(ient) = czero
         atrix(ient,1) = ak1(i,1)
         atrix(ient,2) = czero
         atrix(ient,3) = ak2(i,1)
         atrix(ient,4) = ak4(i,1)
         atrix(ient,5) = ak1(i,2)
         atrix(ient,6) = czero
         atrix(ient,7) = ak2(i,2) - ci*omega
         atrix(ient,8) = ak4(i,2)
         atrix(ient,9) = ak1(i,3)
         atrix(ient,10) = czero
         atrix(ient,11) = ak2(i,3)
         atrix(ient,12) = ak4(i,3)
         atrix(ient,13) = ak1(i,4)
 200  continue
      x(4*np1-1) =-ah1(np1,2)*xnorm
      x(4*np1-2) =-ak1(np1,3)*xnorm
      x(4*np1-3) =-ap1(np1,2)*xnorm
      x(4*np1-4) =-ag1(n,3)*xnorm
      x(4*np1-6) =-ak1(n,4)*xnorm
      call cbmles(atrix,nmax4,1,4*np1-1,13,x)
      err = ag1(n+1,1)*x(4*np1-4) + (ag1(n+1,2)-omsq)*xnorm +
     $ ag3(n+1,1)*x(4*np1-3)+ag4(n+1,1)*x(4*np1-1)+ag2(n+1,1)*x(4*np1-2)
      fomega = err/x(4)
      return
      end
      subroutine amppha(ca,amp,phase)
      implicit real*8(a-h,o-z)
      complex*16 ca
c
c          returns the magnitude of ca in amp and the phase
c       of ca in phase.
c
      amp = cdabs(ca)
      phase = datan2(dimag(ca),dreal(ca))
      return
      end
      function cool(i,qc)
      implicit real*8(a-h,o-z)
c
c      returns the value of luminosity out of a zone as
c   a given function of qc, q1, and q0.
c
c          noburn:= 0, normal nuclear profile.
c                   1, linear ramp.
c                   2, quadratic ramp.
c                   3, square root ramp.
c                   4, inputted luminosity.
c
c          needs in common block coolum:
c                 onemq0:= 1.-q0(inner mass fraction)
c                 onemq1:= 1.-q1(outer mass fraction)
c                 rlumc:= core luminosity
c                 rlums:= surface luminosity
c
      parameter( nmax=512 )
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
c
      noburntst = noburn
      if( noburntst .le. 0 ) then
         return
      elseif( noburntst .eq. 4 ) then
         cool = frft(i)
      elseif( qc .gt. 1.d0-onemq1 ) then
         cool = rlums
      elseif( qc .lt. 1.d0-onemq0 ) then
         cool = rlumc
      else
         f = (qc+onemq0-1.d0)/(onemq0-onemq1)
         if( noburntst .eq. 1 ) then
            cool = rlumc+(rlums-rlumc)*f
         elseif( noburntst .eq. 2 ) then
            cool = rlumc+(rlums-rlumc)*f*f
         elseif( noburntst .eq. 3 ) then
            cool = rlumc+(rlums-rlumc)*sqrt(f)
         endif
      endif
      return
      end
      subroutine zie(t,n,nzmax)
      implicit real*8(a-h,o-z)
c
c          find the maxima and minima in the work per zone function.
c
      parameter ( nmax=512 )
      dimension t(nzmax),imax(nmax),imin(nmax)
c
      nmaxs = 0
      nmin = 0
      ido = n - 2
      do 30 i=3,ido
         max = 0
         if( t(i) .gt. t(i-2) ) max = max + 1
         if( t(i) .gt. t(i-1) ) max = max + 1
         if( t(i) .gt. t(i+1) ) max = max + 1
         if( t(i) .gt. t(i+2) ) max = max + 1
c
         if( max .eq. 0 ) then
c
c          record a minimum
c
            nmin = nmin + 1
            if( nmin .gt. nmax ) goto 30
            imin(nmin) = i
c
         elseif( max .eq. 4 ) then
c
c          record a maximum
c
            nmaxs = nmaxs + 1
            if( nmaxs .gt. nmax ) goto 30
            imax(nmaxs) = i
         endif
  30  continue
      write(11,3000) nmaxs
      if( nmaxs .gt. 0 ) write(11,3001) (imax(i),t(imax(i)),i=1,nmaxs)
c
      write(11,4000) nmin
      if( nmin .gt. 0 ) write(11,3001) (imin(i),t(imin(i)),i=1,nmin)
c
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
      implicit real*8(a-h,o-z)
c
c          find all maxima in normalized work function.
c
      dimension w(nzmax),imax(200)
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
      if( nmax .gt. 0 ) write(11,1001) (imax(i),w(imax(i)),i=1,nmax)
      return
c
 1000 format(/,1x,i4,41h maxima found in normalized work function,
     $  30h...listed below by zone, value)
 1001 format(8(2x,i4,1pe10.3))
      end
      subroutine cptdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
      complex*16 vec
c
c         write the vector vec to the plot file (file12) with appended
c       title ititl. In the version the real and imaginary parts are
c       dumped separately.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(dreal(vec(i)),i=1,n)
      write(12,1001) n,ititl,(dimag(vec(i)),i=1,n)
      return
c
 1000 format(i4,10x,a4,1x,4hreal,/,(1p,6e12.4) )
 1001 format(i4,10x,a4,1x,4himag,/,(1p,6e12.4) )
      end
