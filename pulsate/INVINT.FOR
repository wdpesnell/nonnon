      subroutine invint(omsq, x, finerr, xnorm, n)
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=600 )
      parameter ( nmax3=3*nmax )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $               dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),
     $               rm(nmax),bv(nmax)
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $		     ap1(nmax,2),	     ap3(nmax,3),ap4(nmax)
      common/linear/ atrix(nmax3,7), errvec(nmax3), xold(nmax3)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax3), dr(nmax), dh(nmax)
      parameter ( accur=1.d-10 )
c
c          Inverse interation to reduce the global error in the
c       eigensystem. A description can be found in Press et al.
c       Numerical Recipes, section 11.7.
c
       do 100 iverg=1,10
c
      errm = 0.d0
      errh = 0.d0
      errp = 0.d0
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
	 xoldn = xoldn + x(ix)*x(ix)*dm2(i)
c
c	   horizontal component of motion.
c
	 atrix(idh,1) = zero
	 atrix(idh,2) = ah1(i,1)
	 atrix(idh,3) = ah3(i)
	 atrix(idh,4) = ah4(i) - omsq
	 atrix(idh,5) = ah1(i,2)
	 atrix(idh,6) = zero
	 atrix(idh,7) = zero
	 xoldn = xoldn + x(idh)*x(idh)*dm1(i)
c
c	   poisson equation.
c
	 atrix(igam,1) = ap3(i,1)
	 atrix(igam,2) = zero
	 atrix(igam,3) = ap1(i,1)
	 atrix(igam,4) = ap3(i,2)
	 atrix(igam,5) = ap4(i)
	 atrix(igam,6) = ap1(i,2)
	 atrix(igam,7) = ap3(i,3)
c         xoldn = xoldn + dabs(x(igam))
 200  continue
      l = 0
      ii = 1
      do 210 i=1,3*(n+1)
	 jmin = max0(1, 5-i)
	 jmax = min0(7, 3*(n+1) - i + 4)
	 ereg = 0.0d0
	 l = l + 1
	 if( l .eq. 3 ) then
	    dr(ii) = ereg/x(i)
	    l = 0
	    ii = ii + 1
            do 230 j=jmin, jmax
               ereg = ereg + atrix(i,j)*x(i-4+j)*dm2(ii)
 230        continue
	 elseif( l .eq. 2 ) then
	    dh(ii) = ereg/x(i)
            do 235 j=jmin, jmax
               ereg = ereg + atrix(i,j)*x(i-4+j)*dm1(ii)
 235        continue
	 endif
	 errvec(i) = ereg
 210  continue
      do 25 i=1,3*(n+1)
	 xold(i) = x(i)
  25  continue
      call pltdmp(dr, nmax, n+1, 'drin')
      call pltdmp(dh, nmax, n+1, 'dhin')
c      call rbmles(atrix,nmax3,1,3*(n+1),7,x)
      call rbmles(atrix,nmax3,1,3*(n+1),7,errvec)
c
c          evaluate the error in the new eigenvector.
c
      err = 0.d0
      xoldn = 0.d0
      xnewn = 0.d0
      dlamt = 0.d0
      dlamb = 0.d0
      errmax = 0.d0
      errmax2 = 0.d0
      imax = 0
      imax2 = 0
      do 30 i=1,3*(n+1)
	 xnewn = xnewn + x(i)*x(i)
	 xoldn = xoldn + xold(i)*xold(i)
	 err = err + (errvec(i)/x(i))**2
	 if( dabs(errvec(i)/x(i)) .gt. errmax ) then
	    errmax2 = errmax
	    imax2 = imax
	    errmax = dabs(errvec(i)/x(i))
	    imax = i
	 endif
  30  continue
      write(6,*) ' errmax=',errmax,' zone=',imax
      write(6,*) ' errmax2=',errmax2,' zone2=',imax2
      if( errmax .gt. 5.e-2 ) then
	 uc = 5.d-2/errmax
      else
	 uc = 1.d0
      endif
      do 31 i=1,3*(n+1)
	 x(i) = x(i) - uc*errvec(i)
c	  dlamt = dlamt + x(i)*x(i)
c	  dlamb = dlamb + x(i)*xold(i)
  31  continue
      err = dsqrt(err)
      write(6,*) ' err=',err
c      if( x(3*(n+1))*xold(3*(n+1)) .lt. zero ) xoldn =-xoldn
c      do 35 i=1,3*(n+1)
c	  err = err + dabs((x(i)/xnewn)-(xold(i)/xoldn))
c  35  continue
      write(6,3000) err
      if( err .lt. accur ) then
	 xn = xnorm/x(3*(n+1))
	 do 40 i=1,3*(n+1)
	    x(i) = x(i)*xn
  40	 continue
	 write(6,4000) omsq,domsq
	 finerr = err
	 return
      else
	 do 45 i=1,3*(n+1)
	    x(i) = x(i)/xnewn
  45	 continue
      endif
c	  domsq = dlamt/dlamb*(xoldn/xnewn)
 100  continue
      xn = xnorm/x(3*(n+1))
      do 50 i=1,3*(n+1)
	 x(i) = x(i)*xn
  50  continue
      write(6,5000) omsq
      finerr =-1.d0
c
 2000 format(1x,'INVINT ... residual error is',1pe12.4)
 3000 format(1x,'INVINT ... error =',1pe12.4)
 4000 format(1x,'INVINT ... omsq and domsq are:',1p2e11.4)
 5000 format(1x,'INVINT ... no convergence, omsq is:',1p2e11.4)
      return
      end
