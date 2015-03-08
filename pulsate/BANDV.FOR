      subroutine bandv(nm,n,mbw,a,ie21,m,w,z,ierr,nv,rv,rv6)
      implicit real*8(a-h,o-z)
      real*8 machep
      dimension a(nm,mbw), w(m), z(nm), rv(nv), rv6(n)
c
c	    EISPACK routine for solving a band matrix
c
      machep = 1.d-16
      ierr = 0
      if( m .eq. 0 ) return
      mb = mbw
      if( ie21 .lt. 0 ) mb = (mbw+1)/2
      m1 = mb - 1
      m21 = mb + m1
      iorder = 1 - iabs(ie21)
c
c	   Compute the norm of the matrix
c
      its = 1
      x1 = w(1)
      rnorm = 0.0
c
      do 10 j=1, mb
	 jj = mb + 1 - j
	 kj = jj + m1
	 ij = 1
	 do 15 i=jj, n
	    rnorm = rnorm + abs(a(i,j))
	    if( ie21 .lt. 0 ) then
	       rnorm = rnorm + abs(a(ij,kj))
	       ij = ij + 1
	    endif
  15	 continue
  10  continue
c
      if( ie21 .lt. 0 ) rnorm = rnorm / 2
      if( rnorm .eq. 0.0 ) norm = 1.0
      eps2 = 1.d-3 * rnorm *abs(iorder)
      eps3 = machep*rnorm
      uk = sqrt(float(n))
      eps4 = uk * eps3
  80  igroup = 0
      goto 120
 100  if( abs(x1-x0) .ge. eps2 ) goto 80
      igroup = igroup + 1
      if( iorder*(x1-x0) .le. 0.0 ) x1 = x0 + iorder*eps3
c
c	   Expand matrix, subtract eigenvalue, and initialize vector
c
 120  continue
      do 20 i=1,n
	 ij = i + min0(0, i-m1)*n
	 kj = ij + mb*n
	 ij1 = kj + m1*n
	 if( m1 .ne. 0 ) then
	    do 25 j=1, m1
	       if( ij .gt. m1 ) goto 125
	       if( ij .gt. 0 ) goto 130
		  rv(ij1) = 0.0
		  ij1 = ij1 + n
		  goto 130
 125	       continue
		  rv(i) = a(i,j)
 130	       continue
		  ij = ij + n
		  ii = i + j
		  if( ii .gt. n ) goto 25
		     jj = mb - j
		     if( ie21 .lt. 0 ) then
			ii = i
			jj = mb + j
		     endif
		     rv(kj) = a(ii,jj)
		     kj = kj + n
  25	    continue
c
	    rv(ij) = a(i,mb) - x1
	    rv6(i) = eps4
	    if( iorder .eq. 0 ) rv6(i) = z(i)
       endif
   20  continue
c
c
      if( m1 .eq. 0 ) goto 600
      do 580 i=1, n
	 ii = i + 1
	 maxk = min0(i+m1-1,n)
	 maxj = min0(n-i, m21-2)*n
	 do 360 k=1, maxk
	    kj1 = k
	    j = kj1 + n
	    jj = j + maxj
	    do 340 kj=j, jj, n
	       rv(kj1) = rv(kj)
	       kj1 = kj
 340	    continue
	    rv(kj1) = 0.0
 360	 continue
	 if( i .eq. n ) goto 580
	 u = 0.0
	 maxk = min0(I+m1, n)
	 maxj = min0(n-ii, m21-2)*n
	 do 450 j=1, maxk
	    if( abs(rv(j)) .lt. abs(u) ) goto 450
	    u = rv(j)
	    k = j
 450	 continue
	 j = i + n
	 jj = j + maxj
	 if( k .eq. i ) goto 520
	 kj = k
	 do 500 ij=i, jj, n
	    v = rv(ij)
	    rv(ij) = rv(kj)
	    rv(kj) = v
	    kj = kj + n
 500	 continue
c
	 if( iorder .eq. 0 ) then
	    v = rv6(i)
	    rv6(i) = rv6(k)
	    rv6(k) = v
	 endif
 520	 continue
	 if( u .ne. 0.0 ) then
	    do 560 k=ii, maxk
	       v = rv(k)/u
	       kj = k
	       do 540 ij=j, jj, n
		  kj = kj + n
		  rv(kj) = rv(kj) - v*rv(ij)
 540	       continue
	       if( iorder .eq. 0 ) rv6(k) = rv6(k) - v*rv6(i)
 560	    continue
	 endif
 580  continue
c
c	   Back substitution
c
 600  continue
      do 630 ii=1, n
	 i = n + 1 - ii
	 maxj = min0(ii, m21)
	 if( maxj .eq. 1 ) goto 620
	    ij1 = i
	    j = ij1 + n
	    jj = j + (maxj-2)*n
	    do 610 ij = j, jj, n
	       ij1 = ij1 + 1
	       rv6(i) = rv6(i) - rv(ij)*rv6(ij1)
 610	    continue
 620	    continue
	    v = rv(i)
	    if( abs(v) .lt. eps3 ) then
c
c	   Set error -- nearly singular linear system.
c
	       if( iorder .eq. 0 ) ierr =-1
	       v = sign(eps3, v)
	    endif
	    rv6(i) = rv6(i)/v
 630  continue
c
      xu = 1.0
      if( iorder .eq. 0 ) then
	 do 875 i=1, n
	    z(i) = rv6(i)
 875	 continue
	 return
      else
	 if( group .gt. 0 ) then
	    do 680 jj=1, group
	       j = r - group - 1 + jj
	       xu = 0.0
	       do 640 i=1, n
		  xu = xu + rv6(i)*z(i)
 640	       continue
	       do 660 i=1, n
		  rv6(i) = rv6(i) - xu*z(i)
 660	       continue
 680	    continue
	 endif
	 rnorm = 0.0
	 do 720 i=1, n
	    rnorm = rnorm + abs(rv6(i))
 720	 continue
	 if( rnorm .lt. 0.1 ) then
	    if( its .ge. n ) goto 830
	       its = uts + 1
	       xu = eps4/(uk + 1.0)
	       rv6(1) = eps4
	       do 760 i=2, n
		  rv6(i) = xu
 760	       continue
	       rv6(its) = rv6(its) - eps4*uk
	       goto 600
 830	    continue
	       ierr =-1
	       xu = 0.0
	       goto 870
	 endif
	 u = 0.0
	 do 860 i=1,n
	    u = u + rv6(i)**2
 860	 continue
	 xu = 1.0/sqrt(u)
 870	 continue
	 do 900 i=1,n
	    z(i) = rv6(i)*xu
 900	 continue
	 x0 = x1
      endif
c
      return
      end
