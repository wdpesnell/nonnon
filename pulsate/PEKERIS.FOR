      subroutine pekeris( kin, lin, xin, alfa, npts)
      implicit real*8( a-h, o-z)
      parameter ( nmax=4096 )
      dimension xin(nmax), alfa(nmax), xn(nmax), pgam(nmax)
      integer*2 l, k, i
c
c	   Calculate the Pekeris eigenfunctions for the density
c	variations.
c
      l = lin
      k = iabs( kin )
      B = dfloat( (2*k)*( 2*(k+l) + 5) )
      cj = 1.d0
      rn = xin(npts)**2
      do 10 i=1, npts
	 xn(i) = xin(i)*xin(i)
	 alfa(i) = 1.d0
         pgam(i) = rn/2.d0 - xn(i)*dfloat(2*l+1)/2.d0/dfloat(2*l+3)
  10  continue
      if( k .gt. 0 ) then
	 do 20 j=2, k, 2
	    cj = cj*( dfloat( 2*(j-1)*( 2*(j+l) + 3) ) - B)/
     $		      dfloat( (2*j)*( 2*(j+l) + 1) )
            rn = rn*xin(npts)**2
	    do 25 i=1, npts
	       alfa(i) = alfa(i) + cj*xn(i)
	       xn(i) = xn(i)*xin(i)*xin(i)
               pgam(i) = pgam(i) + cj*( rn/dfloat(j+2) +
     $              xn(i)*( 1.d0/dfloat(2*l+3+j) - 1.d0/dfloat(j+2) ))
  25	    continue
  20	 continue
      endif
      do 30 i=1, npts
	 alfa(i) =-alfa(i)*xin(i)**l
	 pgam(i) =-pgam(i)*xin(i)**l/dfloat(2*l+1)
  30  continue
      call pltdmp( alfa, nmax, npts, 'peka')
      call pltdmp( pgam, nmax, npts, 'pgam')
      return
      end
