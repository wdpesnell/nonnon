      subroutine min_max( y, n )
      implicit real*8(a-h,o-z)
c
c          Find the maxima and minima in the vector Y.
c
      parameter ( nmax=401 )
      dimension y(nmax),imax(nmax),imin(nmax)
c
      nmaxs = 0
      nmin = 0
      ido = n - 2
      do 30 i=3,ido
         max = 0
         if( y(i) .gt. y(i-2) ) max = max + 1
         if( y(i) .gt. y(i-1) ) max = max + 1
         if( y(i) .gt. y(i+1) ) max = max + 1
         if( y(i) .gt. y(i+2) ) max = max + 1
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
      if( nmaxs .gt. 0 ) write(11,3001) (imax(i),y(imax(i)),i=1,nmaxs)
c
      write(11,4000) nmin
      if( nmin .gt. 0 ) write(11,3001) (imin(i),y(imin(i)),i=1,nmin)
c
      return
c
 3000 format(/,1x,i4,' maxima found in potential ',
     $  42hwhere zone(i) .gt. (i-2),(i-1),(i+1),(i+2),
     $  30h...listed below by zone, value)
 4000 format(/,1x,i4,' minima found in potential ',
     $  42hwhere zone(i) .le. (i-2),(i-1),(i+1),(i+2),
     $  30h...listed below by zone, value)
 3001 format(8(1x,i4,1p,e10.3))
      end
