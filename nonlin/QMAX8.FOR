      function qmax8(array,npts)
      implicit real*8(a-h, o-z)
      dimension array(1)
c
c          find the maximum element in array.
c
      nstart = 1
      xmax = array(nstart)
      ido = 6*(npts/6)
      do 10 i=1,ido,6
         xmax = dmax1(xmax,array(i),array(i+1),array(i+2),array(i+3),
     $                  array(i+4),array(i+5) )
  10  continue
      if( ido .ne. npts ) then
         ido = ido + 1
         do 30 i=ido,npts
            xmax = dmax1(xmax,array(i))
  30     continue
      endif
      qmax = xmax
      return
      end
      function qmin8(array,npts)
      implicit real*8(a-h, o-z)
      dimension array(1)
c
c          find the minimum element in array.
c
      xmin = array(1)
      ido = 6*(npts/6)
      do 10 i=1,ido,6
         xmin = dmin1(xmin,array(i),array(i+1),array(i+2),array(i+3),
     $                  array(i+4),array(i+5) )
  10  continue
      if( ido .ne. npts ) then
         ido = ido + 1
         do 30 i=ido,npts
            xmin = dmin1(xmin,array(i))
  30     continue
      endif
      qmin = xmin
      return
      end
