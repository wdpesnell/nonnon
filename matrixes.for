      subroutine rbmles(a,nmax,imin4,imax4,m4,y)
      implicit real*8 (a-h,o-z)
c
c          linear equation solver. (r)eal (b)and (m)atrix (l)inear
c       (e)quation (s)olver. solves a*x = y, destroying a and
c       returns the value of x in y.
c
      dimension a(nmax,m4),y(nmax)
      imin = imin4
      imax = imax4
      m = m4
      id = imax-imin+1
      idm = id-1
      mmid = (m+1)/2
      mmm = mmid-1
      do 20 ii=1,idm
         i = imax+1-ii
         den = 1.0d0/a(i,mmid)
         y(i) = y(i)*den
         kmax = min0(i-imin,mmm)
         do 10 j=1,mmm
            a(i,j) = a(i,j)*den
            do 15 kk=1,kmax
               jk = j+kk
               k = mmid+kk
               ik = i-kk
               a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  15        continue
  10     continue
         do 25 kk=1,kmax
            ik = i-kk
            k = mmid+kk
            y(ik) = y(ik)-y(i)*a(ik,k)
  25     continue
  20  continue
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
         kmax = min0(i-imin,mmm)
         do 35 kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
  35     continue
  30  continue
      return
      end
      subroutine trisol(a,imin4,imax4,x,y,z,nmax,e,f)
      implicit real*8 (a-h,o-z)
c
c          Solves the matrix equation (A-xI)y = z, returning
c       y given A, x, and z. Arrays e and f are scratch matrices
c       for storing the recursion relation.
c
      dimension a(nmax,3), y(nmax), z(nmax)
      dimension e(nmax),f(nmax)
      imin = imin4
      imax = imax4
      id = imax+1-imin
      e(imax) =-a(imax,1)/(a(imax,2)-x)
      f(imax) = y(imax)/(a(imax,2)-x)
      do 10 ii=2,id
         i = imax+1-ii
         den = (a(i,2)-x)+a(i,3)*e(i+1)
         if( ii .ne. id ) then
            e(i) =-a(i,1)/den
         endif
         f(i) = (y(i)-a(i,3)*f(i+1))/den
  10  continue
      imp = imin+1
      z(imin) = f(imin)
      do 20 i=imp,imax
         z(i) = e(i)*z(i-1)+f(i)
  20  continue
      return
      end
      subroutine cbmles(a,nmax,imin,imax,m,y)
      complex*16 a,y,den
c
c          linear equation solver. (c)omplex (b)and (m)atrix (l)inear
c       (e)quation (s)olver. solves a*x = y, destroying a and
c       returns the value of x in y.
c
      dimension a(nmax,m),y(nmax)
      id = imax-imin+1
      idm = id-1
      mmid = (m+1)/2
      mmm = mmid-1
      do 20 ii=1,idm
         i = imax+1-ii
c         den = (1.d0,0.d0)
c         if( dreal(a(i,mmid)).ne.0.d0 .or. dimag(a(i,mmid)).ne.0.d0)
c     $      den = 1.d0/a(i,mmid)
         den = 1.d0/a(i,mmid)
         y(i) = y(i)*den
         kmax = min0(i-imin,mmm)
         do 10 j=1,mmm
            a(i,j) = a(i,j)*den
            do 15 kk=1,kmax
               jk = j+kk
               k = mmid+kk
               ik = i-kk
               a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  15        continue
  10     continue
         do 25 kk=1,kmax
            ik = i-kk
            k = mmid+kk
            y(ik) = y(ik)-y(i)*a(ik,k)
  25     continue
  20  continue
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
         kmax = min0(i-imin,mmm)
         do 35 kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
  35     continue
  30  continue
      return
      end
      subroutine trcsol(a,nmax,imin,imax,x,y,z,e,f)
      implicit real*8 (a-h,o-z)
      complex*16 x,z,e,f,den
      dimension a(nmax,3),y(nmax),z(nmax)
      dimension e(nmax),f(nmax)
c
c          tridiagonal matrix solver for complex arrays.
c
      ido = imax+1 - imin
      e(imax) =-a(imax,1)/(a(imax,2)-x)
      f(imax) = y(imax)/(a(imax,2)-x)
      do 10 ii=2,ido
         i = imax+1 - ii
         den = (a(i,2)-x) + a(i,3)*e(i+1)
         if( ii .ne. ido ) then
            e(i) =-a(i,1)/den
         endif
         f(i) = (y(i)-a(i,3)*f(i+1))/den
  10  continue
      z(imin) = f(imin)
      ido = imin+1
      do 20 i=ido,imax
         z(i) = e(i)*z(i-1) + f(i)
  20  continue
      return
      end
