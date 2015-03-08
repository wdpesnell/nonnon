      subroutine rbmles(a,nmax,imin,imax,m,y)
      implicit real*8(a-h,o-z)
c
c          Linear equation solver. (R)eal (B)and (M)atrix (L)inear
c       (E)quation (S)olver. Solves A*X = Y, destroying A and
c       returns the value of X in Y.
c
      dimension a(nmax,m),y(nmax)
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
            do 10 kk=1,kmax
             jk = j+kk
             k = mmid+kk
             ik = i-kk
             a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  10     continue
         do 20 kk=1,kmax
            ik = i-kk
            k = mmid+kk
            y(ik) = y(ik)-y(i)*a(ik,k)
  20  continue
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
         kmax = min0(i-imin,mmm)
         do 30 kk=1,kmax
            ik = i-kk
            k = mmid-kk
            y(i) = y(i)-a(i,k)*y(ik)
30    continue
      return
      end
