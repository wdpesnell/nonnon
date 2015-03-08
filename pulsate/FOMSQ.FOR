      function fomsq(omsq,x,xnorm,n,lu)
      implicit real*8(a-h,o-z)
      logical lu
c
      parameter ( nmax=400 )
      parameter ( nmax3=3*nmax )
      common/scrtch/ ag1(nmax,3),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),ap4(nmax)
      common/linear/ atrix(nmax3,7)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax3)
      dimension alu(nmax3,nmax3), indx(nmax3)
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega**2.
c
      do 200 i=1,n
c
c      ix is index for radial component, idh for horizontal and
c    igam is for the poisson equation.
c
         ix = 3*i - 2
         idh = ix + 2
         igam = ix + 1
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
      x(3*n) =-ah1(n,2)*xnorm
      x(3*n-1) =-ap1(n,2)*xnorm
      x(3*n-2) =-ag1(n,3)*xnorm
      if( lu ) then
         do 30 i=1,3*n
            jmin = max0(1,5-i)
            jmax = min0(7, 3*n + 4 - i)
            do 35 j=jmin,jmaX
               alu(i, i-4+j) = atrix(i,j)
  35        continue
  30     continue
         call ludcmp(alu,3*n,nmax3,indx,d)
         call lubksb(alu,3*n,nmax3,indx,x)
      else
         call rbmles(atrix,nmax3,1,3*n,7,x)
      endif
      err = ag1(n+1,1)*x(3*n-2) + (ag1(n+1,2)-omsq)*xnorm +
     $      ag3(n+1,1)*x(3*n-1) + ag4(n+1,1)*x(3*n)
      fomsq = err/x(4)
      return
      end
