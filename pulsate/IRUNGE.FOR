      function irunge(n,y,f,x,h,m)
      implicit real*8(a-h,o-z)
      dimension phi(20),savey(20),y(n),f(n)
c
      if( m .eq. 1 ) then
c
c          step 1.
c
	 irunge = 1
      elseif( m .eq. 2 ) then
c
c          step 2.
c
	 do 25 j=1,n
            savey(j) = y(j)
            phi(j) = f(j)
            y(j) = y(j) + h*f(j)/2.d0
  25     continue
         x = x + h/2.d0
         irunge = 1
      elseif( m .eq. 3 ) then
c
c         step 3.
c
	 do 35 j=1,n
            phi(j) = phi(j) + 2.d0*f(j)
            y(j) = savey(j) + f(j)*h/2.d0
  35     continue
         irunge = 1
      elseif( m .eq. 4 ) then
c
c          step 4.
c
	 do 45 j=1,n
            phi(j) = phi(j) + 2.d0*f(j)
            y(j) = savey(j) + f(j)*h
  45     continue
         x = x + h/2.d0
         irunge = 1
      elseif( m .eq. 5 ) then
c
c         final pass.
c
	 do 55 j=1,n
            y(j) = savey(j) + (phi(j) + f(j))*h/6.d0
  55     continue
         irunge = 0
      endif
      return
      end
