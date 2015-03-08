      subroutine lnanon(nin,iout,lin,nomgin,omlow,omhigh)
      implicit real*8(a-h,o-z)
      logical iout0, iout1
c
c          Stability analysis of linear, non-adiabatic nonradial
c      stellar pulsations. First version has the non-adiabatic
c      effects included in the momentum equation but the entropy
c      change will be zero everywhere for now. The notation is
c      based on the lagrangian variations of Castor Ap.J. 162(166)1971,
c      and the actual linearization is based on the technique of
c      Lynden-Bell and Ostriker.
c
c                                                9/15/83 WD Pesnell
c
      parameter ( accur=1.d-12 )
      parameter ( onth=0.333333333333333d0 )
      parameter ( nmax=600 )
      parameter ( nmax3=nmax*3 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	      dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
     $		     bv(nmax)
      common/blkdrs/ drz(nmax),drint(nmax),gor(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk37/  drdm(nmax,2),dl1(nmax,2),dl2(nmax,2)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $		     ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $		     ap1(nmax,2),	     ap3(nmax,3),ap4(nmax)
      common/linear/ dr(nmax),dh(nmax),gam(nmax),y(nmax),z(nmax),
     $		     dp(nmax),adrho(nmax),weight(nmax),cvt(nmax),
     $		     wtherm(nmax),wgrav(nmax),wdiag(nmax),wcross(nmax),
     $		     rkepl(nmax),stwait(nmax),wint(nmax),
     $		     spac(nmax,11)
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      dimension rkl2(nmax)
      dimension x(nmax3),stomeg(50),stomsq(50)
c
      n = nin
      np1 = n + 1
      nomega = nomgin
      lval = lin
      rl = float(lval)
      rl1 = rl*(rl+one)
      iout0 = .false.
      if( iout .gt. 0 ) iout0 = .true.
      iout1 = .false.
      if( iout .gt. 1 ) iout1 = .true.
      rhobar = totmas/(r(np1)**3*pi43)
c
c          The routine matset puts the matrices into order.
c
      call matset(n,iout,lval,rkl2,z)
      rz0 = ahf**onth*r(1)
      rkl20 = rl1/rz0**2
c         drdm01 = zero
c         drdm02 =-(rl+one)/rz0**2
         drdm01 =-ahf*(rl+one)/rz0**2
         drdm02 =-ahf*(rl+one)/rz0**2
c
c          convergence to a normal mode is governed by the momentum
c       equation for the outermost zone. as in the radial case, this
c       function (divided by any x value) is used in a secant method
c       for controlling the iterations. as the equations are not
c       nonlinear (as in non-non) in the frequency, it is hoped the
c       calculation is, in some sense, more stable.
c
c                                               10/4/83 wd pesnell
c
      write(11,3001)
      xnorm = r(np1)**2
      iroot = 0
      dom = (omhigh - omlow)/float(nomega - 1)
      do 30 i=1,nomega
         omsq = omlow + float(i-1)*dom
	 fomeg = fomsq(omsq,x,xnorm,n)
         if( iout0 ) write(11,4400) i,omsq,fomeg
         if( i .gt. 1 ) then
            if( fomeg*fomi .lt. zero ) then
               iroot = iroot + 1
               stomeg(iroot) = omsq - dom*fomeg/(fomeg-fomi)
            endif
         endif
         fomi = fomeg
  30  continue
      if( iroot .eq. 0 ) then
c
c          no roots in this interval, return to main program
c
         write(1,9000) lval,omlow,omhigh
         write(6,9000) lval,omlow,omhigh
         write(11,9000) lval,omlow,omhigh
         return
      else
         write(1,3000) lval,(stomeg(i),i=1,iroot)
         write(6,3000) lval,(stomeg(i),i=1,iroot)
         write(11,3000) lval,(stomeg(i),i=1,iroot)
      endif
c
c          attempt to converge on each of the frequencies found in
c       discriminant search. a secant method to find the zeroes of
c       error function is used.
c
      nomsq = 0
      do 100 nroot=1,iroot
         omsq = stomeg(nroot)
         do 45 icount=1,30
	    err = fomsq(omsq,x,xnorm,n)
            if( icount .eq. 1 ) then
               omsq1 = omsq
               erra  = err
               omsq  = omsq*(one + 1.d-7)
            else
               afac = dabs((omsq-omsq1)/omsq)
               if( iout0 ) write(11,4400) icount,omsq,afac,err,erra
               if( afac .le. accur ) goto 47
               omsq2 = (erra*omsq-err*omsq1)/(erra-err)
               omsq1 = omsq
               erra  = err
               omsq  = omsq2
            endif
  45     continue
         write(11,4500) nroot,omsq
         goto 100
c
c          converged to a root, calculate the eigenvectors.
c
  47     continue
         x(3*(n+1)) = xnorm
         call invint(omsq,x,finerr,xnorm,n)
c
c          do not repeat a root.
c
         nomsq = nomsq + 1
         if( nomsq .gt. 1 ) then
            do 491 ido=1,nomsq-1
               if( dabs((omsq-stomsq(ido))/omsq) .lt. accur ) then
                  nomsq = nomsq - 1
                  goto 100
               endif
 491        continue
         endif
         stomsq(nomsq) = omsq
         do 40 i=1,np1
	    dr(i) = (r(i)/r(np1))**lval*(x(3*i)/r(i)**2)
            if( i .eq. 1 ) then
               dh(i) = (x(3*i-1)/rz0**2)*(rz0/r(np1))**lval
               adrho(i) = ((drdm01+drdm02)*x(3*i) +
     $             rl1*x(3*i-1)/rz0**2 )*(rz0/r(np1))**lval
               dp(i) = g10*adrho(i)
               gam(i) = x(3*i-2)*(rz0/r(np1))**lval
	    else
               dh(i)=(x(3*i-1)/rzone(i-1)**2)*(rzone(i-1)/r(np1))**lval
               adrho(i) = (drdm(i-1,1)*x(3*i-3)+drdm(i-1,2)*x(3*i) +
     $           rl1*x(3*i-1)/rzone(i-1)**2)*(rzone(i-1)/r(np1))**lval
               dp(i) = g1(i-1)*adrho(i)
               gam(i) = x(3*i-2)*(rzone(i-1)/r(np1))**lval
	    endif
  40     continue
         dlnr = ahf*(rl+one)*drz(np1)/r(np1)
         gam(n+2) = gam(n)*(one-dlnr)/(one+dlnr)
         dh(n+2) = zero
	 qdr0 = dr(2)*(r(1)/r(2))**(lval-2)
	 idr = nodes(dr, n+1, 1)
	 idh = nodes(dh, n+1, 1)
c
c          weight function for the adiabatic oscillation.
c
         wtherm(1) = (g10*pc/rhoc)*cormas*adrho(1)**2
         wgrav(1) = dr(1)*dm2(1)*rl*gam(1)
	 wcross(1) = two*rl1*gor(1)*r(1)*dr(1)*dh(1)*dm2(1)
         wdiag(1) = (pi4g*rho(1)-gor(1)*for/r(1) )*
     $			 (dr(1)*r(1))**2*dm2(1)
	 weight(1) = wtherm(1) + wgrav(1) + wcross(1) + wdiag(1)
	 rgwait = weight(1)
	 wint(1) = weight(1)
	 rsum = (dr(1)*r(1))**2*dm2(1)
         hsum = (dh(1)*rz0)**2*cormas
         rkepl(1) = (dr(1)*r(1))**2*dm2(1) +
     $                   (dh(1)*rz0)**2*cormas
         rbar = r(1)*(dr(1)*r(1))**2*dm2(1)
         hbar = rz0*(dh(1)*rz0)**2*cormas
         crsum = zero
         do 48 i=1,n
c
c          weight function for the adiabatic oscillation.
c
            wtherm(i+1) = g1(i)*p(i)*v(i)*dm1(i)*adrho(i+1)**2
            wgrav(i+1) = rl1*dh(i+1)*dm1(i)*gam(i+1) +
     $      (dr(i+1)*r(i+1))*dm2(i+1)*(gam(i+2)-gam(i+1))/drz(i+1)
            wcross(i+1) = two*gor(i+1)*r(i+1)*dr(i+1)*dm2(i+1)*
     $            rl1*ahf*(dh(i+1)+dh(i+2))
            wdiag(i+1) = ( pi4g*ahf*(rho(i)+rho(i-1)) -
     $            for*gor(i+1)/r(i+1) )*(dr(i+1)*r(i+1))**2*dm2(i+1)
            rsum = rsum + (dr(i+1)*r(i+1))**2*dm2(i+1)
            hsum = hsum + (dh(i+1)*rzone(i))**2*dm1(i)
            rkepl(i+1) = (dr(i+1)*r(i+1))**2*dm2(i+1) +
     $                   (dh(i+1)*rzone(i))**2*dm1(i)
            rbar = rbar + r(i+1)*(dr(i+1)*r(i+1))**2*dm2(i+1)
            hbar = hbar + rzone(i)*(dh(i+1)*rzone(i))**2*dm1(i)
            weight(i+1) = wtherm(i+1)+wgrav(i+1)+wcross(i+1)+wdiag(i+1)
	    rgwait = rgwait + weight(i+1)
	    wint(i+1) = wint(i) + weight(i+1)
  48	 continue
c
c          normalize the weight function.
c
         rke = rsum + rl1*hsum
         if( rke .le. zero ) rke = one
         rbar = (rbar + rl1*hbar)/(rke*r(n+1))
         crsum = zero
         dh(n+2) = zero
         do 46 i=1,np1
            crsum = crsum + dr(i)*r(i)**2*dm2(i)*(dh(i)+dh(i+1))
	    rkepl(i) = rkepl(i)/rke
	    wint(i) = wint(i)/(rke*omsq)
	    weight(i) = weight(i)/(rke*omsq)
            wtherm(i) = wtherm(i)/(rke*omsq)
            wgrav(i) = wgrav(i)/(rke*omsq)
            wdiag(i) = wdiag(i)/(rke*omsq)
            wcross(i) = wcross(i)/(rke*omsq)
  46     continue
         rgwait = (rgwait/rke)
	 qch = dabs((omsq-rgwait)/omsq)
	 if( rgwait*omsq .le. 0. ) then
	    qch = abs(abs(omsq)-abs(rgwait))/abs(omsq)
	 endif
c
c          evaluate the rotational spitting coefficients.
c
         crsum = rl1*crsum/rke
         hrsum = hsum/rke
         ckl = crsum + hrsum
         index = knode(dr,dh,nmax,np1,np,ng)
         omega = dsqrt(omsq)
         period = twopi/omega
         call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax,itrans)
         pdays = period/86400.d0
         write(1,5505) nroot,omega,period,pdays,omsq,rgwait,qch
         write(1,5503) index,np,ng,itrans,rbar
         write(1,5506) ckl,crsum,hrsum
         write(6,5505) nroot,omega,period,pdays,omsq,rgwait,qch
         write(6,5503) index,np,ng,itrans,rbar
         write(11,5500) nroot,omega,period,pdays,omsq,rgwait,qch
         write(11,5503) index,np,ng,itrans,rbar
         write(11,5506) ckl,crsum,hrsum
c
         if( iout1 ) then
            write(11,5502)
            write(11,5501) (i,dr(i),dh(i),adrho(i),gam(i),
     $		   dp(i),weight(i),cvt(i),i=1,np1)
         endif
	 call pltdmp(dr,nmax,np1,'dr/r')
	 call pltdmp(dh,nmax,np1,'dh/h')
	 call pltdmp(adrho,nmax,np1,'drho')
	 call pltdmp(rkepl,nmax,np1,'ke  ')
	 call pltdmp(weight,nmax,np1,'wait')
	 call pltdmp(wint  ,nmax,np1,'wint')
	 call pltdmp(wtherm,nmax,np1,'wthr')
	 call pltdmp(wgrav,nmax,np1,'wgrv')
	 call pltdmp(wdiag,nmax,np1,'wdia')
	 call pltdmp(wcross,nmax,np1,'wcrs')
c
	 call orthog(dr,dh,n,omsq,index,1)
	 call eigenf(omsq,n,lval,rke,gor)
c
c	   non-adiabatic loop removed 10/30/89.
c
 100  continue
      call orthog(dr,dh,n,omsq,index,-1)
      return
c
 2200 format(1h1,//,20h ag1,ag2,ag3,ag4,ah4)
 2201 format(1x,i5,1p,10e12.4)
 2300 format(1h1,//,24h ah1,ah2,ah3,ap1,ap3,ap4)
 2600 format(1h1,//,10x,36h the pulsation matrices, ak1,ak2,ak4,/)
c
 3000 format(1x,23h roots for analysis: l=,i3,/,(3x,1p,5e15.6) )
 3001 format(1h1,1x,26hBegin discriminant search.)
c
 4400 format(1x,i5,1p,4e15.6)
 4500 format(1x,27hno convergence root number ,i4,
     $       11h frequency=,1p,2e15.6)
 5500 format(1h1,//,5x,4hroot,i5,/,8h omega =,1pe13.6,9h period =,e13.6,
     $   7h (secs),e15.6,7h (days),
     $   /,2x,33heigenvalue from matrix solution =,e12.4,
     $   /,2x,33heigenvalue from weight function =,e12.4,
     $   8h error =,e10.3)
 5501 format(1x,i5,1p,7e15.6)
 5502 format(/4x,1hi,8x,4hdr/r,10x,4hdh/h,10x,8hdrho/rho,8x,5hgamma,
     $	 11x,4hdp/p,10x,6hweight,10x,4hcvtm/)
 5503 format(1x,24h radial quantum number =,i4,
     $   26h number of nodes, p-type =,i4,9h g-type =,i4,/,
     $   1x,18h transition zone =,i4,/,
     $   1x,28h effective radius of mode = ,1pe11.4)
 5505 format(/,5x,4hroot,i5,/,8h omega =,1pe13.6,9h period =,e13.6,
     $   7h (secs),e15.6,7h (days),
     $   /,2x,33heigenvalue from matrix solution =,e12.4,
     $   /,2x,33heigenvalue from weight function =,e12.4,
     $   8h error =,e10.3)
 5506 format(1x,31h rotational splitting, total = ,1pe11.4,
     $   5h a = ,e11.4,5h b = ,e11.4)
c
 9000 format(1x,44h from lnanon...no adiabatic roots found, l =,i3,/,
     $   1x,17h frequency range:,1p,2e15.6)
      end
      function knode(x,y,nmax,nin,nps,ngs)
      implicit real*8(a-h,o-z)
c
c          returns the number of zero-crossings in the counter-
c       clockwise direction as np and the clockwise crossings
c       in ng.
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax),y(nmax)
      n = nin
      np = 0
      ng = 0
      knode = 0
      if( n .eq. 1 ) return
      nm1 = n - 1
      do 10 i=1,nm1
         if( x(i)*x(i+1) .lt. zero ) then
            if( y(i)*y(i+1) .lt. zero ) write(11,1000)
            if( x(i+1) .gt. x(i) ) then
               if( y(i) .ge. zero ) ng = ng + 1
               if( y(i) .lt. zero ) np = np + 1
            else
               if( y(i+1) .ge. zero ) np = np + 1
               if( y(i+1) .lt. zero ) ng = ng + 1
            endif
         endif
  10  continue
      if( x(n-1)*x(n) .lt. zero ) write(11,1001)
      knode = np - ng
      nps = np
      ngs = ng
      return
 1000 format(1x,50hFrom knode...quadrant jump default rotations used.)
 1001 format(1x,51hFrom knode...node at outer boundary is not counted.)
      end
