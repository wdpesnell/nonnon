      subroutine matset(nin,iout,lin,rkl2,z)
      implicit real*8(a-h,o-z)
      logical rzdone
c
c         Eigen-analysis of linear, adiabatic nonradial
c      stellar pulsations. This version has no non-adiabatic
c      effects included in the momentum equation. The notation is
c      based on the Lagrangian variations of Castor Ap.J. 162(166)1971,
c      and the actual linearization is based on the technique of
c      Lynden-Bell and Ostriker.
c
c                                                9/15/83 WD Pesnell
c
      parameter( onth=1.d0/3.d0 )
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $               dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),
     $               rm(nmax),bv(nmax)
      common/blkdrs/ drz(nmax),drint(nmax),gor(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk37/  drdm(nmax,2),dl1(nmax,2),dl2(nmax,2)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $		     ap1(nmax,2),	     ap3(nmax,3),ap4(nmax)
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/thermo/ rther,a,bk,avagd,ad3
      dimension rkl2(nmax),z(nmax)
      save rzdone
      data rzdone/.false./
c
      n = nin
      np1 = n + 1
      lval = lin
      rl = float(lval)
      rl1 = rl*(rl+one)
      rhobar = totmas/(r(np1)**3*pi43)
c
c          define useful quantities.
c
      rz0 = (ahf**onth)*r(1)
      rkl20 = rl1/rz0**2
      do 10 i=1,n
         rho(i) = one/v(i)
         dkdr(i) = dkdr(i) + g3m1(i)*dkdt(i)
         rzone(i) = (r(i)**3 + ahf*dm1(i)*v(i)/pi43)**onth
         rkl2(i) = rl1/rzone(i)**2
         if( i .eq. 1 ) then
            drz(i) = rzone(i)-rz0
         else
            drz(i) = rzone(i) - rzone(i-1)
         endif
         drint(i) = r(i+1)-r(i)
         gor(i) = g*rm(i)/r(i)**2
         drdm(i,1) = forpi*rzone(i)/(dm1(i)*v(i)) -
     $                  ahf*(rl+one)/rzone(i)**2
         drdm(i,2) =-(forpi*rzone(i)/(dm1(i)*v(i)) +
     $			ahf*(rl+one)/rzone(i)**2)
	 if( .not. rzdone ) z(i) = rzone(i)/r(np1)
  10  continue
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      if( .not. rzdone ) then
	 call pltdmp(z,nmax,n,'rz  ')
	 rzdone = .true.
      endif
      gor(np1) = g*rm(np1)/r(np1)**2
      drz(np1) = p(n)*v(n)/gor(np1)
      pi4g = forpi*g
c
c          initialize the matrices.
c
c          all matrices are prefixed with the letter a, the second
c       letter denotes the equation ( g=radial component of the
c       momentum equation, h=horizontal component of the momentum
c       equation, p=poisson equation and k= thermal equation)
c       in which each matrix is used. the number tells which
c       variable is to be acted on the the matrix: 1= radial
c       component of the motion eigenvector (mechanical
c       eigenvector), 2= entropy variations (to within a factor
c       of the temperature, actuallly the thermal eigenvector),
c       3= gravitation variations multiplied by the square root
c       of the interface radius differences), 4= horizontal
c       component of the mechanical eigenvector. it should be
c       noted that there is no ap2 matrix nor is ak3 present.
c
c                                              11/28/83 wd pesnell
c
      do 20 i=1,np1
      if( i .eq. 1 ) then
c
c          inner boundary conditions.
c
         afac = forpi*r(1)**3/dm2(1)
         vl = g1(1)*p(1)*(ahf*rl*v(1) + afac)
         drdm01 =-ahf*(rl+one)/rz0**2
         drdm02 =-ahf*(rl+one)/rz0**2
         ag11 = drdm01*g10*pc*(ahf*rl/rhoc - afac)
         ag1(1,1) = zero
         ag1(1,2) =-gor(1)/r(1) + drdm(1,1)*vl +
     $             drdm02*g10*pc*(ahf*rl/rhoc - afac)
         ag1(1,3) = drdm(1,2)*vl + ag11
c         ag3(i,1) = ahf*rl - r(i)/drz(i)
c         ag3(i,2) = ahf*rl + r(i)/drz(i)
         ag3(i,1) = zero
         ag3(i,2) = rl
         ag4(1,1) = zero
         ag41 = g10*pc*(ahf*rl/rhoc - afac)*rl1/rz0**2
c         ag4(i,1) = ag41 + ahf*rl1*gor(1)/r(1)
         ag4(1,2) = rkl2(1)*vl + ag41 + rl1*gor(1)/r(1)
c         ag4(i,2) = rkl2(i)*vl + ahf*rl1*gor(i)/r(i)
c
c          horizontal momentum equation.
c
c         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
c     $         ahf*gor(i-1)/r(i-1)
         ah1(i,1) = zero
         ah1(i,2) = (drdm01+drdm02)*g10*pc/rhoc + gor(i)/r(i)
         ah3(i) = one
         ah4(i) = rkl20*g10*pc/rhoc
c         ap1(1,1) =-pi4g*rho(1)*rzone(1)**2*drdm(1,1)
         ap1(1,1) = zero
c         ap1(1,2) =-pi4g*rho(1)*(rzone(1)**2*drdm(1,2)-
c     $     ahf*dlog(rho(2)/rho(1))/dlog(rzone(2)/rzone(1)) )
         ap1(1,2) =-pi4g*rhoc*(rz0**2*drdm02-
     $     dlog(rho(1)/rhoc)/dlog(rzone(1)/rz0) )
         ap311 = r(1)**2/(drint(1)*drz(1))
         ap313 = r(2)**2/(drint(1)*drz(2))
         ap3(1,1) = zero
         ap3(1,2) =-ap313 - ap311*(drz(1)/drz(2))**2
         ap3(1,3) = ap313 + ap311*(drz(1)/drz(2))**2
         ap4(1)   =-pi4g*rl1*rhoc
      elseif( i .eq. np1 ) then
c
c          outer boundary conditions
c
	 dlnr = (rl+ahf)*drz(np1)/r(np1)
	 afac = forpi*r(np1)**3/dm2(np1)
         g1p = g1(n)*p(n)-for*ad3*t(n)**4*g3m1(n)
         g3m1ro = g3m1(n)/v(n)-for*ad3*t(n)**3/cv(n)
	 ag1(np1,1) = drdm(n,1)*g1p*(ahf*rl*v(i-1) - afac)
c
c	  ag1(np1,2) =-for*gor(np1)/r(np1) + pi4g*rho(n) +
c     $       drdm(n,2)*g1p*(ahf*rl*v(i-1) - afac)
c
c          Modified to agree with POLYWD 8/17/90
c
         f = (drz(np1)/r(np1))/(one+dlnr)
	 fdr = ( one + ahf*rl*(drz(np1)/r(np1)) )/(one+ dlnr)
         vl =-g1p*forpi*r(np1)**3/dm2(np1)
         ag1(np1,1) = vl*drdm(n,1)
         ag1(np1,2) =-for*gor(np1)/r(np1) + vl*drdm(n,2) -
     $                pi4g*rho(n)*ahf*fdr*f
         ag1(np1,3) = zero
c         ag3(np1,1) =-(rl+one)/(drz(np1)*drint(n) )
	 ag3(np1,1) =-(rl+one)/(one + dlnr )
         ag3(np1,2) = zero
c         ag4(np1,1) = rl1*gor(np1)/r(np1) +
c     $                rkl2(n)*g1p*(ahf*rl*v(i-1) - afac)
         ag4(np1,1) = rl1*gor(np1)/r(np1) + rkl2(n)*vl
         ag4(np1,2) = zero
c
c          The Poisson equation boundary condition is very messy.
c
         ddr11 =-drz(np1)/(drz(n)*(drz(n)+drz(np1)))
         ddr12 = (drz(np1)-drz(n))/(drz(n)*drz(np1))
         ddr13 = drz(n)/(drz(np1)*(drz(n)+drz(np1)))
         ap3n1 = r(n)**2/(drint(n)*drz(n))
         ap3n3 = r(np1)**2/(drint(n)*drz(np1))
c
         ap1(np1,1) =-pi4g*rho(n)*(rzone(n)**2*drdm(n,1)-
     $      ahf*dlog(rho(n)/rho(n-1))/dlog(rzone(n)/rzone(n-1)) )
         f = (drz(np1)/r(np1))/(one+dlnr)
         ap1(np1,2) =-pi4g*rho(n)*( (ddr13*two*rl*rzone(n) + ap3n3)*f +
     $      ahf*rho(n)*gor(n+1)*r(n+1)/p(n) + rzone(n)**2*drdm(n,2) )
c
	 ap3(np1,1) = ap3n1 + ddr11*two*rl*rzone(n)
c	  fp3 = (rl + ahf)*(drz(np1)/r(np1))
	 fp3 = (rl + ahf)*(drz(np1)/r(np1)) - 1
c	  ap3(np1,2) =-(ap3n1 + ap3n3) + ddr12*two*rl*rzone(n) -
c     $ 	(ddr13*two*rl*rzone(n) + ap3n3)*(rl+one)/(one+dlnr)
         ap3(np1,2) =-(ap3n1 + ap3n3) + ddr12*two*rl*rzone(n) -
     $	       (ddr13*two*rl*rzone(n) + ap3n3)*fp3/(one+dlnr)
         ap3(np1,3) = zero
         ap4(np1) =-pi4g*rl1*rho(n)
c
c          horizontal momentum equation.
c
         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i-1)/r(i-1)
         ah1(i,2) = drdm(i-1,2)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i)/r(i)
         ah3(i) = one
         ah4(i) = rkl2(i-1)*g1(i-1)*p(i-1)*v(i-1)
      else
c
c          horizontal momentum equation.
c
         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i-1)/r(i-1)
         ah1(i,2) = drdm(i-1,2)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i)/r(i)
         ah3(i) = one
         ah4(i) = rkl2(i-1)*g1(i-1)*p(i-1)*v(i-1)
c
c          radial momentum equation.
c
         afac = forpi*r(i)**3/dm2(i)
         ag1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*(ahf*rl*v(i-1)-afac)
         ag1(i,2)=-for*gor(i)/r(i)+pi4g*ahf*(rho(i-1)+rho(i))+
     $     drdm(i,1)*g1(i)*p(i)*(ahf*rl*v(i) + afac) +
     $     drdm(i-1,2)*g1(i-1)*p(i-1)*(ahf*rl*v(i-1) - afac )
         ag1(i,3) = drdm(i,2)*g1(i)*p(i)*(ahf*rl*v(i)+afac)
         ag3(i,1) = ahf*rl - r(i)/drz(i)
         ag3(i,2) = ahf*rl + r(i)/drz(i)
         ag4(i,1) = rkl2(i-1)*(g1(i-1)*p(i-1)*(ahf*rl*v(i-1) - afac)) +
     $     ahf*rl1*gor(i)/r(i)
        ag4(i,2) = rkl2(i)*(g1(i)*p(i)*(ahf*rl*v(i)+afac)) +
     $     ahf*rl1*gor(i)/r(i)
c
c          poisson equation.
c
         ddr11 =-drz(i)/(drz(i-1)*(drz(i-1)+drz(i)))
         ddr12 = (drz(i)-drz(i-1))/(drz(i-1)*drz(i))
         ddr13 = drz(i-1)/(drz(i)*(drz(i-1)+drz(i)))
         if( i .eq. 2 ) then
            ap1(i,1) =-pi4g*rho(i-1)*(rzone(i-1)**2*drdm(i-1,1)-
     $     ahf*dlog(rho(i-1)/rhoc)/dlog(rzone(i-1)/rz0) )
         else
            ap1(i,1) =-pi4g*rho(i-1)*(rzone(i-1)**2*drdm(i-1,1)-
     $     ahf*dlog(rho(i-1)/rho(i-2))/dlog(rzone(i-1)/rzone(i-2)) )
         endif
         ap1(i,2) =-pi4g*rho(i-1)*(rzone(i-1)**2*drdm(i-1,2)-
     $     ahf*dlog(rho(i)/rho(i-1))/dlog(rzone(i)/rzone(i-1)) )
         ap3i1 = r(i-1)**2/( drint(i-1)*drz(i-1))
         ap3i3 = r(i)**2/(drint(i)*drz(i))
         ap3(i,1) = ap3i1 +           ddr11*two*rl*rzone(i-1)
         ap3(i,2) =-(ap3i1 + ap3i3) + ddr12*two*rl*rzone(i-1)
         ap3(i,3) = ap3i3 +           ddr13*two*rl*rzone(i-1)
         ap4(i) =-pi4g*rl1*rho(i-1)
      endif
  20  continue
c
c        write the matrices to file 6.
c
      if( iout .ge. 2 ) then
         write(11,2200)
         do 27 i=1,np1
            write(11,2201) i,ag1(i,1),ag1(i,2),ag1(i,3),
     $         ag3(i,1),ag3(i,2),ag4(i,1),ag4(i,2),ah4(i)
  27     continue
         write(11,2300)
         do 28 i=1,np1
            write(11,2301) i,ah1(i,1),ah1(i,2),ah3(i),ap1(i,1),
     $         ap1(i,2),ap3(i,1),ap3(i,2),ap3(i,3),ap4(i)
  28     continue
      endif
      return
c
 2200 format(1h1,//,1x,'ag1(3), ag3(2), ag4(2) ,ah4(1)')
 2201 format(1x,i5,1p,3e12.4,2x,2e12.4,2x,2e12.4,2x,e12.4)
 2300 format(1h1,//,1x,'ah1(2), ah3(1), ap1(2), ap3(3), ap4(1)')
 2301 format(1x,i5,1p,2e12.4,2x,e12.4,2x,2e12.4,2x,3e12.4,2x,e12.4)
      end
