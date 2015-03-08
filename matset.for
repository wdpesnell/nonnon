      subroutine matset(nin,iout,lin,rkl2,z)
      implicit real*8(a-h,o-z)
      integer*2 n, np1, i
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
c          Modified to full nonadiabatic calculations summer 1989.
c       Uses the r**l transformation on all variables.
c
c                                                8/23/89 WD Pesnell
c
      parameter( onth=1.d0/3.d0 )
      parameter ( nmax=512 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $        dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
     $               bv(nmax)
      common/blkdrs/ drz(nmax),drint(nmax),gor(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk37/  drdm(nmax,2),dl1(nmax,2),dl2(nmax,2)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ag3(nmax,2),ag4(nmax,2),
     $               ah1(nmax,2),ah2(nmax),  ah3(nmax),  ah4(nmax),
     $               ap1(nmax,2),ap3(nmax,3),            ap4(nmax)
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/thermo/ rther,a,bk,avagd,ad3
      common ak1(nmax,4),ak2(nmax,3),ak4(nmax,3)
      dimension rkl2(nmax),z(nmax)
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
c      dkdr0 = dkdr0 + g3m10*dkdt0
      do 10 i=1,n
         rho(i) = one/v(i)
c         dkdr(i) = dkdr(i) + g3m1(i)*dkdt(i)
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
     $                  ahf*(rl+one)/rzone(i)**2)
         z(i) = rzone(i)/r(np1)
  10  continue
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      call pltdmp(z,nmax,nin, 'rz  ')
      gor(np1) = g*rm(np1)/r(np1)**2
      drz(np1) = p(n)*v(n)/gor(np1)
      pi4g = forpi*g
      rlumx = corlum
c
c          Initialize the matrices.
c
c          All matrices are prefixed with the letter A. The second
c       letter denotes the equation (G=radial component of the
c       momentum equation, H=horizontal component of the momentum
c       equation, P=Poisson equation and K=thermal equation)
c       in which each matrix is used. The number tells which
c       variable is to be acted on the the matrix: 1= radial
c       component of the motion eigenvector (mechanical
c       eigenvector), 2= entropy variations (to within a factor
c       of the temperature, actuallly the thermal eigenvector),
c       3= gravitation variations multiplied by the square root
c       of the interface radius differences), 4= horizontal
c       component of the mechanical eigenvector. It should be
c       noted that there is no AP2 matrix nor is AK3 present.
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
c         drdm01 = zero
         drdm02 =-ahf*(rl+one)/rz0**2
c         drdm02 =-(rl+one)/rz0**2
         ag11 = drdm01*g10*pc*(ahf*rl/rhoc - afac)
         ag1(1,1) = zero
         ag1(1,2) =-gor(1)/r(1) + drdm(1,1)*vl +
     $             drdm02*g10*pc*(ahf*rl/rhoc - afac)
         ag1(1,3) = drdm(1,2)*vl + ag11
         ag2(1,1) = g3m10*(ahf*rl - rhoc*afac)
         ag2(1,2) = g3m1(1)*(ahf*rl + rhoc*afac)
         ag3(i,1) = ahf*rl - r(i)/drz(i)
         ag3(i,2) = ahf*rl + r(i)/drz(i)
         ag4(1,1) = zero
         ag41 = g10*pc*(ahf*rl/rhoc - afac)*rl1/rz0**2
         ag4(i,1) = ag41 + ahf*rl1*gor(1)/r(1)
c         ag4(1,2) = rkl2(1)*vl + ag41 + rl1*gor(1)/r(1)
         ag4(i,2) = rkl2(i)*vl + ahf*rl1*gor(i)/r(i)
c
c          horizontal momentum equation.
c
c         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
c     $         ahf*gor(i-1)/r(i-1)
         ah1(i,1) = zero
         ah1(i,2) = (drdm01+drdm02)*g10*pc/rhoc + gor(i)/r(i)
         ah2(i) = g3m10
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
c
c          the thermal equation matrices.
c
         dtl1 = tc**4
         dtl2 = t(i)**4
         t4oki = dtl2/gkap(i)
         t4oki1 = dtl1/opac0
         diff = t4oki-t4oki1
         dl1(i,1) = zero
         dl1(i,2) = (for*g3m10-dkdr0)*(t4oki1/diff)
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rhoc)))
         dl2(i,1) = zero
         dl2(i,2) = ((for-dkdt0)*(t4oki1/diff))/(cv0*tc)
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rhoc)))
c
c          Keep track of the luminosity as a function of mass
c
            glom = zero
            rlumx1 = sorc0*cormas
            if( noburn .ne. 0 ) then
               ipass = i
               rlumx1 = cool(ipass,cormas/totmas)
            endif
            glom1 = frft(i)*(rlumx1/cormas)
     $         *(one+ahf*rl*cormas/(rhoc*forpi*rz0**3))
            dldm = rlumx1/cormas
            aconlm = (one-frft(i))*rlumx1
c            rlumx = rlumx1
c
c          convert the nuclear energy derivatives.
c
         dedr = (-dedv0+g3m10*dedt0)*sorc0
         deds = dedt0*sorc0/(cv0*tc)
           ak1(i,1) = zero
           ak1(i,2) = zero
           ak1(i,3) = - glom1*for/r(i)**2 -
     $        glom1*(dl1(i,1)*drdm02 + dl1(i,2)*drdm(i,1) )
     $        + dedr*drdm02
           ak1(i,4) =-glom1*dl1(i,2)*drdm(i,2)
           ak4(i,1) = zero
           ak4(i,2) =-rkl20*(glom1*dl1(i,1)-dedr)
           ak4(i,3) =-rkl2(i)*glom1*dl1(i,2)
           ak2(i,1) = zero
           ak2(i,2) =-glom1*dl2(i,1) + deds
           ak2(i,3) =-glom1*dl2(i,2)
c
c          add horizontal heat flow terms.
c
         rlom = zero
         rlom1 = frft(i)*rlumx1
         horiz = rkl20*for*ac3*t4oki1/rhoc**2
c         ak1(i,2) = ak1(i,2) - horiz*g3m10*drdm01 -
c     $       rkl20*ahf*rlom/(rhoc*forpi*rz0**3)
         ak1(i,3) = ak1(i,3) - horiz*g3m10*drdm02 -
     $       rkl20*ahf*rlom1/(rhoc*forpi*r(i)**3)
     $       - horiz*g3m10*drdm01 -
     $       rkl20*ahf*rlom/(rhoc*forpi*rz0**3)
         ak2(i,2) = ak2(i,2) - horiz/(cv0*tc)
         ak4(i,2) = ak4(i,2) + rkl20*(dldm-horiz*g3m10 -
     $      aconlm/(rhoc*forpi*rz0**3) )
      elseif( i .eq. np1 ) then
c
c          outer boundary conditions
c
         tau = gkap(n)*dm1(n)/(for*forpi*r(np1)**2)
         tau = 1.5d0*tau/(one+1.5d0*tau)
         fac = (r(np1)-rphoto)/r(np1)
         ef1 = (one+ahf*fac)/(one+fac)
         ef3 = two/(one+fac)
         dl1(np1,1) = ef1*(for*g3m1(n)-tau*dkdr(n))
c     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*rho(i-1)))
         dl1(np1,2) = dl1(np1,1)
         dl1(np1,2) = ef1*(for*g3m1(n)-tau*dkdr(n))
c     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*rho(i-1)))
         dl2(np1,1) = ef1*(for-tau*dkdt(n))/(t(n)*cv(n))
c     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*rho(i-1)))
         dl2(np1,2) = zero
         glom = rlumgv/dm1(n)
     $      *(one-ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
         glom1 = rlumgv/dm1(n)
     $      *(one+ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
         ak1(np1,1) = glom*dl1(n,1)*drdm(n-1,1)
         ak1(np1,2) =-glom1*dl1(np1,1)*drdm(n,1)+glom*for/r(n)**2 +
     $     glom*(dl1(n,1)*drdm(n-1,2)-dl1(n,2)*drdm(n,1) )
         ak1(np1,3) = glom*dl1(n,2)*drdm(n,2) -
     $      glom1*(dl1(np1,1)*drdm(n,2) + (ef3+tau*ef1)/r(np1)**2 )
         ak1(np1,4) = zero
         ak2(np1,1) = glom*dl2(n,1)
         ak2(np1,2) = glom*dl2(n,2) - glom1*dl2(np1,1)
         ak2(np1,3) = zero
         ak4(np1,1) = glom*rkl2(n-1)*dl1(n,1)
         ak4(np1,2) =-rkl2(n)*(glom1*dl1(np1,1) - glom*dl1(n,2))
         ak4(np1,3) = zero
c
c          Add horizontal heat flow terms. Note the absence of the
c       convective luminosity term. The outer boundary zone is 
c       assumed to be radiative, with no incoming luminosity.
c
         rlumx = rlumgv
         horiz = rkl2(n)*for*ac3*t(n)**4*v(n)**2/gkap(n)
         ak1(np1,2) = ak1(np1,2) - horiz*g3m1(n)*drdm(n,1) -
     $          rkl2(n)*ahf*frft(n)*rlumx*v(n)/(forpi*r(n)**3)
         ak1(np1,3) = ak1(np1,3) - horiz*g3m1(n)*drdm(n,2) -
     $          rkl2(n)*ahf*frft(n+1)*rlumx*v(n)/(forpi*r(n+1)**3)
         ak2(np1,2) = ak2(np1,2) - horiz/(cv(n)*t(n))
         ak4(np1,2) = ak4(np1,2) - rkl2(n)*horiz*g3m1(n)
c
c          outer boundary conditions.
c
         afac = forpi*r(np1)**3/dm2(np1)
         g1p = g1(n)*p(n)-for*ad3*t(n)**4*g3m1(n)
         g3m1ro = g3m1(n)/v(n)-for*ad3*t(n)**3/cv(n)
         ag1(np1,1) = drdm(n,1)*g1p*(ahf*rl*v(i-1) - afac)
         ag1(np1,2) =-for*gor(np1)/r(np1) + pi4g*rho(n) +
     $       drdm(n,2)*g1p*(ahf*rl*v(i-1) - afac)
         ag1(np1,3) = zero
         ag2(np1,1) = g3m1ro*(ahf*rl*v(i-1) - afac)
         ag2(np1,2) = zero
c         ag3(np1,1) =-(rl+one)/(drz(np1)*drint(n) )
         dlnr = (rl+ahf)*drz(np1)/r(np1)
         ag3(np1,1) =-(rl+one)/(one + dlnr )
         ag3(np1,2) = zero
         ag4(np1,1) = rl1*gor(np1)/r(np1) +
     $                rkl2(n)*g1p*(ahf*rl*v(i-1) - afac)
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
         ap3(np1,2) =-(ap3n1 + ap3n3) + ddr12*two*rl*rzone(n) -
     $         (ddr13*two*rl*rzone(n) + ap3n3)*(rl+one)/(one+dlnr)
         ap3(np1,3) = zero
         ap4(np1) =-pi4g*rl1*rho(n)
c
c          Horizontal momentum equation.
c
         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i-1)/r(i-1)
         ah1(i,2) = drdm(i-1,2)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i)/r(i)
         ah2(i) = g3m1(i-1)
         ah3(i) = one
         ah4(i) = rkl2(i-1)*g1(i-1)*p(i-1)*v(i-1)
      else
c
c          Horizontal momentum equation.
c
         ah1(i,1) = drdm(i-1,1)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i-1)/r(i-1)
         ah1(i,2) = drdm(i-1,2)*g1(i-1)*p(i-1)*v(i-1)+
     $         ahf*gor(i)/r(i)
         ah2(i) = g3m1(i-1)
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
         ag2(i,1) = g3m1(i-1)*(ahf*rl - rho(i-1)*afac)
         ag2(i,2) = g3m1(i)*(ahf*rl + rho(i)*afac)
         ag3(i,1) = ahf*rl - r(i)/drz(i)
         ag3(i,2) = ahf*rl + r(i)/drz(i)
         ag4(i,1) = rkl2(i-1)*(g1(i-1)*p(i-1)*(ahf*rl*v(i-1) - afac)) +
     $     ahf*rl1*gor(i)/r(i)
        ag4(i,2) = rkl2(i)*(g1(i)*p(i)*(ahf*rl*v(i)+afac)) +
     $     ahf*rl1*gor(i)/r(i)
c
c          Poisson equation.
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
c
c          Thermal equation matrices.
c
         dtl1 = t(i-1)**4
         dtl2 = t(i)**4
         t4oki = dtl2/gkap(i)
         t4oki1 = dtl1/gkap(i-1)
         diff = t4oki-t4oki1
         wiowi1 = dlog(dtl2/dtl1)
         wowsq = wiowi1**2
         gkogk1 = dlog(gkap(i)/gkap(i-1))
         denom = one-gkogk1/wiowi1
         dl1(i,1) = ((-for*g3m1(i-1)+dkdr(i-1))*(t4oki1/diff)
     $    + (-dkdr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom)
     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
         dl1(i,2) = ((for*g3m1(i)-dkdr(i))*(t4oki/diff)
     $    + (dkdr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom)
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
         dl2(i,1) = ((-for+dkdt(i-1))*(t4oki1/diff) -
     $    (dkdt(i-1)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i-1)*t(i-1))
     $      *(one-ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
         dl2(i,2) = ((for-dkdt(i))*(t4oki/diff) +
     $    (dkdt(i)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i)*t(i))
     $      *(one+ahf*rl*dm2(i)/(forpi*r(i)**3*ahf*(rho(i)+rho(i-1))))
c
c          Keep track of the luminosity as a function of mass
c
            rlumx = rlumx1
            glom = frft(i-1)*rlumx/dm1(i-1)
     $         *(one-ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
            rlumx1 = rlumx + sorce(i-1)*dm1(i-1)
            if( noburn .ne. 0 ) then
               ipass = i
               rlumx1 = cool(ipass,rm(ipass)/totmas)
            endif
            glom1 = frft(i)*rlumx1/dm1(i-1)
     $         *(one+ahf*rl*dm1(i-1)*v(i-1)/(forpi*rzone(i-1)**3))
            dldm = (rlumx1-rlumx)/dm1(i-1)
            aconlm = ahf*((one-frft(i-1))*rlumx+(one-frft(i))*rlumx1)
c
c          Convert the nuclear energy derivatives.
c
         dedr = (-dvsorc(i-1)+g3m1(i-1)*dtsorc(i-1))*sorce(i-1)
         deds = dtsorc(i-1)*sorce(i-1)/(cv(i-1)*t(i-1))
         if( i .eq. 2 ) then
            ak1(i,1) = zero
         else
            ak1(i,1) = glom*dl1(i-1,1)*drdm(i-2,1)
         endif
c         if( i .eq. 3 ) ak1(i,1) = zero
         ak1(i,2) =-glom1*dl1(i,1)*drdm(i-1,1)+glom*for/r(i-1)**2 +
     $     glom*(dl1(i-1,1)*drdm(i-2,2)+dl1(i-1,2)*drdm(i-1,1) )
     $     +  dedr*drdm(i-1,1)
         ak1(i,3) = glom*dl1(i-1,2)*drdm(i-1,2) - glom1*for/r(i)**2 -
     $     glom1*(dl1(i,1)*drdm(i-1,2)+dl1(i,2)*drdm(i,1) )
     $     + dedr*drdm(i-1,2)
         ak1(i,4) =-glom1*dl1(i,2)*drdm(i,2)
         if( i .eq. 2 ) then
            ak4(i,1) = zero
         else
            ak4(i,1) = glom*rkl2(i-2)*dl1(i-1,1)
         endif
         ak4(i,2) =-rkl2(i-1)*(glom1*dl1(i,1) -  dedr
     $     - glom*dl1(i-1,2) )
         ak4(i,3) =-glom1*rkl2(i)*dl1(i,2)
         ak2(i,1) = glom*dl2(i-1,1)
         ak2(i,2) = glom*dl2(i-1,2)-glom1*dl2(i,1) + deds
         ak2(i,3) =-glom1*dl2(i,2)
c
c          add horizontal heat flow terms.
c
         rlom = frft(i-1)*rlumx
         rlom1 = frft(i)*rlumx1
         horiz = rkl2(i-1)*for*ac3*t4oki1*v(i-1)**2
         ak1(i,2) = ak1(i,2) - horiz*g3m1(i-1)*drdm(i-1,1) -
     $       rkl2(i-1)*ahf*rlom*v(i-1)/(forpi*r(i-1)**3)
         ak1(i,3) = ak1(i,3) - horiz*g3m1(i-1)*drdm(i-1,2) -
     $       rkl2(i-1)*ahf*rlom1*v(i-1)/(forpi*r(i)**3)
         ak2(i,2) = ak2(i,2) - horiz/(cv(i-1)*t(i-1))
         ak4(i,2) = ak4(i,2) + rkl2(i-1)*(dldm - horiz*g3m1(i-1) -
     $      aconlm*(v(i-1)/(forpi*rzone(i-1)**3)) )
      endif
  20  continue
c
c        Write the matrices to file 11 (NONOUT).
c
      if( iout .ge. 2 ) then
         write(11,2200)
         do 27 i=1,np1
            write(11,2201) i,ag1(i,1),ag1(i,2),ag1(i,3),ag2(i,1),
     $         ag2(i,2),ag3(i,1),ag3(i,2),ag4(i,1),ag4(i,2),ah4(i)
  27     continue
         write(11,2300)
         do 28 i=1,np1
            write(11,2201) i,ah1(i,1),ah1(i,2),ah2(i),ah3(i),ap1(i,1),
     $         ap1(i,2),ap3(i,1),ap3(i,2),ap3(i,3),ap4(i)
  28     continue
         write(11,2600)
         do 26 i=1,np1
            write(11,2201) i,ak1(i,1),ak1(i,2),ak1(i,3),ak1(i,4),
     $          ak2(i,1),ak2(i,2),ak2(i,3),ak4(i,1),ak4(i,2),ak4(i,3)
  26     continue
      endif
      return
c
 2200 format(1h1,//,20h AG1,AG2,AG3,AG4,AH4)
 2201 format(1x,i5,1p,10e12.4)
 2300 format(1h1,//,24h AH1,AH2,AH3,AP1,AP3,AP4)
 2600 format(1h1,//,10x,36h The pulsation matrices, AK1,AK2,AK4,/)
      end
