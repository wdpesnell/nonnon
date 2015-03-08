      subroutine lnanon(nin,iout,lin,nomgin,omlow,omhigh)
      implicit real*8(a-h,o-z)
      complex*16 fomega,cx,comega,dcomega, comeg1,comeg2,cerr,cerra,
     $      cdh,dror,
     $      drho,dtemp,dent,dlum,dpress,dgam,cy,cz,rcont,hcont,domsq,
     $      czero,ci,dlth,dlmech
      integer*2 icount, iroot, nomsq, nroot
      logical iout0, iout1
c
      structure /puls_sum/
         integer*2 mass_flag
         integer*4 ell, k, n_p, n_g
         real*4 omega_ad, omsq_n, period, acc_ad_w, acc_ad_e
         complex*8 del_om_qa
         complex*8 comega
         real*4 kappa, acc_na_k
      end structure
      record /puls_sum/ nonradial
c
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
      parameter ( accur=1.d-10, czero=(0.0d0,0.0d0), ci=(0.0d0,1.0d0) )
      parameter ( onth=1.d0/3.d0 )
      parameter ( nmax=512, nmax3=3*nmax, nmax4=4*nmax )
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
     $               ap1(nmax,2),ap3(nmax,3),           ap4(nmax)
      common/linear/ cdh(nmax),dror(nmax),drho(nmax),dtemp(nmax),
     $      dent(nmax),dlum(nmax),dpress(nmax),dgam(nmax),cy(nmax),
     $      cz(nmax),dlmech(nmax),dlth(nmax),work(nmax),dr(nmax),
     $      dh(nmax),gam(nmax),y(nmax),z(nmax),dp(nmax),dl(nmax),
     $      adrho(nmax),weight(nmax),cvt(nmax),wtherm(nmax),
     $      wgrav(nmax),wdiag(nmax),wcross(nmax),dtph(nmax),dpph(nmax),
     $      dlph(nmax),drhoph(nmax),rkepl(nmax),stwait(nmax),
     $      spac(nmax,59)
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
c
c      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
c      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
c      common/thermo/ rther,a,bk,avagd,ad3
c
      common ak1(nmax,4),ak2(nmax,3),ak4(nmax,3)
      dimension rkl2(nmax)
      dimension cx(nmax4),x(nmax3),stomeg(50),stomsq(50)
      equivalence (cx(1),x(1))
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
      rz0 = ahf**(1./3.)*r(1)
      rkl20 = rl1/rz0**2
c         drdm01 = zero
c         drdm02 =-(rl+one)/rz0**2
         drdm01 =-ahf*(rl+one)/rz0**2
         drdm02 =-ahf*(rl+one)/rz0**2
c
c          Convergence to a normal mode is governed by the momentum
c       equation for the outermost zone. As in the radial case, this
c       function (divided by any x value) is used in a secant method
c       for controlling the iterations. As the equations are not
c       nonlinear (as in NONNON) in the frequency, it is hoped the
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
         fomeg = fomsq(omsq,x,xnorm,nmax,n)
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
            err = fomsq(omsq,x,xnorm,nmax,n)
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
c         call invint(omsq,x,finerr,xnorm,nmax,n)
c
c          Do not repeat a root.
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
               dl(i) = dl1(i,2)*adrho(i)+
     $            for*dr(i) - rl1*dh(i)
            else
               dh(i)=(x(3*i-1)/rzone(i-1)**2)*(rzone(i-1)/r(np1))**lval
               adrho(i) = (drdm(i-1,1)*x(3*i-3)+drdm(i-1,2)*x(3*i) +
     $           rl1*x(3*i-1)/rzone(i-1)**2)*(rzone(i-1)/r(np1))**lval
               dp(i) = g1(i-1)*adrho(i)
               gam(i) = x(3*i-2)*(rzone(i-1)/r(np1))**lval
               dl(i) = dl1(i,1)*adrho(i-1)+dl1(i,2)*adrho(i)+
     $            for*dr(i) - ahf*rl1*(dh(i)+dh(i-1))
            endif
  40     continue
         dlnr = ahf*(rl+one)*drz(np1)/r(np1)
         gam(n+2) = gam(n)*(one-dlnr)/(one+dlnr)
         dh(n+2) = zero
         qdr0 = dr(2)*(r(1)/r(2))**(lval-2)
c
c          weight function for the adiabatic oscillation.
c
         wtherm(1) = (g10*pc/rhoc)*cormas*adrho(1)**2
         wgrav(1) = dr(1)*dm2(1)*rl*gam(1)
         wcross(1) = rl1*gor(1)*r(1)*dr(1)*dh(1)*dm2(1)
         wdiag(1) = (pi4g*rho(1)-gor(1)*for/r(1) )*
     $                   (dr(1)*r(1))**2*dm2(1)
         weight(1) = wtherm(1) + wgrav(1) + wcross(1) + wdiag(1)
         rgwait = weight(1)
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
  48     continue
         dl(np1) = for-ahf*rl1*dh(n)+dl1(np1,1)*adrho(n)
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
            weight(i) = weight(i)/(rke*omsq)
            wtherm(i) = wtherm(i)/(rke*omsq)
            wgrav(i) = wgrav(i)/(rke*omsq)
            wdiag(i) = wdiag(i)/(rke*omsq)
            wcross(i) = wcross(i)/(rke*omsq)
  46     continue
         rgwait = (rgwait/rke)
         qch = dabs((omsq-rgwait)/omsq)
c
c          Evaluate the rotational spitting coefficients.
c
         crsum = rl1*crsum/rke
         hrsum = hsum/rke
         ckl = crsum + hrsum
         index = knode(dr,dh,nmax,np1,np,ng)
         omega = dsqrt(omsq)
         period = twopi/omega
         call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax,itrans)
         pdays = period/86400.d0
c
         nonradial.mass_flag = nint(10.d0*totmas/1.991d33)
         nonradial.ell       = lval
         nonradial.k         = index
         nonradial.n_p       = n_p
         nonradial.n_g       = n_g
         nonradial.omega_ad  = omega
         nonradial.omsq_n    = omsq
         nonradial.period    = period
         nonradial.acc_ad_w  = qch
c
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
     $             dp(i),dl(i),weight(i),cvt(i),i=1,np1)
         endif
	 call pltdmp(dr,nmax,np1,'dr/r')
	 call pltdmp(dh,nmax,np1,'dh/h')
	 call pltdmp(adrho,nmax,np1,'drho')
	 call pltdmp(dl,nmax,np1,'dl/l')
	 call pltdmp(rkepl,nmax,np1,'ke  ')
	 call pltdmp(weight,nmax,np1,'wait')
	 call pltdmp(wtherm,nmax,np1,'wthr')
	 call pltdmp(wgrav,nmax,np1,'wgrv')
	 call pltdmp(wdiag,nmax,np1,'wdia')
	 call pltdmp(wcross,nmax,np1,'wcrs')
c
         call eigenf(nmax,omsq,n,lval,rke,gor, ech)
         nonradial.acc_ad_e  = ech
c
c          Non-adiabatic loop.
c
c          Calculate the quasi-adiabatic stability coefficient.
c      This is essentially the same method as in the radial case.
c
      y(1) =-ak4(1,2)*x(2)-ak4(1,3)*x(5)
      ido = n
      do 60 i=2,ido
         y(i) =-ak4(i,1)*x(3*i-4)-ak4(i,2)*x(3*i-1)-ak4(i,3)*x(3*i+2)
  60  continue
      y(np1) =-ak4(np1,1)*x(3*np1-4)-ak4(np1,2)*x(3*np1-1)
c
c          radial contribution
c
      do 62 i=1,np1
         cz(i) = czero
         jminx = max0(1,4-i)
         jmaxx = 4
         if( i .eq. np1 ) jmaxx = 3
         sy = zero
         do 620 j=jminx,jmaxx
            sy = sy - ak1(i,j)*x(3*(i-3 + j))
 620     continue
         y(i) = y(i) + sy
  62  continue
c
      call trcsol(ak2,nmax,1,np1,ci*omega,y,cz,dlth,dlmech)
      cz(n+2) = czero
      rcont = x(3)*(ag2(1,1)*cz(1)+ag2(1,2)*cz(2))*dm2(1)
     $            *(r(1)/r(np1))**(2*lval)
      hcont = cz(1)*dh(1)*ah2(1)*cormas*(rz0/r(np1))**lval
      do 63 i=1,n
         hcont = hcont + cz(i+1)*dh(i)*ah2(i)*dm1(i)
     $            *(rzone(i)/r(np1))**lval
         rcont = rcont + (ag2(i+1,1)*cz(i+1)+ag2(i+1,2)*cz(i+2))*
     $            x(3*(i+1))*dm2(i+1)*(r(i+1)/r(np1))**(2*lval-1)
  63  continue
      rcont = rcont/xnorm
      domsq = (rcont+rl1*hcont)/rke
      comega = cdsqrt(omsq + domsq)
      growth =-forpi*dimag(comega)/omega
c
c          Limit the quasi-adiabatic growth rate to 0.1
c
      if( abs(growth) .gt. 0.1d0 ) then
         omega_im =-dsign(0.1d0, growth)*omega/forpi
         comega = dcmplx( omega, omega_im)
         growth =-forpi*dimag(comega)/dreal(comega)
      endif
      write(1,6000) comega,growth
      write(1,6001) domsq,rke,rcont,hcont,rsum,hsum
      write(6,6000) comega,growth
      write(11,6000) comega,growth
      write(11,6001) domsq,rke,rcont,hcont,rsum,hsum
      nonradial.del_om_qa = domsq
c
c      if( cdabs(comega) .gt. 1.1d0*omega ) comega = omega*(1.d0,1.d-4)
         do 55 icount=1,30
            cerr = fomega(comega,cx,xnorm,nmax,n)
            if( icount .eq. 1 ) then
               comeg1 = comega
               cerra = cerr
               comega = comega*(one+1.d-7)
            else
               afac = cdabs((comega-comeg1)/comega)
               if( iout0 ) write(11,4400) icount,comega,afac
               if( afac .le. accur ) goto 57
               comeg2 = (cerra*comega-cerr*comeg1)/(cerra-cerr)
               comeg1 = comega
               cerra = cerr
               dcomega = comeg2 - comega
               if( cdabs(dcomega)/cdabs(comega) .gt. 0.05d0 ) then
                  dcomega = dcomega*0.05d0*cdabs(comega)/cdabs(dcomega)
               endif
c               comega = comeg2
               comega = comega + dcomega
            endif
  55     continue
      write(11,4500) nroot,comega
      goto 100
  57  continue
c
c          calculation of the non-adiabatic eigenfunctions.
c
      wgravr = zero
      wgravh = zero
      dlum(1) = czero
      dlth(1) = czero
      dlmech(1) = czero
      dror(1) = czero
      cx(4*np1) = xnorm
      rsum = zero
      hsum = zero
      do 70 i=1,np1
         dror(i) = (cx(4*i)/r(i)**2)*(r(i)/r(np1))**lval
         if( i .eq. 1 ) then
c            write(6,*) tc,cv0,cormas
            cdh(i) = cx(4*i-1)/rz0**2*(rz0/r(np1))**lval
            dgam(i) = cx(4*i-3)*(rz0/r(np1))**lval
            dent(i) = cx(4*i-2)/(tc*cv0)*(rz0/r(np1))**lval
            drho(i) = ((drdm01+drdm02)*cx(4*i)+
     $         rkl20*cx(4*i-1))*(rz0/r(np1))**lval
            dpress(i) = g10*drho(i)+g3m10*cx(4*i-2)*rhoc/pc
            dtemp(i) = g3m10*drho(i)+dent(i)
c
            rsum = conjg(dror(i))*dror(i)*r(i)**2*dm2(i)
            hsum = conjg(cdh(i))*cdh(i)*rz0**2*cormas
            work(i) = dimag(dconjg(drho(i))*dent(i))*tc*cv0*g3m10*cormas
            wgravr = wgravr + dimag(dconjg(dror(i))*(r(i)/r(np1))**lval*
     $        (ag3(i,2)*cx(4*i-3)) )*dm2(i)
            wgravh = wgravh+dimag(dconjg(cdh(i))*dgam(i))*cormas
c
c          various luminosity contributions
c
            dlth(i) = dl2(i,2)*cx(4*i-2)
            dlmech(i) = dl1(i,2)*drho(i) +
     $          for*dror(i) - rl1*cdh(i)
            dlum(i) = dlth(i) + dlmech(i)
         else
            cdh(i) = cx(4*i-1)/rzone(i-1)**2*(rzone(i-1)/r(np1))**lval
            dgam(i) = cx(4*i-3)*(rzone(i-1)/r(np1))**lval
            dent(i) = cx(4*i-2)/(t(i-1)*cv(i-1))*
     $                 (rzone(i-1)/r(np1))**lval
            drho(i) = (drdm(i-1,1)*cx(4*i-4)+drdm(i-1,2)*cx(4*i)+
     $         rkl2(i-1)*cx(4*i-1))*(rzone(i-1)/r(np1))**lval
            dpress(i) = g1(i-1)*drho(i-1)+
     $                   g3m1(i-1)*cx(4*i-2)/(p(i-1)*v(i-1))
            dtemp(i) = g3m1(i-1)*drho(i)+dent(i)
c
            rsum = rsum + dconjg(dror(i))*dror(i)*r(i)**2*dm2(i)
            hsum = hsum + dconjg(cdh(i))*cdh(i)*rzone(i-1)**2*dm1(i-1)
            work(i) = dimag(dconjg(drho(i))*dent(i))*
     $                t(i-1)*cv(i-1)*g3m1(i-1)*dm1(i-1)
            wgravr = wgravr + dimag(dconjg(dror(i))*(r(i)/r(np1))**lval*
     $        (ag3(i,1)*cx(4*i-7)+ag3(i,2)*cx(4*i-3)) )*dm2(i)
            wgravh = wgravh+dimag(dconjg(cdh(i))*dgam(i))*dm1(i-1)
c
c          various luminosity contributions
c
            dlth(i) = dl2(i,1)*cx(4*i-6) + dl2(i,2)*cx(4*i-2)
            dlmech(i) = dl1(i,1)*drho(i-1) + dl1(i,2)*drho(i) +
     $          for*dror(i) - ahf*rl1*(cdh(i)+cdh(i-1))
            dlum(i) = dlth(i) + dlmech(i)
         endif
         call amppha(dlum(i),dum,dlph(i))
         call amppha(dtemp(i),dum,dtph(i))
         call amppha(dpress(i),dum,dpph(i))
         call amppha(drho(i),dum,drhoph(i))
  70  continue
      dlth(np1) = dl2(np1,1)*cx(4*n-2)
      dlmech(np1) = for-ahf*rl1*cdh(n)+dl1(np1,1)*drho(n)
      dlum(np1) = dlth(np1) + dlmech(np1)
      call amppha(dlum(np1),dum,dlph(np1))
c
c          write the eigenvectors to the plot file.
c
      call pltdmp(drhoph,nmax,np1,'drph')
      call pltdmp(dpph,nmax,np1,'dpph')
      call pltdmp(dtph,nmax,np1,'dtph')
      call pltdmp(dlph,nmax,np1,'dlph')
      call cptdmp(dror,nmax,np1,'dr/r')
      call cptdmp(cdh,nmax,np1,'dh/h')
      call cptdmp(drho,nmax,np1,'drho')
      call cptdmp(dent,nmax,np1,'dscv')
      call cptdmp(dlum,nmax,np1,'dl/l')
c
c          calculate the integrated work function.
c
      omega = dreal(comega)
      growth =-forpi*dimag(comega)/omega
      rke = ahf*omega*omega*(rsum + rl1*hsum)
      y(n+2) = zero
      do 75 ii=1,np1
         i = np1 + 1 - ii
         y(i) = y(i+1) - work(i)*pi/rke
         work(i) =-work(i)*(pi/(rke*dabs(growth)))
  75  continue
      call pltdmp(work,nmax,np1,'work')
      call pltdmp(y,nmax,np1,'wint')
      wgravr = wgravr*pi/rke
      wgravh = wgravh*pi/rke
c
      write(11,7000) lval
      period = twopi/omega
      pdays  = period/86400.d0
      qdays  = pdays*dsqrt(rhobar/1.41d0)
      wrkint = y(1)
      qch = dabs((growth-wrkint)/growth)
      write(1,7001) comega,growth,y(1),qch
      write(1,7003) period,pdays,qdays
      write(6,7001) comega,growth,y(1),qch
      write(6,7003) period,pdays,qdays
      write(11,7001) comega,growth,y(1),qch
      write(11,7003) period,pdays,qdays
c
c          sum and display the work integral contributions.
c
      qch = dabs((wgravr+rl1*wgravh)/y(1))
      write(1,7002) rsum,hsum,rke
      write(1,7012) y(1),wgravr,wgravh,qch
      write(11,7002) rsum,hsum,rke
      write(11,7012) y(1),wgravr,wgravh,qch
      nonradial.comega   = comega
      nonradial.kappa    = growth
      nonradial.acc_na_k = qch
      if( iout1 ) then
         write(11,7010)
         write(11,7011) (i,dror(i),cdh(i),dgam(i),dent(i),drho(i),
     $                    i=1,np1)
         write(11,7013)
         write(11,7014) (i,dpress(i),dlum(i),dlth(i),dlmech(i),
     $                    work(i),y(i),i=1,np1)
      endif
c
c          do various and non-obscene things to the work functions.
c
      call workmx(y,np1,nmax)
      call zie(work,np1,nmax)
      call summary( nonradial )
c
 100  continue
      return
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
 5501 format(1x,i5,1p,8e15.6)
 5502 format(/4x,1hi,8x,4hdr/r,10x,4hdh/h,10x,8hdrho/rho,8x,5hgamma,
     $   11x,4hdp/p,11x,4hdl/l,10x,6hweight,10x,4hcvtm/)
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
c          quasi-adiabatic formats.
c
 6000 format(/,1x,27hquasi-adiabatic eigenvalue:,1p,2e11.3,
     $   8h growth=,e11.4)
 6001 format(7h domsq=,1p,2e11.4,5h rke=,e11.4,/,7h rcont=,2e11.4,
     $   7h hcont=,2e11.4,6h rsum=,e11.4,6h hsum=,e11.4)
c
c          non-adiabatic format statements.
c
 7000 format(1h1,//,31h non-adiabatic solution for l =,i3)
 7001 format(//,21h frequency (sec**-1)=,1pe15.6,2x,12h imag. part=,
     $   e15.6,/,9x,24h predicted growth rate =,e14.6,/,
     $   1x,32hgrowth rate from work integral =,e14.6,7h error=,e11.3)
 7002 format(/,23h work integral results:,/,18h kinetic energies:,
     $   1pe15.6,9h (radial),e15.6,9h (horiz.),e15.6,
     $   22h (total), all in ergs.)
 7003 format(/,9h periods:,1pe13.6,7h (secs),e13.6,12h (days), q =,
     $   e11.4)
 7010 format(//3x,1hi,10x,4hdr/r,18x,4hdh/h,17x,5hgamma,15x,
     $   11hd (entropy),13x,8hdrho/rho,8x,4hwork/)
 7011 format(1x,i4,1p,10e11.3 )
 7012 format(1x,29h work integral contributions:,1p,3e12.4,/,3x,
     $   27herror in poisson solution =,e12.4)
 7013 format(//,3x,1hi,9x,6hdpress,18x,4hdlum,14x,13hdlum(thermal),
     $   10x,11hdlum(mech.),14x,4hwork,/6x,4(3x,4hreal,7x,5himag.,3x),
     $   3x,4hzone,6x,8hintegral /)
 7014 format(1x,i4,1p,10e11.3)
 9000 format(1x,44h from lnanon...no adiabatic roots found, l =,i3,/,
     $   1x,17h frequency range:,1p,2e15.6)
      end
      function knode(x,y,nmax,nin,nps,ngs)
      implicit real*8(a-h,o-z)
      integer*2 i, n
c
c          Returns the number of zero-crossings in the counter-
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
 1000 format(1x,51h from knode...quadrant jump default rotations used.)
 1001 format(2x,51hfrom knode...node at outer boundary is not counted.)
      end
