      subroutine cmplna(io,ihmin,ihmax,nin)
      implicit real*8(a-h,o-z)
      integer*2 i, ii, j, icount
c      logical first /.true./
c
c     *****************************************************************
c      Pulsation Equation Solver
c      Basic references -- Castor, Ap. J. 166 (109) 1971.
c                          Pesnell, 1983, Thesis, Univ. of Florida.
c                          Pesnell, PASP, 99 (975), 1987.
c     *****************************************************************
c
c          This version of the linear nonadiabatic radial stability 
c       analysis is suitable for envelope models generated with the 
c       MODEL building code. The inner boundary is valid for a complete 
c       model and models of low central condensation. If the values of 
c       the central ball are not present, or inconsistent, the results 
c       of the stability analysis are questionable.
c
c
c       06/01/94: Moved the conversion of DKDT to DKDS to main routine.
c
      complex*16 x0,y0,y1,z,xs,dror,dlum,omsq1,comega,drho,dpress,
     $  dtemp,amat,oma,erra,err,omb,czero,ci,czi,dlth,dlmech
      parameter ( czero=(0.d0,0.d0), ci=(0.d0,1.d0), accur=1.d-10 )
      type puls_sum
         integer*2 mass_flag
         integer*4 ell, k, n_p, n_g
         real*4 omega_ad, omsq_n, period, acc_ad_w, acc_ad_e
         complex*8 del_om_qa
         complex*8 comega
         real*4 kappa, acc_na_k
      end type
      type(puls_sum) radial
c
c        arrays are rearranged from the original lna code
c
      parameter ( nmax=512 )
      parameter( nmax2=2*nmax)
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	     dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
     $		     drz(nmax),drint(nmax),gor(nmax)
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ak1(nmax,4),ak2(nmax,3),
     $	     bl2(nmax,2),bl1(nmax,3),dr(nmax,2)
      common xs(nmax2),y0(nmax),y1(nmax),x0(nmax),dror(nmax),
     $	 dlth(nmax),dlmech(nmax),dlum(nmax),drho(nmax),dtemp(nmax),
     $   dpress(nmax), xamp(nmax),
     $	 xpha(nmax),rhoamp(nmax),rhopha(nmax),pramp(nmax),prpha(nmax),
     $	 temamp(nmax),tempha(nmax),glamp(nmax),glpha(nmax),
     $	 w(nmax),cvt(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk37/  yo(nmax),xo(nmax),work(nmax),dro(nmax),dp(nmax),
     $	     dt(nmax),dl(nmax)
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      dimension amat(nmax2,7), epstwt(nmax)
      dimension zzero(nmax,36),szero(nmax,19)
      equivalence (zzero(1,1),xs(1)),(szero(1,1),ag1(1,1))
      equivalence (amat(1,1),y0(1))
c
c   Zero common blocks scrtch and (blank).
c
      do 5 j=1,36
         do 5 i=1,nmax
            zzero(i,j) = zero
   5  continue
      do 6 j=1,19
         do 6 i=1,nmax
            szero(i,j) = zero
   6  continue
c
      n = nin
      np1 = n+1
      nm1 = n-1
c
      rlumx = corlum
      afac = forpi*r(1)**2/dsqrt(dm2(1))
      do 60 i=1,np1
         afacp = afac
         if( i .eq. np1 ) goto 40
         afac = forpi*r(i+1)**2/dsqrt(dm2(i+1))
c
c          Change dln(kappa)/dln(rho)<t> to dln(kappa)/dln(rho)<s>
c
c         if( first ) then
c            dkdr(i) = dkdr(i)+g3m1(i)*dkdt(i)
c            first = .false.
c         endif
         dr(i,1) = afacp/v(i)/dm1(i)
         dr(i,2) =-afac/v(i)/dm1(i)
c
c          The inner boundary condition is zero gradient, not zero
c       motion.
c
         if( i .eq. 1 ) then
            dtl1 = tc**4
            dtl2 = t(i)**4
            t4oki = dtl2
            t4oki1 = dtl1
            diff = t4oki-t4oki1
            dkdr0 = dkdr0 + g3m10*dkdt0
            bl1(i,1) = zero
            bl1(i,2) = ((for*g3m10-dkdr0)*(t4oki1/diff)*thre
     $             + for)/(r(i)*dsqrt(dm2(i))) +
     $     dr(i,1)*(for*g3m1(i)-dkdr(i))*(t4oki/diff)
            bl1(i,3) = dr(i,2)*(for*g3m1(i)-dkdr(i))*(t4oki/diff)
            bl2(i,1) = (-for+dkdt(i))*(t4oki1/diff)/(cv0*tc)
            bl2(i,2) = (for-dkdt(i))*(t4oki/diff)/(cv(i)*t(i))
            drdm0 =-thre/(r(1)*dsqrt(dm2(1)))
            qafacp = forpi*r(1)**2/dsqrt(dm2(1))
            qafac = forpi*r(2)**2/dsqrt(dm2(2))
            ag1(1,1) = zero
            ag1(1,2) =-for*g*rm(1)/r(1)**3 + qafacp*(g1(1)*p(1)*dr(1,1)
     $               - drdm0*pc*g10 )
            ag1(1,3) = qafacp*g1(1)*p(1)*dr(1,2)
            ag2(i,1) =-qafacp*g3m10*rhoc
            ag2(i,2) = qafacp*g3m1(1)/v(1)
           goto 50
        endif
c
c          set up the mechanical matrices.
c
         g1p = g1(i)*p(i)
         g3m1ov = g3m1(i)/v(i)
c         if( i .eq. n ) then
c            g1p = g1p-1.00855d-14*t(i)**4*g3m1(i)
c            g3m1ov = g3m1ov-1.00855d-14*t(i)**3/cv(i)
c         endif
         ag1(i,1) =-afacp*g1(i-1)*p(i-1)*dr(i-1,1)
         ag1(i,2) =-for*g*rm(i)/r(i)**3 + afacp*(g1p*dr(i,1)-
     $               g1(i-1)*p(i-1)*dr(i-1,2) )
         ag1(i,3) = afacp*g1p*dr(i,2)
         ag2(i,1) =-afacp*g3m1(i-1)/v(i-1)
         ag2(i,2) = afacp*g3m1ov
c
c          initialize the thermal matrices, if irad=1 use the
c       stellingwerf interpolation formula.
c
         dtl1 = t(i-1)**4
         dtl2 = t(i)**4
         t4oki = dtl2/gkap(i)
         t4oki1 = dtl1/gkap(i-1)
         diff = t4oki-t4oki1
         wiowi1 = dlog(dtl2/dtl1)
         wowsq = wiowi1**2
         gkogk1 = dlog(gkap(i)/gkap(i-1))
         denom = one - gkogk1/wiowi1
       bl1(i,1) = dr(i-1,1)*((-for*g3m1(i-1)+dkdr(i-1))*(t4oki1/diff)
     $  + (-dkdr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom)
c        if( i .eq. 2 ) bl1(2,1) = zero
         bl1(i,2) = for/(r(i)*dsqrt(dm2(i))) +
     $     dr(i,1)*((for*g3m1(i)-dkdr(i))*(t4oki/diff) +
     $    (dkdr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom ) +
     $  dr(i-1,2)*((-for*g3m1(i-1)+dkdr(i-1))*(t4oki1/diff) +
     $    (-dkdr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom )
         bl1(i,3) = dr(i,2)*((for*g3m1(i)-dkdr(i))*(t4oki/diff)
     $     + (dkdr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom)
c
         bl2(i,1) = ((-for+dkdt(i-1))*(t4oki1/diff)
     $  - (dkdt(i-1)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i-1)*t(i-1))
         bl2(i,2) = ((for-dkdt(i))*(t4oki/diff)
     $  + (dkdt(i)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i)*t(i))
         goto 50
  40    continue
         tau = gkap(n)*dm1(n)/(for*forpi*r(np1)*r(np1))
         tau = 1.5d0*tau/(one+1.5d0*tau)
         bbb = (r(np1)-rphoto)/r(np1)
         ef1 = (one+ahf*bbb)/(one+bbb)
         ef3 = two/(one+bbb)
         bl1(np1,1) = ef1*dr(n,1)*( for*g3m1(n) - tau*(one+dkdr(n)) )
         bl1(np1,2) = (ef3-tau*ef1*two)/(r(np1)*dsqrt(dm2(np1))) +
     $      dr(n,2)*ef1*(for*g3m1(n) - tau*(one+dkdr(n)))
         bl1(np1,3) = zero
         bl2(np1,1) = ef1*(for - tau*dkdt(n))/(cv(n)*t(n))
         bl2(np1,2) = zero
         g1p = g1(n)*p(n)
         g3m1ov = g3m1(n)/v(n)
c         g1p = g1p-1.00855d-14*t(n)**4*g3m1(n)
c         g3m1ov = g3m1ov-1.00855d-14*t(n)**3/cv(n)
         ag1(np1,1) =-afacp*g1p*dr(n,1)
         ag1(np1,2) =-for*g*rm(np1)/r(np1)**3-afacp*g1p*dr(n,2)
         ag1(np1,3) = zero
         ag2(np1,1) =-afacp*g3m1ov
         ag2(np1,2) = zero
  50    continue
         if( i .eq. 1 ) then
            glom1 = frft(i)*sorc0
            dedr =-(dedv0-dedt0*g3m10)*sorc0
            deds = dedt0*sorc0/(tc*cv0)
            ak1(i,1) = zero
            ak1(i,2) = zero
            ak1(i,3) =-glom1*bl1(i,2)+dedr*drdm0
            ak1(i,4) =-glom1*bl1(i,3)
            ak2(i,1) = zero
            ak2(i,2) =-glom1*bl2(i,1) + deds
            ak2(i,3) =-glom1*bl2(i,2)
         elseif( i .eq. np1 ) then
            glom = frft(i-1)*rlumx/dm1(i-1)
            glom1 = rlumx/dm1(i-1)
            ak1(i,1) = glom*bl1(i-1,1)
            ak1(i,2) = glom*bl1(i-1,2)-glom1*bl1(i,1)
            ak1(i,3) = glom*bl1(i-1,3)-glom1*bl1(i,2)
            ak1(i,4) = zero
            ak2(i,1) = glom*bl2(i-1,1)
            ak2(i,2) = glom*bl2(i-1,2)-glom1*bl2(i,1)
            ak2(i,3) = zero
         else
            glom = frft(i-1)*rlumx/dm1(i-1)
            rlumx = rlumx + sorce(i-1)*dm1(i-1)
            glom1 = frft(i)*rlumx/dm1(i-1)
            dedr =-(dvsorc(i-1)-dtsorc(i-1)*g3m1(i-1))*sorce(i-1)
            deds = dtsorc(i-1)*sorce(i-1)/(t(i-1)*cv(i-1))
            ak1(i,1) = glom*bl1(i-1,1)
            ak1(i,2) = glom*bl1(i-1,2)-glom1*bl1(i,1)+dedr*dr(i-1,1)
            ak1(i,3) = glom*bl1(i-1,3)-glom1*bl1(i,2)+dedr*dr(i-1,2)
            ak1(i,4) =-glom1*bl1(i,3)
            ak2(i,1) = glom*bl2(i-1,1)
            ak2(i,2) = glom*bl2(i-1,2)-glom1*bl2(i,1) + deds
            ak2(i,3) =-glom1*bl2(i,2)
         endif
  60  continue
c
c          write the matrices to file tape11, if io is less
c       than 3, this write is not done.
c
      if( io .ge. 3 ) then
         write(11,7000)
         do 75 i=1,np1
          write(11,7001) i,ag1(i,1),ag1(i,2),ag1(i,3),ag2(i,1),ag2(i,2),
     $    ak1(i,1),ak1(i,2),ak1(i,3),ak1(i,4),ak2(i,1),ak2(i,2),ak2(i,3)
  75     continue
      endif
c
c          calculate the acoustic tavel time from surface to the
c       innermost zone. this is stored in transt.
c
      transt = zero
      do 77 i=2,np1
         transt = transt+(r(i)-r(i-1))/dsqrt(p(i-1)*v(i-1)*g1(i-1))
  77  continue
      rhom = totmas/(pi43*r(np1)**3)
      termq = dsqrt(rhom/1.41d0)
      tcon = pi43*g*rhom*float(n)**2
      x0np = r(np1)*dsqrt(dm2(np1))
      omsqp = zero
      omsqc = 5.d-4*float(ihmin)*tcon
      iv = 0
c
c          start of adiabatic loop
c
      nmodes = ihmax-ihmin + 1
      do 79 ih=ihmin,ihmax
      write(11,7900) ih
  78  continue
      omsqc = omsqc*two
      do 80 iterad=1,20
c
      if(io.ge.1) write(11,7800) iterad,iv,ih,omsq,omsqc,omsqp
      omsq = ahf*(omsqc+omsqp)
c
c          insure that initial omsqc gives iv .ge. iv-1
c
      if( iterad .eq. 1) omsq = omsqc
      omsqs = omsq
c
c          iteration on the outer boundary condition
c
      do 81 icount=1,30
c
         do 82 i=1,np1
            xo(i) = zero
            yo(i) = zero
  82     continue
         yo(n) =-ag1(n,3)*x0np
         call trisol(ag1,1,n,omsq,yo,xo,nmax,pramp,prpha)
	 if( icount .eq. 1 ) then
c
c          set up initial derivative
c
            qerra = ag1(np1,1)*xo(n)+(ag1(np1,2)-omsq)*x0np
            qerra = qerra/xo(2)
            omsql = omsq
            omsq  = omsq*(one+1.d-7)
         else
            qerr = ag1(np1,1)*xo(n)+(ag1(np1,2)-omsq)*x0np
            afacp = dabs((omsq-omsql)/omsql)
            qerr = qerr/xo(2)
	    if( io .ge. 1 ) write(11,1400) icount,omsq,omsql,afacp
	    if( afacp .lt. 1.d-11 ) goto 86
	       omsq2 = (qerra*omsq-qerr*omsql)/(qerra-qerr)
	       omsql = omsq
	       domsq = omsq2-omsq
	       domsq = dsign(dmin1(dabs(domsq),qrt*dabs(omsq)),domsq)
	       omsq  = omsq + domsq
	       qerra  = qerr
         endif
  81  continue
c
c          no convergence in the adiabatic eigenvalue.
c
      write(6,8100) omsq,iterad,iv,ih
      omsqp = zero
      omsqc = 5.d-6*tcon*float(ih+1)
      goto 79
c
c          converged to omsq value, check if positive and with
c      with proper number of nodes in the displacement eigen-
c      vector. if either is not true, try again.
c
  86  continue
      if( omsq .le. zero ) then
c
c          negative omsq value, stop working on this mode and
c       move on.
c
         write(6,8600) omsq,iv,ih
         omsqp = zero
         omsqc = 5.d-6*tcon*float(ih+1)
         goto 79
      endif
c
      do 90 i=1,n
         yo(i) = zero
  90  continue
      xo(np1) = x0np
      yo(n)  =-ag1(n,3)*x0np
      call trisol(ag1,1,n,omsq,yo,xo,nmax,pramp,prpha)
      iv = nodes(xo,np1,io)
      if( iv .eq. ih-1 ) goto 103
c
c      if initial omsgc gives iv .lt. ih-1, increase omsgc and try again
c
         if( omsqs .eq. omsqc .and. iv .lt. ih-1 ) goto 78
         if( iv .lt. ih-1 ) omsqp = omsqs
         if( iv .gt. ih-1 ) omsqc = omsqs
  80  continue
c
c          not converged to a period with right number of nodes,
c       try for next mode.
c
      omega = dsqrt(omsq)
      period = twopi/(omega*86400.d0)
      write(11,8000) iv,ih,omega,period
      omsqp = zero
      omsqc = 5.d-6*tcon*float(ih)
      goto 79
c
 103  continue
      omega = dsqrt(omsq)
      w(1) = xo(1)/(r(1)*dsqrt(dm2(1)))
      dro(1) = xo(1)*drdm0
      dp(1) = g10*dro(1)
      dt(1) = g3m10*dro(1)
      glamp(1) =-for*(g*rm(1)/r(1))*(xo(1)/r(1))**2 +
     $      pc*g10*dro(1)**2*cormas/rhoc
      stwait = glamp(1)
      rke = xo(1)*xo(1)
      dl(1) = bl1(1,2)*xo(1) + bl1(1,3)*xo(2)
      do 100 i=2,np1
         w(i) = xo(i)/(r(i)*dsqrt(dm2(i)))
         dro(i) = dr(i-1,1)*xo(i-1)+dr(i-1,2)*xo(i)
         dp(i) = g1(i-1)*dro(i)
         dt(i) = g3m1(i-1)*dro(i)
c
c          Set up the weight function calculation. GLAMP will have
c       the weight function per zone as described in the PASP article
c       above. See also Shapiro and Teukolsky? book on Black Holes,
c       White Dwarfs, and Neutron Stars.
c
         glamp(i) =-for*(g*rm(i)/r(i))*(xo(i)/r(i))**2 +
     $      p(i-1)*v(i-1)*g1(i-1)*dro(i)**2*dm1(i-1)
         stwait = stwait + glamp(i)
         rke = rke + xo(i)*xo(i)
         if( i .eq. np1 ) goto 100
         dl(i) = bl1(i,1)*xo(i-1)+bl1(i,2)*xo(i)+bl1(i,3)*xo(i+1)
 100  continue
      dl(np1) = bl1(np1,1)*xo(n)+bl1(np1,2)*x0np
      w(np1) = one
c
c          Given the radial eigenvectors, calculate the Epstein
c       weight functions. See Epstein Ap. J. 112 (6) 1950.
c
      epstwt(1) = (thre*g1(1)-for)*g*rm(1)*dm2(1)/r(1)*w(1)**2
      chwait = epstwt(1)
      do 104 i=2,n
        epstwt(i)=(thre*g1(i)-for)*g*rm(i)*dm2(i)/r(i)*w(i)**2+
     $     p(i)*v(i)*g1(i)*((w(i+1)-w(i))/dlog(r(i+1)/r(i)))**2*dm1(i)
         chwait = chwait + epstwt(i)
 104  continue
      chwait = chwait/rke
      ech = dabs((omsq-chwait)/omsq)
c
      period = twopi/omega
      call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax,itrans)
      radial%period   = period
      period = period/86400.d0
      qvalue = period*termq
      stwait = stwait/rke
      qch = dabs((stwait-omsq)/omsq)
c
      write(6,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch,
     $              chwait,ech
      write(1,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch,
     $              chwait,ech
      write(11,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch,
     $              chwait,ech
      iv = nodes(w,np1,1)
c
c          Set up for summary file
c
      radial%mass_flag = nint(10.d0*(totmas/1.991d33))
      radial%ell      = 0
      radial%k        = iv
      radial%n_p      = iv
      radial%n_g      = 0
      radial%omega_ad = omega
      radial%omsq_n   = omsq
      radial%acc_ad_w = qch
      radial%acc_ad_e = ech
c
c          normalize the weight function per zone to omsq.
c
      do 105 i=1,np1
         epstwt(i) = epstwt(i)/(omsq*rke)
         glamp(i) = glamp(i)/(omsq*rke)
 105  continue
c
      if(io.ge.2) write(11,1000) omega,period,(i,w(i),dl(i),dro(i),
     $ dt(i),dp(i),cvt(i),glamp(i),i=1,np1)
c
c          dump adiabatic eigenvectors to casplt.
c
      call pltdmp(w,     nmax,np1, 'dr/r')
      call pltdmp(dro,   nmax,np1, 'drho')
      call pltdmp(dl,    nmax,np1, 'dl/l')
      call pltdmp(glamp, nmax,np1, 'wait')
      call pltdmp(epstwt,nmax,np1, 'epst')
c
c          Find first guess to the imaginary part of omega by the use
c       of the quasi-adiabatic approximation. i.e.:
c
c               del(omega**2) = xo*ag2*(i*omega-ak2)**-1*ak1*xo.
c
c       Note that the minus signs in yo(i) are compensated for by the
c       routine 'TRCSOL' solving the resolvent as ak2 - i*omega.
c
      yo(1) =-ak1(1,2)*xo(1)-ak1(1,3)*xo(2)-ak1(1,4)*xo(3)
      y1(1) = czero
c      xo(n+2) = zero
      do 110 i=2,n-1
         yo(i) =-ak1(i,1)*xo(i-1)-ak1(i,2)*xo(i)-ak1(i,3)*xo(i+1)-
     $         ak1(i,4)*xo(i+2)
         y1(i) = czero
 110  continue
      yo(n) =-ak1(n,1)*xo(n-1)-ak1(n,2)*xo(n)-ak1(n,3)*x0np
      y1(n) = czero
c      yo(np1) =-ak1(np1,1)*xo(nm1)-ak1(np1,2)*xo(n)-ak1(np1,3)*x0np
      yo(np1) =-ak1(np1,1)*xo(n)-ak1(np1,2)*x0np
      y1(np1) = czero
      call trcsol(ak2,nmax,1,np1,ci*omega,yo,y1,dlth,dlmech)
      czi = czero
      do 120 i=1,n
         czi = czi + xo(i)*(ag2(i,1)*y1(i)+ag2(i,2)*y1(i+1))
 120  continue
      czi = czi + x0np*ag2(np1,1)*y1(np1)
      omsq1 = omsq + czi/rke
c
c          Limit the quasi-adiabatic growth rate to 0.1
c
      comega = cdsqrt(omsq1)
      growth =-forpi*dimag(comega)/dreal(comega)
      if( abs(growth) .gt. 0.1d0 ) then
         omega_im =-dsign(0.1d0, growth)*omega/forpi
         comega = dcmplx( omega, omega_im)
         growth =-forpi*dimag(comega)/dreal(comega)
         omsq1 = comega*comega
      endif
      write(1,1250) rke,comega,growth
      write(11,1250) rke,comega,growth
      radial%del_om_qa = comega
c
c          nonadiabatic iteration on mode iv starts here
c
      do 130 icount=1,20
      do 140 i=1,np1
         iy = 2*i-1
         ix = 2*i
         xs(ix) = czero
         xs(iy) = czero
c
c          AMAT is reset every iteration
c
         amat(iy,1) = ak1(i,1)
         amat(iy,2) = ak2(i,1)
         amat(iy,3) = ak1(i,2)
         amat(iy,4) = ak2(i,2)-ci*comega
         amat(iy,5) = ak1(i,3)
         amat(iy,6) = ak2(i,3)
         amat(iy,7) = ak1(i,4)
         amat(ix,1) = czero
         amat(ix,2) = ag1(i,1)
         amat(ix,3) = ag2(i,1)
         amat(ix,4) = ag1(i,2)-omsq1
         amat(ix,5) = ag2(i,2)
         amat(ix,6) = ag1(i,3)
         amat(ix,7) = czero
 140  continue
      write(11,1400) icount,omsq,comega,omsq1
      xs(2*np1-3) =-ak1(np1-1,4)*x0np
      xs(2*np1-2) =-ag1(n,3)*x0np
      xs(2*np1-1) =-ak1(np1,3)*x0np
      call cbmles(amat,2*nmax,1,2*np1-1,7,xs)
      rke = zero
      do 150 i=1,n
         y1(i) = xs(2*i-1)
         x0(i) = xs(2*i)
         rke = rke + conjg(x0(i))*x0(i)
 150  continue
      y1(np1) = xs(2*np1-1)
      x0(np1) = x0np
      rke = dsqrt(dabs(rke + x0np+x0np))
      if( icount .eq. 1 ) then
         oma = comega
         erra = ag1(np1,1)*x0(n)+(ag1(np1,2)-omsq1)*x0np +
     $             ag2(np1,1)*y1(np1)
c         erra = erra/dsign(rke,dreal(x0(2)))
         erra = erra/x0(2)
         omsq1 = omsq1*(one + 1.d-7)
         comega = cdsqrt(omsq1)
      else
         err = ag1(np1,1)*x0(n)+(ag1(np1,2)-omsq1)*x0np +
     $            ag2(np1,1)*y1(np1)
         afacp = cdabs((comega-oma)/oma)
c         err = err/dsign(rke,dreal(x0(2)))
         err = err/x0(2)
         if( io .ge. 1 ) write(11,1851) icount,afacp,comega,err
         if( afacp .lt. accur ) goto 190
         omb = (erra*comega-err*oma)/(erra-err)
         oma = comega
         comega = omb
         omsq1 = comega*comega
         erra = err
      endif
 130  continue
      write(6,1300) iv,afacp,accur
      write(1,1300) iv,afacp,accur
      write(11,1300) iv,afacp,accur
c
c          Iteration loop ends here, even if omega is not converged,
c       print the eigenvectors.
c
c          converged on mode iv
c
  190 continue
      rke = x0(1)*dconjg(x0(1))
      drho(1) = drdm0*x0(1)
      call amppha ( drho(1),rhoamp(1),rhopha(1))
      dpress(1) = drho(1)*g10+y1(1)*g3m10*rhoc/pc
      call amppha ( dpress(1),pramp(1),prpha(1))
      work(1) = dimag(dconjg(y1(1))*drho(1))*cormas*g3m10
      dtemp(1) = drho(1)*g3m10+y1(1)/(cv0*tc)
      call amppha ( dtemp(1),temamp(1),tempha(1))
      do 200 i=2,np1
         dror(i) = x0(i)/(r(i)*dsqrt(dm2(i)))
         call amppha ( dror(i),xamp(i),xpha(i))
         drho(i) = dr(i-1,1)*x0(i-1)+dr(i-1,2)*x0(i)
         call amppha ( drho(i),rhoamp(i),rhopha(i))
         dpress(i) = drho(i)*g1(i-1)+y1(i)*g3m1(i-1)/v(i-1)/p(i-1)
         call amppha ( dpress(i),pramp(i),prpha(i))
         rke = rke + x0(i)*conjg(x0(i))
         work(i) = dimag(dconjg(y1(i))*drho(i))*dm1(i-1)*g3m1(i-1)
         dtemp(i) = drho(i)*g3m1(i-1)+y1(i)/(cv(i-1)*t(i-1))
         call amppha ( dtemp(i),temamp(i),tempha(i))
c
c          DLTH and DLMECH are the thermal and mechanical
c       contributions to the luminosity variations.
c
         if( i .eq. np1 ) then
            dlmech(i) = bl1(i,1)*x0(i-1)+bl1(i,2)*x0(i)
         else
            dlmech(i) = bl1(i,1)*x0(i-1)+bl1(i,2)*x0(i)+bl1(i,3)*x0(i+1)
         endif
         dlth(i) = bl2(i,1)*y1(i-1) + bl2(i,2)*y1(i)
         dlum(i) = dlth(i) + dlmech(i)
         call amppha ( dlum(i),glamp(i),glpha(i))
         y1(i-1) = y1(i-1)/(cv(i-1)*t(i-1))
 200  continue
      dror(1) = x0(1)/(r(1)*dsqrt(dm2(1)))
      call amppha(dror(1),xamp(1),xpha(1))
      dlum(1) = bl1(1,2)*x0(1)+bl1(1,3)*x0(2)
     $        + bl2(1,1)*y1(1)+bl2(1,2)*y1(2)
      call amppha(dlum(1),glamp(1),glpha(1))
      call cptdmp(dror,nmax,np1,'dr/r')
      call cptdmp(drho,nmax,n,'drho')
      call cptdmp(dlum,nmax,np1,'dl/l')
      call cptdmp(y1,nmax,n,'dscv')
c
      omega = dreal(comega)
      rke = rke*omega*omega*ahf
      omimag = dimag(comega)
      growth =-forpi*dimag(comega)/omega
c
c          Normalized work integral and integrated work function.
c
      work(np1) = zero
      yo(np1) = zero
      do 205 ii=1,n
         i = n - ii + 1
         yo(i) = yo(i+1) + work(i)*pi/rke
         work(i) = work(i)*pi/rke/dabs(growth)
 205  continue
      qch = dabs((growth-yo(1))/growth)
c
      do 215 i=1,np1
         yo(i) = yo(i)/dabs(growth)
 215  continue
      call pltdmp(work,nmax,np1,'work')
      call pltdmp(yo,nmax,np1,'wint')
c
c          define the incarnations of the period.
c
      period = twopi/omega
      phase = transt*two/period
      pdays = period/86400.d0
      radial%comega   = comega
      radial%kappa    = growth
      radial%acc_na_k = qch
c
      s = yo(1)*dabs(growth)
      write(6,3000) comega,growth,s,qch,period,pdays,phase
      write(1,3000) comega,growth,s,qch,period,pdays,phase
      write(11,3000) comega,growth,s,qch,period,pdays,phase
c
      if(io.gt.1) write(11,3100) omega,period,growth,(i,xamp(i+1),
     $  xpha(i+1),glamp(i+1),glpha(i+1),rhoamp(i),rhopha(i),temamp(i),
     $  tempha(i),pramp(i),prpha(i),work(i),yo(i),i=1,n)
      afac = omega**2/(pi43*g*rhom)
      afacp = pdays*termq
      write(11,3200) afac,afacp
      if( io .gt. 1 ) then
         write(11,3300) (i,dror(i+1),dlth(i+1),dlmech(i+1),dlum(i+1),
     $                     drho(i),y1(i),i=1,n)
      endif
      call workmx(yo,np1,nmax)
      write(11,3600) rke
      call zie(work,n,nmax)
c      call summary( radial )
      omsqp = omsq
      omsqc = omsq
      call nldump( n, nmodes, iv, comega)
  79  continue
      return
c
 1250 format(25h quasi-adiabatic results:,/,5x,3h j=,1pe11.3,
     $   7h omega=,2e12.4,14h growth rate =,e12.4)
 1400 format(1x,i4,1p,7e11.4)
 1851 format(1x,10hlna iter =,i3,7h afacp=,1pe10.3,8h comega=,
     $   2e10.3,5h err=,2e10.3)
 3000 format(/,1x,25hlinear non-adiabatic mode,/,8h omega =,1pe16.9,
     $   1x,9h(sec**-1),3x,15h omega (imag) =,e16.9,/,
     $   9x,24h predicted growth rate =,e14.6,/,
     $   1x,32hgrowth rate from work integral =,e14.6,7h error=,e11.3,
     $   //,3x,9hperiods =,e12.4,7h (secs),e13.4,7h (days),/,
     $   1x,23h acoustic travel time =,e13.6,8h periods)
 3100 format(1h1,27h linear, non-adiabatic mode,//,9h omega = ,
     $ 1pe16.8,9h sec(-1),,2x,9hperiod = ,e15.7,6h days,,3x,
     $ 34hfraction energy gain per period = ,e14.6,//,
     $ 3x,1hi,10x,2hdx,20x,2hdl,18x,3hdro,18x,2hdt,19x,2hdp,17x,4hwork,
     $ /,4x,5(5x,3hamp,6x,3hpha,4x),3x,4hzone,4x,8hintegral,/,
     $ (1x,i3,6(1x,2e10.3)) )
 3200 format(/,25h dimensionless frequency=,1pe12.4,
     $ /,17h q-value in days=,e12.5)
 3300 format(1h1,28h cartesian form of dx and dl,/,
     $ 3x,1hi,9x,2hdx,17x,6hdlther,15x,6hdlmech,17x,4hdl/l,17x,
     $ 4hdrho,16x,5hds/cv,/,4x,6(3x,4hreal,6x,5himag.,3x),/,
     $  (1x,i3,6(1x,1p,2e10.3)) )
 3600 format(/,15h kinetic energy,1pe16.8,7h (ergs))
 1300 format(/,43h can only converge on non-adiabatic period ,i2,3h to,
     $ 1pe10.3,7h accur=,e10.3,20h after 20 iterations)
c
 7000 format(1h1,/1x,43hpulsation matrices: ag1(3), ag2(2), ak1(4),,
     $              12h and ak2(3).)
 7001 format(1x,i4,1x,1p,3e10.3,1x,2e10.3,1x,4e10.3,1x,3e10.3)
 7800 format(8h iterad=,i3,4h iv=,i3,4h ih=,i3,6h omsq=,1pe12.4,
     $   7h omsqc=,e12.4,7h omsqp=,e12.4)
 8000 format(40h cant find no. of nodes.eq.ih-1...nodes=,i3,4h ih=,i3,
     $ 7h omega=,1pe13.6,10h period = ,e13.6)
 8600 format(14h adiabatic f2=,1pe17.10,18h is .lt. 0 for iv=,
     $  i3,4h ih=,i3)
 8100 format(14h adiabatic f2=,1pe15.6,20h not converged after,i4,
     $   16h iterations, ih=,i4,17h number of nodes=,i4)
 1000 format(1h1,19h adiabatic solution,/,9h  omega =,1pe17.9,5x,
     $ 8hperiod =,0pf12.6,//,3x,1hi,7x,2hdr,12x,2hdl,11x,4hdrho,11x,
     $  2hdt,12x,2hdp,11x,4hcvtm,10x,6hweight/(1x,i4,1p,7e14.6))
 1001 format(/,22h linear adiabatic mode,4h iv=,i3,4h ih=,i3,/,
     $ 1x,8h omega =,1pe16.8,9h sec(-1),,1x,8hperiod =,e12.5,
     $ 6h days,,1x,9hq-value =,e10.3,/,
     $ 1x,38h eigenfrequency from matrix solution =,e12.4,/,
     $ 1x,38h eigenfrequency from weight function =,e12.4,
     $ 7h error=,e12.4,/,
     $ 1X,42h eigenvalue from Epstein weight function =,e12.5,
     $ 7H error=,e10.3)
 7900 format(1h1,//,20x,28h***** begin search for mode ,i2,6h *****)
      end
      subroutine nldump( npts, nmodes, imode, comega )
      implicit real*8(a-h,o-z)
      integer*2 i
      logical first
c
c          Dump the information for a nonlinear startup file to the
c      file NLDUMP.
c
      parameter ( nmax=512 )
      parameter( nmax2=2*nmax)
      complex*16 x0,y0,y1,xs,dror,dlum,comega,drho,dpress,
     $  dtemp,dlth,dlmech
c
      common xs(nmax2),y0(nmax),y1(nmax),x0(nmax),dror(nmax),
     $ dlth(nmax),dlmech(nmax),dlum(nmax),drho(nmax),dtemp(nmax),
     $   dpress(nmax), xamp(nmax),
     $ xpha(nmax),rhoamp(nmax),rhopha(nmax),pramp(nmax),prpha(nmax),
     $ temamp(nmax),tempha(nmax),glamp(nmax),glpha(nmax),
     $ w(nmax),cvt(nmax)
      data first/.true./
c
      if( first ) then
         open( unit=40, file='nldump', status='unknown')
         write(40, 1000) nmodes
         first = .false.
      endif
c
      write(40,1001) imode
      write(40,1010) comega
      write(40,1010) (dror(i),i=1,npts+1)
      write(40,1010) (dlum(i),i=1,npts+1)
      write(40,1010) (drho(i),i=1,npts)
      write(40,1010) (dtemp(i),i=1,npts)
      write(40,1010) (dpress(i),i=1,npts)
      return
c
 1000 format(i5,10x,'Number of modes dumped to file')
 1001 format(i5,10x,'Node count of current mode')
 1010 format(1p4e20.13)
      end
