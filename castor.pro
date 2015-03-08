function castor, io,ihmin,ihmax, npts, r, t, v, cv, dkdr, $
        dkdt, dm1, kappa, dm2, rm, N2, p, g1, g3m1, rho, rzone, xo
;test = castor(iout, 0, 3, npts-1, data_st.radius, data_st.t, v, data_st.cv, data_st.dkdr, $
;        data_st.dkdt, dm1, data_st.kappa, dm2, data_st.Mass, data_st.N2, data_st.p, $
;        data_st.gamma1, g3m1, data_st.rho, rzone)
;      implicit real*8(a-h,o-z)
;      integer*2 i, ii, j, icount
;c      logical first /.true./
;c
;c     *****************************************************************
;c      Pulsation Equation Solver
;c      Basic references -- Castor, Ap. J. 166 (109) 1971.
;c                          Pesnell, 1983, Thesis, Univ. of Florida.
;c                          Pesnell, PASP, 99 (975), 1987.
;c     *****************************************************************
;c
;c          This version of the linear nonadiabatic radial stability 
;c       analysis is suitable for envelope models generated with the 
;c       MODEL building code. The inner boundary is valid for a complete 
;c       model and models of low central condensation. If the values of 
;c       the central ball are not present, or inconsistent, the results 
;c       of the stability analysis are questionable.
;c
;c
;c       06/01/94: Moved the conversion of DKDT to DKDS to main routine.
;c
;      complex*16 x0,y0,y1,z,xs,dror,dlum,omsq1,comega,drho,dpress,
;     $  dtemp,amat,oma,erra,err,omb,czero,ci,czi,dlth,dlmech
;      parameter ( czero=(0.d0,0.d0), ci=(0.d0,1.d0), accur=1.d-10 )
;      structure /puls_sum/
;         integer*2 mass_flag
;         integer*4 ell, k, n_p, n_g
;         real*4 omega_ad, omsq_n, period, acc_ad_w, acc_ad_e
;         complex*8 del_om_qa
;         complex*8 comega
;         real*4 kappa, acc_na_k
;      end structure
;      record /puls_sum/ radial
;c
;c        arrays are rearranged from the original lna code
;c
;      parameter ( nmax=512 )
;      parameter( nmax2=2*nmax)
;      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
;     $	     dkdt(nmax),dm1(nmax),gkap(nmax),dm2(nmax),rm(nmax),
;     $		     drz(nmax),drint(nmax),gor(nmax)
;      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
;      common/scrtch/ ag1(nmax,3),ag2(nmax,2),ak1(nmax,4),ak2(nmax,3),
;     $	     bl2(nmax,2),bl1(nmax,3),dr(nmax,2)
;      common xs(nmax2),y0(nmax),y1(nmax),x0(nmax),dror(nmax),
;     $	 dlth(nmax),dlmech(nmax),dlum(nmax),drho(nmax),dtemp(nmax),
;     $   dpress(nmax), xamp(nmax),
;     $	 xpha(nmax),rhoamp(nmax),rhopha(nmax),pramp(nmax),prpha(nmax),
;     $	 temamp(nmax),tempha(nmax),glamp(nmax),glpha(nmax),
;     $	 w(nmax),cvt(nmax)
;      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
;      common/blk37/  yo(nmax),xo(nmax),work(nmax),dro(nmax),dp(nmax),
;     $	     dt(nmax),dl(nmax)
;      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
;      common/observ/ teff,rlumgv,totmas,rphoto,corlum
;      common/const/  zero,one,two,thre,for,ten,ahf,qrt
;      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
;     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
;      dimension amat(nmax2,7), epstwt(nmax)
;      dimension zzero(nmax,36),szero(nmax,19)
;      equivalence (zzero(1,1),xs(1)),(szero(1,1),ag1(1,1))
;      equivalence (amat(1,1),y0(1))
;c
;c   Zero common blocks scrtch and (blank).
;c
;      do 5 j=1,36
;         do 5 i=1,nmax
;            zzero(i,j) = zero
;   5  continue
;      do 6 j=1,19
;         do 6 i=1,nmax
;            szero(i,j) = zero
;   6  continue
puls_sum = { mass_flag: 0, ell: 0, k: 0, n_p: 0, n_g: 0, $
             omega_ad: 0., omsq_n: 0., period: 0., acc_ad_w: 0., acc_ad_e: 0., $
             dror: dblarr(npts+1), drho: dblarr(npts+1) }
;         complex*8 del_om_qa
;         complex*8 comega
;         real*4 kappa, acc_na_k
;      end structure
radial = replicate( puls_sum, ihmax-ihmin+1)
;c
   n_win = 2
      n = npts
      np1 = n+1
      nm1 = n-1
;c
;      rlumx = corlum
G = 6.67d-8
      drdm = [ [ (4.d0*!pi)*r[0:n-1]^2/sqrt(dm2[0:n-1])/v/dm1], $
               [-(4.d0*!dpi)*r[1:n]^2/sqrt(dm2[1:n])/v/dm1] ]

      ag1 = dblarr(np1+1,3)
      dr = dblarr(np1+1,2)
      afac = (4.d0*!pi)*r[0]^2/sqrt(dm2[0])
      for i=0,npts do begin
         afacp = afac
         if( i lt n ) then afac = (4.d0*!pi)*r(i+1)^2/sqrt(dm2(i+1))
         if( i lt n ) then begin
            dr[i,0] = afacp/v[i]/dm1[i]
            dr[i,1] =-afac/v[i]/dm1[i]
         endif
         if( i eq n ) then begin
            g1p = g1(n)*p(n)
;            g3m1ov = g3m1[n-1]/v[n-1]
            ag1(n,0) =-afacp*g1p*dr(n,0)
            ag1(n,1) =-4.d0*g*rm[n]/r[n]^3-afacp*g1p*dr[n,1]
            ag1(n,2) = 0.d0
            print, g1p, dr[n,0]
;            ag1[n,0] = ag1[n-1,2]
;            ag1[n,1] = ag1[n-1,1]
;c
;c          The inner boundary condition is zero gradient, not zero
;c       motion.
;c
         endif else if( i eq 0 ) then begin
            drdm0 =-3.d0/(r(1)*sqrt(dm2(1)))
            qafacp = (4.d0*!pi)*r(1)^2/sqrt(dm2(1))
            qafac = (4.d0*!pi)*r(2)^2/sqrt(dm2(2))
            ag1(1,0) = 0.d0
            ag1(1,1) =-4.d0*g*rm[1]/r[1]^3 + qafacp*(g1[1]*p[1]*dr[1,0]) ; - drdm0*pc*g10 )
            ag1(1,2) = qafacp*g1[1]*p[1]*dr[1,1]
         endif else begin
;c
;c          set up the mechanical matrices in a regular zone.
;c
            g1p = g1[i]*p[i]
;c         if( i eq n ) then
;c            g1p = g1p-1.00855d-14*t[i]^4*g3m1[i]
;c         endif
            ag1[i,0] =-afacp*g1[i-1]*p[i-1]*dr(i-1,0)
            ag1[i,1] =-4.d0*g*rm[i]/r[i]^3 + afacp*(g1p*dr[i,0]- g1[i-1]*p[i-1]*dr(i-1,1) )
            ag1[i,2] = afacp*g1p*dr[i,1]
         endelse
       endfor
       for i=n-5, n do print, format='(1x,i4,1x,3e10.3)', i,ag1[i,0],ag1[i,1],ag1[i,2]
;c
;c          write the matrices to file tape11, if io is less
;c       than 3, this write is not done.
;c
      if( io ge 3 ) then begin
         print, 'pulsation matrices: ag1(3)'
;         for i=1,n do print, format='(1x,i4,1x,3e10.3)', i,ag1[i,0],ag1[i,1],ag1[i,2]
      endif
      start_plot, n_win, 'd'
      plot, ag1[*,0], line=1, yrange=[-10,10]
      oplot, ag1[*,1]
      oplot, ag1[*,2], line=5
      end_plot
;c
;c          calculate the acoustic tavel time from surface to the
;c       innermost zone. this is stored in transt.
;c
;      transt = 0.d0
;      for i=0,n-1 do transt = transt+(r[i+1]-r[i])/sqrt(p[i]*v[i]*g1[i])
      transt = total((r[1:n]-r[0:n-1])/sqrt(p*v*g1), /double)
      transt2 = int_tabulated( r, 1./sqrt(p*v*g1), /double)
      rhom = max(rm)/r(n)^3/(4.d0*!pi/3.d0)
      termq = sqrt(rhom/1.41d0)
      tcon = (4.d0*!pi/3.d0)*g*rhom*float(n)^2
      print, transt, transt2, tcon
      x0np = r(n)*sqrt(dm2(n-1))
      omsqp = 0.d0
      omsqc = 5.d-4*float(ihmin)*tcon
      iv = 0
;c
;c          start of adiabatic loop
;c
   nmodes = ihmax-ihmin + 1
   for ih=ihmin,ihmax do begin
      print, ' ***** begin search for mode ', ih, ' *****'
L78: ;  continue
      omsqc = omsqc*2.d0
      for iterad=1,20 do begin
;c
;         if(io ge 1) write(11,7800) iterad,iv,ih,omsq,omsqc,omsqp
         omsq = 0.5d0*(omsqc+omsqp)
;c
;c          insure that initial omsqc gives iv .ge. iv-1
;c
         if( iterad eq 1) then omsq = omsqc
         omsqs = omsq
;c
;c          iteration on the outer boundary condition
;c
         for icount=1,20 do begin
;c
            yo = dblarr(n)
            yo(n-2) =-ag1(n-1,2)*x0np
            B_arr = ag1[1:n,1]-replicate(omsq, n)
            xo = trisol(ag1[1:n,0], B_arr, ag1[1:n,2], yo)
            if( icount eq 1 ) then begin
;c
;c          set up initial derivative
;c
               qerra = ag1(n,0)*xo(n-1)+(ag1(n,1)-omsq)*x0np
               qerra = qerra/xo[1]
               omsql = omsq
               omsq  = omsq*(1.+1.d-7)
               if( io ge 1 ) then print, icount,omsq,omsql,qerra
            endif else begin
               qerr = ag1(n,0)*xo(n-1)+(ag1(n,1)-omsq)*x0np
               afacp = abs((omsq-omsql)/omsql)
               qerr = qerr/xo[1]
               if( io ge 1 ) then print, icount,omsq,omsql,afacp
               if( afacp lt 1.d-11 ) then goto, L86
      	       omsq2 = (qerra*omsq-qerr*omsql)/(qerra-qerr)
	             omsql = omsq
      	       domsq = omsq2-omsq
	             domsq = sign(min([abs(domsq),0.25*abs(omsq)]),domsq)
      	       omsq  = omsq + domsq
	             qerra  = qerr
            endelse
         endfor
;c
;c          no convergence in the adiabatic eigenvalue.
;c
      print, omsq,iterad,iv,ih, format='("adiabatic f2=",e15.6," not converged after",i4," iterations, ih=",i4," number of nodes=",i4)'
      omsqp = 0.d0
      omsqc = 5.d-6*tcon*float(ih+1)
      goto, L79
;c
;c          converged to omsq value, check if positive and with
;c      with proper number of nodes in the displacement eigen-
;c      vector. if either is not true, try again.
;c
L86:
      if( omsq le 0.d0 ) then begin
;c
;c          negative omsq value, stop working on this mode and
;c       move on.
;c
         print, omsq,iv,ih, format='("adiabatic f2=",1pe17.10," is .lt. 0 for iv=",i3," ih=",i3)
         omsqp = 0.d0
         omsqc = 5.d-6*tcon*float(ih+1)
         goto, L79
      endif
;c
      for i=0,n-1 do yo[i] = 0.d0
      xo(n-1) = x0np
      yo(n-1)  =-ag1(n-1,2)*x0np
      xo = trisol(ag1[1:n,0], ag1[1:n,1]-omsq, ag1[1:n,2], yo)
      iv = nodes(xo,io)
      if( iv eq ih-1 ) then goto, L103
;c
;c      if initial omsgc gives iv .lt. ih-1, increase omsgc and try again
;c
         if( omsqs eq omsqc and iv lt ih-1 ) then goto, L78
         if( iv lt ih-1 ) then omsqp = omsqs
         if( iv gt ih-1 ) then omsqc = omsqs
   endfor ;  80  continue
;c
;c          not converged to a period with right number of nodes,
;c       try for next mode.
;c
      omega = sqrt(omsq)
      period = 2.d0*!pi/(omega*86400.d0)
      print, iv,ih,omega,period, format='("cant find no. of nodes.eq.ih-1...nodes=",i3," ih=",i3," omega=",e13.6," period = ",e13.6)'
      
      omsqp = 0.d0
      omsqc = 5.d-6*tcon*float(ih)
      goto, L79
;c
L103: ;  continue
      omega = sqrt(omsq)
      w = [xo/(r*sqrt(dm2)), 1.d0]
;      start_plot, n_win, 'd'
;      plot, r/max(r), w
;      oplot, !x.crange, [0,0], line=1
;      end_plot
      drho = dblarr(n+1)
      drho[0] = xo(1)*drdm0
      for i=1,n-1 do begin
         drho[i] = dr[i-1,0]*xo[i-1]+dr[i-1,1]*xo[i]
      endfor
;      dro(1) = xo(1)*drdm0
;      dp(1) = g10*dro(1)
;      dt(1) = g3m10*dro(1)
;      glamp(1) =-4.d0*(g*rm(1)/r(1))*(xo(1)/r(1))^2 + pc*g10*dro(1)^2*cormas/rhoc
;      stwait = glamp(1)
;      rke = xo(1)*xo(1)
;;      dl(1) = bl1(1,2)*xo(1) + bl1(1,3)*xo(2)
;      for i=2,np1 do begin
;         w[i] = xo[i]/(r[i]*sqrt(dm2[i]))
;         dro[i] = dr(i-1,0)*xo[i-1]+dr(i-1,1)*xo[i]
;         dp[i] = g1[i-1]*dro[i]
;         dt[i] = g3m1[i-1]*dro[i]
;;c
;;c          Set up the weight function calculation. GLAMP will have
;;c       the weight function per zone as described in the PASP article
;;c       above. See also Shapiro and Teukolsky? book on Black Holes,
;;c       White Dwarfs, and Neutron Stars.
;;c
;         glamp[i] =-4.d0*(g*rm[i]/r[i])*(xo[i]/r[i])^2 + p[i-1]*v[i-1]*g1[i-1]*dro[i]^2*dm1[i-1]
;         stwait = stwait + glamp[i]
;         rke = rke + xo[i]*xo[i]
;         if( i lt np1 ) then dl[i] = bl1[i,1]*xo[i-1]+bl1[i,2]*xo[i]+bl1[i,3]*xo(i+1)
;      endfor
 ;     dl(np1) = bl1(np1,1)*xo(n)+bl1(np1,2)*x0np
;      w(np1) = one
;c
;c	    Given the radial eigenvectors, calculate the Epstein
;c	weight functions. See Epstein Ap. J. 112 (6) 1950.
;c
;      epstwt(1) = (3.d0*g1(1)-4.d0)*g*rm(1)*dm2(1)/r(1)*w(1)^2
;      chwait = epstwt(1)
;      for i=2,n do begin
;        epstwt[i]=(3.d0*g1[i]-4.d0)*g*rm[i]*dm2[i]/r[i]*w[i]^2+ $
;          p[i]*v[i]*g1[i]*((w(i+1)-w[i])/dlog(r(i+1)/r[i]))^2*dm1[i]
;         chwait = chwait + epstwt[i]
;      endfor
;      chwait = chwait/rke
;      ech = dabs((omsq-chwait)/omsq)
;c
      period = 2.d0*!pi/omega
;      call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax,itrans)
      radial[ih-1].period   = period
      period = period ;/86400.d0
      qvalue = period*termq
;      stwait = stwait/rke
;      qch = dabs((stwait-omsq)/omsq)
;c
   print, iv,ih,omega,period,qvalue,omsq ;,stwait,qch, chwait,ech
;      write(6,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch, chwait,ech
;      write(1,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch, chwait,ech
;      write(11,1001) iv,ih,omega,period,qvalue,omsq,stwait,qch, chwait,ech
      iv = nodes(w,1)
;c
;c          Set up for summary file
;c
   radial[ih-1].mass_flag = nint(10.d0*(max(rm)/1.991d33))
   radial[ih-1].ell      = 0
   radial[ih-1].k        = iv
   radial[ih-1].n_p      = iv
   radial[ih-1].n_g      = 0
   radial[ih-1].omega_ad = omega
   radial[ih-1].omsq_n   = omsq
   radial[ih-1].dror = w
   radial[ih-1].drho = drho
;      radial.acc_ad_w = qch
;      radial.acc_ad_e = ech
;c
;c          normalize the weight function per zone to omsq.
;c
;      for i=1,np1 do begin
;         epstwt[i] = epstwt[i]/(omsq*rke)
;         glamp[i] = glamp[i]/(omsq*rke)
;      endfor
;c
;      if(io.ge.2) write(11,1000) omega,period,(i,w[i],dl[i],dro[i],dt[i],dp[i],cvt[i],glamp[i],i=1,np1)
;c
;c          dump adiabatic eigenvectors to casplt.
;c
;      call pltdmp(w,     nmax,np1, 'dr/r')
;      call pltdmp(dro,   nmax,np1, 'drho')
;      call pltdmp(dl,    nmax,np1, 'dl/l')
;      call pltdmp(glamp, nmax,np1, 'wait')
;      call pltdmp(epstwt,nmax,np1, 'epst')
;c
;c          Find first guess to the imaginary part of omega by the use
;c       of the quasi-adiabatic approximation. i.e.:
;c
;c               del(omega^2) = xo*ag2*(i*omega-ak2)^-1*ak1*xo.
;c
;c       Note that the minus signs in yo[i] are compensated for by the
;c       routine 'TRCSOL' solving the resolvent as ak2 - i*omega.
;c
;      yo(1) =-ak1(1,2)*xo(1)-ak1(1,3)*xo(2)-ak1(1,4)*xo(3)
;      y1(1) = 0.d0
;;c      xo(n+2) = 0.d0
;      do 110 i=2,n-1
;         yo[i] =-ak1[i,1]*xo[i-1]-ak1[i,2]*xo[i]-ak1[i,3]*xo(i+1)-
;     $         ak1(i,4)*xo(i+2)
;         y1[i] = 0.d0
; 110  continue
;      yo(n) =-ak1(n,1)*xo(n-1)-ak1(n,2)*xo(n)-ak1(n,3)*x0np
;      y1(n) = 0.d0
;;c      yo(np1) =-ak1(np1,1)*xo(nm1)-ak1(np1,2)*xo(n)-ak1(np1,3)*x0np
;      yo(np1) =-ak1(np1,1)*xo(n)-ak1(np1,2)*x0np
;      y1(np1) = 0.d0
;      call trcsol(ak2,nmax,1,np1,ci*omega,yo,y1,dlth,dlmech)
;;
;      czi = 0.d0
;      for i=1,n do czi = czi + xo[i]*(ag2[i,1]*y1[i]+ag2[i,2]*y1(i+1))
;      czi = czi + x0np*ag2(np1,1)*y1(np1)
;;
;      omsq1 = omsq + czi/rke
;c
;c          Limit the quasi-adiabatic growth rate to 0.1
;c
;      comega = sqrt(omsq1)
;      growth =-4.d0*!pi*dimag(comega)/dreal(comega)
;      if( abs(growth) .gt. 0.1d0 ) then
;         omega_im =-dsign(0.1d0, growth)*omega/(4.d0*!pi)
;         comega = dcmplx( omega, omega_im)
;         growth =-(4.d0*!pi)*dimag(comega)/dreal(comega)
;         omsq1 = comega*comega
;      endif
;      write(1,1250) rke,comega,growth
;      write(11,1250) rke,comega,growth
;      radial.del_om_qa = comega
L79: ;  continue
   endfor
return, radial
;c
; 1250 format(25h quasi-adiabatic results:,/,5x,3h j=,1pe11.3,
;     $   7h omega=,2e12.4,14h growth rate =,e12.4)
; 1400 format(1x,i4,1p,7e11.4)
; 1851 format(1x,10hlna iter =,i3,7h afacp=,1pe10.3,8h comega=,
;     $   2e10.3,5h err=,2e10.3)
; 3000 format(/,1x,25hlinear non-adiabatic mode,/,8h omega =,1pe16.9,
;     $   1x,9h(sec**-1),3x,15h omega (imag) =,e16.9,/,
;     $   9x,24h predicted growth rate =,e14.6,/,
;     $   1x,32hgrowth rate from work integral =,e14.6,7h error=,e11.3,
;     $   //,3x,9hperiods =,e12.4,7h (secs),e13.4,7h (days),/,
;     $   1x,23h acoustic travel time =,e13.6,8h periods)
; 3100 format(1h1,27h linear, non-adiabatic mode,//,9h omega = ,
;     $ 1pe16.8,9h sec(-1),,2x,9hperiod = ,e15.7,6h days,,3x,
;     $ 34hfraction energy gain per period = ,e14.6,//,
;     $ 3x,1hi,10x,2hdx,20x,2hdl,18x,3hdro,18x,2hdt,19x,2hdp,17x,4hwork,
;     $ /,4x,5(5x,3hamp,6x,3hpha,4x),3x,4hzone,4x,8hintegral,/,
;     $ (1x,i3,6(1x,2e10.3)) )
; 3200 format(/,25h dimensionless frequency=,1pe12.4,
;     $ /,17h q-value in days=,e12.5)
; 3300 format(1h1,28h cartesian form of dx and dl,/,
;     $ 3x,1hi,9x,2hdx,17x,6hdlther,15x,6hdlmech,17x,4hdl/l,17x,
;     $ 4hdrho,16x,5hds/cv,/,4x,6(3x,4hreal,6x,5himag.,3x),/,
;     $  (1x,i3,6(1x,1p,2e10.3)) )
; 3600 format(/,15h kinetic energy,1pe16.8,7h (ergs))
; 1300 format(/,43h can only converge on non-adiabatic period ,i2,3h to,
;     $ 1pe10.3,7h accur=,e10.3,20h after 20 iterations)
;c
; 7000 format(1h1,/1x,43hpulsation matrices: ag1(3), ag2(2), ak1(4),,
;     $              12h and ak2(3).)
; 7001 format(1x,i4,1x,1p,3e10.3,1x,2e10.3,1x,4e10.3,1x,3e10.3)
; 7800 format(8h iterad=,i3,4h iv=,i3,4h ih=,i3,6h omsq=,1pe12.4,
;     $   7h omsqc=,e12.4,7h omsqp=,e12.4)
; 8000 format(40h cant find no. of nodes.eq.ih-1...nodes=,i3,4h ih=,i3,
;     $ 7h omega=,1pe13.6,10h period = ,e13.6)
; 8600 format(14h adiabatic f2=,1pe17.10,18h is .lt. 0 for iv=,
;     $  i3,4h ih=,i3)
; 8100 format(14h adiabatic f2=,1pe15.6,20h not converged after,i4,
;     $   16h iterations, ih=,i4,17h number of nodes=,i4)
; 1000 format(1h1,19h adiabatic solution,/,9h  omega =,1pe17.9,5x,
;     $ 8hperiod =,0pf12.6,//,3x,1hi,7x,2hdr,12x,2hdl,11x,4hdrho,11x,
;     $  2hdt,12x,2hdp,11x,4hcvtm,10x,6hweight/(1x,i4,1p,7e14.6))
; 1001 format(/,22h linear adiabatic mode,4h iv=,i3,4h ih=,i3,/,
;     $ 1x,8h omega =,1pe16.8,9h sec(-1),,1x,8hperiod =,e12.5,
;     $ 6h days,,1x,9hq-value =,e10.3,/,
;     $ 1x,38h eigenfrequency from matrix solution =,e12.4,/,
;     $ 1x,38h eigenfrequency from weight function =,e12.4,
;     $ 7h error=,e12.4,/,
;     $ 1X,42h eigenvalue from Epstein weight function =,e12.5,
;     $ 7H error=,e10.3)
; 7900 format(1h1,//,20x,28h***** begin search for mode ,i2,6h *****)
end