c*byp3
c*rcft i=%me,x=%mex,lib=(discft,cgscft,intcft)
c      program ladnon(lnrad,tape10=lnrad,nonout,tape6=nonout,plotin,tape8
c     1=plotin,s1,s2,s3,tape11=s1,tape12=s2,tape13=s3,check,tape14=check,
c     2contents,tape15=contents)
      program ladnon
c
c          the lagrangian non-radial nonadiabatic stability analysis
c       for stellar pulsations. the notation is lifted from the model
c       code and the radial stability analysis. written nov. 1983
c       with nuclear burning and (frozen) convection added dec. 1983.
c
c                                                      12/16/83 wd pesnell
c
      common/phypar/ rp(2001),tp(2001),vp(2001),cv(2001),dkdr(2001),
     $   dkdt(2001),dm1(2001),akap(2001),dm2(2001),rm(2001),emdldm(2001)
      common/core/ncca,nccna,termax,err
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,ir,iplot,irp1,isum
      common/blk4/   p(2001),g1(2001),g3m1(2001),brntva(2001)
      common/blk8/ grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/coolum/ noburn,onemq0,onemq1,rlums,rlumc
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/lumins/ frft(2001),sorce(2001),dtsorc(2001),dvsorc(2001)
      common /sumtab/ nrq(100),npn(100),ngn(100),pa(100),pna(100)
     1 ,pgr(100),wgr(100),tke(100),iroot,rl,ititl(3)
      common /rbk/ rc(2001),zone(2001),rn(2001),ltop(3),rcn(2001),
     2 dnsity(2001)
      dimension sound(2001),capa(2001)
      character *100,lxname,lyname,lhead,index
      call filerep
      call libdisp
      call gplot(1hu,6hkidman,6)
      read(10,1001) (ititl(i),i=1,3)
      do 3 i=1,3
3     ltop(i)=ititl(i)
c
c          read in the input variables from file 10.
c
      read(10,*) ir,rlumgv,totmas,teff,rphoto,corlum
      read(10,*) irad,noburn,onemq0,onemq1,rlums,rlumc,fsuba
      irp1 = ir + 1
      read(10,1000) (rp(i),i=1,irp1)
       if(rp(1).eq.0.) rp(1)=1.
      read(10,1000) (tp(i),i=1,ir)
      read(10,1000) (vp(i),i=1,ir)
      read(10,1000) (cv(i),i=1,ir)
      read(10,1000) (dkdr(i),i=1,ir)
      read(10,1000) (dkdt(i),i=1,ir)
      read(10,1000) (dm1(i),i=1,ir)
      read(10,1000) (akap(i),i=1,ir)
      read(10,1000) (dm2(i),i=1,irp1)
      read(10,1000) (rm(i),i=1,irp1)
      read(10,1000) (p(i),i=1,ir)
      read(10,1000) (g1(i),i=1,ir)
      read(10,1000) (g3m1(i),i=1,ir)
      read(10,1000) (frft(i),i=1,irp1)
      read(10,1000) (sorce(i),i=1,ir)
      read(10,1000) (dtsorc(i),i=1,ir)
      read(10,1000) (dvsorc(i),i=1,ir)
      read(10,1000) (brntva(i),i=1,ir),p0
      read(10,*) iout,rl,nomega,omlow,omhigh,iplot,isum,it,ip,ncc,feps,j
      if(j.ne.0) read(10,1000) (emdldm(i),i=1,ir)
c       generate plots of the input data
      do 2 i=1,irp1
      zone(i)=float(i)
      rc(i)=.5*(rp(i)+rp(i+1))
      rn(i)=rp(i)/rp(irp1)
      rcn(i)=rc(i)/rp(irp1)
      if(i.ne.irp1) dnsity(i)=1./vp(i)
 2    continue
      if(iplot.eq.0) go to 19
      lxname="z[one]$"
      lyname="z[one] c[entered] r[adii (cm)$"
      index="zone centered radii vs zone number"
      call russ(zone,rc,ir,0,0,lxname,lyname,ititl,index)
      index="log of zone centered radii vs zone number"
      call russ(zone,rc,ir,0,1,lxname,lyname,ititl,index)
      lyname="t[emperature] (k[elvin)$"
      index="log(temperature) vs zone"
      call russ(zone,tp,ir,0,1,lxname,lyname,ititl,index)
      lyname="s[pecific] v[olume (cc/gm)$"
      index="log(specific volume) vs zone"
      call russ(zone,vp,ir,0,1,lxname,lyname,ititl,index)
      lxname="n[ormalized] r[adius$"
      lyname="d[ensity (gm/cc)$"
      index="log(density) vs normalized radius"
      call russ(rcn,dnsity,ir,0,1,lxname,lyname,ititl,index)
      lxname="z[one]$"
      lyname="s[pecific] h[eat (ergs/gm/kelvin)$"
      index="cv vs zone ... specific heat vs zone"
      call russ(zone,cv,ir,0,0,lxname,lyname,ititl,index)
      lyname="d[ln(opacity)/]d[ln(density)$"
      index="dkdr vs zone ... Dln(opacity)/Dln(density)"
      call russ(zone,dkdr,ir,0,0,lxname,lyname,ititl,index)
      lyname="d[ln(opacity)/]d[ln(temperature)$"
      index="dkdt vs zone ... Dln(opacity)/Dln(temperature)"
      call russ(zone,dkdt,ir,0,0,lxname,lyname,ititl,index)
      lyname="z[one] c[entered] s[hell] m[asses (gm)$"
      index="log(dm1) vs zone ... zone mass vs zone number"
      call russ(zone,dm1,ir,0,1,lxname,lyname,ititl,index)
      lxname="n[ormalized] r[adius$"
      index="log(dm1) vs normalized radius"
      call russ(rn,dm1,ir,0,1,lxname,lyname,ititl,index)
      lxname="z[one]$"
      lyname="o[pacity (cm'eh.5[2'exhx[/gm)$"
      index="akap vs zone ... opacity vs zone"
      call russ(zone,akap,ir,0,0,lxname,lyname,ititl,index)
      lyname="p[ressure (dynes/cm'eh.5[2'exhx[)$"
      index="log(pressure) vs zone "
      call russ(zone,p,ir,0,1,lxname,lyname,ititl,index)
      lyname="d[ln(pressure)/]d[ln(density)$"
      index="g1 vs zone ... Dln(pressure)/Dln(density)"
      call russ(zone,g1,ir,0,0,lxname,lyname,ititl,index)
      lyname="d[ln(temperature)/]d[ln(density)$"
      index="g3m1 vs zone ... Dln(temperature)/Dln(density)"
      call russ(zone,g3m1,ir,0,0,lxname,lyname,ititl,index)
      lyname="l[uminosity] r[adiation] f[raction$"
      index="frft vs zone ... luminosity radiation fraction"
      call russ(zone,frft,irp1,0,0,lxname,lyname,ititl,index)
      lyname="e[nergy] s[ource (erg/gm/sec)$"
      index="log(sorce) vs zone "
c       inner and outer burning zones are n1 and n2
      do 10 i=1,ir
            n2= ir+1-i
            if( sorce(n2) .ne. zero ) goto 15
  10  continue
      go to 9
  15  continue
      do 11 i=1,ir
      n1=i
      if(sorce(i).gt.zero) go to 16
  11  continue
  16  continue
      n=n2-n1+1
8     if(n.le.1) go to 9
      eps=amin1(sorce(n1),sorce(n2))
      do 13 i=n1,n2
      if(sorce(i).eq.0.) go to 13
      eps=amin1(eps,sorce(i))
  13  continue
      do 12 i=n1,n2
      if(sorce(i).le.0.) sorce(i)=eps
  12  continue
      call russ(zone(n1),sorce(n1),n,0,1,lxname,lyname,ititl,index)
      lyname="#q[energy source/#q]t[ (erg/gm/s/kelvin)$"
      index="dtsorc vs zone ... D(energy source)/D(temperature)"
      call russ(zone,dtsorc,n2,0,0,lxname,lyname,ititl,index)
      lyname="#q[energy source/#q]v[ (erg/sec/cc)$"
      index="dvsorc vs zone ... D(energy source)/D(specific volume)"
      call russ(zone,dvsorc,n2,0,0,lxname,lyname,ititl,index)
9     lyname="n'eh.5[2'exhx] (gm/r'eh.5[3'exhx[)$"
      index="Brunt-Vaissala frequency squared vs Zone"
      do 28 i=1,ir
      n2=ir-i
      if(brntva(n2).gt.0..and.brntva(n2+1).le.0.) go to 27
  28  continue
  27  do 26 i=1,ir
      n1=i
      if(brntva(n1).gt.0.) go to 25
  26  continue
  25  n=n2-n1+1
c     call russ(zone(n1),brntva(n1),n,0,0,lxname,lyname,ititl,index)
      con=totmas/rphoto**3*grav
      frl1ocon=rl*(rl+one)/con*for
      do 18 i=1,ir
      brntva(i)=frl1ocon*g1(i)*p(i)*vp(i)/(rp(i)+rp(i+1))**2
  18  continue
      lyname="l[amb] f[requency] s[quared$"
      index="Lamb Frequency squared vs Zone"
c     call russ(zone,brntva,ir,0,1,lxname,lyname,ititl,index)
c
c          iout:= 0 for little output (no eigenvectors)
c                 1 for eigenvectors
c                 2 for eigenvectors, matrices and conver. info.
c
c          rl:= the laditudinal quantum number.
c
c          nomega:= number of points in frequency scan.
c
c          omlow:= starting frequency for scan (rad*sec**-1).
c
c          omhigh:= ending frequency for scan (rad*sec**-1).
c
c          (it is not necessary for omlow to be less than omhigh.)
c
  19  ldeg=ifix(rl*1.01)
      call close(10)
      if(ncc.eq.0) ncc=1
      ncca=3*ncc
      nccna=4*ncc
      if(feps.eq.0.) feps=1.
c     make sure omlow is less than omhigh
      if(omlow.lt.omhigh) go to 1
      temp=omlow
      omlow=omhigh
      omhigh=temp
1     continue
c     fsuba = one
      dm2(irp1) = dm1(ir)/(two*fsuba)
      write(6,2000) (ititl(i),i=1,3),ir
      write(6,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(6,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
      write(6,2005) iout,rl,nomega,omlow,omhigh
      write(6,2003)
      write(6,2004) (i,rp(i+1),tp(i),vp(i),p(i),g1(i),
     $      g3m1(i),akap(i),dkdr(i),dkdt(i),cv(i),
     $      frft(i+1),rm(i+1),i=1,ir)
c
c          check for nuclear burning and print out those
c       zones participating.
c
      if(n2.le.1) write(6,2100)
      if(n2.le.1) goto 20
      write(6,2101) (i,sorce(i),dtsorc(i),dvsorc(i),i=1,n2)
  20  continue
      if(feps.eq.1.) go to 6
      do 5 i=1,n2
      dtsorc(i)=dtsorc(i)*feps
      dvsorc(i)=dvsorc(i)*feps
    5 continue
    6 continue
      call pltdmp(irp1,0,iplot,np,ng,ldeg)
c     write input data to plot tape
      write(8,3001) irp1,8hrp
      write(8,3000) (rp(i),i=1,irp1)
      write(8,3001) ir,8htp
      write(8,3000) (tp(i),i=1,ir)
      write(8,3001) ir,8hvp
      write(8,3000) (vp(i),i=1,ir)
      write(8,3001) ir,8hcv
      write(8,3000) (cv(i),i=1,ir)
      write(8,3001) ir,8hdkdr
      write(8,3000) (dkdr(i),i=1,ir)
      write(8,3001) ir,8hdkdt
      write(8,3000) (dkdt(i),i=1,ir)
      write(8,3001) ir,8hdm1
      write(8,3000) (dm1(i),i=1,ir)
      write(8,3001) ir,8hakap
      write(8,3000) (akap(i),i=1,ir)
      write(8,3001) irp1,8hdm2
      write(8,3000) (dm2(i),i=1,irp1)
      write(8,3001) irp1,8hrm
      write(8,3000) (rm(i),i=1,irp1)
      write(8,3001) ir,8hp
      write(8,3000) (p(i),i=1,ir)
      write(8,3001) ir,8hg1
      write(8,3000) (g1(i),i=1,ir)
      write(8,3001) ir,8hg3m1
      write(8,3000) (g3m1(i),i=1,ir)
      write(8,3001) irp1,8hfrft
      write(8,3000) (frft(i),i=1,irp1)
      write(8,3001) ir,8hsorce
      write(8,3000) (sorce(i),i=1,ir)
      write(8,3001) ir,8hdtsorc
      write(8,3000) (dtsorc(i),i=1,ir)
      write(8,3001) ir,8hdvsorc
      write(8,3000) (dvsorc(i),i=1,ir)
c     write(8,3001) ir,8hbrntva
c     write(8,3000) (brntva(i),i=1,ir)
3000  format(1p6e13.5)
3001  format(i4,1x,10h(1p6e13.5),25x,a8)
c     end of input data to plot tape
      do 999 i=1,ir
      sound(i)=sqrt(g1(i)*p(i)*vp(i))
  999 continue
      gnu0=zero
      ir1=ir-1
      capa(1)=(sound(3)-sound(2))*4./(rp(4)-rp(2))
      do 998 i=1,ir1
      gnu0=(rp(i+2)-rp(i+1))/sound(i+1)+gnu0
      capa(i+1)=(sound(i+1)-sound(i))/rp(i+1)+capa(i)
 998  continue
      gnu0=0.5/gnu0
      capanu0=(sound(ir)/rp(ir)-capa(ir))/pi4/pi
      write(6,199) gnu0,capanu0
  199 format(1p2e11.3)
      call lnanon(iout,rl,nomega,omlow,omhigh)
        call st
      if(ip.gt.0) call ps
      call donepl
      stop
1000  format(1p4e20.13)
1001  format(3a8)
2000  format(50h1 initial model in hydrostatic and thermal balance,//,
     $      3a8,1x,i4,6h zones,/ )
2001  format(4x,18h total luminosity=,1pe11.4,12h total mass=,e11.4,
     $ 17h effective temp.=,e11.4,/,6x,21h photospheric radius=e11.4,
     $ 13h core radius=,e11.4,11h core lum.=e11.4)
2002  format(6h irad=,i2,8h noburn=,i2,8h onemq0=,1pe10.3,
     $ 8h onemq1=,e10.3,7h rlums=,e10.3,7h rlumc=,e10.3)
2003  format(3x,1hi,3x,6hradius,5x,4htemp,5x,6hsp.vol,3x,8hpressure,
     $  5x,2hg1,7x,4hg3m1,5x,7hopacity,3x,6hdlkdlr,4x,6hdlkdlt,6x,
     $  2hcv,7x,4hfrft,4x,8hint.mass)
2004  format(1x,i4,1p12e10.3)
2005  format(6h iout=,i2,13h l-value used,f7.2,
     $ 25h number of points in scan,i4,15h starting freq.,
     $ 1pe10.3,13h ending freq.,e10.3)
2100  format(/,30h model has no nuclear burning.)
2101  format(/,29h following zones are ignited.,/,
     $   (3(1x,i4,1p3e10.3)) )
      end
      block data
      implicit real(a-h,o-z)
c       fundamental constants are from novotny, introduction to
c    stellar atmospheres and interiors, 1973, appendix ii.
c                                                 2/16/83  wdp
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/ grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/thermo/ r,a,bk,avagd,ad3
      data r,a,bk,avagd,ad3 / 8.31434e7,7.56471e-15,8.6170837e-5,
     $     6.02217e23,2.52157e-15 /
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6723e-8,7.5589e-5,
     $  5.66870e-5,3.1415926536e0,6.2831853072e0,12.5663706144e0,
     $   25.1327412288e0,4.1887902048e0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0e0,1.0e0,2.0e0,3.0e0,4.0e0,10.0e0,
     $      0.5e0,0.25e0      /
      end
      subroutine lnanon(iout,rl,nomega,omlow,omhigh)
      complex fomega,cx,comega,comeg1,comeg2,cerr,cerra,ceta,dror,
     $      drho,dtemp,dent,dlum,dpress,dgam,cy,cz,rcont,hcont,domega,
     $      domsq,czero,ci
c
c          stability analysis of linear, non-adiabatic nonradial
c      stellar pulsations. first version has the non-adiabatic
c      effects included in the momentum equation but the entropy
c      change will be zero everywhere for now. the notation is
c      based on the lagrangian variations of castor ap. j. 162(166)1971,
c      and the actual linearization is based on the technique of
c      lynden-bell and ostriker.
c
c                                                9/15/83 wd pesnell
c
      common/phypar/ r(2001),t(2001),v(2001),cv(2001),dlkdlr(2001),
     $ dlkdlt(2001),dm1(2001),gkap(2001),dm2(2001),rm(2001),emdldm(2001)
      common/blk1/   irad,nflag,n,iplot,irp1,isum
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/blk4/   p(2001),g1(2001),g3m1(2001),rzone(2001)
      common/blk37/  drdm(2001,2),dl1(2001,2),dl2(2001,2),
     $      gor(2001),rho(2001),adrho(2001),adtemp(2001)
      common/coolum/ noburn,onemq0,onemq1,rlums,rlumc
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/scrtch/ ag1(2001,3),ag2(2001,2),ag3(2001,2),ag4(2001,2),
     $          ah1(2001,2),ah2(2001),  ah3(2001),  ah4(2001),
     $          ap1(2001,2),ap3(2001,3),ap4(2001)
      common/lumins/ frft(2001),sorce(2001),dtsorc(2001),dvsorc(2001)
      common/linear/ ceta(2001),dror(2001),drho(2001),dtemp(2001),
     $      dent(2001),dlum(2001),dgam(2001),cy(2001),cz(2001),
     $      dpress(2001),work(2001),dr(2001),eta(2001),gam(2001),
     $      y(2001),z(2001),weight(2001),spac(2001,77)
      common ak1(2001,4),ak2(2001,3),ak4(2001,3)
      common /sumtab/ nrq(100),npn(100),ngn(100),pa(100),pna(100)
     1 ,pgr(100),wgr(100),tke(100),iroot
      common /rbk/ rc(2001),zone(2001),rn(2001),ltop(3),rcn(2001),
     2 dnsity(2001)
      common/plot/ plotc(2001),plotd(2001)
      common/core/ncca,nccna,termax,err
      character *100,lxname,lyname,lhead,ides
      dimension drz(2001),drint(2001),fm(30,6),rkl2(2001),cvt(2001)
      dimension cx(8004),x(6003),stomeg(100),omsql(100),omsqu(100)
      data accur/1.e-10/,czero/(0.0e0,0.0e0)/,ci/(0.0e0,1.0e0)/
      data fm/17*0.0,2629.46,2765.03,2899.14,3034.02,3169.40,3303.77,
     1 3439.81,6*0.0,13*0.0,2161.,2293.,2427.,2559.20,2693.68,2828.73,
     2 2963.5,3098.66,3233.43,3369.50,3505.05,3642.47,3779.,3917.,4058.,
     3 2*0.0,12*0.0,2083.,2222.,2352.,2487.,2619.65,2754.96,2890.04,
     4 3024.54,3159.99,3295.23,3433.17,3567.,3703.,3840.,3980.,3*0.0,
     5 14*0.0,2408.,2542.,2677.,2812.,2948.,3082.,3219.,3354.,3491.,
     6 3628.,3763.,3904.,4*0.0,10*0.0,1918.,2052.,0.0,2323.,2459.,0.0,
     7 2728.,2866.,3000.,3136.,3272.,3409.,3545.,3683.,3821.,3962.,
     8 4*0.0,10*0.0,1963.,2100.,2235.,2371.,2505.,2641.,2777.,2914.,
     9 3050.,3187.,3324.,3460.,3597.,3733.,3871.,5*0.0/
      nmax = (2001)
      rlumx = corlum
c
c          define the zone-centered radii and the differences
c       of the two radii (zone and interface centered) for ease
c       of future work. the array gor is the local gravity
c       divided by the radius squared and rho is the zone density.
c
      ll=ifix(rl+.001)
      np1 = n + 1
      rhobar = totmas/(r(np1)**3*pi43)
      onthrd = one/thre
      rl1 = rl*(rl+one)
      do 10 i=1,n
            rzone(i) = (r(i)**3+ahf*dm1(i)*v(i)/pi43)**onthrd
      rkl2(i)=rl1/rzone(i)**2
            if(i.gt.1) drz(i) = rzone(i)-rzone(i-1)
            drint(i) = sqrt(r(i+1)-r(i))
            drdm(i,1) = forpi*r(i)**2/(dm1(i)*v(i))
            drdm(i,2) =-forpi*r(i+1)**2/(dm1(i)*v(i))
            dlkdlr(i) = dlkdlr(i)+g3m1(i)*dlkdlt(i)
            if(i.gt.1) gor(i) = g*rm(i)/r(i)**4
            rho(i) = one/v(i)
  10  continue
      drz(np1)=rzone(np1)-rzone(n)
      drdm(np1,1) = zero
      drdm(np1,2) = zero
      gor(np1) = g*rm(np1)/r(np1)**4
      pi4g = forpi*g
c
c          initialize the matrices.
c
c          all matrices are prefixed with the letter a, the second
c       letter denotes the equation ( g=radial component of the
c       momentum equation, h=horizontal component of the momentum
c       equation, p=poisson equation and k= thermal equation)
c       in which each matrix is used. the number tells which
c       variable is to be acted on by the matrix: 1= radial
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
c
c          inner boundary conditions.
c
      ag1(1,1) = zero
      ag1(1,2) = zero
      ag1(1,3) = zero
      ag2(1,1) = zero
      ag2(1,2) = zero
      ag3(1,1) = zero
      ag3(1,2) = zero
      ag4(1,1) = zero
      ag4(1,2) = zero
      ap1(1,1) = zero
      ap1(1,2) =-drint(1)*pi4g*rzone(1)**2*rho(1)*(drdm(1,2)-
     $  forpi*r(2)**2*ahf*(rho(2)-rho(1))/dm2(2))
      ap3(1,1) = zero
      ap3(1,3) = r(2)**2/(drz(2)*drint(2)*drint(1))
      ap3(1,2) =-(r(2)**2/(drz(2)*drint(2)**2)+rl*rzone(1)/drint(2)**2+
     $          rl1)
      ap4(1)   = -drint(1)*rl1*pi4g*rho(1)
      dl1(1,1) = zero
      dl1(1,2) = zero
      dl2(1,1) = zero
      dl2(1,2) = zero
      do 20 i=1,n
c
c          horizontal momentum equation.
c
            ah1(i,1) = drdm(i,1)*v(i)*g1(i)*p(i) + ahf*gor(i)*r(i)**2
            ah1(i,2) = drdm(i,2)*v(i)*g1(i)*p(i)+ahf*gor(i+1)*r(i+1)**2
            ah2(i) = g3m1(i)
            ah3(i) = one/drint(i)
            ah4(i) = g1(i)*p(i)*rl1/(rho(i)*rzone(i)**2)
            if( i .eq. n ) goto 21
c
c          radial momentum equation.
c
      afac = forpi*r(i+1)**2/dm2(i+1)
            ag1(i+1,1) =-afac*drdm(i,1)*g1(i)*p(i)
            ag1(i+1,2) =-for*gor(i+1)*r(i+1)+pi4g*ahf*(rho(i+1)+rho(i))+
     $ afac*(drdm(i+1,1)*g1(i+1)*p(i+1)-drdm(i,2)*g1(i)*p(i) )
            ag1(i+1,3) = afac*drdm(i+1,2)*g1(i+1)*p(i+1)
            ag2(i+1,1) =-afac*rho(i)*g3m1(i)
            ag2(i+1,2) = afac*rho(i+1)*g3m1(i+1)
            ag3(i+1,1) =-one/(drz(i+1)*drint(i))
            ag3(i+1,2) = one/(drz(i+1)*drint(i+1))
            ag4(i+1,1) = rl1*(gor(i+1)*ahf-afac*g1(i)*p(i)/rzone(i)**2)
            ag4(i+1,2) = rl1*gor(i+1)*ahf +
     $  afac*rl1*g1(i+1)*p(i+1)/rzone(i+1)**2
            if( i .eq. 1 ) goto 20
c
c          poisson equation.
c
            ap1(i,1) = -drint(i)*pi4g*rzone(i)**2*rho(i)*(drdm(i,1)-
     $ forpi*r(i)**2*ahf*(rho(i)-rho(i-1))/dm2(i))
            ap1(i,2) =-drint(i)*pi4g*rzone(i)**2*rho(i)*(drdm(i,2)-
     $  forpi*r(i+1)**2*ahf*(rho(i+1)-rho(i))/dm2(i+1))
            ap3(i,1) = r(i)**2/(drint(i)*drint(i-1)*drz(i))
            ap3(i,3) = r(i+1)**2/(drint(i)*drint(i+1)*drz(i+1))
            ap3(i,2) =-(r(i)**2/(drint(i)**2*drz(i)) + rl1 +
     $          r(i+1)**2/(drint(i)**2*drz(i+1)) )
            ap4(i)   =-drint(i)*pi4g*rl1*rho(i)
  21  continue
c
c          the thermal equation matrices.
c
      dtl1 = t(i-1)**4
      dtl2 = t(i)**4
      t4oki = dtl2/gkap(i)
      t4oki1 = dtl1/gkap(i-1)
      diff = t4oki-t4oki1
      wiowi1 = alog(dtl2/dtl1)
      wowsq = wiowi1**2
      gkogk1 = alog(gkap(i)/gkap(i-1))
      denom = one-gkogk1/wiowi1
      dl1(i,1) = (-for*g3m1(i-1)+dlkdlr(i-1))*(t4oki1/diff)
     $  + (-dlkdlr(i-1)/wiowi1+for*g3m1(i-1)*gkogk1/wowsq)/denom
      dl1(i,2) = (for*g3m1(i)-dlkdlr(i))*(t4oki/diff)
     $  + (dlkdlr(i)/wiowi1-for*g3m1(i)*gkogk1/wowsq)/denom
      dl2(i,1) = ((-for+dlkdlt(i-1))*(t4oki1/diff) -
     $  (dlkdlt(i-1)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i-1)*t(i-1))
      dl2(i,2) = ((for-dlkdlt(i))*(t4oki/diff) +
     $  (dlkdlt(i)/wiowi1-for*gkogk1/wowsq)/denom)/(cv(i)*t(i))
      glom = rlumx*frft(i-1)/dm1(i-1)
      rlumx1= rlumx + sorce(i-1)*dm1(i-1)
      if( noburn .eq. 1 ) rlumx1= cool(rm(i)/totmas)
      glom1 = rlumx1*frft(i)/dm1(i-1)
       dldm=(rlumx1-rlumx)/dm1(i-1)
       aconlm=ahf*((one-frft(i-1))*rlumx+(one-frft(i))*rlumx1)
      rlumx=rlumx1
c
c      convert the nuclear energy derivatives
c
      dedr=(-dvsorc(i-1)*v(i-1)+g3m1(i-1)*dtsorc(i-1)*t(i-1))
      deds=dtsorc(i-1)/cv(i-1)
      if(i.gt.2) go to 22
      ak1(i-1,1)=zero
      ak1(i-1,2)=zero
      ak1(i-1,3)=glom*dl1(i-1,2)*drdm(i-1,2)-glom1*for/r(i)-
     1glom1*(dl1(i,1)*drdm(i-1,2)+dl1(i,2)*drdm(i,1))
     2+dedr*drdm(i-1,2)
      ak1(i-1,4)=-glom1*dl1(i,2)*drdm(i,2)
      ak4(i-1,1)=zero
      ak4(i-1,2)=-rkl2(i-1)*(glom1*dl1(i,1)-dedr)
      ak4(i-1,3)=-rkl2(i)*glom1*dl1(i,2)
      ak2(i-1,1)=zero
      ak2(i-1,2)=-glom1*dl2(i,1)+deds
      ak2(i-1,3)=-glom1*dl2(i,2)
c
c         add horizontal heat flow terms
c
      horiz=rl1*for*ac3*t4oki1*(v(i-1)/rzone(i-1))**2
      ak1(i-1,2)=ak1(i-1,2)-horiz*g3m1(i-1)*drdm(i-1,1)-
     1rkl2(i-1)*ahf*glom*dm1(i-1)*v(i-1)/(forpi*r(i-1)**2)
      ak1(i-1,3)=ak1(i-1,3)-horiz*g3m1(i-1)*drdm(i-1,2)-
     1rkl2(i-1)*ahf*glom1*dm1(i-1)*v(i-1)/(forpi*r(i)**2)
      ak2(i-1,2)=ak2(i-1,2)-horiz/(cv(i-1)*t(i-1))
      ak4(i-1,2)=ak4(i-1,2)+rkl2(i-1)*(dldm-horiz*g3m1(i-1)
     1-aconlm*(v(i-1)/(forpi*rzone(i-1)**3)))
c     ak4(i-1,2)=ak4(i-1,2)-rkl2(i-1)*(horiz*g3m1(i-1)+
c    1(glom1-glom))
      go to 20
  22  continue
      ak1(i-1,1)=glom*dl1(i-1,1)*drdm(i-2,1)
      if(i.eq.3) ak1(i-1,1)=zero
      ak1(i-1,2)=-glom1*dl1(i,1)*drdm(i-1,1)+glom*for/r(i-1)+
     1glom*(dl1(i-1,1)*drdm(i-2,2)+dl1(i-1,2)*drdm(i-1,1))
     2+dedr*drdm(i-1,1)
      ak1(i-1,3)=glom*dl1(i-1,2)*drdm(i-1,2)-glom1*for/r(i)-
     1glom1*(dl1(i,1)*drdm(i-1,2)+dl1(i,2)*drdm(i,1))
     2+dedr*drdm(i-1,2)
      ak1(i-1,4)=-glom1*dl1(i,2)*drdm(i,2)
      ak4(i-1,1)=glom*rkl2(i-2)*dl1(i-1,1)
      ak4(i-1,2)=-rkl2(i-1)*(glom1*dl1(i,1)-dedr
     1-glom*dl1(i-1,2))
      ak4(i-1,3)=-glom1*rkl2(i)*dl1(i,2)
      ak2(i-1,1)=glom*dl2(i-1,1)
      ak2(i-1,2)=glom*dl2(i-1,2)-glom1*dl2(i,1)+deds
      ak2(i-1,3)=-glom1*dl2(i,2)
c
c         add horizontal heat flow terms
c
      horiz=rkl2(i-1)*for*ac3*t4oki1*v(i-1)**2
      ak1(i-1,2)=ak1(i-1,2)-horiz*g3m1(i-1)*drdm(i-1,1)-
     1rkl2(i-1)*ahf*glom*dm1(i-1)*v(i-1)/(forpi*r(i-1)**2)
      ak1(i-1,3)=ak1(i-1,3)-horiz*g3m1(i-1)*drdm(i-1,2)-
     1rkl2(i-1)*ahf*glom1*dm1(i-1)*v(i-1)/(forpi*r(i)**2)
      ak2(i-1,2)=ak2(i-1,2)-horiz/(cv(i-1)*t(i-1))
      ak4(i-1,2)=ak4(i-1,2)+rkl2(i-1)*(dldm-horiz*g3m1(i-1)
     1-aconlm*(v(i-1)/(forpi*rzone(i-1)**3)))
c     ak4(i-1,2)=ak4(i-1,2)-rkl2(i-1)*(horiz*g3m1(i-1)+
c    1(glom1-glom))
  20  continue
      tau = gkap(n)*dm1(n)/(for*forpi*r(np1)**2)
      tau = 1.5e0*tau/(one+1.5e0*tau)
      fac = (r(np1)-rphoto)/r(np1)
      ef1 = (one+ahf*fac)/(one+fac)
      ef3 = two/(one+fac)
      dl1(np1,1) = ef1*(for*g3m1(n)-tau*dlkdlr(n))
      dl1(np1,2) = dl1(np1,1)
      dl2(np1,1) = ef1*(for-tau*dlkdlt(n))/(t(n)*cv(n))
      dl2(np1,2) = zero
      glom = rlumgv/dm1(n)
      ak1(n,1) = glom*dl1(n,1)*drdm(n-1,1)
      ak1(n,2) =-glom*(dl1(np1,1)*drdm(n,1)-for/r(n) -
     $  dl1(n,1)*drdm(n-1,2)-dl1(n,2)*drdm(n,1) )
      ak1(n,3) = glom*(dl1(n,2)*drdm(n,2)-dl1(np1,1)*drdm(n,2)-
     $  (ef3+tau*ef1)/r(np1) )
      ak1(n,4) = zero
      ak2(n,1) = glom*dl2(n,1)
      ak2(n,2) = glom*(dl2(n,2)-dl2(np1,1))
      ak2(n,3) = zero
      ak4(n,1) = glom*rkl2(n-1)*dl1(n,1)
      ak4(n,2) =-glom*rkl2(n)*(dl1(np1,1)-dl1(n,2))
      ak4(n,3) = zero
c
c        add horizontal heat flow terms
c
      horiz=rkl2(n)*for*ac3*t(n)**4*v(n)**2/gkap(n)
      ak1(n,2)=ak1(n,2)-horiz*g3m1(n)*drdm(n,1)-
     1rkl2(n)*ahf*frft(n)*rlumx*v(n)/(forpi*r(n)**2)
      ak1(n,3)=ak1(n,3)-horiz*g3m1(n)*drdm(n,2)-
     1rkl2(n)*ahf*frft(n+1)*rlumx*v(n)/(forpi*r(n+1)**2)
      ak2(n,2)=ak2(n,2)-horiz/(cv(n)*t(n))
      ak4(n,2)=ak4(n,2)-rkl2(n)*horiz*g3m1(n)
c
2600  format(1h1,//,10x,36h the pulsation matrices, ak1,ak2,ak4,/)
c
c          outer boundary conditions.
c
      afac = forpi*r(np1)**2/dm2(np1)
      g1p=g1(n)*p(n)-for*ad3*t(n)**4*g3m1(n)
      g3m1ro=g3m1(n)/v(n)-for*ad3*t(n)**3/cv(n)
      ag1(np1,1) =-afac*drdm(n,1)*g1p
      ag1(np1,2) =-for*gor(np1)*r(np1) + pi4g*ahf*rho(n) -
     $    afac*drdm(n,2)*g1p
      ag1(np1,3) = zero
      ag2(np1,1) =-afac*g3m1ro
      ag2(np1,2) = zero
      ag3(np1,1) =-(rl+one)/(drz(np1)*drint(n))
      ag3(np1,2) = zero
      ag4(np1,1) = rl1*gor(np1)*(r(np1)/rzone(n))**2*
     1(one-afac*g1p/(gor(np1)*r(np1)**2))
      ag4(np1,2) = zero
      ap1(n,1) = -drint(n)*pi4g*rzone(n)**2*rho(n)*(drdm(n,1)-
     $      forpi*r(n)**2*ahf*(rho(n)-rho(n-1))/dm2(n))
      ap1(n,2) =-drint(n)*pi4g*rzone(n)**2*rho(n)*(drdm(n,2)+
     $      forpi*r(np1)**2*rho(n)/(dm2(n)+dm2(np1)))
      ap3(n,1) = r(n)**2/(drint(n)*drint(n-1)*drz(n))
      ap3(n,2) =-(r(n)**2/(drint(n)**2*drz(n)) + rl1 +
     $      r(np1)*(rl+one)/drint(n)**2)
      ap3(n,3) = zero
      ap4(n) =-drint(n)*pi4g*rl1*rho(n)
      ah4(np1) = zero
c
c        write the matrices to file 6.
c
      if( iout .lt. 2 ) goto 29
      write(6,2200)
      do 27 i=1,np1
            write(6,2201) i,ag1(i,1),ag1(i,2),ag1(i,3),ag2(i,1),
     $   ag2(i,2),ag3(i,1),ag3(i,2),ag4(i,1),ag4(i,2),ah4(i)
  27  continue
      write(6,2300)
      do 28 i=1,n
            write(6,2201) i,ah1(i,1),ah1(i,2),ah2(i),ah3(i),ap1(i,1),
     $      ap1(i,2),ap3(i,1),ap3(i,2),ap3(i,3),ap4(i)
  28  continue
      write(6,2600)
      do 26 i=1,n
      write(6,2201) i,ak1(i,1),ak1(i,2),ak1(i,3),ak1(i,4),ak2(i,1),
     $  ak2(i,2),ak2(i,3),ak4(i,1),ak4(i,2),ak4(i,3)
  26  continue
  29  continue
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
      xnorm = r(np1)
      fomp=zero
      fomi = zero
      iroot = 0
      dom = (omhigh - omlow)/float(nomega - 1)
      do 30 i=1,nomega
         omega = omlow + float(i-1)*dom
         omsq  = omega*omega
         fomeg = fomsq(omsq,x,xnorm,nmax,n,ind)
         if( fomeg*fomi .ge. zero ) goto 35
         if((i.gt.2).and.(abs(fomi).gt.abs(fomp))) go to 34
         iroot = iroot + 1
         stomeg(iroot) = omega - dom*fomeg/(fomeg-fomi)
         omsql(iroot)=(omega-dom)**2
         omsqu(iroot)=omsq
      go to 35
  34  write(6,2700) omsq
 2700 format(27h pole in period scan, omsq=,1pe15.6)
  35        continue
         fomp=fomi
         fomi = fomeg
         if( iout .gt. 0 ) write(6,3500) i,omsq,fomeg
  30  continue
      if( iroot .eq. 0 ) goto 900
      write(6,3000) rl,(stomeg(i),i=1,iroot)
c
c          attempt to converge on each of the frequencies found in
c       discriminant search. a secant method to find the zeroes of
c       error function is used.
c
      omsq1 = zero
      erra  = zero
      if(isum.ne.0) write(11,7014) ll
      if(isum.ne.0) write(12,7014) ll
      if(isum.ne.0) write(13,7014) ll
      do 40 nroot=1,iroot
            omega = stomeg(nroot)
            omsq  = omega*omega
      accuru=accur
   39       do 45 icount=1,20
      if(omsq.gt.1.1*omsqu(nroot)) omsq=.5*(omsq1+omsqu(nroot))
      if(omsq.lt.0.9*omsql(nroot)) omsq=.5*(omsq1+omsql(nroot))
                  err = fomsq(omsq,x,xnorm,nmax,n,ind)
                  if( icount .le. 1 ) goto 44
                  afac = abs((omsq-omsq1)/omsq)
      eotm=err/termax
           if(iout.gt.0) write(6,4401) icount,omsq,afac,err,eotm,ind
                  if( afac .le. accuru.and.afacp.le.accuru) goto 47
      if(abs(eotm).lt.1.e-14) go to 47
      if(abs(err-erra)/abs(err).lt.1.e-13) go to 46
c                 omsq2 = (erra*omsq-err*omsq1)/(erra-err)
      delomsq=omsq-omsq1
      delerr=err-erra
      if(abs(delomsq/omsq).lt.1.e-8.or.abs(delerr/err).lt.1.e-8)go to 42
      ratio=delomsq/delerr
  42  domsq=err*ratio
      omsq2=omsq-domsq
                  omsq1 = omsq
                  erra  = err
                  omsq  = omsq2
      afacp=afac
                  goto 45
  44              continue
                  omsq1 = omsq
                  erra  = err
                  omsq  = omsq*(one + 2.e-7)
      afacp=1.
  45        continue
            if( omsq .ge. zero ) goto 46
43          write(6,4501)
            omsq =-omsq
  46        continue
            omega = sqrt(omsq)
      accuru=accuru*ten
            write(6,4500) nroot,omega,accuru
c           goto 40
      go to 39
  47        continue
      if(omsq.lt.zero) go to 43
c
c          calculate the adiabatic eigenvectors and weight
c       function associated with the mode. the various
c       pieces of the weight function are stored in local
c       variables wtherm, wgrav, wcross and wdiag for ease
c       of coding and the weight function as a function of
c       position is printed on file 6.
c
      wtherm = zero
      wcross = zero
      wdiag = zero
      rsum = zero
      hsum = zero
      stwait=zero
         dr(1) = zero
         x(3*n) = r(np1)
         do 48 i=1,n
             dr(i+1)= x(3*i)/r(i+1)
             eta(i) = x(3*i-1)/rzone(i)**2
         adrho(i) = drdm(i,1)*r(i)*dr(i)+drdm(i,2)*r(i+1)*dr(i+1) +
     $         rl1*eta(i)
             if(i.gt.1) gam(i) = x(3*i-2)*(r(i)/(g*rm(i)*drint(i)))
      adtemp(i)=g3m1(i)*adrho(i)
c
c          weight function.
c
          wtherm = g1(i)*p(i)*v(i)*adrho(i)**2*dm1(i)
          wgrav = rl1*eta(i)*x(3*i-2)*dm1(i)/drint(i) + x(3*i)*dm2(i+1)*
     $            (ag3(i+1,1)*x(3*i-2)+ag3(i+1,2)*x(3*i+1))
          wcross = rl1*gor(i+1)*x(3*i)*dm2(i+1)*(x(3*i-1)+x(3*i+2))
          wdiag = x(3*i)**2*dm2(i+1)*(pi4g*ahf*(rho(i)+rho(i+1))-
     $             for*gor(i+1)*r(i+1))
          rsum = rsum + x(3*i)**2*dm2(i+1)
          hsum = hsum + x(3*i-1)*eta(i)*dm1(i)
          weight(i) = wtherm + wgrav + wcross + wdiag
      stwait=stwait+weight(i)
            if( i .eq. 1 ) goto 48
          drz(i) = dl1(i,1)*adrho(i-1)+dl1(i,2)*adrho(i) + for*dr(i) -
     $         ahf*rl1*(eta(i)+eta(i-1))
  48     continue
        drz(np1) =+for-ahf*rl1*eta(n) + dl1(np1,1)*adrho(n)
      rke = rsum + rl1*hsum
      do 49 i=1,n
         weight(i) = weight(i)/(rke*omsq)
  49  continue
      stwait=stwait/rke
      qch=abs(omsq-stwait)/omsq
         omega = sqrt(omsq)
         index = knode(dr,eta,nmax,n,np,ng)
         period = twopi/omega
      call cvtm(n,period,rlumgv,cv,t,dm1,cvt,nmax)
         pdays = period/86400.e0
         fuhza=1.e6/period
         write(6,5500) nroot,omega,period,pdays,fuhza
         write(6,5503) index,np,ng
      write(6,5504) omsq,stwait,qch
         nrq(nroot)=index
         npn(nroot) = np
         ngn(nroot)=ng
         pa(nroot)=period
         if(iout .gt. 0 ) write(6,5502)
      if(iout.gt.0)write(6,5501) (i,dr(i+1),eta(i),adrho(i),adtemp(i),
     $  drz(i+1),weight(i),cvt(i),i=1,n)
      write(lyname,"(11hw[eight (p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      lxname="z[one$"
      write(ides,"(17hweight vs zone p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call russ(zone(1),weight(1),n,0,0,lxname,lyname,ltop,ides)
      call pltdmp(irp1,1,iplot,np,ng,ll)
c
c          non-adiabatic loop.
c
c
c          calculate the quasi-adiabatic stability coefficient.
c      this is essentially the same method as in the radial case.
c
c          the quasi-adiabatic entropy contributions, first
c       sum the horizontal components.
c
      y(1) =-ak4(1,2)*x(2)-ak4(1,3)*x(5)
      ido = n-1
      do 60 i=2,ido
         y(i) =-ak4(i,1)*x(3*i-4)-ak4(i,2)*x(3*i-1)-ak4(i,3)*x(3*i+2)
  60  continue
      y(n) =-ak4(n,1)*x(3*n-4)-ak4(n,2)*x(3*n-1)
c
c          radial contribution
c
      do 62 i=1,n
      cz(i)=czero
            jmin = max0(1,4-i)
            jmax = min0(4,n+3-i)
      sy=zero
            do 620 j=jmin,jmax
                  sy   = sy   - ak1(i,j)*x(3*(i+j-3))
  620 continue
      y(i)=y(i)+sy
  62  continue
      call trcsol(ak2,nmax,1,n,ci*omega,y,cz,dlum,cy)
      cz(n+1) = czero
c
c          the quasi-adiabatic stability coefficient can be
c       written as in castor's paper on radial pulsations with
c       addition of the horizontal integrals.
c
      rcont = czero
      hcont = czero
      hsum = zero
      rsum = zero
      do 63 i=1,n
       rcont=rcont+x(3*i)*(ag2(i+1,1)*cz(i)+ag2(i+1,2)*cz(i+1))*dm2(i+1)
            hcont = hcont + eta(i)*ah2(i)*cz(i)*dm1(i)
            hsum = hsum + x(3*i-1)*eta(i)*dm1(i)
            rsum = rsum + x(3*i)*x(3*i)*dm2(i+1)
  63  continue
      rke = rsum + rl1*hsum
      domsq = (rcont+rl1*hcont)/rke
      comega = csqrt(omsq + domsq)
      growth =-forpi*aimag(comega)/omega
      write(6,6000) comega,growth
      write(14,6000) comega,growth
      write(6,6001) domsq,rke,rcont,hcont,rsum,hsum
      write(14,6001) domsq,rke,rcont,hcont,rsum,hsum
c
c          if the quasi-adiabatic guess is more than 10 percent,
c       than use the old first guess of an arbitrary multple of
c       the adiabatic eigenfrequency.
c
      if( cabs(comega) .gt. 1.1e0*omega ) comega = omega*(1.e0,1.e-4)
c      comega = omega*(one+ci*1.e-7)
      comeg1 = czero
      cerra  = czero
      if(real(omsq).gt.omsqu(nroot))   omsqu(nroot)=real(omsq)
      if(real(omsq).lt.omsql(nroot)) omsql(nroot)=real(omsq)
      sru=sqrt(omsqu(nroot))*1.1
      srl=sqrt(omsql(nroot))*0.9
  50  do 55 icount=1,20
      romega=real(comega)
      criti=abs(1.e-3*accuru*romega/aimag(comega))
      if(romega.gt.sru) comega=.5*(comeg1+sru*(1.,.0001))
      if(romega.lt.srl) comega=.5*(comeg1+srl*(1.,.0001))
            cerr = fomega(comega,cx,xnorm,nmax,n)
            if( icount .le. 1 ) goto 54
c     if(real(cerr)*real(cerra).lt.0.) go to 51
c 52  continue
            afac = cabs((comega-comeg1)/comega)
          afaci=abs((aimag(comega)-aimag(comeg1))/aimag(comega))
            if(iout.gt.0) write(6,4400) icount,comega,afac,cerr,cerra
       if(afac.le.accuru.and.afacp.le.accuru.and.afaci.le.criti) goto 57
      if(cabs(cerr/termax).lt.1.e-14) go to 57
      if(cabs(cerr-cerra)/cabs(cerr).lt.1.e-13) go to 56
c           comeg2 = (cerra*comega-cerr*comeg1)/(cerra-cerr)
      comeg2=comega-cerr*(comega-comeg1)/(cerr-cerra)
            comeg1 = comega
            cerra = cerr
            comega = comeg2
      afacp=afac
            goto 55
c 51  fac1=abs(aimag(comega)/real(comega))+1.
c     sru=amax1(real(comega),real(comega1))*fac1
c     srl=amin1(real(comega),real(comega1))/fac1
c     go to 52
  54           continue
            comeg1 = comega
            cerra = cerr
            comega = comega*(one+2.e-7)
      afacp=1.
  55  continue
 56   accuru=accuru*ten
      write(6,4500) nroot,real(comega),accuru
c     goto 40
      go to 50
  57  continue
c
c          calculation of the non-adiabatic eigenfunctions.
c
      twork = zero
      dror(1) = czero
      cx(4*n) = xnorm
      rsum = zero
      hsum = zero
      d=zero
      do 70 i=1,n
         dror(i+1) = cx(4*i)/r(i+1)
         ceta(i) = cx(4*i-1)/rzone(i)**2
         if(i.ne.1) dgam(i) = cx(4*i-3)*(r(i)/(g*rm(i)*drint(i)))
         dent(i) = cx(4*i-2)/(t(i)*cv(i))
         drho(i) = drdm(i,1)*dror(i)*r(i)+drdm(i,2)*dror(i+1)*r(i+1)
     $            + rl1*ceta(i)
         dpress(i) = g1(i)*drho(i)+g3m1(i)*cx(4*i-2)/(p(i)*v(i))
         dtemp(i) = g3m1(i)*drho(i)+dent(i)
         rsum = rsum + conjg(cx(4*i))*cx(4*i)*dm2(i+1)
         hsum = hsum + conjg(cx(4*i-1))*ceta(i)*dm1(i)
         work(i) = -aimag(conjg(drho(i))*cx(4*i-2))*g3m1(i)*dm1(i)
         twork = twork + work(i)
         if( i .eq. 1 ) goto 70
          dlum(i) = dl1(i,1)*drho(i-1)+dl1(i,2)*drho(i) + for*dror(i) -
     $       ahf*rl1*(ceta(i)+ceta(i-1)) + dl2(i,1)*cx(4*i-6) +
     $       dl2(i,2)*cx(4*i-2)
      d=d+(g1(i)-one)*g3m1(i)*emdldm(i)*conjg(drho(i))*drho(i)*dm1(i)
  70  continue
      call pltdmp(irp1,2,iplot,np,ng,ll)
      call pltdmp(irp1,3,iplot,np,ng,ll)
c
c          vertical gravity work.
c
 
      write(6,7000) rl
      omega = real(comega)
      growth =-forpi*aimag(comega)/omega
      if(abs(-(growth/forpi)-1.e-4).gt.1.e-5) go to 71
      write(6,7004)
      write(6,4400) icount,sru,srl
  71  period = twopi/omega
      pdays  = period/86400.e0
      qdays  = pdays*sqrt(rhobar/1.41e0)
      write(6,7001) comega,growth
      pgr(nroot)=growth
      gp=pdays/growth
      rke = ahf*omega*omega*(rsum+rl1*hsum)
      d=-d/(rke*for)*forpi/omega
      write(6,7012) twork
      growth =+pi*twork/rke
      write(6,7002) rsum,hsum,rke,growth,d
c
c          work function niceties.
c
      y(np1) = zero
      do 75 ii=1,n
         i = np1 - ii
         work(i) = work(i)*pi/rke
         y(i) = y(i+1) + work(i)
c        use array dgam to store computed diagnostic ddel
         if(ii.eq.1.or.ii.eq.n) go to 75
      del=alog(t(ii+1)/t(ii))/alog(p(ii+1)/p(ii))
         dgam(ii)=-((dtemp(ii+1)-dtemp(ii))-del*(dpress(ii+1)-
     1dpress(ii)))/
     1alog(p(ii+1)/p(ii))
  75  continue
      write(lyname,"(29hr[e -{d[del (nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      lxname="z[one$"
      write(ides,"(39hRe -{d[del vs zone ... nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call seper(dgam,plotd,plotc,1,n)
      call russ(zone(2),plotd(2),n-2,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(29hi[m -{d[del (nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      write(ides,"(39hIm -{d[del vs zone ... nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call russ(zone(2),plotc(2),n-2,0,0,lxname,lyname,ltop,ides)
c     do 76 i=2,n
c      use array dgam to store computed diagnostic 4*dr/r/dl/l
c     dgam(i)=dror(i)/dlum(i)*4.
c  76 continue
c     write(lyname,"(35hr[e {d[r/r/{d[l/l (nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
c     write(ides,"(42hRe Dr/r/Dl/l  vs zone ... nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
c     call seper(dgam,plotd,plotc,1,n)
c     call russ(zone(2),plotd(2),n-1,0,0,lxname,lyname,ltop,ides)
c     write(lyname,"(35hi[m {d[r/r/{d[l/l (nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
c     write(ides,"(42hIm Dr/r/Dl/l  vs zone ... nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
c     call russ(zone(2),plotc(2),n-1,0,0,lxname,lyname,ltop,ides)
c     do 77 i=2,n
c      use array dgam to store computed diagnostic gradient
c     dgam(i)=(t(i)**4*dtemp(i)-t(i-1)**4*dtemp(i-1))/(t(i)**4
c    1-t(i-1)**4)/dlum(i)*4.
c 77  continue
c     write(lyname,"(35hr[e {d[t/t/{d[l/l (nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
c     write(ides,"(42hRe Dt/t/Dl/l  vs zone ... nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
c     call seper(dgam,plotd,plotc,1,n)
c     call russ(zone(2),plotd(2),n-1,0,0,lxname,lyname,ltop,ides)
c     write(lyname,"(35hi[m {d[t/t/{d[l/l (nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
c     write(ides,"(42hIm Dt/t/Dl/l  vs zone ... nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
c     call russ(zone(2),plotc(2),n-1,0,0,lxname,lyname,ltop,ides)
c     do 78 i=2,n
c      use array dgam to store horizontal density effect on dl/l
c     dgam(i)=rl1*(ceta(i)+ceta(i-1))/dlum(i)*ahf
c  78 continue
c     write(lyname,"(39hr[e {d[rho/rho/{d[l/l (nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
c     write(ides,"(46hRe Drho/rho/Dl/l  vs zone ... nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
c     call seper(dgam,plotd,plotc,1,n)
c     call russ(zone(2),plotd(2),n-1,0,0,lxname,lyname,ltop,ides)
c     write(lyname,"(39hi[m {d[rho/rho/{d[l/l (nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
c     write(ides,"(46hIm Drho/rho/Dl/l  vs zone ... nonadiabatic, p=,
c    1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
c     call russ(zone(2),plotc(2),n-1,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(29hr[e {d[temp (nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      write(ides,"(38hRe Dtemp  vs zone ... nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call seper(dtemp,plotd,plotc,1,n)
      call russ(zone(2),plotd(2),n-1,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(29hi[m {d[temp (nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      write(ides,"(38hIm Dtemp  vs zone ... nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call russ(zone(2),plotc(2),n-1,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(32hr[e {d[lum/lum (nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      write(ides,"(41hRe Dlum/lum  vs zone ... nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call seper(dlum,plotd,plotc,1,n)
      call russ(zone(2),plotd(2),n-1,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(32hi[m {d[lum/lum (nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2,2h)$)") np,ng,ll
      write(ides,"(41hIm Dlum/lum  vs zone ... nonadiabatic, p=,
     1i2,4h, g=,i2,4h, l=,i2)") np,ng,ll
      call russ(zone(2),plotc(2),n-1,0,0,lxname,lyname,ltop,ides)
      call workmax(y,n,nmax)
      call zie(work,n,nmax)
      tke(nroot)=rke
      wgr(nroot)=growth
      gw=pdays/growth
      fuhz=1.e6/period
      write(6,7003) period,pdays,qdays,fuhz
      pna(nroot)=period
      if(isum.ne.0) write(11,7013) np,period,gp,gw-gp,fuhz,fuhza-fuhz,
     1 fm(np,ll+1),fuhz-fm(np,ll+1)
      if(isum.ne.0) write(12,7013) np-1,period,gp,gw-gp,fuhz,fuhza-fuhz,
     1 fm(np-1,ll+1),fuhz-fm(np-1,ll+1)
      if(isum.ne.0) write(13,7013) np+1,period,gp,gw-gp,fuhz,fuhza-fuhz,
     1 fm(np+1,ll+1),fuhz-fm(np+1,ll+1)
      if( iout .le. 0 ) goto 40
      write(6,7010)
      write(6,7011) (i,dror(i+1),ceta(i),dgam(i),dent(i),drho(i),
     $                 work(i),i=1,n)
      write(6,7016) (i,dpress(i),dlum(i+1),dtemp(i),y(i),i=1,n)
  40  continue
      return
c
c
c          error recovery.
c
 900  continue
      write(6,9000) rl,omlow,omhigh
      stop
c
2200  format(1h1,//,20h ag1,ag2,ag3,ag4,ah4)
2201  format(1x,i5,1p10e12.4)
2300  format(1h1,//,24h ah1,ah2,ah3,ap1,ap2,ap4)
3000  format(23h roots for analysis: l=,f7.2,/,(3x,1p5e15.6) )
3500  format(1x,i5,1p2e15.6)
4400  format(i5,1p7e15.6)
 4401 format(i5,1p4e15.6,i5)
4500  format(28h no convergence root number ,i4,
     $       11h frequency=,1pe15.6,8h accuru=,e10.3)
4501  format(35h from lnanon...omsq less than zero.)
5500  format(1h1,//,5x,4hroot,i5,7h omega=,1pe15.6,8h period=,2e15.6,
     $ 13h,  frequency=,e12.5,6h (uhz))
 5501 format(1x,i5,1p7e15.6)
 5502 format(4x,1hi,6x,4hdr/r,12x,3heta,9x,8hdrho/rho,9x,11hdtemp/temp,
     $   5x,4hdl/l,7x,6hweight,5x,4hcvtm,/)
 5503 format(23h radial quantum number=,i5,25h number of nodes, p-type=,
     $   i5,8h g-type=,i5)
 5504 format(/,2x,32heigenvalue from matrix solution=,1pe12.3,
     1/,2x,32heigenvalue from weight function=,1pe12.3,7h error=,e12.3)
c
c          quasi-adiabatic formats.
c
 6000 format(28h quasi-adiabatic eigenvalue:,1p2e11.4,8h growth=,e11.4)
 6001 format(7h domsq=,1p2e11.4,5h rke=,e11.4,/,7h rcont=,2e11.4,
     $  7h hcont=,2e11.4,6h rsum=,e11.4,6h hsum=,e11.4)
c
c          non-adiabatic format statements.
c
7000  format(1h1,//,30h non-adiabatic solution for l=,f7.2)
7001  format(//,21h frequency (sec**-1)=,1pe15.6,2x,12h imag. part=,
     $       e15.6,2x,23h predicted growth rate=,e15.6)
7002  format(//,23h work integral results:,/,18h kinetic energies:,
     $       1pe15.6,9h (radial),e15.6,9h (horiz.),e15.6,
     $       22h (total), all in ergs.,/,
     $       32h growth rate from work integral=,e15.6,3h d=,e15.6)
7003  format(/,9h periods:,1pe15.6,7h (secs),e15.6,11h (days), q=,e15.6,
     1 13h,  frequency=,e12.5,6h (uhz))
7004  format(33h period solution may be defective)
7010  format(3x,1hi,10x,4hdr/r,19x,3heta,17x,5hgamma,15x,
     $       11hd (entropy),13x,8hdrho/rho,8x,4hwork)
7011  format( 1x,i4,1p11e11.3 )
7012  format(20h work contributions:,1p4e12.4)
7013  format(i3,f11.2,2f11.3,4f11.2)
7014  format(1h1,30x,17hsolar summary (l=,i2,1h),//,
     1 32x,4hwork,/,28x,8hintegral,34x,10hcalculated,/,
     2 16x,9hpredicted,5x,6hgrowth,13x,9hadiabatic,13x,9hfrequency,/,
     3 19x,6hgrowth,7x,4hrate,11h calculated,11h  frequency,
     3 11h   measured,
     4 11h  -measured,/,8x,6hperiod,7x,4hrate,2x,9hincrement,2x,
     4 9hfrequency,
     5 2x,9hincrement,2x,9hfrequency,2x,9hfrequency,/,2x,1hn,6x,5h(sec),
     6 5x,6h(days),5x,6h(days),6x,5h(uhz),6x,5h(uhz),6x,5h(uhz),6x,
     6 5h(uhz),/)
 7016 format(1x,i4,1p7e10.2)
 9000 format(43h from lnanon...no adiabatic roots found, l=,f7.2,/,
     $      17h frequency range:,1p2e15.6)
      end
      function cool(qc)
c
c    returns the value of luminosity out of a zone
c  as a linear ramp from q1 to q0
c
c     needs in common block coolum:
c            onemq0:=1.-q0 (inner mass fraction)
c            onemq1:=1.-q1 (outer mass fraction)
c            rlumc:=core luminosity
c            rlums:=surface luminosity
c      q is present mass fraction
c          cool:=f*rlums+(1-f)*rlumc
c               f between 0 and 1
c
      common/coolum/  noburn,onemq0,onemq1,rlums,rlumc
      if( qc .gt. 1.0e0-onemq1 ) goto 30
      if( qc .lt. 1.0e0-onemq0 ) goto 20
      f = (qc+onemq0-1.0e0)/(onemq0-onemq1)
      cool = f*rlums+(1.0e0-f)*rlumc
      return
  20  cool = rlumc
      return
  30  cool = rlums
      return
      end
      function knode(x,y,nmax,n,np,ng)
c
c          returns the number of zero-crossings in the counter-
c       clockwise direction as np and the clockwise crossings
c       in ng.
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(nmax),y(nmax)
      np = 0
      ng = 0
      knode = 0
      if( n .eq. 1 ) return
      nm1 = n - 1
      do 10 i=1,nm1
            if( x(i)*x(i+1) .ge. zero ) goto 10
      if((i.ne.1).and.(x(i-1)*x(i).lt.zero)) go to 10
      if((i.ne.nm1).and.(x(i+1)*x(i+2).lt.zero)) go to 10
            if( y(i)*y(i+1) .lt. zero ) write(6,1000)
            if( x(i+1) .lt. x(i) ) goto 15
            if( y(i) .ge. zero ) ng = ng + 1
            if( y(i) .lt. zero ) np = np + 1
            goto 10
  15        if( y(i+1) .ge. zero ) np = np + 1
            if( y(i+1) .lt. zero ) ng = ng + 1
  10  continue
      if( x(n-1)*x(n) .lt. zero ) write(6,1001)
      knode = np - ng
      return
1000  format(51h from knode...quadrant jump default rotations used.)
1001  format(52h from knode...node at outer boundary is not counted.)
      end
      function fomsq(omsq,x,xnorm,nmax,n,ind)
      common/scrtch/ ag1(2001,3),ag2(2001,2),ag3(2001,2),ag4(2001,2),
     $      ah1(2001,2),ah2(2001),ah3(2001),ah4(2001),ap1(2001,2),
     $      ap3(2001,3),ap4(2001)
      common/linear/atrix(6003,10),work(6003,11),iwork(6003)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      dimension x(6003)
      dimension ipvt(6003)
      common/core/ncca,nccna,termax,err
c
c       dimensions of spac in /linear/ are 83*nmax
c
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega**2.
c
      do 200 i=1,n
c
c      ix is index for radial component, ieta for horizontal and
c    igam is for the poisson equation.
c
            ix = 3*i
            ieta = ix - 1
            igam = ix - 2
c
c          radial component of motion.
c
            x(ix) = zero
            atrix(ix,1) = ag1(i+1,1)
            atrix(ix,2) = ag3(i+1,1)
            atrix(ix,3) = ag4(i+1,1)
            atrix(ix,4) = ag1(i+1,2) - omsq
            atrix(ix,5) = ag3(i+1,2)
            atrix(ix,6) = ag4(i+1,2)
            atrix(ix,7) = ag1(i+1,3)
c
c          horizontal component of motion.
c
            x(ieta) = zero
            atrix(ieta,1) = zero
            atrix(ieta,2) = ah1(i,1)
            atrix(ieta,3) = ah3(i)
            atrix(ieta,4) = ah4(i) - omsq
            atrix(ieta,5) = ah1(i,2)
            atrix(ieta,6) = zero
            atrix(ieta,7) = zero
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
      x(3*n-1) =-ah1(n,2)*xnorm
      x(3*n-2) =-ap1(n,2)*xnorm
      x(3*n-3) =-ag1(n,3)*xnorm
c     call rbmles(atrix,3*nmax,1,3*n-1,7,x)
      call snbir(atrix,3*nmax,3*n-1,3,3,x,1,ind,work,iwork)
c     err = ag1(n+1,1)*x(3*n-3) + (ag1(n+1,2)-omsq)*xnorm +
c    $      ag3(n+1,1)*x(3*n-2) + ag4(n+1,1)*x(3*n-1)
 14   term1=ag1(n+1,1)*x(3*n-3)
      term2=(ag1(n+1,2)-omsq)*xnorm
      term3=ag3(n+1,1)*x(3*n-2)
      term4=ag4(n+1,1)*x(3*n-1)
      err=term1+term2+term3+term4
      xncca=x(ncca)+x(3)
      termax=amax1(abs(term1),abs(term2),abs(term3),abs(term4))
     1/xncca
      fomsq = err/xncca
      return
      end
      subroutine rbmles(a,nmax,imin,imax,m,y)
      implicit real(a-h,o-z)
      dimension a(nmax,7),y(1)
      id = imax-imin+1
      idm = id-1
      mmid=(m+1)/2
      mmm=mmid-1
      do 20 ii=1,idm
      i = imax+1-ii
      den = 1.e0/a(i,mmid)
      y(i) = y(i)*den
      kmax=min0(i-imin,mmm)
      do 10 j=1,mmm
      a(i,j) = a(i,j)*den
      do 10 kk=1,kmax
      jk=j+kk
      k=mmid+kk
      ik = i-kk
      a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  10  continue
      do 20 kk=1,kmax
      ik = i-kk
      k=mmid+kk
20    y(ik) = y(ik)-y(i)*a(ik,k)
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
      kmax=min0(i-imin,mmm)
      do 30 kk=1,kmax
      ik = i-kk
      k=mmid-kk
      y(i) = y(i)-a(i,k)*y(ik)
30    continue
      return
      end
      subroutine trcsol(a,nmax,imin,imax,x,y,z,e,f)
      complex x,z,e,f,den
      dimension a(nmax,3),y(nmax),z(nmax),e(nmax),f(nmax)
c
c          tridiagonal matrix solver for complex arrays.
c
      ido = imax+1 - imin
      e(imax) =-a(imax,1)/(a(imax,2)-x)
      f(imax) = y(imax)/(a(imax,2)-x)
      do 10 ii=2,ido
            i = imax+1 - ii
            den = a(i,2)-x + a(i,3)*e(i+1)
            if( ii .eq. ido ) goto 10
            e(i) =-a(i,1)/den
  10        f(i) = (y(i)-a(i,3)*f(i+1))/den
      z(imin) = f(imin)
      ido = imin+1
      do 20 i=ido,imax
            z(i) = e(i)*z(i-1) + f(i)
  20  continue
      return
      end
      complex function fomega(omega,x,xnorm,nmax,n)
      complex atrix,x,omsq,omega,err,czero,term,term1,term2,term3,term4,
     1term5
      common/scrtch/ ag1(2001,3),ag2(2001,2),ag3(2001,2),ag4(2001,2),
     $      ah1(2001,2),ah2(2001),ah3(2001),ah4(2001),ap1(2001,2),
     $      ap3(2001,3),ap4(2001)
      common/linear/ atrix(8004,13)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/ grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/lumins/ frft(2001),sorce(2001),dtsorc(2001),dvsorc(2001)
      common/phypar/ r(2001),t(2001),v(2001),cv(2001),dlkdlr(2001),
     $ dlkdlt(2001),dm1(2001),gkap(2001),dm2(2001),rm(2001),emdldm(2001)
      common/blk4/   p(2001),g1(2001),g3m1(2001),rzone(2001)
      common/blk37/  drdm(2001,2),dl1(2001,2),dl2(2001,2),
     $      gor(2001),rho(2001),adrho(2001),adtemp(2001)
      common /sumtab/ nrq(100),npn(100),ngn(100),pa(100),pna(100)
     1 ,pgr(100),wgr(100),tke(100),iroot,rl,ititl(3)
      common ak1(2001,4),ak2(2001,3),ak4(2001,3)
      dimension x(8004)
      data czero/(0.0e0,0.0e0)/
      common/core/ncca,nccna,termax,err
c
c          evaluate the error in the outer boundary equation
c       for an arbitrary omega.
c
      omsq = omega*omega
c     term=(0.,1.)*870*real(omega)/pi2
      term=0.
      rl1=rl*(rl+one)
      do 200 i=1,n
c
c      ix is index for radial component, ieta for horizontal,
c    igam is for the poisson equation and ient is the entropy index.
c
            ix = 4*i
            ieta = ix - 1
            igam = ix - 3
            ient = ix - 2
c
c          radial component of motion.
c
            x(ix) = czero
            atrix(ix,1) = czero
            atrix(ix,2) = czero
            atrix(ix,3) = ag1(i+1,1)
            atrix(ix,4) = ag3(i+1,1)
            atrix(ix,5) = ag2(i+1,1)
            atrix(ix,6) = ag4(i+1,1)
            atrix(ix,7) = ag1(i+1,2) - omsq
            atrix(ix,8) = ag3(i+1,2)
            atrix(ix,9) = ag2(i+1,2)
            atrix(ix,10) = ag4(i+1,2)
            atrix(ix,11) = ag1(i+1,3)
            atrix(ix,12) = czero
            atrix(ix,13) = czero
c
c          horizontal component of motion.
c
            x(ieta) = czero
            atrix(ieta,1) = czero
            atrix(ieta,2) = czero
            atrix(ieta,3) = czero
            atrix(ieta,4) = ah1(i,1)
            atrix(ieta,5) = ah3(i)
            atrix(ieta,6) = ah2(i)
            atrix(ieta,7) = ah4(i) - omsq
            atrix(ieta,8) = ah1(i,2)
            atrix(ieta,9) = czero
            atrix(ieta,10) = czero
            atrix(ieta,11) = czero
            atrix(ieta,12) = czero
            atrix(ieta,13) = czero
c
c          poisson equation.
c
            x(igam) = czero
            atrix(igam,1) = czero
            atrix(igam,2) = czero
            atrix(igam,3) = ap3(i,1)
            atrix(igam,4) = czero
            atrix(igam,5) = czero
            atrix(igam,6) = ap1(i,1)
            atrix(igam,7) = ap3(i,2)
            atrix(igam,8) = czero
            atrix(igam,9) = ap4(i)
            atrix(igam,10) = ap1(i,2)
            atrix(igam,11) = ap3(i,3)
            atrix(igam,12) = czero
            atrix(igam,13) = czero
c
c          thermal equations.
c
      dedr=-dvsorc(i)*v(i)+g3m1(i)*dtsorc(i)*t(i)
      dedv12=dedr*drdm(i,1)
      dedv13=dedr*drdm(i,2)
      dedv4=rl1*dedr/rzone(i)**2
      deds=dtsorc(i)/cv(i)
            x(ient) = czero
            atrix(ient,1) = ak1(i,1)
            atrix(ient,2) = czero
            atrix(ient,3) = ak2(i,1)
            atrix(ient,4) = ak4(i,1)
            atrix(ient,5) = ak1(i,2)-term*dedv12
            atrix(ient,6) = czero
            atrix(ient,7) = ak2(i,2) - (0.0e0,1.0e0)*omega-term*deds
            atrix(ient,8) = ak4(i,2)-term*dedv4
            atrix(ient,9) = ak1(i,3)-term*dedv13
            atrix(ient,10) = czero
            atrix(ient,11) = ak2(i,3)
            atrix(ient,12) = ak4(i,3)
            atrix(ient,13) = ak1(i,4)
 200  continue
      x(4*n-1) =-ah1(n,2)*xnorm
      x(4*n-2) =-ak1(n,3)*xnorm
      x(4*n-3) =-ap1(n,2)*xnorm
      x(4*n-4) =-ag1(n,3)*xnorm
      x(4*n-6) =-ak1(n-1,4)*xnorm
      call cbmles(atrix,4*nmax,1,4*n-1,13,x)
c     err = ag1(n+1,1)*x(4*n-4) + (ag1(n+1,2)-omsq)*xnorm +
c    $ ag3(n+1,1)*x(4*n-3)+ag4(n+1,1)*x(4*n-1)+ag2(n+1,1)*x(4*n-2)
      term1=ag1(n+1,1)*x(4*n-4)
      term2=(ag1(n+1,2)-omsq)*xnorm
      term3=ag3(n+1,1)*x(4*n-3)
      term4=ag4(n+1,1)*x(4*n-1)
      term5=ag2(n+1,1)*x(4*n-2)
      err=term1+term2+term3+term4+term5
      termax=amax1(cabs(term1),cabs(term2),cabs(term3),cabs(term4),
     1cabs(term5))/cabs(x(nccna))
      fomega = err/x(nccna)
      return
      end
      subroutine cbmles(a,nmax,imin,imax,m,y)
      implicit real(a-h,o-z)
      complex a,y,den
      dimension a(nmax,13),y(1)
      id = imax-imin+1
      idm = id-1
      mmid=(m+1)/2
      mmm=mmid-1
      do 20 ii=1,idm
      i = imax+1-ii
      den = 1.e0/a(i,mmid)
      y(i) = y(i)*den
      kmax=min0(i-imin,mmm)
      do 10 j=1,mmm
      a(i,j) = a(i,j)*den
      do 10 kk=1,kmax
      jk=j+kk
      k=mmid+kk
      ik = i-kk
      a(ik,jk) = a(ik,jk)-a(i,j)*a(ik,k)
  10  continue
      do 20 kk=1,kmax
      ik = i-kk
      k=mmid+kk
20    y(ik) = y(ik)-y(i)*a(ik,k)
      y(imin) = y(imin)/a(imin,mmid)
      imp = imin+1
      do 30 i=imp,imax
      kmax=min0(i-imin,mmm)
      do 30 kk=1,kmax
      ik = i-kk
      k=mmid-kk
      y(i) = y(i)-a(i,k)*y(ik)
30    continue
      return
      end
      subroutine cvtm(n,p,rl,cv,t,dm,x,nmax)
c
c        compute x=sum(cv*t*dm/(rl*p)) from outside in
c
      dimension cv(nmax),t(nmax),dm(nmax),x(nmax)
      x(n)=cv(n)*t(n)*dm(n)/(rl*p)
      ido=n-1
      do 10 i=1,ido
      j=n-i
      x(j)=x(j+1)+cv(j)*t(j)*dm(j)/(rl*p)
  10  continue
      return
      end
      subroutine zie(t,n,nzmax)
      dimension t(nzmax)
      dimension imax(2000), imin(2000)
c
c find minima and maxima in work per zone function
c
      nmax = 0
      nmin = 0
      nm2 = n-2
      do 50 i=3,nm2
         max = 0
         if( t(i) .gt. t(i-2) ) max = max+1
         if( t(i) .gt. t(i-1) ) max = max+1
         if( t(i) .gt. t(i+1) ) max = max+1
         if( t(i) .gt. t(i+2) ) max = max+1
         if( max .eq. 0 ) goto 40
         if( max .lt. 4 ) goto 50
c maximum found
      nmax=nmax+1
         if( nmax .gt. 1000 ) goto 50
      imax(nmax) = i
      goto 50
c minimum found
   40 nmin=nmin+1
         if( nmin .gt. 1000 ) goto 50
      imin(nmin) = i
   50 continue
      write(11,60) nmax
      if( nmax .le. 0 ) goto 55
      write(11,70) (imax(i),t(imax(i)),i=1,nmax)
   55 write(11,65) nmin
      if( nmin .le. 0 ) return
      write(11,70) (imin(i),t(imin(i)),i=1,nmin)
      return
c
   60 format(/,1x,i2,40h maxima found in work per zone function ,
     $ 42hwhere zone(i) .gt. (i-2),(i-1),(i+1),(i+2),
     $ 29h...listed below by zone,value)
   65 format(/,1x,i2,40h minima found in work per zone function ,
     $ 42hwhere zone(i) .le. (i-2),(i-1),(i+1),(i+2),
     $ 29h...listed below by zone,value)
   70 format(8(1x,i4,1pe10.3))
      end
      subroutine workmax(w,nz,nzmax)
c
c
c find and print all maxima in normalized work function
c
      dimension w(1),imax(2000)
      nmax=0
      nz1=nz-1
      do 10 i=2,nz1
         if( w(i).le.w(i-1) .or. w(i).le.w(i+1) ) goto 10
c
c maximum found... record it
c
         nmax = nmax+1
         imax(nmax) = i
   10 continue
      write(11,20) nmax
      if(nmax.le.0) return
      write(11,30) (imax(i),w(imax(i)),i=1,nmax)
      return
   20 format(/,1x,i2,41h maxima found in normalized work function,
     $  29h...listed below by zone,value)
   30 format(8(1x,i3,1pe10.3))
      end
      subroutine pltdmp(n,nonfl,iplot,np,ng,ldeg)
c
c        does various writes to the plot file.
c             iplot= 0 no plotting
c                   1 r(i)/r(n),y1,cy1(real),cy1(imag)
c                   2 r(i)/r(n),y1,y2,y3,y4,cy1(r&i)
c                   3 r(i)/r(n),y1,cy1(r&i)++++cy6(r&i)
c                   4 everybody
c
c             nonfl:= 0 write r(i)/r(n2)
c                     1 write adiabatic eigenvectors
c                     2 write non-adiabatic eigenvectors
c                     3 write work per zone
c
      complex fomega,cx,comega,comeg1,comeg2,cerr,cerra,ceta,dror,
     $      drho,dtemp,dent,dlum,dpress,dgam,cy,cz,rcont,hcont,domega,
     $      czero,ci
      common/phypar/ rp(2001),tp(2001),vp(2001),cv(2001),dkdr(2001),
     $   dkdt(2001),dm1(2001),akap(2001),dm2(2001),rm(2001),emdldm(2001)
      common/linear/ ceta(2001),dror(2001),drho(2001),dtemp(2001),
     $      dent(2001),dlum(2001),dgam(2001),cy(2001),cz(2001),
     $      dpress(2001),work(2001),dr(2001),eta(2001),gam(2001),
     $      y(2001),z(2001),spac(2001,78)
      common/blk37/  drdm(2001,2),dl1(2001,2),dl2(2001,2),
     $      gor(2001),rho(2001),adrho(2001),adtemp(2001)
      common/plot/ plotc(2001),plotd(2001)
      common /rbk/ rc(2001),zone(2001),rn(2001),ltop(3),rcn(2001),
     2 dnsity(2001)
      character *100,lxname,lyname,lhead,ides
      character *10 ititl(20)
c
c          titles for dumps
c
      data ititl/10h norm. rad ,10h dr/r      ,10h eta       ,
     &           10h drho/rho  ,10h gamma     ,10h dr/r   re ,
     &           10h dr/r   im ,10h eta    re ,10h eta    im ,
     &           10h gamma  re ,10h gamma  im ,10h dentro re ,
     &           10h dentro im ,10h drho/r re ,10h drho/r im ,
     &           10h cy6(i) re ,10h cy6(i) im ,10h work     ,
     &           10hcum. work  ,10hzone        /
      if( iplot .eq. 0 ) return
      n1=n-1
      if( nonfl .eq. 1 ) go to 10
      if( nonfl .eq. 2 ) go to 20
      if( nonfl .eq. 3 ) go to 30
c
c          normalized radius
c
      do 1 i=1,n
            plotd(i) = rp(i)/rp(n)
            plotc(i) = float(i)
  1   continue
      write(8,1001) n,ititl(1)
      write(8,1000) (plotd(i),i=1,n)
      write(8,1001) n,ititl(20)
      write(8,1000) (plotc(i),i=1,n)
      return
c
c          adiabatic eigenvectors
c
  10  continue
      write(8,1001) n,ititl(2)
      write(8,1000) (dr(i) ,i=1,n)
      write(15,21)
21    format(1h )
      write(lyname,"(21h{d[r/r (adiabatic, p=,i2,4h, g=,i2,4h, l=,i2,2h)
     1$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(30hdr/r vs zone ... adiabatic, p=,i2,4h, g=,i2,4h, l=
     1,i1)") np,ng,ldeg
      call russ(zone,dr,n,0,0,lxname,lyname,ltop,ides)
      lxname="n[ormalized] r[adius$"
      write(ides,"(43hdr/r vs normalized radius ... adiabatic, p=,i2,4h,
     1 g=,i2,4h, l=,i2)") np,ng,ldeg
      call russ(rn,dr,n,0,0,lxname,lyname,ltop,ides)
      if( iplot .eq. 1 ) return
      if( iplot .eq. 3 ) return
c
c          remainder of adiabatic eigenvectors
c
            write(8,1001) n1,ititl(3)
            write(8,1000) (eta(i),i=1,n1)
      write(8,1001) n1,ititl(4)
      write(8,1000) (adrho(i), i=1,n1)
      write(8,1001) n1,ititl(5)
      write(8,1000) (gam(i), i=1,n1)
      write(lyname,"(21h{d[h/r (adiabatic, p=,i2,4h, g=,i2,4h, l=,i2,2h)
     1$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(30hdh/r vs zone ... adiabatic, p=,i2,4h, g=,i2,4h, l=
     1,i1)") np,ng,ldeg
      call russ(zone,eta,n1,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(25h{d[rho/rho (adiabatic, p=,i2,4h, g=,i2,4h, l=,i2
     1,2h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(34hdrho/rho vs zone ... adiabatic, p=,i2,4h, g=,i2,4h
     1, l=,i2)") np,ng,ldeg
      call russ(zone,adrho,n1,0,0,lxname,lyname,ltop,ides)
      write(lyname,"(21hg[amma (adiabatic, p=,i2,4h, g=,i2,4h, l=,i2,2h)
     1$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(31hgamma vs zone ... adiabatic, p=,i2,4h, g=,i2,4h, l
     1=,i1)") np,ng,ldeg
      call russ(zone,gam,n1,0,0,lxname,lyname,ltop,ides)
      return
c
c
c        non-adiabatic plot dump
c
  20  continue
      call seper(dror,plotd,plotc,1,n)
      write(8,1001) n,ititl(6)
      write(8,1000) (plotd(i),i=1,n)
      write(lyname,"(24h{d[r/r (nonadiabatic, p=,i2,4h, g=,i2,4h, l=,i2,
     12h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(33hdr/r vs zone ... nonadiabatic, p=,i2,4h, g=,i2,4h,
     1 l=,i2)") np,ng,ldeg
      call russ(zone,plotd,n,0,0,lxname,lyname,ltop,ides)
      write(8,1001) n,ititl(7)
      write(8,1000) (plotc(i),i=1,n)
      if( iplot .eq. 1 ) return
      if( iplot .eq. 2 ) return
      m = 8
            call seper(ceta,plotd,plotc,2,n1)
            write(8,1001) n1,ititl(m)
            write(8,1000) (plotd(i),i=1,n1)
      write(lyname,"(24h{d[h/r (nonadiabatic, p=,i2,4h, g=,i2,4h, l=,i2,
     12h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(33hdh/r vs zone ... nonadiabatic, p=,i2,4h, g=,i2,4h,
     1 l=,i2)") np,ng,ldeg
      call russ(zone,plotd,n1,0,0,lxname,lyname,ltop,ides)
            write(8,1001) n1,ititl(m+1)
            write(8,1000) (plotc(i),i=1,n1)
            m = 10
      call seper(dgam, plotd,plotc,3,n1)
      write(8,1001) n1,ititl(m)
      write(8,1000) (plotd(i),i=1,n1)
      write(lyname,"(24hg[amma (nonadiabatic, p=,i2,4h, g=,i2,4h, l=,i2,
     12h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(34hgamma vs zone ... nonadiabatic, p=,i2,4h, g=,i2,4h
     1, l=,i2)") np,ng,ldeg
      call russ(zone,plotd,n1,0,0,lxname,lyname,ltop,ides)
      write(8,1001) n1,ititl(m+1)
      write(8,1000) (plotc(i),i=1,n1)
      m=12
      call seper(dent, plotd,plotc,4,n1)
      write(8,1001) n1,ititl(m)
      write(8,1000) (plotd(i),i=1,n1)
      write(lyname,"(25h{d[S/cv (nonadiabatic, p=,i2,4h, g=,i2,4h, l=,
     1i1,
     12h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(34hdS/cv vs zone ... nonadiabatic, p=,i2,4h, g=,i2,4h
     1, l=,i2)") np,ng,ldeg
      call russ(zone,plotd,n1,0,0,lxname,lyname,ltop,ides)
      write(8,1001) n1,ititl(m+1)
      write(8,1000) (plotc(i),i=1,n1)
      m=14
      call seper(drho,plotd,plotc,5,n1)
      write(8,1001) n1,ititl(m)
      write(8,1000) (plotd(i),i=1,n1)
      write(lyname,"(24h{dr/r[ (nonadiabatic, p=,i2,4h, g=,i2,4h, l=,i2,
     12h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(37hdrho/rho vs zone ... nonadiabatic, p=,i2,4h, g=,i2
     1,4h, l=,i2)") np,ng,ldeg
      call russ(zone,plotd,n1,0,0,lxname,lyname,ltop,ides)
      write(8,1001) n1,ititl(m+1)
      write(8,1000) (plotc(i),i=1,n1)
      return
c
c          work per zone
c
  30  continue
      write(8,1001) n1,ititl(18)
      write(8,1000) (work(i),i=1,n1)
      write(lyname,"(23hw[ork]/z[one (ergs) (p=,i2,4h, g=,i2,4h, l=,i2,2
     1h)$)") np,ng,ldeg
      lxname="z[one$"
      write(ides,"(33hwork vs zone ... nonadiabatic, p=,i2,4h, g=,i2,4h,
     1 l=,i2)") np,ng,ldeg
      call russ(zone,work,n1,0,0,lxname,lyname,ltop,ides)
      return
1000  format(1p6e13.5)
1001  format(i4,1x,10h(1p6e13.5),25x,a10)
      end
      subroutine seper(cyi,re,ai,j,n)
c
c         returns real and imag. part of cyi in re and ai.
c
      complex cyi
      dimension cyi(1),re(1),ai(1)
      do 10 i=1,n
            re(i) = real(cyi(i))
            ai(i) = aimag(cyi(i))
  10  continue
      return
      end
      subroutine ps
c     perform and print a period spacing analysis. rbk
      common /sumtab/ nrq(100),npn(100),ngn(100),pa(100),pna(100)
     1 ,pgr(100),wgr(100),tke(100),iroot,rl,ititl(3)
      common/phypar/ rp(1)
      common/blk1/   irad,nflag,ir,iplot,irp1,isum
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      if(pa(iroot).le.0.0) iroot=iroot-1
      if(iroot.lt.2) return
      l=ifix(rl+.1)
      write(6,100) l
100   format(1h1/////30x,22hPERIOD SPACINGS FOR l=,i2///
     1 9x6hRadial6x9hNumber of6x9hNumber of3x12h   Adiabatic9x6hPeriod/
     2 8x7hQuantum9x6hP-type9x6hG-type9x6hPeriod8x7hSpacing/
     3 9x6hNumber10x5hNodes10x5hNodes6x9h(minutes)6x9h(minutes)/)
      write(6,200) nrq(1),npn(1),ngn(1),pa(1)/60.,0.0
200   format(3i15,2f15.2)
      pp=pa(1)
      do 10 i=2,iroot
      if(pa(i).le.0.0) go to 10
      del=abs(pp-pa(i))
      pp=pa(i)
      write(6,200) nrq(i),npn(i),ngn(i),pa(i)/60.,del/60.
10    continue
      avedel=abs(pa(1)-pa(iroot))/float(iroot-1)
      pzero=avedel*sqrt(rl*(rl+1.))
      write(6,300) avedel/60.,pzero/60.
300   format(//36x24hAverage Period Spacing =f15.2//
     1 26x34h(Average Spacing)*(SQRT(L(L+1))) =f15.2)
      write(6,3000) (ititl(i),i=1,3),ir,rlumgv,totmas,
     1 teff,rphoto,rp(1),corlum
3000  format(///1x,3a8,i4,7h zones://
     1 1x,24hTotal Luminosity (erg/s)1x1pe11.4/
     2 1x,15hTotal Mass (gm)10x1pe11.4/
     3 1x,25hEffective Temperature (K)1pe11.4/
     4 1x,24hPhotospheric Radius (cm)1x1pe11.4/
     5 1x,16hCore Radius (cm)9x1pe11.4/
     6 1x,23hCore Luminosity (erg/s)2x1pe11.4)
      end
      subroutine st
c     print a summary table of important parameters. rbk
      common /sumtab/ nrq(100),npn(100),ngn(100),pa(100),pna(100)
     1 ,pgr(100),wgr(100),tke(100),iroot,rl,ititl(3)
      common/phypar/ rp(1)
      common/blk1/   irad,nflag,ir,iplot,irp1,isum
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      if(pa(iroot).le.0.0) iroot=iroot-1
      if(iroot.lt.1) return
      l=ifix(rl+.1)
      write(6,100) l
100   format(1h1/////30x20hSUMMARY TABLE FOR l=i1///
     1 64x4hWork/
     2 8x2(2x6hNumber)18x4hNon-2x9hPredicted3x8hIntegral6x5hTotal/
     3 2x6hRadial2(6x2hof)2(2x9hAdiabatic)2(5x6hGrowth)4x7hKinetic/
     4 8h Quantum2x6hP-type2x6hG-type2(5x6hPeriod)2(7x4hRate)5x6hEnergy/
     *2x6hNumber2(3x5hNodes)2(2x9h(seconds))2(2x9h(e/cycle))5x6h(ergs)/)
      do 10 i=1,iroot
      if(pa(i).le.0.0) go to 10
      write(6,200) nrq(i),npn(i),ngn(i),pa(i),pna(i),pgr(i),
     1 wgr(i),tke(i)
200   format(3i8,2f11.2,1p3e11.2)
10    continue
      write(6,3000) (ititl(i),i=1,3),ir,rlumgv,totmas,
     1 teff,rphoto,rp(1),corlum
3000  format(///1x,3a8,i4,7h zones://
     1 1x,24hTotal Luminosity (erg/s)1x1pe11.4/
     2 1x,15hTotal Mass (gm)10x1pe11.4/
     3 1x,25hEffective Temperature (K)1pe11.4/
     4 1x,24hPhotospheric Radius (cm)1x1pe11.4/
     5 1x,16hCore Radius (cm)9x1pe11.4/
     6 1x,23hCore Luminosity (erg/s)2x1pe11.4)
      end
      subroutine russ(xx,yy,n,lx,ly,lxname,lyname,ititl,index)
      dimension xx(1),yy(1)
      dimension x(2001),y(2001),ititl(1)
      character *(*),lxname,lyname,index
      do 1 i=1,n
      x(i)=xx(i)
      y(i)=yy(i)
 1    continue
      pagex=14.
      pagey=11.
      xaxis=11.
      yaxis=9.
      xphys=2.25
      yphys=1.
      call mx1alf("standard","]")
      call mx2alf("l/cstd","[")
      call mx3alf("math","#")
      call mx4alf("instru","'")
      call mx5alf("l/cgreek","{")
      call nobrdr
      call page(pagex,pagey)
      call physor(xphys,yphys)
      call area2d(xaxis,yaxis)
      call height(.28)
      alen=xmess(ititl,24)
      call messag(ititl,24,.5*(xaxis-alen),yaxis*1.025)
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      do 10 i=2,n
      xmin=amin1(xmin,x(i))
      xmax=amax1(xmax,x(i))
      ymin=amin1(ymin,y(i))
10    ymax=amax1(ymax,y(i))
      if(lxname(1:1).eq."z") call xintax
      call xrevtk
      call yrevtk
      call yaxang
      call graf(xmin,"scale",xmax,ymin,ymax-ymin,ymax)
      if(lx.eq.1) go to 11
      call xgraxs(xmin,"scale",xmax,xaxis,lxname,100,.0,.0)
      go to 12
11    continue
      call algplt(xmin,xmax,xaxis,xorg,xcyc)
      call xlgaxs(xorg,xcyc,xaxis,lxname,100,.0,.0)
12    continue
      if(ly.eq.1) go to 13
      ym=amax1(abs(ymax),abs(ymin))
      iexp=ifix(alog10(ym))
      ymaxt=ymax*10.**(-iexp)
      ymint=ymin*10.**(-iexp)
      call axsplt(ymint,ymaxt,yaxis,yorg,ystep,axisy)
      ymaxt=(yorg+ystep*axisy)*10.**iexp
      yorg=yorg*10.**iexp
      ystep=ystep*10.**iexp
      call ygraxs(yorg,ystep,ymaxt,yaxis,lyname,100,.0,.0)
      go to 14
13    continue
      call algplt(ymin,ymax,yaxis,yorg,ycyc)
      call ylgaxs(yorg,ycyc,yaxis,lyname,100,.0,.0)
14    continue
      call xnonum
      call ynonum
      if(lx.eq.0) call xgraxs(xmin,"scale",xmax,xaxis," ",-1,.0,yaxis)
      if(lx.eq.1) call xlgaxs(xorg,xcyc,xaxis," ",-1,.0,yaxis)
      if(ly.eq.0) call ygraxs(yorg,ystep,ymaxt,yaxis," ",-1,xaxis,.0)
      if(ly.eq.1) call ylgaxs(yorg,ycyc,yaxis," ",-1,xaxis,.0)
      call reset("xintax")
      call reset("xnonum")
      call reset("ynonum")
      call height(.14)
      call curve(x,y,n,0)
      call endpl(0)
      iplot=iplot+1
      if(iplot.eq.1) write(15,16)
      write(15,20) iplot,index
20    format(i5,3x,a)
16    format(" For your convenience, a PLOT file has been generated"/
     1" which contains many informative graphs of data and results."//
     2" If you imagine that you do not have time to sit at a graphics"/
     3" terminal and look at selected pictures (but you can wait"/
     4" awhile before you have to see the pictures), you can type:"/
     5" PESP . . . which sends them all to the electrostatic printer"/
     6" PFILM . . . which sends them all to microfiche"/
     7" PFILM FORMAT=35 . . . which sends them all to 35mm film."//
     8" On the other hand you can begin to look at your picture"/
     9" survey right away by typing PSCAN, wait for a prompt, and"/
     1" then type Dm where m is the number of the plot you want to"/
     2" display (see table of contents below).  While in PSCAN you"/
     3" can generate a smaller plot file (NPLOT) to keep or send to"/
     4" PESP (or PFILM) by first typing CREATE=NPLOT and then typing"/
     5" Wm for each plot you want to go to NPLOT (or Wm p i which"/
     6" writes plots m through p with increment i)."////
     7"        PLOT FILE CONTENTS"//)
      return
      end
