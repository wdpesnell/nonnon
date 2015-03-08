      program nonad
      implicit real*8(a-h,o-z)
      character*80 ititl
      logical bigout
c
c          the lagrangian non-radial adiabatic eigenvalue analysis
c       for stellar pulsations. the notation is lifted from the model
c       code and radial stability analysis. written summer 1986.
c          this version of the code treats the "realistic" stellar
c       models with the r**l transformation in the eigenvectors.
c
c                                                12/20/83 wd pesnell
c
c          Adapted to do nonadiabatic pulsations, with central 
c       ball included, July 9-20, 1989, Los Alamos Nat'l. Lab.
c
c                                                 7/30/89 WD Pesnell
c
      parameter ( nmax=600 )
      common/phypar/ rp(nmax),tp(nmax),vp(nmax),cv(nmax),dkdr(nmax),
     $   dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),rho(nmax),rzone(nmax)
      common/blk8/   g,ac3,acrad,pi,twopi,forpi,pi8,pi43
      common/coretc/ pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,g10,g3m10,
     $               cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/coolum/ onemq0,onemq1,rlums,rlumc,noburn
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/stnmod/ nmode
c
c          read in the input variables from file 10.
c
c      open (unit=1,file='sumout',status='unknown',acess='append')
c      open (unit=5,file='sys$input',status='unknown')
c      open (unit=6,file='sys$output',status='unknown')
      open (unit=10,file='lnrad',status='old')
cPC      open (unit=10,file=' ',status='old')
c
      open (unit=1,file='sumout.adn',status='unknown')
      open (unit=11,file='adnout',status='unknown')
c
c      call filerep
c
      read(10,1001,end=100) ititl
c
c          read in the input variables from file 10.
c
      read(10,*) npts,rlumgv,totmas,teff,rphoto,corlum
      read(10,*) irad,noburn,onemq0,onemq1,rlums,rlumc
      np1 = npts+1
      read(10,1000) (rp(i),i=1,np1)
      read(10,1000) (tp(i),i=1,npts)
      read(10,1000) (vp(i),i=1,npts)
      read(10,1000) (cv(i),i=1,npts)
      read(10,1000) (dkdr(i),i=1,npts)
      read(10,1000) (dkdt(i),i=1,npts)
      read(10,1000) (dm1(i),i=1,npts)
      read(10,1000) (akap(i),i=1,npts)
      read(10,1000) (dm2(i),i=1,np1)
      read(10,1000) (rm(i),i=1,np1)
      read(10,1000) (p(i),i=1,npts)
      read(10,1000) (g1(i),i=1,npts)
      read(10,1000) (g3m1(i),i=1,npts)
      read(10,1000) (frft(i),i=1,np1)
      read(10,1000) (sorce(i),i=1,npts)
      read(10,1000) (dtsorc(i),i=1,npts)
      read(10,1000) (dvsorc(i),i=1,npts)
      read(10,1000) (bv(i),i=1,np1)
      read(10,1000,end=900) pc,rhoc,tc,cormas,rl0,chit0,chr0,q0,
     $   g10,g3m10,cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
      close (unit=10)
      totmas = rm(npts+1)
c
c	   Does user want to write the initial model?
c
      bigout = .false.
      write(6,1100)
      read(5,1101,end=110) iyorn
      if( iyorn .eq. 'y' .or. iyorn .eq. 'Y' ) bigout = .true.
 110  continue
c
c          output the initial model data
c
      write(11,2000) ititl,npts
      write(11,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(11,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
      write(11,2004)
      write(11,2005) npts,rp(1),tc,1.e0/rhoc,pc,g10,g3m10,opac0,
     $	    dkdr0,dkdt0,cv0,frft(1),rm(1)
      if( bigout ) then
	 write(11,2005) (i,rp(i+1),tp(i),vp(i),p(i),g1(i),g3m1(i),
     $	    akap(i),dkdr(i),dkdt(i),cv(i),frft(i+1),rm(i+1),i=1,npts)
c
c          check for nuclear burning and print out those
c       zones participating.
c
	 if( tp(1) .ge. 3.d6 ) then
	    do 10 i=1,npts
	       imax = np1-i
	       if( sorce(imax) .ne. zero ) goto 15
  10	    continue
	    write(11,2100)
	    goto 20
  15	    continue
	    do 16 i=1,imax
	       imin = i
	       if( sorce(i) .ne. zero ) goto 17
  16	    continue
	    write(11,2100)
	    goto 20
  17	    continue
	    write(11,2101)
	    write(11,2102) (i,sorce(i),dtsorc(i),dvsorc(i),i=imin,imax)
  20	    continue
	 endif
      endif
      call pltint(npts,ititl)
      write(1,2000) ititl,npts
      write(1,2001) rlumgv,totmas,teff,rphoto,rp(1),corlum
      write(1,2002) irad,noburn,onemq0,onemq1,rlums,rlumc
c
c          The pulsation control loop.
c
c          iout:= 0 for little output (no eigenvectors)
c                 1 for convergence information
c                 2 for eigenvectors, matrices and conver. info.
c
c          lin:= laditudinal quantum number.
c
c          nomega:= number of points in frequency scan. (lin>0)
c                   number of nodes -1 in first mode (lin=0)
c
c          ihmax:= not used (lin>0), enter a 0.
c                   number of nodes - 1 in last mode (lin=0)
c
c          permin:= starting period for scan (sec).
c
c          permax:= ending period for scan (sec).
c
c          (it is not necessary for permin to be less than permax.)
c
c          Example input line: 0 0 1 3 0 0
c
c          This would calculate the first three radial modes (fundamental
c       and first and second overtone, with little output.
c
c          Example input line: 3 2 10 0 4000 4500
c
c          This would search 10 periods between 4000 and 4500 seconds
c       for l=2 modes, producing a great deal of output.
c
c          In every case, a plot file is produced with the eigenvectors,
c       weight functions, and global quantities such as pressure.
c
 200  continue
         write(6,1003)
         read(5,*,end=900) iout,lin,nomega,ihmax,permin,permax
         if( lin .lt. 0 ) then
            goto 900
         elseif( lin .eq. 0 ) then
	    write(6,3000)
         else
            omlow = (twopi/permax)**2
            omhigh = (twopi/permin)**2
            write(11,2003) iout,lin,nomega,omlow,omhigh
            lval = lin
            nmode = 0
            call lnanon(npts,iout,lval,nomega,omlow,omhigh)
         endif
      goto 200
c
c          if end of file, try again.
c
 100  continue
      write(1,1002)
      stop
 900  continue
      stop
c
 1000 format(1p,4e20.13)
 1001 format(a80)
 1002 format(1x,'End-of-File on title, try again.')
 1003 format(1x,27henter pulsation parameters:,/,
     $   1x,38h io, l, nom, ihmax, permin, and permx.)
 1100 format(1x,'Write initial model data? (y/n) ',$)
 1101 format(a)
c
 2000 format(50h1 initial model in hydrostatic and thermal balance,//,
     $      a80,1x,i5,6h zones,/)
 2001 format(4x,18h total luminosity=,1pe11.4,12h total mass=,e11.4,
     $ 17h effective temp.=,e11.4,/,6x,21h photospheric radius=,e11.4,
     $ 13h core radius=,e11.4,11h core lum.=,e11.4)
 2002 format(6h irad=,i2,8h noburn=,i2,8h onemq0=,1pe10.3,
     $ 8h onemq1=,e10.3,7h rlums=,e10.3,7h rlumc=,e10.3)
 2003 format(6h iout=,i2,13h l-value used,i4,
     $ 25h number of points in scan,i4,15h starting freq.,1pe10.3,
     $ 13h ending freq.,e10.3)
 2004 format(3x,1hi,3x,6hradius,5x,4htemp,5x,6hsp.vol,3x,8hpressure,
     $ 5x,2hg1,7x,4hg3m1,5x,7hopacity,3x,6hdlkdlr,4x,6hdlkdlt,6x,
     $ 2hcv,7x,4hfrft,4x,8hint.mass)
 2005 format(1x,i4,1p,12e10.3)
 2100 format(/,30h model has no nuclear burning.)
 2101 format(/,29h following zones are ignited.)
 2102 format(3(1x,i4,1p,3e10.3))
 3000 format(1x,'No radial analysis available at this time.')
      end
      block data
      implicit real*8(a-h,o-z)
c
c       fundamental constants are from novotny, introduction to
c    stellar atmospheres and interiors, 1973, appendix ii.
c                                                 2/16/83  wdp
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
      common/thermo/ r,a,bk,avagd,ad3
      data r,a,bk,avagd,ad3 / 8.31434d7,7.56471d-15,8.6170837d-5,
     $     6.02217d23,2.52157d-15 /
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726d-8,7.5595d-5,
     $  5.66961d-5,3.1415926536d0,6.2831853072d0,12.5663706144d0,
     $   25.1327412288d0,4.1887902048d0 /
      data zero,one,two,thre,for,ten,ahf,qrt  /
     $      0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,10.0d0,
     $      0.5d0,0.25d0      /
      end
      subroutine pltint(n,ititl)
      implicit real*8(a-h,o-z)
      character*80 ititl
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=600 )
      common/phypar/ r(nmax),t(nmax),v(nmax),cv(nmax),dkdr(nmax),
     $	 dkdt(nmax),dm1(nmax),akap(nmax),dm2(nmax),rm(nmax),bv(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk4/   p(nmax),g1(nmax),g3m1(nmax),omq(nmax),tlog(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      dimension plog(nmax),xrad(nmax)
c
      open (unit=12,file='nonplt',status='unknown')
      write(12,1000) ititl
      np1 = n + 1
      call pltdmp(r,nmax,n+1,'r   ')
      sfact = one - 5.d-14
      do 10 i=1,n
         omq(i) =-dlog10(one - sfact*rm(i)/totmas)
         tlog(i) = dlog10(t(i))
         plog(i) = dlog10(p(i))
         xrad(i) = r(i)/r(np1)
  10  continue
      xrad(np1) = one
      omq(n+1) =-dlog10(one-sfact)
      call pltdmp(tlog,nmax,n,'t   ')
      call pltdmp(plog,nmax,n,'p   ')
      call pltdmp(omq,nmax,n+1,'1-q ')
      call pltdmp(xrad,nmax,n+1,'x   ')
      call pltdmp(g1,nmax,n,'g1   ')
      call pltdmp(g3m1,nmax,n,'g3m1')
      do 20 i=1,n
         xrad(i) = dfloat(i)
  20  continue
      call pltdmp(xrad,nmax,n,'zone')
      do 30 i=1,n
         xrad(i) = 0.
	 if( bv(i) .ne. 0. ) xrad(i) = dlog10(dabs(bv(i)))
  30  continue
      call pltdmp(xrad,nmax,n,'bv  ')
      return
 1000 format(a80)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
      character*4 ititl
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(12,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p,6e12.4) )
      end
