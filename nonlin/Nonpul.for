      PROGRAM NONPUL
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 iyorn
      integer*2 np1
      LOGICAL IPRT, linmod
      complex*16 dr, dl, dro, dt, dp, comega, eiomt, ci
      parameter (rlsol=3.826d33, rmsol=1.991d33, rsol=6.959d10 )
      parameter ( ci=(0.d0, 1.d0) )
c
c          The Lagrangian radial nonadiabatic nonlinear analysis
c       for stellar pulsations. The notation is lifted from the MODEL
c       code and radial stability analysis. Begun March 1988, NMSU.
c
c                                                3/17/88 WD Pesnell
c       12/03/91 Modifications begun at NASA/GSFC to move code to
c                time dependent convection.
c       02/28/92 Modified extensively during move to Building 22. There
c                was nothing else to do. Included the development of the
c                STRUCTURE/RECORD for plotting and linear eigenvector
c                stuff.
c       03/05/92 Removed references to file 9 (the summary file). Changed
c                the error limits on \Delta T/T to prevent an inversion.
c                Basically replicated the radius limitations. Modified 
c                the STRUCTURE PLOTVEC to prevent overwriting of other
c                variables. Gave time-step to time-step oscillations in
c                all variables.
c
      parameter ( nmax=256, nmaxtd=5000, nmodmx=3 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),p(nmax),up(nmax),
     $               ep(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/scrtch/ dr(nmax,nmodmx), dro(nmax,nmodmx), dp(nmax,nmodmx),
     $               dt(nmax,nmodmx), dl(nmax,nmodmx)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/blk8/   grav,pi,pi2,pi4,pi43
      common/inhom/  x,y,z,rmuc,pargrm
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/stlums/ rlumn1(nmax),rlumn(nmax),rlcon(nmax)
      parameter ( nzone=5 )
      type plotvec
         real*8 t_vals(nmaxtd,nzone)
         real*8 int_val(nzone)
      end type
      type(plotvec):: rt, tmax, umax, lmax
      dimension comega(nmodmx), sttime(nmaxtd)
C
C          FILES IN THIS PROGRAM.
C  UNIT   ROUTINE    NAME                PURPOSE
C    5      MAIN    SYS$INPUT           TERMINAL INPUT FILE
C    6      MAIN    SYS$OUTPUT          TERMINAL OUTPUT FILE
C
C   10      MAIN    LNRAD               CONTAINS THE MODEL QUANTITIES
C                                          FOR THE STABILITY ANALYSES.
C   11      MAIN    PULOUT              MAIN OUTPUT FILE
C
C   12     PLTINT   PULPLT              HAS VECTORS OF EIGENVECTORS
C                                          QUANTITIES FOR PLTCMD PROGRAM
C
cVAX      OPEN (UNIT=11,FILE='PULOUT',STATUS='NEW')
      open (unit=11, file='pulout.dat',status='unknown')
c
c          Get the physical parameters of the initial model (from MODEL)
c       and store in the common blocks.
c
      call readin
c
c          Read the linear, nonadiabatic eigenvectors from NLDUMP.
c
      IOUT = 0
      call linread( npts, nmodes, comega, ampmx)
      omega = dreal(comega(1))
      PERIOD = PI2/OMEGA
      PDAYS  = period/86400.D0
      if( npts .le. 0 ) stop
c
      np1 = npts + 1
      WRITE(6,1000) ampmx
      READ(5,*,END=900) NTIME,NPRT,AMP,DTFAC, phi_0
      IF( NTIME .LE. 0 ) GOTO 900
      IF( NTIME .GT. NMAXTD ) NTIME = NMAXTD
      IF( AMP .le. 0.d0 ) AMP = 10.D9/(OMEGA*RN1(NP1))
c
      linmod = .false.
      write(6,1010)
      read(5,1011,end=910) iyorn
      if( iyorn .eq. 'N' .or. iyorn .eq. 'n' ) then
         linmod = .true.
      else
         linmod = .false.
      endif
 910  continue
c
c          Set up a first guess to variables
c
      do 20 ii=1,nzone
         i = (ii*npts)/nzone
         rt%int_val(ii) = RN1(i)
         umax%int_val(ii) = UN1(i)*1.d-5
         tmax%int_val(ii) = TN1(i)
         lmax%int_val(ii) = RLUMN1(i)
  20  continue
c
      delta_0 = phi_0*pi/180.d0
      eiomt = cdexp( ci*delta_0 )
      UP(1) = ZERO
      DO 90 I=1,NPTS
         CALL STATE(I,1)
         CALL RADLUM(I,1)
         CALL CONLUM(I, 1, RLUMC)
         UN1(I+1)=-RP(I+1)*AMP*dimag(comega(1)*DR(I+1,1)*eiomt)
         CS(I) = DSQRT(G1(I)*PN1(I)*VN1(I))
  90  CONTINUE
      IOUT = 3
      DTMAX = PERIOD/1.D2
      ELTIME = ZERO
      DO 100 ITIME=1,NTIME
         IPRT = .FALSE.
         IF( MOD(ITIME,NPRT) .EQ. 0 ) IPRT = .TRUE.
c
         DTMIN = QTIME(RN1,CS,NMAX,NPTS,IPRT)*DTFAC
         IF( DTMIN .LE. ZERO ) DTMIN =-DTMIN
         DELT = DMIN1(DTMAX,DTMIN)
         ELTIME = ELTIME + DELT
         IF( linmod ) THEN
            eiomt = cdexp( ci*(comega(1)*eltime+delta_0) )
            if( itime .eq. 1 ) then
C
C          TRANSFER THE INITIAL VALUES INTO TO STORAGE.
C
               RP(1) = RN1(1)
               UP(1) = 0.d0
               RLUMN(1) = RLUMN1(1)
c
               DO 120 I=1,NPTS
                  RP(I+1) = RN1(I+1)
                  UP(I+1) = UN1(I+1)
                  VP(I) = VN1(I)
                  TP(I) = TN1(I)
                  P(I) = PN1(I)
                  RLUMN(I+1) = RLUMN1(I+1)
 120          continue
            endif
C
C          UPDATE THE TIME STEP VALUE WITH THE LINEAR EIGENVECTORS
C
            RP(1) = RN1(1)
            RLUMN(1) = RLUMN1(1)
            do 130 i=1, npts
               RN1(I+1) = RP(I+1)*(ONE + AMP*dreal(DR(I+1,1)*eiomt))
               UN1(I+1)=-RP(I+1)*AMP*dimag(comega(1)*DR(I+1,1)*eiomt)
               VN1(I) = VP(I)*(ONE - AMP*dreal(DRO(I,1)*eiomt))
               TN1(I) = TP(I)*(ONE + AMP*dreal(DT(I,1)*eiomt))
               PN1(I) = P(I)*(ONE + AMP*dreal(DP(I,1)*eiomt))
              RLUMN1(I+1)=RLUMN(I+1)*(ONE+AMP*dreal(DL(I+1,1)*eiomt))
 130        continue
         else
C
C          TRANSFER THE LAST TIME STEP TO STORAGE.
C
            RP(1) = RN1(1)
            RLUMN(1) = RLUMN1(1)
c
            DO 110 I=1,NPTS
               RP(I+1) = RN1(I+1)
               UP(I+1) = UN1(I+1)
               VP(I) = VN1(I)
               TP(I) = TN1(I)
               P(I) = PN1(I)
               EP(I) = EN(I)
               RLUMN(I+1) = RLUMN1(I+1)
               if( ITIME .EQ. 1 ) THEN
                  eiomt = cdexp( ci*(comega(1)*eltime+delta_0) )
                  RN1(I+1) = RN1(I+1)*(ONE+AMP*dreal(DR(I+1,1)*eiomt))
                  UN1(I+1)=-RP(I+1)*AMP*dimag(comega(1)*DR(I+1,1)*eiomt)
                  VN1(I) = VN1(I)*(ONE - AMP*dreal(DRO(I,1)*eiomt))
                  TN1(I) = TN1(I)*(ONE + AMP*dreal(DT(I,1)*eiomt))
                  PN1(I) = PN1(I)*(ONE + AMP*dreal(DP(I,1)*eiomt))
                  CALL STATE(I,1)
                  CALL RADLUM(I,1)
                  CALL CONLUM(I, 1, RLUMC)
                RLUMN1(I+1)=RLUMN1(I+1)*(ONE+AMP*dreal(DL(I+1,1)*eiomt))
               ELSE
                  RN1(I+1) = RN1(I+1) + UN1(I+1)*DELT*QRT
C                  VN1(I) = PI43*(RN1(I+1)**3-RN1(I)**3)/DM1(I)
C                  CALL STATE(I,1)
C                  CALL RADLUM(I,1)
C                  CALL CONLUM(I, 1, RLUMC)
C                  RLUMN1(I+1) = RLUMN1(I+1)*(ONE + AMP*DL(I+1))
               ENDIF
 110        CONTINUE
            CALL HYDRO(IOUT,DELT,ITIME,IERRHY)
            IF( IERRHY .LE. 0 ) GOTO 200
         endif
c
         PHASE = ELTIME/PERIOD
         STTIME(ITIME) = PHASE
         do 300 ii=1,nzone
            i = (ii*npts)/nzone
            rt%t_vals(ITIME,ii) = RN1(i)/rt%int_val(ii) - 1.d0
            umax%t_vals(ITIME,ii) = UN1(i)*1.d-5
            tmax%t_vals(ITIME,ii) = TN1(i)/tmax%int_val(ii) - 1.d0
            lmax%t_vals(ITIME,ii) = RLUMN1(i)/lmax%int_val(ii)-1.d0
 300     continue
         IF( IPRT ) THEN
            WRITE(11,2006) (I,RN1(I+1),UN1(I+1),TN1(I),VN1(I),
     $       PN1(I),I=NPTS-5,NPTS)
            WRITE(6,2006) ITIME,PHASE,RN1(NPTS),RN1(NP1),UN1(NPTS),
     $       UN1(NP1)
         ENDIF
 100  CONTINUE
 200  CONTINUE
      call pltdmp(sttime,nmaxtd,itime-1,'time')
      do 400 i=1, nzone
         call pltdmp(rt%t_vals(1:itime,i),nmaxtd,itime-1,'rt  ')
         call pltdmp(umax%t_vals(1:itime,i),nmaxtd,itime-1,'umax')
         call pltdmp(tmax%t_vals(1:itime,i),nmaxtd,itime-1,'tmax')
         call pltdmp(lmax%t_vals(1:itime,i),nmaxtd,itime-1,'rl  ')
 400  CONTINUE
 900  CONTINUE
      STOP
C
 1000 FORMAT(1X,'Enter ntime, nprt, amplitude (<',1pe8.1,
     $          '), dtfac, and initial phase: ')
 1010 format(1x,'Do you wish to run a nonlinear model? [Y]: ',$)
c 1010 format(1x,
c     $   '[L]inear, [N]on-adiabatic, or [A]diabatic model? [N]: ',$)
 1011 format(a)
 2006 FORMAT(1X,I4,1P,5E15.8)
      END
      BLOCK DATA
      IMPLICIT REAL*8(A-H,O-Z)
C
C       FUNDAMENTAL CONSTANTS ARE FROM NOVOTNY, INTRODUCTION TO
C    STELLAR ATMOSPHERES AND INTERIORS, 1973, APPENDIX II.
C                                                 2/16/83  WDP
C
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      COMMON/BLK8/   GRAV,PI,PI2,PI4,PI43
      COMMON/THERMO/ R,A,BK,AVAGD,AD3
      DATA R,A,BK,AVAGD,AD3 / 8.31434D7,7.56471D-15,8.6170837D-5,
     $     6.02217D23,2.52157D-15 /
      DATA GRAV,PI,PI2,PI4,PI43/6.6726D-8,
     $   3.1415926536D0,6.2831853072D0,12.5663706144D0,
     $   4.1887902048D0 /
      DATA ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT/
     $      0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,10.0D0,
     $      0.5D0,0.25D0/
      END
      subroutine readin
      implicit real*8(a-h, o-z)
      character*80 ititl
      integer*2 i, np1
c
c          Read the initial model information from LNRAD
c
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/blk8/   grav,pi,pi2,pi4,pi43
      common/inhom/  x,y,z,rmuc,pargrm
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/stlums/ rlumn1(nmax),rlumn(nmax), rlcon(nmax)
C
C          READ IN THE INPUT VARIABLES FROM FILE 10.
C
c      OPEN (UNIT=10,FILE='LNRAD',STATUS='OLD',READONLY,err=900)
      OPEN (UNIT=10,FILE='evnondmp.dat',STATUS='OLD',action='READ',
     $      err=900)
      READ(10,1001) ITITL
      READ(10,*) NPTS,RLUMGV,TOTMAS,TEFF,RPHOTO,CORLUM
      READ(10,*) IRAD,NOBURN,ONEMQ0,ONEMQ1,RLUMS,RLUMC
      NP1 = NPTS+1
      READ(10,1000) (RN1(I),I=1,NP1)
      READ(10,1000) (TN1(I),I=1,NPTS)
      READ(10,1000) (VN1(I),I=1,NPTS)
      READ(10,1000) (CV(I),I=1,NPTS)
      READ(10,1000) (DKDR(I),I=1,NPTS)
      READ(10,1000) (DKDT(I),I=1,NPTS)
      READ(10,1000) (DM1(I),I=1,NPTS)
      READ(10,1000) (AKAP(I),I=1,NPTS)
      READ(10,1000) (DM2(I),I=1,NP1)
      READ(10,1000) (RM(I),I=1,NP1)
      READ(10,1000) (PN1(I),I=1,NPTS)
      READ(10,1000) (G1(I),I=1,NPTS)
      READ(10,1000) (G3M1(I),I=1,NPTS)
      READ(10,1000) (FRFT(I),I=1,NP1)
      READ(10,1000) (SORCE(I),I=1,NPTS)
      READ(10,1000) (DTSORC(I),I=1,NPTS)
      READ(10,1000) (DVSORC(I),I=1,NPTS)
      CLOSE (UNIT=10)
c
      iout = 0
      ihmin = 0
      ihmax = 0
      WRITE(11,2000) ITITL,NPTS
      WRITE(11,2001) RLUMGV,TOTMAS,TEFF,RPHOTO,RN1(1),CORLUM
      WRITE(11,2002) IRAD,NOBURN,ONEMQ0,ONEMQ1,RLUMS,RLUMC
      WRITE(11,2003) IOUT,IHMIN,IHMAX
      WRITE(11,2004)
      WRITE(11,2005) (I,RN1(I+1),TN1(I),VN1(I),PN1(I),G1(I),G3M1(I),
     $   AKAP(I),DKDR(I),DKDT(I),CV(I),FRFT(I+1),RM(I+1),I=1,NPTS)
c
      NFLAG = 3
      X = 0.7D0
      Y = 0.28D0
      Z = 0.02D0
      PARGRM = X/1.00797D0 + Y/4.0026D0 + Z/21.02447D0
      RMUC = ONE/PARGRM
      CALL PLTINT(NPTS)
c
c          Find the initial luminosity profile
c
      rltot = rlumgv
      rlumn1(np1) = rlumgv
      rlcon(np1) = zero
      do 10 i=npts, 1, -1
         rltot = rltot - sorce(i)*dm1(i)
         rlumn1(i) = frft(i)*rltot
         rlcon(i) = (one-frft(i))*rltot
  10  continue
      return
c
 900  continue
      write(6,9000)
      stop 901
C
 1000 FORMAT(1P,4E20.13)
 1001 FORMAT(A)
 1010 FORMAT(1P,4E20.13)
C
 2000 FORMAT(1H1,/,1X,19HOUTPUT FROM NONPUL:,/,
     $      50H  INITIAL MODEL IN HYDROSTATIC AND THERMAL BALANCE,//,
     $      2X,A,1X,I5,6H ZONES,/)
 2001 FORMAT(4X,18H TOTAL LUMINOSITY=,1PE11.4,12H TOTAL MASS=,E11.4,
     $ 17H EFFECTIVE TEMP.=,E11.4,/,6X,21H PHOTOSPHERIC RADIUS=,E11.4,
     $ 13H CORE RADIUS=,E11.4,11H CORE LUM.=,E11.4)
 2002 FORMAT(6H IRAD=,I2,8H NOBURN=,I2,8H ONEMQ0=,1PE10.3,
     $ 8H ONEMQ1=,E10.3,7H RLUMS=,E10.3,7H RLUMC=,E10.3)
 2003 FORMAT(10X,6H IOUT=,I2,7H IHMIN=,I3,7H IHMAX=,I3)
 2004 FORMAT(3X,1HI,3X,6HRADIUS,5X,4HTEMP,5X,6HSP.VOL,3X,8HPRESSURE,
     $ 5X,2HG1,7X,4HG3M1,5X,7HOPACITY,3X,6HDLKDLR,4X,6HDLKDLT,6X,
     $ 2HCV,7X,4HFRFT,4X,8HINT.MASS)
 2005 FORMAT(1X,I4,1P,12E10.3)
c
 9000 format(/1x,'Model file from initial model builder was not found.',
     $ /,1x,'Please make sure file is available before running again.'/)
      end
      SUBROUTINE LINREAD( npts, nmodes, comega, ampmx)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 dr, dl, dt, dp, dro, comega
      integer*2 imode, i
C
C          Read the linear eigenvectors calculated by CMPLNA.
C
      parameter ( nmax=256, nmodmx=3 )
      common/scrtch/ dr(nmax,nmodmx), dro(nmax,nmodmx), dp(nmax,nmodmx),
     $               dt(nmax,nmodmx), dl(nmax,nmodmx)
      dimension reg(nmodmx), comega(nmodmx), dlabs(nmax)
c
c      open(unit=40,file='nldump',status='old',err=900,readonly)
      open(unit=40,file='nldump',status='old',err=900,action='read')
      read(40,1000) nmodes
      if( nmodes .gt. nmodmx ) nmodes = nmodmx
      do 10 imode=1, nmodes
         read(40,1000) imode1
         read(40,1010) comega(imode)
         write(6,2000) imode, comega(imode)
         read(40,1010) (dr(i,imode),i=1,npts+1)
         read(40,1010) (dl(i,imode),i=1,npts+1)
         read(40,1010) (dro(i,imode),i=1,npts)
         read(40,1010) (dt(i,imode),i=1,npts)
         read(40,1010) (dp(i,imode),i=1,npts)
  10  continue
      close(unit=40)
c
c          Find the maximum amplitude in the DL/L eigenvector. This is
c       usually the largest and limits the amplitude the most.
c
      do 20 imode=1, nmodes
         do 25 i=1,npts+1
            dlabs(i) = cdabs( dl(i,imode) )
  25     continue
         reg(imode) = qmax8( dlabs, npts+1 )
  20  continue
      ampmx = qmax8( reg, nmodes )
      ampmx = 0.5d0/ampmx
      return
c
 900  continue
      write(6,9000)
      npts = -1
      return
c
 1000 format(i5)
 1010 format(1p,4e20.12)
 2000 format(1x,'Linear mode ',i2,' has frequency ',1p,2e11.3)
 9000 format(1x,'File with linear eigenvectors not found ... ')
      end
      SUBROUTINE PLTINT(N)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (rlsol=3.826d33, rmsol=1.991d33, rsol=6.9599d10)
C
C          INITIALIZE THE PLOT FILE AND WRITE THE TITLE LINE.
C
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/blk37/  x(nmax),tlog(nmax),rholog(nmax),onemq(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
C
cVAX      OPEN (UNIT=12,FILE='PULPLT',STATUS='NEW')
      open (unit=12, file='pulplt.dat', status='unknown')
C
      RMS = TOTMAS/RMSOL
      RL = DLOG10(RLUMGV/RLSOL)
      WRITE(12,1000) RMS,RL,TEFF
      QM = AHF*DM1(N)
      DO 10 I=1,N
         X(I+1) = RN1(I+1)/RN1(N+1)
         RHOLOG(I) =-DLOG10(VN1(I))
         ONEMQ(N+2-I) =-DLOG10(QM/TOTMAS)
         QM = QM + DM1(N+1-I)
  10  CONTINUE
      ONEMQ(1) =-DLOG10(QM/TOTMAS)
      X(1) = RN1(1)/RN1(N+1)
c
      call pltdmp(x,nmax,n+1, 'x   ')
      call pltdmp(tn1,nmax,n, 't   ')
      call pltdmp(pn1,nmax,n, 'p   ')
      call pltdmp(rholog,nmax,n,'rho ')
      call pltdmp(onemq,nmax,n+1,'1-q ')
      call pltdmp(g1,nmax,n,'g1  ')
c
      RETURN
 1000 FORMAT(1X,F7.2,'M!dsun!n L=',F7.3,'L!dsun!n T!deff!n=',f8.1)
      END
      SUBROUTINE PLTDMP(VEC,NMAX,N,ITITL)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*4 ITITL
C
C          WRITE THE VECTOR VEC TO THE PLOT FILE (TAPE12) WITH APPENDED
C       TITLE ITITL.
C
      DIMENSION VEC(NMAX)
      WRITE(12,1000) N,ITITL,(VEC(I),I=1,N)
      RETURN
 1000 FORMAT(I4,10X,A4,/,(1P,6E12.4) )
      END
      SUBROUTINE HYDRO(IOUT,DTIN,ITIME,IERRHY)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter ( onth=1.d0/3.d0, toths=2.d0/3.d0 )
C      PARAMETER ( accur=1.D-6 ) ! USED 12/2/91 2000 STEPS 22.75 HRS DE614
      PARAMETER ( accur=1.D-3 ) ! USED 12/4/91 1000 STEPS 4.0 HRS DE614
      PARAMETER ( ITRYMX=40 )
c
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),p(nmax),up(nmax),
     $               ep(nmax)
      common/stlums/ rlumn1(nmax),rlumn(nmax), rlcon(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/blk8/   g,pi,pi2,pi4,pi43
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/thermo/ r,a,bk,avagd,ad3
c      dimension x(nmax),tlog(nmax),plog(nmax),rholog(nmax)
      dimension gn(nmax)
      dimension pbar(nmax),fn(nmax)
C
      IERRHY = 0
      NP1 = NPTS + 1
C
      DT = DTIN
      DO 100 ITRY=1,ITRYMX
         PBAR(1) = AHF*(PN1(1) + P(1))
         GN(1) = ZERO
         FN(1) = EP(1) - PBAR(1)*(VN1(1)-VP(1)) -
     $     DT*((TOTHS*(RLUMN1(2)+RLCON(2)-CORLUM) +
     $           ONTH*(RLUMN(2)+RLCON(2)-CORLUM))/DM1(1) -
     $           SORCE(1) )
         DO 10 I=2,NPTS
            PBAR(I) = AHF*(PN1(I) + P(I))
            RBAR2 = (RN1(I)**2 + RN1(I)*RP(I) + RP(I)**2)*ONTH
            RBARM2 = 1.D0/(RN1(I)*RP(I))
            GN(I) = RP(I) + DT*UP(I) - AHF*DT*DT*
     $         (G*RM(I)*RBARM2 + PI4*RBAR2*(PBAR(I)-PBAR(I-1))/DM2(I))
            FN(I) = EP(I) - PBAR(I)*(VN1(I)-VP(I)) -
     $      DT*((TOTHS*(RLUMN1(I+1)-RLUMN1(I)) +
     $           ONTH*(RLUMN(I+1)-RLUMN(I)) )/DM1(I) +
     $          (TOTHS*(RLCON(I+1)-RLCON(I)) +
     $           ONTH*(RLCON(I+1)-RLCON(I)) )/DM1(I) -
     $           SORCE(I) )
  10     CONTINUE
         RBAR2 = (RN1(NP1)**2 + RN1(NP1)*RP(NP1) +RP(NP1)**2)*ONTH
         RBARM2 = 1.D0/(RN1(NP1)*RP(NP1))
         PNPTS = AHF*(PN1(NPTS) - AD3*TN1(NPTS)**4 +
     $                P(NPTS) - AD3*TP(NPTS)**4)
         GN(NP1) = RP(NP1) + DT*UP(NP1) - AHF*DT*DT*
     $      (G*RM(NP1)*RBARM2 - PI4*RBAR2*PNPTS/DM2(NP1))
         IF( ITRY .GT. ITRYMX/2 ) THEN
            IOUT = 3
         ELSE
            IOUT = 0
         ENDIF
         CALL NDERIV(IOUT,IERR,RETERR,DT,PBAR,GN,FN,nmax)
         IF( IERR .LE. 0 ) THEN
C
C          RADII CROSSED IN CONVERGENCE ROUTINE, DECREASE TIME STEP
C       AND TRY AGAIN.
C
            DT = DT/2.D0
            RN1(1) = RP(1) + UP(1)*DT
            RLUMN1(1) = RLUMN(1)
            UN1(1) = UP(1)
            DO 15 K=1,NPTS
               UN1(I+1) = UP(I+1)
               RN1(I+1) = RP(I+1) + UP(I+1)*DT
               VN1(I) = VP(I)
               TN1(I) = TP(I)
               PN1(I) = P(I)
               EN(I) = EP(I)
               RLUMN1(I+1) = RLUMN(I+1)
  15        CONTINUE
            GOTO 100
         ENDIF
c         WRITE(11,1000) RETERR,ACCUR
         DO 20 I=2,NPTS+1
            UN1(I) = TWO*(RN1(I) - RP(I))/DT - UP(I)
  20     CONTINUE
         DO 30 I=1,NPTS
            IF( I .LT. NPTS ) THEN
c               IF( TN1(I+1) .GT. TN1(I) ) CALL MODDMP(I,DT)
               IF( TN1(I+1) .GT. TN1(I) ) write(11,3010) I
            ENDIF
            CALL STATE (I,1)
            CALL RADLUM(I,0)
            CALL CONLUM(I, 0, RLUMC)
            CS(I) = DSQRT(G1(I)*VN1(I)*PN1(I))
  30     CONTINUE
         DTIN = DT
         IF( RETERR .LE. ACCUR ) THEN
            WRITE(11,3000) ITRY,RETERR,ACCUR
            IERRHY = 1
            RETURN
         ENDIF
 100  CONTINUE
      CALL MODDMP(I,DT)
      WRITE(6,1010) ITIME,DT
      IERRHY =-1
      RETURN
C
 1000 FORMAT(1X,20HFROM HYDRO...ERROR =,1PE11.3,11H ACCURACY =,E11.3)
 1010 FORMAT(1X,19HNO GO IN TIME STEP ,I3,8H DTIME =,1PE11.3)
 3000 FORMAT(1X,'Converged after ',i3,' steps to ',1PE10.3,
     $          ' accuracy =',e10.3)
 3010 format(1x,'Temperature inversion in zone ',i4)
      END
      FUNCTION QTIME(R,CS,NMAX,N,IPRTIN)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL IPRT,IPRTIN
C
C          RETURNS THE MINIMUM TIME STEP ALLOWED USING AN
C       ACOUSTIC PROPAGATION VELOCITY ACROSS EACH ZONE.
C
      DIMENSION R(NMAX),CS(NMAX)
      IPRT = IPRTIN
      REG = (R(3)-R(2))/CS(2)
      IS = 2
      DO 10 I=3,N
         IF( (R(I+1)-R(I))/CS(I) .GT. REG ) GOTO 10
         REG = (R(I+1)-R(I))/CS(I)
         IS = I
  10  CONTINUE
      QTIME = REG
C      WRITE(6,1000) IS,REG
      IF( IPRT ) WRITE(11,1000) IS,REG,R(IS),R(IS-1),CS(IS)
      RETURN
 1000 FORMAT(21H MINIMUM TIME IN ZONE,I4,11H TIME STEP=,1PE10.3,
     $   3X,1P,3E11.4)
      END
      SUBROUTINE NDERIV(IO,IERR,RETERR,DTIN,PBAR,GN,FN,MMAX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     *****************************************************************
C      PULSATION EQUATION SOLVER
C      BASIC REFERENCES -- CASTOR, AP.J. 166 (109) 1971.
C                          PESNELL, 1983, THESIS, UNIV. OF FLORIDA.
C                          PESNELL, PASP, 99 (975), 1987.
C     *****************************************************************
C
C
C        ARRAYS ARE REARRANGED FROM THE ORIGINAL LNA CODE
C
      parameter ( onth=1.d0/3.d0, toths=2.d0/3.d0 )
      parameter ( nmax=256, nmaxvc=2*nmax )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),p(nmax),up(nmax),
     $               ep(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common xs(nmaxvc)
      common/blk1/   irad,nflag,n
      common/blk8/   g,pi,pi2,pi4,pi43
      common/blk37/  useno(nmax,2), delr(nmax),delt(nmax)
      common/lumins/ frft(nmax),sorce(nmax),dtsorc(nmax),dvsorc(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/stlums/ rlumn1(nmax),rlumn(nmax), rlcon(nmax)
      common/thermo/ r,a,bk,avagd,ad3
c
      dimension ag1(nmax,3),ag2(nmax,2),ak1(nmax,4),ak2(nmax,3),
     $       bl2(nmax,2),bl1(nmax,3),dr(nmax,2)
      dimension pbar(mmax),gn(mmax),fn(mmax)
      dimension amat(nmaxvc,7)
      dimension szero(nmax,19)
      equivalence (szero(1,1),ag1(1,1))
C
C   ZERO COMMON BLOCKS SCRTCH AND (BLANK).
C
      DO 5 J=1,NMAXVC
         XS(J) = ZERO
   5  CONTINUE
      DO 6 J=1,19
         DO 6 I=1,NMAX
            SZERO(I,J) = ZERO
   6  CONTINUE
      DT = DTIN
      NP = N+1
      NM1 = N-1
C
      RLUMX = CORLUM
      AFAC = PI43*(RN1(1)**2+RP(1)*RN1(1)+RP(1)**2)/DM2(1)
      DO 60 I=1,NP
         AFACP = AFAC
         IF( I .EQ. NP ) GOTO 40
         AFAC = PI43*(RN1(I+1)**2+RP(I+1)*RN1(I+1)+RP(I+1)**2)/DM2(I+1)
         DR(I,1) = PI4*RN1(I)**2/VN1(I)/DM1(I)
         IF( I .EQ. 1 ) DR(I,1) = ZERO
         DR(I,2) =-PI4*RN1(I+1)**2/VN1(I)/DM1(I)
         IF( I.GT.1 ) GOTO 20
         BL1(I,1) = ZERO
         BL1(I,2) = ZERO
         BL1(I,3) = ZERO
         BL2(I,1) = ZERO
         BL2(I,2) = ZERO
         GOTO 60
  20    CONTINUE
C
C          SET UP THE MECHANICAL MATRICES.
C
         AG1(I,1) = AFACP*AHF*CHR(I-1)*PN1(I-1)*DR(I-1,1)
         IF( I .EQ. 2 ) AG1(2,1) = ZERO
         RBARM2 = 1.D0/(RN1(I)*RP(I))
         AG1(I,2) = G*RM(I)*RBARM2/RN1(I) -PI43*(TWO*RN1(I)+RP(I))*
     $      (PBAR(I)-PBAR(I-1))/DM2(I) -
     $    AFACP*(CHR(I)*PN1(I)*DR(I,1) -
     $               CHR(I-1)*PN1(I-1)*DR(I-1,2) )*AHF
         AG1(I,3) =-AFACP*AHF*CHR(I)*PN1(I)*DR(I,2)
         AG2(I,1) = AFACP*AHF*CHIT(I-1)*PN1(I-1)/TN1(I-1)
         AG2(I,2) =-AFACP*AHF*CHIT(I)*PN1(I)/TN1(I)
C
C          INITIALIZE THE THERMAL MATRICES, IF IRAD=1 USE THE
C       STELLINGWERF INTERPOLATION FORMULA.
C
         DTL1 = TN1(I-1)**4
         DTL2 = TN1(I)**4
         T4OKI = DTL2/AKAP(I)
         T4OKI1 = DTL1/AKAP(I-1)
         DIFF = T4OKI-T4OKI1
         WIOWI1 = DLOG(DTL2/DTL1)
         WOWSQ = WIOWI1**2
         GKOGK1 = DLOG(AKAP(I)/AKAP(I-1))
         DENOM = ONE - GKOGK1/WIOWI1
         BL1(I,1) = DR(I-1,1)*( DKDR(I-1)*(T4OKI1/DIFF)
     $  - (DKDR(I-1)/WIOWI1)/DENOM)
         IF( I .EQ. 2 ) BL1(2,1) = ZERO
         BL1(I,2) = FOR/RN1(I) +
     $     DR(I,1)*( -DKDR(I)*(T4OKI/DIFF) + (DKDR(I)/WIOWI1)/DENOM) +
     $  DR(I-1,2)*( DKDR(I-1)*(T4OKI1/DIFF) - (DKDR(I-1)/WIOWI1)/DENOM )
         BL1(I,3) = DR(I,2)*( -DKDR(I)*(T4OKI/DIFF)
     $     + (DKDR(I)/WIOWI1)/DENOM)
         BL2(I,1) = ((-FOR+DKDT(I-1))*(T4OKI1/DIFF)
     $  - (DKDT(I-1)/WIOWI1-FOR*GKOGK1/WOWSQ)/DENOM)/TN1(I-1)
         BL2(I,2) = ((FOR-DKDT(I))*(T4OKI/DIFF)
     $  + (DKDT(I)/WIOWI1-FOR*GKOGK1/WOWSQ)/DENOM)/TN1(I)
         GOTO 50
  40    CONTINUE
         TAU = AKAP(N)*DM1(N)/(FOR*PI4*RN1(NP)*RN1(NP))
         TAU = 1.5D0*TAU/(ONE+1.5D0*TAU)
         BL1(NP,1) = DR(N,1)*( -TAU*(ONE+DKDR(N)) )
         BL1(NP,2) = TWO*(ONE-TAU)/RN1(NP) +
     $      DR(N,2)*( -TAU*(ONE+DKDR(N)))
         BL1(NP,3) = ZERO
         BL2(NP,1) = (FOR - TAU*DKDT(N))/TN1(N)
         BL2(NP,2) = ZERO
         PNPTS = AHF*(PN1(N)-AD3*TN1(N)**4+P(N)-AD3*TP(N)**4)
         AFACP = PI43*(RN1(NP)**2+RP(NP)*RN1(NP)+RP(NP)**2)/DM2(NP)
         RBARM2 = 1.D0/(RN1(NP)*RP(NP))
C
         AG1(NP,1) = AFACP*AHF*CHR(N)*PN1(N)*DR(N,1)
         AG1(NP,2) = G*RM(NP)*RBARM2/RN1(NP) +
     $   AFACP*CHR(N)*PN1(N)*DR(N,2)*AHF +
     $  PI43*(TWO*RN1(NP)+RP(NP))*PNPTS/DM2(NP)
         AG1(NP,3) = ZERO
C
         PNCHIT = CHIT(N)*PN1(N)/TN1(N)
         PNCHIT = PNCHIT - 4.d0*AD3*TN1(N)**3
         AG2(NP,1) = AFACP*AHF*PNCHIT
         AG2(NP,2) = ZERO
  50    CONTINUE
         GLOM = RLUMN1(I-1)/DM1(I-1)
         GLOM1 = RLUMN1(I)/DM1(I-1)
         DEDR = DVSORC(I-1)*SORCE(I-1)
         DEDS = DTSORC(I-1)*SORCE(I-1)/TN1(I-1)
         AK1(I-1,1) = GLOM*BL1(I-1,1)*TOTHS
         AK1(I-1,2) = (GLOM*BL1(I-1,2)-GLOM1*BL1(I,1))*TOTHS
c     $                - DEDR*DR(I-1,1)
         AK1(I-1,3) = (GLOM*BL1(I-1,3)-GLOM1*BL1(I,2))*TOTHS
c     $                - DEDR*DR(I-1,2)
         AK1(I-1,4) =-GLOM1*BL1(I,3)*TOTHS
         AK2(I-1,1) = GLOM*BL2(I-1,1)*TOTHS
         AK2(I-1,2) = (GLOM*BL2(I-1,2)-GLOM1*BL2(I,1))*TOTHS
c     $               + DEDS
         AK2(I-1,3) =-GLOM1*BL2(I,2)*TOTHS
  60  CONTINUE
      GLOM = RLUMN1(NP)/(AHF*DM1(N))
C     AK1(NP,1) = GLOM*BL1(NP,1)
      AK1(NP,1) = ZERO
C     AK1(NP,2) = GLOM*BL1(NP,2)
      AK1(NP,2) = ZERO
      AK1(NP,3) = ZERO
      AK1(NP,4) = ZERO
C     AK2(NP,1) = GLOM*BL2(NP,1)
      AK2(NP,1) = ZERO
      AK2(NP,2) = ZERO
      AK2(NP,3) = ZERO
C
C          WRITE THE MATRICES TO FILE TAPE11, IF IO IS LESS
C       THAN 3, THIS WRITE IS NOT DONE.
C
      IF( IO .GE. 3 ) THEN
         WRITE(11,7000)
         DO 75 I=1,NP
         WRITE(11,7001) I,AG1(I,1),AG1(I,2),AG1(I,3),AG2(I,1),AG2(I,2),
     $    AK1(I,1),AK1(I,2),AK1(I,3),AK1(I,4),AK2(I,1),AK2(I,2),AK2(I,3)
  75     CONTINUE
      ENDIF
C
C          LOAD MATRIX TO CALCULATE THE CORRECTIONS TO RADIUS AND
C       TEMPERATURE
C
      XS(1) = ZERO
      DO 140 I=1,N
         IY = 2*I-1
         IX = 2*I
         XS(IX) = RN1(I+1) - GN(I+1)
         XS(IY) = EN(I) - FN(I)
C
C          AMAT RESET EVERY ITERATION
C
         AMAT(IY,1) =-DT*AK1(I,1)
         AMAT(IY,2) =-DT*AK2(I,1)
         AMAT(IY,3) = (PBAR(I)*VN1(I)+PN1(I)*((CHIT(I)-ONE)*VN1(I) -
     $                 AHF*CHR(I)*(VN1(I)-VP(I))) )*DR(I,1)
     $               - DT*AK1(I,2)
         AMAT(IY,4) = CV(I) + AHF*PN1(I)*CHIT(I)*(VN1(I)-VP(I))/TN1(I) -
     $                 DT*AK2(I,2)
         AMAT(IY,5) = (PBAR(I)*VN1(I)+PN1(I)*((CHIT(I)-ONE)*VN1(I) -
     $                 AHF*CHR(I)*(VN1(I)-VP(I))) )*DR(I,2)
     $               - DT*AK1(I,3)
         AMAT(IY,6) =-DT*AK2(I,3)
         AMAT(IY,7) =-DT*AK1(I,4)
         AMAT(IX,1) = ZERO
         AMAT(IX,2) =-AHF*DT**2*AG1(I+1,1)
         AMAT(IX,3) =-AHF*DT**2*AG2(I+1,1)
         AMAT(IX,4) = ONE - AHF*DT**2*AG1(I+1,2)
         AMAT(IX,5) =-AHF*DT**2*AG2(I+1,2)
         AMAT(IX,6) =-AHF*DT**2*AG1(I+1,3)
         AMAT(IX,7) = ZERO
 140  CONTINUE
      CALL RBMLES(AMAT,NMAXVC,1,2*N,7,XS)
c
c          Guarantee that the corrections are too small to allow
c       radii to cross or a temperature inversion.
c
      DTMAX = ZERO
      DTERR = ZERO
      DRMAX = ZERO
      DRERR = ZERO
      DO 180 I=1,N
         IF( I .EQ. 1 ) THEN
            DELT(I) = XS(2*I-1)/TN1(I)
            DELR(I+1) = XS(2*I)/(RN1(I+1)-RN1(I))
         ELSE
            DELT(I) = (XS(2*I-1)-XS(2*I-3))/(TN1(I)-TN1(I-1))
            DELR(I+1) = (XS(2*I)-XS(2*I-2))/(RN1(I+1)-RN1(I))
         ENDIF
         DTMAX = DMAX1(DTMAX,DABS(XS(2*I-1)/TN1(I)))
         DTERR = DTERR + (XS(2*I-1)/TN1(I))**2
         DRMAX = DMAX1(DRMAX,DABS(DELR(I+1)))
         DRERR = DRERR + DELR(I+1)**2
 180  CONTINUE
C
C          LIMIT MAXIMUM RELATIVE CORRECTION TO FMAX
C
      FMAX = 0.10D0
      FFAC = ONE
      IF( DTMAX .GT. FMAX .OR. DRMAX .GT. FMAX )
     $   FFAC = FMAX/DMAX1(DTMAX,DRMAX)
C
C          Find the new radii and temperatures.
C
      DO 185 I=1,N
         TN1(I) = TN1(I) - FFAC*XS(2*I-1)
         RN1(I+1) = RN1(I+1) - FFAC*XS(2*I)
 185  CONTINUE
C
C          Finish up by finding the new specific volume and thermodynamic
c       quantities.
C
      DO 190 I=1,N
         IF( RN1(I+1) .LE. RN1(I) ) THEN
            WRITE(6,'(1x,i4,1p2e15.5)') I,RN1(I+1),RN1(I)
            GOTO 900
         ENDIF
         VN1(I) = PI43*(RN1(I+1)**3 - RN1(I)**3)/DM1(I)
         CALL STATE(I,1)
         CALL RADLUM(I,1)
         CALL CONLUM(I, 1, RLUMC)
         IF( IO .GT. 2 )
     $      WRITE(11,1400) I,RN1(I+1),XS(2*I),VN1(I),TN1(I),XS(2*I-1)
 190  CONTINUE
      RETERR = DSQRT(DTERR)
      IERR = 1
C      WRITE(11,1900) DTMAX,DRMAX
C
      RETURN
C
C          ERRORS DETECTED, PROGRAM REDIRECTED
C
 900  CONTINUE
      IERR =-1
      WRITE(6,9000) DT
      RETURN
C
 1400 FORMAT(1X,I4,1P,7E11.4)
 1900 FORMAT(1X,21HFROM NDERIV...DTMAX =,1PE11.3,8H DRMAX =,E11.3)
 7000 FORMAT(1H1,/1X,43HPULSATION MATRICES: AG1(3), AG2(2), AK1(4),,
     $              12H AND AK2(3).)
 7001 FORMAT(1X,I4,1X,1P,3E10.3,1X,2E10.3,1X,4E10.3,1X,3E10.3)
 9000 FORMAT(1X,26HRADII CROSSED AT TIME STEP,1PE11.4,
     $       17H RETURN TO HYDRO.)
      END
      SUBROUTINE RADLUM(I,IOPT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     *****************************************************************
C      RADLUM CALCULATES THE RADIATIVE LUMINOSITY AND ROSSELAND OPACITY
C     *****************************************************************
C
      parameter ( SIGMA=5.66961D-5 )
      parameter ( AC3=7.5595D-5 )
      parameter ( onth=1.d0/3.d0, toths=2.d0/3.d0 )
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),p(nmax),un1(nmax),en(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/blk1/   irad,nflag,npts
      common/stlums/ rlumn1(nmax),rlumn(nmax), rlcon(nmax)
      common/blk8/   grav,pi,pi2,pi4,pi43
      common/inhom/  x,y,z,rmuc,pargrm
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      VPI = VN1(I)
      TPI = TN1(I)
C
C  GET OPACITY ONLY IF IOPT=-1
C
      CALL STLOPC(TPI,VPI,AKAPI,DKDT(I),DKDR(I),IOPT)
      AKAP(I) = AKAPI
      IF( IOPT .NE. 0 ) RETURN
C
C          RETURN THE LUMINOSITY EXITING ZONE I.
C
      IF( I .EQ. NPTS ) THEN
C
C          OUTER ZONE IS CONVERGED TO THE DESIRED OPTICAL DEPTH.
C
         TAU = AKAPI*DM1(NPTS)/(FOR*PI4*RN1(NPTS+1)**2)
         RLUM = PI4*RN1(NPTS+1)**2*SIGMA*TPI**4/(0.75D0*(TAU+TOTHS))
         RLUMN1(I+1) = RLUM
      ELSE
C
C       RATIO:=AC/3*(4PI*R(EFF)**2)**2/(KAPPA(EFF)*MASS(EFF))
C      SUCH THAT RATIO*(TP(I)**4-TP(I+1)**4) = LUM.(I+1)
C        SEE STELLINGWERF AP. J. 195(441)1975
C
         W  = TPI**4
         WP = TN1(I+1)**4
         XLNW = DLOG(WP/W)
         XLNK = DLOG(AKAP(I+1)/AKAPI)
         C1 =-(PI4*RN1(I+1)**2)**2/DM2(I+1)
         C2 = WP/AKAP(I+1) - W/AKAPI
         C2 = C2/(ONE-XLNK/XLNW)
         RLUMN1(I+1) = AC3*C1*C2
      ENDIF
      RETURN
      END
      SUBROUTINE CONLUM(I, IOPT, RLUMC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C      CONLUM CALCULATES THE CONVECTIVE LUMINOSITY
C
      parameter ( onth=1.d0/3.d0, toths=2.d0/3.d0 )
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),p(nmax),un1(nmax),en(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/blk1/   irad,nflag,npts
      common/stlums/ rlumn1(nmax),rlumn(nmax), rlcon(nmax)
      common/blk8/   grav,pi,pi2,pi4,pi43
      common/inhom/  x,y,z,rmuc,pargrm
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
c
      VPI = VN1(I)
      TPI = TN1(I)
C
      RLUMC = RLCON(I+1)
C
      RETURN
      END
      SUBROUTINE STATE(I,IOPT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C    ******************************************************************
C            STATE CALCULATES THE EQUATION OF STATE
C    ******************************************************************
C            REQUIRES THE ROUTINE NEOS FOR THE ANALYTIC EOS
C            REQUIRES THE ROUTINE INTERP FOR THE TABULAR EOS
C
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/blk17/  pfi,dtp,dvp,pe,pet,pev
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/thermo/ r,a1,bk,avagd,ad3
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/inhom/  x,y,z,rmuc,pargrm
c
      TPI = TN1(I)
      PRAD = AD3*TPI**4
      VPI = VN1(I)
C
C          THE ANALYTIC EOS
C
      KUSE = 5
      IF( IOPT .EQ. 1 ) KUSE = 2
      CALL NEOS(TPI,VPI,PFI,DVP,DTP,E,DVE,DTE,BP,DBDV,DBDT,PE,PET,
     $       PEV,FRE,ENT,KUSE)
      PN1(I) = PFI
      EN(I) = E
      IF( IOPT .EQ. 1 ) THEN
         DBDV = DBDV*VPI/BP
         DBDT = DBDT*TPI/BP
         BETAI = ONE - PRAD/PFI
         CV(I) = DTE
C
C          FIND THE ADIABATIC EXPONENTS
C
         G3M1(I) = DTP*VPI/DTE
         G1(I) = (TPI*DTP*G3M1(I)-VPI*DVP)/PFI
         DELAD(I) = G3M1(I)/G1(I)
C
C          EQUATION 27.28 COX + GIULI
C
         CHITI = FOR-THRE*BETAI+BETAI*DBDT
         CHIRHO = BETAI*(ONE - DBDV)
         CHR(I) = CHIRHO
         CHIT(I) = CHITI
         Q(I) = CHITI/CHIRHO
         CON = PFI*VPI/TPI
C          EQUATION 9.86 COX + GIULI
         CP(I) = DTE+CON*CHITI*Q(I)
C         ENTPRY(I) = ENT
      ENDIF
      RETURN
      END
      SUBROUTINE NEOS(TIN,VIN,P,PV,PT,E,EV,ET,BP,BPV,BPT,PE,PET,PEV,
     $          FRE,ENT,KUSE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     EQUATION OF STATE, OPACITY
C        FIRST ORDER ONLY
C     ARGUMENTS...TIN (DEGREES K), VIN=1/RHO (CM**3/GM)
C     METALS...
C        NA,AL ALWAYS IONIZED
C        MG,SI,FE INCLUDED AS SINGLE ELEMENT
C        ALL OTHERS IGNORED
C      KUSE=0, PRESSURE, INTERNAL ENERGY AND DERIVATIVES
C      KUSE=1, PRESSURE, INTERNAL ENERGY, OPACITY AND DERIVATIVES
C      KUSE=2, PRESSURE, INTERNAL ENERGY, ENTROPY AND DERIVATIVES
C      KUSE=3, PRESSURE, INTERNAL ENERGY, ENTROPY (WITH DERIVATIVES)
C      KUSE=5, PRESSURE AND DERIVATIVES
C
C
C          TABLE OF RETURNED QUANTITIES.
C
C      FUNCTION       NAME   DERIVATIVE WITH RESPECT TO
C                               TEMP.       SP. VOL.
C-------------------------------------------------------
C      PRESSURE     I    P I     PT      I     PV      I
C      INT. ENERGY  I    E I     ET      I     EV      I
C      ELEC. PRS.   I   PE I     PET     I     PEV     I
C      MOLAR DENSITYI   BP I     BPT     I     BPV     I
C      ENTROPY      I  ENT I     ---     I     ---     I
C      FREE ENERGY  I  FRE I     ---     I     ---     I
C-------------------------------------------------------
C
C          X = HYDROGEN MASS FRACTION
C          Y = HELIUM MASS FRACTION
C          Z = METALLIC MASS FRACTION
C          PARGRM = MEAN MOLECULAR MOLAR DENSITY WITHOUT ELECTRONS
C
      COMMON/INHOM/  X,Y,Z,RMUC,PARGRM
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
C
C          R = GAS CONSTANT 8.31434E7
C          A = STEFAN-BOLTZMAN CONSTANT 7.56471E-15
C          BK = BOLTZMAN S CONSTANT 8.6170837E-5
C          AVAGD = AVAGADRO S NUMBER 6.02217E23
C          AD3 = A/3
C
      COMMON/THERMO/ R,A,BK,AVAGD,AD3
      DIMENSION QMASP(10),PARTN(10)
      DATA PARTN/ 10*0.0D0 /
      DATA QMASP/ 2.023958D0,1.011979D0,8.0078D0,16.015D0,
     $  8.0078D0,124.768D0,310.026D0,73.961D0,2.56984D-5,0.0D0 /
      DATA T3OUT,T4OUT/1.665795163D-25,3.802592017D-28/
      DATA T2OUT,T5OUT/ 5.347896D-35,6.614536D-34/
C
C          IONIZATION POTENTIALS FOR HYDROGEN AND HELIUM
C
      PARAMETER (XH=13.595D0, XHE=24.581D0, XHE2=54.403D0)
      PARAMETER ( PREC=1.D-12, THREHLF=1.5D0 )
      PARAMETER (C1=4.0092926D-9,C2=1.00797D0,C3=4.0026D0,C4=21.02447D0)
      PARAMETER ( XM=7.9D0, CM=0.7D0, ZPZP=0.12014D0 )
C
      V = VIN
      T = TIN
      IF( V .LE. ZERO ) GO TO 11
      IF( T .LE. ZERO ) GO TO 10
      FRE = ZERO
      ENT = ZERO
      RT = R*T
      TT4 = T**4
      TK = ONE/(T*BK)
      SQT = DSQRT(T)
C    C1=ORIGINAL(C1(0.33334622))/R
      T1 = V*SQT**3*C1
      T2 = T2OUT
      IF( T .GT. 2.D3 ) T2 = DEXP(-XH*TK)
      T3 = T3OUT
      IF( T .GT. 5.D3 ) T3 = DEXP(-XHE*TK)
      T4 = T4OUT
      IF( T .GT. 1.D4 ) T4 = DEXP(-XHE2*TK)
      T5 = T5OUT
      IF( T .GT. 1.2D3 ) T5 = DEXP(-XM*TK)
      D = T1*T2
      B = FOR*T1*T3
      C = B*T1*T4
      DD = TWO*CM*T1*T5
      ZNA = Z*2.48D-3/24.969D0
      ZMG = Z*ZPZP/45.807D0
C
C          CONVERGE ON ELECTRON DENSITY USING THE SAHA EQUATION.
C
C          GES IS THE MOLAR DENSITY OF ELECTRONS.
C
      GES = (X+Y*AHF)/(ONE+Y/(FOR*C))
      IF(GES.LT.X) GES = AHF*(DSQRT(D*(D+FOR*X))-D)
      IF( GES .LT. 1.D-6*Z ) GES = 1.D-6*Z
      XC2 = X/C2
      YC3 = Y/C3
C
C
C          NEWTON METHOD FOR ELECTRON DENSITY.
C
      DO 1 I=1,25
         T2 = C/GES+GES+B
         GEP = XC2*D/(GES+D)+YC3*(B+TWO*C/GES)/T2
     $        + ZMG*DD/(GES+DD) + ZNA
         T1 = ONE+XC2*D/(D+GES)**2+YC3/T2*
     $     (TWO*C/GES**2+(B+TWO*C/GES)*(ONE-C/GES**2)/T2)
     $     + ZMG*DD/(GES+DD)**2
         DGES = (GEP-GES)/T1
         GES = GES+DGES
         IF( DABS(DGES)/GES .LT. PREC ) GOTO 3
   1  CONTINUE
      GOTO 12
   3  CONTINUE
C
C  ELECTRON PRESSURE
C
      PE = RT*GES/V
C
C      TOTLN = 1/MU = X/C2+Y/C4+Z/C3+GES
C
      TOTLN = PARGRM+GES
      XX = D/(GES+D)
      T2 = GES+B+C/GES
      YY = B/T2
      ZZ = C/(GES*T2)
      WW = DD/(GES+DD)
C
C          DERIVATIVES OF THE SAHA EQUATION FOR THE PRESSURE AND
C       INTERNAL ENERGY TEMPERATURE AND DENSITY DERIVATIVES.
C
      T1 = YC3*(B+TWO*C/GES)
      QC0 = ONE+XC2*XX/(GES+D)+ZMG*WW/(GES+DD)+YC3/T2*
     $      (TWO*C/GES**2+(B+TWO*C/GES)*(ONE-C/GES**2)/T2)
      QC1 = XC2*(ONE-XX)/(GES+D)
      QC4 = ZMG*(ONE-WW)/(GES+DD)
      QC2 = (YC3-T1/T2)/T2
      QC3 = (YC3*TWO-T1/T2)/(GES*T2)
      QGV = (QC1*D+QC2*B+QC3*TWO*C+QC4*DD)/(QC0*V)
      QP1 = D*(THREHLF+XH*TK)/T
      QP2 = B*(THREHLF+XHE*TK)/T
      QP3 = C*(THRE+(XHE+XHE2)*TK)/T
      QP4 = DD*(THREHLF+XM*TK)/T
      QGT = (QC1*QP1+QC2*QP2+QC3*QP3+QC4* QP4)/QC0
C
C          ELECTRON PRESSURE DERIVATIVES.
C
      PET = ONE + T*QGT/GES
      PEV =-ONE + V*QGV/GES
C
C          PRESSURE DUE TO THE IDEAL GAS
C
      P = RT*TOTLN/V
      PT = P/T+RT*QGT/V
      PV = RT*QGV/V-P/V
C
C          BP IS R/MU
C
      BP = R*TOTLN
      BPV = R*QGV
      BPT = R*QGT
C
C           ADD THE RADIATION PRESSURE
C
      P = P+AD3*TT4
      PT = PT+FOR*AD3*TT4/T
      IF( KUSE .EQ. 5 ) RETURN
C
C          IONIZATION ENERGY
C
      EI = (R/BK)*(XH*XX*XC2+YC3*(XHE*YY+(XHE+XHE2)*ZZ)
     $      +ZMG*XM*WW + ZNA*5.524D0 )
C          TOTAL INTERNAL ENERGY
      E = THREHLF*RT*TOTLN + A*V*TT4 + EI
      EV = T*PT-P
      QXT = ((ONE-XX)*QP1-XX*QGT)/(GES+D)
      DT2 = QGT*(ONE-C/GES**2)+QP2+QP3/GES
      QYT = (QP2-B*DT2/T2)/T2
      QZT = (QP3-C*QGT/GES-C*DT2/T2)/(T2*GES)
      QWT = ((ONE-WW)*QP4-WW*QGT)/(GES+DD)
      EIT = (R/BK)*(XH*QXT*XC2+YC3*(XHE*QYT+(XHE+XHE2)*QZT)+
     $     ZMG*XM*QWT)
      ET = THREHLF*R*(TOTLN+T*QGT)+FOR*A*V*TT4/T+EIT
      IF( KUSE .EQ. 2 ) GOTO 30
      RETURN
C
C            ENTROPY AND FREE ENERGY
C     SEE ZELDOVICH & RAZIER,PHYSICS OF SHOCK WAVES AND
C     HIGH TEMP. PHENOMENA,1966,PP. 192-197
C
  30  CONTINUE
C      7.76E11=SQRT(AVAGADRO S NUMBER)
      T32 = SQT**3*V/7.76026D11
      T32 = DLOG(T32)*TOTLN
      PARTN(2) = XC2*XX
      PARTN(4) = YC3*YY
      PARTN(5) = YC3*ZZ
      PARTN(6) = ZNA
      PARTN(7) = ZMG*WW
      PARTN(1) = XC2 - PARTN(2)
      PARTN(3) = YC3-PARTN(4)-PARTN(5)
      PARTN(8) = Z/C4-ZMG*WW-ZNA
      PARTN(9) = GES
      ENTX = 8.38508D0*TOTLN
      ETA = ONE
      DO 40 J=1,10
         IF( PARTN(J) .LE. ZERO ) GOTO 40
         ETA = ETA * (PARTN(J)/QMASP(J))**PARTN(J)
  40  CONTINUE
      ETA = DLOG(ETA)
      ENT = (2.5D0*TOTLN+ENTX-ETA+T32)*R
     $      +(FOR*AD3)*V*T**3
C      FRE = (ETA-T32-ENTX-TOTLN)*R*T-AD3*V*T**4
      RETURN
C
C     DIAGNOSTICS
C
  10  WRITE(6,1000) T,V
 1000 FORMAT(/,1X,'NEGATIVE TEMP. IN NEOS,   T= ',1PE10.3,2X,2HV=,E10.3)
      STOP
  11  WRITE(6,1010) T,V
 1010 FORMAT(/,1X,'NEGATIVE DENSITY IN NEOS  T= ',1PE10.3,2X,2HV=,E10.3)
      STOP
  12  WRITE(6,1020) T,V
 1020 FORMAT(/,1X,'NO CONVERGENCE IN NEOS    T= ',1PE10.3,2X,2HV=,E10.3)
      STOP
      END
      SUBROUTINE STLOPC(TIN,VIN,AKAP,DKDT,DKDR,IOPT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C          STELLINGWERF KING4A OPACITY FIT.
C          IF IOPT = 0, RETURNS ONLY THE OPACITY.
C          IF IOPT = 1, RETURN THE OPACITY AND DERIVATIVES.
C
C      ARGUMENTS...TIN (DEGREES K), VIN=1/RHO (CM**3/GM)
C          X = HYDROGEN MASS FRACTION
C          Y = HELIUM MASS FRACTION
C          Z = METALLIC MASS FRACTION
C          PARGRM = MEAN MOLECULAR MOLAR DENSITY WITHOUT ELECTRONS
C
      COMMON/INHOM/  X,Y,Z,RMUC,PARGRM
      COMMON/BLK17/  PFI,DTP,DVP,PE,PET,PEV
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
C
C          OPACITY IS FROM A FIT BY STELLINGWERF AP. J. 195(441)1975
C       AND AP. J. 195(705)1975.
C
      ZVAR1 = 21.D0*Z + 0.979D0
      ZVAR2 = 105.D0*Z + 0.895D0
      YVAR1 =-6.D-5*Y + 6.294D-5
      YVAR2 = 4.53D6*Y - 3.0447D5
      V = VIN
      T = TIN
      OP = OPAC1(T,V,ZVAR1,ZVAR2,YVAR1,YVAR2)
      AKAP = PE*OP
C     DKDR = ZERO
C     DKDT = ZERO
      IF( IOPT .EQ. 0 ) RETURN
      V = V*(ONE+1.D-8)
C
C    DKDR=-D(LOG(KAPPA))/D(LOG(V))=D(LOG(KAPPA))/D(LOG(RHO))
C
      OP1 = OPAC1(T,V,ZVAR1,ZVAR2,YVAR1,YVAR2)
      DKDR =-PEV + (OP-OP1)/OP*1.D8
      V = VIN
      T = T*(ONE+1.D-8)
C
C    DKDT=D(LOG(KAPPA))/D(LOG(T))
C
      OP1 = OPAC1(T,V,ZVAR1,ZVAR2,YVAR1,YVAR2)
      DKDT = PET + (OP1-OP)/OP*1.D8
      RETURN
      END
      FUNCTION OPAC1(T,V,ZVAR1,ZVAR2,YVAR1,YVAR2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
C
C     GENERAL FIT TO KING OPACITY TABLES
C
      U = T*1.D-4
      U2 = U*U
      U4 = U2*U2
      U6 = U4*U2
      SQU = DSQRT(U)
      UM8 = ONE/(U4*U4)
      UM10 = UM8/U2
      U25 = SQU*U2
      U35 = U*U25
      UM45 = ONE/(U*U35)
      U5 = U6/U
      V1 = V**0.35D0
      V2 = DSQRT(V1)
      T1 = 760.D0*U5+316.D0/V2
      T1 = YVAR1*V1*U35+ONE/T1
      T1 = ONE/(TEN*U6+ONE/T1)
      T2 = ZVAR1*YVAR2*UM10+2.13D-3*V2*ZVAR2*UM45
c
C          OLD STELLINGWERF FORMULA
C      T2 = 1780.E0*U25/ZVAR1+ONE/T2
c
c           NEW STELLINGWERF FORMULA
c
      T2 = 513.4D0*U4/(ZVAR1*U) + ONE/T2
      T2 = 47.3D0*UM8+ONE/T2
      T2 = ONE/(4.D3+ONE/T2)
c
c       OPAC1 is the total opacity normalized by the electron pressure.
c
      OPAC1 = 4.819D-13*V/U+T1+T2
      RETURN
      END
      SUBROUTINE MODDMP(IIN,DTIME)
      IMPLICIT REAL*8(A-H,O-Z)
      integer*2 i
c
c          Dumps the extant model and the prior time step to the
c       printer file for diagnosis.
c
c                                                8/6/88 WD PESNELL
c
      parameter ( nmax=256 )
      common/phypar/ rn1(nmax),tn1(nmax),vn1(nmax),dm1(nmax),dm2(nmax),
     $               rm(nmax),pn1(nmax),un1(nmax),en(nmax)
      common/phypm1/ rp(nmax),tp(nmax),vp(nmax),p(nmax),up(nmax),
     $               ep(nmax)
      common/thrder/ chr(nmax),chit(nmax),q(nmax),cp(nmax),cv(nmax),
     $               g1(nmax),g3m1(nmax),delad(nmax),cs(nmax),
     $               akap(nmax),dkdr(nmax),dkdt(nmax)
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk1/   irad,nflag,npts
      common/observ/ teff,rlumgv,totmas,rphoto,corlum
      common/stlums/ rlumn1(nmax),rlumn(nmax), rlcon(nmax)
C
      DT = DTIME
      WRITE(11,2000) IIN,DT
      WRITE(11,2001) RLUMGV,TOTMAS,TEFF,RPHOTO,RN1(1),CORLUM
      WRITE(11,2004)
      WRITE(11,2005) (I,RN1(I+1),TN1(I),VN1(I),PN1(I),G1(I),RLUMN1(I+1),
     $   AKAP(I),DKDR(I),DKDT(I),CV(I),UN1(I+1),EN(I),I=1,NPTS)
      WRITE(11,2100)
      WRITE(11,2105) (I,RP(I+1),TP(I),VP(I),P(I),UP(I),EP(I),RLUMN(I+1),
     $   I=1,NPTS)
C
      RETURN
C
 2000 FORMAT(1H1,/,1X,19HOUTPUT FROM MODDMP.,/,
     $      1X,48H DIAGNOSTIC MODEL, CURRENT TIME THEN PRIOR STEP.,//,
     $      2X,11HBAD ZONE IS,I4,14H TIME STEP WAS,1PE11.4)
 2001 FORMAT(4X,18H TOTAL LUMINOSITY=,1PE11.4,12H TOTAL MASS=,E11.4,
     $ 17H EFFECTIVE TEMP.=,E11.4,/,6X,21H PHOTOSPHERIC RADIUS=,E11.4,
     $ 13H CORE RADIUS=,E11.4,11H CORE LUM.=,E11.4)
 2004 FORMAT(3X,1HI,3X,6HRADIUS,5X,4HTEMP,5X,6HSP.VOL,3X,8HPRESSURE,
     $ 5X,2HG1,7X,4HRLUM,5X,7HOPACITY,3X,6HDLKDLR,4X,6HDLKDLT,6X,
     $ 2HCV,7X,3HUN1,5X,8HINT.ENER)
 2005 FORMAT(1X,I4,1P,12E10.3)
 2006 FORMAT(1X,I4,1P,5E15.8)
 2100 FORMAT(1H1,/,1X,19HOUTPUT FROM MODDMP.,/,
     $      1X,44H DIAGNOSTIC MODEL, LAST CONVERGED TIME STEP.,//,
     $      3X,1HI,3X,6HRADIUS,5X,4HTEMP,5X,6HSP.VOL,3X,8HPRESSURE,
     $ 5X,2HUP,7X,2HEP,7X,4HRLUM)
 2105 FORMAT(1X,I4,1P,7E10.3)
C
      END
