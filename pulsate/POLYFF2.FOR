      program polyff
      implicit real*8(a-h,o-z)
c
c          Polytrope Integrator. Integrates the second order equation
c       governing polytropes using a Runge-Kutta method and calculates
c       the linear, adiabatic eigenvalues of the resulting model. The
c       adiabatic index is set to 5/3 in this version.
c
c          This version will calculate models with the gravity term
c       modified to include the Eotvos contribution. The Yukawa type
c       modification requires two input parameters.
c
c       03/09/95: Rederived the analytic solution, there was an error
c                 in ApJ, 344, 851-855, 1989. The corrected solution
c                 goes to 1 as the radius goes to 0 and gives the same
c                 external radius.
c
      parameter ( nmax=500 )
      COMMON/PHYPAR/ R(nmax),THETA(nmax),DTHDR(nmax),V(nmax),RM(nmax),
     $      GOR(nmax),XI(nmax),DM1(nmax),DM2(nmax),BV(nmax)
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      COMMON/BLK4/   P(nmax),G1(nmax),RHO(nmax),RZONE(nmax)
      COMMON/BLK8/   G,AC3,ACRAD,PI,TWOPI,FORPI,PI8,PI43
      COMMON/COREVL/ RZ0,P0,RHO0
      COMMON/STNMOD/ NMODE
      COMMON/EOTVOS/ DBET(nmax),BETA(nmax),FFALFA,FFLAMB
      DIMENSION C2(nmax),SL2(nmax)
C
      OPEN (1,FILE='SUMOUT.DAT',STATUS='UNKNOWN',access='append')
      OPEN (11,FILE='POLFF.OUT',STATUS='NEW')
      OPEN (12,FILE='POLFF.PLT',STATUS='NEW')
C
C          READ IN THE INPUT VARIABLES FROM FILE 5.
C
C          IOMOD := OUTPUT FLAG. IOUT GREATER THAN UNITY PRODUCES
C                 VOLUMINOUS LINES OF OUTPUT.
C
C          NPTS := NUMBER OF POINTS IN MODEL. NPTS MUST BE LESS THAN
C               NMAX-2, WHICH IS THE DIMENSION OF THE MATRICES
C               USED.
C
C          FNPOL := POLYTROPIC INDEX OF THE INITIAL MODEL. DEFAULT
C                  VALUE FOR FNPOL .LT. 0 IS 3.
C
C          R0 := INITIAL GUESS AT THE RADIUS OF THE MODEL. THIS
C               PARAMETER IS FOUND AS AN EIGENVALUE IN THE SOLUTION
C               AND A GOOD GUESS SPEEDS THE CONVERGENCE. DEFAULT
C               VALUE IS 5 (USED IF THE INPUTTED VALUE IS .LE. 0).
C
C          FFALFA := PARAMETER DESCRIBING THE STRENGTH OF THE EOTVOS
C                   FORCE. TYPICAL VALUE IS 0.01 TO 0.5, DEFAULT IS
C                   0, USED IF FFALFA < 0.0.
C
C          FFLAMB := PARAMETER DESCRIBING THE LENGTH SCALE OF THE
C                   EOTVOS FORCE. MUST BE GREATER THAN ZERO FOR THE
C                   PROGRAM TO CONTINUE. NO DEFAULT, IF FFLAMB <= 0,
C                   THE PROGRAM STOPS.
C
      IDONPL = 1
   1  CONTINUE
      WRITE(*,1000)
      READ(*,*,END=910) IOMOD,NPTS,FNPOL,R0,FFALFA,FFLAMB
      IF( NPTS .LE. 0 ) GOTO 910
      IF( R0 .LE. ZERO ) R0 = 5.D0
      IF( FFALFA .LT. ZERO ) FFALFA = ZERO
      IF( NPTS+2 .GT. NMAX ) NPTS = NMAX - 2
C
      CALL POLINT(IOMOD,NPTS,FNPOL,R0,IDONPL,FFALFA,FFLAMB)
      if( npts .le. 0 ) stop
      IF( FNPOL .LT. 0 ) FNPOL = 0.D0
      BETANP = ONE-(ONE+XI(NPTS)/FFLAMB)*FFALFA*DEXP(-XI(NPTS)/FFLAMB)
      RHOBAR = 1.D0/PI43
      RLAM = RHOBAR*(XI(NPTS)**3/(THRE*DTHDR(NPTS)))*BETANP
      ALFA = 1.D0/XI(NPTS)
      Pc = FORPI*G*(ALFA*RLAM)**2/(ONE+FNPOL)
      NPULS = NPTS-2
      DO 5 I=1,NPTS
         THETA(I) = AHF*(THETA(I)+THETA(I+1))
   5  CONTINUE
      do 10 i=1,npuls
         if( theta(i+1) .gt. 0 ) then
            rho(i) = rlam*theta(i+1)**fnpol
            v(i) = one/rho(i)
         else
            rho(i) = 0.d0
            v(i) = 0.d0
         endif
         p(i) = theta(i+1)*rho(i)*pc/rlam
         g1(i) = 5.d0/3.d0
c         c2(i) = p(i)*v(i)*g1(i)
  10  continue
      bvfac = fnpol - (fnpol+one)/g1(1)
      om02 = pi*g*rhobar
      DO 15 I=1,NPULS+1
         R(I) = ALFA*XI(I+1)
         zeta = xi(i+1)/fflamb
         beta(i)=one-(one+zeta)*ffalfa*dexp(-zeta)
         rm(i) = forpi*alfa**3*rlam*dthdr(i+1)/beta(i)
         gor(i) = g*forpi*alfa*rlam*dthdr(i+1)/xi(i+1)**2
         dbet(i) = gor(i)*ffalfa*fflamb*r(i)**2*dexp(-zeta)
     $                        /beta(i)
        if( i .ne. npuls + 1 .and. abs(bvfac) .le. 1.d-12 ) then
C            RAI = BVFAC*DTHDR(I+1)/(THETA(I)*XI(I+1))
C            BV(I) = DLOG10((DABS(RAI)*GOR(I)/R(I))/OM02)
         endif
  15  CONTINUE
      IGMAX = 0
      DO 16 I=1,NPULS
         IF( GOR(I+1) .LT. GOR(I) ) GOTO 17
          IGMAX = I+1
  16  CONTINUE
  17  CONTINUE
      WRITE(*,1700) IGMAX,R(IGMAX)/R(NPULS+1),RM(IGMAX)/RM(NPULS+1)
C      CALL PLTDMP(BV,1000,NPULS,'BV  ')
C      CALL PLTDMP(GOR,1000,NPULS,'G   ')
C      CALL PLTDMP(C2,1000,NPULS,'C2  ')
      RZ0 = AHF*R(1)
      RHO0 = RLAM
      P0 = PC
C      ASYMFR = ZERO
      DO 20 I=1,NPULS
         RZONE(I) = AHF*(R(I)+R(I+1))
C         ASYMFR = ASYMFR + (R(I+1)-R(I))/DSQRT(C2(I))
C         C2(I) = DLOG10(C2(I))
C         SL2(I) = DLOG10(6.D0/RZONE(I)**2) + C2(I)
  20  CONTINUE
C      ASYMFR = ONE/ASYMFR
C      WRITE(*,1099) ASYMFR
C      CALL PLTDMP(SL2,1000,NPULS,'SL2 ')
      R(NPULS+1) = ALFA*R0
      TMAS = PI43*R(1)**3*RLAM
      DO 30 I=1,NPULS
         DM1(I) = (R(I+1)**3-R(I)**3)*PI43*RHO(I)
         TMAS = TMAS + DM1(I)
  30  CONTINUE
      DM2(1) = AHF*(RM(1)+DM1(1))
      DO 40 I=2,NPULS
         DM2(I) = AHF*(DM1(I)+DM1(I-1))
  40  CONTINUE
      DM2(NPULS+1) = AHF*DM1(NPULS)
      RMCH = (RM(NPULS+1)-TMAS)/RM(NPULS+1)
      WRITE(11,2000)
      WRITE(11,2003) FNPOL,R0,ALFA,RLAM,RHOBAR,RM(NPULS+1),TMAS,RMCH
      WRITE(11,2006) RZ0,P0,RHO0
      WRITE(11,1700) IGMAX,R(IGMAX)/R(NPULS+1),RM(IGMAX)/RM(NPULS+1)
      WRITE(11,1095) ASYMFR,FFALFA,FFLAMB
      WRITE(11,2004)
      WRITE(11,2005) (I,XI(I),THETA(I),DTHDR(I),R(I),P(I),RHO(I),
     $   GOR(I),RM(I),DM1(I),RZONE(I),I=1,NPTS)
C
      QCH = RLAM/RHOBAR
      WRITE(1,2007) NPULS,FNPOL,R0,ALFA,RLAM,RHOBAR,RM(NPULS+1),
     $      TMAS,RMCH
      WRITE(1,2008) RZ0,P0,RHO0,QCH
      WRITE(1,1700) IGMAX,R(IGMAX)/R(NPULS+1),RM(IGMAX)/RM(NPULS+1)
      WRITE(1,1095) ASYMFR,FFALFA,FFLAMB
c
      goto 1
 910  continue
      stop
c
 1000 format(1x,'Enter iout, no. zones, fnpol, r0 and',/,
     $       3x,'Eotvos parameters (alpha and lambda): ',$)
 1010 format(1x,27hEnter pulsation parameters:,/,
     $  1x,38h io, l, nom, ihmax, omlow, and omhigh.)
 1700 format(1x,15hmax. g in zone ,i4,9h radius =,1pe10.3,
     $      8h mass = ,e10.3)
 1095 format(1x,36hAsymptotic frequency for p modes is ,1pe11.4,/,
     $   3x,22hEotvos parameters are ,e11.3,0pf7.3)
 1099 format(1x,36hAsymptotic frequency for p modes is ,1pe11.4)
 2000 format(41h1 polytropic model in hydrostatic balance,//,
     $      1x,i5,6h zones,/)
 2003 format(7h fnpol=,f6.2,26h guess for surface radius=,1pe10.3,
     $   7h alpha=,e10.3,/,2x,8h lambda=,e10.3,8h rhobar=,e10.3,
     $   8h totmas=,e10.3,11h int. mass=,e10.3,7h diff.=,e10.3)
 2004 format(3x,1hi,5x,2hxi,7x,5htheta,5x,5hdthdr,4x,6hradius,6x,1hp,
     $   8x,3hrho,8x,1hg,9x,2hrm,7x,3hdm1,7x,5hrzone)
 2005 format(1x,i4,1p,10e10.3)
 2006 format(1x,5h rz0=,1pe10.3,4h p0=,e10.3,6h rho0=,e10.3)
 2007 format(1h1,/,1x,18houtput from polyrk,i5,6h zones,/,
     $   7h fnpol=,f6.2,26h guess for surface radius=,1pe10.3,
     $   7h alpha=,e10.3,/,2x,8h lambda=,e10.3,8h rhobar=,e10.3,/,
     $   8h totmas=,e10.3,11h int. mass=,e10.3,7h diff.=,e10.3)
 2008 format(1x,5h rz0=,1pe10.3,4h p0=,e10.3,6h rho0=,e10.3,/
     $   1x,29h central/average densities = ,e10.3)
      end
      SUBROUTINE POLINT(IOUT,NIN,FNPOL,R0,IDONPL,FFALFA,FFLAMB)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LOUT
C
      parameter ( nmax=500 )
      COMMON/PHYPAR/ R(nmax),THETA(nmax),DTHDR(nmax),V(nmax),RM(nmax),
     $      GOR(nmax),XI(nmax),DM1(nmax),DM2(nmax),BV(nmax)
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      DIMENSION X(2),F(2),Y(2)
      DATA ACCUR/1.D-4/
C
      N = NIN
      LOUT = .FALSE.
      IF( IOUT .GE. 1 ) LOUT = .TRUE.
C
      IF( FNPOL .LT. 0 ) then
C
C          ANALYTIC SOLUTION FOR N=0, REACHED BY ENTERING N<0.
C
         R0 = DSQRT(6.D0)
         DO 65 IVERG=1,50
c
c          3/4/95
c
c         FR0 = ONE - R0**2/6.D0 - FFALFA*
c     $      DEXP(-R0/FFLAMB)*(R0**2/THRE + (ONE+R0/FFLAMB)*FFLAMB**2)
c         DFR0 = R0*(ONE+FFALFA*DEXP(-R0/FFLAMB)*(ONE-R0/FFLAMB))/THRE
            zeta = R0/fflamb
            FR0 = ONE - R0**2/6.D0 + FFALFA*fflamb**2*
     $        (one - dexp(-zeta)*(1.d0 + zeta*(one+zeta/3.d0)) )
            DFR0 = R0*(ONE-FFALFA*DEXP(-zeta)*(ONE+zeta))/THRE
            DR0 = FR0/DFR0
            R0 = DR0 + R0
            WRITE(6,6501) IVERG,R0,DR0,FR0,DFR0
            IF( DABS(DR0/R0) .LT. ACCUR ) then
               N = N+1
               DR = R0/FLOAT(N+1)
               DO 67 I=1,N+2
                  XI(I) = FLOAT(I-1)*DR
                  zeta = xi(i)/fflamb
                  theta(i) = one - xi(i)**2/6.d0 + ffalfa*fflamb**2*
     $                (one - dexp(-zeta)*(1.d0 + zeta*(1.d0+zeta/thre)))
                  dthdr(i) = xi(i)**3*(one-ffalfa*dexp(-zeta)*
     $               (one+zeta))/thre
  67           CONTINUE
               NIN = N+1
               N = N+1
               call pltint(n,fnpol,xi(np),idonpl)
               call pltdmp(theta,nmax,n,'thet')
               call pltdmp(dthdr,nmax,n,'dtdr')
               return
            endif
  65     CONTINUE
         WRITE(6,6500) R0
         nin = -1
      else
c
         ITRY =-1
         ITERMX = 20
         DO 40 ITER=1,ITERMX
C
            DR = R0/FLOAT(N+1)
C
            NP = 0
            XI(1) = ZERO
            CENVAL = ONE 
            THETA(1) = CENVAL
            DTHDR(1) = ZERO
            XI(2) = DR
            RFAC = ONE - FFALFA
            zeta = DR/FFLAMB
         theta(2) = cenval-dr**2*((rfac-fnpol*(rfac*dr)**2/20.d0)/6.d0
     $      + 23.d0*ffalfa*(zeta)**2/120.d0)
            dthdr(2) = dr**3*((rfac-fnpol*(rfac*dr)**2/10.d0)/3.d0
     $      + 23.d0*ffalfa*(zeta)**2/30.d0)
            x(1) = dr
            y(1) = theta(2)
            y(2) =-dthdr(2)/dr**2
            do 30 ic=2,n+1
               do 35 m=1,5
                  kst = irunge(2,y,f,x(1),dr,m)
                  if( kst .eq. 0 ) goto 31
                  if( y(1) .lt. zero ) goto 31
                     zeta = x(1)/fflamb
                     ffbeta = one-(one+zeta)*ffalfa*dexp(-zeta)
                     dffbet = ffalfa*zeta/fflamb*dexp(-zeta)
                     f(1) = y(2)
                     f(2) =-two*y(2)/x(1) - ffbeta*y(1)**fnpol +
     $                 y(2)*dffbet/ffbeta
  35           continue
  31          continue
               NP = IC + 1
               XI(NP) = X(1)
               THETA(NP) = Y(1)
               DTHDR(NP) =-X(1)*X(1)*Y(2)
               IF( Y(1) .LT. ZERO ) GOTO 45
  30        CONTINUE
            WRITE(11,4500) FNPOL,R0
            IF( ITRY .LE. 0 ) then
               R0 = R0*1.25E0
               WRITE(*,4700) NP,FNPOL,R0,R0
C
C          CONVERGED TO A SOLUTION, PUT THE SURFACE
C       AT R(N).
C
            endif
            goto 40
  45     CONTINUE
               IF( NP .LE. 1 ) then
C
C          UNKNOWN PROBLEM IN THE INTEGRATION. RUN AWAY, RUN AWAY...
C
                  WRITE(*,7000) K,NP,X(1),Y(1),Y(2),F(1),F(2)
                  nin =-1
                  return
               endif
               ITRY = 1
               DR0 = THETA(NP)*XI(NP)*XI(NP)/DTHDR(NP)
               R0 = XI(NP) + DR0
               IF( DABS(DR0/R0) .LT. ACCUR ) GOTO 50
               WRITE(*,4700) NP,FNPOL,R0,DR0
               IF( LOUT ) WRITE(11,5501) (I,XI(I),THETA(I),DTHDR(I),
     $                   I=1,NP)
  40     CONTINUE
C
C          NOT CONVERGED AFTER ITERMX TRIES, WRITE ERROR MESSAGE AND
C       STOP THE CODE.
C
         WRITE(11,4000) ITERMX,FNPOL,R0
         WRITE(*,4000) ITERMX,FNPOL,R0
         nin =-1
         return
      endif
  50  CONTINUE
C
C          CONVERGED TO A SOLUTION, SET NUMBER OF ZONES AND OUTER
C       BOUNDARY VALUE OF THETA.
C
      theta(np) = zero
      write(11,5500) fnpol,xi(np)
      if( lout ) write(11,5501) (i,xi(i),theta(i),dthdr(i),i=1,np)
      n = np
      nin = np
      call pltint(n,fnpol,xi(np),idonpl)
      call pltdmp(theta,nmax,n,'thet')
      call pltdmp(dthdr,nmax,n,'dtdr')
      return
C
 4000 FORMAT(1X,20HNO CONVERGENCE AFTER,I4,18H TRIES FOR FNPOL =,F7.3,
     $      9H RADIUS =,1PE11.4)
 4500 FORMAT(1X,22HNO CONVERGENCE FNPOL =,F7.3,9H RADIUS =,1PE11.4)
 4700 FORMAT(1X,4HNP =,I4,6H FNPOL,F7.3,8H RADIUS=,1PE11.3,5H DR0=,
     $    E11.3)
 5500 FORMAT(/,5X,5HFNPOL,F7.3,8H RADIUS=,1PE15.6)
 5501 FORMAT(1X,I5,1P,3E15.6)
 6500 FORMAT(1X,43HNO CONVERGENCE IN R0 (ANALYTIC MODEL), R0 =,1PE11.4)
 6501 FORMAT(1X,I3,1P,4E11.4)
 7000 FORMAT(1X,2I3,1P,5E11.3)
C
      END
      FUNCTION IRUNGE(N,Y,F,X,H,M)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PHI(20),SAVEY(20),Y(N),F(N)
C
C          STEP 1.
C
      if( m .eq. 1 ) then
         IRUNGE = 1
C
C          STEP 2.
C
      elseif( m .eq. 2 ) then
         DO 25 J=1,N
            SAVEY(J) = Y(J)
            PHI(J) = F(J)
            Y(J) = Y(J) + H*F(J)/2.D0
  25     CONTINUE
         X = X + H/2.D0
         IRUNGE = 1
C
C         STEP 3.
C
      elseif( m .eq. 3 ) then
         DO 35 J=1,N
            PHI(J) = PHI(J) + 2.D0*F(J)
            Y(J) = SAVEY(J) + F(J)*H/2.D0
  35     CONTINUE
         IRUNGE = 1
C
C          STEP 4.
C
      elseif( m .eq. 4 ) then
         DO 45 J=1,N
            PHI(J) = PHI(J) + 2.D0*F(J)
            Y(J) = SAVEY(J) + F(J)*H
  45     CONTINUE
         X = X + H/2.D0
         IRUNGE = 1
C
C         FINAL PASS.
C
      elseif( m .eq. 5 ) then
         DO 55 J=1,N
            Y(J) = SAVEY(J) + (PHI(J) + F(J))*H/6.D0
  55     CONTINUE
         IRUNGE = 0
      endif
      return
      end
      block data
      implicit real*8(a-h,o-z)
c
c       Fundamental constants are from Novotny, Introduction to
c    Stellar Atmospheres and Interiors, 1973, Appendix II.
c                                                 2/16/83  WDP
c
      common/const/  zero,one,two,thre,for,ten,ahf,qrt
      common/blk8/   grav,ac3,sigma,pi,pi2,pi4,pi8,pi43
c      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/6.6726d-8,7.5595d-5,
      data grav,ac3,sigma,pi,pi2,pi4,pi8,pi43/ 1.d0, 7.5595d-5,
     $ 5.66961d-5, 3.1415926536d0, 6.2831853072d0, 12.5663706144d0,
     $ 25.1327412288d0, 4.1887902048d0 /
      data zero,one,two,thre,for,ten,ahf,qrt /0.d0, 1.d0, 2.d0, 3.d0,
     $       4.d0, 10.d0, 0.5D0, 0.25d0 /
      end
      subroutine pltint(n,fnpol,r0,idonpl)
      implicit real*8(a-h,o-z)
      character*4 cxi,cx
C
C          INITIALIZE THE PLOT FILE AND WRITE THE TITLE LINE.
C
      parameter ( nmax=500 )
      COMMON/PHYPAR/ R(nmax),THETA(nmax),DTHDR(nmax),V(nmax),RM(nmax),
     $      GOR(nmax),XI(nmax),DM1(nmax),DM2(nmax),BV(nmax)
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      COMMON/BLK4/   P(nmax),G1(nmax),RHO(nmax),RZONE(nmax)
c
      data cxi/'xi  '/,cx/'x   '/
c
      if( idonpl .gt. 0 ) then
         idonpl =-1
         if( fnpol .lt. zero ) then
            write(12,1000) zero,r0
         else
            write(12,1000) fnpol,r0
         endif
      endif
      call pltdmp(xi(2),nmax,n,cxi)
      do  20 i=2,n
         r(i) = xi(i)/xi(n)
  20  continue
      call pltdmp(r(2),nmax,n-1,cx   )
      return
 1000 FORMAT(1X,'Eotvos sphere with n =',F5.2,15H outer radius =,F7.4)
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
