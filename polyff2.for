      PROGRAM POLYFF
      IMPLICIT REAL*8(A-H,O-Z)
C
C          POLYTROPE INTEGRATOR. INTEGRATES THE SECOND ORDER EQUATION
C       GOVERNING POLYTROPES USING A RUNGE-KUTTA METHOD AND CALCULATES
C       THE LINEAR, ADIABATIC EIGENVALUES OF THE RESULTING MODEL. THE
C       ADIABATIC INDEX IS SET TO 5/3 IN THIS VERSION.
C
C
C          THIS VERSION WILL CALCULATE MODELS WITH THE GRAVITY TERM
C       MODIFIED TO INCLUDE THE EOTVOS CONTRIBUTION. THE YUKAWA TYPE
C      MODIFICATION REQUIRES THREE INPUT PARAMETERS.
C
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
      CALL INTDAT
      OPEN (1,FILE='SUMOUT.DAT',STATUS='UNKNOWN',access='append')
C      OPEN (UNIT=5,FILE='SYS$INPUT',STATUS='UNKNOWN')
C      OPEN (UNIT=6,FILE='SYS$OUTPUT',STATUS='UNKNOWN')
      OPEN (11,FILE='POLOUT.DAT',STATUS='UNKNOWN')
      OPEN (12,FILE='POLPLT.DAT',STATUS='UNKNOWN')
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
C      IF( FNPOL .LT. ZERO ) FNPOL = THRE
      IF( FFALFA .LT. ZERO ) FFALFA = ZERO
      IF( NPTS+2 .GT. NMAX ) NPTS = NMAX - 2
C
      CALL POLINT(IOMOD,NPTS,FNPOL,R0,IDONPL,FFALFA,FFLAMB)
      IF( FNPOL .LT. 0 ) FNPOL = 0.0D0
      BETANP = ONE-(ONE+XI(NPTS)/FFLAMB)*FFALFA*DEXP(-XI(NPTS)/FFLAMB)
      RHOBAR = 1.D0/PI43
      RLAM = RHOBAR*(XI(NPTS)**3/(THRE*DTHDR(NPTS)))*BETANP
      ALFA = 1.D0/XI(NPTS)
      PC = FORPI*G*(ALFA*RLAM)**2/(ONE+FNPOL)
      NPULS = NPTS-2
      DO 5 I=1,NPTS
         THETA(I) = AHF*(THETA(I)+THETA(I+1))
   5  CONTINUE
      DO 10 I=1,NPULS
         RHO(I) = RLAM*THETA(I+1)**FNPOL
         V(I) = ONE/RHO(I)
         P(I) = THETA(I+1)*RHO(I)*PC/RLAM
         G1(I) = 5.D0/3.D0
C         C2(I) = P(I)*V(I)*G1(I)
  10  CONTINUE
      BVFAC = FNPOL - (FNPOL+ONE)/G1(1)
      OM02 = PI*G*RHOBAR
      DO 15 I=1,NPULS+1
         R(I) = ALFA*XI(I+1)
         BETA(I)=ONE-(ONE+XI(I+1)/FFLAMB)*FFALFA*DEXP(-XI(I+1)/FFLAMB)
         RM(I) = FORPI*ALFA**3*RLAM*DTHDR(I+1)/BETA(I)
         GOR(I) = G*FORPI*ALFA*RLAM*DTHDR(I+1)/XI(I+1)**2
         DBET(I) = GOR(I)*FFALFA*FFLAMB*R(I)**2*DEXP(-XI(I+1)/FFLAMB)
     $                        /BETA(I)
        IF( I .EQ. NPULS + 1 ) GOTO 15
        IF( DABS(BVFAC) .LE. 1.D-12 ) GOTO 15
C         RAI = BVFAC*DTHDR(I+1)/(THETA(I)*XI(I+1))
C         BV(I) = DLOG10((DABS(RAI)*GOR(I)/R(I))/OM02)
  15  CONTINUE
      IGMAX = 0
      DO 16 I=1,NPULS
         IF( GOR(I+1) .LT. GOR(I) ) GOTO 17
          IGMAX = I+1
  16  CONTINUE
  17  CONTINUE
      WRITE(*,1700) IGMAX,R(IGMAX)/R(NPULS+1),RM(IGMAX)/RM(NPULS+1)
C      CALL PLTDMP(BV,1000,NPULS,4HBV  )
C      CALL PLTDMP(GOR,1000,NPULS,4HG   )
C      CALL PLTDMP(C2,1000,NPULS,4HC2  )
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
C      CALL PLTDMP(SL2,1000,NPULS,4HSL2 )
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
C
      GOTO 1
 910  CONTINUE
      STOP
C
 1000 FORMAT(1X,36HENTER IOUT, NO. ZONES, FNPOL, R0 AND,/,
     $       3X,35HEOTVOS PARAMETERS (ALPHA AND BETA).)
 1010 FORMAT(1X,27HENTER PULSATION PARAMETERS:,/,
     $  1X,38H IO, L, NOM, IHMAX, OMLOW, AND OMHIGH.)
 1700 FORMAT(1X,15HMAX. G IN ZONE ,I4,9H RADIUS =,1PE10.3,
     $      8H MASS = ,E10.3)
 1095 FORMAT(1X,36HASYMPTOTIC FREQUENCY FOR P MODES IS ,1PE11.4,/,
     $   3X,22HEOTVOS PARAMETERS ARE ,E11.3,0PF7.3)
 1099 FORMAT(1X,36HASYMPTOTIC FREQUENCY FOR p MODES IS ,1PE11.4)
 2000 FORMAT(41H1 POLYTROPIC MODEL IN HYDROSTATIC BALANCE,//,
     $      1X,I5,6H ZONES,/)
 2003 FORMAT(7H FNPOL=,F6.2,26H GUESS FOR SURFACE RADIUS=,1PE10.3,
     $   7H ALPHA=,E10.3,/,2X,8H LAMBDA=,E10.3,8H RHOBAR=,E10.3,
     $   8H TOTMAS=,E10.3,11H INT. MASS=,E10.3,7H DIFF.=,E10.3)
 2004 FORMAT(3X,1HI,5X,2HXI,7X,5HTHETA,5X,5HDTHDR,4X,6HRADIUS,6X,1HP,
     $   8X,3HRHO,8X,1HG,9X,2HRM,7X,3HDM1,7X,5HRZONE)
 2005 FORMAT(1X,I4,1P,10E10.3)
 2006 FORMAT(1X,5H RZ0=,1PE10.3,4H P0=,E10.3,6H RHO0=,E10.3)
 2007 FORMAT(1H1,/,1X,18HOUTPUT FROM POLYRK,I5,6H ZONES,/,
     $   7H FNPOL=,F6.2,26H GUESS FOR SURFACE RADIUS=,1PE10.3,
     $   7H ALPHA=,E10.3,/,2X,8H LAMBDA=,E10.3,8H RHOBAR=,E10.3,/,
     $   8H TOTMAS=,E10.3,11H INT. MASS=,E10.3,7H DIFF.=,E10.3)
 2008 FORMAT(1X,5H RZ0=,1PE10.3,4H P0=,E10.3,6H RHO0=,E10.3,/
     $   1X,29H CENTRAL/AVERAGE DENSITIES = ,E10.3)
      END
      SUBROUTINE POLINT(IOUT,NIN,FNPOL,R0,IDONPL,FFALFA,FFLAMB)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LOUT
C
      parameter ( nmax=500 )
      COMMON/PHYPAR/ R(nmax),THETA(nmax),DTHDR(nmax),V(nmax),RM(nmax),
     $      GOR(nmax),XI(nmax),DM1(nmax),DM2(nmax),BV(nmax)
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      DIMENSION X(2),F(2),Y(2)
      DATA ACCUR/1.D-6/
C
      N = NIN
      LOUT = .FALSE.
      IF( IOUT .GE. 1 ) LOUT = .TRUE.
C
      IF( FNPOL .LT. 0 ) GOTO 60
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
         THETA(2) = CENVAL-DR**2*((RFAC-FNPOL*(RFAC*DR)**2/20.D0)/6.D0
     $      + 23.D0*FFALFA*(DR/FFLAMB)**2/120.D0)
         DTHDR(2) = DR**3*((RFAC-FNPOL*(RFAC*DR)**2/10.D0)/3.D0
     $      + 23.D0*FFALFA*(DR/FFLAMB)**2/30.D0)
         X(1) = DR
         Y(1) = THETA(2)
         Y(2) =-DTHDR(2)/DR**2
         DO 30 IC=2,N+1
            DO 35 M=1,5
               KST = IRUNGE(2,Y,F,X(1),DR,M)
               IF( KST .EQ. 0 ) GOTO 31
               IF( Y(1) .LT. ZERO ) GOTO 31
               FFBETA = ONE-(ONE+X(1)/FFLAMB)*FFALFA*DEXP(-X(1)/FFLAMB)
               DFFBET = FFALFA*X(1)/FFLAMB/FFLAMB*DEXP(-X(1)/FFLAMB)
                F(1) = Y(2)
                F(2) =-TWO*Y(2)/X(1) - FFBETA*Y(1)**FNPOL +
     $                 Y(2)*DFFBET/FFBETA
  35        CONTINUE
  31       CONTINUE
            NP = IC + 1
            XI(NP) = X(1)
            THETA(NP) = Y(1)
            DTHDR(NP) =-X(1)*X(1)*Y(2)
            IF( Y(1) .LT. ZERO ) GOTO 45
  30     CONTINUE
         WRITE(11,4500) FNPOL,R0
         IF( ITRY .LE. 0 ) GOTO 49
C
C          CONVERGED TO A SOLUTION, PUT THE SURFACE
C       AT R(N).
C
  45     CONTINUE
         IF( NP .LE. 1 ) GOTO 70
         ITRY = 1
         DR0 = THETA(NP)*XI(NP)*XI(NP)/DTHDR(NP)
         R0 = XI(NP) + DR0
         IF( DABS(DR0/R0) .LT. ACCUR ) GOTO 50
         WRITE(*,4700) NP,FNPOL,R0,DR0
         IF( LOUT ) WRITE(11,5501) (I,XI(I),THETA(I),DTHDR(I),I=1,NP)
         GOTO 40
  49    CONTINUE
         R0 = R0*1.25E0
         WRITE(*,4700) NP,FNPOL,R0,R0
  40  CONTINUE
C
C          NOT CONVERGED AFTER ITERMX TRIES, WRITE ERROR MESSAGE AND
C       STOP THE CODE.
C
      WRITE(11,4000) ITERMX,FNPOL,R0
      WRITE(*,4000) ITERMX,FNPOL,R0
      STOP
  50  CONTINUE
C
C          CONVERGED TO A SOLUTION, SET NUMBER OF ZONES AND OUTER
C       BOUNDARY VALUE OF THETA.
C
      THETA(NP) = ZERO
      WRITE(11,5500) FNPOL,XI(NP)
      IF( LOUT ) WRITE(11,5501) (I,XI(I),THETA(I),DTHDR(I),I=1,NP)
      N = NP
      NIN = NP
      CALL PLTINT(N,FNPOL,XI(NP),IDONPL)
C      CALL PLTDMP(THETA,NMAX,N,4HTHET)
C      CALL PLTDMP(DTHDR,NMAX,N,4HDTDR)
      RETURN
C
C          ANALYTIC SOLUTION FOR N=0, REACHED BY ENTERING N<0.
C
  60  CONTINUE
      R0 = DSQRT(6.D0)
      DO 65 IVERG=1,50
c         FR0 = ONE - R0**2/6.D0 - FFALFA*
c     $      DEXP(-R0/FFLAMB)*(R0**2/THRE + (ONE+R0/FFLAMB)*FFLAMB**2)
c         DFR0 = R0*(ONE+FFALFA*DEXP(-R0/FFLAMB)*(ONE-R0/FFLAMB))/THRE
         zeta = R0/fflamb
         FR0 = ONE - R0**2/6.D0 - FFALFA*DEXP(-zeta)*
     $      fflamb**2*(1.d0 + zeta*(ONE+zeta))/3.d0
         DFR0 = R0*(ONE+FFALFA*DEXP(-zeta)*(ONE-zeta))/THRE
         DR0 = FR0/DFR0
         R0 = DR0 + R0
         WRITE(*,6501) IVERG,R0,DR0,FR0,DFR0
         IF( DABS(DR0/R0) .LT. ACCUR ) GOTO 66
  65  CONTINUE
 6501 FORMAT(1X,I3,1P,4E11.4)
      WRITE(*,6500) R0
      STOP
  66  CONTINUE
      N = N+1
      DR = R0/FLOAT(N+1)
      DO 67 I=1,N+2
         XI(I) = FLOAT(I-1)*DR
         THETA(I) = ONE - XI(I)**2/6.D0 -
     $      FFALFA*DEXP(-XI(I)/FFLAMB)*(XI(I)**2/THRE +
     $      (ONE+XI(I)/FFLAMB)*FFLAMB**2)
         DTHDR(I) = XI(I)**3*(ONE+FFALFA*DEXP(-XI(I)/FFLAMB)*
     $      (ONE-XI(I)/FFLAMB))/THRE
  67  CONTINUE
      NIN = N+1
      N = N+1
      RETURN
  70  CONTINUE
C
C          UNKNOWN PROBLEM IN THE INTEGRATION. RUN AWAY, RUN AWAY...
C
      WRITE(*,7000) K,NP,X(1),Y(1),Y(2),F(1),F(2)
      STOP
C
 4000 FORMAT(1X,20HNO CONVERGENCE AFTER,I4,18H TRIES FOR FNPOL =,F7.3,
     $      9H RADIUS =,1PE11.4)
 4500 FORMAT(1X,22HNO CONVERGENCE FNPOL =,F7.3,9H RADIUS =,1PE11.4)
 4700 FORMAT(1X,4HNP =,I4,6H FNPOL,F7.3,8H RADIUS=,1PE11.3,5H DR0=,
     $    E11.3)
 5500 FORMAT(/,5X,5HFNPOL,F7.3,8H RADIUS=,1PE15.6)
 5501 FORMAT(1X,I5,1P,3E15.6)
 6500 FORMAT(1X,43HNO CONVERGENCE IN R0 (ANALYTIC MODEL), R0 =,1PE11.4)
 7000 FORMAT(1X,2I3,1P,5E11.3)
C
      END
      FUNCTION IRUNGE(N,Y,F,X,H,M)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PHI(20),SAVEY(20),Y(N),F(N)
C
      GOTO (10,20,30,40,50),M
C
C          STEP 1.
C
  10  CONTINUE
      IRUNGE = 1
      RETURN
C
C          STEP 2.
C
  20  CONTINUE
      DO 25 J=1,N
         SAVEY(J) = Y(J)
         PHI(J) = F(J)
         Y(J) = Y(J) + H*F(J)/2.D0
  25  CONTINUE
      X = X + H/2.D0
      IRUNGE = 1
      RETURN
C
C         STEP 3.
C
  30  CONTINUE
      DO 35 J=1,N
         PHI(J) = PHI(J) + 2.D0*F(J)
         Y(J) = SAVEY(J) + F(J)*H/2.D0
  35  CONTINUE
      IRUNGE = 1
      RETURN
C
C          STEP 4.
C
  40  CONTINUE
      DO 45 J=1,N
         PHI(J) = PHI(J) + 2.D0*F(J)
         Y(J) = SAVEY(J) + F(J)*H
  45  CONTINUE
      X = X + H/2.D0
      IRUNGE = 1
      RETURN
C
C         FINAL PASS.
C
  50  CONTINUE
      DO 55 J=1,N
         Y(J) = SAVEY(J) + (PHI(J) + F(J))*H/6.D0
  55  CONTINUE
      IRUNGE = 0
      RETURN
      END
      SUBROUTINE INTDAT
      IMPLICIT REAL*8(A-H,O-Z)
C
C       FUNDAMENTAL CONSTANTS ARE FROM NOVOTNY, INTRODUCTION TO
C    STELLAR ATMOSPHERES AND INTERIORS, 1973, APPENDIX II.
C                                                 2/16/83  WDP
C
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      COMMON/BLK8/   GRAV,AC3,SIGMA,PI,PI2,PI4,PI8,PI43
C      DATA GRAV,AC3,SIGMA,PI,PI2,PI4,PI8,PI43/6.6726D-8,7.5595D-5,
      GRAV = 1.D0
      AC3 = 7.5595D-5
      SIGMA = 5.66961D-5
      PI = 3.1415926536D0
      PI2 = 6.2831853072D0
      PI4 = 12.5663706144D0
      PI8 = 25.1327412288D0
      PI43 = 4.1887902048D0
      ZERO = 0.0D0
      ONE = 1.D0
      TWO = 2.D0
      THRE = 3.D0
      FOR = 4.D0
      TEN = 10.D0
      AHF = 0.5
      QRT = 0.25D0
      END
      SUBROUTINE PLTINT(N,FNPOL,R0,IDONPL)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*4 CXI,CX
C
C          INITIALIZE THE PLOT FILE AND WRITE THE TITLE LINE.
C
      parameter ( nmax=500 )
      COMMON/PHYPAR/ R(nmax),THETA(nmax),DTHDR(nmax),V(nmax),RM(nmax),
     $      GOR(nmax),XI(nmax),DM1(nmax),DM2(nmax),BV(nmax)
      COMMON/CONST/  ZERO,ONE,TWO,THRE,FOR,TEN,AHF,QRT
      COMMON/BLK4/   P(nmax),G1(nmax),RHO(nmax),RZONE(nmax)
C
      DATA CXI/'XI  '/,CX/'X   '/
C
      IF( IDONPL .LE. 0 ) GOTO 10
      IDONPL =-1
      WRITE(12,1000) FNPOL,R0
  10  CONTINUE
      CALL PLTDMP(XI(2),NMAX,N,CXI)
      DO  20 I=2,N
         R(I) = XI(I)/XI(N)
  20  CONTINUE
      CALL PLTDMP(R(2),NMAX,N-1,CX   )
      RETURN
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