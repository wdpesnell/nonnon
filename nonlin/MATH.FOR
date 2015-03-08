      SUBROUTINE PWDERIV(NP,X,Y,D)
      DIMENSION X(1),Y(1),D(1)
c
c          Returns the derivative of Y with respect to X by
c       simple subtraction. The vector X is assumed to be 
c       monotonically increasing with index.
c
      DO I=1,NP-1
         D(I) = (Y(I+1)-Y(I))/(X(I+1)-X(I))
      ENDDO
C
      RETURN
      END
      FUNCTION FLINTER(X,X1,Y1,X2,Y2)
c
c          Interpolate a linear value between the points (X1, Y1)
c       and (X2, Y2).
c
      XDIF=X2-X1
      IF(XDIF.EQ.0.) THEN
         FLINTER=Y1
      ELSE
         SLOP=(Y2-Y1)/XDIF
         FLINTER=SLOP*(X-X1) + Y1
      ENDIF
C
      RETURN
      END
      FUNCTION EXPINTER(X,X1,Y1,X2,Y2)
c
c           Interpolate a function between (X1, Y1) and (X2, Y2) that
c       has a logarithmic dependence.
c
      XDIF=X2-X1
      IF( XDIF .EQ. 0. ) THEN
         EXPINTER=Y1
      ELSEIF( Y1 .LE. 0. ) THEN
         EXPINTER=Y2
      ELSEIF( Y2 .LE. 0. ) THEN
         EXPINTER=Y1
      ELSE
         SLOP=LOG(Y2/Y1)/(X2-X1)
         EXPINTER = Y1*EXP(SLOP*(X-X1))
      ENDIF
C
      RETURN
      END
C***************************

	SUBROUTINE SURFINTERP( NX, X, IMAX, XI, LX, 
	1		       NY, Y, JMAX, YJ, LY,
	2		       NDIM, FXY, IDIM, ZF, FINTERP)

	DIMENSION FXY(NDIM,1), ZF(IDIM,1)
	1		,X(1), XI(1), LX(1)
	2		,Y(1), YJ(1), LY(1)

c   ZF(I,J) is the result of linear interpolation at the points
c           XI(I=1:IMAX), YJ(J=1:JMAX)
c   of the input surface defined by matrix
c      FXY(N,M)=suface( X(N=1:NX), Y(M=1:NY) )

c   LX,LY are the location indices of arrays XI,YJ in arrays X,Y

	DO J=1,JMAX

		YY = YJ (J)		
		LJ2 = LY (J)
		IF(LJ2.GT.NY) THEN
			LJ2=NY
			LJ1=NY
		ELSE IF(LJ2.LE.1) THEN
			LJ2=1
			LJ1=1
		ELSE
			LJ1 = LJ2-1
		END IF
		Y1 = Y (LJ1)
		Y2 = Y (LJ2)

		DO I=1,IMAX

			LI2 = LX (I)
			IF(LI2.GT.NX) THEN
				LI2=NX
				LI1=NX
			ELSE IF(LI2.LE.1) THEN
				LI2=1
				LI1=1
			ELSE
				LI1 = LI2-1
			END IF

			F1 = FINTERP( YY, Y1, FXY (LI1,LJ1), 
	1				Y2, FXY (LI1,LJ2) )

			F2 = FINTERP( YY, Y1, FXY (LI2,LJ1), 
	1				Y2, FXY (LI2,LJ2) )

			ZF (I,J) = 
	1		FINTERP( XI(I), X(LI1), F1, X(LI2), F2 )

		END DO

	END DO

	RETURN
	END
C***************************

      SUBROUTINE INTERP(X1,Y1,X2,Y2,INVERS,X,Y)
c
c          Interpolate a linear function between (X1, Y1) and (X2, Y2).
c       If INVERS .lt. 0, return the value of X that would give Y.
c
      SLOP=(Y2-Y1)/(X2-X1)
      IF(INVERS.GE.0) THEN
         Y=SLOP*(X-X1)+Y1
      ELSE
         X=X1+(Y-Y1)/SLOP
      ENDIF
C
      RETURN
      END
      FUNCTION VNORM2 (N,U,V)
      DIMENSION U(1),V(1)
c
c          Unnormed difference of two vectors
c
      SUM=0.
      DO I=1,N
         DIFF=U(I)-V(I)
         SUM = SUM+DIFF*DIFF
      ENDDO
      VNORM2=SUM
C
      RETURN
      END
      FUNCTION FNORM2(F,G,X,N)
      DIMENSION F(1),G(1),X(1),D(500)
C
c          Intergrated vector differnece. Uses a trapezoidal rule
c       for the integral instead of simple sum.
c
      DO I=1,N
         DIFF=F(I)-G(I)
         D(I) = DIFF*DIFF
      ENDDO
      FNORM2=TRAPEZ(D,X,N)
C
      RETURN
      END
      SUBROUTINE POLCART(NP,RADIUS,ANGLE,X,Y)
      DIMENSION RADIUS(1),ANGLE(1),X(1),Y(1)
c
      DATA TORAD/.0174533/
c
      IF(NP.LT.0) THEN
         N=-NP
         DO I=1,N
            R=RADIUS(I)
            A=ANGLE(I)
            X(I) = R*COS(TORAD*A)
            Y(I)=R*SIN(TORAD*A)
         ENDDO
      ELSE
         DO I=1,N
            R=RADIUS(I)
            A=ANGLE(I)
            X(I)=R*COS(A)
            Y(I)=R*SIN(A)
         ENDDO
      ENDIF
C
      RETURN
      END
C******************************

	SUBROUTINE LLSQ1 (N,X,Y,D,C,CHI)

	DIMENSION X(1), Y(1), WORK(900), W2(900)

	SX = VECSUM (N,X)
	SY = VECSUM (N,Y)
	SXY = VECDOT (N,X,Y)
	SXX = VECDOT (N,X,X)

	D = ( SXY - SY*SX/N ) / ( SXX - SX*SX/N )
	C = ( SY - D*SX ) / N

	CALL VECSCAL (D,N,X,WORK)
	CALL VECSET (C,W2,1,N)
	CALL VECADD (N,WORK,W2,WORK)
	CALL VECSUB (N,WORK,Y,W2)
	CHI = SQRT ( VECDOT (N,W2,W2) ) / (N-2)

	RETURN
	END
