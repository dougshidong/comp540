      PROGRAM LUMAIN

      USE PREC
      IMPLICIT NONE

      INTEGER, PARAMETER        :: N=30
      INTEGER                   :: I,J,K
      REAL(DP), DIMENSION (N,N) :: A0, A, L, U
      REAL(DP), DIMENSION (N)   :: B, E, X, PERM

      CALL ONESMATRIX(A0,N)
      DO I=1,N
        A0(I,I) = A0(I,I) + 1.0E-8
      ENDDO

      E = 1.

      B = MATMUL(A0,E)

      DO I=1,N
        WRITE(*,*) A(I,:)
      ENDDO

      CALL LUFAC(N,A0,L,U,PERM)

      
      CALL SOLVELU(N,L,A0,B,X)
      print *, x
      
      END PROGRAM

      SUBROUTINE ONESMATRIX(A,N)

      USE PREC
      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: N
      REAL(DP), INTENT(OUT)   :: A(N,N)
      INTEGER                 :: I,J

      A=0.
      DO I=1,N
      DO J=1,N
        IF(I.EQ.J) THEN 
          A(I,J) = 1.
        ELSE IF(I.GT.J) THEN 
          A(I,J) = -1.
        ELSE IF(J.EQ.N) THEN 
          A(I,J) = 1.
        ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE
      
      SUBROUTINE LUFAC(N, A, L, U, PERM)

      USE PREC
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)        :: N
      REAL(DP), INTENT(INOUT)    :: A(N,N)
      REAL(DP), INTENT(OUT)      :: L(N,N), U(N,N), PERM(N)
      real(dp)                   :: LU(N,N)
      INTEGER                    :: I,J,K
      INTEGER                    :: IMAX
      REAL(DP)                   :: AMAX

      L = 0.
      U = 0.
      
      DO J=1,N-1

C       PARTIAL PIVOTING
        
        AMAX = ABS(A(J,J))
        IMAX = J
        DO I=J+1,N
          IF(AMAX.LT.ABS(A(I,J))) THEN
            AMAX = ABS(A(I,J))
            IMAX = I
          ENDIF
        ENDDO
        PERM(J) = IMAX

        IF(IMAX.NE.J) THEN
          CALL SWAPARRAY(A(J,J:N), A(IMAX,J:N), N)
          IF(J.GT.1) CALL SWAPARRAY(L(J,1:J-1), L(IMAX,1:J-1), N)
        ENDIF

C       FACTORIZATION

        DO I=J+1,N

          L(I,J) = A(I,J) / A(J,J)
          A(I,J) = 0
          DO K=J+1,N
            A(I,K) = A(I,K) - L(I,J) * A(J,K)
          ENDDO
        ENDDO

        L(J,J) = 1.
    
      ENDDO

      END SUBROUTINE

      SUBROUTINE SWAPARRAY(A,B,NN)

      USE PREC

      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: NN
      REAL(DP), INTENT(INOUT) :: A(NN),B(NN)
      INTEGER                 :: I
      REAL(DP)                :: TEMP

      DO I=1,NN
        TEMP = A(I)
        A(I) = B(I)
        B(I) = TEMP
      ENDDO

      END
      
      SUBROUTINE SOLVELU(N,L,U,B,X)
C     SOLVES LUX=B SYSTEM
C     INPUT:
C     N   -   SIZE OF MATRIX/VECTOR
C     L   -   LOWER TRIANGULAR MATRIX
C     U   -   UPPER TRIANGULAR MATRIX
C     B   -   N X 1 VECTOR
C     OUTPUT:
C     X   -   N X 1 SOLUTION
      USE PREC
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: L(N,N), U(N,N), B(N)
      REAL(DP), INTENT(OUT):: X(N)
      REAL(DP)             :: Y(N)

C     FORWARD SUBSTITUTION LY=B
      DO I=1,N
        Y(I) = B(I)
        DO J=1,I-1
          Y(I) = Y(I) - L(I,J) * Y(J)
        ENDDO
        Y(I) = Y(I) / L(I,I)
      ENDDO

C     BACK SUBSTITUTION UX=Y
      DO I=N,1,-1
        X(I) = Y(I)
        DO J=I+1,N
          X(I) = X(I) - U(I,J) * X(J)
        ENDDO
        X(I) = X(I) / U(I,I)
      ENDDO

      END SUBROUTINE
