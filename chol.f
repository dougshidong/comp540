      PROGRAM MAIN

      USE PREC

      IMPLICIT NONE


      INTEGER                 :: I, IERR
      INTEGER, PARAMETER      :: N = 8
      REAL(P), DIMENSION(N,N) :: A, L, LLTR
      REAL(P), DIMENSION(N)   :: B, E, X
      REAL(P)                 :: FROB

C     EXACT SOLUTION E
      E = 1.
C     INITIALIZE A AND B MATRICES
C     SUCH THAT AE=B
      CALL INIT(A,E,B,N)

      WRITE(*,*) 'A MATRIX'
      DO I = 1,N
        WRITE(*,*) A(I,:)
      ENDDO


      WRITE(*,*) 'B MATRIX'
      DO I = 1,N
        WRITE(*,*) B(I)
      ENDDO

      CALL CHOL(A,N,L,IERR)
      IF(IERR.LT.0) STOP

      WRITE(*,*) 'L MATRIX'
      DO I = 1,N
        WRITE(*,*) L(I,:)
      ENDDO
        
      LLTR = MATMUL(L,TRANSPOSE(L))

      WRITE(*,*) 'LL^T MATRIX'
      DO I = 1,N
        WRITE(*,*) LLTR(I,:)
      ENDDO

      X = 0.
      CALL SOLVELU(N,L,TRANSPOSE(L),B,X)
      WRITE(*,*) 'SOLUTION X: ', X

      WRITE(*,*) 'RELATIVE ERROR IN THE SOLUTION:'
      WRITE(*,*) NORM2(X-E)/NORM2(X)

      WRITE(*,*) 'RELATIVE RESIDUAL:'
      WRITE(*,*) (NORM2(B-MATMUL(A,X)) / (FROB(A,N,N)/NORM2(X)))

      WRITE(*,*) 'RELATIVE MATRIX RESIDUAL:'
      WRITE(*,*) FROB(A-LLTR,N,N)/FROB(A,N,N)

      END PROGRAM


C***********************************************************************
      SUBROUTINE CHOL(A,N,L,IERR)
C     CHOLESKY DECOMPOSITION OF SYMMETRIC POSITIVE DEFINITE MATRIX A
      USE PREC

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(OUT) :: IERR
      REAL(P), INTENT(IN)  :: A(N,N)
      REAL(P), INTENT(OUT) :: L(N,N)
      REAL(P), PARAMETER   :: TOL = 2e-7
      INTEGER              :: I,J

      L = 0.
      IERR = 0
      DO 10 J = 1,N
C       DIAGONAL VALUES     
        L(J,J) = A(J,J)
        IF(J.GE.2) L(J,J) = L(J,J) - DOT_PRODUCT(L(J,1:J-1),L(J,1:J-1))

C       CATCH ERROR IF SQUARE-ROOTING NEGATIVE NUMBER
        IF(L(J,J).LT.0.) THEN
          WRITE(*,*) 'NEGATIVE SQUARE ROOT'
          WRITE(*,*) 'I, L(I,I): ', J, L(J,J)
          IERR = -1
          EXIT
        ENDIF

        L(J,J)=SQRT(L(J,J))

        IF(ABS(L(J,J)).LT.TOL)THEN
          WRITE(*,*) 'L(J,J) IS CLOSE TO SINGLE PREC 0:'
          WRITE(*,*) 'I, L(I,I): ', J, L(J,J)
          IERR = -2
          EXIT
        ENDIF


C       OFF-DIAGONAL VALUES
        DO 20 I = J+1,N
          L(I,J) =( A(I,J) - 
     .              DOT_PRODUCT(L(I,1:J-1),L(J,1:J-1)) ) / L(J,J)     
   20   CONTINUE

   10 CONTINUE



      END SUBROUTINE
C***********************************************************************
      SUBROUTINE INIT(A,E,B,N)
C     ****************************************
C     INITIALIZES A AND B MATRICES
C     A IS A HILBERT MATRIX N X N
C     B IS A N X 1 VECTOR DEFINED BY B=A.E
C     E IS THE (1,1,1,..)^T VECTOR
C     ****************************************
      USE PREC

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: N
      REAL(P), INTENT(IN)    :: E(N)
      REAL(P), INTENT(OUT)   :: A(N,N), B(N)
      INTEGER                :: I, J

      DO J = 1,N
      DO I = 1,N
        A(I,J) = 1./(I+J-1)
      ENDDO
      ENDDO
c      A = MATMUL(C,TRANSPOSE(C))

      B = MATMUL(A,E)

      RETURN
      END SUBROUTINE

C***********************************************************************
      REAL(P) FUNCTION FROB(A,M,N)
C     FROBENIUS NORM OF A, DIMENSION M X N
C     RETURNS FNORM
      USE PREC

      IMPLICIT NONE

      INTEGER             :: M,N,I,J
      REAL(P), INTENT(IN) :: A(M,N)
      REAL(P)             :: FNORM
      
      FNORM = 0.
      DO J = 1,M  
      DO I = 1,N
        FNORM = FNORM + A(I,J)*A(I,J)
      ENDDO
      ENDDO

      FROB = SQRT(FNORM)

      RETURN

      END FUNCTION

C***********************************************************************
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
      REAL(P), INTENT(IN) :: L(N,N), U(N,N), B(N)
      REAL(P), INTENT(OUT):: X(N)
      REAL(P)             :: Y(N)

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
c      SUBROUTINE RELERR(X,E,)
     

c      END SUBROUTINE
