      PROGRAM MAIN

      USE PRES

      IMPLICIT NONE


      INTEGER :: I 
      INTEGER, PARAMETER :: N = 6
      REAL(P), DIMENSION(N,N) :: A
      REAL(P), DIMENSION(N)   :: B

C     INITIALIZE A AND B MATRICES      
      CALL INIT(A,B,N)

      WRITE(*,*) 'A MATRIX'
      DO I = 1,N
         WRITE(*,*) A(I,:)
      ENDDO


      WRITE(*,*) 'B MATRIX'
      DO I = 1,N
         WRITE(*,*) B(I)
      ENDDO



      END PROGRAM


      SUBROUTINE INIT(A,B,N)
C     ****************************************
C     INITIALIZES A AND B MATRICES
C     A IS A HILBERT MATRIX N X N
C     B IS A N X 1 VECTOR DEFINED BY B=A.E
C     E IS THE (1,1,1,..)^T VECTOR
C     ****************************************
      USE PRES

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: N
      REAL(P), INTENT(OUT)   :: A(N,N), B(N)
      INTEGER                :: E(N), I, J

      E = 1
      DO I = 1,N
      DO J = 1,N
        A(I,J) = 1./(I+J-1)
      ENDDO
      ENDDO

      B = MATMUL(A,E)

      RETURN
      END SUBROUTINE

      SUBROUTINE CHOL


      IMPLICIT NONE


      END SUBROUTINE

C     FROBENIUS NORM OF A, DIMENSION M X N
C     RETURNS FNORM
      SUBROUTINE FROB(A,M,N,FNORM)

      USE PRES

      IMPLICIT NONE

      INTEGER             :: M,N,I,J
      REAL(P), INTENT(IN) :: A(M,N)
      REAL(P)             :: FNORM
      

      DO I = 1,M  
      DO J = 1,N
        FNORM = FNORM + A(I,J)**2
      ENDDO
      ENDDO

      FNORM = SQRT(FNORM)

      END SUBROUTINE
