      PROGRAM LUMAIN

      USE PREC
      IMPLICIT NONE

      INTEGER, PARAMETER        :: N=2
      INTEGER                   :: I,J,K
      REAL(DP), DIMENSION (N,N) :: A0, A, L, U
      REAL(DP), DIMENSION (N)   :: B, E, X, PERM

c      CALL ONESMATRIX(A0,N)

      DO I=1,N
c        WRITE(*,*) A0(I,:)
      ENDDO

      a = reshape((/ 0.001,1.,1.,0. /), shape(a))
      DO I=1,N
        WRITE(*,*) A(I,:)
      ENDDO

      CALL LUFAC(N,A,L,U,PERM)
      
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
    
          print *,j
          DO I=1,N
            WRITE(*,*) A(I,:)
          ENDDO

      ENDDO
      L(N,N) = 1.
      print *,'a'
      DO I=1,N
        WRITE(*,*) A(I,:)
      ENDDO
      print *,'l'
      DO I=1,N
        WRITE(*,*) L(I,:)
      ENDDO
      print *,'u'
      DO I=1,N
        WRITE(*,*) U(I,:)
      ENDDO
      print *,'lu'
      LU=MATMUL(L,A)
      DO I=1,N
        WRITE(*,*) LU(I,:)
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
