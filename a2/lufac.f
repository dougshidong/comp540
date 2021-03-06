      PROGRAM LUMAIN

      USE PREC
      IMPLICIT NONE

      INTEGER, PARAMETER        :: N=30
      INTEGER                   :: I
      REAL(P2), DIMENSION (N,N) :: A, L, U, LU
      REAL(P2), DIMENSION (N)   :: B, E, X, R, D
      REAL(P2)                  :: GFAC, TOL, IMPROV
      INTEGER, DIMENSION (N)    :: IPERM
 
      CALL ONESMATRIX(A,N)
      DO I=1,N
        A(I,I) = A(I,I) + 1.0E-8
      ENDDO
      E = 1.

      B = MATMUL(A,E)

c      CALL PRINTARR(A,N)

      CALL LUFAC(N,A,L,U,IPERM)

C     GROWTH FACTOR
      GFAC = MAXVAL(U)/MAXVAL(A)
      WRITE(*,*) 'GROWTH FACTOR: ', GFAC
      WRITE(*,*) 

      LU = MATMUL(L,U)
      GOTO 10
      WRITE(*,*) 'A'
      CALL PRINTARR(A,N)
      WRITE(*,*) 'L'
      CALL PRINTARR(L,N)
      WRITE(*,*) 'U'
      CALL PRINTARR(U,N)
      WRITE(*,*) 'LU'
      CALL PRINTARR(LU,N)
      WRITE(*,*) 'PERMUTATION'
      WRITE(*,*) IPERM
   10 CONTINUE

      CALL PERMB(B,IPERM,N)
      CALL SOLVELU(N,L,U,B,X)
      
      WRITE(*,*) 'COMPUTED X'
      DO I=1,N
        WRITE(*,*) X(I)
      ENDDO
      WRITE(*,*) 

      WRITE(*,*) '||E-X|| / ||X||'
      WRITE(*,*) MAXVAL(ABS(E(1:N)-X(1:N)))/MAXVAL(X)
      WRITE(*,*) 

C     ITERATIVE REFINEMENT
      TOL = 1E-15
      DO I=1,99
        WRITE(*,'(A,I2)') 'ITERATIVE REFINEMENT #', I
        R = B - MATMUL(A, X)
        CALL SOLVELU(N,L,U,R,D)
        X = X + D
        IMPROV = NORM2(D) / NORM2(X)

        WRITE(*,*) '||E-X|| / ||X||'
        WRITE(*,*) MAXVAL(ABS(E(1:N)-X(1:N)))/MAXVAL(X)
        WRITE(*,*) '||D|| / ||X||'
        WRITE(*,*) IMPROV
        WRITE(*,*) 
        IF(IMPROV.LE.TOL) EXIT

      ENDDO
      
      WRITE(*,'(A,I2,A)') 'X AFTER ',I,' REFINEMENTS:'
      DO I=1,N
        WRITE(*,*) X(I)
      ENDDO

      END PROGRAM

C***********************************************************************
      SUBROUTINE ONESMATRIX(A,N)

      USE PREC
      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: N
      REAL(P2), INTENT(OUT)   :: A(N,N)
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
      
C***********************************************************************
      SUBROUTINE LUFAC(N, A, L, U, IPERM)
C     LU FACTORIZATION WITH PARTIAL PIVOTING
C     LU = PA
C     INPUT:
C     N     -    SIZE OF MATRIX/VECTOR
C     A     -    N X N MATRIX
C     OUTPUT:
C     L     -    LOWER TRIANGULAR MATRIX
C     U     -    UPPER TRIANGULAR MATRIX
C     IPERM -    PERMUTATION INDICES OF PERMUTATION MATRIX
      USE PREC
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)        :: N
      INTEGER, INTENT(OUT)       :: IPERM(N)
      REAL(P2), INTENT(IN)       :: A(N,N)
      REAL(P2), INTENT(OUT)      :: L(N,N), U(N,N)
      INTEGER                    :: I,J,K
      INTEGER                    :: IMAX
      REAL(P2)                   :: AMAX

      IPERM = 0
      L = 0.
      U = A

C     LOOP THROUGH COLUMNS   
      DO J=1,N-1

C       PARTIAL PIVOTING
        AMAX = ABS(U(J,J))
        IMAX = J
C       LOOP THROUGH ROWS TO FIND LARGEST PIVOT        
        DO I=J+1,N
          IF(AMAX.LT.ABS(U(I,J))) THEN
            AMAX = ABS(U(I,J))
            IMAX = I
          ENDIF
        ENDDO
        IPERM(J) = IMAX

C       SWAPPING ROWS
        IF(IMAX.NE.J) THEN
          CALL SWAPARRAY(U(J,J:N), U(IMAX,J:N), N-J+1)
          IF(J.GT.1) CALL SWAPARRAY(L(J,1:J-1), L(IMAX,1:J-1), J-1)
        ENDIF

C       FACTORIZATION
        DO I=J+1,N

          L(I,J) = U(I,J) / U(J,J)
          U(I,J) = 0.
          DO K=J+1,N
            U(I,K) = U(I,K) - L(I,J) * U(J,K)
          ENDDO
        ENDDO

      ENDDO

C     DIAGONAL ENTRIES OF L ARE 1.0
      DO I=1,N
        L(I,I) = 1.0
      ENDDO

      END SUBROUTINE

C***********************************************************************
      SUBROUTINE SWAPARRAY(A,B,NN)
C     SWAP ELEMENTS OF ARRAY A AND B OF SIZE NN
      USE PREC

      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: NN
      REAL(P2), INTENT(INOUT) :: A(NN),B(NN)
      INTEGER                 :: I
      REAL(P2)                :: TEMP

      DO I=1,NN
        TEMP = A(I)
        A(I) = B(I)
        B(I) = TEMP
      ENDDO

      END
      
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
      REAL(P2), INTENT(IN) :: L(N,N), U(N,N), B(N)
      REAL(P2), INTENT(OUT):: X(N)
      REAL(P2)             :: Y(N)

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
     
C***********************************************************************
      SUBROUTINE PERMB(B,IPERM,N)
C     PERMUTE THE ROWS OF B WITH IPERM      
      USE PREC 
      IMPLICIT REAL(P2)(A-H,O-Z)
      DIMENSION  B(N), IPERM(N)
      DO I=1,N-1
        IF(IPERM(I).NE.I) THEN
          TEMP = B(I)
          B(I) = B(IPERM(I))
          B(IPERM(I)) = TEMP
        ENDIF
      ENDDO
      END SUBROUTINE

C***********************************************************************
      SUBROUTINE PRINTARR(A,N)
C     PRINT ARRAY A WITH N ROWS 
      USE PREC 
      IMPLICIT REAL(P2)(A-H,O-Z)
      DIMENSION  A(N,N)
      DO I=1,N
        WRITE(*,*) A(I,:)
      ENDDO
      WRITE(*,*)
      ENDSUBROUTINE

