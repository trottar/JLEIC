*CMZ :  2.06/32 03/05/98  10.20.39  by  Hannes Jung
*-- Author :    Hannes Jung   06/04/98
      SUBROUTINE DFRIDR(FUNC,X,H,DFDX,ERR)
c to calculate the derivative of function FUNC
c taken from: Numerical Recepies in Fortran
c             W. H. Press, S.A. Teulolsky, W.T. Vetterling, B.P. Flannery
c             Cambridge University Press 1992
      IMPLICIT NONE
      INTEGER  NTAB
      DOUBLE PRECISION DFDX,ERR,H,X,FUNC,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4D0,CON2=CON*CON,BIG=1.D30,NTAB=10,SAFE=2.D0)
      EXTERNAL FUNC
      INTEGER I,J
      DOUBLE PRECISION ERRT,FAC,HH,A(NTAB,NTAB)
      IF(H.EQ.0.D0) THEN
         write(6,*) 'DFRIFR: h must be nonzero: program stopped'
         STOP
      ENDIF
      HH=H
      A(1,1)=(FUNC(X+HH)-FUNC(X-HH))/(2.D0*HH)
      ERR=BIG
      DO I=2,NTAB
         HH=HH/CON
         A(1,I)=(FUNC(X+HH)-FUNC(X-HH))/(2.D0*HH)
         FAC = CON2
         DO J=2,I
            A(J,I)=(A(J-1,I)*FAC-A(J-1,I-1))/(FAC-1.D0)
            FAC = CON2*FAC
            ERRT=MAX(ABS(A(J,I)-A(J-1,I)),ABS(A(J,I)-A(J-1,I-1)))
            IF(ERRT.LE.ERR) THEN
               ERR=ERRT
               DFDX=A(J,I)
            ENDIF
         ENDDO
         IF(ABS(A(I,I)-A(I-1,I-1)).GE.SAFE*ERR) RETURN
      ENDDO
      RETURN
      END
