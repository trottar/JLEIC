*CMZ :  2.07/03 15/05/99  14.33.38  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      FUNCTION DFUN(N,X)
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION DFUN,FXN1,WEIGHT,X(20)
      EXTERNAL FXN1
      DFUN=0
      WEIGHT=0.D0
      DO 10 I=1,N
10    IF(X(I).EQ.0..OR.X(I).EQ.1.) RETURN
      DFUN = FXN1(X,WEIGHT)
      RETURN
      END
