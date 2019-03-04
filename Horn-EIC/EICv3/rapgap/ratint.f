*CMZ :  2.08/04 14/07/99  18.40.56  by  Hannes Jung
*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 17/03/99  10*-- Author :
C===================================================================
C===================================================================

      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 10 I=1,N
         H=ABS(X-XA(I))
         IF (H.EQ.0.)THEN
            Y=YA(I)
            DY=0.0
            RETURN
         ELSE IF (H.LT.HH) THEN
            NS=I
            HH=H
         ENDIF
         C(I)=YA(I)
         D(I)=YA(I)+TINY
   10 CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 30 M=1,N-1
         DO 20 I=1,N-M
            W=C(I+1)-D(I)
            H=XA(I+M)-X
            T=(XA(I)-X)*D(I)/H
            DD=T-C(I+1)
            IF(DD.EQ.0.)PAUSE
            DD=W/DD
            D(I)=C(I+1)*DD
            C(I)=T*DD
   20    CONTINUE
         IF (2*NS.LT.N-M)THEN
            DY=C(NS+1)
         ELSE
            DY=D(NS)
            NS=NS-1
         ENDIF
         Y=Y+DY
   30 CONTINUE
      RETURN
      END
