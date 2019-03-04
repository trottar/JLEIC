*CMZ :  2.08/04 14/07/99  18.40.56  by  Hannes Jung
*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 17/03/99  10*-- Author :
C===================================================================
C===================================================================

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 10 I=1,N
         DIFT=ABS(X-XA(I))
         IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
         ENDIF
         C(I)=YA(I)
         D(I)=YA(I)
   10 CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 30 M=1,N-1
         DO 20 I=1,N-M
            HO=XA(I)-X
            HP=XA(I+M)-X
            W=C(I+1)-D(I)
            DEN=HO-HP
            IF(DEN.EQ.0.)PAUSE
            DEN=W/DEN
            D(I)=HP*DEN
            C(I)=HO*DEN
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
