*CMZ :  2.07/03 15/05/99  14.45.03  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94

      FUNCTION DOT1(I,J)
	IMPLICIT NONE
      INTEGER N,K
      REAL SP,V
      DOUBLE PRECISION P
*KEEP,RGLUPARM.
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


*KEND.
      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
      Double Precision DOT1
	INTEGER I,J
C+++++++++++
C    DOT PRODUCT OF FOUR VECTOR IN MINKOWSKI METRIK
C    WITH VECTORS FROM LUJETS
C+++++++++++
      DOT1= P(I,4)*P(J,4)-(P(I,1)*P(J,1))-(P(I,2)*P(J,2))-
     .     (P(I,3)*P(J,3))
      RETURN
      END
