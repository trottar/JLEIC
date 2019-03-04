
*CMZ :  2.08/05 27/03/2000  16.08.25  by  Hannes Jung
*CMZ :  2.08/04 21/12/99  11.48.06  by  Hannes Jung
*CMZ :  2.08/00 27/05/99  14.56.13  by  Hannes Jung
*-- Author :
C----------------------------------------------------------------------
      FUNCTION FourDot(P,Q)
      DOUBLE PRECISION FourDot
C---LORENTZ 4-VECTOR DOT PRODUCT
      DOUBLE PRECISION P(5),Q(5)
      FourDot=P(4)*Q(4)-(P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3))
      END
