*CMZ :  2.07/03 15/05/99  14.21.38  by  Hannes Jung
*CMZ :  2.01/09 27/02/96  09.50.59  by  Hannes Jung
*CMZ :  1.04/01 23/02/95  11.50.45  by  Hannes Jung
*CMZ :  1.03/01 03/04/94  16.45.55  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      FUNCTION AALAM(A,B,C)
      IMPLICIT NONE
	double precision A,B,C,AALAM
      AALAM = A**2 + B**2 + C**2 - 2.D0*A*B - 2.D0*B*C - 2.D0*C*A
      RETURN
      END
