*CMZ :  2.07/03 15/05/99  14.45.03  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      FUNCTION DOT(A,B)
      IMPLICIT NONE
C+++++++++++
C    DOT PRODUCT OF FOUR VECTOR IN MINKOWSKI METRIK
C+++++++++++
      Double Precision A(4),B(4),DOT
      DOT = A(4)*B(4)-A(1)*B(1)-A(2)*B(2)-A(3)*B(3)
      RETURN
      END
