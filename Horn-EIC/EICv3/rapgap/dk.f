      DOUBLE PRECISION FUNCTION DK(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      DOUBLE PRECISION a1ma1,Q2_c
      COMMON /QQG_CON/a1ma1,Q2_c
      DK = a1ma1*q2_c + X
      IF(DK.LT.0.) THEN
c        write(6,*) ' DK,a1ma1,q2_c,X ',DK,a1ma1,q2_c,X
      ENDIF
      RETURN
      END
