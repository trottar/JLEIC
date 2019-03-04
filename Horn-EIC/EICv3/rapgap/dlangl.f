*CMZ :  2.07/03 15/05/99  14.46.46  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      FUNCTION DLANGL(X,Y)
      IMPLICIT NONE
	Double Precision PI,DLANGL,R,X,Y
      PI = 4.D0*DATAN(1.D0)
      DLANGL=0.D0
      R=DSQRT(X**2+Y**2)
      IF(R.LT.1D-20) RETURN
      IF(DABS(X)/R.LT.0.8D0) THEN
         DLANGL=DSIGN(DACOS(X/R),Y)
      ELSE
         DLANGL=DASIN(Y/R)
         IF(X.LT.0.d0.AND.DLANGL.GE.0.d0) THEN
            DLANGL=PI-DLANGL
         ELSEIF(X.LT.0.) THEN
            DLANGL=-PI-DLANGL
         ENDIF
      ENDIF

      RETURN
      END
