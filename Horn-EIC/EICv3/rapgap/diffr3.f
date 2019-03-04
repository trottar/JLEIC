*CMZ :  2.08/00 06/06/99  15.53.52  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  14.38.44  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE DIFFR3(X,WDIF)
      IMPLICIT NONE
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
      Double Precision GEV2NB,WPART,WDIF
	Double Precision X(20)
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDF(X,WPART)
      ELSE
         CALL PARTDFHS(X,WPART)
      ENDIF
      IF(WPART.LE.0.D0) WPART = 0.D0
c matrix element in partdf included
      WDIF=WPART*GEV2NB
      RETURN
      END
