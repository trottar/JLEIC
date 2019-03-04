*CMZ :  2.08/00 06/06/99  15.47.21  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  14.41.38  by  Hannes Jung
*-- Author :    Hannes Jung   12/01/95
      SUBROUTINE DIR3(X,WDIF)
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
         CALL PARTDI(X,WPART)
      ELSE
         CALL PARTDIHS(X,WPART)
      ENDIF
c matrix element in partdi included
      IF(WPART.LE.0.D0) WPART=0.D0
      WDIF=WPART*GEV2NB
c      write(6,*) ' dir3 ',wpart,wdif
      RETURN
      END
