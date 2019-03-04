*CMZ :  2.08/00 14/06/99  19.13.57  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  14.42.27  by  Hannes Jung
*-- Author :    Hannes Jung   02/03/97
      SUBROUTINE DIR99(X,WDIF)
      IMPLICIT NONE
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
      Double Precision GEV2NB,WPART,WDIF,WT1
      Double Precision X(20)
      INTEGER IDEBUG
      COMMON/ INTERN/IDEBUG
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0

      CALL PARTDI99(X,WPART)

      IF(WPART.GT.0.) THEN
c         CALL ELERES(WT1)
         WT1 = 1.D0
      ENDIF
cc      write(6,*) wpart,wt1
      WDIF=WPART*WT1*GEV2NB
      IF(WDIF.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' DIR99 WPART,WT1 ',WPART,WT1
      ENDIF
      RETURN
      END
