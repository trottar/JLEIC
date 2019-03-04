*CMZ :  2.07/03 24/05/99  17.28.48  by  Hannes Jung
*-- Author :    Hannes Jung   02/03/97
      SUBROUTINE DIFFR7(X,WDIF)
      IMPLICIT None
	Double Precision X,WDIF,WPART,WT1,GEV2NB
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
      DIMENSION X(20)
	Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDF(X,WPART)
      ELSE
         write(6,*) ' resloved photon with HERACLES not implemented '
      ENDIF
      IF(WPART.GT.0.) THEN
         CALL ELERES(WT1)
      ENDIF
      WDIF=WPART*WT1*GEV2NB
      IF(WDIF.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' DIFFR7 WPART,WT1 ',WPART,WT1
      ENDIF
c      IF(WDIF.LE.0.D0) THEN
c         write(6,*) ' DIFFR7 WPART,WT1 ',WPART,WT1
c      ENDIF
      RETURN
      END
