      SUBROUTINE RAPGEN(NDIMEN,XG)
      IMPLICIT NONE
      REAL ACC1,ACC2
	Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
	Double Precision XG
	Double Precision XGS,FXNB
	Integer NDIMEN,MXTRY,I
      DIMENSION XG(20)
      COMMON/XFXNB/XGS(20)
      EXTERNAL FXNB
      IF(IINT.EQ.1) THEN
         CALL RANGEN(NDIMEN,XG)
      ELSE
         MXTRY = 500
         CALL SPRING( FXNB, MXTRY )
         DO I=1,NDIMEN
            XG(I) =XGS(I)
         ENDDO
      ENDIF
      RETURN
      END
