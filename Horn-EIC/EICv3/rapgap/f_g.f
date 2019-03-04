
      DOUBLE PRECISION FUNCTION F_G(Y)
*     IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      implicit none
      DOUBLE PRECISION XPOM,Y,FT,DGRV_NL
      REAL DGRV_NLN
      DOUBLE PRECISION YOLD,XPOLD,FTOLD
      COMMON/BARTELS/XPOM
      EXTERNAL DGRV_NLN,DGRV_NL
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         YOLD = 0
         XPOLD = 0
         FTOLD = 0
         FIRST=.FALSE.
      ENDIF
      IF(Y.EQ.YOLD.AND.XPOM.EQ.XPOLD) THEN
         FT = FTOLD
      ELSE
         YOLD = Y
         XPOLD = XPOM
         IF(XPOM.GT.0.095D0) THEN
            FT = 0.0D0
         ELSE
            FT = DBLE(DGRV_NLN(XPOM,Y))
            IF(FT.LE.0) THEN
c            write(6,*) ' F_G,xpom,Y :',FT,XPOM,Y
               FT = 0.D0
            ENDIF
         ENDIF
         FTOLD = FT
      ENDIF
      F_G = FT
      RETURN
      END
