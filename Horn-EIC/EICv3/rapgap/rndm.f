*CMZ :  2.08/04 22/12/99  15.39.28  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  16.58.07  by  Hannes Jung
*-- Author :    Hannes Jung   19/08/97
      FUNCTION RNDM()
	Implicit None
	Double Precision RNDM
      Double Precision draprn
      EXTERNAL draprn
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         write(6,*) 'call rndm = draprn'
         FIRST = .FALSE.
      ENDIF
      RNDM = draprn()
      RETURN
      END
