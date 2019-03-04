*CMZ :  2.08/04 22/12/99  15.39.36  by  Hannes Jung
*CMZ :  2.06/39 07/07/98  15.33.52  by  Hannes Jung
*CMZ :  2.06/02 19/08/97  16.26.26  by  Hannes Jung
*-- Author :    Hannes Jung   19/08/97
      DOUBLE PRECISION FUNCTION HSRNDM()
      Double Precision draprn
      EXTERNAL draprn
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         write(6,*) ' change of random number generator:               '
     +   //'              call hsrndm = draprn'
         FIRST = .FALSE.
      ENDIF
      HSRNDM = draprn()
      RETURN
      END
