*CMZ :  2.08/04 22/12/99  15.39.28  by  Hannes Jung
*CMZ :  2.06/32 14/04/98  17.02.34  by  Hannes Jung
*-- Author :    Hannes Jung   14/04/98
      FUNCTION RLU()
      Double Precision draprn
      EXTERNAL draprn
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         write(6,*) 'call rlu = draprn'
         FIRST = .FALSE.
      ENDIF
      RLU =  draprn()
      RETURN
      END
