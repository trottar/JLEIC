*CMZ :  2.08/04 14/07/99  18.40.56  by  Hannes Jung
*CMZ :  2.07/03 24/03/99  17.03.03  by  Hannes Jung
*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 17/03/99  10
*-- Author :
C================================================================
C================================================================

      FUNCTION NEXTUN()
C                                    Returns an unallocated FORTRAN i/o unit.
      LOGICAL EX
C
      DO 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         IF (.NOT. EX) then
            nextun = n
            RETURN
         end if
   10 CONTINUE
      RETURN
C               *************************
      END
