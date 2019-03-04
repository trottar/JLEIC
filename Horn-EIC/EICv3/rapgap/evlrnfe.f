*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 17/03/99
*-- Author :
C===================================================================
C===================================================================

      SUBROUTINE EVLRNFE (FLNM, HEADER, IRET)
C         Read evolved parton density from named formatted file,
C         and set HEADER to the parton density name.
C         Parton density will be put in the standard common blocks
C         Error code returned.   IRET .EQ. 0 => no problem.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      character *(*) flnm, header, myname
      parameter (myname = 'EVLRNFE')

      iret = 0
      nunit = nextun()
      OPEN (Nunit, FILE=flnm, FORM='FORMATTED', STATUS='OLD',
     +      iostat = iret)
      IF (Iret .NE. 0) THEN
         PRINT *, 'I could not open the parton density file ', flnm
         PRINT *, 'Did you forget to make it?'
         PRINT *, 'I am ', myname
         return
      ENDIF

      CALL EVLRD (NUNIT, HEADER, IRET)
      IF (Iret .NE. 0) THEN
         PRINT *, 'Error reading from file in ', myname
      ENDIF

      END
