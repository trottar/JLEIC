*CMZ :  6.12/02 22/03/95  12.29.13  by  Martin Hampel
*CMZ :  6.00/00 10/05/94  13.50.16  by  Martin Hampel
*-- Author :
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ck DJlepto6vh:
Ck 14/04/93 routines from LEPTO6.1 used by DJANGO6 but unmodified:

      SUBROUTINE LTIMEX(TIME)
C...Interface routine to transfer a call to some machine-dependent
C...routine to get the execution time used since job started.
C...Nice, but not necessary information. Can also be called by user.

      TIME=0.
C...Use of CERN library routine Z007, replace/delete if not available.
ck      CALL TIMEX(TIME)
C-MH activated...
      CALL TIMEX(TIME)
      RETURN
      END
