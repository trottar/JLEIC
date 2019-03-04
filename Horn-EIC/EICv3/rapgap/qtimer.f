*CMZ :  2.06/39 07/07/98  15.33.52  by  Hannes Jung
*CMZ :  6.22/00 20/07/95  17.30.40  by  Martin Hampel
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     QTIMER(TIME)
C     TIMER ROUTINE, RETURNS THE TIME OF A CPU-CLOCK COUNTING IN
C     MICRO SECONDS (MACHINE DEPENDENT)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE QTIMER(RTIME)
      REAL*8 RTIME
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      IF(FIRST) THEN
         FIRST=.FALSE.
         NCL=0
      ENDIF
      NCL=NCL+1
      RTIME=NCL*1D6
      RETURN
      END
