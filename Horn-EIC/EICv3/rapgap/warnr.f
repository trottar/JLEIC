*CMZ :  2.07/02 24/03/99  16.21.10  by  Hannes Jung
*-- Author :
C===================================================================
C===================================================================

      SUBROUTINE WARNR (IWRN, NWRT, MSG, NMVAR, VARIAB,
     +                  VMIN, VMAX, IACT)

C   Subroutine to handle warning messages.  Writes the (warning) message
C   and prints out the name and value of an offending variable to SYS$OUT
C   the first time, and to output file unit # NWRT in subsequent times.
C
C   The switch IACT decides whether the limits (VMIN, VMAX) are active or
C   not.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      CHARACTER*(*) MSG, NMVAR
      Data Nmax / 10 /

      IW = IWRN
      VR = VARIAB

      IF  (IW .EQ. 0) THEN
         PRINT '(1X, A/1X,2A,1PD16.7/A,I4)', MSG, NMVAR, ' = ', VR,
     +         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
            PRINT '(A/2(1PE15.4))', ' The limits are: ', VMIN, VMAX
            WRITE (NWRT,'(A/2(1PE15.4))') ' The limits are: ', VMIN,
     +      VMAX
         ENDIF
      ENDIF

      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(I5, 2A/1X,2A,1PD16.7)') IW, '   ', MSG,
     +                  NMVAR, ' = ', VR
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     +     '!! Too many warnings, Message suppressed from now on !!'
      Endif

      IWRN = IW + 1

      RETURN
C           ****************************
      END
