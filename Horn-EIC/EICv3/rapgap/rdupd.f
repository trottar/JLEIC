*CMZ :  2.07/01 17/03/99  10*-- Author :
C================================================================
C================================================================

      SUBROUTINE RDUPD (UPD, NTL, NDAT, IRET)
C
C                       The I/O operation is made a stand-alone subprogram
C                       so that fast execution is achieved with block
C                       data transfer for the actual size of the array
C                       (rather than the declared size in the main program).
      DOUBLE PRECISION UPD
      DIMENSION UPD (NTL)

      READ (NDAT, *, IOSTAT=IRET) UPD
C                             To check read-in result against saved file.
C      Print '(1pE13.5, 5E13.5)', UPD
C

      END
