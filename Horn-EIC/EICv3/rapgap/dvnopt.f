*CMZ :  2.08/04 22/12/99  15.39.24  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.08.11  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  21.04.02  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94

C***********************************************************************
C         Change DIVON default options                                 C
C                                                                      C
C***********************************************************************
      SUBROUTINE DVNOPT
C

      COMMON /D151DT/ IDATE
      DOUBLE PRECISION IDATE
      COMMON /PRINT/ IPRINT
      COMMON /ISTRGE/ MXRGNS , ISTOR(12000)
      COMMON /RSTRGE/ RSTSZE,RSTOR(18001)
      INTEGER RSTSZE
      COMMON /QUADRE/ IDEG
      COMMON /START/ ISTART
      COMMON /EXFILE/ NFILE
      COMMON /DISPOS/ IDISP
      COMMON /DEPTHS/ ISTDPH , INCDPH
      COMMON /SAMPLE/ NPOINT
      COMMON /CUTOLS/ BNDTOL, FRACT, REGNTL, FNLTOL
      COMMON /BNDLMT/ FLOBD,FUPBD
      COMMON /PRSTOP/ NSTOP
      COMMON /ZEETRM/ ITRMF
C
      IPRINT= 10
      IDEG=1
      FLOBD=0.0
C
      PRINT 10000
10000 FORMAT(44H0DIVON4 default options are altered by user.  )
C     WRITE(6,*) ' NPOINT = ',NPOINT,' IDEG = ',IDEG
      RETURN
      END
