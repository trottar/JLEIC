*CMZ :  2.08/00 06/06/99  17.48.16  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  14.37.42  by  Hannes Jung
*-- Author :    Hannes Jung   30/10/94
      SUBROUTINE DIFFR30(X,WDIF)
      IMPLICIT NONE
*KEEP,RGPARAS.
      DOUBLE PRECISION PT2CUT,THEMA,THEMI,Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPRO
      COMMON/RAPA /IPRO,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA
      COMMON/PTCUT/ PT2CUT(100)
      COMMON/ELECT/ THEMA,THEMI
      REAL ULALPS,ULALEM
      EXTERNAL ULALPS,ULALEM
C     SAVE

*KEND.
      Double Precision WPART,WDIF
	Double Precision X(20)
      WPART=0.D0
      WDIF=0.D0
c      write(6,*) 'diffr30 WPART before rgsatrap'
      CALL RGSATRAP(X,WPART)
c      write(6,*) 'diffr30 WPART after rgsatrap'
c      WDIF=WPART*GEV2NB
      WDIF=WPART
c      write(6,*) 'diffr30: WDIF ',WDIF
c      call lulist(1)
      RETURN
      END
