C...  END INITIALIZE
*CMZ :  2.08/04 22/12/99  15.39.24  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  20.59.48  by  Hannes Jung
*CMZ :  2.03/08 03/10/96  20.08.56  by  Julian Phillips
*-- Author :    Julian Phillips   01/10/96
*#**********************************************************************
*#
*#    SUBROUTINE GRAWQ2
*#
*# PURPOSE: Decide whether to write event for weighted files
*#
*# INPUT: RAPGAP Common RGPARAS
*#
*# OUTPUT: IWRITE = 0 If event to be rejected
*#                = 1 If event to be kept
*#
*# CALLED BY: GRAPGA
*#
*# CALLING: draprn()
*#
*# AUTHOR: Julian Phillips                 CREATED AT: 96/10/01
*#**********************************************************************
      SUBROUTINE GRAWQ2(IWRITE,Q2EG)
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

*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEND.
      Integer NCHK
      PARAMETER (NCHK=1000)
      INTEGER  IWRITE
      Double Precision draprn
	Real Q2EG
      EXTERNAL draprn
      DOUBLE PRECISION WTEST

      INTEGER ILOOP,ICALL,IKEEP
      DATA ILOOP/0/, ICALL/0/, IKEEP/0/

      IWRITE = 0
      ICALL = ICALL + 1
      ILOOP = ILOOP + 1

* Contruct weight
      W_Q2 = 1.D0
      IF(Q2.GT.0.D0) W_Q2 = DMAX1(1.D0,Q2START/DBLE(Q2EG))

* Keep fraction 1/W_Q2 of sample, weighted by W_Q2
      WTEST = draprn()
      IF(WTEST.LT.1.D0/W_Q2) THEN
         IWRITE=1
         IKEEP = IKEEP + 1
      ENDIF

      IF(ILOOP.EQ.NCHK) THEN
         ILOOP = 0
         WRITE(6,*)'Weighting: from ',ICALL,' kept ',IKEEP
      ENDIF

      RETURN
      END
