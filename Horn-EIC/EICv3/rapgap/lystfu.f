*CMZ :  2.08/02 14/07/99  18.40.56  by  Hannes Jung
*CMZ :  2.08/00 15/06/99  17.22.43  by  Hannes Jung
*CMZ :  2.06/39 07/07/98  15.33.52  by  Hannes Jung
*-- Author :
C######################################################################
C
C   Various routines to give structure function parametrizations.
C
C ********************************************************************

c...hs taken from LEPTO 6.5 (modified)
      SUBROUTINE LYSTFU(KF,X,Q2,XPQ)
      COMMON /ARSTRF/ KFSAVE(2),XSAVE(2),XQ2SAV(2),XPQSAV(2,-6:6)
      REAL*4          PYSTOP,PYSLAM
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMOD,NPYMAX,NPYMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12),IPART
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      DOUBLE PRECISION POLARI
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      DOUBLE PRECISION SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      DOUBLE PRECISION MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION XPQ(-6:6),XPQS(-6:6)
      DOUBLE PRECISION HSLOQS,DQ2,DX
*KEEP,RGLQ2.
      DOUBLE PRECISION Q2SUPP
      COMMON/LOWQ2S/Q2SUPP
*KEND.
C...hs
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./

      IF (LFIRST) THEN
         LFIRST=.FALSE.
         write(6,*) ' RAPGAP version for LYSTFU '
      ENDIF

C...Reset arrays etc.
      DO 10  KFL=-6,6
         XPQ(KFL)=0.0
   10 XPQSAV(1,KFL)=0.
      XSAVE(1)=X
      XQ2SAV(1)=Q2
      KFSAVE(1)=KF
C...Check x and particle species.
      IF(X.LE.0..OR.X.GE.1.) THEN
         WRITE(6,10000) X
         RETURN
      ENDIF

      CALL PYSTFU(KF,X,Q2,XPQS)
      DO 20  KFL=-NPYMAX,NPYMAX
   20 XPQ(KFL)=XPQS(KFL)

C...H.S.
C...APPLY LOW Q2 SUPPRESSION
      IF (Q2SUPP.NE.0) THEN
         DCORRP=1.-EXP(-Q2SUPP*Q2)
         DO 30 IF=1,6
            XPQ(IF)=XPQ(IF)*DCORRP
   30    XPQ(-IF)=XPQ(-IF)*DCORRP
      ENDIF

      DO 40  KFL=-6,6
   40 XPQSAV(1,KFL)=XPQ(KFL)

C...Formats for error printouts.
10000 FORMAT(' Error in LYSTFU: x =',1P,E12.4,' outside physical range')

      RETURN
      END
