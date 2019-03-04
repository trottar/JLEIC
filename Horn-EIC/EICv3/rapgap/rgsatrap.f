*CMZ :  2.08/06 31/03/2000  11.20.05  by  Hannes Jung
*CMZ :  2.08/05 27/03/2000  15.57.41  by  Hannes Jung
*CMZ :  2.08/04 08/01/2000  16.54.02  by  Hannes Jung
*CMZ :  2.08/02 04/11/99  16.39.09  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.58  by  Hannes Jung
*-- Author :    Hannes Jung   29/05/99
      Subroutine RGSATRAP(X,WDIF)
      Implicit NONE
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.
      Double Precision X(20)
      Double Precision WDIF

      Double Precision  Xlow, Xup

      Double Precision  Stot, Q2s, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2s, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision    DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT
      Double Precision    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT
      COMMON /CCOMPCROS/  DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT,
     &                    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT

      Integer         NQUAD
      COMMON /CNQUAD/ NQUAD

      Double Precision  DSIGMax,  Mx, W, DUM
	Data DUM/0.d0/
      Double Precision  EPS
      Common  /mysatrap/EPS
      Double Precision  Kt2Min, Kt2Max
      Double Precision  ZWMin, ZWMax
      Double Precision  THREEDOT
      Double Precision  Pi
      Integer NMaxWe, Iev
      Double Precision W02,W12
      DOUBLE PRECISION ME,MP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get Parameters of the Wuesthoff Model

      Call SATURPAR
C -----------------------------------------------------

      WDIF = 0.
      ME =0.511e-3
      MP =0.938


      Stot = SSS

c calculate kinematically alowed y range:
      W02=(1.+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))


      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
c     write(6,*) 'rgsatrap ',ymin,ymax

c Mx limits moved to    RGQ2W2Mx

      NQUAD = 24

c     write(6,*) ' in RGSATRAP '

c Compute diffractive sigma's by MC computation

      Call RGQ2W2Mx(X)
      IF(AJAC.LE.0) RETURN
      Call INPUTKINE

      Call COMPCROSS

      EPS = 2.*(1.-YBj)/(1.+(1.-YBj)**2)

      DLTOT = DLQQT+DLQQGT + EPS*(DLQQL)
      DLCTOT= DCQQT+DCQQGT + EPS*(DCQQL+DCQQGL) + DLTOT

      IF(DLCTOT.LE.0.0) DLCTOT=0.
      WDIF = DLCTOT

      Return
      End
