*CMZ :  2.08/05 27/03/2000  16.15.24  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      FUNCTION SWAVEG(Lt2)
C----------------------------------------------------------------------
      Implicit NONE

      Double Precision  Lt2, PHIFUNC, SWAVEG

      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision    KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT

      Double Precision  WAVEWUG
      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Double Precision    MQU2
      MQU2 = MQUARK**2

c...
c  QQG Leading Twist Wave Function Convoluted with the
c  proton vituality function
c
c  factor 2 used for R02 because of gluon charge
      PHIFUNC = 3.*SIGM0*R02*DEXP(-R02*Lt2)
      PHIFUNC = PHIFUNC/4./9.8696

      WAVEWUG =  ZW**2+(1.-ZW)**2 + Lt2/KW2 -
     &   (((1.-2.*ZW)*KW2-Lt2)**2 + 2.*ZW*(1.-ZW)*KW2**2)/KW2/
     &   dsqrt((KW2+Lt2)**2-4.*(1.-ZW)*Lt2*KW2)

      SWAVEG = WAVEWUG*PHIFUNC

c      write(*,1000)  SWAVEL,PHIFUNC, Lt2
 1000 Format('SWAVET ',7E12.3)
 999  RETURN
      END
