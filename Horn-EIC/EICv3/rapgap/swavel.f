*CMZ :  2.08/05 27/03/2000  16.14.59  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      FUNCTION SWAVEL(Lt2)
C----------------------------------------------------------------------
      Implicit NONE

      Double Precision  Lt2, PHIFUNC, SWAVEL

      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2

      Double Precision  WAVEWUL, QW2
      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Double Precision    MQU2
      MQU2 = MQUARK**2

c...
c  QQBar Longitudinal Wave Function Convoluted with the
c  proton vituality function
c

      PHIFUNC = 3.*SIGM0*R02*DEXP(-R02*Lt2)
      PHIFUNC = PHIFUNC/4./9.8696

c HK 19.12.99
c old      QW2 = Kt2/(1.-BETA)
      QW2 = (Kt2+MQU2)/(1.-BETA)

      WAVEWUL = 1.-
     &      QW2/sqrt((Lt2+QW2)**2 -4.*Kt2*Lt2)

      SWAVEL = WAVEWUL*PHIFUNC

c      write(*,1000)  SWAVEL,PHIFUNC, Lt2
 1000 Format('SWAVET ',7E12.3)
 999  RETURN
      END
