*CMZ :  2.08/05 27/03/2000  16.14.35  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      FUNCTION SWAVET(Lt2)
C----------------------------------------------------------------------
      Implicit NONE

      Double Precision  Lt2, PHIFUNC, SWAVET

      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2

      Double Precision  WAVEWUT, QW2
      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Double Precision    MQU2
      MQU2 = MQUARK**2
c...
c  QQBar Transverse Wave Function Convoluted with the
c  proton vituality function
c

      PHIFUNC = 3.*SIGM0*R02*DEXP(-R02*Lt2)
      PHIFUNC = PHIFUNC/4./9.8696

      QW2 = (Kt2+MQU2)/(1.-BETA)
      WAVEWUT = 1.-2.*BETA -2.*MQU2/QW2+(Lt2 -(1.-2.*BETA)*QW2+2.*MQU2)
     &        /sqrt((Lt2+QW2)**2 - 4.*Kt2*Lt2)
      SWAVET = WAVEWUT*PHIFUNC

c      write(*,1000)  SWAVET,PHIFUNC, Lt2,BETA,Kt2
 1000 Format('SWAVET ',7E12.3)
 999  RETURN
      END
