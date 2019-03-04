*CMZ :  2.08/00 27/05/99  14.56.13  by  Hannes Jung
*-- Author :
C----------------------------------------------------------------------

      Function PHIFUNC(Lt2)
      Implicit NONE

      Double Precision  PHIFUNC,Lt2
      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS
      PHIFUNC = 3.*SIGM0*R02*DEXP(-R02*Lt2)
      PHIFUNC = PHIFUNC/4./9.8696
      END
