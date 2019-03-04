*CMZ :  2.08/05 27/03/2000  16.12.12  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      Subroutine  COMPCROSS

      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF
      Double Precision      CHARFAC
      COMMON /CCHARFAC/     CHARFAC

      Double Precision  DSQQT, DSQQL, DSQQGT, DSQQGL
      COMMON /CSWUEST/  DSQQT, DSQQL, DSQQGT, DSQQGL

      Double Precision    DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT
      Double Precision    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT
      COMMON /CCOMPCROS/  DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT,
     &                    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT

      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2
      Double Precision    KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT

      Double Precision    Mx

      Mx = dsqrt(Mx2)
c Compute cross sections for light quarks

      MQUARK = 0.0001
      CHARFAC = 6./9.
      if(Mx.lt.1.) CHARFAC = 5./9.

      DSQQT  = 0.
      DSQQL  = 0.
      DSQQGT = 0.
      DSQQGL = 0.

       Call RAN2BODY
       Call RAN3BODY

       Call INTWAVE
       Call SIGDIF

      DLQQT  = DSQQT
      DLQQL  = DSQQL
      DLQQGT = DSQQGT
      DLQQGL = DSQQGL
      DLQQGL = 0.

cc Compute cross section with Charm

      DSQQT  = 0.
      DSQQL  = 0.
      DSQQGT = 0.
      DSQQGL = 0.

      if(Mx.gt.(2.*MCHARM+0.05)) Then
        MQUARK = MCHARM
        CHARFAC = 4./9.

        Call RAN2BODY
        Call RAN3BODY

        Call  INTWAVE
        Call  SIGDIF
      Endif

      DCQQT  = DSQQT
      DCQQL  = DSQQL
      DCQQGT = DSQQGT
      DCQQGL = DSQQGL

      Return
      END
