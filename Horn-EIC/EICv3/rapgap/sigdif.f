
*CMZ :  2.08/05 27/03/2000  16.14.13  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      Subroutine  SIGDIF
c
c author H.Kowalski
c COMPCROSS computes the diffractive X-section dsig/Mx in nB for the
c QQb  transverse   process - stored in DSQQT
c QQb  longitudinal process - stored in DSQQL
c QQbG transverse   process - stored in DSQQGT
c QQbG longitudinal process - stored in DSQQGL
c
      Implicit NONE
      Double Precision  FSIGT, FSIGL, FSIGGT, FSIGGL, DUM
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF
      Double Precision      CHARFAC
      COMMON /CCHARFAC/     CHARFAC

      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2
      Double Precision    KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT

      Double Precision   QW2, Mx, WTPhS, MQU2, PI
      Double Precision   FacT, FACL, FACGT, FACGL, FACQM

      Double Precision  IntWUT, INTWUL, INTWUG
      COMMON /CINTWAVE/ IntWUT, INTWUL, INTWUG

      Double Precision  DSQQT, DSQQL, DSQQGT, DSQQGL, DSTOT
      COMMON /CSWUEST/  DSQQT, DSQQL, DSQQGT, DSQQGL, DSTOT

      MQU2 = MQUARK**2
      Mx = dsqrt(Mx2)
      PI = 3.1415927

      QW2 = (Kt2+MQU2)/(1.-BETA)
      WTPhS = (Kt2+MQU2)/Kt2/dsqrt(1.-4.*BETA*QW2/Q2)/Mx2
      WTPhS = WTPhS*BETA/(1.-BETA)/Q2
      FacT = (1.-2.*BETA*QW2/Q2)*(IntWUT**2) +
     &        4.*Kt2*MQU2/QW2**2*(IntWUL**2)
      FacT = FacT * (1./137.)*(PI**2/12.)*1000000./2.568

      FSIGT = WTPhS*FacT*2.*Mx*CHARFAC*(1./6.)
      DSQQT = FSIGT*WT2*AJAC

c      write (*,*) 'SIGDIF',FacT,WTPhs,Wt2,AJAC

      WTPhS = 1./dsqrt(1.-4.*BETA*QW2/Q2)/Mx2
      FacL = (QW2*BETA**3/Q2/Q2)*(IntWUL**2)
      FacL = FacL *(1./137.)*(4.*PI**2/3.)*1000000./2.568
      FSIGL = WTPHS*FacL*2.*Mx*CHARFAC*(1./6.)
      DSQQL = FSIGL*WT2*AJAC

      DSQQGT  = 0.
      DSQQGL  = 0.
      IF(Mx.lt.1.) Return
cH.K 19.12.99 introduce charm mass
      FACQM = (QtG2+MQU2)
c instead of
c old      FACG = (1./QtG2)*((BETA/ZW)**2+(1.-BETA/ZW)**2)
c  new

      FACGT = ((1.-2.*FACQM/SM2)/FACQM)*((BETA/ZW)**2+(1.-BETA/ZW)**2)
c      write (*,*) 'SIGDIF',FacGT,FACQM,SM2,BETA,ZW
      FACGT = FACGT+4./FACQM*(BETA/ZW)*(1.-2.*BETA/ZW)*MQU2/Q2
      FACGT = FACGT+2.*SM2/FACQM**2*(BETA/ZW)**2*MQU2/Q2*
     1                 (1.-2.*MQU2/Q2)
cH.K 19.12.99
      FACGT = FACGT*((BETA/ZW)**2)/(1-ZW)**3
      FACGT = FACGT*(INTWUG**2)

      FACGT  =  FACGT*(9.*PI/32.)*(1./137.)*ALPHAS/Q2/Q2

      FSIGGT = (FACGT *1000000./2.568)*2.*Mx*CHARFAC*(1./6.)

c      write (*,*) 'SIGDIF',FacGT,INTWUG

c  new H.K 19.12.99
      FACGL = 8./SM2*(BETA/ZW)*(1.-BETA/ZW)
      FACGL = FACGL-8./FACQM*(BETA/ZW)**2*MQU2/Q2

      FACGL = FACGL*((BETA/ZW)**2)/(1-ZW)**3
      FACGL = FACGL*(INTWUG**2)

      FACGL = FACGL*(9.*PI/32.)*(1./137.)*ALPHAS/Q2/Q2

      FSIGGL = (FACGL *1000000./2.568)*2.*Mx*CHARFAC*(1./6.)

      if(WT3.gt.0.0) THEN
      DSQQGT  = FSIGGT* WT3*AJAC
      DSQQGL  = FSIGGL* WT3*AJAC
      Endif
      Return
      END
