*CMZ :  2.08/04 22/12/99  15.39.29  by  Hannes Jung
*CMZ :  2.06/34 20/05/98  09.12.56  by  Hannes Jung
*CMZ :  2.06/31 25/03/98  16.02.23  by  Hannes Jung
*CMZ :  2.06/28 12/03/98  14.44.31  by  Hannes Jung
*CMZ :  2.05/05 20/03/97  12.01.08  by  Hannes Jung
*CMZ :  2.04/00 23/12/96  11.32.31  by  Hannes Jung
*CMZ :  2.03/07 16/09/96  08.58.52  by  Hannes Jung
*CMZ :  2.03/05 06/09/96  08.15.27  by  Hannes Jung
*CMZ :  2.03/04 27/08/96  17.16.20  by  Hannes Jung
*CMZ :  2.03/03 21/08/96  09.28.27  by  Hannes Jung
*CMZ :  2.03/02 18/08/96  17.43.17  by  Hannes Jung
*CMZ :  2.02/03 01/08/96  09.40.25  by  Hannes Jung
*CMZ :  2.02/02 25/07/96  12.26.13  by  Hannes Jung
*CMZ :  2.02/00 02/07/96  17.40.39  by  Hannes Jung
*CMZ :  2.01/19 24/06/96  12.34.52  by  Hannes Jung
*CMZ :  2.01/18 20/06/96  11.21.16  by  Hannes Jung
*CMZ :  2.01/17 05/06/96  15.53.05  by  Hannes Jung
*CMZ :  2.01/15 19/05/96  13.20.55  by  Hannes Jung
*-- Author :    Hans Lotter   19/05/96
      SUBROUTINE F2BLW(BET,X_POM,Q2,T2,F2QT,FLQ)
* calculate  F_2^D in the hard approach of J. Bartels, H. Lotter, M. Wuesthoff
*                                          (DESY 96-026)
* the subroutines come from H. Lotter
* THANKS A LOT
*

      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      common   /phi/       phi
      common   /parameter/ s,q,beta,xp,t
      INTEGER IPHI
      COMMON /SEL/IPHI
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGPARA.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
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

      REAL * 4 SLO,SHI,FBLW,XEPS,RESULT,ULMASS
      EXTERNAL FBLW,ULMASS
      DATA SMALL/1.D-6/
      IPHI = 0
      XEPS = 0.01
      SLO = SNGL(SMALL)
      SHI = 1. - SLO
      q = Q2
      beta = bet
      xp = x_pom
      t = t2
      phi = 0.d0
      aem=1.d0/137.d0
      pi=2.d0*dasin(1.d0)
      s = SSS
      pt2min=pt2cut(12)
      xm=q2/beta-q2 + t
      pt2max=xm/4.d0 - 4.d0*DBLE(ULMASS(3))**2
      F2QT = 0.D0
      FLQ = 0.D0
      if(pt2min.gt.pt2max) return
      NI = 25

c      write(6,*) 'F2BLW: beta,xpom,q2,t2',beta,x_pom,q2,t2
      CALL INTGA(SLO,SHI,FBLW,XEPS,RESULT)
c      write(6,*) ' after INTGA : ',RESULT
c          write(6,*) ' f2blw  xng,xpom',,xng,xpom
c      write(6,*) 'f2blw f2qt ',f2qt
      F2QT=DBLE(RESULT)
      FLQ=0.0D0

   10 CONTINUE
c         write(6,*) ' f2blw iwei =',iwei
      if(iwei.eq.1) then
c         write(6,*) ' f2blw iwei =1'
         F2QMAX = 0.D0
         DO 20 I1 = 0,NI
            rn1 = dfloat(I1)/dfloat(NI)
            if(rn1.eq.0.d0) rn1 = small
            if(rn1.eq.1.d0) rn1 = 1.d0 - small
c weighting with 1/pt**4
            pt2 = pt2min*pt2max/(rn1*pt2min + pt2max*(1.D0 - rn1))
            wpt2 = pt2**2*(pt2max-pt2min)/pt2min/pt2max
            if(rn1.gt.0.998d0) goto 20
            CALL SIGBLW(pt2,sigt,sigl)
            F2TEST = sigt*Q2/(4.D0*PI**2*aem)*wpt2
c         write(6,*) ' f2blw 1st: F2TEST,wpt2 ',F2TEST,wpt2,sigt,pt2
            IF(F2TEST.GE.F2QMAX) F2QMAX = F2TEST

   20    CONTINUE
         F2QMAX = 2.D0*F2QMAX
   30    CONTINUE
         iphi = 1
         rn1 = draprn()
c weighting with 1/pt**4
         pt2 = pt2min*pt2max/(rn1*pt2min + pt2max*(1.D0 - rn1))
         wpt2 =  pt2**2*(pt2max-pt2min)/pt2min/pt2max
         if(rn1.gt.0.998d0) goto 30
         rn2 = draprn()
c chi integration (0 - 2*pi)
         PHI = rn2*2.d0*pi

         CALL SIGBLW(pt2,sigt,sigl)
         F2GEN =  sigt*Q2/(4.D0*PI**2*aem)*wpt2
         F2RN = draprn()
         IF(F2GEN.GE.F2QMAX) THEN
            write(6,*) ' F2BLW : F2GEN > F2MAX ',F2GEN,F2QMAX
         ELSEIF(F2GEN.LT.0.D0.OR.F2QMAX.LT.0.D0) THEN
            write(6,*) ' F2BLW: F2GEN OR F2MAX < 0 ',F2GEN,F2QMAX
            write(6,*) ' F2BLW: beta,x_pom,q2,t2',beta,x_pom,q2,t2
            write(6,*) ' F2BLW: sigt,wpt2 ',sigt,wpt2
         ENDIF
c rotation with pt2gen is done in pyremn similar to primordial pt
         IF(F2QMAX*F2RN.GT.F2GEN) THEN
c            write(6,*) 'F2QMAX*F2RN.GT.F2GEN',F2QMAX,F2RN,F2GEN
            GOTO 30
         ENDIF
         PT2GEN = pt2
         PHIGEN = PHI
c         write(6,*) ' F2BLW final : PT2GEN',pt2gen,' phigen ',phigen
      endif
      RETURN
      END
