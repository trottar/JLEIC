*CMZ :  2.08/04 22/12/99  15.39.29  by  Hannes Jung
*CMZ :  2.06/34 20/05/98  09.12.56  by  Hannes Jung
*CMZ :  2.06/29 15/03/98  14.20.39  by  Hannes Jung
*CMZ :  2.06/27 01/02/98  12.55.36  by  Hannes Jung
*CMZ :  2.06/26 21/01/98  17.46.26  by  Hannes Jung
*CMZ :  2.06/24 14/01/98  19.40.49  by  Hannes Jung
*CMZ :  2.05/05 20/03/97  12.01.07  by  Hannes Jung
*CMZ :  2.04/00 03/12/96  09.39.16  by  Hannes Jung
*CMZ :  2.03/04 08/09/96  14.54.40  by  Hannes Jung
*CMZ :  2.03/03 19/08/96  17.22.58  by  Hannes Jung
*CMZ :  2.02/03 01/08/96  09.33.47  by  Hannes Jung
*CMZ :  2.02/01 02/07/96  18.53.24  by  Hannes Jung
*CMZ :  2.01/18 19/06/96  11.40.39  by  Hannes Jung
*CMZ :  2.01/17 05/06/96  10.36.55  by  Hannes Jung
*CMZ :  2.01/16 23/05/96  12.18.53  by  Hannes Jung
*CMZ :  2.01/15 14/05/96  16.21.03  by  Hannes Jung
*CMZ :  2.01/14 13/05/96  15.59.42  by  Hannes Jung
*CMZ :  2.01/09 07/03/96  13.50.54  by  Hannes Jung
*-- Author :
      SUBROUTINE F2DHW(BETA,X_POM,Q2,T2,F2QTF,FLQF,XGX)
* calculate  F_2^D in the hard approach of M. Wuesthoff
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
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
      EXTERNAL POM,FG
      DATA SMALL/1.D-6/
      q02=1.0d0
      NTRY = 100
      NI = 5
      F2QTF = 0.D0
      FLQF = 0.D0
      XGX = 0.D0
c      write(6,*) 'F2DHW: beta,xpom,q2,t2',beta,x_pom,q2,t2
c      write(6,*) ' f2dhw iwei',iwei
      xi0=log(4.0d0*q02/q2*beta)
      xi1=12.d0
      F2QMAX = 0.D0
      DO 10 IRN = 1,NTRY
         rn1 = draprn()
         v=exp(xi0+rn1*(xi1-xi0))
         wv=(xi1-xi0)*v
c
         rn2 = draprn()
         z=beta/(beta+rn2*(1.d0-beta))
         wz = z**2 *(1.d0 -beta)/beta
         CALL XNTOTAL(F2QT,F2G,F2QL,Q2,x_pom,BETA,V,Z,T2)
         F2QTF = F2QTF + F2QT * wv/NTRY + F2G * wv * wz/NTRY
         FLQF =FLQF + F2QL * wv/NTRY
c get the gluon
         xng = 0.d0
c
c
         xpom=pom(x_pom,-t2,v)
         xpom = xpom/x_pom
         if (v.lt.q02/q2) then
            xng=0.0d0
         else
            xng=fg(beta,v)/v
            xng=xng**2/q02
            xng=xng*9.0d0/4.0d0
c
            xng = xng/16.d0 * xpom
         endif
         xgx = xgx + xng*wv/ntry

   10 CONTINUE
c          write(6,*) ' f2dhw f2qtf xng,xpom',f2qtf,xng,xpom

c      write(6,*) 'f2dhw : ',f2gf
      DO 20 I1 = 0,NI
         DO 20 I2 = 0,NI
            rn1 = dfloat(I1)/dfloat(NI)
            if(rn1.eq.0.d0) rn1 = small
            if(rn1.eq.1.d0) rn1 = 1.d0 - small
            v=exp(xi0+rn1*(xi1-xi0))
            wv=(xi1-xi0)*v
c
            rn2 = dfloat(I2)/dfloat(NI)
            if(rn2.eq.0.d0) rn2 = small
            if(rn2.eq.1.d0) rn2 = 1.d0 - small
            z=beta/(beta+rn2*(1.d0-beta))
            wz = z**2 *(1.d0 -beta)/beta
            CALL XNTOTAL(F2QT,F2G,F2QL,Q2,x_pom,BETA,V,Z,T2)
            F2TEST = (F2QT+F2QL)*wv + F2G*wv*wz
c         write(6,*) ' f2dhw 1st: (F2g) ',F2g,z,v,wv,wz
c         write(6,*) ' f2dhw 1st: (F2QT+F2QL),wv ',F2TEST,(F2QT+F2QL)*wv
            IF(F2TEST.GE.F2QMAX) F2QMAX = F2TEST

   20 CONTINUE
      F2QMAX = 5.D0*F2QMAX
   30 CONTINUE
      if(iwei.eq.1) then
         rn1 = draprn()
         v=exp(xi0+rn1*(xi1-xi0))
         wv=(xi1-xi0)*v
         rn2 = draprn()
         z=beta/(beta+rn2*(1.d0-beta))
         wz = z**2 *(1.d0 -beta)/beta
         CALL XNTOTAL(F2QT,F2G,F2QL,Q2,x_pom,BETA,V,Z,T2)
         F2GEN = (F2QT+F2QL)*wv + F2G*wv*wz
         F2RN = draprn()
         IF(F2GEN.GE.F2QMAX) THEN
            write(6,*) ' F2DHW : F2GEN > F2MAX ',F2GEN,F2QMAX
         ELSEIF(F2GEN.LE.0.D0.OR.F2QMAX.LE.0.D0) THEN
            write(6,*) ' F2DHW: F2GEN OR F2MAX < 0 ',F2GEN,F2QMAX
            write(6,*) ' F2DHW: beta,x_pom,q2,t2',beta,x_pom,q2,t2
         ENDIF
c rotation with pt2gen is done in pyremn similar to primordial pt
         IF(F2QMAX*F2RN.GT.F2GEN) GOTO 30
         FT = F2G*wv*wz/F2GEN
         PT2TT = (1.d0 - beta)*q02/v
         IF(FT.GE.draprn()) THEN
            PT2GEN = (1.d0 - z)*q02/v
         ELSE
            PT2GEN = (1.d0 - beta)*q02/v
         ENDIF
c      write(6,*) ' F2DHW final : PT2GEN,v,beta,q2,z',pt2gen,v,beta,q2,z
c      write(6,*) ' F2DHW final : PT2GEN(beta),PT2GEN',pt2tt,pt2gen
c      write(6,*) ' F2DHW final :F2GEN,F2G*wv*wz',F2GEN,F2G*wv*wz
      endif
      RETURN
      END
