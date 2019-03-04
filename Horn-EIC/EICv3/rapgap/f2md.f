*CMZ :  2.08/04 22/12/99  15.39.30  by  Hannes Jung
*CMZ :  2.06/34 20/05/98  09.12.56  by  Hannes Jung
*CMZ :  2.06/26 23/01/98  22.07.10  by  Hannes Jung
*CMZ :  2.06/10 16/10/97  10.32.53  by  Hannes Jung
*CMZ :  2.06/02 29/07/97  09.22.42  by  Hannes Jung
*CMZ :  2.05/05 20/03/97  12.36.52  by  Hannes Jung
*CMZ :  2.03/04 05/09/96  17.11.59  by  Hannes Jung
*CMZ :  2.03/03 22/08/96  19.41.59  by  Hannes Jung
*CMZ :  2.03/02 16/08/96  15.14.32  by  Hannes Jung
*CMZ :  2.02/03 01/08/96  09.49.34  by  Hannes Jung
*CMZ :  2.02/02 25/07/96  12.26.14  by  Hannes Jung
*CMZ :  2.02/01 02/07/96  18.57.50  by  Hannes Jung
*CMZ :  2.02/00 02/07/96  17.58.21  by  Hannes Jung
*CMZ :  2.01/19 25/06/96  12.14.36  by  Hannes Jung
*CMZ :  2.01/18 20/06/96  11.21.16  by  Hannes Jung
*CMZ :  2.01/15 20/05/96  07.27.09  by  Hannes Jung
*-- Author :    Hannes Jung   19/05/96

      SUBROUTINE F2MD(BETA,X_POM,Q2,T2,XDX,FLQ)
* calculate  F_2^D in the approach of M. Diehl
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
*KEEP,RGPARAS.
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

*KEEP,RGPARAM.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
*KEEP,RGLUCO.
      REAL PLEPIN,PPIN
      INTEGER KE,KP,KEB,KPH,KGL,KPA,NFRAG,ILEPTO,IFPS,IHF,IALMKT
      INTEGER INTER,ISEMIH
      INTEGER NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT,NFLAV,NFLQCDC
      COMMON/LUCO  /KE,KP,KEB,KPH,KGL,KPA,NFLAV,NFLQCDC
      COMMON/INPU  /PLEPIN,PPIN,NFRAG,ILEPTO,IFPS,IHF,IALMKT,INTER,
     +              ISEMIH
      COMMON/HARD/ NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT
      INTEGER IHFLA
      COMMON/HFLAV/ IHFLA
C      SAVE

*KEND.
      DOUBLE PRECISION mq,mp,mx2,beta0,mu0
      DOUBLE PRECISION PT2GEN,PHIGEN
      DOUBLE PRECISION XDX(-4:4)
      COMMON/HARDPOM/PT2GEN,PHIGEN
      DOUBLE PRECISION PT2MD,PHIMD
      COMMON/PT/PT2MD,PHIMD
      COMMON /PARMD/ bet,xpom,q2t,t2t,mq,eq2
      COMMON /SEL/IPHI
      REAL ULMASS
      EXTERNAL ULMASS,LUCHGE
      REAL SLO,SHI,FMD,XEPS,RESULT
*KEEP,RGLUDAT2.
      REAL PMAS,PARF,VCKM
      INTEGER KCHG
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
C      SAVE

*KEND.
      EXTERNAL FMD
      DATA SMALL/1.D-6/
      IPHI = 0
      pi=2.d0*dasin(1.d0)
      XEPS = 0.01
      SLO = SNGL(SMALL)
      SHI = 1. - SLO
      q2t = q2
      t2t = t2
      bet = beta
      xpom = x_pom
c      write(6,*) ' xpom',x_pom,beta,q2,t2
      mp = 0.938d0
      mx2 = Q2*(1.d0-beta)/beta + t2
c      write(6,*) ' f2md mx2,t2 ',mx2,t2
      eq2 = 2.d0/3.d0
      beta0 = 2.d0
      mu0 = 1.2d0
      mq = 0.3d0
      alphas0 = 1.d0
      NI = 25
      DO 10 I=-4,4
   10 XDX(I)=0.d0
      F2GF = 0.D0
      pt2min=pt2cut(12)
c      write(6,*) ' F2MD ptmin,pt2cut',pt2min,pt2cut(12)
c      write(6,*) 'F2MD: beta,xpom,q2,t2,f1t',beta,x_pom,q2,t2
      F2QMAX = 0.D0
      DO 20 IP=1,NFLAV
c current quark masses
c         mq=dble(ulmass(ip))
c constituent quark masses
         mq=dble(PARF(100+IP))
         eq2 = dble(LUCHGE(IP))**2/9.0d0
c         write(6,*) ' mass ',mq,eq2,ip
         F2QTF =0.d0
         CALL INTGA(SLO,SHI,FMD,XEPS,RESULT)
c        write(6,*) ' after INTGA : ',RESULT
         F2QTF=DBLE(RESULT)
         XDX(IP)=F2QTF/eq2/2.D0
   20 CONTINUE
      XDX(-4)=XDX(4)
      XDX(-3)=XDX(3)
      XDX(-2)=XDX(2)
      XDX(-1)=XDX(1)
      XDX(0) = 0.d0
c      write(6,*) ' F2MD:F2QTF ',F2QTF,XDX
c      write(6,*) ' F2MD final : PT2GEN,v,beta,q2',pt2gen,v,beta,q2
      if(iwei.eq.1) then
         F2QMAX = 0.D0
         q2t = q2
         t2t = t2
         bet = beta
         xpom = x_pom
c      write(6,*) ' xpom',x_pom,beta,q2,t2
         mp = 0.938d0
         mx2 = Q2*(1.d0-beta)/beta + t2
         DO 30 I1 = 0,NI
            ip = KPA
c current quark masses
c           mq=dble(ulmass(ip))
c constituent quark masses
            mq=dble(PARF(100+IP))
            eq2 = dble(LUCHGE(IP))**2/9.0d0
            rn2 = dfloat(I1)/dfloat(NI)
            if(rn2.eq.0.d0) rn2 = small
            if(rn2.eq.1.d0) rn2 = 1.d0 - small


            F2QT = dble(FMD(SNGL(rn2)))
c         FLQ =  sigl * wpt2
            F2TEST = F2QT
c         write(6,*) ' f2md 1st: F2QT FLQ ',F2QT,FLQ,rn2
c         write(6,*) ' pt2,wpt2,sigt',pt2,wpt2,sigt
            IF(F2TEST.GE.F2QMAX) F2QMAX = F2TEST
   30    CONTINUE
         F2QMAX = 2.0D0*F2QMAX

   40    CONTINUE
         IPHI = 1
         rn1 = draprn()
         rn2 = draprn()
c chi integration (0 - 2*pi)
         PHIMD = rn1*2.d0*pi
         if(rn2.gt.0.998d0) goto 40
         F2GEN = dble(FMD(SNGL(rn2)))
         PT2GEN = PT2MD
         PHIGEN = PHIMD
         F2RN = draprn()
c            write(6,*) ' F2MD: KPA = ',IP
c            write(6,*) ' F2MD : F2GEN F2MAX ',F2GEN,F2QMAX
c            write(6,*) ' f2MD: pt2,pl2,chi ',pt2,pl2,chi
         IF(F2GEN.GE.F2QMAX) THEN
            write(6,*) ' F2MD : F2GEN > F2MAX ',F2GEN,F2QMAX
         ELSEIF(F2GEN.LT.0.D0.OR.F2QMAX.LT.0.D0) THEN
            write(6,*) ' F2MD : F2GEN OR F2MAX < 0 ',F2GEN,F2QMAX
            write(6,*) ' F2MD : beta,x_pom,q2,t2',beta,x_pom,q2,t2

         ENDIF
c rotation with pt2gen is done in pyremn similar to primordial pt
         IF(F2QMAX*F2RN.GT.F2GEN) GOTO 40
c        write(6,*) ' F2MD final : PT2GEN,PHIGEN',pt2gen,phigen
      endif
      RETURN
      END
