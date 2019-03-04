*CMZ :  2.08/04 22/12/99  15.39.28  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.58  by  Hannes Jung
*CMZ :  2.07/03 16/05/99  09.40.42  by  Hannes Jung
*-- Author :    Hannes Jung   12/01/95
      SUBROUTINE partdi99(X,WPART)
      IMPLICIT NONE
*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGLUJETS.
      INTEGER N,K
      REAL SP,V
      DOUBLE PRECISION P
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
      REAL ULMASS
      DOUBLE PRECISION DOT1,DPLU,DLANGL
      EXTERNAL DLANGL,ULMASS,DPLU,LUCHGE,DOT1
C      SAVE

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


*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGSTRU.
      REAL SCAL1,XPD1,SCAL2,XPD2
      COMMON/STRU/ SCAL1,XPD1,SCAL2,XPD2
C     SAVE
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
      Double Precision PHI
      COMMON/DIFFA/ PHI
      Double Precision draprn
      DOUBLE PRECISION ME,MP
      Double Precision X,XY,WPART,WEIGHT
      COMMON/XVAR/ XY(10)
      DOUBLE PRECISION STHETA,SPHI
      DIMENSION X(20)
      REAL XPQ(-6:6)
      REAL F2CC(-6:6),XF3CC(-6:6)
      Double Precision WMAX
      Integer NDIMC,IMIX,IDEBUG,NDIM,NPOIN,LST,IRES,NTHRE,NVMRE,IGENFL
      COMMON /DIMEN/ NDIMC
      COMMON /OALPHAS/ WMAX,IMIX
      COMMON/ INTERN/IDEBUG
      COMMON/DIVO/ NDIM,NPOIN
      COMMON/EPPARA/LST(30),IRES(2)
      COMMON/HERTHET/ NTHRE,NVMRE
      COMMON/GENWEI/IGENFL
      DOUBLE PRECISION QG2
      COMMON/SEMIH/ QG2
      REAL SNGL
      EXTERNAL draprn
      Double Precision QF,QFT,XR12,SMALL,ALPH_EM,GF,DOT,WTGLU
      Integer ICHECK,NFLAVP,LUCHGE,KPHF,NACC
      Double Precision  WT,W02,W12,YX,SMIN,XP2Q,XG1,XP1MIN
      Double Precision QG2MAX,QG2MIN,QG20,XZOHRN,WQG2
      Double Precision  XP2,SHAT1,PIO,FGAM,FWEI,FLUX
      Double Precision THETE,ECM,PT2,PT,WTG15,XP1,PR,t,at
      Double Precision COSTP,PHIP,PHIO,SPHP,CPHP,STHP,CTHP,BOCHCK
      Double Precision PEP,PEG,PEGZ,POMDGA,PEZ,EN,PZC,WMAT,XPINT,XP
      Integer I,J,IN,KI,MSTJ24,NP,NPP,KPA2,KPAT,KPFL,KPAO,KPAO2
      Integer NDIMEN,NRN,IST,NPFIN,NB2,nafl,KPF,KFLCH,KFLSP


      DATA QF/0.0D0/,W12/0.0D0/
      NDIMEN = NDIM
c      IDEBUG=1
c      write(6,*) ' in partdi99 ',x
      SCAL1 = -99999.
      XPD1 = -99999.
      SCAL2 = -99999.
      XPD2 = -99999.
      DO 10  IN=1,20
         K(IN,1) = 0
         K(IN,2) = 0
   10 CONTINUE

C...  GIVE BEAM  FOUR VECTORS
      DO 20 I=1,2
         DO 20 J=1,5
            K(I,J)=KBEAM(I,J)
   20 P(I,J) = PBEAM(I,J)
      ME =P(1,5)
      MP =P(2,5)
      NP=2
      NPP=3*NP
      N=2

      QF = 0.D0
      NFLAVP = 0

      IF(IPRO.EQ.10.OR.IPRO.EQ.12.OR.IPRO.EQ.13.
     +OR.IPRO.EQ.15.OR.IPRO.EQ.18) THEN
         AM(1) = 0.0D0
         AM(2) = 0.0D0
         IF(IHFLA.GE.4.AND.IPRO.EQ.18) THEN
            KPA=IHFLA
            AM(1) = DBLE(ULMASS(KPA))
            AM(2) = DBLE(ULMASS(KPA))
         ENDIF
      ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
         IF(IWEI.EQ.0) THEN
            KPA = 4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         AM(1) = DBLE(ULMASS(KPA))
         AM(2) = DBLE(ULMASS(KPA))
      ELSEIF(IPRO.EQ.99) THEN
         AM(1) = 0.0D0
         AM(2) = 0.0D0
         AM(3) = 0.0D0
         NP=3
         NPP=3*NP

      ELSE
         WRITE(6,*) ' partdi99: wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
      IF(IPRO.EQ.14.AND.IMIX.EQ.1) THEN
      ELSE
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
         W02=(AM(1)+AM(2)+MP)**2
         IF(IPRO.EQ.12) W02=(AM(1)+MP)**2
         IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
         IF(AM(1).LT.1.0D0.and.pt2cut(ipro).ne.0.)
     +     W02=(DSQRT(PT2CUT(IPRO)) + MP)**2
         W12=W02-MP*MP
         YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
         YMAX=YMAX/(2.D0*(SSS+ME*ME))
         YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
         YMIN=YMIN/(2.D0*(SSS+ME*ME))
         IF(YMI.GT.YMIN) YMIN=YMI
         IF(YMA.LT.YMAX) YMAX=YMA
      ENDIF
c         write(6,*) w02
c         WRITE(6,10500) YMIN,YMAX
c10500 FORMAT(' limits on y ',/,
c     +' YMIN = ',E10.5,' YMAX = ',E10.5)

C ... select particle code for light flavour production according
C ... to charges
      IF(IPRO.EQ.10.OR.IPRO.EQ.13.OR.IPRO.EQ.18) THEN

      ELSEIF(IPRO.EQ.99) THEN
         NFLAVP = NFLAV
         KPA = 1
      ENDIF
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSEIF(IGENFL.EQ.1) THEN
         KPA = KPAO
         KPA2 = KPAO2
         IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            AM(1)=DBLE(ULMASS(KPA))
            AM(2)=DBLE(ULMASS(KPA2))
         ENDIF
      ENDIF

C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... XP2 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,XP2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NRN = 0
c      IRES(1)=0
      XMAX=0.999d0
      IF(IRES(1).EQ.0.OR.IPRO.EQ.12) THEN
c.......
         IF(KE.NE.22) THEN
            NRN = NRN + 1

            YX = YMIN*((YMAX/YMIN)**X(NRN))
c            write(6,*) ' y,ymax,ymin,x(NRN),NRN',yx,ymax,ymin,x(nrn),nrn
            IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 70
            Q2MIN=ME*ME*YX*YX/(1.D0-YX)
            IF(QMI.GT.Q2MIN) Q2MIN = QMI
            Q2MAX=YX*SSS - W12
            IF(QMA.LT.Q2MAX) Q2MAX = QMA
            IF(Q2MAX.LT.Q2MIN) GOTO 80
            NRN = NRN + 1
            IF(QMI.EQ.0.D0) THEN
               Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
            ELSE
c new try weighting with 1/q**4
               Q2 = Q2MIN*Q2MAX/(X(NRN)*Q2MIN + Q2MAX*(1.D0 - X(NRN)))
            ENDIF
            XMIN = 0.0D0
            XP2Q=Q2/YX/SSS
            QG2 = 0.D0
            WQG2 = 1.D0
            IF(IPRO.NE.12) THEN
               SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
               XMIN=(SMIN+Q2+QG2)/(YX*SSS)
               IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
                  XMIN = DMAX1(XMIN,Q2/YX/SSS)
c                    XMAX = XR12
               ENDIF
            ELSE
               XMIN=Q2/YX/SSS
            ENDIF
         ELSEIF(KE.EQ.22) THEN
            Q2 = 0.D0
            YX = 1.D0
            QG2 = 0.D0
            WQG2 = 1.D0
c               write(6,*) ' partdi99 yx,q2 ',yx,q2
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            XMIN= (SMIN+QG2)/(YX*SSS)
         ELSE
            write(6,*) ' wrong KF selected; program stop '
            YX = 0.D0
         ENDIF
      ELSE
         YX = 0.D0
         XG1 = 0.D0
      ENDIF
      XG1=YX
      XMIN=Q2/YX/SSS
      IF(XMIN.GE.XMAX) GOTO 90
C ... XP2 = E_parton /E_proton

      IF(IPRO.NE.12.AND.IMIX.EQ.0) THEN
         NRN = NRN + 1
         XP2= XMIN*((XMAX/XMIN)**X(NRN))
      ELSEIF(IPRO.NE.12.AND.IMIX.EQ.1.AND.IWEI.EQ.1) THEN
         XP2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
      ELSEIF(IPRO.EQ.12.OR.IPRO.EQ.100) THEN
         XP2 = Q2/YX/SSS
      ELSE
         XP2 = 0.D0
      ENDIF
c      write(6,*) ' xp2 ',xp2
      IF(SNGL(XP2).GE.0.9999.OR.XP2.LT.XMIN) GOTO 100
      SHAT1 = SSS*XG1*XP2
      IF(IPRO.NE.12)  SHAT1=SSS*XG1*XP2 - Q2 - QG2
c      write(6,*) ' sss,xg1,xp2,q2,qg2 ',sss,xg1,xp2,q2,qg2
      IF(SHAT1.LT.(AM(1)+AM(2))**2) GOTO110

      IST = 0
      YY=SNGL(YX)
      XEL=SNGL(XG1)
      XPR=SNGL(XP2)
      PHI = 99999.D0
      NRN=NRN+1
      PHI = 2.D0*PI*X(NRN)
      CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      XP=0.D0
      IRES(2)=1
C here do some gymnastics to make cut on theta angle of electron
      IF(THEMA.LT.180.D0.OR.THEMI.GE.0.D0) THEN
C go first to ep LAB system
         CALL DUDBRB(3,3,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
         THETE = DPLU(3,14)
         CALL DUDBRB(3,3,0.D0,0.D0,-CM(1)/CM(4),
     +                -CM(2)/CM(4),-CM(3)/CM(4))
c         CALL DULIST(1)
         IF(THETE .GT. 180.01D0) THEN
            WRITE(6,*) ' FATAL: theta_electron > 180 deg ',THETE
         ENDIF
         IF(THETE .GE. THEMA) GOTO 120
         IF(THETE .LE. THEMI) GOTO 120
      ENDIF
C end of these gymnastics

C FINAL STATE PROTON
      NPFIN=NIA1+4
      IF(IPRO.EQ.99) THEN
         NPFIN=NIA1+5
      ENDIF
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      K(NPFIN,3)=2
      P(NPFIN,5)=DBLE(ULMASS(KP))
      PR = P(2,3)*(1.D0 - XP2)
      P(NPFIN,4) = PR
      CPHP=1.D0
      SPHP=DSQRT(1.D0 - CPHP**2)
	STHP=1.D0
      IF(DABS(STHP).GT.1.D0) goto 180
      CTHP=DSQRT(1.D0 - STHP**2)
      P(NPFIN,1)= PR*STHP*CPHP
      P(NPFIN,2)= PR*STHP*SPHP
      P(NPFIN,3)= PR*CTHP
c      CALL DULIST(1)
C MOMENTA OF PARTON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KGL
      IF(IPRO.EQ.12) K(NIA1+1,2) = KPA
      K(NIA1+1,3)=2
      DO 30  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   30 CONTINUE
      IF(ISEMIH.EQ.0) THEN
         P(NIA1+1,4)=ABS(P(NIA1+1,3))
         P(NIA1+1,5)=0.0D0
      ELSE
         P(NIA1+1,4)=DSQRT(P(NIA1+1,1)**2+P(NIA1+1,2)**2+P(NIA1+1,3)**2
     +               - QG2)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
      ENDIF
      NIA2 = NIA1+1
c      write(6,*) ' partdi99 p(nia2,5)**2',DOT1(NIA2,NIA2)
c      write(6,*) ' partdi99 p(npfin,5)**2',DOT1(NPFIN,NPFIN)
c      write(6,*)  ' partdi99 '
c      call dulist(1)
c      pause
      IF(IPRO.NE.12) THEN
         NF1=NIA1+2
         NF2=NIA1+3
         K(NF1,1)=2

         IF(IPRO.EQ.15) THEN

         ELSE
            K(NF2,2)=KPA2
         ENDIF
         IF(IPRO.NE.15.AND.IPRO.NE.12) KPF=KPA
         K(NF1,2)=KPF
         K(NF1,3)=NIA1
         K(NF2,1)=1
         K(NF2,3)=NIA1
         P(NF1,5)=AM(1)
         P(NF2,5)=AM(2)
         NB2 = NIA2
         IF(IPRO.EQ.99) THEN
            NF1=NIA1+2
            NF2=NIA1+4
            K(NF1,1)=2
            K(NF1+1,1)=2
            K(NF2,1)=1
            K(NF1,2)=KPF
            K(NF1+1,2)=21
            K(NF2,2)=-KPF
            K(NF1,3)=NIA1
            K(NF1+1,3)=NIA1
            K(NF2,3)=NIA1
         ENDIF
C...   VECTOR OF GAMMA GLUON CMS SYSTEM

      ELSE
         NB2 = 0
      ENDIF
      DBCMS(1)=  P(NIA1,1) + P(NB2,1)
      DBCMS(2)=  P(NIA1,2) + P(NB2,2)
      DBCMS(3)=  P(NIA1,3) + P(NB2,3)
      DBCMS(4)=  P(NIA1,4) + P(NB2,4)
      DO 40  I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
   40 CONTINUE
      SHAT=DOT(DBCMS,DBCMS)
      IF(SHAT.LE.0.0) THEN
         GOTO 140
      ENDIF
c          call dulist(1)
C NOW BOOST TO GAMMA GLUON
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2 +
     +(DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 130
      CALL DUDBRB(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +-DBCMS(3)/DBCMS(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      call DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      call DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
      IF(IPRO.NE.12) THEN
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
c         write(6,*) ECM
         ECM =DSQRT(SHAT)
C      IF(ECM.LE.(AM(1)+AM(2))) WRITE(6,*) ' ECM LE MASS',ECM
         IF(ECM.LE.(AM(1)+AM(2))) GOTO 150
c         write(6,*) ' IPRO = 13 ,NRN= ',NRN
         IF(IPRO.EQ.99) THEN
            NRN=NRN+1
            XY(1)=X(NRN)
            NRN=NRN+1
            XY(2)=X(NRN)
            NRN=NRN+1
            XY(3)=X(NRN)
            NRN=NRN+1
            XY(4)=X(NRN)
            NRN=NRN+1
            XY(5)=X(NRN)
         ELSE
            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XY(1)=X(NDIMEN+3)
               XY(2)=X(NDIMEN+4)
            ELSE
               NRN = NRN + 1
               XY(1) = X(NRN)
               NRN = NRN + 1
               XY(2) = X(NRN)
            ENDIF
         ENDIF
c         write(6,*) am(1),am(2),ipro
c         write(6,*) 'partdi99 ',XY(1),XY(2),NDIM,NDIMEN,IMIX,IWEI
         CALL PHASE(NP,ECM,AM,PCM,WT)
         IF(WT.LE.0.D0) GOTO 160
         IF(IPRO.EQ.99) THEN
            DO 50 I=1,4
               P(NF1,I)=PCM(I,1)
               P(NF1+1,I)=PCM(I,2)
               P(NF2,I)=PCM(I,3)
               P(NF1,5)=AM(1)
               P(NF1+1,5)=AM(2)
               P(NF2,5)=AM(3)
   50       CONTINUE
         ELSE
            DO 60 I=1,4
               P(NF1,I)=PCM(I,1)
               P(NF2,I)=PCM(I,2)
   60       CONTINUE
         ENDIF
c         write(6,*) ' partdi99 ', KPA,KPA2,K(NF1,2),K(NF2,2)
         PT2 = DPLU(NF1,9)
         CALL CUTG(PT2,NACC)
         IF(NACC.EQ.0) GOTO 170
cscale        IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
cscale             Q2Q = (AM(1)+AM(2))**2
cscale            Q2Q = (AM(1)+AM(2))**2 + Q2
cscale         ELSE
         IF(IQ2.EQ.1) THEN
            Q2Q = (2.D0*AM(1))**2
         ELSEIF(IQ2.EQ.2) THEN
            Q2Q = SHAT
         ELSEIF(IQ2.EQ.3) THEN
            Q2Q = (2.D0*AM(1))**2 + PT2
         ELSEIF(IQ2.EQ.4) THEN
            Q2Q = Q2
         ELSEIF(IQ2.EQ.5) THEN
            Q2Q = Q2 + PT2 + (2.D0*AM(1))**2
c this is for testing
         ELSEIF(IQ2.EQ.6) THEN
c torbjorn suggestion
            Q2Q = PT2 + 8.D0*PT2**2/SHAT
C kt**2 as from Zeppenfeld/Mirkes
            t = AM(1)**2 - 2.D0 * DOT1(NIA1,NF1) - P(NIA1,5)**2
            at=DSQRT((t+Q2)**2+4.D0*Q2*PT2)
            Q2Q = at*(at+t+Q2)/2.D0/Q2
c         write(6,*) 'partdi99 : t,PT2,Q2,Q2Q ',t,PT2,Q2,Q2Q
         ELSEIF(IQ2.EQ.7) THEN
c torbjorn suggestion
            Q2Q = PT2 + 8.D0*PT2**2/SHAT
         ELSE
            WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
            STOP
         ENDIF
cscale         ENDIF
         IF(IQ2.EQ.5) THEN
            Q2Q = Q2 + PT2*SCALFA + (2.D0*AM(1))**2
         ELSE
            Q2Q = SCALFA*Q2Q
         ENDIF

C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
         call DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
         call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
         XPR = SNGL(XP2)
c         IDIRO = IDIR
c         IDIR = 0
         if(xp2.GT.1.d0) write(6,*) 'partdi99 xp2>1',xp2,ipro
c         write(6,*) 'partdi99 xp2,xg1',xp2,xg1
         CALL PYSTFU(K(2,2),SNGL(XP2),SNGL(Q2Q),XPQ)

c         IDIR = IDIRO
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
c for resolved photon move parton densities into ELERES
         IF(IRES(1).EQ.1) THEN
            XP = 1.D0/XP2
         ELSE
c            write(6,*) ' partdi99 ires=0'
            XP = DBLE(XPQ(0))/XP2
         ENDIF
         SCAL2 = SNGL(Q2Q)
         XPD2 = XPQ(0)
         XPINT=DLOG(XMAX/XMIN)
         WTGLU = XP2 * XP * XPINT

         IF(ISEMIH.EQ.1.AND.IPRO.NE.12.AND.IPRO.NE.15)
     +      WTGLU = WTGLU*WQG2
         WMAT=WT *WTGLU*FGAM*FWEI
         WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
         FLUX=(2.D0*PI)**4
         FLUX=FLUX/(2.D0*(SHAT+Q2))
         WPART = WMAT *FLUX
         IF(IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.99) THEN
            WPART = WTGLU * FWEI * WT
            IF(IPRO.EQ.15) THEN
               WPART = WPART*FLUX/((2.D0*PI)**NPP)
            ENDIF

         ENDIF
c         IF(WPART.LE.0) THEN
c           write(6,*) ' partdi99: SHAT,Q2 ',SHAT,Q2
c           write(6,*) ' partdi99: WTGLU,FWEI,WT,FGAM',WTGLU,FWEI,WT,FGAM
c           write(6,*) ' partdi99: XP2,XPINT,XP,WQG2 ',XP2,XPINT,XP,WQG2
c           call dulist(1)
c           ENDIF
         PT2H = SNGL(PT2)
         SHH = SNGL(SHAT)
c           call dulist(1)


      ENDIF
      NDIMC = NRN
c      write(6,*) ' end of partdi99 ',NDIMC,NRN
c      CALL DULIST(1)
c      write(6,*)  'partdi99: WPART ',WPART
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSE
         KPA = KPAO
         KPA2 = KPAO2
      ENDIF
c      write(6,*) ' partdi99 end ,kpa,kpa1,kpa2,kpao2',kpa,kpa1,kpa2,kpao2
      RETURN
   70 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      RETURN
   80 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
         NDIMC = 9999
      ENDIF
      RETURN
   90 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: xlimit ; RETURN xmin,xmax',xmin,xmax
      ENDIF
      NDIMC = 9999
      RETURN
  100 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: xp2imit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  110 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: slimit ; RETURN ',shat1,am(1),am(2),IPRO
         write(6,*) ' sss,xg1,xp2,q2,qg2 ',sss,xg1,xp2,q2,qg2
         write(6,*) ' partdi99 q2,y,sss',q2,yx,sss
      ENDIF
      NDIMC = 9999
      RETURN
  120 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: theta limit ; RETURN '
      ENDIF
      NDIMC = 9999
      IF(IHERAC.EQ.1) THEN
         NTHRE = NTHRE +1
      ENDIF
      RETURN
  130 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: bochck limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  140 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: shat ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  150 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: ECM limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  160 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: PHASE WT=0 ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN

  170 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: PTCUT limit ; RETURN ',PT2,KPA
      ENDIF
      NDIMC = 9999
      RETURN
  180 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi99: QG2 limit ; RETURN ',DSQRT(QG2),PR,STHP
      ENDIF
      NDIMC = 9999
      RETURN


      END
