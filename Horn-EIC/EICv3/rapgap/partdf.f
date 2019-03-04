*CMZ :  2.08/05 10/03/2000  11.32.27  by  Hannes Jung
*CMZ :  2.08/04 22/12/99  15.39.25  by  Hannes Jung
*CMZ :  2.08/02 12/08/99  16.34.55  by  Hannes Jung
*CMZ :  2.08/01 22/06/99  16.58.44  by  Hannes Jung
*CMZ :          02/06/99  12.39.58  by  Hannes Jung
*CMZ :  2.07/03 16/05/99  09.27.29  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94


      SUBROUTINE PARTDF(X,WPART)
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

*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
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
*KEEP,RGFULL.
      INTEGER IFULL,IQCDGRID
      COMMON /OALPINI/ IFULL,IQCDGRID
*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEEP,RGRAHER.
      REAL XPQDIF,XPQPI
	Integer IHERPYS
	Integer NBQ2,NBX
      PARAMETER (NBQ2=20)
      PARAMETER (NBX=20)
      COMMON /RAHER/ IHERPYS,XPQDIF(-6:6,NBX,NBQ2),XPQPI(-6:6,NBX,NBQ2)
*KEND.
      Double Precision PHI
      COMMON/DIFFA/ PHI
      Double Precision draprn
      DOUBLE PRECISION ME,MP
      Double Precision XY,X,WPART,WEIGHT
      COMMON/XVAR/ XY(10)
      DOUBLE PRECISION STHETA,SPHI
      DIMENSION X(20)
      Integer IJOIN,NJOIN
      DIMENSION IJOIN(10)
      REAL XPQ(-6:6)
      REAL F2CC(-6:6),XF3CC(-6:6)
      Integer NDIMC
      COMMON /DIMEN/ NDIMC
      Integer IMIX
      Double Precision WMAX
      COMMON /OALPHAS/ WMAX,IMIX
      Integer NDIM,NPOIN
      COMMON/DIVO/ NDIM,NPOIN
      Integer LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
      Integer IVM
      COMMON/VMESON/IVM
      Integer IGENFL
      COMMON/GENWEI/IGENFL
      Integer NTHRE,NVMRE
      COMMON/HERTHET/ NTHRE,NVMRE
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      EXTERNAL draprn
      Double Precision QF,QFT,XR12,SMALL,ALPH_EM,GF,DOT
      Integer INRNCH,ICHECK,NFLAVP,LUCHGE,NACC
      Double Precision  WT,W02,W12,YX,SMIN,XP2Q,XG1,XP1MIN,QG2
      Double Precision  XP2,SHAT1,XMCHECK,PIO,FGAM,FWEI,FLUX
      Double Precision THETE,ECM,PT2,PT,WTG15
      Double Precision XP,XR,XX,VMAX,VMIN,XMAXV,XMINV,XMINHF
      Double Precision T2,T2MIN,T2MN,T2MAX1,XR1,WTDIST,WTD12
      Double Precision COSTP,PHIP,PHIO,SPHP,CPHP,BOCHCK
      Double Precision PEP,PEG,PEGZ,POMDGA,PEZ,EN,PZC,WMAT,XPINT
      Double Precision SPOM,WTGLU,tdraprn,XP1
      Double Precision xmax1
      Integer I,J,IN,KI,MSTJ24,NP,NPP,NDFF,KPA2,KPAT,KPFL,KPAO,KPAO2
      Integer NDIMEN,NRN,IST,NPFIN,NB2,nafl,KPF
      DATA QF/0.0D0/,W12/0.0D0/,NFLAVP/0/,XR12/-999./
      DATA SMALL/1.D-3/
      DATA ICHECK/0/
      Data INRNCH/0/
c      IDEBUG = 1
      SCAL1 = -99999.
      XPD1 = -99999.
      SCAL2 = -99999.
      XPD2 = -99999.
      MSTJ24 = MSTJ(24)
      MSTJ(24) = 0
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
      NDFF = 0


      IF(IPRO.EQ.10.OR.IPRO.EQ.12.OR.IPRO.EQ.13.OR.
     +   IPRO.EQ.15.OR.IPRO.EQ.18) THEN
         AM(1) = 0.0D0
         AM(2) = 0.0D0
         IF(IHFLA.GE.4.AND.IPRO.EQ.18) THEN
            KPA=IHFLA
            AM(1) = DBLE(ULMASS(KPA))
            AM(2) = DBLE(ULMASS(KPA))
         ENDIF
         IF(IPRO.EQ.12.OR.IPRO.EQ.15) THEN
            NFLAVP = NFLAV
            IF(NPOM.EQ.2.AND.NG.EQ.11) THEN
c this is for DL pomeron and DL parton density
               NFLAVP=3
            ENDIF
         ENDIF
      ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
         IF(IWEI.EQ.0) THEN
            KPA = 4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         AM(1) = DBLE(ULMASS(KPA))
         AM(2) = DBLE(ULMASS(KPA))

      ELSEIF(IPRO.EQ.100) THEN
         AM(1)=1.D0
         AM(2)=2.D0*DBLE(ULMASS(211))

ccc            KPA = 213
      ELSE
         WRITE(6,*) 'wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
      IF(IPRO.EQ.14.AND.IMIX.EQ.1) THEN
      ELSE
         W02=(AM(1)+AM(2)+MP)**2
         IF(IPRO.EQ.12) W02=(AM(1)+MP)**2
         IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
         W12=W02-MP*MP
         YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
         YMAX=YMAX/(2.D0*(SSS+ME*ME))
         YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
         YMIN=YMIN/(2.D0*(SSS+ME*ME))
         IF(YMI.GT.YMIN) YMIN=YMI
         IF(YMA.LT.YMAX) YMAX=YMA
      ENDIF
c         WRITE(6,10000) YMIN,YMAX
c 10000 FORMAT(' limits on y ',/,' YMIN = ',E10.5,' YMAX = ',E10.5)

C ... select particle code for light flavour production according
C ... to charges
      IF(IPRO.EQ.10.OR.IPRO.EQ.13.OR.IPRO.EQ.18) THEN
         IF(IHFLA.GE.4.AND.IPRO.EQ.18) THEN
            IF(IWEI.EQ.0) THEN
               KPA=4
               IF(IHFLA.GE.4) KPA=IHFLA
            ELSEIF(IWEI.EQ.1) THEN
               IF(draprn().GT.0.5) KPA = -KPA
            ENDIF
            KPA2 = - KPA

         ELSEIF(INTER.LT.2) THEN
            QF=DFLOAT(LUCHGE(1))**2 + DFLOAT(LUCHGE(2))**2 + DFLOAT(
     +      LUCHGE(3))**2
            QF = 2.D0*QF
            KPA = -4
            QFT = - draprn()*QF
   30       KPA=KPA+1
            QFT = QFT + DBLE(LUCHGE(KPA))**2
            IF(QFT.LT.0.0D0) GOTO 30
            IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!! ',
     +      KPA
            KPA2 = - KPA
         ELSEIF(INTER.EQ.2) THEN
            KPAT = 2

            IF(draprn().GT.0.5) THEN
               KPA = -ISIGN(1,K(1,2))*IABS(KPAT)
               KPA2 = ISIGN(1,K(1,2))*(IABS(KPAT)-1)
            ELSE
               KPA = ISIGN(1,K(1,2))*(IABS(KPAT)-1)
               KPA2 = -ISIGN(1,K(1,2))*IABS(KPAT)
            ENDIF
         ENDIF
      ELSEIF(IPRO.EQ.14) THEN
         IF(IWEI.EQ.0) THEN
            KPA = 4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         KPFL = IABS(KPA)
         KPA2=-KPA
         IF(INTER.EQ.2) THEN
            KPA = -ISIGN(1,K(1,2))*KPFL
            KPA2 = ISIGN(1,K(1,2))*(KPFL-1)
            IF(draprn().LT.0.5) THEN
               KPA = ISIGN(1,K(1,2))*(KPFL-1)
               KPA2 = -ISIGN(1,K(1,2))*KPFL
            ENDIF
         ENDIF
         AM(1)=DBLE(ULMASS(KPA))
         AM(2)=DBLE(ULMASS(KPA2))
      ELSEIF(IPRO.EQ.12.OR.IPRO.EQ.15) THEN
         NFLAVP = NFLAV
         KPA = 1
      ENDIF
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
c          write(6,*) '0th',kpa,kpa2
      ELSEIF(IGENFL.EQ.1) THEN
         KPA = KPAO
         KPA2 = KPAO2
c         write(6,*) '1st',kpa,kpa2
         IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            AM(1)=DBLE(ULMASS(KPA))
            AM(2)=DBLE(ULMASS(KPA2))
         ENDIF
      ENDIF
c         IF(IDEBUG.EQ.1) write(6,*) 'partdf KPA = ',KPA
      IF(IPRO.EQ.100) THEN
         AM(1)=1.D0
         AM(2)=2.D0*DBLE(ULMASS(211))
c         write(6,*) ' am ',am(1),am(2)
      ENDIF
c      write(6,*) ' ULMASS ',ULMASS(113),ULMASS(211)
C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... XP2 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,XP2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NDIMEN = NDIM
      IF(IDISDIF.GE.1) NDIMEN=NDIM+2
      NRN = 0
      XMAX=1.d0-XF
      IF(IRES(1).EQ.0) THEN
c.......
         IF(KE.NE.22) THEN
            NRN = NRN + 1
            YX = YMIN*((YMAX/YMIN)**X(NRN))
            IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 260
            Q2MIN=ME*ME*YX*YX/(1.D0-YX)
            IF(QMI.GT.Q2MIN) Q2MIN = QMI
            Q2MAX=YX*SSS - W12
            IF(QMA.LT.Q2MAX) Q2MAX = QMA
            IF(Q2MAX.LT.Q2MIN) GOTO 270
            NRN = NRN + 1
            IF(QMI.EQ.0.D0) THEN
               Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
            ELSE
c new try weighting with 1/q**4
               Q2 = Q2MIN*Q2MAX/(X(NRN)*Q2MIN + Q2MAX*(1.D0 - X(NRN)))
            ENDIF
c            write(6,*) 'PARTDF: Q2MIN,Q2MAX,Q2 ',Q2MIN,Q2MAX,Q2
            IF(INRNCH.EQ.1) Then
               write(6,*) ' after y,q2: nrn = ',nrn
            endif
            XMIN = 0.0D0
            XP2Q=Q2/YX/SSS
            IF(IPRO.NE.12) THEN
               SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
               IF(IPRO.EQ.100) SMIN = 2.D0*(AM(1)**2+AM(2)**2)
               XMIN=(SMIN+Q2)/(YX*SSS)
               IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
                  XMIN = DMAX1(XP2Q,XMIN)
cnew
                  XMAX = XR12
cnew
               ENDIF
            ELSE
               XMIN=Q2/YX/SSS
            ENDIF
         ELSEIF(KE.EQ.22) THEN
            Q2 = 1.D0
            YX = 1.D0
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            IF(IPRO.EQ.100) SMIN = 2.D0*(AM(1)**2+AM(2)**2)
            XMIN=SMIN/(YX*SSS)
         ELSE
            write(6,*) ' wrong KF selected; program stop '
            YX = 0.D0
         ENDIF
         XG1=YX
      ELSEIF(IRES(1).EQ.1) THEN
C XG1 = E_gluon/E_electron
C YX  = E_photon/E_electron
C XP1 = E_gluon/E_photon
         NRN = NRN + 1
         YX = YMIN*((YMAX/YMIN)**X(NRN))
         NRN = NRN + 1
         XP1MIN =DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
     +           /YX/SSS
         XP1 = XP1MIN*(1.D0/XP1MIN)**X(NRN)
         XG1 = XP1*YX
         Q2MIN = ME*ME*YX*YX/(1.D0 - YX)
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         Q2MAX=YX*(SSS-MP**2) - W12
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
         IF(Q2MAX.LT.Q2MIN) GOTO 270
         NRN = NRN + 1
         IF(QMI.EQ.0.D0) THEN
            Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         ELSE
c new try weighting with 1/q**4
            Q2 = Q2MIN*Q2MAX/(X(NRN)*Q2MIN + Q2MAX*(1.D0 - X(NRN)))
         ENDIF
c         write(6,*) ' partdf yx =',yx,' xel = ',xg1
         XMIN = 0.0D0
         XP2Q=Q2/XG1/SSS
         QG2 = 0.D0
         IF(IPRO.NE.12) THEN
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            XMIN=(SMIN+Q2+QG2)/(XG1*SSS)
c         write(6,*) yx,xg1,xp1,smin,q2,xmin
            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XMIN = DMAX1(XMIN,Q2/YX/SSS)
c                    XMAX = XR12
            ENDIF
         ELSE
            XMIN=Q2/YX/SSS
         ENDIF
      ELSE
         YX = 0.D0
         XG1 = 0.D0
      ENDIF
C ... XP2 = E_glue /E_proton
C ... XX  = E_glue/E_pomeron = xg_pom
C ... XR  = E_pomeron/E_proton = x_pom
      IF(IPRO.NE.12.AND.IMIX.EQ.0) THEN
         NRN = NRN + 1
         XP2= XMIN*((XMAX/XMIN)**X(NRN))
c       write(6,*) ' partdf xp2,xmin,xmax',xp2,xmin,xmax
      ELSEIF(IPRO.NE.12.AND.IMIX.EQ.1.AND.IWEI.EQ.1) THEN
         XP2= XMIN*((XMAX/XMIN)**X(NDIMEN+3))
         IF(ICHECK.EQ.1)  THEN
            write(6,*) 'partdf x(ndimen+3),ndimen,xp2,xmax,xmin '
     +      ,x(ndimen+3),ndimen,xp2,xmax,xmin
         ENDIF
      ELSEIF(IPRO.EQ.12.OR.IPRO.EQ.100) THEN
         XP2 = Q2/YX/SSS
      ELSE
         XP2 = 0.D0
      ENDIF
      IF(INRNCH.EQ.1) Then
         write(6,*) ' after xp2: nrn = ',nrn
      endif

c check for w range
c      ww = dsqrt(-q2 + yx*sss)
c      if(ww.lt.50.or.ww.gt.220) goto 456



      IF(SNGL(XP2).GT.SNGL(XMAX).OR.SNGL(XP2).LT.SNGL(XMIN)) GOTO 290
      SHAT1=SSS*XG1*XP2
      XMCHECK = DMAX1((AM(1)+AM(2))**2,DBLE(ULMASS(113)))
      IF(IPRO.EQ.12) XMCHECK = AM(1)
      IF(SHAT1.LT.XMCHECK) GOTO 300
      ENTRY PARTDFHS(X,WPART)
      IF(IHERAC.EQ.1) THEN
c test
         MSTJ24 = MSTJ(24)
         MSTJ(24) = 0
         DO 40    IN=1,20
            K(IN,1) = 0
            K(IN,2) = 0
   40    CONTINUE

C...  GIVE BEAM  FOUR VECTORS
         DO 50   I=1,2
            DO 50   J=1,5
               K(I,J)=KBEAM(I,J)
   50    P(I,J) = PBEAM(I,J)
         ME =P(1,5)
         MP =P(2,5)
         NP=2
         NPP=3*NP
         N=2
         NFLAVP = NFLAV
         NRN = 2
         MP =P(2,5)
         XMAX= 1.d0 -XF
         IST = 0
         YX = YHS
         XP2 = XHS
         Q2 = Q2HS
c        write(6,*) '1st partdf IPRO,IHFLA,KPA',IPRO,IHFLA,KPA,IGENFL
         IF(IPRO.EQ.12.OR.IPRO.EQ.15) THEN
            KPA = KPAHS
         ENDIF
         XG1 = YX
         NDIMEN = NDIM
C ... select particle code for light flavour production according
C ... to charges
         IF(IPRO.EQ.10.OR.IPRO.EQ.13) THEN
            AM(1) = 0.0D0
            AM(2) = 0.0D0
            IF(INTER.LT.2) THEN
               QF=DFLOAT(LUCHGE(1))**2 + DFLOAT(LUCHGE(2))**2 +
     +         DFLOAT( LUCHGE(3))**2
               QF = 2.D0*QF
               KPA = -4
               QFT = - draprn()*QF
   60          KPA=KPA+1
               QFT = QFT + DBLE(LUCHGE(KPA))**2
               IF(QFT.LT.0.0D0) GOTO 60
               IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!!'
     +         //' ', KPA
               KPA2 = - KPA
            ELSEIF(INTER.EQ.2) THEN
               KPAT = 2
               IF(draprn().GT.0.5) THEN
                  KPA = -ISIGN(1,K(1,2))*IABS(KPAT)
                  KPA2 = ISIGN(1,K(1,2))*(IABS(KPAT)-1)

               ELSE
                  KPA = ISIGN(1,K(1,2))*(IABS(KPAT)-1)
                  KPA2 = -ISIGN(1,K(1,2))*IABS(KPAT)
               ENDIF
            ENDIF
         ENDIF
         IF(IPRO.EQ.14) THEN
            IF(IWEI.EQ.0) THEN
               KPA=4
               IF(IHFLA.GE.4) KPA=IHFLA
               KPA2=-KPA
            ELSEIF(IWEI.EQ.1) THEN
               IF(draprn().GT.0.5) KPA = -KPA
               KPA2 = -KPA
            ENDIF
            KPFL = IABS(KPA)
            IF(INTER.EQ.2) THEN
               KPA = -ISIGN(1,K(1,2))*KPFL
               KPA2 = ISIGN(1,K(1,2))*(KPFL-1)
               IF(draprn().LT.0.5) THEN
                  KPA = ISIGN(1,K(1,2))*(KPFL-1)
                  KPA2 = -ISIGN(1,K(1,2))*KPFL
               ENDIF
            ENDIF
            AM(1)=DBLE(ULMASS(KPA))
            AM(2)=DBLE(ULMASS(KPA2))
            IF(AM(1).LT.1.D-6) AM(1)=0.0D0
            IF(AM(2).LT.1.D-6) AM(2)=0.0D0
         ENDIF
c         IF(IDEBUG.EQ.1) write(6,*) 'partdf KPA = ',KPA
c         IF(IPRO.EQ.14) write(6,*) ' 2 partdf ',KPA,AM(1)
c         write(6,*) 'parths ',Q2
c         write(6,*) 'PARTDF SSS ',SSS
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
         IF(IPRO.EQ.14.AND.IMIX.EQ.1) THEN
         ELSE
            W02=(1.D0 + MP)**2
            W12=W02-MP*MP
            YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
            YMAX=YMAX/(2.D0*(SSS+ME*ME))
            YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
            YMIN=YMIN/(2.D0*(SSS+ME*ME))
c         write(6,*) ' YX,YHS,YMIN,YMAX ',YX,YHS,ymin,ymax
            IF(YMI.GT.YMIN) YMIN=YMI
            IF(YMA.LT.YMAX) YMAX=YMA
c         IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 260
         ENDIF
         IF(IGENFL.EQ.0) THEN
            KPAO = KPA
            KPAO2 = KPA2
         ELSE
            KPA=KPAO
            KPA2 = KPAO2
         ENDIF

         Q2MIN=ME*ME*YX*YX/(1.D0-YX)
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         Q2MAX=YX*SSS - W12
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
         IF(Q2MAX.LT.Q2MIN) GOTO 270
         XP2Q=Q2/YX/SSS
         IF(IPRO.NE.12) THEN
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            IF(IPRO.EQ.100) SMIN = 2.D0*(AM(1)**2+AM(2)**2)
            XMIN=(SMIN+Q2)/(YX*SSS)
            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XMIN = DMAX1(XP2Q,XMIN)
               XMAX = XR12
               NDIMEN = NDIM + 3
               XP2= XMIN*((XMAX/XMIN)**X(NDIMEN+3))
c               write(6,*) 'partdf x',x
c               write(6,*) 'partdf: xp2 ',x(ndimen+3),ndimen
            ENDIF
         ENDIF
c         write(6,*) ' PARTDHS: KPA,AM',KPA,AM(1),AM(2)
c         write(6,*) ' here in heracl partdf ',IHERAC
         IF(XP2.GT.XMAX) GOTO 280
c        write(6,*) '2nd partdf IPRO,IHFLA,KPA',IPRO,IHFLA,KPA,IGENFL

      ENDIF

      IST = 0
      YY=SNGL(YX)
      XEL=SNGL(XG1)
      XPR=SNGL(XP2)
      PHI = 99999.D0
      IF(IPRO.EQ.15) NFLAVP = NFLQCDC
      IF((IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.18)
     +  .AND.IMIX.EQ.0) THEN
         NRN=NRN+1
         PHI = 2.D0*PI*X(NRN)
         IF(ICHECK.EQ.1)  THEN
            write(6,*) 'partdf  phi random nrn,x(nrn),phi',nrn,x(nrn),
     +      phi
         ENDIF

      ELSEIF((IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.18)
     +  .AND.IMIX.EQ.1.
     +        AND.IWEI.EQ.1) THEN
         PHI = 2.D0*PI*X(NDIMEN+4)
         IF(ICHECK.EQ.1)  THEN
            write(6,*) 'partdf  phi imix/iwei ndimen,x(ndimen+4),phi',
     +      ndimen,x(ndimen+4),phi
         ENDIF
c         write(6,*) 'partdf phi',x(ndimen+4),ndimen
      ELSEIF(IHERAC.EQ.0) THEN
         IF(IGENFL.EQ.0) THEN
            PHI = 2.D0*PI*draprn()
            PIO = PHI
         ELSEIF(IGENFL.EQ.1) THEN
            PHI = PIO
         ENDIF
         IF(ICHECK.EQ.1)  THEN
            write(6,*) 'partdf  phi herac phi,pio',phi,pio
         ENDIF

      ENDIF
      IF(INRNCH.EQ.1) Then
         write(6,*) ' after phi: nrn = ',nrn
      endif

      IF(IRES(1).EQ.0) THEN
         CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      ELSEIF(IRES(1).EQ.1) THEN
         CALL PARTI(KE,YX,FGAM,WEIGHT,1,IST)
c factor xg1/yx cancelled because parton densities of photon taken
c xgx instead of g(x)
         FWEI =  DLOG(1.D0/XP1MIN) * WEIGHT
c       write(6,*) ' partdf fwei ',fwei
      ELSE
         WRITE(6,*) ' IRES(1) > 1 not implemented '
         FGAM = 0.D0
         FWEI = 0.D0
      ENDIF

c         write(6,*) 'parths Q2_gen,dot ',Q2,dot1(NIA1,NIA1)
c      write(6,*) ' partdf: iherpys',iherpys
      XP=0.D0
      XX=0.D0
      XR=0.D0
      IRES(2) = 1
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
         IF(THETE .GE. THEMA) GOTO 310
         IF(THETE .LE. THEMI) GOTO 310
      ENDIF
C end of these gymnastics

      IF(IPRO.EQ.10.OR.IPRO.EQ.11.OR.IPRO.EQ.12.OR.IPRO.EQ.13.
     +   OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.18.OR.IPRO.EQ.100) THEN
         IF(IPRO.NE.100) THEN
            NRN = NRN + 1
            IF(IVM.EQ.0) THEN
               XR = XP2*(XMAX/XP2)**X(NRN)
            ELSE
               IF(IVM.GE.1.AND.IVM.LT.443) THEN
                  VMAX = DBLE(1.020 + 2.*ULMASS(211))
                  VMIN = DBLE(0.780 - 2.*ULMASS(211))
               ELSEIF(IVM.EQ.443) THEN
                  VMAX = 2.D0*DBLE(ULMASS(421))
                  VMIN = 2.D0*DBLE(ULMASS(4))
               ELSEIF(IVM.EQ.553) THEN
                  VMAX = 2.D0*DBLE(ULMASS(521))
                  VMIN = 2.D0*DBLE(ULMASS(5))
               ENDIF
               VMIN = VMIN+SMALL
               VMAX = VMAX-SMALL
c                  write(6,*) 'partdf 1:',XMAXV,XMINV
               XMAXV = XP2*(1 + VMAX**2/Q2)
               XMINV = XP2*(1 + VMIN**2/Q2)
c                  write(6,*) 'partdf 2:',XMAXV,XMINV
               XR = XMINV*(XMAXV/XMINV)**X(NRN)
               IF(XMAXV.GT.XMAX) GOTO 290
            ENDIF
            IF(IHF.EQ.1.and.ipro.eq.12) THEN
               XMINHF = XP2*(1.d0 + 4.d0*DBLE(ULMASS(4))**2/Q2)
               IF(XMINHF.GT.XMAX) GOTO 290
               XR = XMINHF*(XMAX/XMINHF)**X(NRN)
            ENDIF
c          write(6,*) ' PARTDF: XR,X(NRN),NRN ',XR,X(NRN),NRN

            IF(IMIX.EQ.1.AND.IPRO.EQ.12.AND.IWEI.EQ.1) THEN
               XR12=XR
            ELSEIF(IMIX.EQ.1.AND.IPRO.NE.12.AND.IWEI.EQ.1) THEN
               XR=XR12
            ENDIF
c      IF(IWEI.EQ.1) THEN
c           write(6,*) ' PARTDF: XR,X(NRN),NRN ',XR,X(NRN),NRN
c           ENDIF
cnew
c            IF(XP2.GT.XR) THEN
c          write(6,*) ' partdf xp2 > xr ',xp2,xr
c          goto 345
c          endif
cnew
            IF((XP2/XR).LT.1.D-20) GOTO 320
            XX = XP2/XR
c            IF(XX.GT.1.0.OR.XX.LT.0.0)
c     +      write(6,*) 'PARTDF xx,xpr,xr',xx,xpr,xr
            IF(ICHECK.EQ.1) THEN
               write(6,*) 'partdf XX,xr,xp2,x(NRN) ',xx,xr,xp2,x(NRN)
            ENDIF
c         IF(XX.GT.0.9999) THEN
c         write(6,*) 'XX>.9999 ',xx,xr,xp2,x(NRN)
c         ENDIF
         ELSEIF(IPRO.EQ.100) THEN
            XR=XP2
            NDFF=4
         ENDIF
         IF(INRNCH.EQ.1) Then
            write(6,*) ' after x_pom: nrn = ',nrn
         endif

         T2MIN=MP*MP*XR*XR/(1.D0-XR)
c         T2MN = -(Q2  - XR*YX*SSS  + 4.D0*DBLE(ULMASS(211))**2)
         T2MN = -(Q2  - XR*XG1*SSS  + 4.D0*DBLE(ULMASS(211))**2)
c         write(6,*) 'partdf : t2mn,t2min,t2max',t2mn,t2min,t2max
         IF(T2MN.LT.T2MAX) THEN
            T2MAX1 = T2MN
         ELSE
            T2MAX1 = T2MAX
         ENDIF

         NDFF = 6
      ELSE
         WRITE(6,*) ' no t_min, t_max calculated '
         T2MIN = 0.D0
         NDFF = 0
      ENDIF
      IF(T2MIN.GE.T2MAX1) THEN
CCC            WRITE(6,*) 'T2MIN=',T2MIN,'.GT.T2MAX1=',T2MAX1
         GOTO 330
      ENDIF
      NRN = NRN + 1
      T2 = T2MIN * ((T2MAX1/T2MIN)**X(NRN))
      IF(ICHECK.EQ.1) THEN
         write(6,*) ' PARTDF: T2,X(NRN),NRN ',T2,X(NRN),NRN
      ENDIF

c      IF(IWEI.EQ.1) THEN
c          write(6,*) ' PARTDF: T2,X(NRN),NRN ',T2,X(NRN),NRN
c          ENDIF
      XR1= T2/XR/SSS
      T2=-T2
      T2GKI = SNGL(T2)
      XFGKI = SNGL(XR)
      IF(INRNCH.EQ.1) Then
         write(6,*) ' after t: nrn = ',nrn
      endif

      CALL RAT2DI(KINT(2,2),XR,T2,WTDIST)
C weight for t distribution
      WTDIST = WTDIST * (-T2) *DLOG(T2MAX1/T2MIN)
      IF(IPRO.NE.100) THEN
         IF(IVM.EQ.0.AND.IHF.EQ.0) THEN
            WTDIST=WTDIST*XR*DLOG((1.D0 - XF)/XP2)
         ELSEIF(IVM.EQ.0.AND.IHF.EQ.1.and.ipro.eq.12) THEN
            WTDIST=WTDIST*XR*DLOG(XMAX/XMINHF)
         ELSEIF(IVM.EQ.0.AND.IHF.EQ.1.and.ipro.ne.12) THEN
            WTDIST=WTDIST*XR*DLOG((1.D0 - XF)/XP2)
         ELSEIF(IVM.GE.1) THEN
            WTDIST=WTDIST*XR*DLOG(XMAXV/XMINV)
         ENDIF
      ELSE
         WTDIST=WTDIST*XR*DLOG(XMAX/XMIN)
      ENDIF
      IF(IMIX.EQ.1.AND.IPRO.EQ.12.AND.IWEI.EQ.1) THEN
         WTD12=WTDIST
         IF(ICHECK.EQ.1)  THEN
            write(6,*) 'partdf wtdist 1st ',wtdist,wtd12
         ENDIF

      ELSEIF(IMIX.EQ.1.AND.IPRO.NE.12.AND.IWEI.EQ.1) THEN
         WTDIST=WTD12
         IF(ICHECK.EQ.1)  THEN
            write(6,*) 'partdf wtdist 2nd ',wtdist,wtd12
         ENDIF

      ENDIF
C end weight for t distribution
      IF(DABS(WTDIST).LE.1.D-45) GOTO 340
      COSTP=(1.D0 - XR - XR1*XR)/(1.D0 - XR + XR1*XR)
      IF(COSTP.GT.1.D0) COSTP=1.0D0
      IF(COSTP.LT.-1.D0) COSTP=-1.0D0
C...  COSTP IS SCATTERING ANGLE OF PROTON
      IF(IDISDIF.EQ.0.OR.IMIX.EQ.0) THEN
         NRN = NRN + 1
         PHIP=2.D0*PI*X(NRN)
	   IF(IHERAC.EQ.1) PHIP=2.D0*PI*draprn()
      ELSE
         PHIP=2.D0*PI*draprn()
      ENDIF
      IF(IGENFL.EQ.0) THEN
         PHIO = PHIP
      ELSEIF(IGENFL.EQ.1) THEN
         PHIP = PHIO
      ENDIF
      IF(INRNCH.EQ.1) Then
         write(6,*) ' after phip: nrn = ',nrn
      endif

c      PHIP = 0.D0
c      IF(IDEBUG.EQ.1) write(6,*) 'partdf: PHIP=',PHIP
      SPHP=DSIN(PHIP)
      CPHP=DCOS(PHIP)

C FINAL STATE PROTON
      NPFIN=NIA1+NDFF
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      IF(IABS(KINT(2,2)).EQ.211) K(NPFIN,2) = 2112
      K(NPFIN,3)=2
      P(NPFIN,5)=DBLE(ULMASS(K(NPFIN,2)))
C NEW gymnastics because of numerical precision problems.....

c for construction of  pomeron goto gamma p system
      DBCMS(1)=  P(NIA1,1) + P(2,1)
      DBCMS(2)=  P(NIA1,2) + P(2,2)
      DBCMS(3)=  P(NIA1,3) + P(2,3)
      DBCMS(4)=  P(NIA1,4) + P(2,4)
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2
     +              + (DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 360

      CALL DUDBRB(0,0,0.D0,0.D0,
     +  -DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +  -DBCMS(3)/DBCMS(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c new
      PEP = DSQRT(P(2,3)**2 + DBLE(ULMASS(KP)**2))
      P(2,4) = PEP
      PEZ = ABS(P(2,3))
      PEG = P(NIA1,4)
      PEGZ = DABS(P(NIA1,3))
cres      POMDGA = XR * YX * SSS/2.D0
      POMDGA = XR * XG1 * SSS/2.D0
      EN= PEZ*(PEGZ*PEZ+PEG*PEP)/(PEZ*PEG+PEP*PEGZ) - ((T2-P(2,5)**2-
     +P(NPFIN,5)**2)*PEGZ/2.D0 + POMDGA*PEZ)/(PEZ*PEG+PEP*PEGZ)
      P(NPFIN,4) = EN
      PZC =       - (PEP*(PEGZ*PEZ+PEG*PEP)/(PEZ*PEG+PEP*PEGZ) +
     +             ((T2-P(2,5)**2-P(NPFIN,5)**2)*PEG/2.D0
     +             - POMDGA*PEP)/(PEZ*PEG+PEP*PEGZ))
      P(NPFIN,3) = PZC
      PT = (-(T2-P(NPFIN,5)**2)*(PEGZ*PEZ+PEG*PEP)**2 + (T2-P(2,5)**2-
     +P(NPFIN,5)**2)**2*Q2/4.D0 -P(NPFIN,5)**2*POMDGA**2 +(T2+P(2,5)**
     +2-P(NPFIN,5)**2)*POMDGA *(PEGZ*PEZ+PEG*PEP) )/(PEZ*PEG+PEGZ*PEP)*
     +*2 - P(2,5)**2
      P(NPFIN,1) = DSQRT(DMAX1(0.D0,PT))*CPHP
      P(NPFIN,2) = DSQRT(DMAX1(0.D0,PT))*SPHP
      CALL DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +  DBCMS(3)/DBCMS(4))
c end new

C end of these gymnastics

C MOMENTA OF POMERON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KINT(2,2)
      K(NIA1+1,3)=2
      DO 70  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   70 CONTINUE
      P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
      T2GKI=-SNGL(P(NIA1+1,5)**2)
c      IF(IDEBUG.EQ.1) THEN
c         write(6,*) ' PARTDF: T2GKI,T2 ',T2GKI,T2
c         write(6,*) ' partdf: yx=',yx,dot1(nia1,2)/dot1(1,2)
c         write(6,*) ' partdf: xr=',xr,dot1(nia1,nia1+1)/dot1(nia1,2)
c      ENDIF
c for construction of parton in pomeron goto gamma pomeron system

      DBCMS(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBCMS(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBCMS(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBCMS(4)=  P(NIA1,4) + P(NIA1+1,4)
      spom = dot(dbcms,dbcms)
c      IF(IWEI.EQ.1) write(6,*) ' PARTDF: SPOM= ',spom
      IF(SPOM.LT.0.0D0) GOTO 380
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2
     +              + (DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 360

      CALL DUDBRB(NIA1,NIA1+4,0.D0,0.D0,
     +  -DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +  -DBCMS(3)/DBCMS(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(NIA1,NIA1+4,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(NIA1,NIA1+4,-STHETA,0.d0,0.d0,0.d0,0.d0)
c      write(6,*) ' now in parton pom system partdf NIA1 ',NIA1
c      CALL DULIST(1)
      IF(IPRO.EQ.100) THEN
         NIA2= NIA1+1
         NF1=NIA1+2
         NF2=NIA1+3
         K(NF1,1)=1
         K(NF1,2)=KPA
         K(NF1,3)=NIA1
         K(NF2,1)=1
         K(NF2,2)=KINT(2,2)
         K(NF2,2)=83
         K(NF2,3)=NIA1
         NJOIN=0
         NB2 = NIA1 + 1
      ELSE
         K(NIA1+2,1)=21
         K(NIA1+2,2)=KGL
         K(NIA1+2,3)=NIA1+1
         IF(IPRO.NE.12) THEN
c            IF(Q2.LE.2) THEN
            IF(IRES(1).EQ.1) THEN
               P(NIA1+2,3) = XX*P(NIA1+1,3)
               P(NIA1+2,4) = DABS(P(NIA1+2,3))
            ELSE
               IF(Q2.LE.2.D0) THEN
                  P(NIA1+2,4)= XX*(P(NIA1+1,4)*P(NIA1,4)+ ABS(P(NIA1+1,
     +            3)* P(NIA1,3)) + DBLE(T2GKI)) /(P(NIA1,4)+ABS(P(NIA1,
     +            3)))
                  P(NIA1+2,3)= -P(NIA1+2,4)
               ELSE
                  P(NIA1+2,4)= -XX*(P(NIA1+1,4)*P(NIA1,4)+ ABS(P(NIA1+
     +            1,3) * P(NIA1,3)) + DBLE(T2GKI))*(P(NIA1,4)-P(NIA1,3)
     +            )/Q2
                  P(NIA1+2,3)= -P(NIA1+2,4)
               ENDIF
            ENDIF
c            write(6,*) 'construct NIA1+2: P,Q2,XX ',P(NIA1+2,4),Q2,XX
c            write(6,*) 'NIA1  ',NIA1

         ELSE
            DO 80 I=-6,6
   80       XPQ(I) = 0.0

            Q2Q = Q2
            IF(XX.LT.0.999D0)
     +      CALL RASTFU(KINT(2,2),SNGL(XX),SNGL(Q2Q),XPQ)
            IF(XX.GE.1.D0) THEN
               write(6,*) ' partdf: XX',XX
               write(6,*) ' X(I) ',X
            ENDIF
c            if(xpr.ne.xx*xfgki) write(6,*) xpr,xfgki*xx,xfgki,xx
c            if(iwei.eq.1) call HF1(301,sngl(xx),1.)
c            if(iwei.eq.1) write(6,*) 'partdf: beta = ',xx
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
            WTGLU = 0.D0
c           write(6,*) ' partdf: NFLAV ',NFLAVP,XPQ
            nafl =0
            DO 90  I=1,NFLAVP
               IF(INTER.LT.2) THEN
                  IF(SPOM.LT.(4.*ULMASS(I)**2)) GOTO 90
               ELSEIF(INTER.EQ.2) THEN
                  IF(I.GE.(NFLAVP+1)) GOTO 90
                  IF(SPOM.LT.((ULMASS(I)+ULMASS(I+1))**2)) GOTO 90
               ENDIF
               IF(XPQ(I).LE.0.) GOTO 90
               NAFL = I
   90       CONTINUE
            IF(INTER.LT.2) THEN
               DO 100 I=-NAFL,NAFL
  100          WTGLU =WTGLU+DBLE(XPQ(I))*DFLOAT(LUCHGE(I))**2/9.D0
               IF(WTGLU.LT.1.D-20) GOTO 350
               IF(IGENFL.EQ.1) THEN
                  KPA = KPAO
                  WTGLU = DBLE(XPQ(KPA))*DFLOAT(LUCHGE(KPA))**2/9.D0
                  KPF = KPA
               ELSE
  110             tdraprn =  draprn()
                  QFT = - tdraprn*WTGLU
                  KPA=-NAFL-1
  120             KPA=KPA+1
                  QFT = QFT + DBLE(LUCHGE(KPA))**2/9.D0*DBLE(XPQ(KPA))
                  IF(QFT.LT.0.0D0) GOTO 120
                  KPF = KPA
               ENDIF
               if(wtglu.le.0) then
                  write(6,*) ' wtglu < 0 '
                  write(6,*) 'check 1 ',nafl,wtglu,kpa,xpq,xx,q2q,
     +            igenfl
                  CALL RASTFU(KINT(2,2),SNGL(XX),SNGL(Q2Q),XPQ)
                  write(6,*) 'check 2 ',nafl,wtglu,kpa,xpq
               endif
               if(iabs(kpa).gt.nafl) then
                  write(6,*) 'partdf kpa,nafl',kpa,nafl
                  write(6,*) ' partdf wtglu = ',wtglu,tdraprn
                  write(6,*) ' partdf xpq ',xpq
                  write(6,*) ' luchge(0) ',luchge(0),igenfl
               endif
               IF(IVM.GT.0) THEN
                  IF(IVM.GE.1.AND.IVM.LT.443) THEN
                     VMAX = DBLE(1.020 + 2.*ULMASS(211))
                     VMIN = DBLE(0.780 - 2.*ULMASS(211))
                  ELSEIF(IVM.EQ.443) THEN
                     VMAX = 2.D0*DBLE(ULMASS(421))
                     VMIN = 2.D0*DBLE(ULMASS(4))
                     KPA = 4
                     WTGLU = 2.D0*DBLE(XPQ(4))*DFLOAT(LUCHGE(4))**2/
     +               9.D0
                  ELSEIF(IVM.EQ.553) THEN
                     VMAX = 2.D0*DBLE(ULMASS(521))
                     VMIN = 2.D0*DBLE(ULMASS(5))
                     KPA = 5
                     WTGLU = 2.D0*DBLE(XPQ(5))*DFLOAT(LUCHGE(5))**2/
     +               9.D0
                  ENDIF
                  VMIN = VMIN+SMALL
                  VMAX = VMAX-SMALL
                  KPF = KPA
c                  write(6,*) 'partdf spom =',spom,vmin,vmax
                  IF(SNGL(SPOM).LT.VMIN**2.OR.SPOM.GT.VMAX**2) GOTO
     +            440
               ENDIF
            ELSEIF(INTER.EQ.2) THEN
c            write(6,*) ISIGN(1,K(1,2))
               DO 130 I=-6,6
                  F2CC(I) = 0.0
  130          XF3CC(I) = 0.0
               DO 140 I=1,NAFL-1,2
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     F2CC(I+1) = 2.*XPQ(I+1)
                     F2CC(-I) = 2.*XPQ(-I)
                     XF3CC(I+1) = 2.*XPQ(I+1)
                     XF3CC(-I) = -2.*XPQ(-I)
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(-I)+F2CC(I+1)) + (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(-I)+XF3CC(I+1))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     F2CC(I) = 2.*XPQ(I)
                     F2CC(-I-1) = 2.*XPQ(-I-1)
                     XF3CC(I) = 2.*XPQ(I)
                     XF3CC(-I-1) = -2.*XPQ(-I-1)
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(I)+F2CC(-I-1)) - (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(I)+XF3CC(-I-1))
                  ENDIF
  140          CONTINUE
               IF(IGENFL.EQ.1) THEN
                  KPA = KPAO
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) + (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) - (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ENDIF
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     KPF = (KPA - 1)
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     KPF = (KPA + 1)
                  ENDIF
               ELSE
                  QFT = - draprn()*WTGLU
                  KPA=-NFLAVP-1
  150             KPA=KPA+1
                  IF(KPA.EQ.0) GOTO 150
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     QFT = QFT + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) + (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     QFT = QFT + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) - (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ENDIF

                  IF(QFT.LT.0.0D0) GOTO 150
c                  write(6,*) KPA, QFT
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     KPF = (KPA - 1)
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     KPF = (KPA + 1)
                  ENDIF
               ENDIF
            ELSE
               WRITE(6,*) ' wrong interaction INTER = ',INTER
            ENDIF
            K(NIA1+2,2) = KPA
            SCAL2 = SNGL(Q2Q)
            XPD2 = XPQ(KPA)
            IF(Q2.LE.2) THEN
               P(NIA1+2,3)= -(DBLE(ULMASS(KPF))**2+Q2)
     +          /(P(NIA1,4)+P(NIA1,3))/2.D0
c               P(NIA1+2,3)= -(Q2)/(P(NIA1,4)+P(NIA1,3))/
c     +         2.
               P(NIA1+2,4)= ABS(P(NIA1+2,3))
            ELSE
               P(NIA1+2,3)= (DBLE(ULMASS(KPF))**2+Q2)
     +           *(P(NIA1,4)-P(NIA1,3))/Q2/2.D0
c               P(NIA1+2,3)= (Q2) *(P(NIA1,4)-P(NIA1,3))/
c     +         Q2/2.
               P(NIA1+2,4)= ABS(P(NIA1+2,3))
            ENDIF
         ENDIF
         P(NIA1+2,1)= 0.0D0
         P(NIA1+2,2)= 0.0D0
         P(NIA1+2,5)= DOT1(NIA1+2,NIA1+2)
         IF(P(NIA1+2,5).LT.0)
     + P(NIA1+2,5)=-DSQRT(DMAX1(0.D0,P(NIA1+2,5)))
         IF(P(NIA1+2,5).GE.0)
     + P(NIA1+2,5)= DSQRT(DMAX1(0.D0,P(NIA1+2,5)))
         NIA2 = NIA1+2

         IF(IPRO.NE.12) THEN
            NF1=NIA1+3
            NF2=NIA1+5
            K(NF1,1)=2
            K(NF1,2)=KPA
            K(NF1,3)=NIA1
            K(NF2,1)=1
            K(NF2,2)=KPA2
            K(NF2,3)=NIA1
            IF(IPRO.EQ.15) THEN
               K(NF1,2)=KPF
               K(NF2,2) = 21

            ENDIF
            IJOIN(1)=NF1
            IF(IPRO.NE.15) THEN
               IJOIN(2)=NF1+1
               IJOIN(3)=NF2
            ELSE
               IJOIN(2) = NF2
               IJOIN(3) = NF1+1
            ENDIF
            NJOIN=3
C POMERON Remnant
            K(NF1+1,1)=2
            K(NF1+1,2)=KGL
            IF(IPRO.EQ.15) THEN
               K(NF1+1,2) = -KPA
            ENDIF
            K(NF1+1,3)=NIA1+1
            P(NF1+1,1)=P(NIA1+1,1) - P(NIA1+2,1)
            P(NF1+1,2)=P(NIA1+1,2) - P(NIA1+2,2)
            P(NF1+1,3)=P(NIA1+1,3) - P(NIA1+2,3)
            P(NF1+1,4)=P(NIA1+1,4) - P(NIA1+2,4)
            P(NF1+1,5)=DBLE(ULMASS(K(NF1+1,2)))
            P(NF1+1,5)= DOT1(NF1+1,NF1+1)
            IF(P(NF1+1,5).LT.0)
     +      P(NF1+1,5)=-DSQRT(DMAX1(0.D0,P(NF1+1,5)))
            IF(P(NF1+1,5).GE.0)
     +      P(NF1+1,5)= DSQRT(DMAX1(0.D0,P(NF1+1,5)))
            P(NF1,5)=AM(1)
            P(NF2,5)=AM(2)
            NB2 = NIA2
C...   VECTOR OF GAMMA GLUON CMS SYSTEM or gamma pomeron system
         ELSEIF(IPRO.EQ.12) THEN

c            P(NIA1+2,5)=ULMASS(KPA)
c            P(NIA1+2,5) = 0.0
c            P(NIA1+2,4)= SQRT(P(NIA1+2,3)**2 + P(NIA1+2,5)**2)
            NF1=NIA1+3
            NF2=NF1+2
            K(NF1,1)=2
            K(NF1,2)=KPF
            K(NF1,3)=NIA1
            K(NF2,1)= 0
            IJOIN(1)=NF1
            IJOIN(2)=NF1+1
            NJOIN=2
            NB2 = NIA1 + 1
c            write(6,*) 'partdf nf1,nia1,nf2',nf1,nia1,nf2
         ELSE
            NB2 = 0
         ENDIF
      ENDIF

      CALL DUDBRB(NIA1,NIA1+4,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(NIA1,NIA1+4,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(NIA1,NIA1+4,0.D0,0.D0,DBCMS(1)/DBCMS(4),
     +  DBCMS(2)/DBCMS(4),DBCMS(3)/DBCMS(4))
      DBCMS(1)=  P(NIA1,1) + P(NB2,1)
      DBCMS(2)=  P(NIA1,2) + P(NB2,2)
      DBCMS(3)=  P(NIA1,3) + P(NB2,3)
      DBCMS(4)=  P(NIA1,4) + P(NB2,4)
      DO 160 I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
  160 CONTINUE
c      write(6,*) ' before shat calc: NB2',NB2
c      Call Dulist(1)
      SHAT=DOT(DBCMS,DBCMS)
      IF(SHAT.LE.0.0) THEN
         GOTO 390
      ENDIF

C NOW BOOST TO GAMMA GLUON or gamma pomeron CMS
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2 +
     +(DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 370
      CALL DUDBRB(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +-DBCMS(3)/DBCMS(4))

      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      call DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      call DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
      IF(IPRO.EQ.12) THEN
         IF((SHAT-SMALL).LT.4.*ULMASS(211)**2) GOTO 430
         K(NIA2,2)=KPA
         IF((SHAT-SMALL).LT.4.*ULMASS(KPA)**2) GOTO 430
         IF(INTER.EQ.2) THEN
            IF((SHAT-SMALL).LT.(ULMASS(KPA)+ULMASS(KPF))**2) GOTO 430
         ENDIF
         WTGLU=WTGLU*WTDIST

         IF(SHAT.LE.0) THEN
            DO 170 I=1,4
  170       P(NF1,I) = P(NIA1,I) + P(NIA1+1,I)
            P(NF1,5)=DSQRT(SHAT)
            N=NPFIN
            K(NF1,1)=1
            K(NF1,2)=100*IABS(KPA) + 10*IABS(KPA) + 3
            K(NF1,3)=NIA1
            K(NF1+1,1)=0
         ELSE
            CALL DU2ENT(NF1,KPF,-KPA,SNGL(DSQRT(SHAT)))
c            write(6,*) 'partdf shat,KPA',SHAT,KPA
c            CALL DULIST(1)
            N = NPFIN
            K(NF1,3) = NIA1+1
            K(NF1+1,3) = NIA1 + 1
            K(NF1+1,1) = 1
            K(NF1+1,2) = -KPA
         ENDIF
         CALL DUEDIT(12)
c         IF(IWEI.EQ.1) write(6,*) 'PARTDF: spom,shat',spom,shat
c         DO 190 I=1,4
c  190    P(NIA1+2,I) = P(NF1,I) - P(NIA1,I)
c        CALL DULIST(1)
         NF2 = NF1 + 1
         ALPH_EM = ALPH
         IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
         IF(INTER.LT.2) THEN
            WMAT = 2.D0 * ALPH_EM*ALPH_EM*PI/YX/Q2**2
            WMAT = WMAT *(1.D0+(1.D0-YX)**2)
         ELSEIF(INTER.EQ.2) THEN
            GF = PI*ALPH_EM/(SIN2W*XMW2*DSQRT(2.D0))
            WMAT = GF**2/2.D0/PI /(1.D0 + Q2/XMW2)**2 /YX
capply additional factor 0.5 because only right/or left handed electrons contr.
            WMAT = WMAT * 0.5D0
         ENDIF
         WPART = WMAT * WTGLU *FWEI

C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
         call DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
         call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
         CALL DUDBRB(0,N,0.D0,0.D0,DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +   DBCMS(3)/DBCMS(4))

         PT2H = 99999.
         SHH = 99999.
C here do gymnastics for IPRO = 100 gamma pomeron --> rho pomeron
      ELSEIF(IPRO.EQ.100) THEN
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
         ECM =DSQRT(SHAT)
c         write(6,*) ' 100: am(1),am(2) ',am(1),am(2)
c       IF(ECM.LE.(AM(1)+AM(2))) WRITE(6,*) 'ECM LE MASS',ECM,dsqrt(spom)
         IF(ECM.LE.(AM(1)+AM(2))) GOTO 400
c    we set scale to Q2Q = 99999. because it's irrelevant here
         Q2Q = 99999.D0
c here we set WT=1 because in ELERHO the momenta are seperately generated
         WT = 1.D-4
C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
         call DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
         call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
         WMAT=WT *WTDIST*FGAM*FWEI
         WPART = WMAT
      ELSE
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
         ECM =DSQRT(SHAT)
         IF(ECM.LE.(AM(1)+AM(2))) GOTO 400
         NRN = NRN + 1
         XY(1) = X(NRN)
         NRN = NRN + 1
         XY(2) = X(NRN)
         IF(INRNCH.EQ.1) Then
            write(6,*) ' after phase: nrn = ',nrn
         endif

c          write(6,*) 'partdf:',x(nrn-1),x(nrn),nrn,ndim,ndimen,ipro
         CALL PHASE(NP,ECM,AM,PCM,WT)
         IF(WT.LE.0.D0) GOTO 410
c         write(6,*) ' PARTDF phase 1sr :',WT,XY(2)*DLOG(1.D0/SIMIN)
c         write(6,*) ' PARTDF phase AM:',am(1),am(2)
         DO 180 I=1,4
            P(NF1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,2)
  180    CONTINUE
         PT2 = DPLU(NF1,9)
         CALL CUTG(PT2,NACC)
c         IF(IDEBUG.EQ.1) write(6,*) 'PARTDF: PT2 = ',PT2
         IF(NACC.EQ.0) GOTO 420

cscale         IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
cscale            Q2Q = (AM(1)+AM(2))**2
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
         ELSE
            WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
            STOP
         ENDIF
cscale         ENDIF
         IF(IQ2.EQ.5) THEN
            Q2Q=Q2 + PT2*SCALFA + (2.D0*AM(1))**2
         ELSE
            Q2Q=SCALFA*Q2Q
         ENDIF
c do qcd compton stuff here
         IF(IPRO.EQ.15) THEN
c this should not be here               Q2Q = Q2
            IF(XX.LT.0.999D0) CALL RASTFU(KINT(2,2),SNGL(XX),SNGL(Q2Q),
     +      XPQ)
            IF(XX.GE.1.D0) THEN
c                  write(6,*) ' partdf: XX',XX
c                  write(6,*) ' X(I) ',X
               DO 190 I=-6,6
  190          XPQ(I) = 0.0
            ENDIF
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
            DO 200 I=1,NFLAVP
               IF(SPOM.LT.(4.*ULMASS(I)**2)) GOTO 200
               NAFL = I
  200       CONTINUE
            WTGLU = 0.D0
c           write(6,*) ' partdf: NFLAV ',NFLAVP,XPQ
            IF(INTER.LT.2) THEN
               DO 210 I=-NAFL,NAFL
  210          WTGLU =WTGLU+DBLE(XPQ(I))*DFLOAT(LUCHGE(I))**2/9.D0
               IF(WTGLU.LT.1.D-20) GOTO 350
               IF(IGENFL.EQ.1) THEN
                  KPA = KPAO
                  WTGLU = DBLE(XPQ(KPA))*DFLOAT(LUCHGE(KPA))**2/
     +            9.D0
               ELSE
                  QFT = - draprn()*WTGLU
                  KPA=-NAFL-1
  220             KPA=KPA+1
                  QFT = QFT + DBLE(LUCHGE(KPA))**2/9.D0* DBLE(XPQ(KPA))
                  IF(QFT.LT.0.0D0) GOTO 220
               ENDIF
               KPF = KPA
            ELSEIF(INTER.EQ.2) THEN
c            write(6,*) ISIGN(1,K(1,2))
               DO 230 I=-6,6
                  F2CC(I) = 0.0
  230          XF3CC(I) = 0.0
               DO 240 I=1,NAFL-1,2
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     F2CC(I+1) = 2.*XPQ(I+1)
                     F2CC(-I) = 2.*XPQ(-I)
                     XF3CC(I+1) = 2.*XPQ(I+1)
                     XF3CC(-I) = - 2.*XPQ(-I)
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)* DBLE(F2C
     +               C(-I)+F2CC(I+1)) + (YX - YX**2/2.D0)* DBLE(XF3CC(-
     +               I)+XF3CC(I+1))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     F2CC(I) = 2.*XPQ(I)
                     F2CC(-I-1) = 2.*XPQ(-I-1)
                     XF3CC(I) = 2.*XPQ(I)
                     XF3CC(-I-1) = - 2.*XPQ(-I-1)
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)* DBLE(F2C
     +               C(I)+F2CC(-I-1)) - (YX - YX**2/2.D0)* DBLE(XF3CC(I
     +               )+XF3CC(-I-1))
                  ENDIF
  240          CONTINUE
               IF(IGENFL.EQ.1) THEN
                  KPA = KPAO
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)* DBLE(F2C
     +               C(KPA)) + (YX - YX**2/2.D0)* DBLE(XF3CC(KPA))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)* DBLE(F2C
     +               C(KPA)) - (YX - YX**2/2.D0)* DBLE(XF3CC(KPA))
                  ENDIF
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     KPF = (KPA - 1)
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     KPF = (KPA + 1)
                  ENDIF
               ELSE
                  QFT = - draprn()*WTGLU
                  KPA=-NFLAVP-1
  250             KPA=KPA+1
                  IF(KPA.EQ.0) GOTO 250
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     QFT = QFT + (1.D0 - YX + YX**2/2.D0)* DBLE(F2CC(KP
     +               A)) + (YX - YX**2/2.D0)* DBLE(XF3CC(KPA))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     QFT = QFT + (1.D0 - YX + YX**2/2.D0)* DBLE(F2CC(KP
     +               A)) - (YX - YX**2/2.D0)* DBLE(XF3CC(KPA))
                  ENDIF

                  IF(QFT.LT.0.0D0) GOTO 250
c               write(6,*) KPA, QFT
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     KPF = (KPA - 1)
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     KPF = (KPA + 1)
                  ENDIF
               ENDIF
            ELSE
               WRITE(6,*) ' wrong interaction INTER = ',INTER
            ENDIF
            WTG15 = WTGLU
            K(NIA1+2,2) = KPA
            AM(1) = 0.D0
            IF(IABS(KPA).GE.4) AM(1) = DBLE(ULMASS(KPA))
            SCAL2 = SNGL(Q2Q)
            XPD2 = XPQ(KPA)
            K(NF1,2)=KPF
            K(NF2,2) = 21
            IJOIN(1)=NF1
            IJOIN(2) = NF2
            IJOIN(3) = NF1+1
            NJOIN=3
            K(NF1+1,2) = -KPA

         ENDIF
c end qcd compton

C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
         call DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
         call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
         XPR = SNGL(XP2)
         IF(XX.LT.0.999D0)
     +   CALL RASTFU(KINT(2,2),SNGL(XX),SNGL(Q2Q),XPQ)
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
ccccc         write(6,*) ' partdf XPQ ',IPRO,xpq

c for resolved photon move parton densities into ELERES
         IF(IRES(1).EQ.1) THEN
            XP = 1.D0/XX
         ELSE
c            write(6,*) ' partdf ires=0'
c            if(xx.gt.0.999999) write(6,*) ' partdf xx>1. ',xx,xpq(0)
            XP = DBLE(XPQ(0))/XX
         ENDIF
         SCAL2 = SNGL(Q2Q)
         XPD2 = XPQ(0)
         IF(IPRO.EQ.15) THEN
            XP = WTG15/XX
            XPD2 = XPQ(KPA)
         ENDIF
         XPINT=DLOG(XMAX/XMIN)
cnew
c         IF(IMIX.EQ.1.AND.IPRO.NE.12.AND.IWEI.EQ.1) THEN
c         XPINT=DLOG(XR/XMIN)
c        endif
cnew
         WTGLU = XX * XP * XPINT
         WTGLU=WTGLU*WTDIST
         WMAT=WT *WTGLU*FGAM*FWEI
         WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
         FLUX=(2.D0*PI)**4
         FLUX=FLUX/(2.D0*(SHAT+Q2))
         WPART = WMAT *FLUX
         IF(IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15) THEN
            WPART = WTGLU * FWEI * WT
            IF(IPRO.EQ.15) THEN
               WPART = WPART*FLUX/((2.D0*PI)**NPP)
            ENDIF
         ENDIF
         PT2H = SNGL(PT2)
         SHH = SNGL(SHAT)

      ENDIF
      IF(IFPS.NE.10.AND.NPOM.NE.20.AND.NPOM.NE.21.AND.IPRO.NE.100) THEN
c         CALL DULIST(1)
c         IF(IWEI.EQ.1) THEN
c            write(6,*) ' partdf ', kpa,kpa2
c            call dulist(1)
c            endif
         CALL LUJOIN(NJOIN,IJOIN)
c         CALL DULIST(1)
      ENDIF
      NDIMC = NRN
      IF(INRNCH.EQ.1) Then
         write(6,*) ' end: nrn = ',nrn
      endif

      IF(WPART.LT.0.D0) THEN
         WRITE(6,*) ' WPART < 0 '
         WRITE(6,*) ' WMAT = ',WMAT,' WTGLU = ',WTGLU,' FWEI = ',FWEI
         WRITE(6,*) ' WT = ',WT, 'WTDIST = ',WTDIST
         WRITE(6,*) ' NAFL =',nafl,' x = ',xpr,' xpom =',xr,' xx=',xx
      ENDIF
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
c         write(6,*) '2nd ', kpa, kpa2
      ELSE
         KPA = KPAO
         KPA2 = KPAO2
c         write(6,*) '3rd ', kpa, kpa2
      ENDIF
c      write(6,*) 'partdf wpart,q2,xp2,yy,ipro',wpart,q2,xp2,q2q,yy,ipro
c         WRITE(6,*) 'partdf igenfl,imix,iwei ',igenfl,imix,iwei
      MSTJ(24) = MSTJ24
c      if(iwei.eq.1)      write(6,*) kpa,kpa2,kpf
c      call lulist(1)
      RETURN
  260 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  270 CONTINUE
      WPART = 0.D0
      IF(IWEI.EQ.1)
     +   write(6,*) ' partdf: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  280 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdf: xlimit ; RETURN xmin,xmax',xmin,xmax,IPRO
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: xlimit ; RETURN xmin,xmax',xmin,xmax,
     +   IPRO
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  290 CONTINUE
      WPART = 0.D0
c      write(6,*) ' partdf: xp2 limit ; RETURN ',xp2,xmax,xmin,IPRO
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: xp2 limit ; RETURN ',xp2,xmax,xmin,IPRO
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  300 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdf: slimit ; RETURN '
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: slimit ; RETURN '
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  310 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdf: theta limit ; RETURN '
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: theta limit ; RETURN '
      ENDIF
      IF(IHERAC.EQ.1) THEN
         NTHRE = NTHRE +1
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  320 CONTINUE
      WPART = 0.D0
c      write(6,*) ' partdf: xp2/xx limit ; RETURN ',IPRO
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: xp2/xx limit ; RETURN ',IPRO
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  330 CONTINUE
      WPART = 0.D0
c      write(6,*) ' partdf: tlimit ; RETURN '
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: tlimit ; RETURN ',T2,T2MIN,T2MAX1
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  340 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdf: wtdist limit ; RETURN ',WTDIST
      IF(IDEBUG.EQ.1) THEN
      write(6,*) ' partdf: wtdist limit; RETURN ',WTDIST,t2,T2MIN,T2MAX1
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  350 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdf: wtglu limit ; RETURN ',WTGLU
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: wtglu limit ; RETURN '
         write(6,*) ' partdf: WTGLU,xpq(4),nafl,nflavp,spom',
     +   WTGLU,xpq(4),nafl,nflavp,spom
         write(6,*) ' partdf: xr,q2,xbj,t2',xr,q2,xp2,t2
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  360 CONTINUE
      WPART = 0.D0

      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: bochck limit,pom ; RETURN ',BOCHCK,NDIM,
     +   NRN
         write(6,*) ' bochk:   P(NIA1,I)',(P(NIA1,I),I=1,4)
         write(6,*) ' bochk: P(NIA1+1,I)',(P(NIA1+1,I),I=1,4)
         write(6,*) ' bochk: DBCMS ',DBCMS
         write(6,*) ' bochk: X ',X
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  370 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: bochck limit, phase ; RETURN '
         write(6,*) ' bochk:   P(NIA1,I)',(P(NIA1,I),I=1,4)
         write(6,*) ' bochk: P(NIA1+1,I)',(P(NIA1+1,I),I=1,4)
         write(6,*) ' bochk: X ',X
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  380 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: spom=',spom,' test =',
     +   (-Q2+xr*yx*sss+t2)
c         write(6,*) ' partdf: nia1,nia1+1 ',nia1,nia1+1
         write(6,*) ' partdf: yx=',yx,dot1(nia1,2)/dot1(1,2)
         write(6,*) ' partdf: xr=',xr,dot1(nia1,nia1+1)/dot1(nia1,2)
         write(6,*) ' partdf: x_bj = ',xpr
c         write(6,*) 'partdf: 2gamma pom ',2*dot1(nia1,nia1+1),
c     +    xr*yx*sss
         write(6,*) ' partdf: Q2 ',Q2,dot1(nia1,nia1)
         write(6,*) ' partdf: T2 ',T2,T2min,T2max1
c         write(6,*) ' partdf:  test',(-Q2+dot1(nia1,nia1+1)/dot1(nia1,2)
c     +    *yx*sss+t2)
c        call lulist(1)
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  390 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: shat ; RETURN :x_bj',xpr,'x_pom',xr
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  400 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: ECM limit ; RETURN ',ECM,AM(1),AM(2)
         write(6,*) ' partdf: ecm_min ',dsqrt(-q2+xmin*yx*sss),xmin
      ENDIF
      MSTJ(24) = MSTJ24
      NDIMC = 9999
      RETURN
  410 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: PHASE wt=0 ; RETURN '
      ENDIF
      MSTJ(24) = MSTJ24
      NDIMC = 9999
      RETURN

  420 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: PTCUT limit ;RETURN ',PT2,PT2CUT(IPRO),
     +   IPRO
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  430 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) '350 partdf: shat;RETURN ',shat,ulmass(211),
     +   ulmass(kpa),nflav,nflavp
         write(6,*) '350 partdf: shat ; RETURN :x_bj',xpr,'x_pom',xr
         write(6,*) '350 partdf: shat,spom',shat,spom,kpa,ipro
      endif
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  440 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdf: VM_mass check ',IVM,SPOM,VMIN**2,VMAX**2
      ENDIF
      IF(IHERAC.EQ.1) THEN
         NVMRE = NVMRE +1
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  450 CONTINUE
      WPART = -1.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: xp2 > xr ',xp2,xr
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
      END
