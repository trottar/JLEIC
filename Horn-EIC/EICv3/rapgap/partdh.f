*CMZ :  2.08/04 22/12/99  15.39.28  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.56  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  17.25.59  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94


      SUBROUTINE partdh(X,WPART)
      IMPLICIT None
      Double Precision X,WPART,XY
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
*KEND.
      Double Precision PHI
      COMMON/DIFFA/ PHI
      Double Precision draprn
      DOUBLE PRECISION ME,MP
      COMMON/XVAR/ XY(10)
      DOUBLE PRECISION STHETA,SPHI,WMAX
      DIMENSION X(20)
      Integer IJOIN,NJOIN,NDIMC,IMIX
      DIMENSION IJOIN(10)
      COMMON /DIMEN/ NDIMC
      COMMON /OALPHAS/ WMAX,IMIX
*KEEP,RGPQCDPOM.
      Integer Iqqg
	COMMON/pqcdpom/ Iqqg
C     SAVE


*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEND.
      Integer NDIM,NPOIN,LST,IRES
      COMMON/DIVO/ NDIM,NPOIN
      COMMON/EPPARA/LST(30),IRES(2)
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
      Integer IVM,IGENFL,NTHRE,NVMRE,IDEBUG
      COMMON/VMESON/IVM
      COMMON/GENWEI/IGENFL
      COMMON/HERTHET/ NTHRE,NVMRE
*KEEP,RGRAHER.
      REAL XPQDIF,XPQPI
	Integer IHERPYS
	Integer NBQ2,NBX
      PARAMETER (NBQ2=20)
      PARAMETER (NBX=20)
      COMMON /RAHER/ IHERPYS,XPQDIF(-6:6,NBX,NBQ2),XPQPI(-6:6,NBX,NBQ2)
*KEND.
      COMMON/ INTERN/IDEBUG
      Double PRecision XPOM
      COMMON/BARTELS/XPOM
      DOUBLE PRECISION   sst,betat,xppt,t,qqt
      common   /parameter/ sst,qqt,betat,xppt,t
      Double Precision DHFORMF,W12,W02,QF,QFT,YX,SMIN
      Double Precision XG1,ww,XP2,FGAM,FWEI,XX,THETE
      Double Precision WTDIST,xbj,XR,djac,beta,T2MN,T2MAX1,T2MIN,T2
      Double Precision PEGZ
      Double Precision COSTP,PHIP,PHIO,SPHP,CPHP,BOCHCK,PEP,PEZ,PEG
      Double Precision POMDGA,EN,PZC,PT,spom,WTGLU,ECM,WT,PT2,XR1,DOT
      Integer I,J,In,NFLAVP,ICHECK,ITCH,IFIXB,MSTJ24,NP,NDFF
      Integer LUCHGE
      Integer KPA2,KPAO,KPAO2,NRN,IPSR,IST,NPFIN,KI,nafl
      Integer iph,NACC
      EXTERNAL draprn,DHFORMF
      DATA QF/0.0D0/,W12/0.0D0/,NFLAVP/0/
      DATA ICHECK/0/
c      IDEBUG = 0
      itch=0
      ifixb=0
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
      NP=3
      N=2
      QF = 0.D0
      NFLAVP = 0
      NDFF = 0
      AM(1) = 0.0D0
      AM(2) = 0.0D0
      NFLAVP = NFLAV
      IF(IHFLA.GE.4) THEN
         KPA = IHFLA
         AM(1) = DBLE(ULMASS(KPA))
         AM(2) = AM(1)
         NFLAVP = IHFLA
      ENDIF

C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
      W02=(AM(1)+AM(2)+MP)**2
      IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
c         WRITE(6,10000) YMIN,YMAX
c 10000 FORMAT(' limits on y ',/,' YMIN = ',E10.5,' YMAX = ',E10.5)

C ... select particle code for light flavour production according
C ... to charges
      IF(IHFLA.GE.4) THEN
         QF=DFLOAT(LUCHGE(4))**2
         KPA = IHFLA
      ELSE
         QF=DFLOAT(LUCHGE(1))**2 + DFLOAT(LUCHGE(2))**2 + DFLOAT(
     +   LUCHGE(3 ))**2
         QF = 2.D0*QF
         KPA = -4
         QFT = - draprn()*QF
   30    KPA=KPA+1
         QFT = QFT + DBLE(LUCHGE(KPA))**2
         IF(QFT.LT.0.0D0) GOTO 30
         IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!! ',
     +   KPA
      ENDIF
      KPA2 = - KPA
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
c          write(6,*) '0th',kpa,kpa2
      ELSEIF(IGENFL.EQ.1) THEN
         KPA = KPAO
         KPA2 = KPAO2
c         write(6,*) '1st',kpa,kpa2
      ENDIF
      AM(1) = DBLE(ULMASS(KPA))
      AM(2) = AM(1)
c         IF(IDEBUG.EQ.1) write(6,*) 'partdh KPA = ',KPA
c      write(6,*) ' ULMASS ',ULMASS(113),ULMASS(211)
C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... XP2 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,XP2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NRN = 0
      IPSR=0
      XMAX=1.d0-XF
      IF(IPSR.EQ.0) THEN
c.......
         IF(KE.NE.22) THEN
            NRN = NRN + 1
            YX = YMIN*((YMAX/YMIN)**X(NRN))
            IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 110
            Q2MIN=ME*ME*YX*YX/(1.D0-YX)
            IF(QMI.GT.Q2MIN) Q2MIN = QMI
            Q2MAX=YX*SSS - W12
            IF(QMA.LT.Q2MAX) Q2MAX = QMA
            IF(Q2MAX.LT.Q2MIN) GOTO 120
            NRN = NRN + 1
            IF(QMI.EQ.0.D0) THEN
               Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
            ELSE
c new try weighting with 1/q**4
               Q2 = Q2MIN*Q2MAX/(X(NRN)*Q2MIN + Q2MAX*(1.D0 - X(NRN)))
            ENDIF
c            write(6,*) 'partdh: Q2MIN,Q2MAX,Q2 ',Q2MIN,Q2MAX,Q2
            XMIN = 0.0D0
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),
     +        2.D0*(AM(1)**2+AM(2)**2+AM(3)**2))
            XMIN=(Q2+SMIN)/YX/SSS
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
      ELSE
         YX = 0.D0
         XG1 = 0.D0
      ENDIF
c check for w range
      ww = dsqrt(-q2 + yx*sss)
c      if(ww.lt.50.or.ww.gt.220) goto 300


C ... XP2 = E_glue /E_proton
C ... XX  = E_glue/E_pomeron
C ... XR  = E_pomeron/E_proton
      XP2 = Q2/YX/SSS
      IF(SNGL(XP2).GT.SNGL(XMAX)) GOTO 140
      IF(SNGL(XMIN).GT.SNGL(XMAX)) GOTO 130
c      SHAT1=SSS*XG1*XP2
c      XMCHECK = DMAX1((AM(1)+AM(2))**2,DBLE(ULMASS(113)))
c      IF(SHAT1.LT.XMCHECK) GOTO 150
      ENTRY partdhHS(X,WPART)
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
         N=2
         NFLAVP = NFLAV
         NRN = 2
         MP =P(2,5)
         IPSR=0
         XMAX= 1.d0 -XF
         IST = 0
         YX = YHS
         XP2 = XHS
         Q2 = Q2HS
         XG1 = YX
C ... select particle code for light flavour production according
C ... to charges
         QF=DFLOAT(LUCHGE(1))**2 + DFLOAT(LUCHGE(2))**2 + DFLOAT(
     +   LUCHGE(3))**2
         QF = 2.D0*QF
         KPA = -4
         QFT = - draprn()*QF
   60    KPA=KPA+1
         QFT = QFT + DBLE(LUCHGE(KPA))**2
         IF(QFT.LT.0.0D0) GOTO 60
         IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!!'
     +   //' ', KPA
         KPA2 = - KPA
c         IF(IDEBUG.EQ.1) write(6,*) 'partdh KPA = ',KPA
c         write(6,*) 'parths ',Q2
c         write(6,*) 'PARTDI SSS ',SSS
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
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
         IF(Q2MAX.LT.Q2MIN) GOTO 120
c         write(6,*) ' PARTDHS: KPA,AM',KPA,AM(1),AM(2)
c         write(6,*) ' here in heracl partdh ',IHERAC
         IF(XP2.GT.XMAX) GOTO 130

      ENDIF

      IST = 0
      YY=SNGL(YX)
      XEL=SNGL(XG1)
      XPR=SNGL(XP2)
      PHI = 99999.D0
      NRN=NRN+1
      PHI = 2.D0*PI*X(NRN)
c      PHI = 2.D0*PI*draprn()
      IF(ICHECK.EQ.1) THEN
         write(6,*) 'partdh  phi random nrn,x(nrn),phi',nrn,x(nrn),
     +   phi
      ENDIF
      IF(ICHECK.EQ.1) THEN
         write(6,*) 'partdh  phi herac phi',phi
      ENDIF

      IF(IPSR.EQ.0) THEN
         CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      ELSEIF(IPSR.EQ.1) THEN
         WRITE(6,*) ' IPSR = 1 not implemented '
         FGAM = 0.D0
         FWEI = 0.D0
      ELSE
         WRITE(6,*) ' IPSR > 1 not implemented '
         FGAM = 0.D0
         FWEI = 0.D0
      ENDIF
c         write(6,*) 'parths Q2_gen,dot ',Q2,dot1(NIA1,NIA1)
c      write(6,*) ' partdh: iherpys',iherpys
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
         IF(THETE .GE. THEMA) GOTO 160
         IF(THETE .LE. THEMI) GOTO 160
      ENDIF
      WTDIST = 1.D0
C end of these gymnastics
cc here for fixed beta
      xbj = q2/yx/sss
      IF(ifixb.eq.0) then
         NRN = NRN + 1
         XR = XMIN*(XMAX/XMIN)**X(NRN)
         WTDIST=WTDIST*XR*DLOG(XMAX/XMIN)
c jacobian dM**2 --> d x_pom
         djac = q2/xbj
         WTDIST = WTDIST*DJAC
      else
         beta = 0.1D0
         xr = xp2/beta
      endif
check beta
      beta = xp2/xr
      XPOM = XR
c     write(6,*) beta,xp2,xr
c      if(beta.gt.0.1) goto 300


c end of fixed beta
c      IF(IWEI.EQ.1) THEN
c           write(6,*) ' partdh: XR,X(NRN),NRN ',XR,X(NRN),NRN
c           ENDIF
      IF((XP2/XR).LT.1.D-20) GOTO 170
      XX = XP2/XR
      IF(XX.GT.1.0.OR.XX.LT.0.0) write(6,*) 'partdh xx,xpr,xr',xx,xpr,
     +xr
c      write(6,*) 'partdh xx,xpr,xr',xx,xpr,xr
      IF(ICHECK.EQ.1) THEN
         write(6,*) 'partdh XX,xr,xp2,x(NRN) ',xx,xr,xp2,x(NRN)
      ENDIF
c         IF(XX.GT.0.9999) THEN
c         write(6,*) 'XX>.9999 ',xx,xr,xp2,x(NRN)
c         ENDIF
      T2MIN=MP*MP*XR*XR/(1.D0-XR)
      T2MN = -(Q2 - XR*YX*SSS + 4.D0*DBLE(ULMASS(211))**2)
c         write(6,*) 'partdh : t2mn,t2min,t2max',t2mn,t2min,t2max
      IF(T2MN.LT.T2MAX) THEN
         T2MAX1 = T2MN
      ELSE
         T2MAX1 = T2MAX
      ENDIF

      NDFF = 5
      IF(T2MIN.GE.T2MAX1) THEN
CCC            WRITE(6,*) 'T2MIN=',T2MIN,'.GT.T2MAX1=',T2MAX1
         GOTO 180
      ENDIF
      if(itch.eq.0) then
         NRN = NRN + 1
         T2 = T2MIN * ((T2MAX1/T2MIN)**X(NRN))
C weight for t distribution
         if(wtdist.le.0.) then
            write(6,*) ' partdh: wtdist',wtdist
         endif
         WTDIST=WTDIST*dhformf(-T2)
         WTDIST = WTDIST * T2 *DLOG(T2MAX1/T2MIN)
         if(wtdist.le.0.) then
            write(6,*) ' partdh: t2,T2max1,t2min,form',t2,t2max1,t2min,
     +      dhformf(-T2)
         endif
c       write(6,*) ' partdh: wtdist ',wtdist,t2
      else
         t2 = 0d0
      endif
      IF(ICHECK.EQ.1) THEN
         write(6,*) ' partdh: T2,X(NRN),NRN ',T2,X(NRN),NRN
      ENDIF

c      IF(IWEI.EQ.1) THEN
c          write(6,*) ' partdh: T2,X(NRN),NRN ',T2,X(NRN),NRN
c          ENDIF
      XR1= T2/XR/SSS
      T2=-T2
      T2GKI = SNGL(T2)
      XFGKI = SNGL(XR)
c      write(6,*) ' partdh xpom,xr ',xpom,xr
      IF(DABS(WTDIST).LE.1.D-45) GOTO 190
      COSTP=(1.D0 - XR - XR1*XR)/(1.D0 - XR + XR1*XR)
      IF(COSTP.GT.1.D0) COSTP=1.0D0
      IF(COSTP.LT.-1.D0) COSTP=-1.0D0
C...  COSTP IS SCATTERING ANGLE OF PROTON
      if(itch.eq.0) then
         NRN = NRN + 1
         PHIP=2.D0*PI*X(NRN)
         IF(IGENFL.EQ.0) THEN
c         PHIP=2.D0*PI*X(NRN)
            PHIO = PHIP
         ELSEIF(IGENFL.EQ.1) THEN
            PHIP = PHIO
         ENDIF
      else
         PHIP = 0.D0
         PHIO = PHIP
      endif
c      PHIP = 0.D0
c      IF(IDEBUG.EQ.1) write(6,*) 'partdh: PHIP=',PHIP
      SPHP=DSIN(PHIP)
      CPHP=DCOS(PHIP)

C FINAL STATE PROTON
      NPFIN=NIA1+NDFF
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
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
      IF(BOCHCK.GT.0.99999999D0) goto 210

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
      POMDGA = XR * YX * SSS/2.D0
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
c         write(6,*) ' partdh: T2GKI,T2 ',T2GKI,T2
c         write(6,*) ' partdh: yx=',yx,dot1(nia1,2)/dot1(1,2)
c         write(6,*) ' partdh: xr=',xr,dot1(nia1,nia1+1)/dot1(nia1,2)
c      ENDIF
c for construction of parton in pomeron goto gamma pomeron system
      DBCMS(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBCMS(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBCMS(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBCMS(4)=  P(NIA1,4) + P(NIA1+1,4)
      spom = dot(dbcms,dbcms)
c      IF(IWEI.EQ.1) write(6,*) ' partdh: SPOM= ',spom
      IF(SPOM.LT.0.0D0) GOTO 230
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2
     +              + (DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 210

      CALL DUDBRB(NIA1,NIA1+4,0.D0,0.D0,
     +  -DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +  -DBCMS(3)/DBCMS(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(NIA1,NIA1+4,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(NIA1,NIA1+4,-STHETA,0.d0,0.d0,0.d0,0.d0)
c            if(xpr.ne.xx*xfgki) write(6,*) xpr,xfgki*xx,xfgki,xx
c            if(iwei.eq.1) call HF1(301,sngl(xx),1.)
c            if(iwei.eq.1) write(6,*) 'partdh: beta = ',xx
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
      WTGLU = 0.D0
c           write(6,*) ' partdh: NFLAV ',NFLAVP
      nafl =0
      DO 80 I=1,NFLAVP
         IF(INTER.LT.2) THEN
            IF(SPOM.LT.(4.*ULMASS(I)**2)) GOTO 80
         ELSEIF(INTER.EQ.2) THEN
            IF(I.GE.(NFLAVP+1)) GOTO 80
            IF(SPOM.LT.((ULMASS(I)+ULMASS(I+1))**2)) GOTO 80
         ENDIF
         NAFL = I
   80 CONTINUE
      if(nafl.le.0) then
         write(6,*) ' partdh: spom ',spom,nafl,nflavp,nflav
      endif
      KPA = KPAO
      SCAL2 = SNGL(Q2Q)
      IF(IPRO.EQ.20) THEN
         NF1=NIA1+2
         NF2=NIA1+4
         K(NF1,1)=2
         K(NF1,2)=KPA
         K(NF1,3)=NIA1
         K(NF2,1)=1
         K(NF2,2)=KPA2
         K(NF2,3)=NIA1
         IJOIN(1)=NF1
         IJOIN(3) = NF2
         IJOIN(2) = NF1+1
         NJOIN=3
C
         K(NF1+1,1)=2
         K(NF1+1,2)=KGL
         K(NF1+1,3)=NIA1
      ELSEIF(IPRO.EQ.21) THEN
         NF1=NIA1+2
         NF2=NIA1+3
         K(NF1,1)=2
         K(NF1,2)=KPA
         K(NF1,3)=NIA1
         K(NF2,1)=1
         K(NF2,2)=KPA2
         K(NF2,3)=NIA1
         IJOIN(1)=NF1
         IJOIN(2) = NF2
         NJOIN=2
      ENDIF
      CALL LUJOIN(NJOIN,IJOIN)

      DO 90  I=1,4
         P(NF1,I)=0.0D0
         P(NF1+1,I)=0.0D0
         P(NF2,I)=0.0D0
   90 CONTINUE
      SHAT=DOT(DBCMS,DBCMS)
      IF(SHAT.LE.0.0) THEN
         GOTO 240
      ENDIF

C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
      ECM =DSQRT(SHAT)
      IF(ECM.LE.(AM(1)+AM(2)+AM(3))) GOTO 250
c      write(6,*) ' ECM  ',ECM,AM(1),AM(2)
c          write(6,*) 'partdh:',x(nrn-1),x(nrn),nrn,ndim
      IF(IPRO.EQ.20) THEN
         iph=1
         if(iph.eq.0) then
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
         else
c            XY(1)=draprn()
c            XY(2)=draprn()
c            XY(3)=draprn()
c            XY(4)=draprn()
c            XY(5)=draprn()
            NRN=NRN+1
            XY(1)=X(NRN)
c this is for alpha2
            NRN=NRN+1
            XY(2)=X(NRN)
c this is for mass m**2
            NRN=NRN+1
            XY(3)=X(NRN)
            NRN=NRN+1
            XY(4)=X(NRN)
            NRN=NRN+1
            XY(5)=X(NRN)

         endif
         NP = 3
      ELSEIF(IPRO.EQ.21) THEN
         NRN=NRN+1
         XY(1)=X(NRN)
         NRN=NRN+1
         XY(2)=X(NRN)
         NP = 2
      ENDIF
ccc          write(6,*) 'partdh:',nrn,ndim
      CALL PHASE(NP,ECM,AM,PCM,WT)
      IF(WT.LE.0.D0) GOTO 260
c         write(6,*) ' partdh phase 1sr :',WT,XY(2)*DLOG(1.D0/SIMIN)
c         write(6,*) ' partdh phase AM:',am(1),am(2)
      DO 100 I=1,4
         IF(IPRO.EQ.20) THEN
            P(NF1,I)=PCM(I,2)
            P(NF1+1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,3)
            P(NF1,5)=AM(1)
            P(NF2,5)=AM(2)
         ELSEIF(IPRO.EQ.21) THEN
            P(NF1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,2)
            P(NF1,5)=AM(1)
            P(NF2,5)=AM(2)
         ENDIF

  100 CONTINUE
      if(pcm(1,2).ne.pcm(1,2)) then
         call dulist(1)
      endif
      IF(IPRO.EQ.21) THEN
         PT2 = DPLU(NF2,9)
         CALL CUTG(PT2,NACC)
         IF(NACC.EQ.0) GOTO 270
c get weight correction from d cost --> d k_T**2
c include factor 1/2 because from cost integration two possible states
c count for pt**2 i.e. cost and -cost
c         write(6,*) ' partdh WT,1/8pi ',WT,1/(8*PI)
         WT = (8.D0*PI*DABS(P(NF2,3))*DSQRT(PT2+P(NF2,3)**2))/2.D0
      ELSEIF(IPRO.EQ.20) THEN
         IF(Iqqg.eq.0) Then
            PT2 = DPLU(NF1+1,9)
            CALL CUTG(PT2,NACC)
            IF(NACC.EQ.0) GOTO 270
         Endif
c now the other cuts are moved to eleqqg
      ENDIF
      CALL DUDBRB(NIA1,NIA1+4,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(NIA1,NIA1+4,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(NIA1,NIA1+4,0.D0,0.D0,DBCMS(1)/DBCMS(4),
     +  DBCMS(2)/DBCMS(4),DBCMS(3)/DBCMS(4))
c      write(6,*) ' partdh before last boost'
c       call dulist(1)

C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
c      call DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
c      call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
c      write(6,*) ' end of partdh '
c       call dulist(1)
      XPR = SNGL(XP2)
      WPART=WTDIST*FWEI*WT
      PT2H = SNGL(PT2)
      SHH = SNGL(SHAT)

      NDIMC = NRN
      CALL DUEDIT(12)
      IF(WPART.LT.0.D0) THEN
         WRITE(6,*) ' WPART < 0 '
         WRITE(6,*) ' WTGLU = ',WTGLU,' FWEI = ',FWEI
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
c      write(6,*) 'partdh wpart,q2,xp2,yy,ipro',wpart,q2,xp2,q2q,yy,ipro
c         WRITE(6,*) 'partdh igenfl,imix,iwei ',igenfl,imix,iwei
      MSTJ(24) = MSTJ24
c      if(iwei.eq.1)      write(6,*) kpa,kpa2
c      call lulist(1)
      RETURN
  110 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  120 CONTINUE
      WPART = 0.D0
      IF(IWEI.EQ.1)
     +   write(6,*) ' partdh: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  130 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdh: xlimit ; RETURN xmin,xmax',xmin,xmax,IPRO
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: xlimit ; RETURN xmin,xmax',xmin,xmax,
     +   IPRO
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  140 CONTINUE
      WPART = 0.D0
c      write(6,*) ' partdh: xp2 limit ; RETURN ',xp2,xmax,xmin,IPRO
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: xp2 limit ; RETURN ',xp2,xmax,xmin,IPRO
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  150 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdh: slimit ; RETURN '
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: slimit ; RETURN ',xg1,xp2
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  160 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdh: theta limit ; RETURN '
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: theta limit ; RETURN '
      ENDIF
      IF(IHERAC.EQ.1) THEN
         NTHRE = NTHRE +1
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  170 CONTINUE
      WPART = 0.D0
c      write(6,*) ' partdh: xp2/xx limit ; RETURN ',IPRO
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: xp2/xx limit ; RETURN ',IPRO
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  180 CONTINUE
      WPART = 0.D0
c      write(6,*) ' partdh: tlimit ; RETURN '
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: tlimit ; RETURN ',T2,T2MIN,T2MAX1
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  190 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdh: wtdist limit ; RETURN ',WTDIST
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: wtdist limit ; RETURN ',WTDIST
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  200 CONTINUE
      WPART = 0.D0
c      IF(IWEI.EQ.1)
c     +write(6,*) ' partdh: wtglu limit ; RETURN ',WTGLU
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: wtglu limit ; RETURN '
         write(6,*) ' partdh: WTGLU,nafl,nflavp,spom',
     +   WTGLU,nafl,nflavp,spom
         write(6,*) ' partdh: xr,q2,xbj,t2',xr,q2,xp2,t2
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  210 CONTINUE
      WPART = 0.D0

      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: bochck limit,pom ; RETURN ',BOCHCK,NDIM,
     +   NRN
         write(6,*) ' bochk:   P(NIA1,I)',(P(NIA1,I),I=1,4)
         write(6,*) ' bochk: P(NIA1+1,I)',(P(NIA1+1,I),I=1,4)
         write(6,*) ' bochk: DBCMS ',DBCMS
         write(6,*) ' bochk: X ',X
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  220 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: bochck limit, phase ; RETURN '
         write(6,*) ' bochk:   P(NIA1,I)',(P(NIA1,I),I=1,4)
         write(6,*) ' bochk: P(NIA1+1,I)',(P(NIA1+1,I),I=1,4)
         write(6,*) ' bochk: X ',X
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  230 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: spom=',spom,' test =',
     +   (-Q2+xr*yx*sss+t2)
c         write(6,*) ' partdh: nia1,nia1+1 ',nia1,nia1+1
         write(6,*) ' partdh: yx=',yx,dot1(nia1,2)/dot1(1,2)
         write(6,*) ' partdh: xr=',xr,dot1(nia1,nia1+1)/dot1(nia1,2)
         write(6,*) ' partdh: x_bj = ',xpr
c         write(6,*) 'partdh: 2gamma pom ',2*dot1(nia1,nia1+1),
c     +    xr*yx*sss
         write(6,*) ' partdh: Q2 ',Q2,dot1(nia1,nia1)
         write(6,*) ' partdh: T2 ',T2,T2min,T2max1
c         write(6,*) ' partdh:  test',(-Q2+dot1(nia1,nia1+1)/dot1(nia1,2)
c     +    *yx*sss+t2)
c        call lulist(1)
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  240 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: shat ; RETURN :x_bj',xpr,'x_pom',xr
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  250 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: ECM limit ; RETURN ',ECM,AM(1),AM(2)
         write(6,*) ' partdh: ecm_min ',dsqrt(-q2+xmin*yx*sss),xmin
      ENDIF
      MSTJ(24) = MSTJ24
      NDIMC = 9999
      RETURN
  260 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: PHASE wt=0 ; RETURN '
      ENDIF
      MSTJ(24) = MSTJ24
      NDIMC = 9999
      RETURN

  270 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: PTCUT limit ;RETURN ',PT2,PT2CUT(IPRO),
     +   IPRO
      ENDIF

      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  280 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) '350 partdh: shat;RETURN ',shat,ulmass(211),
     +   ulmass(kpa),nflav,nflavp
         write(6,*) '350 partdh: shat ; RETURN :x_bj',xpr,'x_pom',xr
         write(6,*) '350 partdh: shat,spom',shat,spom,kpa,ipro
      endif
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  290 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: VM_mass check ',IVM,SPOM
      ENDIF
      IF(IHERAC.EQ.1) THEN
         NVMRE = NVMRE +1
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
  300 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdh: w chech '
      ENDIF
      NDIMC = 9999
      MSTJ(24) = MSTJ24
      RETURN
      END
