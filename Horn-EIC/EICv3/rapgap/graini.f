*CMZ :  2.08/04 10/01/2000  10.19.35  by  Hannes Jung
*CMZ :  2.08/02 20/10/99  08.36.03  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.06.21  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  20.57.12  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE GRAINI

*#**********************************************************************
*#
*#    SUBROUTINE GRAINI
*#
*# PURPOSE: INITIALISE RAPGAP PARAMETER COMMONS RAPINI,
*#
*#
*# CHANGED BY:                              AT:
*# REASON:
*#
*#**********************************************************************
      IMPLICIT NONE
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

*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGFULL.
      INTEGER IFULL,IQCDGRID
      COMMON /OALPINI/ IFULL,IQCDGRID
*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEEP,RGSCQ2.
      DOUBLE PRECISION SCALQ2
	Integer ISET,IP2
      COMMON/RESGAM/ISET,IP2,SCALQ2
*KEEP,RGLQ2.
      DOUBLE PRECISION Q2SUPP
      COMMON/LOWQ2S/Q2SUPP
*KEEP,PSHWR.
      INTEGER IORDER,IALPS,ITIMSHR,ISOFTG,ICCFM
      COMMON /PSHWR/IORDER,IALPS,ITIMSHR,ISOFTG,ICCFM
*KEEP,RGPQCDPOM.
      Integer Iqqg
	COMMON/pqcdpom/ Iqqg
C     SAVE


*KEND.
      REAL ALAM
	Integer NF
      COMMON/ALPHA/NF,ALAM
      CHARACTER CHIN*100
	INTEGER IVM
      COMMON/VMESON/IVM
      REAL PYPAR,PYVAR
	Integer IPY
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
	Integer IDEBUG
      COMMON/ INTERN/IDEBUG
	Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      REAL ACC1,ACC2
	Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      DOUBLE PRECISION C1,Cg
      COMMON/BUCHMUE/C1,Cg
*KEEP,RGPRKT.
      DOUBLE PRECISION PRKT1,PRKT2,pktm
	Integer IGAMKT
      COMMON /RGPRKT/PRKT1,PRKT2,pktm,IGAMKT

*KEND.
      DOUBLE PRECISION  XMW
      CALL BLEPIN
      IDEBUG = 0
c select integration procedure
c IINT = 1 DIVON
C IINT = 0 BASES/SPRING
      IINT = 0
      NCB = 60000
      ACC1 = 0.75
      ACC2 = 0.5
c... call block data for heracles
      CALL HS46INI
      NF=4
C FINAL STATE PARTON SHOWER
      IFPS = 3
      ILEPTO = 1
C LEPIN KE = 11
      KE = 11
      KGL = 21
      KPH=22
      KEB = KPH
      KPA = 1
      NFLAV = 5
      NFLQCDC = 3
      IHF = 0
      IALMKT = 0
c for photon structure function
      MSTP(12) = 2
c select proton strcuture function
      MSTP(51)=9
c select pion strcuture function
      MSTP(52)=1
c select photon strcuture function
      MSTP(56)=1
c virtual photon supression for resolved photons
      OMEG2 = 0.01D0
C select interaction type
C INTER = 0 --> photon exchange
C INTER = 2 --> CC (W) exchange
      INTER = 0
C select BGF in semihard approach a la Zotov et al. ISEMIH = 1
      ISEMIH = 0
C RUNNING ALPHA_S IRUNA = 1
      IRUNA=1
C RUNNING ALPHA_EM IRUNAEM = 1
      IRUNAEM=0
C SET ELECTROWEAK PARAMETERS
      SIN2W = 0.23D0
      XMW = 80.D0
      XMW2 = XMW**2
C SELECT SCALE FOR STRUCTURE FKT AND ALPHA_S
C IQ2 =1 SCALE Q2 = MASS**2
C IQ2 =2 SCALE Q2 = SHAT
C IQ2 =3 SCALE Q2 = MASS**2 + PT**2
C IQ2 =4 SCALE Q2 = Q2
C IQ2 =5 SCALE Q2 = Q2 + pt**2
      IQ2=5
C SELECT PROCESS
      IPRO= 12
c qpm and martix elements
      IFULL = 1
c for matrix elements: grid (IQCDGRID = 1)
c       or calculated event by event (IQCDGRID = 0)
      IQCDGRID = 1
c select vm production
      IVM = 0
c select p remnant treatment
      IREM = 1
c select pomeron distribution
c npom = 0 (a la streng hera proc. 1987)
c      = 1 ( a la Ingelman)
      NPOM = -10
C...  GLUON DENSITY IN POMERON : ALA STRENG HERA PROCEEDINGS 1987
C...  IF NG = 0 THEN G(X) = 6(1-X)
C...  FOR 0 <= N <= 5 USE NG=...
      NG=-14
C... PARAMETERS FOR POMERON DISTRIBUTION ALA STRENG HERA PROCEEDINGS
C... R_N **2 = 4.7 GEV**-2
      RN2=4.7D0
C... ALPHA_P = 0.25 GEV**-2
      ALPHP=0.25D0
C... EPSP = 0.085
      EPSP = 0.085D0
      XF=0.90D0
C   TMAX =  ACCEPTANCE THETA = 0.5 DEG
      T2MAX=5.D0
C Minimum Q^2 of electron to be generated
      QMI = 5.d0
C Maximum Q^2 of electron to be generated
      QMA = 10D8
C Minimum y of electron to be generated
      YMI=0.0d0
C Minimum y of electron to be generated
      YMA=1.0d0
C pt^2_hat cut for light quark Matrix Elements
      PT2CUT(10)=5.D0
      PT2CUT(13)=5.D0
      PT2CUT(15)=5.D0
      PT2CUT(18)=5.D0
C Maximium theta angle of scattered electron
      THEMA = 180.0D0
C Minimum  theta angle of scattered electron
      THEMI =   0.0D0
C ELECTRON MOMENTUM
      PLEPIN =-30.
C PROTON MOMENTUM
      PPIN   = 820.
C PERFORM FRAGMENTATION NFRAG=0/1
      NFRAG = 1
c safety for fragmentation
      PYPAR(12) = 1.5
      MSTU(12) = 0
      CHIN='CHAF(C100)=pomeron;'
      CALL LUGIVE(CHIN)
      CHIN='PMAS(4,1)=1.5'
      CALL LUGIVE(CHIN)
c      CALL LUGIVE('PMAS(1,1)=0.450')
c      CALL LUGIVE('PMAS(2,1)=0.450')
c      CALL LUGIVE('PMAS(2101,1)=0.650')
c      CALL LUGIVE('PMAS(2103,1)=0.650')
c      CALL LUGIVE('PMAS(2203,1)=0.650')
c for HERACLES low Q2 supreesion
      Q2SUPP=3.37D0
c for IPRO =14,18 ,1400 which flavor is produced
      IHFLA = 3
c set scale factor for QCD scale:
      SCALFA = 1.D0
c Buchmueller et al parameters
      C1 = 1.D0
      Cg = 1.D0
c for DIS resolved photons .... scale q2
      SCALQ2=1.D0

c
C     Matrix elements for resolved photon processes
C     gg --> qq_bar
      IRPA = 1
C     g + g --> g + g
      IRPB = 1
C     g + q --> g + q
      IRPC = 1
C     qq_bar --> gg
      IRPD = 1
C     q + q_bar --> q + q_bar
      IRPE = 1
C     qq --> qq
      IRPF = 1
C     BFKL qq --> qq
      IRPG = 0

c Schuler Sjostrand
      ISET = 3
      IP2 = 0
c initial state parton shower parameters
      IORDER=1
      IALPS=1
      ITIMSHR=1
      ISOFTG=1
	ICCFM = 0
c intrinsic kt in proton
      IGAMKT = 1 ! for gaussian
c      IGAMKT = 2 ! for dkt^2/(kt0^2 + kt^2)
c      IGAMKT = 3 ! for dkt^2/(kt0^2 + kt^2)^2
c gaussian width:
      PRKT2=0.44D0
c power like: for IGAMKT=2,3
	pktm=2.d0	! maximum kt
c intrinsic kt in photon
      PRKT1=0.77D0
c full/approximate integration for diffr. qqg
      Iqqg = 0
      WRITE(6,*) ' GRAINI: initialisation of RAPGAP '

      RETURN
      END
