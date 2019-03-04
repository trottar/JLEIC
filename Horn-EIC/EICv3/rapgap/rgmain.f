*CMZ :  2.08/04 26/01/2000  15.07.40  by  Hannes Jung
*CMZ :  2.08/03 07/12/99  18.22.41  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  10.08.43  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      PROGRAM RGMAIN
	Implicit None
*KEEP,RGFULL.
      INTEGER IFULL,IQCDGRID
      COMMON /OALPINI/ IFULL,IQCDGRID
*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEND.
C---initialise JETSET 7.3 parameters
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


*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGSCQ2.
      DOUBLE PRECISION SCALQ2
	Integer ISET,IP2
      COMMON/RESGAM/ISET,IP2,SCALQ2
*KEEP,RGLQ2.
      DOUBLE PRECISION Q2SUPP
      COMMON/LOWQ2S/Q2SUPP
*KEND.
      Integer IVM
      COMMON/VMESON/IVM
      REAL PYPAR,PYVAR
      INTEGER IPY
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
C SEQUENCES FOR HERACLES.
      Double Precision THEMIN,CTHMIN,CTHCON
      COMMON /HSTCUT/ THEMIN,CTHMIN,CTHCON
*KEEP,RGHSUNTS.
      Integer LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
*KEEP,RGHSOPTN.
      Integer INT2,INT3,ISAM2,ISAM3,IOPLOT,IPRINT,ICUT
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     +                IOPLOT,IPRINT,ICUT

*KEEP,RGHSCUTS.
      DOUBLE PRECISION XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN


*KEEP,RGHSVGLP.
      Integer         NPOVEG,NUMINT,NPHYP
      COMMON /HSVGLP/ NPOVEG,NUMINT,NPHYP

*KEND.
      Integer LPAR,LPARIN,IPART
      COMMON /HSPARL/ LPAR(20),LPARIN(12),IPART
      REAL*4 PAR111,PAR112,PARL11,PARL19
      INTEGER MST111,MST115
	COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
C-MH weighting COMMONs:

      Integer I,ISEED
C---initialise JETSET 7.3 parameters
      CALL GARINI
      CALL GJEINI
C     initialize random number generator
      ISEED = 213123
      CALL H1RNIN(ISEED)
C---initialise RAPGAP parameters
      CALL GRAINI
C-- change standard parameters
C**********************************************************************
C     1ST INCOMING PARTICLE  (KE=11 ELECTRON)
C     1ST INCOMING PARTICLE  (KE=22 PHOTON)
      KE =  -11
c      KE =  22
C     LEPTON MOMENTIM (D=-30)
      PLEPIN = -27.5
C**********************************************************************
C     2ND INCOMING PARTICLE (KP = 2212 PROTON)
      KP = 2212
C     proton momentum (D=820)
      PPIN   = 820.
C**********************************************************************
C     PARTON SHOWER OFF/ON (D=10)
C     (off = 0, =1 initial state, =2 final state, =3 inital + final state
C       = 10 ariadne)

      IFPS   = 1
C**********************************************************************
C     fragmentation on/off (D=1)
C     NFRAG=1         (D=1) FRAGMENATION
C         =10 FRAGMENATION+ p dissocia
      NFRAG  = 0
C     IREM=1   p dissociation treatment
C              (D=1): P(beta') = 2(1-beta')
C         =2        :          = (a+1)(1-beta')**a
      IREM = 1
C**********************************************************************
C     fixed/running alpha_s (D=1)
      IRUNA  = 1
C     alphas parameters
C     nr of flavours wrt lambda_QCD
C     MSTU(112) = 5
C     lambda QCD
C     PARU(112) = 0.25
C     lambda QCD for MRSR2
C     PARU(112) = 0.239
C     min nr of flavours for alphas
      MSTU(113) = 3
C     max nr of flavours for alphas
      MSTU(114) = 5
C**********************************************************************
C      scale for alpha_s
C      1:  q2 = m**2
C      2:  q2 = shat
C      3:  q2 = m**2 + pt**2
C      4:  q2 = Q2
C      5:  q2 = Q2 + pt**2
      IQ2    = 2
C**********************************************************************
C      scale factor q2 = SCALFA * q2
      SCALFA = 1.D0
C**********************************************************************
C     select process to be generated
c      WRITE(6,*) ' select process IPRO and IFULL,IDIR '
c      READ(5,*) IPRO,IFULL,IDIR
c      WRITE(6,*) ' IPRO = ',IPRO,' and event mixing IFULL = ',IFULL
c      write(6,*) ' IDIR = ',IDIR,' selected '
      IDIR = 1
      IPRO = 12
      IFULL = 0
C**********************************************************************
C  Processes gamma gluon fusion using EPA and gamma glu Martix Element
C    10: gamma pommeron_gluon --> q qbar (light quarks:u,d,s)
C    11:  gamma pommeron_gluon --> Q Qbar (charm quarks)
C
C**********************************************************************
C  Process e q --> e' q'
C    12:  e q --> e' q'
C    1200:  e q --> e' q' using HERACLES
C**********************************************************************
C  Processes gamma gluon fusion using full e glu Matrix Element
C    13:  e gluon --> e' q qbar (light quarks:u,d,s)
C    14:  e gluon --> e' Q Qbar (charm quarks)
C    1400:  e gluon --> e' Q Qbar (charm quarks) using HERACLES
C    18:  resolved gamma processes
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
c Parameters for SaS (Schuler Sjostrand) virtual photon structure function
      ISET = 3
      IP2 = 0
C**********************************************************************
C  Process gamma* pomeron --> rho pomeron
C   100: e pomeron --> e' rho pomeron
C
C**********************************************************************
C    parameters for pomeron distribution
C     gluon density in pomeron (D=0)
C     (D=0) gluon density in pomeron
C     <0 user supplied via SUBROUTINE USDIFFR
C     20 = pi+- exchange
C     21 = pi0 exchange
C     30 = Nikolaev Zakharov model
C     40 = hard pomeron M.Wuestoff
C     41 = hard pomeron  Bartels,Lotter,Wuesthoff
C     42 = 2 - glu pomeron(soft)  M.Diehl
C     45 = Buchmueller/McDermott/Hebecker

      NG   =  0
C     if MSTP(52) < 10 use inbuild structure functionsof PYTHIAfor pi
C if PDFLIB is used MSTP(51) gives
C     10^6*NPTYPE+1000*NGROUP+NSET    is used as coding scheme
C example 2001001 for NPTYPE 2 Ngroup 1 Nset 1
C      MSTP(52) = 2001001
C proton structure function
      MSTP(51) = 9
C
C
C     Q2SUPP=3.37  exponential supression factor for small Q2
C               in parton densities for HERACLES
C               =5 for Q2 > 1
C               =10 for Q2 >0.5
      Q2SUPP = 3.37

C
C     NFLAV=5 nr of flavours used in str.fct.
      NFLAV = 4
C
c photon structure function
c MSTP(56) = 1 GRS parametristaion
c          = 2 SaSgam (Schuler-Sjostrand)
c          > 1000 PDFLIB with Drees Godbole virtual gamma suppression
      MSTP(56) = 1
C
C     SCALQ2=1. scale/Q2 for resolved gamma in DIS
      SCALQ2 = 1.0
C
C select pomeron distribution
C    (D=0) pomeron density
C  =0 streng density ,
C  =1 Ingelman density,
C  =2 Donnachie Landshoff density,
C  =20 pion exchange
C  =30 Nikolaev Zakharov model
C  =40 hard pomeron M.Wuestoff
C  =41 hard pomeron  Bartels,Lotter,Wuesthoff
C  =42 2 - glu pomeron(soft)  M.Diehl
C  =45 Buchmueller/McDermott/Hebecker
C  <0 user supplied via SUBROUTINE USDIFFR
      NPOM = 0
C epsilon for rising of x Section in Streng and DL pomeron distribution
      EPSP = 0.085
      RN2  = 4.0
      ALPHP = 0.25
C cut of x_f
      XF  = 0.9
c Maximum -t allowed in generation
      T2MAX = 1.
C
C     IALMKT=0   ( no primordial kt for partons in diffraction)
C           =1 ( primordial kt for partons for IPRO=12
C                 acc. aligned jet model)
      IALMKT=0
C
c select Vector meson production
c IVM = 1 light VM
C IVM = 443 J/psi
C IVM = 553 Upsilon
C IVM = 0 no special selection
      IVM = 0
C Minimum Q^2 of electron to be generated
      QMI = 5.0d0
C Minimum y of electron to be generated
      YMI=0.04d0
C pt^2_hat cut for light quark Matrix Elements
      PT2CUT(10)=4.
      PT2CUT(13)=4.
      PT2CUT(15)=4.
      PT2CUT(18)=4.
C Maximium theta angle of scattered electron
      THEMA = 180.0D0
C Minimium theta angle of scattered electron
      THEMI =   0.0D0
c select heavy flavour code for IPRO=14,1400
c      IHFLA=5
*                                        1=cuts in x and lower Q2
*                                        2=cuts in x, lower Q2, W
*                                        3=cuts in x, y, Q2, W
*                                        must be 3 for the use of THMI!
*! (D=3) definition of cuts applied
      ICUT=3
*! (D=1) upper cut in x
      XMAX=1.0
*! (D=1E-4) lower cut in x
      XMIN=0.00001
*! (D=0) integration of the non-
      INT2(1) = 1
*                                        radiative contribution (NC)
*                                        0=off, 1=on
*! (D=0) integration of the non-
      INT2(2) = 0
*                                        radiative contribution (CC)
*                                        0=off, 1=on
*! (D=0) integration of initial
      INT3(2) = 0
*                                        state bremsstr. contr. NC
*                                        number of iterations for VEGAS
*                                        <100:VEGAS, >100:VEGAS1/2
*! (D=0) as INT3(1) but for final
      INT3(1) = 0
*                                        state bremsstr. NC
*! (D=0) as INT3(1) but for Compton
      INT3(2) = 0
*                                        contribution NC
*! (D=0) integration of initial
      INT3(3) = 0
*                                        state bremsstr. contr. CC
*                                        number of iterations for VEGAS
*                                        <100:VEGAS, >100:VEGAS1/2
*! (D=1000) number of integration
      NPOVEG=50000
*                                        points used for VEGAS
*! (D=1) only Born x-section term
      LPARIN(2) = 0
*                                        or 1=including corrections
*! (D=0) sampling of the non-
      ISAM2(1) = 1
*                                        radiative contribution (NC)
*                                        0=off, 1=on
*      ISAM2(2) = 0                  ! (D=0) sampling of the non-
*                                        radiative contribution (CC)
*                                        0=off, 1=on
*! (D=0) sampling of initial
      ISAM3(1) = 0
*                                        state bremsstrahlung  NC
*! (D=0) sampling of final
      ISAM3(2) = 0
*                                        state bremsstr. NC
*! (D=0) sampling of Compton
      ISAM3(3) = 0
*                                        contribution NC
*      ISAM3(7) = 0                   ! (D=0) sampling of initial
*                                        state bremsstrahlung  CC

C**********************************************************************
C Initialize ARIADNE
      CALL ARINIT('RAPGAP')
C--- CALCULATE X SECTION
      CALL RAPGAP
C--- print x section
      CALL RAEND(1)
C--- event generation
      DO 10 I=1,50000
         CALL EVENT
C--- user analysis routine
         CALL ANALYS
C---
   10 CONTINUE
C---PRINT NR OF GENERATED EVENTS
      CALL RAEND(20)
      STOP
      END
