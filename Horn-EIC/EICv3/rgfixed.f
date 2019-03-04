      PROGRAM RGMAIN

      IMPLICIT NONE

CPERvvvvvvvvvvvvvvvvvvvvvv
      INTEGER nPrintedA, nPrintedB, luntmp
CPER^^^^^^^^^^^^^^^^^^^^^^

      INTEGER mxEvent
C      PARAMETER (mxEvent = 20000000)
C      PARAMETER (mxEvent = 500000000)
      PARAMETER (mxEvent = 200)

      INTEGER nhbmem, lrecl, istat, icycle
      PARAMETER (nhbmem=10000000)
      REAL hmem(nhbmem)
      COMMON /pawc/ hmem

      CHARACTER ntFileName*80, histFileName*80, fileName*70
      INTEGER lenocc
      EXTERNAL lenocc

      INTEGER iE
      REAL phiE

      CHARACTER chtmp*80

      REAL fixed
      EXTERNAL fixed

      INTEGER iquest(100)
      COMMON /quest/ iquest

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

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
C-MH COMMONs included for weighting:
*KEEP,RGHSUNTS.
	Integer LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND,LUNPD6,LUNPD7
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

C**** Begin ****

      call hlimit(nhbmem)

      call hcdir('//PAWC', ' ')


C---initialise JETSET 7.3 parameters
      CALL GARINI
      CALL GJEINI
C     initialize random number generator
      ISEED = 213123
      CALL H1RNIN(ISEED)
      call readSeed
C---initialise RAPGAP parameters
C *PER* note that some/many of these parameters are changed 
C *PER* below, including, e.g. plepin
      CALL GRAINI
C-- change standard parameters
C**********************************************************************
C     1ST INCOMING PARTICLE  (KE=11 ELECTRON)
C     1ST INCOMING PARTICLE  (KE=-11 positron)
C     1ST INCOMING PARTICLE  (KE=22 PHOTON)
C      KE =  -11
      KE =  11
c      KE =  22
C     LEPTON MOMENTIM (D=-30)
CPER      PLEPIN = -27.5
      PLEPIN = -4.0
C**********************************************************************
C     2ND INCOMING PARTICLE (KP = 2212 PROTON)
      KP = 2212
C     proton momentum (D=820)
CPER      PPIN   = 820.
C*boost*      PPIN   = 10.
C*boost*
C*boost*      call getenv('PLEPIN', chtmp)
C*boost*      read(chtmp,*,err=100) plepin
C*boost*      plepin = -abs(plepin)
C*boost*      GOTO 101
C*boost* 100  CONTINUE
C*boost*      plepin = -4.0
C*boost* 101  CONTINUE
C*boost*
C*boost*      call getenv('PPIN', chtmp)
C*boost*      read(chtmp,*,err=102) ppin
C*boost*      GOTO 103
C*boost* 102  CONTINUE
C*boost*      ppin = 10.0
C*boost* 103  CONTINUE

      plepin = -1.0
      ppin = fixed(plepin)

C**********************************************************************
C     PARTON SHOWER OFF/ON (D=3)
C     (off = 0, =1 initial state, =2 final state, =3 inital + final state
C       = 10 ariadne)
      ifps = 3
C**********************************************************************
C     fragmentation on/off (D=1)
C     NFRAG=1         (D=1) FRAGMENATION
C         =10 FRAGMENATION+ p dissocia
      NFRAG  = 1
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
      PARU(112) = 0.25
C     lambda QCD for MRSR2
C     PARU(112) = 0.239
C#PER#NoInfo     min nr of flavours for alphas
C#PER#NoInfo      MSTU(113) = 3
C#PER#NoInfoC     max nr of flavours for alphas
C#PER#NoInfo      MSTU(114) = 5
C**********************************************************************
C      scale for alpha_s
C      1:  q2 = m**2
C      2:  q2 = shat
C      3:  q2 = m**2 + pt**2
C      4:  q2 = Q2
C      5:  q2 = Q2 + pt**2
      IQ2    = 5
C**********************************************************************
C      scale factor q2 = SCALFA * q2
      SCALFA = 1.D0
C**********************************************************************
C     select process to be generated
c      WRITE(6,*) ' select process IPRO and IFULL,IDIR '
c      READ(5,*) IPRO,IFULL,IDIR
c      WRITE(6,*) ' IPRO = ',IPRO,' and event mixing IFULL = ',IFULL
c      write(6,*) ' IDIR = ',IDIR,' selected '
      IDIR = 0 ! PER diffractive and pion process => idir = 0
      IPRO = 12
      IFULL = 0
C PER perhaps we should set idisdif here as well.  The default, 
C PER idisdif = 0 forces only  process chosen by idir to be generated.
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
C#PER#NoInfoc Parameters for SaS (Schuler Sjostrand) virtual photon structure function
C#PER#NoInfo      ISET = 3
C#PER#NoInfo      IP2 = 0
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
C *PER* This is important.  it determines the internal structure 
C *PER* of the pion being probed
      NG   =  20
C     if MSTP(52) < 10 use inbuild structure functionsof PYTHIAfor pi
C if PDFLIB is used MSTP(51) gives
C     10^6*NPTYPE+1000*NGROUP+NSET    is used as coding scheme
C example 2001001 for NPTYPE 2 Ngroup 1 Nset 1
C      MSTP(52) = 2001001
C proton structure function
      MSTP(51) = 9

C     Q2SUPP=3.37  exponential supression factor for small Q2
C               in parton densities for HERACLES
C               =5 for Q2 > 1
C               =10 for Q2 >0.5
      Q2SUPP = 10.0

C     NFLAV=5 nr of flavours used in str.fct.
      NFLAV = 4
C
c photon structure function
c MSTP(56) = 1 GRS parametristaion
c          = 2 SaSgam (Schuler-Sjostrand)
c          > 1000 PDFLIB with Drees Godbole virtual gamma suppression
      MSTP(56) = 1

C     SCALQ2=1. scale/Q2 for resolved gamma in DIS
      SCALQ2 = 1.0
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
      npom = 0
      IF (ng .eq. 20) npom = 20
      IF (ng .eq. 21) npom = 21
C#PER#C epsilon for rising of x Section in Streng and DL pomeron distribution
C#PER#      EPSP = 0.085
C#PER#      RN2  = 4.0
C#PER#      ALPHP = 0.25
C cut of x_f
CPER      XF  = 0.9
      XF  = 0.01
C#PER#c Maximum -t allowed in generation
C#PER#      T2MAX = 1.
C#PER#C
C#PER#C     IALMKT=0   ( no primordial kt for partons in diffraction)
C#PER#C           =1 ( primordial kt for partons for IPRO=12
C#PER#C                 acc. aligned jet model)
C#PER#      IALMKT=0
C#PER#C
C#PER#c select Vector meson production
C#PER#c IVM = 1 light VM
C#PER#C IVM = 443 J/psi
C#PER#C IVM = 553 Upsilon
C#PER#C IVM = 0 no special selection
C#PER#      IVM = 0
C Minimum Q^2 of electron to be generated
      QMI = 1.0d0
C Minimum y of electron to be generated
      YMI=0.004d0
C#PER#C pt^2_hat cut for light quark Matrix Elements
C#PER#      PT2CUT(10)=4.
C#PER#      PT2CUT(13)=4.
C#PER#      PT2CUT(15)=4.
C#PER#      PT2CUT(18)=4.
C Maximium theta angle of scattered electron
      THEMA = 180.0D0
C Minimium theta angle of scattered electron
      THEMI =   0.0D0
C#PER#c select heavy flavour code for IPRO=14,1400
C#PER#c      IHFLA=5
C#PER#*                                        1=cuts in x and lower Q2
C#PER#*                                        2=cuts in x, lower Q2, W
C#PER#*                                        3=cuts in x, y, Q2, W
C#PER#*                                        must be 3 for the use of THMI!
C#PER#NoInfo*! (D=3) definition of cuts applied
C#PER#NoInfo      ICUT=3
*! (D=1) upper cut in x
      XMAX=1.0
*! (D=1E-4) lower cut in x
      XMIN=0.00001
C#PER#NoInfo*! (D=0) integration of the non-
C#PER#NoInfo      INT2(1) = 1
C#PER#NoInfo*                                        radiative contribution (NC)
C#PER#NoInfo*                                        0=off, 1=on
C#PER#NoInfo*! (D=0) integration of the non-
C#PER#NoInfo      INT2(2) = 0
C#PER#NoInfo*                                        radiative contribution (CC)
C#PER#NoInfo*                                        0=off, 1=on
C#PER#NoInfo*! (D=0) integration of initial
C#PER#NoInfo      INT3(2) = 0
C#PER#NoInfo*                                        state bremsstr. contr. NC
C#PER#NoInfo*                                        number of iterations for VEGAS
C#PER#NoInfo*                                        <100:VEGAS, >100:VEGAS1/2
C#PER#NoInfo*! (D=0) as INT3(1) but for final
C#PER#NoInfo      INT3(1) = 0
C#PER#NoInfo*                                        state bremsstr. NC
C#PER#NoInfo*! (D=0) as INT3(1) but for Compton
C#PER#NoInfo      INT3(2) = 0
C#PER#NoInfo*                                        contribution NC
C#PER#NoInfo*! (D=0) integration of initial
C#PER#NoInfo      INT3(3) = 0
C#PER#NoInfo*                                        state bremsstr. contr. CC
C#PER#NoInfo*                                        number of iterations for VEGAS
C#PER#NoInfo*                                        <100:VEGAS, >100:VEGAS1/2
*! (D=1000) number of integration
      NPOVEG=mxEvent
      npoveg=10
C#PER#*                                        points used for VEGAS
C#PER#*! (D=1) only Born x-section term
C#PER#      LPARIN(2) = 0
C#PER#*                                        or 1=including corrections
C#PER#NoInfo*! (D=0) sampling of the non-
C#PER#NoInfo      ISAM2(1) = 1
C#PER#NoInfo*                                        radiative contribution (NC)
C#PER#NoInfo*                                        0=off, 1=on
C#PER#NoInfo*      ISAM2(2) = 0                  ! (D=0) sampling of the non-
C#PER#NoInfo*                                        radiative contribution (CC)
C#PER#NoInfo*                                        0=off, 1=on
C#PER#NoInfo*! (D=0) sampling of initial
C#PER#NoInfo      ISAM3(1) = 0
C#PER#NoInfo*                                        state bremsstrahlung  NC
C#PER#NoInfo*! (D=0) sampling of final
C#PER#NoInfo      ISAM3(2) = 0
C#PER#NoInfo*                                        state bremsstr. NC
C#PER#NoInfo*! (D=0) sampling of Compton
C#PER#NoInfo      ISAM3(3) = 0
C#PER#NoInfo*                                        contribution NC
C#PER#NoInfo*      ISAM3(7) = 0                   ! (D=0) sampling of initial
C#PER#NoInfo*                                        state bremsstrahlung  CC

C**********************************************************************
C Initialize ARIADNE
      CALL ARINITRAPGAP('RAPGAP')
C--- CALCULATE X SECTION
      CALL RAPGAP
C--- print x section
      CALL RAEND(1)

CPER      call histCross

C--- event generation

      filename = ' '
      call getenv('RGFIXED_OUTPUT', filename)
      if (filename(1:1) .eq. ' ') then
        filename = 'rapgap'
      endif
      ntFileName = filename(1:lenocc(filename))//'.nt'

C      iquest(10)=1024000
C      lrecl = 8192
C      call hropen(11, 'rgnt', ntFileName, 'n', lrecl, istat)

      iE = 4
      DO 10 I=1,mxEvent

         write(6,*) 'beginning of event loop', plepin

C         call saveSeed

         CALL EVENT

C        boost events to fixed proton frame
         call fixedBoost

C--- user analysis routine

         CALL ANALYS
C---
         write(6,*) 'at end of event loop'

   10 CONTINUE

      call reject(-1)

C      call hrout(1, icycle, ' ')
C      call hrendc('rgnt')
C      close(11)

      call hdelet(1)

      histFileName = filename(1:lenocc(filename))//'.his'
      lrecl = 1024
      iquest(10) =  64000
      call hropen(10, 'rapgap', histFileName, 'n', lrecl, istat)
      call hrout(0, icycle, ' ')
      call hrendc('rapgap')
      close(10)

C---PRINT NR OF GENERATED EVENTS
      CALL RAEND(20)

      call saveSeed

      STOP
      END

C
      SUBROUTINE readSeed

      IMPLICIT NONE

      CHARACTER mcSeedFile*1024
      CHARACTER command*1024
      INTEGER ivec(25)
      INTEGER iscan
      EXTERNAL iscan
      INTEGER index
      INTEGER lunsee
      PARAMETER (lunsee=9)

C**** Begin ****

      call getenv('RGFIXED_MCSEED',mcSeedFile)
      index = iscan(mcSeedFile, ' ')-1
      IF ((0 .lt. index) .and. (index .lt. 256)) THEN
      ELSE
        mcSeedFile = 'mcseed.dat'
      ENDIF

      OPEN(LUNSEE,FILE=mcSeedFile,STATUS='OLD', err=100)
      call h1rnv(ivec,1)
      READ(LUNSEE,9877) ivec
 9877 FORMAT(20X,25Z9)
      CLOSE(UNIT=LUNSEE)
      call h1rniv(ivec)

      RETURN

 100  continue
      index = iscan(mcSeedFile, ' ')-1
      write(6,*) 'Seed file, ',
     &            mcSeedFile(1:index), ' does not exist'

      return

      END
C
      SUBROUTINE saveSeed

      IMPLICIT NONE

      CHARACTER mcSeedFile*1024
      CHARACTER command*1024
      INTEGER ivec(25)
      INTEGER iscan
      EXTERNAL iscan
      INTEGER index
      INTEGER lunsee
      PARAMETER (lunsee=9)

C**** Begin ****

      call getenv('RGFIXED_MCSEED',mcSeedFile)
      index = iscan(mcSeedFile, ' ')-1
      IF ((0 .lt. index) .and. (index .lt. 256)) THEN
        command = 'mv '//mcSeedFile(1:index)//' '//
     &            mcSeedFile(1:index)//'.bak'
        call system(command)
      ELSE
        mcSeedFile = 'mcseed.dat'
      ENDIF
      call h1rnsv(ivec)
      open(lunsee,file=mcSeedFile)
      write(lunsee,9876) ivec
 9876 format(' Random number seed ',25z9)
      close(unit=lunsee)

      RETURN

      END
C
      SUBROUTINE histCross

      IMPLICIT none

      REAL vector(1)

*KEEP,RGEFFIC.
      DOUBLE PRECISION AVGI,SD
      INTEGER NIN,NOUT
      COMMON/EFFIC/AVGI,SD,NIN,NOUT
C     SAVE

C**** Begin ****

C     Book histogram to store cross section in
      call hbook1(2, 'Cross Section', 1, 0.0, 1.0, 0.0)

      vector(1) = sngl(avgi)
      call hpak(2, vector)

      vector(1) = sngl(sd)
      call hpake(2,sd)

      RETURN

      END
C
      SUBROUTINE fixedBoost

      IMPLICIT none

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

      INTEGER iPart, iDir

      REAL*8 s(4), pIn(4), pOut(4)

      INTEGER iTest, mxTest
      DATA iTest /0/
C      PARAMETER (mxTest = 5)
      PARAMETER (mxTest = 0)
      SAVE iTest

C**** Begin ****

      DO iDir = 1, 4
        s(iDir) = dp(2,iDir)
      ENDDO

      DO iPart = 1, n

        DO iDir = 1, 4
          pIn(iDir) = dp(iPart,iDir)
        ENDDO

        call dloren4(s, pIn, pOut)

        DO iDir = 1, 4
          dp(iPart,iDir) = pOut(iDir)
          p(iPart,iDir) = sngl(pOut(iDir))
        ENDDO

      ENDDO

      IF (iTest .le. mxTest) THEN
        iTest = iTest + 1
        call lulist(1)
      ENDIF

      RETURN

      END
C
      SUBROUTINE DLOREN4  (DIR,P4IN,P4OUT)
C 2 Mar 2000 PER modified for complete double precision

      DOUBLE PRECISION PCM2, ONMCM, EPBETA, PROD

      DOUBLE PRECISION DIR(4),P4IN(4),P4OUT(4)
C
C--                VN(A) MEANS N-VECTOR A
C--                GAMMA=ECM/MCM
C--                EPBETA=ECM*V3(PCM)*V3(BETA)
C--                V3(BETA)=V3(PCM)/ECM
C--

      PCM2=DIR(1)*DIR(1)+DIR(2)*DIR(2)+DIR(3)*DIR(3)
      ONMCM=1.D0/ dSQRT (DIR(4)*DIR(4)-PCM2)

      EPBETA=P4IN(1)*DIR(1)+P4IN(2)*DIR(2)+P4IN(3)*DIR(3)
      PROD=EPBETA*(DIR(4)*ONMCM-1.D0)/PCM2-P4IN(4)*ONMCM
      P4OUT(4)=ONMCM*(P4IN(4)*DIR(4)-EPBETA)
         DO 50 I=1,3
   50 P4OUT(I)=P4IN(I)+DIR(I)*PROD
      RETURN
      END
C
      REAL FUNCTION fixed(pE)

C     This routine calculates what energy the proton needs to be in the
C     reaction p + e given pE momentum of the electron to have, when 
C     boosted to the proton rest frame, and electron momentum of 
C     pFixed.  Rather than doing the calculation explicitly, I 
C     iteratively find the solution (easier on me to save the math 
C     by hand).

      IMPLICIT none

      REAL pE

      REAL pFixed
C      PARAMETER (pFixed = -27.5)
C      PARAMETER (pFixed = -12.0)
C      PARAMETER (pFixed = -6.0)


      REAL*8 limit
      PARAMETER (limit = 1.0e-12)

      REAL*8 pP, eP, pEIn, eEin, pEOut, eEOut

      REAL*8 aMassE, aMassP
      PARAMETER (aMassE = 0.000510998902)
      PARAMETER (aMassP = 0.93827200)

      REAL*8 step, dir

      INTEGER i

      REAL*8 beta, gamma

      INTEGER iterations, maxIteration
      PARAMETER (maxIteration = 5000)

      CHARACTER chtmp*80

C**** Begin ****

      call getenv('RGFIXED_PE', chtmp)
      read(chtmp,*,err=100) pFixed
      pFixed = -abs(pFixed)
      GOTO 101
 100  CONTINUE
      pFixed = -6.0
 101  CONTINUE

      pEIn = pE
      eEIn = sqrt(aMassE**2 + pE**2)

      pP = 0.0
      eP = sqrt(aMassP**2 + pP**2)

      beta = pP / eP
      gamma = 1 / sqrt(1 - beta**2)

      pEOut = gamma * (pEIn - beta * eEIn)
      eEOut = gamma * (eEIn - beta * pEIn)

      step = 1.0
      dir = -1.0

      iterations = 0

      DO WHILE ( (abs(pEOut - dble(pFixed)) .ge. limit) .and.
     &           (iterations .le. maxIteration) )
        IF ( (pEOut .le. pFixed) .and. 
     &       (dir .ge. 0.0) ) THEN
          step = step / 2.0
          dir = -1.0
        ELSE IF ( (pEOut .ge. pFixed) .and. 
     &       (dir .le. 0.0) ) THEN
          step = step / 2.0
          dir = +1.0
        ENDIF

        pP = pP + dir * step
        eP = sqrt(aMassP**2 + pP**2)

        beta = pP / eP
        gamma = 1 / sqrt(1 - beta**2)

        pEOut = gamma * (pEIn - beta * eEIn)
        eEOut = gamma * (eEIn - beta * pEIn)

C        iterations = iterations + 1
C        write(6,*) iterations, pFixed, pE, pEOut, pP, eP

      ENDDO

      IF (iterations .ge. maxIteration) THEN
        write(6,*) 'Maximum iterations reached, pE = ', pE
        STOP
      ENDIF

      fixed = sngl(pP)

      RETURN

      END


*CMZ :  2.06/51 10/10/98  10.29.56  by  Hannes Jung
*CMZ :  2.06/34 20/05/98  09.12.48  by  Hannes Jung
*CMZ :  2.05/05 20/03/97  11.13.39  by  Hannes Jung
*CMZ :  2.04/00 10/12/96  11.15.20  by  Hannes Jung
*CMZ :  2.03/03 22/08/96  19.18.54  by  Hannes Jung
*CMZ :  2.01/09 10/02/96  19.20.34  by  Hannes Jung
*CMZ :  2.00/30 29/12/95  14.41.47  by  Hannes Jung
*CMZ :  2.00/23 19/11/95  20.13.16  by  Hannes Jung
*CMZ :  2.00/19 06/11/95  08.35.01  by  Hannes Jung
*CMZ :  1.04/00 11/01/95  14.38.48  by  Hannes Jung
*CMZ :  1.03/01 03/04/94  16.45.58  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE ARINITRAPGAP(MODE)

C...ARiadne subroutine INITialize
C...Initializes Ariadne to run with other (Lund) MC programs
      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     +                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     +                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     +                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     +                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     +                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

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

*KEND.

*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.

*KEEP,RGLUDAT2.
      REAL PMAS,PARF,VCKM
      INTEGER KCHG
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
C      SAVE

*KEND.
      DOUBLE PRECISION PARL
      COMMON /RAPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,XQ2,U
      SAVE /RAPTOU/
*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEND.

      COMMON /PYSUBS/ MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      SAVE /PYSUBS/
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
      CHARACTER MODE*(*)


C...Set output files if not already done
      IF(MSTA(7).LT.0) MSTA(7)=MSTU(11)
      IF(MSTA(8).LT.0) MSTA(8)=MSTU(11)

C...Write out header
      WRITE(MSTA(7),10000)
      MSTA(2)=1

C...If Ariadne mode, do nothing special
      IF(MODE(1:6).EQ.'ARIADN') THEN
         MSTA(1)=0

C...If JETSET mode, switch off cascade and fragmentation in JETSET
      ELSEIF(MODE(1:6).EQ.'JETSET') THEN
         MSTA(1)=1
         MSTA(5)=MIN(MAX(MSTJ(105),0),1)
         MSTJ(101)=5
         MSTJ(41)=0
         MSTJ(105)=0
         WRITE(MSTA(7),10100)

C...If PYTHIA mode, switch off cascades and fragmentation. Check that
C...Ariadne can handle selected processes
      ELSEIF(MODE(1:6).EQ.'PYTHIA') THEN

         MSTA(1)=2
         WRITE(MSTA(7),10200)
         MSTA(5)=MIN(MAX(MSTP(111),0),1)
         MSTP(61)=0
         MSTP(71)=0
         MSTP(111)=0

C...If LEPTO mode, switch off cascades and fragmentation.
      ELSEIF(MODE(1:5).EQ.'LEPTO') THEN
         MSTA(1)=3
         WRITE(MSTA(7),10300)
         LST(8)=9
         MSTA(5)=MIN(MAX(LST(7),0),1)
         LST(7)=0

C...If RAPGAP mode, switch off cascades and fragmentation.
      ELSEIF(MODE(1:6).EQ.'RAPGAP') THEN
         MSTA(1)=4
         WRITE(MSTA(7),10400)
         MSTA(5)=MIN(MAX(LST(7),0),1)

C...Warn if mode is none of the above
      ELSE
         WRITE(MSTA(7),10500) MODE
         MSTA(1)=0
      ENDIF

C...Set quark masses
      IF(MSTA(24).GT.0) THEN
         DO 10 I=1,8
            PQMAS(I)=PMAS(I,1)
   10    CONTINUE
      ENDIF

      IF(MSTA(24).GE.2) THEN
         DO 20 I=1,5
            PQMAS(I)=PARF(100+I)
   20    CONTINUE
      ENDIF

      IF(MSTA(3).EQ.1) CALL ARTUNE('EMC')

10000 FORMAT(/,14X,
     +     'The Lund Monte Carlo - Ariadne version 4 revision 03',/,
     +     23X,'Latest date of change: 16 Jul 1992')
10100 FORMAT(18X,'Initialization done for running with JETSET')
10200 FORMAT(18X,'Initialization done for running with PYTHIA')
10300 FORMAT(18X,'Initialization done for running with LEPTO')
10400 FORMAT(18X,'Initialization done for running with RAPGAP')
10500 FORMAT(/,15X,'WARNING: Ariadne cannot be initialized for "',A,
     +     '".',/,21X,'Using default initialization instead.')
      RETURN

C**** END OF ARINIT ****************************************************
      END
