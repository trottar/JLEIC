*CMZ :  2.08/05 28/03/2000  17.43.49  by  Hannes Jung
*CMZ :  2.08/04 10/01/2000  11.03.01  by  Hannes Jung
*CMZ :  2.08/02 13/08/99  12.44.05  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.12  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  09.38.43  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE RAPGAP
      IMPLICIT None
      Double Precision XL,XU,ACC,XLB,XUB
      Integer NDI,NCALL,ITMX,NPRN,MXDIM,NDIMB,NWILDB,IG,NCALLB
      COMMON/BVEG1/XL(10),XU(10),ACC,NDI,NCALL,ITMX,NPRN
c common for bases/spring 5.1
      PARAMETER (MXDIM = 50 )
      COMMON /BPARM1/ XLB(MXDIM),XUB(MXDIM),NDIMB,NWILDB,
     +   IG(MXDIM),NCALLB
      Integer ITMX1B,ITMX2B
      Double Precision ACC1B,ACC2B
      COMMON /BPARM2/ ACC1B,ACC2B,ITMX1B,ITMX2B
*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

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


*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEEP,RGEFFIC.
      DOUBLE PRECISION AVGI,SD
      INTEGER NIN,NOUT
      COMMON/EFFIC/AVGI,SD,NIN,NOUT
C     SAVE

*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
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
*KEEP,RGPRKT.
      DOUBLE PRECISION PRKT1,PRKT2,pktm
	Integer IGAMKT
      COMMON /RGPRKT/PRKT1,PRKT2,pktm,IGAMKT

*KEEP,RGPQCDPOM.
      Integer Iqqg
	COMMON/pqcdpom/ Iqqg
C     SAVE


*KEEP,RGRAHER.
      REAL XPQDIF,XPQPI
	Integer IHERPYS
	Integer NBQ2,NBX
      PARAMETER (NBQ2=20)
      PARAMETER (NBX=20)
      COMMON /RAHER/ IHERPYS,XPQDIF(-6:6,NBX,NBQ2),XPQPI(-6:6,NBX,NBQ2)
*KEEP,RGRGRIDF2.
      DOUBLE PRECISION XX,Q2X
      COMMON /RGRIDF2/ XX(NBX),Q2X(NBQ2)
      REAL F2_DIS,F2_DIF,F2_PI
      COMMON /F2VAL/ F2_DIS(NBX,NBQ2),F2_DIF(NBX,NBQ2),F2_PI(NBX,NBQ2)
cNEW
      DOUBLE PRECISION F2DIS,F2DIF,F2PI
      COMMON /F2INT/ F2DIS,F2DIF,F2PI
cNEW
*KEEP,RGFULL.
      INTEGER IFULL,IQCDGRID
      COMMON /OALPINI/ IFULL,IQCDGRID
*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEEP,RGLUDAT2.
      REAL PMAS,PARF,VCKM
      INTEGER KCHG
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
C      SAVE

*KEND.
      REAL ACC1,ACC2
      Integer IINT,NCB,NDIM,NPOIN,IDEG,NPOINT,ISTART
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB

      COMMON/DIVO/ NDIM,NPOIN
      REAL XMINUS(10),XPLUS(10),ANS,ERROR,PRE
      COMMON/QUADRE/IDEG
      COMMON/SAMPLE/NPOINT
      COMMON/START/ISTART
      INTEGER LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
      REAL PYPAR,PYVAR
      Integer IPY
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      Integer NDIMC
      COMMON /DIMEN/ NDIMC
      Double Precision XD
      DIMENSION XD(20)
      Double Precision XG(20)
      CHARACTER CHIN*100
      CHARACTER * 16 CNAM
      CHARACTER * 7 CINT
      CHARACTER * 4 CRHO
      CHARACTER*8 VERSQQ
      REAL PARA
      Integer MSTA
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      Integer IVM
      COMMON/VMESON/IVM
      REAL XPQ_DIS(-6:6),XPQ_DIF(-6:6),XPQ_PI(-6:6)
      DOUBLE PRECISION POMAX,PIMAX,WTMA
      COMMON/WEIGHT/ POMAX(NBX,NBQ2),PIMAX(NBX,NBQ2),WTMA
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      REAL XPT,Q2T
      Integer LERR,IGENFL
      COMMON/ERR/LERR(100)
      COMMON/GENWEI/IGENFL
      REAL XXI(NBX),Q2I(NBQ2)
      REAL XHF,Q2HF
      Integer KPHF
      COMMON /HEAVYF/XHF,Q2HF,KPHF
      CHARACTER CSTRU*38
      DIMENSION CSTRU(0:10)
      REAL SNGL
      DOUBLE PRECISION C1,Cg
      COMMON/BUCHMUE/C1,Cg
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      Double precision DFUN,DVNOPT,USRINT,FXN1,FXNB
      Integer II,I,J,KPCH,IPROF,KPAO,Int,LU,IFIL,NFLT,IRR,idirf
      Integer idisdiff
      Integer nfltf,ngf,npomf,IDIRO,NGO,NPOMO,IHFO,IVERSQ
      Integer mstp52f,ifullf,IHERPYSO,LUCHGE
      Double Precision DOT,SQRTS,EC,DDUM,qmif,sssf,xff,t2maxf
      Double Precision CTIME
      Integer IT1,IT2
	Double Precision DQ,DX
	COMMON /qxbin/DQ,DX
      EXTERNAL DFUN,DVNOPT,USRINT,FXN1,FXNB
      LOGICAL ex
	
      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF
	
*KEEP,VERSQQ.
      VERSQQ = ' 2.08/06'
      IVERSQ =  20806
*KEND.
      DATA CSTRU/' Simple scaling Function      ',
     +           ' EHLQ set 1                   ',
     +           ' EHLQ set 2                   ',
     +           ' Duke-Owens set 1             ',
     +           ' Duke-Owens set 2             ',
     +           ' Morfin-Tung set 1 (S1)       ',
     +           ' Morfin-Tung set 2 (B1)       ',
     +           ' Morfin-Tung set 3 (B2)       ',
     +           ' Morfin-Tung set 4 (E1)       ',
     +           ' Gluck-Reya-Vogt LO set       ',
     +           ' Gluck-Reya-Vogt HO set       '/
      IGENFL = 0
      DO 10  I=1,100
   10 LERR(I) = 0
      ERROR = 0.0
      CHIN='PARU(11)=0.010;'
      CALL LUGIVE(CHIN)
      MCHARM = PMAS(4,1)
      ALPH=7.299D-3
      ALPHS=0.3D0
      PI=4.D0*DATAN(1.D0)
      IWEI = 0
C...  GIVE ELECTRON FOUR VECTOR
      K(1,1)=21
      K(1,2)=KE
      P(1,1) = 0.0D0
      P(1,2) = 0.0D0
      P(1,3) = DBLE(PLEPIN)
      P(1,5)=DBLE(ULMASS(KE))
      P(1,4)=DSQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
C...  GIVE PROTON FOUR VECTOR
      KP = 2212
      K(2,1)=21
      K(2,2)=KP
      P(2,1)= 0.0D0
      P(2,2)= 0.0D0
      P(2,3)=DBLE(PPIN)
      P(2,5)=DBLE(ULMASS(KP))
      P(2,4)=DSQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
      N=2
      CALL DULIST(1)
C... CALCULATE CMS ENERGY
      DO 20 I=1,4
         CM(I) =P(1,I)+P(2,I)
   20 CONTINUE

C BOOST TO EP CMS
      CALL DUDBRB(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4))
      DO 30 I=1,2
         DO 30 J=1,5
            PBEAM(I,J)=P(I,J)
   30 KBEAM(I,J)=K(I,J)
      KINT(1,2) = 22
      IF(IPRO.EQ.20.OR.IPRO.EQ.21) IDIR = 0
      IF(IDIR.EQ.1) THEN
         KINT(2,2)=2212
      ELSE
         IF(NPOM.EQ.20) THEN
            KINT(2,2) = 211
         ELSEIF(NPOM.EQ.21) THEN
            KINT(2,2) = 111
         ELSE
            KINT(2,2) = 100
         ENDIF
      ENDIF
      IF(INTER.EQ.0) THEN
         KEB = KPH
      ELSEIF(INTER.EQ.2) THEN
         KEB=-24*ISIGN(1,K(1,2))
      ELSE
         WRITE(6,*) ' interaction not implemented: INTER = ',INTER
         STOP
      ENDIF
      SSS = DOT(CM,CM)
      SQRTS=DSQRT(SSS)
      EC= 0.5D0*SQRTS
      MSTA(34)=0
      IF(IPRO.EQ.10.OR.IPRO.EQ.11.OR.IPRO.EQ.12
     +.OR.IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15
     +.OR.IPRO.EQ.18
     +.OR.IPRO.EQ.20.OR.IPRO.EQ.21
     +.OR.IPRO.EQ.100
     +.OR.IPRO.EQ.1200.OR.IPRO.EQ.1400
     +.OR.IPRO.EQ.30.OR.IPRO.EQ.3000
     +.OR.IPRO.EQ.99) THEN
      ELSE
         WRITE(6,*) ' wrong subprocess selected ', IPRO
         WRITE(6,*) ' PROGRAM STOPPED '
         STOP
      ENDIF
	IF(IPRO.EQ.30.OR.IPRO.EQ.3000) IDIR=0
      IF(IPRO.EQ.12.OR.IPRO.EQ.1200.OR.IPRO.EQ.1400) THEN
         IF(IFULL.EQ.0) IRES(1) = 0
      ELSEIF(IPRO.EQ.18) THEN
         IRES(1)=1
         IFULL = 0
      ELSEIF(IPRO.EQ.20) THEN
         IRES(1)=0
         IFULL = 0
      ELSE
         IFULL = 0
      ENDIF
      If (IFULL.EQ.1.OR.IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.20) Then
         MSTA(32) = 0
      EndIf
      IF(NG.GE.30.AND.NG.LT.100) NPOM=NG
      CALL LUNAME(KINT(2,2),CNAM)
      CINT = CNAM
      KPCH = 213
      CALL LUNAME(KPCH,CNAM)
      CRHO=CNAM
      WRITE(6,10000)

      WRITE(6,10100)

      WRITE(6,10200)

      WRITE(6,10300) VERSQQ
      IF(INTER.EQ.0) WRITE(6,10400)
      IF(INTER.EQ.2) WRITE(6,10500)
      IF(IPRO.NE.21) THEN
         IF(NG.EQ.41.OR.NG.EQ.42.AND.IPRO.LT.1000) IPRO=12
      ENDIF
      IF(IPRO.EQ.10.OR.IPRO.EQ.13) THEN
         IF(INTER.EQ.0) WRITE(6,10600) CINT
         IF(INTER.EQ.2) WRITE(6,10700) CINT
      ENDIF
      IF(IPRO.EQ.10)
     +WRITE(6,10800)

      IF(IPRO.EQ.15)
     +WRITE(6,10900)

      IF(IPRO.EQ.13)
     +WRITE(6,11000)
      IF(IPRO.EQ.18) THEN
         IF(IHFLA.GE.4) THEN
            IRPB = 0
            IF(NFLQCDC.LT.4) IRPC = 0
            IRPD = 0
c            IF(NFLQCDC.LT.4) IRPF = 0
            IRPF = 0
            IRPG = 0
         ENDIF
         WRITE(6,11100)
         IF(IRPA.NE.0) THEN
            IF(IHFLA.LT.4) write(6,11200)
            IF(IHFLA.GE.4) write(6,11300)
         ENDIF
         IF(IRPB.NE.0) write(6,11400)
         IF(IRPC.NE.0) THEN
            IF(IHFLA.LT.4) write(6,11500)
            IF(IHFLA.GE.4) write(6,11600)
         ENDIF
         IF(IRPD.NE.0) write(6,11700)
         IF(IRPE.NE.0) THEN
            IF(IHFLA.LT.4) write(6,11800)
            IF(IHFLA.GE.4) write(6,11900)
         ENDIF
         IF(IRPF.NE.0) THEN
            IF(IHFLA.LT.4) write(6,12000)
            IF(IHFLA.GE.4) write(6,12100)
         ENDIF
         IF(IRPG.EQ.1) write(6,12200)
         IF(IRPG.EQ.2) write(6,12300)

         WRITE(6,12400) SCALQ2
c         IF(IQ2.EQ.4) THEN
c            write(6,*) ' DIS resolved photons must not have Q2 as a '
c     +      //'scale'
c            write(6,*) ' PROGRAM STOPPED '
c            STOP
c         ENDIF
c         IF(IDIR.EQ.0) THEN
c            write(6,*) ' DIS resolved photons not possible for '
c     +      //'diffraction'
c            write(6,*) ' PROGRAM STOPPED '
c            STOP
c         ENDIF
      ENDIF
      IF(IPRO.EQ.20) THEN
         IF(IHFLA.GE.4) THEN
            write(6,14400) IHFLA
         ELSE
            write(6,14300) IHFLA
         ENDIF
         IF(Iqqg.eq.0) Then
            write(6,14500)
         elseif(Iqqg.eq.1) Then
            write(6,14600)
         elseif(Iqqg.eq.2) Then
            write(6,14700)
         else
            write(6,*) ' wrong Iqqg selected '
            write(6,*) ' valid values are Iqqg=0,1,2 '
            write(6,*) ' program stopped '
            stop
         endif
         IF(IFPS.EQ.3) THEN
            IFPS = 2
         ELSEIF(IFPS.EQ.1) THEN
            IFPS = 0
         ENDIF
      ELSEIF(IPRO.EQ.21) THEN
         IF(IHFLA.GE.4) THEN
            write(6,14900) IHFLA
         ELSE
            write(6,14800)IHFLA
         ENDIF
         IF(IFPS.EQ.3) THEN
            IFPS = 2
         ELSEIF(IFPS.EQ.1) THEN
            IFPS = 0
         ENDIF
      ENDIF
      IF(IPRO.EQ.30) THEN

         WRITE(6,15200)
         NPOM = 999
         NG = 999
         IFULL = 0
         PT2CUT(30) = 0.
         IF(IFPS.EQ.3) THEN
            IFPS = 2
         ELSEIF(IFPS.EQ.1) THEN
            IFPS = 0
         ENDIF

      ENDIF
      IF(IPRO.EQ.3000) THEN

         WRITE(6,15300)
         NPOM = 999
         NG =999
         IFULL = 0
         PT2CUT(30) = 0.
         IF(IFPS.EQ.3) THEN
            IFPS = 2
         ELSEIF(IFPS.EQ.1) THEN
            IFPS = 0
         ENDIF

      ENDIF
      IF(ISEMIH.EQ.1) THEN
         IF(IPRO.NE.12.AND.IPRO.NE.100) write(6,12500)
         IF(IFULL.EQ.1) write(6,12500)
      ENDIF
      IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
         IF(INTER.EQ.0) THEN
            WRITE(6,12600) CINT
            WRITE(6,12700) MAX(4,IHFLA)
         ENDIF
         IF(INTER.EQ.2) WRITE(6,12800) CINT
         KPA=4
      ENDIF
      IF(IPRO.EQ.100) THEN
         WRITE(6,12900) CRHO

         WRITE(6,13000) CRHO
         IDIR = 0

      ENDIF
      IF(IPRO.EQ.11)
     +WRITE(6,13100)

      IF(IPRO.EQ.14) THEN
         WRITE(6,13200)
      ENDIF

      IF(IPRO.EQ.1200) THEN
         WRITE(6,13300)
         WRITE(6,13500) Q2SUPP

      ENDIF
      IHF = 0
      IF(IPRO.EQ.1400) THEN
         WRITE(6,13400)
         WRITE(6,13500) Q2SUPP
         IHFLA=MAX(4,IHFLA)
         WRITE(6,12700) IHFLA
         IHF = 1
         IPRO = 1200
         IF(IHFLA.GT.NFLAV) NFLAV=IHFLA
      ENDIF
      IF((IPRO.EQ.12.OR.IPRO.EQ.1200).AND.IDIR.EQ.0) THEN
         IF(INTER.EQ.0) WRITE(6,13600) CINT
         IF(INTER.EQ.2) WRITE(6,13700) CINT
         KPA=1
         IF(IVM.GT.0) THEN
            WRITE(6,13800) IVM

            PARJ(11)=1.
            PARJ(12)=1.
            PARJ(13)=1.
         ENDIF
         IF(IFULL.EQ.1.AND.NG.EQ.30.AND.NPOM.EQ.30) IFULL=0
         IF(IFULL.EQ.1.AND.NG.EQ.40.AND.NPOM.EQ.40) IFULL=0
         IF(IFULL.EQ.1.AND.NG.EQ.41.AND.NPOM.EQ.41) IFULL=0
         IF(IFULL.EQ.1.AND.NG.EQ.42.AND.NPOM.EQ.42) IFULL=0
cGB
         IF(IFULL.EQ.1.AND.NG.EQ.13) IFULL=0
cGB
         IF(NG.GT.30.AND.NG.LT.100) THEN
            WRITE(6,13900)
            IF(NG.EQ.30) WRITE(6,14000)
            IF(NG.EQ.40) THEN
               IF(IPRO.LT.1000) THEN
                  IPRO = 12
               ELSE
                  IPRO = 1200
               ENDIF
               WRITE(6,14100)
            ELSEIF(NG.EQ.41) THEN
               IFULL = 0
               IF(IPRO.LT.1000) THEN
                  IPRO = 12
               ELSE
                  IPRO = 1200
               ENDIF
               WRITE(6,14200)
            ELSEIF(NG.EQ.42) THEN
               IFULL = 0
               IF(IPRO.LT.1000) THEN
                  IPRO = 12
               ELSE
                  IPRO = 1200
               ENDIF
               WRITE(6,15000)
            ENDIF
         ENDIF

         IF(NG.EQ.45) THEN
            IFULL = 0
            WRITE(6,15100)
            WRITE(6,15500) C1,Cg
         ENDIF
         IF(IFULL.EQ.1.AND.NG.EQ.40) THEN
            WRITE(6,15600)

         ELSEIF(IFULL.EQ.1) THEN
            WRITE(6,15700)

            IF(IQCDGRID.EQ.0) THEN
               WRITE(6,15800)

            ELSE
               WRITE(6,15900)

            ENDIF
         ENDIF
      ELSEIF((IPRO.EQ.12.OR.IPRO.EQ.1200).AND.IDIR.EQ.1) THEN
         IF(INTER.EQ.0) WRITE(6,16000) CINT
         IF(INTER.EQ.2) WRITE(6,16100) CINT
         KPA=1
         IF(IFULL.EQ.1) THEN
            WRITE(6,16200)

            IF(IQCDGRID.EQ.0) THEN
               WRITE(6,15800)
            ELSE
               WRITE(6,15900)
            ENDIF
         ENDIF
      ENDIF
      IF(NG.EQ.45) THEN
         IFULL = 0
         WRITE(6,15100)
         WRITE(6,15500) C1,Cg
      ENDIF
      IF(IPRO.EQ.100) THEN

         WRITE(6,*)
     +    '*  EPA + gamma -->',CRHO,'  + pomeron                   *'
         IF(NPOM.GE.20) WRITE(6,*) ' Wrong NPOM selected. NPOM = ',NPOM
      ENDIF

      IF(IPRO.EQ.99) THEN

         WRITE(6,15400)
      ENDIF


      IF(THEMA.NE.180.0D0) THEN
         WRITE(6,16300) THEMA

      ELSE
         WRITE(6,*) '*  no cut on max angle of scattered electron      '
     +   //'  *'
      ENDIF
      IF(THEMI.NE.0.0D0) THEN
         WRITE(6,16400) THEMI

      ELSE
         WRITE(6,*) '*  no cut on min angle of scattered electron      '
     +   //'  *'
      ENDIF
      IF(QMI.GT.0.0) THEN
         WRITE(6,16500) QMI

      ELSE
         WRITE(6,*) '*  Q^2 _min according to kinematics               '
     +   //'  *'
      ENDIF
      IF(QMA.LT.10D8) THEN
         WRITE(6,16600) QMA

      ELSE
         WRITE(6,*) '*  Q^2 _max according to kinematics               '
     +   //'  *'
      ENDIF
      IF(YMI.GT.0.0) THEN
         WRITE(6,16700) YMI

      ELSE
         WRITE(6,*) '*  y_min according to kinematics                  '
     +   //'  *'
      ENDIF
      IF(YMA.LT.1.0) THEN
         WRITE(6,16800) YMA

      ELSE
         WRITE(6,*) '*  y_max according to kinematics                  '
     +   //'  *'
      ENDIF
      IF(IPRO.EQ.20.OR.IPRO.EQ.21.OR.IPRO.EQ.30.OR.IPRO.EQ.3000) THEN
      ELSE
         WRITE(6,16900) NFLAV
         WRITE(6,17000) NFLQCDC
      ENDIF
      IPY(8) = NFLAV
      IPROF = IPRO
      IF(IPRO.GT.1000) IPROF = IPRO/100
      IF(IPROF.EQ.12) THEN
         IF(NG.GE.41.AND.NG.LT.50) PT2CUT(IPROF) = PT2CUT(13)
      ENDIF
      IF(PT2CUT(IPROF).GT.0.0.AND.IPRO.NE.100) THEN
         WRITE(6,17100) PT2CUT(IPROF),IPRO

      ELSE
         WRITE(6,17200) IPRO

         IF(IFULL.EQ.1) THEN
            WRITE(6,17300) PT2CUT(15),PT2CUT(13)

         ENDIF
      ENDIF
cc      write(6,*) IDIR,IDISDIF
      IF(IDIR.EQ.0) IDISDIF=0
      IF(IPRO.EQ.12.OR.IPRO.EQ.1200) THEN
      ELSE
         WRITE(6,*) '* IDISDIF=0 forced, since IPRO.ne.12           '
         IDISDIF = 0
      ENDIF
      IF(NG.GE.30.and.IDISDIF.GT.0) THEN
         WRITE(6,*) '* IDISDIF=0 forced, since NG>30                '
         IDISDIF = 0
      ENDIF
      IF(IPRO.EQ.18) THEN
	   IF(IGAMKT.EQ.1) THEN
         write(6,17400) prkt1
	   ELSEIF(IGAMKT.EQ.2) THEN
         write(6,17401) prkt1,pktm
	   ELSEIF(IGAMKT.EQ.3) THEN
         write(6,17402) prkt1,pktm
	   ELSE
	   write(6,*) ' no valid distr. for prim. kt. ',
     &  'in photon selected: IGAMKT=',IGAMKT
	
	   ENDIF
      endif
      IF(IDIR.EQ.1) THEN
         write(6,17500) prkt2
      endif

      IF(IPRO.EQ.20.OR.IPRO.EQ.21) THEN
         WRITE(6,*) '* cuts for diffractive events:               '
     +   //'       *'
         WRITE(6,18000) T2MAX
         WRITE(6,18200) XF
      ELSE
         IF(IDIR.EQ.0.OR.IDISDIF.GE.1) THEN
            IF(IPRO.NE.100) THEN
               WRITE(6,*) '***         diffractive hard scattering     '
     +         //'      ***'
               IF(NG.EQ.0) THEN
                  WRITE(6,*) '* parton density in pomeron xf_0(x)=6x(1-'
     +            //'x)         *'
               ELSEIF(NG.GT.0.AND.NG.LE.5) THEN
                  WRITE(6,18400) NG,NG+1,NG
               ELSEIF(NG.EQ.10) THEN
                  WRITE(6,*) '*  parton density', 'in pomeron xf_10(x)='
     +            //'           **   (0.18+5.46x)(1-x)            '
     +            //'                   *'
               ELSEIF(NG.EQ.11) THEN
                  WRITE(6,*) '* parton density in pomeron              '
     +            //'           *'
                  WRITE(6,*) ' acc. to Donnchie Landshoff '
               ELSEIF(NG.EQ.12) THEN
                  WRITE(6,*) '* parton density in pomeron              '
     +            //'           *'
                  WRITE(6,*) '* acc. to Kniehl, Kohrs, Kramer DESY 94-'
     +            //'140         *'
                  WRITE(6,*) '* with pointlike contribution            '
     +            //'           *'
cGB
               ELSEIF(NG.EQ.13) THEN
                  Write(6,*)' '
                  WRITE(6,*) '* ==> Parton density in pomeron xf(x)=   '
     +            //'           **  (3/4){x(1-x)+(0.57/2.)(1-x)^2}'
     +            //'                   *'
                  Write(6,*)' '
cGB
               ELSEIF(NG.LT.0) THEN
                  IF(NG.ge.-15.and.NG.le.-10) Then
                  ELSEIF(NG.ge.-5.and.NG.le.-1) Then
			ELSEIF(NG.EQ.-20) THEN
                  ELSE
                     WRITE(6,*) '* user has HOPEFULLY provided gluon '
     +               //'density         *'
                     WRITE(6,*) '* via FUNCTION USDIFFR                '
     +               //'              *'
                  ENDIF
               ELSEIF(NG.EQ.30) THEN
                  WRITE(6,17900)

               ENDIF
               IF(NPOM.LT.20.OR.NPOM.GE.40) THEN

                  IF(NPOM.LT.20.AND.NPOM.GE.0) THEN
                     WRITE(6,18100) RN2,ALPHP,XF,EPSP
                  ENDIF
                  WRITE(6,*) '* cuts for diffractive events:           '
     +            //'           *'
                  WRITE(6,18000) T2MAX
                  WRITE(6,18200) XF

               ELSEIF(NPOM.EQ.20.OR.NPOM.EQ.21) THEN
                  WRITE(6,*) '* parameters for pion distribution:      '
     +            //'           *'
                  WRITE(6,18000) T2MAX
                  IF(XF.GT.0.8) THEN
                     WRITE(6,*) '* selected x_f too small for pion '
     +               //'exchange          *'
                     WRITE(6,*) '* make sure that x_f is what you want '
     +               //'              *'
                  ENDIF
                  WRITE(6,18200) XF

               ENDIF
               IF(NPOM.EQ.0) THEN
                  WRITE(6,*) '* pomeron distribution a la Streng is '
     +            //'used          *'
               ELSEIF(NPOM.EQ.1) THEN
                  WRITE(6,*) '* pomeron distribution a la Ingelman is '
     +            //'used        *'
               ELSEIF(NPOM.LT.0) THEN
                  IF(NPOM.ge.-12.and.NPOM.le.-10) Then
                  ELSEIF(NPOM.ge.-3.and.NPOM.le.-1) Then
                  ELSEIF(NPOM.eq.-20) Then
                  ELSE
                     WRITE(6,*) '* user has HOPEFULLY provided pomeron '
     +               //'distributions *'
                     WRITE(6,*) '* via FUNCTION USDIFFR                '
     +               //'              *'
                  ENDIF
               ELSEIF(NPOM.EQ.30) THEN
                  WRITE(6,18300)

               ENDIF
            ENDIF
         ENDIF
      ENDIF
      WRITE(6,*) '*****************************************************'
      IF(IFPS.NE.0.AND.IPRO.EQ.100) THEN
         WRITE(6,*) '##   no parton shower possible for IPRO = 100   ##'
         IFPS = 0
      ENDIF
      IF(IFPS.NE.0.AND.IPRO.EQ.12.AND.NG.EQ.30.AND.NPOM.EQ.30) THEN
         IF(IFPS.EQ.1.OR.IFPS.EQ.3) THEN
            WRITE(6,17600)

         ENDIF
         IF(IFPS.EQ.3) IFPS = 2
         IF(IFPS.EQ.1) IFPS = 0
      ENDIF
      IF(IFPS.EQ.0) THEN
         IPY(13) = 0
         IPY(14) = 0
      ENDIF
      IF(IFPS.NE.0) THEN
         WRITE(6,*) '##################################################'
     +   //'###'
         WRITE(6,*) '#     parton shower selection:                    '
     +   //'  #'
         IF(IFPS.EQ.0) THEN
            IPY(13) = 0
            IPY(14) = 0
         ELSEIF(IFPS.EQ.1) THEN
            IPY(13) = 0
         ELSEIF(IFPS.EQ.2) THEN
            IPY(14) = 0
         ELSEIF(IFPS.EQ.10) THEN
            IPY(13) = 0
            IPY(14) = 0
         ENDIF
         IF(IPY(14).GE.1) THEN
c            IF(ILEPTO.EQ.1) WRITE(6,*) '#     initial state parton '
c     +      //'shower a la LEPTO        #'
c            IF(ILEPTO.EQ.0) WRITE(6,*) '#     initial state parton '
c     +      //'shower a la PYTHIA       #'
            write(6,*) '#     initial state parton '
     +      //'parameters:              #'
            IF(IORDER.EQ.0) write(6,*) '#     no ordering of emitted '
     +      //'partons                #'
            IF(IORDER.EQ.1) write(6,*) '#     ordering in Q2 of '
     +      //'emitted partons             #'
            IF(IORDER.EQ.2) write(6,*) '#     ordering in Q2 and angle'
     +      //' of emitted partons   #'
            IF(IALPS.EQ.1) write(6,*) '#     alphas first order with s'
     +      //'cale Q2              #'
            IF(IALPS.EQ.2) write(6,*) '#     alphas first order with s'
     +      //'cale k_t**2=(1-z)*Q2 #'
            IF(ITIMSHR.EQ.0) write(6,*) '#     no shower of time like '
     +      //'partons                #'
            IF(ITIMSHR.EQ.1) write(6,*) '#     time like partons may '
     +      //'shower                  #'
            IF(ISOFTG.EQ.0) write(6,*) '#     soft gluons are entire'
     +      //'ly neglected          #'
            IF(ISOFTG.EQ.1) write(6,*) '#     soft gluons are resummed'
     +      //'                      #'
         ENDIF
         IF(IPY(13).EQ.1) WRITE(6,*) '#     final state parton shower  '
     +   //'                   #'
         IF(IFPS.EQ.10) Then
            WRITE(6,*) '#     parton shower with ARIADNE          '
     +      //'          #'
            write(6,*) '# BGF treatment (0/1=BGF from RAPGAP/ARIADNE)',
     +      MSTA(32)
            write(6,*) '# pomeron treatment (0/1=pom from RAPGAP/'
     +      //'ARIADNE)', MSTA(34)
            write(6,*) '# treatment of struck quark and remnant   '
     +      //'          #'
            write(6,*) '#               MSTA(30)= ',MSTA(30)
            write(6,*) '#               PARA(11)= ',PARA(11)
         ENDIF
         WRITE(6,*) '##################################################'
     +   //'###'
      ENDIF
      IF(IPRO.EQ.10.OR.
     +   IPRO.EQ.11.OR.IPRO.EQ.13.OR.IPRO.EQ.12.OR.
     +   IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.100) NFT = 2
      IF(IPRO.EQ.20) NFT=3
      IF(IPRO.EQ.21) NFT=2
c      IRES(1)=0
      IRES(2)=0
      WRITE(6,17700) SQRTS
      IF(IPRO.GE.1000) THEN
         IF(IRUNAEM.NE.0)
     +    write(6,*) ' alpha_em forced to be fixed for HERACLES'
         IRUNAEM=0
      ENDIF
      IF(IRUNAEM.EQ.0) THEN
         WRITE(6,*) ' alpha_em fixed; alpha_em = ',ALPH
      ELSE
         WRITE(6,*) ' running alpha_em selected '
      ENDIF
      IF(IRUNA.EQ.0) THEN
         MSTU(111) = 0
         PARU(111) = SNGL(ALPHS)
         WRITE(6,*) ' alpha_s fixed; alpha_s = ',ALPHS
      ENDIF
      IF((NG.EQ.41.OR.NG.EQ.42).AND.IDIR.EQ.0) THEN
         WRITE(6,*) ' scale for alpha_s: (p_t**2+m_q**2)/(1-beta) '
      ELSEIF(IPRO.EQ.21.AND.IDIR.EQ.0) THEN
         WRITE(6,*) ' scale for alpha_s: (p_t**2+m_q**2)/(1-beta) '
      ELSEIF(IPRO.EQ.20.AND.IDIR.EQ.0) THEN
         WRITE(6,*) ' fixed alpha_s=0.25 '
      ELSE
         IF(IQ2.EQ.1) WRITE(6,*) ' scale for alpha_s: 4*m_q**2 '
         IF(IQ2.EQ.2) WRITE(6,*) ' scale for alpha_s: shat '
         IF(IQ2.EQ.3) WRITE(6,*) ' scale for alpha_s: 4*m_q**2 + p_t **'
     +   //'2 '
         IF(IQ2.EQ.4) WRITE(6,*) ' scale for alpha_s: Q2 '
         IF(IQ2.EQ.5) WRITE(6,*) ' scale for alpha_s: Q2 + p_t **2'
         IF(IQ2.EQ.6) WRITE(6,*) ' scale for alpha_s: kt**2 (Breit fr.)'
         IF(IQ2.EQ.7) WRITE(6,*) ' scale for alpha_s: pt**2+8*pt**2/s'
         IF(IQ2.LT.1.OR.IQ2.GT.7) THEN
            WRITE(6,*) ' no valid scale selected. PROGRAM STOPS '
            STOP
         ENDIF
      ENDIF
      WRITE(6,*) ' scale is multiplied by: ',SCALFA
      WRITE(6,17800) PARU(112),MSTU(112)
      PYPAR(21) = PARU(112)

      IF(IDIR.EQ.1.OR.IDISDIF.GE.1) THEN
         WRITE(6,*) ' *** Deep Inelastic Scattering      ***'
         IF(MSTP(51).GT.10) THEN
            WRITE(6,*) ' p structure function from PDFLIB ',MSTP(51)
         ELSE
            WRITE(6,*) ' p structure function',MSTP(51),CSTRU(MSTP(51))
         ENDIF
      ENDIF


      KPAO = KPA

      IF(IPRO.LT.1000) THEN

         IF(IPRO.EQ.10.OR.IPRO.EQ.11.OR.IPRO.EQ.15) NDIM = 7
         IF(IPRO.EQ.13.OR.IPRO.EQ.14) NDIM=8
         IF(IPRO.EQ.12) NDIM=4
         IF(IPRO.EQ.100) NDIM = 4
         IF(IPRO.EQ.30) NDIM = 4
         NDIMC = 9999
   40    CALL draprnV(XG,20)
         DO 50 I=1,20
   50    XD(I)=XG(I)
         DDUM=DFUN(20,XD)
*         WRITE(6,*) ' NDIM dynamically calculated ',NDIM,NDIMC
         IF(NDIMC.GT.20) goto 40
         NDIM = NDIMC
         INT=0
         NIN=0
         NOUT=0

         IF(IINT.EQ.1) THEN
            IF(NDIM.GT.10) THEN
               write(6,*) ' more than 10 dimensions requested '
               write(6,*) ' DIVON can handle only 10 '
               write(6,*) ' user must switch to BASES integration '
               write(6,*) ' program stopped '
               stop
            endif
            DO 60 I=1,10
               XMINUS(I)=1.E-11
   60       XPLUS(I)=1.E0-XMINUS(I)
            PRE =1.75
c         PRE =0.75
            IDEG=1
            NPOINT=200
            NPOIN= NPOINT
            WRITE(6,*) ' precision parameter in DIVON = ',PRE
            WRITE(6,*) ' npoint = ',NPOINT,IDEG
            CALL PARTN(NDIM,XMINUS,XPLUS,PRE,750000)
            if(NDIM.GT.6) then
               istart = 2
               pre = 0.75
               CALL PARTN(NDIM,XMINUS,XPLUS,PRE,750000)
            endif
            CALL INTGRL(NDIM,0,400,ANS,ERROR)
c         CALL INTGRL(NDIM,1,200,ANS,ERROR)
         ELSE
********************************************************************
*     Initialization of BASES/SPRING V5.1
********************************************************************
*===> Initialization of BASES by calling BSINIT

            CALL BSINIT
            NDIMB=NDIM
            NWILDB = NDIM
            NCALLB = NCB
            ITMX1B = 30
            ACC1B = DBLE(ACC1)
            ITMX2B = 100
            ACC2B = DBLE(ACC2)
            DO 70 I=1,NDIMB
               IG(I) = 1
               XLB(I)=1.0D-12
   70       XUB(I)=1.D0
********************************************************************
*     Nimerical Integration by BASES V5.1
********************************************************************

            CALL BASES( FXNB, AVGI, SD, CTIME, IT1, IT2 )

            LU = 6
            ANS = SNGL(AVGI)
            ERROR = SNGL(SD)
         ENDIF
         SD=DBLE(ERROR)
         AVGI=DBLE(ANS)

c check for calculating QCDweights in grid
         ifil = 0
         NFLT = 6
         YMIN = DMAX1(YMI,0.00001D0)
         YMAX = DMIN1(YMA,1.D0)
         Q2MIN = DMAX1(QMI,0.5D0)
         Q2MAX = DMIN1(QMA,SSS)
         IF((IFULL.EQ.1.AND.IDIR.EQ.0).OR.(IDISDIF.GE.1)) THEN
            write(6,*) ' calculating xpq_dif in grid IFULL/IDIR/IDISDIF'
            write(6,*) ' check existence of file with pdf"s '
            inquire(FILE='pdf.dat',EXIST=ex)
            if(ex) then
               open(31,FILE='pdf.dat', FORM='unformatted',STATUS= 'OLD',
     +         IOSTAT=IRR,ERR=220)
               write(6,*) ' file 31 ',irr
               read(31) idirf,idisdiff,nfltf,ngf,npomf,mstp52f
               read(31) ifullf,qmif,sssf,xff,t2maxf
c               write(6,*) idirf,idisdiff,nfltf
c               write(6,*) ifullf,qmif,sssf
c               write(6,*) idir,idisdif,nflt
c               write(6,*) ifull,qmi,sss
               if(idirf.eq.idir.and.idisdiff.eq.idisdif.and.nfltf.eq.nf
     +         lt .and.ifullf.eq.ifull.and.ngf.eq.ng.and.npomf.eq.npom
     +         .and.mstp52f.eq.mstp(52) .and.qmif.eq.qmi.and.sssf.eq.ss
     +         s .and.xff.eq.xf.and.t2maxf.eq.t2max) then
c               if(idirf.eq.idir.and.idisdiff.eq.idisdif.and.nfltf.eq.nf
c     +         lt .and.ifullf.eq.ifull) then
                  write(6,*) ' file found with consistent parameters '
                  ifil = 1

               else
c                  write(6,*) idirf,idisdiff,nfltf
c                  write(6,*) ifullf,qmif,sssf
c                  write(6,*) idir,idisdif,nflt
c                  write(6,*) ifull,qmi,sss
                  write(6,*) ' inconsistent parameters on file on unit '
     +            //'31'
                  write(6,*) ' cold start necessary '
                  close(31)
               endif
            endif
            If(ifil.eq.0) then
               open(31,FILE='pdf.dat', FORM='unformatted',STATUS=
     +         'NEW',IOSTAT=IRR,ERR=210)
               write(6,*) ' no file found on unit 31 '
               write(6,*) ' cold start necessary '
               write(31) idir,idisdif,nflt,ng,npom,MSTP(52)
               write(31) ifull,qmi,sss,xf,t2max
            endif
            IDIRO = IDIR
            IDIR = 0
            DQ = (LOG(SSS) - LOG(QMI))/DBLE(NBQ2-1)
            DX = (LOG10(1.D0) - LOG10(QMI/SSS))/DBLE(NBX-1)
            DO 90  I=1,NBX
               POM = 0.D0
               PIM = 0.D0
               DO 90  J=1,NBQ2
                  XXI(I) = SNGL(10**(LOG10(QMI/SSS) + DX*DFLOAT(I-1)))
                  IF(XXI(I).GE.1.D0) XXI(I) = 0.999
                  XPT = XXI(I)
                  Q2I(J) = SNGL(EXP(LOG(QMI) + DQ*DFLOAT(J-1)))
                  Q2T = Q2I(J)
                  XX(I) = DBLE(XXI(I))
                  Q2X(J) = DBLE(Q2I(J))
                  XPR = XPT
                  Q2 = DBLE(Q2T)
                  If(ifil.eq.0) then
                     CALL PYSTFU(2212,XPT,Q2T,XPQ_DIF)
                     write(31) pom
                  else
                     read(31) pom
                  endif
                  IF(POM.GT.POMAX(I,J)) POMAX(I,J) = POM
                  DO 80  II = -NFLT,NFLT
                     if(ifil.eq.0) then
                        write(31) XPQ_DIF(II)
                     else
                        read(31) XPQ_DIF(II)
                     endif
                     XPQDIF(II,I,J) = XPQ_DIF(II)
                     XPQPI(II,I,J) = XPQ_DIF(II)
cc                  write(6,*)'rapgap x,q2,pom', xx(i),q2x(j),pom
cc                   write(6,*)'rapgap XOQ ',XPQ_DIF(II)
   80             CONTINUE
   90       CONTINUE
            write(6,*) ' xpq for DIF done '
            IF(IDISDIF.EQ.2) THEN
               write(6,*) ' calculating xpq_pi in grid '
               NPOMO = NPOM
               NGO = NG
               IDIR = 0
               NPOM = 20
               NG = 20
               DO 110 I=1,NBX
                  PIM = 0.D0
                  DO 110 J=1,NBQ2
                     XXI(I) = SNGL(10**(LOG10(QMI/SSS) + DX*DFLOAT(I-1)
     +               ))
                     IF(XXI(I).GE.1.D0) XXI(I) = 0.999
                     XPT = XXI(I)
                     Q2I(J) = SNGL(EXP(LOG(QMI) + DQ*DFLOAT(J-1)))
                     Q2T = Q2I(J)
                     XX(I) = DBLE(XXI(I))
                     Q2X(J) = DBLE(Q2I(J))
                     XPR = XPT
                     Q2 = DBLE(Q2T)
                     If(ifil.eq.0) then
                        CALL PYSTFU(2212,XPT,Q2T,XPQ_PI)
                        write(31) pim
                     else
                        read(31) pim
                     endif
                     IF(PIM.GT.PIMAX(I,J)) PIMAX(I,J) = PIM
                     DO 100 II = -NFLT,NFLT
                        if(ifil.eq.0) then
                           write(31) XPQ_PI(II)
                        else
                           read(31) XPQ_PI(II)
                        endif
                        XPQPI(II,I,J) = XPQ_PI(II)
c                  write(6,*)'rapgap x,q2,pom', xx(i),q2x(j),pom
  100                CONTINUE
  110          CONTINUE
               write(6,*) ' xpq for PI done '
               NG = NGO
               NPOM = NPOMO
            ENDIF
            IHERPYS = 1
c            write(6,*) ' RAPGAP: IHERPYS ',IHERPYS
            IDIR = IDIRO
         ENDIF
      ELSEIF(IPRO.EQ.1200.OR.IPRO.EQ.3000) THEN
c      write(6,*) IDIR,IDISDIF,IFULL
c      write(6,*) NG,NPOM,MSTP(51)
         YMIN = DMAX1(YMI,0.00001D0)
         YMAX = DMIN1(YMA,1.D0)
         Q2MIN = DMAX1(QMI,0.5D0)
         Q2MAX = DMIN1(QMA,SSS)
         IHERPYS = 0
         ifil = 0
         NFLT = 6
         IHFO = IHF
         IHF = 0
         IF(IDIR.EQ.0.OR.IDISDIF.GE.1) THEN
            write(6,*) ' calculating xpq_dif in grid HERACLES'
            IDIRO = IDIR
            NGO = NG
            NPOMO = NPOM
            write(6,*) ' check existence of file with pdf"s '
            inquire(FILE='pdf.dat',EXIST=ex)
            if(ex) then
               open(31,FILE='pdf.dat', FORM='unformatted',STATUS= 'OLD',
     +         IOSTAT=IRR,ERR=220)
               write(6,*) ' file 31 ',irr
               read(31) idirf,idisdiff,nfltf,ngf,npomf,mstp52f
               read(31) ifullf,qmif,sssf,xff,t2maxf
               write(6,*) idirf,idisdiff,nfltf,mstp52f
               write(6,*) ifullf,qmif,sssf,xff,t2maxf
               write(6,*) idir,idisdif,nflt,mstp(52)
               write(6,*) ifull,qmi,sss,xf,t2max
               if(idirf.eq.idir.and.idisdiff.eq.idisdif.and.nfltf.eq.nf
     +         lt .and.ifullf.eq.ifull.and.ngf.eq.ng.and.npomf.eq.npom
     +         .and.mstp52f.eq.mstp(52) .and.qmif.eq.qmi.and.sssf.eq.ss
     +         s .and.xff.eq.xf.and.t2maxf.eq.t2max) then
c               if(idirf.eq.idir.and.idisdiff.eq.idisdif.and.nfltf.eq.nf
c     +         lt .and.ifullf.eq.ifull) then
                  write(6,*) ' file found with consistent parameters '
                  ifil = 1

               else
c                  write(6,*) idirf,idisdiff,nfltf
c                  write(6,*) ifullf,qmif,sssf
c                  write(6,*) idir,idisdif,nflt
c                  write(6,*) ifull,qmi,sss
                  write(6,*) ' inconsistent parameters on file on unit '
     +            //'31'
                  write(6,*) ' cold start necessary '
                  close(31)
               endif
            endif
            If(ifil.eq.0) then
               open(31,FILE='pdf.dat', FORM='unformatted',STATUS=
     +         'NEW',IOSTAT=IRR,ERR=210)
               write(6,*) ' no file found on unit 31 '
               write(6,*) ' cold start necessary '
               write(31) idir,idisdif,nflt,ng,npom,MSTP(52)
               write(31) ifull,qmi,sss,xf,t2max
            endif
            IDIR = 0
            DQ = (LOG(SSS) - LOG(QMI))/DBLE(NBQ2-1)
            DX = (LOG10(1.D0) - LOG10(QMI/SSS))/DBLE(NBX-1)
            DO 130 I=1,NBX
               DO 130 J=1,NBQ2
                  XXI(I) = SNGL(10**(LOG10(QMI/SSS) + DX*DFLOAT(I-1)))
                  IF(XXI(I).GE.1.D0) XXI(I) = 0.999
                  XPT = XXI(I)
c                  if(xpt.lt.1.e-3) goto 120
                  Q2I(J) = SNGL(EXP(LOG(QMI) + DQ*DFLOAT(J-1)))
                  Q2T = Q2I(J)
                  XX(I) = DBLE(XXI(I))
                  Q2X(J) = DBLE(Q2I(J))
                  XPR = XPT
                  Q2 = DBLE(Q2T)
                  POM = 0.D0
                  PIM = 0.D0
c                  write(6,*) ' RAPGAP integrate pom/pi pdf: X,Q2',xpt,q2t
                  If(ifil.eq.0) then
                     CALL PYSTFU(2212,XPT,Q2T,XPQ_DIF)
c              write(6,*) 'XPT,Q2T,XPQ_DIF',XPT,Q2T,XPQ_DIF
                     IF(NG.EQ.20.or.NG.eq.21) THEN
                     write(31) pim
			   ELSE
                     write(31) pom
			   ENDIF
                  else
                     IF(NG.EQ.20.or.NG.eq.21) THEN
                     read(31) pim
			   ELSE
                     read(31) pom
			   ENDIF
                  endif
                     IF(NG.EQ.20.or.NG.eq.21) THEN
                  IF(PIM.GT.PIMAX(I,J)) PIMAX(I,J) = PIM
			   ELSE
                  IF(POM.GT.POMAX(I,J)) POMAX(I,J) = POM
                     ENDIF

                  DO 120 II = -NFLT,NFLT
                     if(ifil.eq.0) then
                        write(31) XPQ_DIF(II)
                     else
                        read(31) XPQ_DIF(II)
                     endif
                     XPQDIF(II,I,J) = XPQ_DIF(II)
                     XPQPI(II,I,J) = XPQ_DIF(II)
c                  write(6,*) xpq_dif(II),xx(i),q2x(j)
  120             CONTINUE
  130       CONTINUE
            write(6,*) ' xpq  done '
            IF(IDISDIF.EQ.2) THEN
               NPOMO = NPOM
               NGO = NG
               IDIR = 0
               NPOM = 20
               NG = 20
               PIM = 0.D0
               DO 150 I=1,NBX
                  DO 150 J=1,NBQ2
                     XXI(I) = SNGL(10**(LOG10(QMI/SSS) + DX*DFLOAT(I-1)
     +               ))
                     IF(XXI(I).GE.1.D0) XXI(I) = 0.999
                     XPT = XXI(I)
                     Q2I(J) = SNGL(EXP(LOG(QMI) + DQ*DFLOAT(J-1)))
                     Q2T = Q2I(J)
                     XX(I) = DBLE(XXI(I))
                     Q2X(J) = DBLE(Q2I(J))
                     XPR = XPT
                     Q2 = DBLE(Q2T)
                     POM = 0.D0
                     PIM = 0.D0
                     If(ifil.eq.0) then
                        CALL PYSTFU(2212,XPT,Q2T,XPQ_PI)
                        write(31) pim
                     else
                        read(31) pim
                     endif
                     IF(PIM.GT.PIMAX(I,J)) PIMAX(I,J) = PIM
			
                     DO 140 II = -NFLT,NFLT
                        if(ifil.eq.0) then
                           write(31) XPQ_PI(II)
                        else
                           read(31) XPQ_PI(II)
                        endif
                        XPQPI(II,I,J) = XPQ_PI(II)
c                  write(6,*) xpq_dif(II),xx(i),q2x(j)
  140                CONTINUE
  150          CONTINUE
               write(6,*) ' xpq for pi done '
            ENDIF
            IDIR = IDIRO
            NG = NGO
            NPOM = NPOMO
            IF(IDIR.EQ.0) IHERPYS = 1
         ENDIF
         IHF = IHFO
         CALL HERACL(AVGI,SD,1)
      ELSE
         write(6,*) 'process IPRO = ',IPRO,' not implemented  '
      ENDIF

      write(6,*) IDIR,IDISDIF,IFULL
      write(6,*) NG,NPOM,MSTP(51)
c check whether mixing of DIS, DIF, and pi exchange is wanted
      IF(IDISDIF.GE.1.AND.IDIR.EQ.1) THEN
         DQ = (LOG(SSS) - LOG(QMI))/DBLE(NBQ2-1)
         DX = (LOG10(1.D0) - LOG10(QMI/SSS))/DBLE(NBX-1)
         IHERPYSO = IHERPYS
         DO 190 I=1,NBX
            POM = 0.D0
            DO 190 J=1,NBQ2
               XXI(I) = SNGL(10**(LOG10(QMI/SSS) + DX*DFLOAT(I-1)))
               IF(XXI(I).GE.1.D0) XXI(I) = 0.999
               XPT = XXI(I)
               Q2I(J) = SNGL(EXP(LOG(QMI) + DQ*DFLOAT(J-1)))
               Q2T = Q2I(J)
               XX(I) = DBLE(XXI(I))
               Q2X(J) = DBLE(Q2I(J))
               XPR = XPT
               NPOMO = NPOM
               IDIRO = IDIR
               NGO = NG
               IDIR = 1
               IHERPYS = 0
               CALL PYSTFU(2212,XPT,Q2T,XPQ_DIS)
               F2_DIS(I,J) = 0.
               DO 160 II=-6,6
c               write(6,*) 'RAPGAP',XPQ_DIS(II),LUCHGE(II),II
  160          F2_DIS(I,J) = F2_DIS(I,J) + FLOAT(LUCHGE(II))**2 /9.0 *
     +         XPQ_DIS(II)
c               write(6,*) ' F2_DIS ',F2_DIS(I,J),XX(I),Q2X(J)
c               write(6,*) iherpys,iherpyso
               IHERPYS = IHERPYSO
               IHERPYS = 1
               IDIR = 0
               NPOM = NPOMO
               POM = 0.D0
               CALL PYSTFU(2212,XPT,Q2T,XPQ_DIF)
               IF(POM.GT.POMAX(I,J)) POMAX(I,J)=POM
               F2_DIF(I,J) = 0.
               DO 170 II=-6,6
c                  write(6,*) 'RAPGAP',XPQ_DIF(II),LUCHGE(II),II
  170          F2_DIF(I,J) = F2_DIF(I,J) + FLOAT(LUCHGE(II))**2 /9. *
     +         XPQ_DIF(II)
c               write(6,*) ' F2_DIF ',F2_DIF(I,J),XX(I),Q2X(J)
               F2_PI(I,J) = 0.
               IF(IDISDIF.EQ.2) THEN
                  IDIR = 0
                  NPOM = 20
                  NG = 20
                  PIM = 0.D0
                  CALL PYSTFU(2212,XPT,Q2T,XPQ_PI)
                  IF(PIM.GT.PIMAX(I,J)) PIMAX(I,J)=PIM
                  DO 180 II=-6,6
c               write(6,*) 'RAPGAP',XPQ_PI(II),LUCHGE(II),II
  180             F2_PI(I,J) = F2_PI(I,J) + FLOAT(LUCHGE(II))**2 /9. *
     +            XPQ_PI(II)
c               write(6,*) ' F2_PI ',F2_PI(I,J),XX(I),Q2X(J)
               ENDIF
               IDIR =IDIRO
               NPOM =NPOMO
               NG = NGO
  190    CONTINUE
         IHERPYS = IHERPYSO
         WRITE(6,18500)

         DO 200 J=1,NBQ2
            DO 200 I=1,NBX
               WRITE(6,18600) XX(I),Q2X(J),F2_DIS(I,J),F2_DIF(I,J),
     +         F2_PI(I,J)

  200    CONTINUE
      ENDIF
      IF(IFULL.EQ.1.AND.IQCDGRID.EQ.1) THEN
         CALL QCDGRID
      ENDIF

      KPA = KPAO
      RETURN
  210 write(6,*) ' error opening new file for pdf"s '
      stop
  220 write(6,*) ' error opening old file for pdf"s '
      STOP
10000 FORMAT(' *****************************************************')
10100 FORMAT(' *                                                   *')
10200 FORMAT(' *     You are using the RAPGAP MC generator         *')
10300 FORMAT(' *           version  ',A8,'                       *')
10400 FORMAT(' *    neutral current interaction selected           *')
10500 FORMAT(' *    charged current interaction selected           *')
10600 FORMAT(' *  gamma + gluon_',A6,' --> q q_bar(light quarks)   *')
10700 FORMAT(' *      W + gluon_',A6,' --> q q_bar(light quarks)   *')
10800 FORMAT(' *  EPA + gamma* gluon --> q q_bar (massless) used   *')
10900 FORMAT(' *  full ME    e q --> q gluon (massless) used       *')
11000 FORMAT(' *  full matrixelement e gluon --> e" q q_bar        *')
11100 FORMAT(' *  resolved photon processes                        *')
11200 FORMAT(' *               g + g     --> q + q_bar             *')
11300 FORMAT(' *               g + g     --> Q + Q_bar             *')
11400 FORMAT(' *               g + g     --> g + g                 *')
11500 FORMAT(' *               g + q     --> g + q                 *')
11600 FORMAT(' *               g + Q     --> g + Q                 *')
11700 FORMAT(' *               q + q_bar --> g + g                 *')
11800 FORMAT(' *               q + q_bar --> q + q_bar             *')
11900 FORMAT(' *               q + q_bar --> Q + Q_bar             *')
12000 FORMAT(' *               q + q     --> q + q                 *')
12100 FORMAT(' *               q + q     --> Q + Q                 *')
12200 FORMAT(' *  Mueller/Tang q + q     --> q + q                 *')
12300 FORMAT(' *  massive glu  q + q     --> q + q                 *')
12400 FORMAT(' *  for DIS:  scale > ',F7.3,' Q2                     *')
12500 FORMAT(' *  semihard approach for BGF a la Zotov et al       *')
12600 FORMAT(' *  gamma + gluon_',A4,
     +   ' --> Q Q_bar                   *')
12700 FORMAT(' *  heavy flavor produced is :',I3,'                   *')
12800 FORMAT(' *      W + gluon_',A4,
     +   ' --> c c_bar                   *')
12900 FORMAT(' *      gamma + pomeron --> ',A4,' + pomeron           *')
13000 FORMAT(' *  M.E.(gamma + pomeron --> ',A4,' + pomeron) = 1.E-4 *')
13100 FORMAT(' *  EPA + gamma gluon --> Q Q_bar (massive)   used   *')
13200 FORMAT(' *  full matrixelement e gluon --> e" Q Q_bar        *')
13300 FORMAT(' *  using HERACLES 4.6 for DIS                       *')
13400 FORMAT(' *  using HERACLES 4.6 for for heavy quark BGF        *')
13500 FORMAT(' *  HERACLES low q2 supp.: 1.-EXP(-',F4.2,
     +    '*Q2)          *')
13600 FORMAT(' *  gamma + q_',A7, ' --> q  (DIS diffractive)       *')
13700 FORMAT(' *      W + q_',A7, ' --> q  (DIS diffractive)       *')
13800 FORMAT(' *  vector meson production selected for IVM = ',I4)
13900 FORMAT(' *  QCD calc. for hard diffractive scattering used   *')
14000 FORMAT(' *  Nikolaev Zakharov model                          *')
*13600 FORMAT(' *  parametrisation of M. Wuesthoff  (DESY 95-166)   *')
14100 FORMAT(' *  M. Wuesthoff (ANL-HEP-PR 97-03)                  *')
14200 FORMAT(' *  Bartels, Lotter, Wuestoff (DESY 96-026)          *')
14300 FORMAT(' *  Bartels/Jung/Wuesthoff qqg for DIF DESY 99-027   *',/,
     +       ' *  gamma* p --> q q_bar g p (via 2 gluon exchange)  *',/,
     +       ' *  heaviest flavor produced:',I2,
     +       '                      *')
14400 FORMAT(' *  Bartels/Jung/Wuesthoff qqg for DIF DESY 99-027   *',/,
     +       ' *  gamma* p --> Q Q_bar g p (via 2 gluon exchange)  *',/,
     +       ' *  heaviest flavor produced:',I2,
     +       '                      *')
14500 FORMAT(' *   full integration selceted                       *')
14600 FORMAT(' *   approximate formula acc. eq.(4.5) DESY 99-027   *',/
     +       ' *   ordering in k_t: k_t_quark > k_t_gluon          *')
14700 FORMAT(' *   approximate formula acc. eq (4.5) DESY 99-027   *',/
     +       ' *   ordering in k_t: k_t_quark > k_t_gluon          *',/
     +       ' *   and suppression at large beta                   *')
14800 FORMAT(' *  Bartels, Lotter, Wuestoff (DESY 96-026)          *',/
     +       ' *  gamma* p --> q q_bar p (via 2 gluon exchange)    *',/,
     +       ' *  heaviest flavor produced:',I2,
     +       '                      *')
14900 FORMAT(' *  Bartels, Lotter, Wuestoff (DESY 96-026)          *',/
     +       ' *  gamma* p --> Q Q_bar p (via 2 gluon exchange)    *',/,
     +       ' *  heaviest flavor produced:',I2,
     +       '                      *')

15000 FORMAT(' *  M. Diehl with Donnachie/Landshoff 2 glu pomeron  *')
15100 FORMAT(' *  Buchmueller/McDermott/Hebecker DESY 97-099       *')
15200 FORMAT(' *  saturation model Golec-Biernat and Wuesthoff     *',/
     +       ' *  coded by H. Kowalski and M. Wuesthoff            *')
15300 FORMAT(' *  saturation model Golec-Biernat and Wuesthoff     *',/
     +       ' *  coded by H. Kowalski and M. Wuesthoff            *',/
     +       ' *  with rad. correction from HERACLES               *')
15400 FORMAT(' *  only phase space for gamma glu --> q q_bar g     *')
15500 FORMAT(' *  C1 = ',F6.2,' Cg = ',F6.2,
     + '                           *')
15600 FORMAT(' *  gamma glu glu --> q q_bar                        *',
     +     /,' *  gamma glu glu --> q q_bar glu                    *')
15700 FORMAT(' *  mixing of QPM and O(alpha_s) processes selected  *')
15800 FORMAT(' *      QCD weights calculated event by event        *')
15900 FORMAT(' *      QCD weights calculated in grid               *')
16000 FORMAT(' *  gamma + q_',A6, ' --> q  (DIS )                   *')
16100 FORMAT(' *      W + q_',A6, ' --> q  (DIS )                   *')
16200 FORMAT(' *  mixing of QPM and O(alpha_s) processes selected  *')
16300 FORMAT(' *  scattered electron theta_max = ',F7.3,'           *')
16400 FORMAT(' *  scattered electron theta_min = ',F7.3,'           *')
16500 FORMAT(' *  Q^2 _min = ',F7.3,'                               *')
16600 FORMAT(' *  Q^2 _max = ',F7.3,'                               *')
16700 FORMAT(' *  y _min = ',F7.3,'                                 *')
16800 FORMAT(' *  y _max = ',F7.3,'                                 *')
16900 FORMAT(' *  Nr of flavours in target sea  = ',I4,'             *')
17000 FORMAT(' *  Nr of flavours for QCDC  = ',I4,'                  *')
17100 FORMAT(' *  cut on p_t ^2  = ',f5.2,' for IPRO = ',I5,
     +       '          *')
17200 FORMAT(' *  no p_t ^2 cut imposed for IPRO = ',I4,
     +       '            *')
17300 FORMAT(' *  p_t^2 cut:',F5.2,' for IPRO=15 ',F5.2,
     +       ' for IPRO=13    *')
17400 FORMAT(' * gaussian intrinsic kt in photon:',F5.2,
     +'             *')
17401 FORMAT(' * dkt^2/(kt0^2 + kt^2) intrinsic kt in photon:      *',
     +/,' *  kt0 = ',F5.2,' ktm = ',F5.2,'                          *')
17402 FORMAT(' * dkt^2/(kt0^2 + kt^2)^2 intrinsic kt in photon:    *',
     +/,' *  kt0 = ',F5.2,' ktm = ',F5.2,'                          *')
17500 FORMAT(' *  intrinsic kt in proton:',F5.2,
     +'                     *')
17600 FORMAT(' ##    no initial state PS possible for NG,NPOM=30   ##',/
     +      ,' ##    switch to final state PS                      ##')
17700 FORMAT('  cm energy of ep system ',F10.3,' GEV')
17800 FORMAT('  Lambda_qcd = ',F5.3,' at NF= ',I2,' flavours')
17900 FORMAT(' * Nikolaev Zakharov process selected ')
18000 FORMAT(' * t_max  = ',F7.3,'                                  *')
18100 FORMAT(' * Rn**2=',F5.3,/,'   alpha_Po =',F5.3
     +    ,/,' * cut on x_f:  ',F5.3,/,'   epsilon =',F6.4)
18200 FORMAT(' * cut on x_f:  ',F5.3,
     +    '                                *')
18300 FORMAT(' * Nikolaev Zakharov process selected ')
18400 FORMAT(' * parton density in pomeron xf_',
     +             I1,' (x) = ',I1,'*(1-x)**',I1)
18500 FORMAT('    X     |   Q^2   |',
     +      ' F_2 (DIS) | F_2 (DIF) | F_2 (PI)   | ')
18600 FORMAT(' ',E8.2,' |',F8.2,' |  ',F8.4,' |   '
     +                   ,F8.4,' |   ',F8.4,' |')

      END
