*CMZ :  2.08/05 10/03/2000  12.33.04  by  Hannes Jung
*CMZ :  2.08/04 11/01/2000  19.47.19  by  Hannes Jung
*CMZ :  2.08/02 04/11/99  16.39.08  by  Hannes Jung
*CMZ :  2.08/01 22/06/99  16.58.44  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.55  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  19.54.45  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE EVENT
      IMPLICIT NONE
      Integer NDIM,NPOIN
      COMMON/DIVO/ NDIM,NPOIN
*KEEP,RGEFFIC.
      DOUBLE PRECISION AVGI,SD
      INTEGER NIN,NOUT
      COMMON/EFFIC/AVGI,SD,NIN,NOUT
C     SAVE

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


*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGLUPARM.
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
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
*KEEP,RGRGRIDF2.
      DOUBLE PRECISION XX,Q2X
      COMMON /RGRIDF2/ XX(NBX),Q2X(NBQ2)
      REAL F2_DIS,F2_DIF,F2_PI
      COMMON /F2VAL/ F2_DIS(NBX,NBQ2),F2_DIF(NBX,NBQ2),F2_PI(NBX,NBQ2)
cNEW
      DOUBLE PRECISION F2DIS,F2DIF,F2PI
      COMMON /F2INT/ F2DIS,F2DIF,F2PI
cNEW
*KEND.
c henry's commons
      Double Precision  Stot, Q2_henry, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2_henry, W2, Mx2, AJAC
      Double Precision    KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
c end
      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)
      REAL ULMASS
      EXTERNAL ULMASS

      REAL XPQ(-6:6)
      Double Precision X(20)
      Double Precision XG(20)
      Integer LST,IRES,MSTA
      COMMON /EPPARA/ LST(30),IRES(2)
      REAL PARA
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      Double Precision WMAX
      Integer IMIX
      COMMON /OALPHAS/ WMAX,IMIX

      Integer NQPM,NQQB,NQQBC,NQQBB,NQCDC
      COMMON /QCDEV/ NQPM,NQQB,NQQBC,NQQBB,NQCDC
      Double Precision XV
      Integer NDIMEN
      COMMON /XVAL/ XV(20),NDIMEN
      Double Precision draprn
      DOUBLE PRECISION MP
      DOUBLE PRECISION POMAX,PIMAX,WTMA
      COMMON/WEIGHT/ POMAX(NBX,NBQ2),PIMAX(NBX,NBQ2),WTMA
      Integer IGENFL,IVM
      COMMON/GENWEI/IGENFL
      COMMON/VMESON/IVM
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      Integer NDIS,NQPMDS,NQQBDS,NQQBCDS,NQQBBDS,NQCDCDS,
     +              NDIF,NQPMDF,NQQBDF,NQQBCDF,NQQBBDF,NQCDCDF,
     +              NPI,NQPMPI,NQQBPI,NQQBCPI,NQQBBPI,NQCDCPI,
     +              NRPA,NRPB,NRPC,NRPD,NRPE,NRPF,NRPG,NPRT
      COMMON/NEVOUT/NDIS,NQPMDS,NQQBDS,NQQBCDS,NQQBBDS,NQCDCDS,
     +              NDIF,NQPMDF,NQQBDF,NQQBCDF,NQQBBDF,NQCDCDF,
     +              NPI,NQPMPI,NQQBPI,NQQBCPI,NQQBBPI,NQCDCPI,
     +              NRPA,NRPB,NRPC,NRPD,NRPE,NRPF,NRPG,NPRT
      REAL PYPAR,PYVAR
      Integer IPY
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      REAL QMAX
      Integer LERR
      COMMON/ERR/LERR(100)
C JPP Add variables for calculating Q2EG
      DOUBLE PRECISION KEP(4)
      REAL Q2EG
      REAL XHF,Q2HF
      Integer KPHF
      COMMON /HEAVYF/XHF,Q2HF,KPHF
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      REAL SNGL
      Real     remPARJ32
      Real     remPARJ(5)
      LOGICAL FIRST

      Integer IQPML1,IQQBL1,IQQBCL1,IQCDCL1,IQPML0
      Integer IQQBL0,IQCDCL0,IQPMH1,IQQBH1,IQQBCH1
      Integer IQCDCH1,IQPMH0,IQQBH0,IQCDCH0
      Integer NTRMAX,NPRIN,IDEBUGO,IPROO,IDIRO,NGO,NPOMO,KPAO,KINT22
      Integer IHERPYSO,IHFO
      Integer I,J,NTRY,IX,IQ,IERR,IKEEP
      Integer NFI1,NFI2
      Integer IRAD
      Double Precision WDUM,SMALL,T2MAXO,XBJ,WPA,XPRT,XRMIN,XR
      Double Precision VMAX,VMIN,XMAXV,XMINV,xmch
      Double Precision XMINHF,T2MIN,XP2T,Q2T,T2,T2MN,WTDIST,WTT,WTRN
      Double Precision FXN1,DCORR,XD,QD,X1P,X2P,RF2DIS,RF2DIF,RF2PI
      Double Precision RNTEST,DOT1,RNSAT
      EXTERNAL DOT1
      EXTERNAL draprn
      DATA WDUM/0.D0/
      DATA IQPML1/0/,IQQBL1/0/,IQQBCL1/0/,IQCDCL1/0/,IQPML0/0/
      DATA IQQBL0/0/,IQCDCL0/0/,IQPMH1/0/,IQQBH1/0/,IQQBCH1/0/
      DATA IQCDCH1/0/,IQPMH0/0/,IQQBH0/0/,IQCDCH0/0/
      DATA FIRST/.TRUE./
      DATA NTRMAX/8000/
      DATA NPRIN/5/
      DATA  SMALL/1.D-3/
      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2

CPER
      integer nPrint
      parameter (nPrint = 100000)
CPER

      IDEBUGO = IDEBUG
   10 CONTINUE
      NIN=NIN+1
c      IF(MOD(NOUT,1).EQ.0) write(6,*) ' start: Nr of events ',NOUT
      IF(FIRST) THEN
         first = .FALSE.
         remPARJ32 = PARJ(32)
         IPROO = IPRO
         IDIRO = IDIR
         NGO = NG
         NPOMO = NPOM
         KPAO = KPA
         KINT22 = KINT(2,2)
         FIRST=.FALSE.
         MP = DP(2,5)
         T2MAXO = T2MAX
         NDIS=0
         NQPMDS=0
         NQQBDS=0
         NQQBCDS=0
         NQQBBDS=0
         NQCDCDS=0
         NDIF=0
         NQPMDF=0
         NQQBDF=0
         NQQBCDF=0
         NQQBBDF=0
         NQCDCDF=0
         NPI=0
         NQPMPI=0
         NQQBPI=0
         NQQBCPI=0
         NQQBBPI=0
         NQCDCPI=0
         NRPA=0
         NRPB=0
         NRPC=0
         NRPD=0
         NRPE=0
         NRPF=0
         NRPG=0
         NPRT=0
         IHERPYSO = IHERPYS
         IHFO = IHF
      ELSE
         PARJ(32)=remPARJ32
         IHERPYS = IHERPYSO
         IHF = IHFO
         IPRO = IPROO
         IDIR=IDIRO
         NG = NGO
         NPOM = NPOMO
         KPA = KPAO
         KINT(2,2) = KINT22
      ENDIF
      T2GKI = 0.
      XFGKI = 0.
      PT2H = 0.
      SHH = 0.
      IWEI = 0
      IGENFL = 0
      IMIX = 0
      DO 20 I=1,20
         XV(I)= 0.D0
   20 X(I) = 0.D0
      NDIMEN = NDIM
      IF(IDISDIF.GE.1) IHERPYS = 0
      IF(IPRO.EQ.1200.AND.IDIR.EQ.1) IHERPYS=0
      IF(IPRO.EQ.1200.AND.IDIR.EQ.0) IHERPYS=1
      IF(IPRO.EQ.3000.AND.IDIR.EQ.0) IHERPYS=1
      IF(IPRO.LT.1000) THEN
         CALL RAPGEN(NDIM,XG)
         DO 30 I=1,NDIM
            X(I) = XG(I)
   30    XV(I) = X(I)
c      write(6,*) ' check kt2 after',kt2
      ELSE
         CALL HERACL(AVGI,SD,2)
         IHERPYS = 0
c       IF(MOD(NOUT,1).EQ.0) write(6,*) ' heracl: Nr of events ',NOUT

         XBJ = XHS
         Q2 = Q2HS
         IF(IHF.EQ.1) THEN
            KPHF = IHFLA
            KPA = KPHF
         ENDIF
         IF(IDIR.EQ.1) THEN
            WPA = FXN1(X,WDUM)
         ELSEIF(IDIR.EQ.0.AND.NPOM.NE.20.AND.NPOM.NE.21) THEN
            XPRT = XBJ
		Q2T = Q2
            IF(XPRT.LT.XX(1).OR.XPRT.GT.XX(NBX)) THEN
               WRITE(6,*) ' EVENT: X values outside grid '
               WRITE(6,*) ' X_min ',XX(1),' X_max ',XX(NBX),
     +          ' actual X ', XPRT
               IF(XPRT.LT.XX(1)) XPRT=XX(1)
               IF(XPRT.GT.XX(NBX)) XPRT=XX(NBX)
            ENDIF
            IX = 0
   40       IX = IX + 1
            IF(XPRT.GT.XX(IX+1)) GOTO 40
            IF(Q2T.LT.Q2X(1).OR.Q2T.GT.Q2X(NBQ2)) THEN
               WRITE(6,*) ' EVENT: Q2 values outside grid '
               WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NBQ2),
     +          ' actual Q2 ', Q2T
               IF(Q2T.LT.Q2X(1)) Q2T=Q2X(1)
               IF(Q2T.GT.Q2X(NBQ2)) Q2T=Q2X(NBQ2)
            ENDIF
            IQ = 0
   41       IQ = IQ + 1
            IF(Q2T.GT.Q2X(IQ+1)) GOTO 41
            NTRY = 0
            NG = NGO
            NPOM =NPOMO
            IDIR = 0
            KINT(2,2) = 100
            WTMA = POMAX(IX,IQ)

            XD = (XPRT - XX(IX))/(XX(IX+1)-XX(IX))
            QD = (Q2T - Q2X(IQ))/(Q2X(IQ+1) - Q2X(IQ))
            X1P=DBLE(
     +      (POMAX(IX+1,IQ)-POMAX(IX,IQ))*SNGL(XD)+POMAX(IX,IQ))
            X2P=DBLE((POMAX(IX+1,IQ+1)-POMAX(IX,IQ+1))
     +                *SNGL(XD)+POMAX(IX,IQ+1))
            WTMA = (X2P-X1P)*QD + X1P


c            NDIMEN = NDIM + 2
c            write(6,*) ' event: IPRO = ',IPRO,IHERAC
            IF(IPRO.EQ.30) THEN
c               NTRMAX = 15000
   50          CALL draprnV(XG,2)
               IF(NTRY.GT.NTRMAX) THEN
                  write(6,*) 'too many trials for DIF; event rejected,'
     +            //'q2 = ',Q2,' x = ',xbj
                  GOTO 10
               ENDIF
               NTRY = NTRY + 1
               DO  I=1,2
                  XV(NDIM+I)=XG(I)
                  X(NDIM+I)=XG(I)
               ENDDO
c           write(6,*) ' event x ',x(ndim+1),x(ndim+2),xg(1),xg(2)
               call RGSATRAP(X,WPA)
               RNSAT = draprn()
c               write(6,*) ' EVENT :Wpa,WTMA ',WPA,WTMA,RNSAT,ntry
               IF(WPA.GT.1.0*WTMA) THEN
                  write(6,*) ' EVENT: WPA > WTMA from RGSATRAP',WPA,
     +            WTMA,' q2,x ',q2,xbj,ix,iq
               endif

               IF(RNSAT.GT.WPA/(1.0*WTMA)) GOTO 50
c             write(6,*) ' EVENT accepted:Wpa,WTMA ',WPA,WTMA,RNSAT,ntry
            ELSE
c generate xr and t2
   60          CALL draprnV(XG,2)
               NTRY = NTRY + 1
               DO 70 I=1,2
                  XV(NDIM+I)=XG(I)
   70          X(NDIM+I)=XG(I)
               XMAX = 1.D0 - XF
               XP2T = XBJ
               IF(XMAX.LE.XP2T) THEN
                  write(6,*) ' 1. xmax<xp2t ',xmax,xp2t
                  GOTO 10
               endif
               IF(NTRY.GT.NTRMAX) THEN
                  write(6,*) 'too many trials for DIF; event rejected,'
     +            //'t,x'
                  GOTO 10
               ENDIF
               XRMIN = XP2T
cnew
c            IF(IHF.EQ.1) THEN
c               XMINHF = XP2T*(1.d0 + 4.d0*DBLE(ULMASS(KPA))**2/Q2)
c           Xrmin = MAX(XrMIN,XMINHF)
c           Xrmin = MIN(Xrmin,XMAX)
cc          write(6,*) ' 1st xminhf ',xminhf,xrmin,xmax
c           ENDIF
cnew
               IF(IVM.EQ.0) THEN
                  XR = XRMIN*(XMAX/XRMIN)**X(NDIM+1)
c           write(6,*) ' event xr,xmax,xrmin ',xr,xmax,xrmin
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
c                  write(6,*) 'event:',XMAX,XMIN
                  XMAXV = XP2T*(1.d0 + VMAX**2/Q2)
                  XMINV = XP2T*(1.d0 + VMIN**2/Q2)
                  XR = XMINV*(XMAXV/XMINV)**X(NDIM+1)
c                  write(6,*) 'event:',XMAXV,XMINV,XR
                  IF(XMAXV.GT.XMAX) THEN
                     write(6,*) 'xmaxv>xmax'
                     GOTO 10
                  ENDIF
               ENDIF
               IF(IHF.EQ.1) THEN
                  XMINHF = XP2T*(1.d0 + 4.d0*DBLE(ULMASS(KPA))**2/Q2)
c keep XR generated in full range, and only reject those values which fall
c outside the phase space, otherwise POMAX calculated in POMSTR will not match
C maximum calculated here
                  IF(XMINHF.GT.XMAX) THEN
                     write(6,*) ' event xminhf,xrmin',xminhf,xrmin
                     write(6,*) ' event xp2t,q2',xp2t,q2,ulmass(4)
                     write(6,*) ' goto 10 '
                     GOTO 10
                  ENDIF
                  IF(XMINHF.GT.XR) THEN
c             write(6,*) 'xminhf > xr ',xminhf,xr,xmin,xmax
c             write(6,*) ' goto 50 '
                     GOTO 60
                  ENDIF

c               XR = XMINHF*(XMAX/XMINHF)**X(NDIM+1)
               ENDIF
               T2MIN = MP*MP*XR*XR/(1.D0-XR)
               T2MN = -(Q2 *(1.D0 - XR/XBJ) + 4.D0*DBLE(ULMASS(211))**
     +         2)
               IF(T2MN.LT.T2MAXO) THEN
                  T2MAX = T2MN
               ELSE
                  T2MAX = T2MAXO
               ENDIF
               IF(T2MAX.LE.T2MIN) THEN
c               write(6,*) ' event: T2min T2max',T2MIN,T2MAX
c           write(6,*) ' goto 50 '
                  GOTO 60
               ENDIF
               T2 = T2MIN*((T2MAX/T2MIN)**X(NDIM+2))
               WTDIST = 0.D0
               CALL RAT2DI(KINT(2,2),XR,-T2,WTDIST)
c        write(6,*) ' event WTDIST',WTDIST
               WTDIST = WTDIST * T2*DLOG(T2MAX/T2MIN)
               IF(IVM.EQ.0.AND.IHF.EQ.0) THEN
                  WTDIST=WTDIST*XR*DLOG(XMAX/XRMIN)
               ELSEIF(IVM.EQ.0.AND.IHF.EQ.1) THEN
ccc               WTDIST=WTDIST*XR*DLOG(XMAX/XMINHF)
                  WTDIST=WTDIST*XR*DLOG(XMAX/XRMIN)
                  xmch = q2*(xr-xp2t)/xp2t - T2
                  if(xmch.lt.4.d0*dble(ulmass(kphf))**2) Then
c              write(6,*) ' xmch failed xmch= ',xmch,t2,kphf,kpa
c              write(6,*) ' goto 50 '
                     goto 60
                  endif
c        write(6,*) ' event WTDIST',WTDIST,XR,XR*DLOG(XMAX/XRMIN)
               ELSE
                  WTDIST=WTDIST*XR*DLOG(XMAXV/XMINV)
               ENDIF
               IF(DABS(WTDIST).LE.1.D-45) GOTO 60
               XPR = SNGL(XP2T)
               T2GKI = SNGL(-T2)
               XFGKI = SNGL(XR)
               CALL RASTFU(100,SNGL(XP2T/XR),SNGL(Q2),XPQ)
               WTDIST = WTDIST*DBLE(XPQ(KPA))
               IF(WTDIST.LE.0.0) GOTO 60
c add safety margin 1.5 for rejection
               WTT = WTDIST/(1.5D0*WTMA)
               WTRN = draprn()
C JPP Turn debugging off
c            IDEBUG = 1
               WPA = FXN1(X,WDUM)
c            IDEBUG =IDEBUGO
c               write(6,*) ' T2,Xr,xbj ',t2,xr,xbj
c               write(6,*) ' x(NDIM+1),X(NDIM+2),NDIM ',X(NDIM+1),
c     +       X(NDIM+2),NDIM
               IF(WTDIST.GT.1.5D0*WTMA) THEN
                  write(6,*) ' WPA,NTRY ',WPA,NTRY
                  write(6,*) ' EVENT: WTDIST,WTMA ',WTDIST,WTMA
                  write(6,*) ' WTRN,WTT ',WTRN,WTT
                  write(6,*) ' T2,Xr,xbj ',t2,xr,xbj
                  write(6,*) ' x(NDIM+1),X(NDIM+2),NDIM ',X(NDIM+1),
     +            X(NDIM+2),NDIM
               endif
               IF(WPA.LE.0.) THEN
c               write(6,*) '1st WPA = 0 '
                  GOTO 10
               ENDIF
               IF(WTRN.GT.WTT) THEN
c               write(6,*) ' WTRN>WTT ',wtrn,wtt
                  GOTO 60
               ENDIF
c            write(6,*) ntry
            ENDIF
         ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
            NTRY = 0
            IDIR = 0
            XPRT = XBJ
		Q2T = Q2
            IF(NPOM.EQ.20) THEN
               NG = 20
               KINT(2,2)=211
            ELSEIF(NPOM.EQ.21) THEN
               NG = 21
               KINT(2,2)=111
            ENDIF
            XPRT = XBJ
            IF(XPRT.LT.XX(1).OR.XPRT.GT.XX(NBX)) THEN
               WRITE(6,*) ' EVENT: X  values outside grid '
               WRITE(6,*) ' X_min ',XX(1),' X_max ',XX(NBX),
     +         ' actual X ', XPRT
               IF(XPRT.LT.XX(1)) XPRT=XX(1)
               IF(XPRT.GT.XX(NBX)) XPRT=XX(NBX)
            ENDIF
            IX = 0
   80       IX = IX + 1
            IF(XPRT.GT.XX(IX+1)) GOTO 80
            IF(Q2T.LT.Q2X(1).OR.Q2T.GT.Q2X(NBQ2)) THEN
               WRITE(6,*) ' EVENT: Q2 values outside grid '
               WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NBQ2),
     +          ' actual Q2 ', Q2T
               IF(Q2T.LT.Q2X(1)) Q2T=Q2X(1)
               IF(Q2T.GT.Q2X(NBQ2)) Q2T=Q2X(NBQ2)
            ENDIF
            IQ = 0
   81       IQ = IQ + 1
            IF(Q2T.GT.Q2X(IQ+1)) GOTO 81

            WTMA = PIMAX(IX,IQ)
            NDIMEN = NDIM + 2
c generate xr and t2
   90       CALL draprnV(XG,2)
            DO 100  I=1,2
               XV(NDIM+I)=XG(I)
  100       X(NDIM+I)=XG(I)
            NTRY = NTRY + 1
            XMAX = 1.D0 - XF
            XP2T = XBJ
            IF(XMAX.LE.XP2T) THEN
               write(6,*) ' 2 xmax<x2pt ',xmax,xp2t
               GOTO 10
            endif
            IF(NTRY.GT.NTRMAX) THEN
               write(6,*) ' too many trials for PI; event rejected,x,t'
               GOTO 10
c               RETURN
            ENDIF
            XR = XP2T*(XMAX/XP2T)**X(NDIM+1)
            T2MIN = MP*MP*XR*XR/(1.D0-XR)
            T2MN = -(Q2 *(1.D0 - XR/XBJ) + 4.D0*DBLE(ULMASS(211))**2)
c            write(6,*) T2MN,T2MAXO,T2min
c            write(6,*) Q2,Xr,XBJ,ULMASS(211)
            IF(T2MN.LT.T2MAXO) THEN
               T2MAX = T2MN
            ELSE
               T2MAX = T2MAXO
            ENDIF
            IF(T2MAX.LE.T2MIN) GOTO 90

            T2 = T2MIN*((T2MAX/T2MIN)**X(NDIM+2))
            WTDIST = 0.D0
            CALL RAT2DI(KINT(2,2),XR,-T2,WTDIST)
            WTDIST = WTDIST *T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)
            IF(DABS(WTDIST).LE.1.D-45) GOTO 90
cnew
            XPR = SNGL(XP2T)
            T2GKI = SNGL(-T2)
            XFGKI = SNGL(XR)
            CALL RASTFU(KINT(2,2),SNGL(XP2T/XR),SNGL(Q2),XPQ)
            WTDIST = WTDIST*DBLE(XPQ(KPA))
cnew
c add safety margin 1.5 for rejection
            WTT = WTDIST/(1.5D0*WTMA)
            WTRN = draprn()
c		IDEBUGO = IDEBUG
c            IDEBUG = 1
               WPA = FXN1(X,WDUM)
c            IDEBUG =IDEBUGO
            IF(WTDIST.GT.1.5D0*WTMA) THEN
               write(6,*) ' WPA,NTRY ',WPA,NTRY
               write(6,*) ' EVENT: WTDIST,WTMA ',WTDIST,WTMA
               write(6,*) ' WTRN,WTT ',WTRN,WTT
               write(6,*) ' T2,Xr,xbj ',t2,xr,xbj
               write(6,*) ' x(NDIM+1),X(NDIM+2),NDIM ',X(NDIM+1),
     +       X(NDIM+2),NDIM
            endif
            IF(WPA.LE.0.0) THEN
c               write(6,*) ' 3rd WPA = 0'
c               write(6,*) T2MN,T2MAXO,T2min
c            write(6,*) Q2,Xr,XBJ,ULMASS(211),DBLE(ULMASS(211))
c            write(6,*)  Q2*(1.D0 - XR/XBJ),4.D0*ULMASS(211)**2
               GOTO 10
            ENDIF
            IF(WTRN.GT.WTT) THEN
               GOTO 90
            ENDIF

         ENDIF
      ENDIF

c      CALL DULIST(1)

      IGENFL = 0
      IDEBUG = 0
c       write(6,*) 'EVENT: before 1st FXN1, WPA= ',WPA,IPRO,KPA,idir,x
      IF(IPRO.EQ.30) THEN
      ELSEIF(ISEMIH.EQ.0) THEN
         WPA = FXN1(X,WDUM)
         IF(WPA.LE.0.0) THEN
c         write(6,*) ' 2nd after rangen :',(x(i),i=1,ndim)
c         write(6,*) ' 4th WPA = 0',KPA
c         call dulist(1)
c         pause
            GOTO 10
         ENDIF
      ENDIF
      IDEBUG = IDEBUGO
c       write(6,*) 'EVENT: after 1st FXN1, WPA= ',WPA,IPRO,KPA,idir,x

      IF(IPRO.EQ.100) GOTO 180
      XBJ = Q2/DBLE(YY)/SSS
      IF(IDISDIF.GE.1) THEN
         IGENFL = 0
         XPRT = XBJ
         Q2T = Q2
         IF(XPRT.LT.XX(1).OR.XPRT.GT.XX(NBX) .OR.Q2T.LT.Q2X(1)
     +   .OR.Q2T.GT.Q2X(NBQ2)) THEN
            IF(IHERAC.EQ.0) THEN
               WRITE(6,*) ' EVENT: X or Q2 values outside grid '
               WRITE(6,*) ' X_min ',XX(1),' X_max ',XX(NBX),
     +         ' actual X ', XPRT
               WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NBQ2), ' '
     +         //'actual Q2 ',Q2T
            ENDIF
            IF(XPRT.LT.XX(1)) XPRT=XX(1)
            IF(XPRT.GT.XX(NBX)) XPRT=XX(NBX)
            IF(Q2T.LT.Q2X(1)) Q2T = Q2X(1)
            IF(Q2T.GT.Q2X(NBQ2)) Q2T = Q2X(NBQ2)
         ENDIF
         IX = 0
  110    IX = IX + 1
         IF(XPRT.GT.XX(IX+1)) GOTO 110
         IQ = 0
  120    IQ = IQ + 1
         IF(Q2T.GT.Q2X(IQ+1)) GOTO 120
         DCORR = 1.D0
         XD = (XPRT - XX(IX))/(XX(IX+1)-XX(IX))
         QD = (Q2T - Q2X(IQ))/(Q2X(IQ+1) - Q2X(IQ))
         X1P=DBLE(
     +    (F2_DIS(IX+1,IQ)-F2_DIS(IX,IQ))*SNGL(XD)+F2_DIS(IX,IQ))
         X2P=DBLE(
     +    (F2_DIS(IX+1,IQ+1)-F2_DIS(IX,IQ+1))*SNGL(XD)+F2_DIS(IX,IQ+1))
         F2DIS = (X2P-X1P)*QD + X1P
         F2DIS = F2DIS*DCORR
         X1P=DBLE(
     +    (F2_DIF(IX+1,IQ)-F2_DIF(IX,IQ))*SNGL(XD)+F2_DIF(IX,IQ))
         X2P=DBLE(
     +    (F2_DIF(IX+1,IQ+1)-F2_DIF(IX,IQ+1))*SNGL(XD)+F2_DIF(IX,IQ+1))
         F2DIF = (X2P-X1P)*QD + X1P
         F2DIF = F2DIF*DCORR
         X1P=DBLE(
     +    (F2_PI(IX+1,IQ)-F2_PI(IX,IQ))*SNGL(XD)+F2_PI(IX,IQ))
         X2P=DBLE(
     +    (F2_PI(IX+1,IQ+1)-F2_PI(IX,IQ+1))*SNGL(XD)+F2_PI(IX,IQ+1))
         F2PI = (X2P-X1P)*QD + X1P
         F2PI = F2PI*DCORR
         IF((1.D0 - XF).LE.XBJ) THEN
            F2DIF = 0.D0
            F2PI = 0.D0
         ENDIF
         RF2DIS = (F2DIS-F2DIF-F2PI)/F2DIS
         IF(RF2DIS.LE.0.0D0) THEN
            write(6,*) ' pomeron or pion contribution larger than '
     +      //'total'
            write(6,*) ' F2TOT= ',F2DIS,' F2_DIF= ',F2DIF,' F2_PI= ',
     +      F2PI
            write(6,*) ' event rejected ..... mixing might be wrong '
            RETURN
         ENDIF
         RF2DIF = F2DIF/F2DIS
         RF2PI = F2PI/F2DIS
         KPA = KPAO
         RNTEST = draprn()
         IF(RNTEST.LT.RF2DIS) THEN
            IDIR = 1
            KINT(2,2)=2212
            NDIMEN = NDIM
            WTMA = 999999.D0
         ELSEIF(RNTEST.LT.(RF2DIS+RF2DIF)) THEN
            NTRY = 0
            NG = NGO
            NPOM =NPOMO
            IDIR = 0
            KINT(2,2) = 100
            WTMA = POMAX(IX,IQ)
            NDIMEN = NDIM + 2
c generate xr and t2
  130       CALL draprnV(XG,2)
            NTRY = NTRY + 1
            DO 140 I=1,2
               XV(NDIM+I)=XG(I)
  140       X(NDIM+I)=XG(I)
            XMAX = 1.D0 - XF
            XP2T = XBJ
            IF(XMAX.LE.XP2T) THEN
               write(6,*) '3. xmax<xp2t ',xmax,xp2t
               GOTO 10
            endif
            IF(NTRY.GT.NTRMAX) THEN
               write(6,*) ' too many trials for DIF;',
     +          ' event rejected,mix,po'
               GOTO 10
c               RETURN
            ENDIF
            XR = XP2T*(XMAX/XP2T)**X(NDIM+1)
            T2MIN = MP*MP*XR*XR/(1.D0-XR)
            T2MN = -(Q2 *(1.D0 - XR/XBJ) + 4.D0*DBLE(ULMASS(211))**2)
            IF(T2MN.LT.T2MAXO) THEN
               T2MAX = T2MN
            ELSE
               T2MAX = T2MAXO
            ENDIF

            T2 = T2MIN*((T2MAX/T2MIN)**X(NDIM+2))
            WTDIST = 0.D0
            CALL RAT2DI(KINT(2,2),XR,-T2,WTDIST)
            WTDIST = WTDIST * T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)
c include that when ng <0 or ng>30 because weight is only included in RASTFU,
c in RAT2DI weight=1.
cnew
            XPR = SNGL(XP2T)
            T2GKI = SNGL(-T2)
            XFGKI = SNGL(XR)
            CALL RASTFU(100,SNGL(XP2T/XR),SNGL(Q2),XPQ)
            WTDIST = WTDIST*DBLE(XPQ(KPA))
cnew

c add safety margin 1.5 for rejection
            WTT = WTDIST/(1.5D0*WTMA)
            WTRN = draprn()
            IF(WTRN.GT.WTT) THEN
               GOTO 130
            ENDIF
            WPA = FXN1(X,WDUM)
            IF(WPA.LE.0.0) GOTO 130
         ELSEIF(RNTEST.LT.(RF2DIS+RF2DIF+RF2PI)) THEN
            NTRY = 0
            IDIR = 0
            NG = 20
            NPOM = 20
            KINT(2,2)=211
            WTMA = PIMAX(IX,IQ)
            NDIMEN = NDIM + 2
c generate xr and t2
  150       CALL draprnV(XG,2)
            DO 160 I=1,2
               XV(NDIM+I)=XG(I)
  160       X(NDIM+I)=XG(I)
            NTRY = NTRY + 1
            XMAX = 1.D0 - XF
            XP2T = XBJ
            IF(XMAX.LE.XP2T) THEN
               write(6,*) '4. xmax<xp2t ',xmax,xp2t
               GOTO 10
            endif
            IF(NTRY.GT.NTRMAX) THEN
               write(6,*) ' too many trials for PI;',
     +         '  event rejected,mix,pi'
               GOTO 10
c               RETURN
            ENDIF
            XR = XP2T*(XMAX/XP2T)**X(NDIM+1)
            T2MIN = MP*MP*XR*XR/(1.D0-XR)
            T2MN = -(Q2 *(1.D0 - XR/XBJ) + 4.D0*DBLE(ULMASS(211))**2)
            IF(T2MN.LT.T2MAXO) THEN
               T2MAX = T2MN
            ELSE
               T2MAX = T2MAXO
            ENDIF

            T2 = T2MIN*((T2MAX/T2MIN)**X(NDIM+2))
            WTDIST = 0.D0
            CALL RAT2DI(KINT(2,2),XR,-T2,WTDIST)
            WTDIST = WTDIST * T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)
c include that when ng <0 or ng>30 because weight is only included in RASTFU,
c in RAT2DI weight=1.
cnew
            XPR = SNGL(XP2T)
            T2GKI = SNGL(-T2)
            XFGKI = SNGL(XR)
            CALL RASTFU(100,SNGL(XP2T/XR),SNGL(Q2),XPQ)
            WTDIST = WTDIST*DBLE(XPQ(KPA))
cnew
c add safety margin 1.5 for rejection
            WTT = WTDIST/(1.5D0*WTMA)
            WTRN = draprn()
            IF(WTRN.GT.WTT) THEN
               GOTO 150
            ENDIF
            WPA = FXN1(X,WDUM)
            IF(WPA.LE.0.0) GOTO 150
         ELSE
            WRITE(6,*) ' something"s wrong here '
            WRITE(6,*) ' RNTEST = ',RNTEST,' RF2DIS = ',RF2DIS, ' '
     +      //'RF2DIF = ',RF2DIF,' RF2PI = ',RF2PI
         ENDIF
         WPA = FXN1(X,WDUM)
      ENDIF


      IGENFL = 0
c      write(6,*)'here before MEPS mixing',WPA,KPA
      DO 170 I=1,NDIMEN
  170 XV(I) = X(I)
c check whether event mixing wanted by user
      IHERPYS = 0
      IMIX = IFULL
      IWEI = 1
c      write(6,*) ' 1. event xpr,xfgki',xpr,xfgki
c      write(6,*) ' event ',imix,ndim,ndimen
      IF(IMIX.EQ.1) THEN
c         WPA = FXN1(X,WDUM)
c         write(6,*) ' wpa before qcdmix ',wpa,x
         CALL QCDMIX(X,IERR)
c check on QCDMIX
         IF(LST(21).NE.0) THEN
            LERR(LST(21)) = LERR(LST(21)) + 1
         ENDIF
         IF(IERR.NE.0) THEN
            write(6,*) ' qcdmix error ',ierr
            GOTO 10
         ENDIF
      ENDIF
c      write(6,*) ' 2. event xpr,xfgki',xpr,xfgki
      IGENFL = 1
c zero MEINFO
      ZQGKI = 99999.
      XPGKI = 99999.
      PHITGKI = 99999.
c............
c      IF(IPRO.NE.20) THEN
      IF(IPRO.EQ.20.OR.IPRO.EQ.30) THEN
      ELSE
         IDEBUG = 0
c       write(6,*) 'EVENT: before 2nd FXN1, WPA= ',WPA,IPRO,KPA,idir,x
         IF(ISEMIH.EQ.0) THEN
            WPA = FXN1(X,WDUM)
            IDEBUG=IDEBUGO
c       write(6,*) 'EVENT: after 2nd FXN1, WPA= ',WPA,IPRO,KPA,idir,x
            IF(WPA.LE.0.0D0) THEN
c           write(6,*) 'EVENT: after 2nd FXN1, WPA= ',WPA,IPRO,KPA,idir,x
               GOTO 10
            ENDIF
         ENDIF
      ENDIF
  180 CONTINUE
c      write(6,*) ' 3. event xpr,xfgki',xpr,xfgki
      IWEI = 0
c      write(6,*) ' event NIA1 ',nia1
c      IF(IPRO.EQ.20.OR.IPRO.EQ.21) THEN
      IF(IPRO.EQ.30) THEN
c      write(6,*) ' check kt2 before satrev',kt2
         CALL RGSATREV
c        write(6,*) '  event before reshuffling'
c        write(6,*) ' event nia1,nia2,nf1,nf2 ',nia1,nia2,nf1,nf2
c        call dulist(1)
         IRAD = 0
         IF(K(4,2).eq.22.and.K(4,1).eq.1) IRAD=1
         DO I=1,5
            DP(N+1,I) = DP(3,I)
            K(N+1,I) = K(3,I)
            K(N+2,I) = K(4,I)
            DP(N+2,I) = DP(4,I)
            DP(3,I) = DP(NIA1,I)
            K(3,I) = K(NIA1,I)
            IF(IRAD.eq.1) then
               K(5,I) = K(N+2,I)
               DP(5,I) = DP(N+2,I)
            endif
            DP(4,I) = DP(N+1,I)
            K(4,I) = K(N+1,I)
            DP(N+1,I) = 0
            K(N+1,I) = 0
            DP(N+2,I) = 0
            K(N+2,I) = 0
         ENDDO
         nia1 = nia1 - 1 - irad
c        write(6,*) ' irad ',irad
c         write(6,*) ' event nf1,nf2,shh',nf1,nf2,mx2,sm2
c         call dulist(1)
         IF(IPY(13).EQ.1) THEN
            QMAX = SQRT(Mx2)/2.
            IF(NF2.EQ.NF1+1) THEN
               NFI2 = N
               NFI1 = N-1
            ELSE
               QMAX = SQRT(SM2)
               NFI2 = N
               NFI1 = N-2
            ENDIF
c        write(6,*) ' qmax ',qmax,mx2,sm2
            CALL DUSHOW(NFI1,NFI2,QMAX)
         ENDIF
         call rghaprep

         CALL DUEDIT(14)
      ELSEIF(IPRO.EQ.20.OR.IPRO.EQ.21) THEN
         DO I=1,5
            DP(N+1,I) = DP(3,I)
            K(N+1,I) = K(3,I)
            DP(3,I) = DP(4,I)
            K(3,I) = K(4,I)
            DP(4,I) = DP(N+1,I)
            K(4,I) = K(N+1,I)
            DP(N+1,I) = 0
            K(N+1,I) = 0
         ENDDO
         IF(IPY(13).EQ.1) THEN
c            write(6,*) ' event nf1,nf2,shh',nf1,nf2,shh
            IF(IPRO.EQ.21) THEN
               QMAX = 0.5*SQRT(SHH)
            ELSEIF(IPRO.EQ.20) THEN
               QMAX=SNGL(DOT1(NF1,NF2))
            ENDIF
            CALL DUSHOW(NF1,NF2,QMAX)
         ENDIF
         call rghaprep
         K(NF1,1)=21
         K(NF1+1,1)=21
         K(NF2,1)=21
         CALL DUEDIT(14)
      ELSE

         IF(IFPS.NE.10) THEN
            IF(IPRO.NE.100) THEN
c             call dulist(1)
c            write(6,*) ' before LMEPS '
               IHF=0
               IF(IPRO.EQ.99) THEN
                  CALL LMEPS99
               ELSE
                  CALL LMEPS
               ENDIF
               IHF = IHFO
c            write(6,*) ' after LMEPS '
c check on lmeps and pyremn
               IF(LST(21).NE.0) THEN
                  LERR(LST(21)) = LERR(LST(21)) + 1
                  GOTO 10
               ENDIF
            ELSEIF(IPRO.EQ.100) THEN
C...Delete empty lines
               CALL DUEDIT(12)
            ENDIF
         ELSE
            CALL LMEPS
c check on lmeps and pyremn
            IF(LST(21).NE.0) THEN
               LERR(LST(21)) = LERR(LST(21)) + 1
               GOTO 10
            ENDIF
            MSTA(5)=0
c         write(6,*) ' event: NF1,NF2 ',NF1,NF2
            DO 210 I=1,N
               DO 210 J=1,5
  210       P(I,J) = SNGL(DP(I,J))
c         write(6,*) ' before arexec '
            CALL AREXEC
            DO 220 I=1,N
               DO 220 J=1,5
  220       DP(I,J) = DBLE(P(I,J))
         ENDIF
      ENDIF
c     write(6,*) ' before prodiff'
      IF(NFRAG.EQ.10.AND.NPOM.NE.20.AND.NPOM.NE.21
     +           .AND.IDIR.EQ.0) CALL PRODIFF
c change interacting gamma to z0 on request from experiment
      DO 230 I=3,N
         IF(K(I,2).EQ.22.AND.K(I,1).EQ.21) K(I,2)=23
  230 CONTINUE
c     write(6,*) ' end of event nia1,nia2,nf1,nf2 ',nia1,nia2,nf1,nf2
c      call lulist(1)

      CALL DUDBRB(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
      DO 240  I=1,N
         DO 240  J=1,5
  240 P(I,J) = SNGL(DP(I,J))
      CALL LUDBRB(0,N,0.,0.,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4))
      IF(NFRAG.GE.1) CALL LUEXEC
      CALL LUDBRB(0,N,0.,0.,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
c only for decay of  p - diss
      IF(NFRAG.EQ.10.AND.NPOM.NE.20.AND.NPOM.NE.21
     +           .AND.IDIR.EQ.0) CALL PDISDC
c....
      IF(MSTU(24).NE.0) THEN
         WRITE(6,*) 'MSTU(24)= ',MSTU(24)
         CALL LULIST(1)
      ENDIF
      NOUT = NOUT + 1
      IF(MOD(NOUT,nPrint).EQ.0) write(6,*) ' Nr of events ',NOUT
      IF(XBJ.LT.0.001) THEN
         IF(IDIR.EQ.1.AND.IPRO.EQ.12) IQPML1 = IQPML1+1
         IF(IDIR.EQ.1.AND.IPRO.EQ.13) IQQBL1 = IQQBL1+1
         IF(IDIR.EQ.1.AND.IPRO.EQ.14) IQQBCL1 = IQQBCL1+1
         IF(IDIR.EQ.1.AND.IPRO.EQ.15) IQCDCL1 = IQCDCL1+1
         IF(IDIR.EQ.0.AND.IPRO.EQ.12) IQPML0 = IQPML0+1
         IF(IDIR.EQ.0.AND.IPRO.EQ.13) IQQBL0 = IQQBL0+1
         IF(IDIR.EQ.0.AND.IPRO.EQ.15) IQCDCL0 = IQCDCL0+1
      ELSE
         IF(IDIR.EQ.1.AND.IPRO.EQ.12) IQPMH1 = IQPMH1+1
         IF(IDIR.EQ.1.AND.IPRO.EQ.13) IQQBH1 = IQQBH1+1
         IF(IDIR.EQ.1.AND.IPRO.EQ.14) IQQBCH1 = IQQBCH1+1
         IF(IDIR.EQ.1.AND.IPRO.EQ.15) IQCDCH1 = IQCDCH1+1
         IF(IDIR.EQ.0.AND.IPRO.EQ.12) IQPMH0 = IQPMH0+1
         IF(IDIR.EQ.0.AND.IPRO.EQ.13) IQQBH0 = IQQBH0+1
         IF(IDIR.EQ.0.AND.IPRO.EQ.15) IQCDCH0 = IQCDCH0+1
      ENDIF
      IF(IMIX.EQ.1) THEN
         IF(IPRO.EQ.12) THEN
            NQPM = NQPM + 1
         ELSEIF(IPRO.EQ.13) THEN
            NQQB = NQQB + 1
         ELSEIF(IPRO.EQ.14) THEN
            IF(AM(1).LT.2.D0) NQQBC = NQQBC + 1
            IF(AM(1).GT.2.D0) NQQBB = NQQBB + 1
         ELSEIF(IPRO.EQ.15) THEN
            NQCDC = NQCDC + 1
         ENDIF
      ENDIF
      IF(IPRO.EQ.18) THEN
         NPRT = NPRT + 1
         IF(IRESPRO.EQ.1) NRPA=NRPA+1
         IF(IRESPRO.EQ.2) NRPB=NRPB+1
         IF(IRESPRO.EQ.3) NRPC=NRPC+1
         IF(IRESPRO.EQ.4) NRPD=NRPD+1
         IF(IRESPRO.EQ.5) NRPE=NRPE+1
         IF(IRESPRO.EQ.6) NRPF=NRPF+1
         IF(IRESPRO.EQ.7) NRPG=NRPG+1
      ENDIF
      IF(IDISDIF.GE.1) THEN
c         write(6,*) IDIR,NG,NPOM
         IF(NG.EQ.20.AND.NPOM.EQ.20) IDIR =2
         IF(IDIR.EQ.1) THEN
            NDIS = NDIS + 1
            IF(IPRO.EQ.12) THEN
               NQPMDS = NQPMDS + 1
            ELSEIF(IPRO.EQ.13) THEN
               NQQBDS = NQQBDS + 1
            ELSEIF(IPRO.EQ.14) THEN
               IF(AM(1).LT.2.D0) NQQBCDS = NQQBCDS + 1
               IF(AM(1).GT.2.D0) NQQBBDS = NQQBBDS + 1
            ELSEIF(IPRO.EQ.15) THEN
               NQCDCDS = NQCDCDS + 1
            ENDIF
         ELSEIF(IDIR.EQ.0) THEN
            NDIF = NDIF + 1
            IF(IPRO.EQ.12) THEN
               NQPMDF = NQPMDF + 1
            ELSEIF(IPRO.EQ.13) THEN
               NQQBDF = NQQBDF + 1
            ELSEIF(IPRO.EQ.14) THEN
               IF(AM(1).LT.2.D0) NQQBCDF = NQQBCDF + 1
               IF(AM(1).GT.2.D0) NQQBBDF = NQQBBDF + 1
            ELSEIF(IPRO.EQ.15) THEN
               NQCDCDF = NQCDCDF + 1
            ENDIF
         ELSEIF(IDIR.EQ.2) THEN
            NPI = NPI + 1
            IF(IPRO.EQ.12) THEN
               NQPMPI = NQPMPI + 1
            ELSEIF(IPRO.EQ.13) THEN
               NQQBPI = NQQBPI + 1
            ELSEIF(IPRO.EQ.14) THEN
               IF(AM(1).LT.2.D0) NQQBCPI = NQQBCPI + 1
               IF(AM(1).GT.2.D0) NQQBBPI = NQQBBPI + 1
            ELSEIF(IPRO.EQ.15) THEN
               NQCDCPI = NQCDCPI + 1
            ENDIF
         ENDIF
      ENDIF
      IF(MOD(NOUT,nPrint*5).EQ.0) THEN
         write(6,*) ' Nr of events ',NOUT,NQPM,NQQB
         write(6,*) ' Nr of DIS events x<10**-3 ',IQPML1,IQQBL1,IQCDCL1
         write(6,*) ' Nr of DIS events x>10**-3 ',IQPMH1,IQQBH1,IQCDCH1
         write(6,*) ' Nr of DIF events x<10**-3 ',IQPML0,IQQBL0,IQCDCL0
         write(6,*) ' Nr of DIF events x>10**-3 ',IQPMH0,IQQBH0,IQCDCH0
         IF(IDISDIF.GE.1) THEN
            write(6,*) ' all qpm qqb ccb qcdc '
            write(6,*) ' Nr of DIS events ',NDIS,NQPMDS,NQQBDS,NQQBCDS,
     +       NQQBBDS,NQCDCDS
            write(6,*) ' Nr of DIF events ',NDIF,NQPMDF,NQQBDF,NQQBCDF,
     +       NQQBBDF,NQCDCDF
            write(6,*) ' Nr of pi  events ',NPI,NQPMPI,NQQBPI,NQQBCPI,
     +       NQQBBPI,NQCDCPI
         ENDIF
      ENDIF
c      write(6,*) ' EVENT: x,Q2,XR',XBJ,Q2,XR
c      write(6,*) ' event: IPRO = ',IPRO,KPA
      IF(NOUT.LE.NPRIN) THEN
         write(6,*) ' event: IPRO = ',IPRO
         CALL LULIST(1)
      ENDIF
      IF(MOD(NOUT,nPrint*5).EQ.0) THEN
         write(6,*) 'Nr qpm ',NQPM
         write(6,*) 'Nr bgf ',NQQB
         write(6,*) 'Nr ccbar',NQQBC
         write(6,*) 'Nr bbbar',NQQBB
         write(6,*) 'Nr qcdc ',NQCDC
      endif
      NG = NGO
      NPOM = NPOMO
      RETURN
      END
