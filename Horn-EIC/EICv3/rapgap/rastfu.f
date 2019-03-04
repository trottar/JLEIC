*CMZ :  2.08/05 28/03/2000  14.06.17  by  Hannes Jung
*CMZ :  2.08/00 23/06/99  08.08.05  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  09.56.21  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94

      SUBROUTINE RASTFU(KF,XRX,SCALA,XPQ)
      Implicit None
c note that xpq is parton density of pomeron for NG<20
c but diffractive parton density for NG>20: xpq = xpq(beta,x_pom,t2,scala)
c        also for user defined function USDIFFR(beta,scala,xpq,xpom,t2)
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
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
*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

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
cGB
      Real  XPGB
cGB
C...MODIFIED FOR USE IN HERACLES 4.2 (BY H.SPIESBERGER, 8.2.93)
      Real PYSTOP,PYSLAM
      Integer NPYMOD,NPYMAX,NPYMIN
      COMMON /PYSTFUC/ PYSTOP,PYSLAM,NPYMOD,NPYMAX,NPYMIN
C Input:
C   PYSTOP = top mass
C   NPYMAX = maximal flavour
C   NPYMIN = minimal flavour
C   NPYMOD = Choice of parametrization
      Integer LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      Integer LPAR,LPARIN,IPART
      COMMON /HSPARL/ LPAR(20),LPARIN(12),IPART
*KEEP,RGRAHER.
      REAL XPQDIF,XPQPI
	Integer IHERPYS
	Integer NBQ2,NBX
      PARAMETER (NBQ2=20)
      PARAMETER (NBX=20)
      COMMON /RAHER/ IHERPYS,XPQDIF(-6:6,NBX,NBQ2),XPQPI(-6:6,NBX,NBQ2)
*KEND.
      REAL COW(3,5,4,2),XQ(9),TS(6)
      DOUBLE PRECISION XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
C coded is xQ(X,SCALE) and xG(x,SCALE), not G(x) or Q(x)
      REAL XRX,X,SCALE,SCALA
      INTEGER KF
      REAL F2SIG
      LOGICAL FIRST,FIRS
      LOGICAL PDFFIRST
      COMMON /W50516/PDFFIRST

      REAL XPQ(-6:6)
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION BETA,F2QT,FLQ,X_POM,T2,XGX,XDX(-4:4)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      REAL PYGAMM
      EXTERNAL PYGAMM
      Real GEV2NB,DNC,DC,XP,EULBET,XPQU,DL1,DL2,DL3,Q2IN,ALPH_EM,SD
      Real POINT,Y,ALAM
      Real D2S_VL,D2S_VS,D2S_VC,D2S_SL,D2S_SS,D2S_SC
      Integer I,IP,IS,INIT,NPDF,NSET,KFL,KPHF
      DATA GEV2NB/.3893857E+6/
      DATA DNC/3./,DC/1./
      DATA XP/0.0/,INIT/0/,FIRST/.TRUE./,FIRS/.TRUE./
      DATA NPDF/0/
C...The following data lines are coefficients needed in the
C...Owens pion structure function parametrizations, see below.
C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1,1),IS=1,5),IP=1,3)/
     +  4.0000E-01,  7.0000E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     + -6.2120E-02,  6.4780E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     + -7.1090E-03,  1.3350E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/
      DATA ((COW(IP,IS,1,2),IS=1,5),IP=1,3)/
     +  4.0000E-01,  6.2800E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     + -5.9090E-02,  6.4360E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     + -6.5240E-03,  1.4510E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2,1),IS=1,5),IP=1,3)/
     +  8.8800E-01,  0.0000E+00,  3.1100E+00,  6.0000E+00,  0.0000E+00,
     + -1.8020E+00, -1.5760E+00, -1.3170E-01,  2.8010E+00, -1.7280E+01,
     +  1.8120E+00,  1.2000E+00,  5.0680E-01, -1.2160E+01,  2.0490E+01/
      DATA ((COW(IP,IS,2,2),IS=1,5),IP=1,3)/
     +  7.9400E-01,  0.0000E+00,  2.8900E+00,  6.0000E+00,  0.0000E+00,
     + -9.1440E-01, -1.2370E+00,  5.9660E-01, -3.6710E+00, -8.1910E+00,
     +  5.9660E-01,  6.5820E-01, -2.5500E-01, -2.3040E+00,  7.7580E+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3,1),IS=1,5),IP=1,3)/
     +  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00,
     + -2.4280E-01, -2.1200E-01,  8.6730E-01,  1.2660E+00,  2.3820E+00,
     +  1.3860E-01,  3.6710E-03,  4.7470E-02, -2.2150E+00,  3.4820E-01/
      DATA ((COW(IP,IS,3,2),IS=1,5),IP=1,3)/
     +  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00,
     + -1.4170E-01, -1.6970E-01, -2.4740E+00, -2.5340E+00,  5.6210E-01,
     + -1.7400E-01, -9.6230E-02,  1.5750E+00,  1.3780E+00, -2.7010E-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4,1),IS=1,5),IP=1,3)/
     +  0.0000E+00, -2.2120E-02,  2.8940E+00,  0.0000E+00,  0.0000E+00,
     +  7.9280E-02, -3.7850E-01,  9.4330E+00,  5.2480E+00,  8.3880E+00,
     + -6.1340E-02, -1.0880E-01, -1.0852E+01, -7.1870E+00, -1.1610E+01/
      DATA ((COW(IP,IS,4,2),IS=1,5),IP=1,3)/
     +  0.0000E+00, -8.8200E-02,  1.9240E+00,  0.0000E+00,  0.0000E+00,
     +  6.2290E-02, -2.8920E-01,  2.4240E-01, -4.4630E+00, -8.3670E-01,
     + -4.0990E-02, -1.0820E-01,  2.0360E+00,  5.2090E+00, -4.8400E-02/


C...Euler's beta function, requires ordinary Gamma function
      EULBET(X,Y)=PYGAMM(X)*PYGAMM(Y)/PYGAMM(X+Y)
      DO 10 I=-6,6
   10 XPQ(I) = 0.0
      IF(NG.LT.40) PT2GEN = 0.D0
      SCALE = SCALA
c      IF(IWEI.EQ.1) THEN
c         write(6,*) ' RASTFU: NG = ',NG,' NPOM = ',NPOM,' IPRO ',IPRO
c         write(6,*) ' RASTFU: KF = ',KF,' X ',X
c         ENDIF
      IF(FIRS) THEN
         write(6,*) ' RASTFU: NG = ',NG,' NPOM = ',NPOM
         FIRS=.FALSE.
      ENDIF
      X=XRX
      IF(X.LT.0.0.OR.X.GE.1) THEN
         write(6,*) ' RASTFU error: X = ',X
         IF(X.LT.0) X=0.
         IF(X.GE.1.) X=0.99999
         RETURN
      ENDIF
      IF(SCALE.LT.0.0) THEN
         write(6,*) ' RASTFU error: SCALE = ',SCALE
         SCALE = 0.0001
         RETURN
      ENDIF
      IF(NG.LT.0) THEN
         IF(NG.GE.-5.AND.NG.LE.-1) THEN
c ACTW Fit
            CALL ACTWFIT(X,SCALE,XPQ,XFGKI,T2GKI)
            IF(FIRS) THEN
               WRITE(6,*) ' ACTW pdf"s are used '
            ENDIF
            GOTO 120
         ELSEIF(NG.GE.-15.AND.NG.LE.-10) THEN
c H1 QCD Fit
            CALL H1QCDFIT(X,SCALE,XPQ,XFGKI,T2GKI)
            GOTO 120
         ELSEIF(NG.EQ.-20.AND.NG.EQ.-20) THEN
c Soft Colour Interaction approach from Buchmueller/Hebecker
            CALL SCAPDF(X,SCALE,XPQ,XFGKI,T2GKI)
            GOTO 120
         ELSE
            IF(FIRS) THEN
               WRITE(6,*) ' User defined routine USDIFFR used '
            ENDIF
            FIRST = .FALSE.
c      write(6,*) ' RASTFU: beta,scale,x_pom,t2',beta,scale,xfgki,t2gki
c      write(6,*) ' RASTFU: beta,x',beta,x
            CALL USDIFFR(X,SCALE,XPQ,XFGKI,T2GKI)
            GOTO 120
         ENDIF
      ENDIF
c         write(6,*) 'rastfu: x,scale',x,scale
      IF(NG.EQ.12) THEN

         XPQ(0) = 9./2 *X*(1.-X)
         XPQU = 3./10. *X*(1.-X)
C add pointlike contribution
c         write(6,*) ' log arg = ',SCALE*(1.-X)/0.3/0.3/X,x,scale
         DL1 = ALOG(SCALE*(1.-X)/0.3/0.3/X)
         IF(DL1.LE.0.) DL1=0.0
         DL2 = ALOG(SCALE*(1.-X)/0.3/0.3/X)
         IF(DL2.LE.0.) DL2=0.0
         DL3 = ALOG(SCALE*(1.-X)/0.5/0.5/X)
         IF(DL3.LE.0.) DL3=0.0
         POINT = DNC/8.0/SNGL(PI)**2 * DC**2*(X**2 + (1. - X)**2)
c          IF(POINT.LE.0) write(6,*) DL1,POINT,XPQ(1),X
         XPQ(1) = XPQU + X * POINT*DL1
         XPQ(-1)= XPQ(1)
         XPQ(2) = XPQU + X * POINT*DL2
         XPQ(-2)= XPQ(1)
         XPQ(3) = 3./20. *X*(1.-X) + X * POINT * DL3
         XPQ(-3) = XPQ(3)
c charm here is only for testing
c         xpq(4) = xpqu/10.
c         xpq(-4)= xpq(4)
c end testing...........................
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: 12 NG,NPOM',NG,NPOM
cGB
      ELSEIF(NG.EQ.13) THEN
         XPGB = (3./4.)*(X*(1.-X)+(0.57/2.)*(1.-X)**2)
         XPQ(0)  = 0.
         XPQ(1)  = XPGB
         XPQ(2)  = XPGB
         XPQ(3)  = XPGB
         XPQ(-1) = XPGB
         XPQ(-2) = XPGB
         XPQ(-3) = XPGB
cGB
      ELSEIF(NG.GT.0.AND.NG.LE.5) THEN
C.... xG_N (X)= (N+1)(1-X)**N
         XPQ(0) = FLOAT(NG+1)*(1.-X)**NG
         XPQ(1) = XPQ(0)/4.
         XPQ(2) = XPQ(0)/4.
         XPQ(-1) = XPQ(1)
         XPQ(-2) = XPQ(2)
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: <5 NG,NPOM',NG,NPOM
      ELSEIF(NG.EQ.10) THEN
         XPQ(0) = (0.18 + 5.46*X)*(1.-X)
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: 10 NG,NPOM',NG,NPOM
      ELSEIF(NG.EQ.0) THEN
C.... xG_0 (X)= 6x(1-X)
         XPQ(0) = 6.*(1.-X)*X
         XPQ(1) = XPQ(0)/4.
         XPQ(2) = XPQ(0)/4.
         XPQ(-1) = XPQ(1)
         XPQ(-2) = XPQ(2)
c        write(6,*) ' RASTFU: NG,NPOM',xpq(0)
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: 10 NG,NPOM',NG,NPOM
      ELSEIF(NG.EQ.11) THEN
         XP = SNGL(PI) * 0.16*x*(1-x)/3.
         DO 20 I=-2,2
   20    XPQ(I) = XP
         XPQ(0) = 0.
         XPQ(-3) = XP/2.
         XPQ(3) = XP/2.
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: 11 NG,NPOM',NG,NPOM
      ELSEIF(NG.GE.15.AND.NG.LT.20) THEN
c call parametrisation from file acc. to H1 data
         IF(NG.LE.17) THEN
            CALL RES_LS(X,SCALE,XPQ)
         ENDIF
      ELSEIF(NG.EQ.100) THEN
         IF(INIT.EQ.0) WRITE(6,*) ' SIMPLE SCALING FUNCTION USED'
         INIT=1
C...Simple scaling structure functions, valence quarks and gluon only.
C. this applies for protons....
         XPQ(0)=3.*(1.-X)**5
         XPQ(3)=0.125*(1.-X)**7
         XPQ(2)=70./32 *(1.-X)**3 *SQRT(X)
CCCC         XPQ(2)=(1.-X)**3*(1.274+.589*(1.-X)-1.675*(1.-X)**2)
         XPQ(1)=2.*XPQ(2)
         XPQ(-1)=XPQ(3)
         XPQ(-2)=XPQ(3)
         XPQ(-3)=XPQ(3)
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: 100 NG,NPOM',NG,NPOM
      ELSEIF(NG.EQ.20.OR.NG.EQ.21) THEN
C... xG(x) for pion exchange
         IF(FIRST) THEN
            IF(MSTP(52).LE.10) THEN
               WRITE(6,*) ' pi structure function in RASTFU '
            ELSE
               WRITE(6,*) ' pi structure function in PDFLIB'
            ENDIF
            FIRST = .FALSE.
         ENDIF
c        write(6,*) 'MSTP(52) = ',MSTP(52),KF
C...Call user-supplied structure function.
C...Proton structure function call.
C.. modify to call PDFLIB

         IF(MSTP(52).GE.10) THEN
C...Call PDFLIB structure functions.
            XX=DBLE(X)
            QQ=DBLE(SQRT(MAX(0.,SCALE)))
C if you use PDFLIB version 3.  then uncoment the following 2 lines
C and comment the lines with 'cold pdf'
cold pdf  PARM(1)='MODE'
cold pdf  VAL(1)=MSTP(52)
C if you use PDFLIB 4.   then use the next lines.
CNEW
            PARM(1) = 'NPTYPE'
            VAL(1) = DBLE(MSTP(52)/1000000.)
            PARM(2) = 'NGROUP'
            VAL(2) = DBLE(MOD(MSTP(52),1000000)/1000)
            PARM(3) = 'NSET'
            VAL(3) = DBLE(MOD(MOD(MSTP(52),1000000),1000))
CNEW
c            write(6,*) ' RATSFU : PDFLIB CALL ',NG,NPOM
            NPDF=NPDF+1
            PDFFIRST = .FALSE.
            IF(NPDF.LE.1) PDFFIRST=.TRUE.
c call PDFSET each time, because when DIS with normal p structure function is
c also selected one would get in confusion....
            CALL PDFSET(PARM,VAL)
            CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
            XPQ(0)=SNGL(GLU)
            XPQ(1)=SNGL(DNV+DSEA)
            XPQ(-1)=SNGL(DSEA)
            XPQ(2)=SNGL(UPV+USEA)
            XPQ(-2)=SNGL(USEA)
            XPQ(3)=SNGL(STR)
            XPQ(-3)=SNGL(STR)
            XPQ(4)=SNGL(CHM)
            XPQ(-4)=SNGL(CHM)
            XPQ(5)=SNGL(BOT)
            XPQ(-5)=SNGL(BOT)
            XPQ(6)=SNGL(TOP)
            XPQ(-6)=SNGL(TOP)
c            write(6,*) ' pion pdf ',(xpq(i),i=-3,3)
            IF(NG.EQ.21) THEN
c construct partons in pi0
               DO 30 I=1,6
   30          XPQ(I) = (XPQ(I) + XPQ(-I))/2.
            ENDIF
            GOTO 120
         ENDIF
C end modify to call PDFLIB


c        write(6,*) ' rastfu :',x,SCALE
C copied from pystfu-------------------------------------------------------
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
         NSET=1
c         IF(MSTP(51).EQ.2.OR.MSTP(51).EQ.4.OR.MSTP(51).EQ.13) NSET=2
         IF(NSET.EQ.1) ALAM=0.2
         IF(NSET.EQ.2) ALAM=0.4
c        IF(MSTP(52).LE.0) THEN
c           SD=0.
c        ELSE
         Q2IN=MIN(2E3,MAX(4.,SCALE))
         SD=LOG(LOG(Q2IN/ALAM**2)/LOG(4./ALAM**2))
c         ENDIF
C...Calculate structure functions.
         DO 50 KFL=1,4
            DO 40 IS=1,5
   40       TS(IS)=COW(1,IS,KFL,NSET)+COW(2,IS,KFL,NSET)*SD+ COW(3,IS,
     +      KFL,NSET)*SD**2
            IF(KFL.EQ.1) THEN
               XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/EULBET(TS(1),TS(2)+1.)
            ELSE
               XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*
     +         X**2)
            ENDIF
   50    CONTINUE

C...Put into output arrays.
         XPQ(0)=XQ(2)
         XPQ(1)=XQ(3)/6.
         XPQ(2)=XQ(1)+XQ(3)/6.
         XPQ(3)=XQ(3)/6.
         XPQ(4)=XQ(4)
         XPQ(-1)=XQ(1)+XQ(3)/6.
         XPQ(-2)=XQ(3)/6.
         XPQ(-3)=XQ(3)/6.
         XPQ(-4)=XQ(4)
c         IF(IWEI.EQ.1) write(6,*) ' RASTFU: pi NG,NPOM',NG,NPOM
         IF(NG.EQ.21) THEN
c construct partons in pi0
            DO 60 I=1,6
   60       XPQ(I) = (XPQ(I) + XPQ(-I))/2.
         ENDIF
      ELSEIF(NG.EQ.30) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' RASTFU: Nikolaev Zakharov process selected '
         ENDIF
         FIRST = .FALSE.
         DO 70  I=-6,6
            XPQ(I) = 0.0
   70    CONTINUE
*************************************************************************
*
* author: M.Grothe , 15.12.95
*
* D^2 SIGMA/(DT DXP) (GAMMA*P->XP) ACCORDING TO
* Genovese, Nikolaev, Zakharov, Cern/Th-95/13
*
* input : kinematic variables x_bj, x_pomeron, Q_sqrt, ABS(t)
* output: D^2 SIGMA/(DT DX_pomeron) in mbarn/GeV^2
*         D2S_VL  valence part with light quarks in final state
*         D2S_VS  valence part with strange quarks in final state
*         D2S_VC  valence part with charm quarks in final state
*         D2S_SL  sea part with light quarks in final state
*         D2S_SS  sea part with strange quarks in final state
*         D2S_SC  sea part with charm quarks in final state
*
*
*************************************************************************

         CALL NIKZAK(XPR,XFGKI,SNGL(Q2),ABS(T2GKI),
     +           D2S_VL,D2S_VS,D2S_VC,D2S_SL,D2S_SS,D2S_SC)
         ALPH_EM = SNGL(ALPH)
         IF(IRUNAEM.EQ.1) ALPH_EM = ULALEM(SNGL(Q2))

         F2SIG = SNGL(Q2/(4.D0*PI**2*DBLE(ALPH_EM)))/GEV2NB*1.E6
c sum_u,d e_q**2 = 10/9 for u and d quarks
         XPQ(1) = F2SIG*(D2S_VL+D2S_SL)*9./10.
         XPQ(-1) = XPQ(1)
         XPQ(2) = XPQ(1)
         XPQ(-2) = XPQ(1)
         XPQ(3) = F2SIG*(D2S_VS+D2S_SS)*9./2.
         XPQ(-3) = XPQ(3)
         XPQ(4) = F2SIG*(D2S_VC+D2S_SC)*9./8.
         XPQ(-4) = XPQ(4)
c         write(6,*) 'RASTFU: XPR,XFGKI,T2GKI,Q2',XPR,XFGKI,T2GKI,Q2
c         write(6,*) 'XPQ ',(XPQ(I),I=-4,4)
      ELSEIF(NG.EQ.40) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' RASTFU: M. Wuesthoff model used '
         ENDIF
         FIRST = .FALSE.
         beta = dble(xpr/xfgki)
         x_pom = dble(xfgki)
         T2 = dble(t2gki)
         if(beta.lt.0.0.or.beta.gt.1.) then
            write(6,*) ' RASTFU: beta,x,xpr,xfgki ',beta,x,xpr,xfgki
         endif
         CALL F2DHW(BETA,x_pom,Q2,T2,F2QT,FLQ,XGX)
         DO 80  I=-3,3
            XPQ(I) = sngl(F2QT)*9./6.
   80    CONTINUE
         XPQ(0) = SNGL(XGX)
c         write(6,*) 'rastfu: xpq',xpq
      ELSEIF(NG.EQ.41) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' RASTFU: Bartels, Lotter, Wuesthoff DESY 96-'
     +      //'026 '
         ENDIF
         FIRST = .FALSE.
         beta = dble(xpr/xfgki)
         x_pom = dble(xfgki)
         T2 = dble(t2gki)
         if(beta.lt.0.0.or.beta.gt.1.) then
            write(6,*) ' RASTFU: beta,x,xpr,xfgki ',beta,x,xpr,xfgki
         endif
         CALL F2BLW(BETA,x_pom,Q2,T2,F2QT,FLQ)
         DO 90  I=-3,3
c divide by 2 because of quark and antiquark include in calc but explicitely her
            XPQ(I) = sngl(F2QT)*9./6./2.
   90    CONTINUE
         XPQ(0) = 0.0
c         write(6,*) 'rastfu: xpq',xpq
      ELSEIF(NG.EQ.42) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' RASTFU: M. Diehl calc. used '
         ENDIF
         FIRST = .FALSE.
         beta = dble(xpr/xfgki)
         x_pom = dble(xfgki)
         T2 = dble(t2gki)
         if(beta.lt.0.0.or.beta.gt.1.) then
            write(6,*) ' RASTFU: beta,x,xpr,xfgki ',beta,x,xpr,xfgki
         endif
         CALL F2MD(BETA,x_pom,Q2,T2,XDX,FLQ)
         DO 100 I=-4,4
            XPQ(I) = sngl(XDX(I))
  100    CONTINUE
         XPQ(0) = 0.0
c         write(6,*) 'rastfu: xpq',xpq
      ELSEIF(NG.EQ.45) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' RASTFU: Buchmueller/McDermott/Hebecker
     +        calc. used '
         ENDIF
         FIRST = .FALSE.
         beta = dble(xpr/xfgki)
         x_pom = dble(xfgki)
         T2 = dble(t2gki)
         if(beta.lt.0.0.or.beta.gt.1.) then
            write(6,*) ' RASTFU: beta,x,xpr,xfgki ',beta,x,xpr,xfgki
         endif
         CALL F2MCD(BETA,x_pom,Q2,T2,XDX,FLQ)
         DO 110 I=-4,4
            XPQ(I) = sngl(XDX(I))
  110    CONTINUE
c         write(6,*) 'rastfu: xpq',xpq
      ELSE
         WRITE(6,*) ' requested parton distribution not implemented'
         WRITE(6,*) ' NG = ',NG,' NPOM = ',NPOM
      ENDIF
  120 CONTINUE
      DO 130 I=-6,6
         IF(XPQ(I).LT.0.0) XPQ(I) = 0.0
  130 CONTINUE
      IF(IHF.EQ.1) THEN
         DO 140 I=1,6
            KPHF = 4
            IF(I.EQ.KPHF) GOTO 140
            XPQ(I) = 0.0
            XPQ(-I) = 0.0
  140    CONTINUE
      ENDIF
      RETURN
      END
