

*CMZ :  2.08/00 14/06/99  19.13.57  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  21.48.57  by  Hannes Jung
*CMZ :          15/05/99  21.37.38  by  Hannes Jung
* changes done for low mass diffraction by G. Briskin
*CMZ :          17/12/96  12.37.25  by  G.M. Briskin (TAU)
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE LMEPS99
* for low mass systems changes introduced by gena briskin 2.12.1996
      IMPLICIT NONE
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

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
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

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
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

*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEND.
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
c this is the modified but original LMEPS routine LEPTO 6.1
      DOUBLE PRECISION DETOT
      DOUBLE PRECISION DROBO(5)
      REAL CUT,X,Y,W2,Q2LP,U
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /RAPTOU/ CUT(14),LLST(40),PARL(30),X,Y,W2,Q2LP,U
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      DOUBLE PRECISION DBETA(2,3),DBETAR(2,3)
      DOUBLE PRECISION STHETA(2),SPHI(2),STHETAR(2),SPHIR(2)
      REAL PYPAR,PYVAR
      Integer IPY
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      REAL XPY,SH,TH,UH,Q2PY
      Integer ISUB,KFL
      COMMON /PYPROC/ ISUB,KFL(3,2),XPY(2),SH,TH,UH,Q2PY
      REAL XQPY
      COMMON /MYINT1/ XQPY(2,-6:6)
c.hju DIMENSION KS(9,5),PS(9,5),ROBO(5),XPQ(-6:6)
      REAL XPQ1,XPQ2
      DOUBLE PRECISION ROBO
      DIMENSION ROBO(5),XPQ1(-6:6),XPQ2(-6:6)
      DOUBLE PRECISION PS(40,5)
      REAL QMAX,XPRO
c..hju
      Integer LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
      DOUBLE PRECISION DELTAP(4), DPLONG, DBTOT, DGAMMA
c..hju
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      Integer KS,NS
      COMMON /COLR/ KS(40,5),NS
      REAL SNGL
cGB
      Real     remPARJ32
      Integer I,J,NCALL,NPRIN,NCHECK,NIPH,NB2,NSTB,IFL,IT,IPU1,IPU2
      Integer NIEL,IRFL,NIA1T,IPOM,MSTJO
      Double Precision SHRAP,Q2LM,UMAS2,XR,DETOTR,XGTEST,Q12
      Double Precision PC1,PC2,PC3,PC4,CHEC,P1,P2,P3,P4
      Double Precision phielt,PHIROT,PZT,STHARD
*PER
        integer iLMEPSPrint, MxLMEPSPrint
        parameter (MxLMEPSPrint = 10)
        save iLMEPSPrint
*PER

      logical first
      data    first /.TRUE./
c
c      SAVE KS,PS
cGB

      DATA NCALL/0/
      DATA NPRIN/0/
      DATA NCHECK/0/
      If(first) Then
         first = .FALSE.
         remPARJ32 = PARJ(32)
*PER
         iLMEPSPrint = 0
*PER
      ELSE
         PARJ(32)=remPARJ32
      EndIf
cGB
      DO I=1,40
         DO J=1,5
            KS(I,J) = 0
            PS(I,J) = 0.D0
         ENDDO
      ENDDO

c..hju
      NCALL = NCALL + 1
      IPY(48)=0
      LST(21)=0
C ILEPTO=1 parton shower a la LEPTO
C ILEPTO=0 parton shower a la PYTHIA
      SHRAP=DBLE(SHH)
      Q2PY = SNGL(Q2Q)
      Q2LP = SNGL(Q2)
      Q2LM = Q2
      IF(IRES(1).EQ.1) THEN
         DO 10 I=1,N
            IF(K(I,1).EQ.21.AND.K(I,2).EQ.22.AND.K(I,3).EQ.1) NIPH = I
            Q2LM = 0.0d0
   10    CONTINUE
      ELSE
         NIPH = NIA1
      ENDIF
c change to e in + z direction
      DO 20  I=1,N
         P(I,3) = - P(I,3)
   20 continue
      NB2 = 2
      UMAS2 = DBLE(ULMASS(K(2,2)))**2
c.hju parl(21) = 2P.k = invariant mass s
c.hjutest      PARL(21) = 2.D0*DOT1(1,NB2)
      PARL(21) = 2.D0*DOT1(1,2)
      IF(ILEPTO.EQ.0) PARL(21) =PARL(21)+DBLE(ULMASS(K(1,2)))**2 +UMAS2
c.hju parl(22) = 2P.q P = proton(or pomeron) q = photon
      PARL(22) = 2.D0*DOT1(NB2,NIPH)
      IF(IABS(K(1,2)).EQ.22.OR.IABS(K(1,2)).EQ.23) Q2LM=0.0D0

      XPRO = XPR
      PYVAR(32) = XPRO
      XR = DBLE(XPRO)
      PYVAR(31) = XEL
      IF(NCALL.LE.NPRIN) CALL DULIST(1)
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' LMEPS: 1st '
         CALL DULIST(1)
      ENDIF
      IF(ILEPTO.EQ.1) THEN
C...Transform to gamma p(or pomeron) cms, boost parameters in double precision.
         DETOT=P(NIPH,4)+P(NB2,4)
         DBETA(2,1)=(P(NIPH,1)+P(NB2,1))/DETOT
         DBETA(2,2)=(P(NIPH,2)+P(NB2,2))/DETOT
         DBETA(2,3)=(P(NIPH,3)+P(NB2,3))/DETOT
         CALL DUDBRB(0,0,0.D0,0.D0,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
         SPHI(2)=DLANGL(P(NIPH,1),P(NIPH,2))
         CALL DUDBRB(0,0,0.D0,-SPHI(2),0.D0,0.D0,0.D0)
         STHETA(2)=DLANGL(P(NIPH,3),P(NIPH,1))
         CALL DUDBRB(0,0,-STHETA(2),0.D0,0.D0,0.D0,0.D0)
      ENDIF
c         write(6,*) ' LMEPS: in gamma p ',NB2
c         CALL DULIST(1)
      IF(ILEPTO.EQ.1.and.IRES(1).EQ.1) THEN
         NSTB = 5
         DETOTR=P(NIA1,4)+P(NIA2,4)
         DBETAR(2,1)=(P(NIA1,1)+P(NIA2,1))/DETOTR
         DBETAR(2,2)=(P(NIA1,2)+P(NIA2,2))/DETOTR
         DBETAR(2,3)=(P(NIA1,3)+P(NIA2,3))/DETOTR
         CALL DUDBRB(nstb,0,0.D0,0.D0,-DBETAR(2,1),-DBETAR(2,2),
     +   -DBETAR(2,3))
         SPHIR(2)=DLANGL(P(NIA1,1),P(NIA1,2))
         CALL DUDBRB(nstb,0,0.D0,-SPHIR(2),0.D0,0.D0,0.D0)
         STHETAR(2)=DLANGL(P(NIA1,3),P(NIA1,1))
         CALL DUDBRB(nstb,0,-STHETAR(2),0.D0,0.D0,0.D0,0.D0)

         if(ncall.le.nprin) call dulist(1)
      ENDIF
c      IF(NCHECK.EQ.1) THEN
c         write(6,*) ' LMEPS: 2nd '
c         CALL DULIST(1)
c      ENDIF

      IF(NCALL.LE.NPRIN) CALL DULIST(1)
C...Save event record in  cms
      DO 30 I=1,N
         DO 30 J=1,5
            KS(I,J)=K(I,J)
   30 PS(I,J)=P(I,J)
C...Rearrange event record to PYSSPA standard
      DO 40 J=1,5
         K(3,J)=0
         P(3,J)=0.D0
         K(4,J)=0
         P(4,J)=0.D0
         K(5,J)=KS(NIA1,J)
         P(5,J)=PS(NIA1,J)
         IF(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13) THEN
            K(7,J)=KS(3,J)
            P(7,J)=PS(3,J)
            K(9,J)=KS(3,J)
            P(9,J)=PS(3,J)
         ENDIF
         K(8,J)=KS(NF1,J)
         P(8,J)=PS(NF1,J)

         K(11,J) = 0
         K(12,J) = 0

C check for radiative gamma
         IF(KS(4,2).EQ.22.AND.KS(4,1).EQ.1) THEN
            K(13,J) = KS(4,J)
            P(13,J) = PS(4,J)
         ENDIF
         K(11,J)=KS(NF1+1,J)
         P(11,J)=PS(NF1+1,J)
         K(10,J)=KS(NF2,J)
   40 P(10,J)=PS(NF2,J)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(9,3)=5
      K(10,3)=6
      DO 50 I=5,11
   50 K(I,1)=21
      K(9,1)=0
C...Incoming parton
      DO 60 J=1,4
         IF(ILEPTO.EQ.1) THEN
            IF(IPRO.NE.12) THEN

               P(6,J)=P(8,J) + P(10,J) + P(11,J)  - P(5,J)
            ENDIF
         ELSE
            P(6,J)=PS(NIA2,J)
         ENDIF
c this redinition of incoming and outgoing parton is needed
c because in PYREMN there is a boost to correct for the long. frame in case
c there are different masses for outgoing parton and remnant.

   60 CONTINUE
c end of redefintion
      P(6,5)=0.0D0
      K(6,2)=KS(NIA2,2)
      N=15
      IF(NCALL.LE.NPRIN) THEN
         CALL DULIST(1)
      ENDIF
C...Partons with colour information in hadronic cms frame.
      DO 70  I=16,MSTU(4)
         DO 70 J=1,5
            K(I,J)=0
            P(I,J)=0.D0
   70 V(I,J)=0.
      NS=20

      DO 80  J=1,5
         N=NS+7
         K(NS+1,J)=K(5,J)
         P(NS+1,J)=P(5,J)
         K(NS+3,J)=K(6,J)
         P(NS+3,J)=P(6,J)
         K(NS+5,J)=K(8,J)
         P(NS+5,J)=P(8,J)
         K(NS+6,J)=K(11,J)
         P(NS+6,J)=P(11,J)
         K(NS+7,J)=K(10,J)
         P(NS+7,J)=P(10,J)
   80 CONTINUE
C...Old standard continuation lines
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      P(NS+4,3)=28
      P(NS+4,4)=28
C...Origin and Colour info for incoming parton
      K(NS+3,1)=13
      K(NS+3,3)=2
      K(NS+3,4)=28
      K(NS+3,5)=28
C...Colour info for two outgoing partons
      K(NS+5,1)=3
      K(NS+6,1)=3
      K(NS+7,1)=3
C...Effective outgoing parton = sum of both outgoing partons
      K(NS+8,1)=14
      K(NS+8,3)=3
C... q qbar gluon event
      K(NS+5,4)=(NS+6)*MSTU(5)
      K(NS+5,5)=(NS+8)*MSTU(5)
      K(NS+6,4)=(NS+7)*MSTU(5)
      K(NS+6,5)=(NS+5)*MSTU(5)
      K(NS+7,4)=(NS+8)*MSTU(5)
      K(NS+7,5)=(NS+6)*MSTU(5)


      K(NS+5,4)=(NS+6)*MSTU(5)
      K(NS+5,5)=0
      K(NS+6,4)=(NS+8)*MSTU(5)
      K(NS+6,5)=(NS+5)*MSTU(5)
      K(NS+7,4)=0
      K(NS+7,5)=(NS+8)*MSTU(5)
      K(NS+8,2)=21
      IF(K(NS+5,2).GT.0) THEN
         K(NS+8,4)=(NS+3)*MSTU(5)+26
         K(NS+8,5)=(NS+3)*MSTU(5)+27
      ELSE
         K(NS+8,4)=(NS+3)*MSTU(5)+27
         K(NS+8,5)=(NS+3)*MSTU(5)+26
      ENDIF
      DO 90 J=1,4
         P(NS+8,J)=P(8,J)+P(10,J)+P(11,J)
   90 CONTINUE
      P(NS+8,5)=DSQRT(DMAX1(0.D0, P(NS+8,4)**2-P(NS+8,1)**2-P(NS+
     +8,2)**2- P(NS+8,3)**2))
      N=NS+8

CPHI
C...Scale for bremsstrahlung etc.

      SH = SNGL(SHRAP)
      IPY(40)=10
      IF(IPRO.EQ.12) IPY(40)=8
      IPY(47)=N
C...Save quantities for later use.
      IPY(41)=K(1,2)
      IF(IRES(1).EQ.1) IPY(41)=22
      IPY(42)=KINT(2,2)
      XPY(1)=1.
      XPY(2)=SNGL(XR)
      IF(IRES(1).EQ.1) THEN
         IPY(41)=22
         PYVAR(2) = -Q2LP +SNGL(PARL(22))
         PYVAR(31)=SNGL(SHRAP)/PYVAR(2)/PYVAR(32)
         XPY(1)=PYVAR(31)
         IF(XPY(1).GT.1.1) THEN
            write(6,*) ' LMEPS : XPY(1) > 1. ',XPY(1)
            write(6,*) ' LMEPS :s_hat, W**2, x_p',
     +        SHRAP,PYVAR(2),PYVAR(32)
            write(6,*) ' LMEPS :y, sss, q2lp, xel, xel/yy ',
     +        yy,sss,q2lp,xel,xel/yy
         ENDIF
         XPY(1)=AMIN1(1.0,XPY(1))
         CALL PYSTFU(IPY(41),XPY(1),Q2PY,XPQ1)
         IF(NCALL.LE.NPRIN) THEN
            write(6,*) ' IPY(41),XPY(1),Q2PY,XPQ1', IPY(41),XPY(1),
     +      Q2PY,XPQ1
         ENDIF
      ENDIF

c.hju CALL PYSTFU(K(2,2),XR,Q2,XPQ)
      IF(XPY(2).GT.1.) write(6,*) 'LMEPS : XPY(2) > 1. ',XPY(2)

      CALL PYSTFU(IPY(42),XPY(2),Q2PY,XPQ2)
      DO 100 I=-6,6
         IF(XPQ2(I).EQ.0.0) XPQ2(I)=1.E-11
  100 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' IPY(42),XPY(2),Q2PY,XPQ2',
     +    IPY(42),XPY(2),Q2PY,XPQ2
      ENDIF
      DO 110 IFL=-6,6
         XQPY(1,IFL)=XPQ1(IFL)
  110 XQPY(2,IFL)=XPQ2(IFL)
      ISUB=39
      IPY(11)=1
      KFL(2,1)=KS(NIA1,2)
      KFL(2,2)=K(6,2)
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      KFL(3,1)=K(1,2)
      KFL(3,2)=K(28,2)
      IF(IPRO.EQ.12) KFL(3,2)=K(8,2)
      PYVAR(2)=SNGL(PARL(21))
      IF(IRES(1).EQ.1) PYVAR(2) = -Q2LP +SNGL(PARL(22))
      PYVAR(3)=SNGL(P(1,5))
      PYVAR(4)=SNGL(P(2,5))
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(5)=PYVAR(1)/2.


C...Generate timelike parton shower (if required)
      IF(IPY(13).EQ.1) THEN
         CALL LSCALE(1,QMAX)
         IF(IPRO.NE.12) THEN
            CALL DUSHOW(25,27,QMAX)
c        write(6,*) ' qmax,shat ',qmax,sqrt(shh)
         ELSEIF(IPRO.EQ.12) THEN
            QMAX = MIN(QMAX,SNGL(P(25,4)))
            CALL DUSHOW(25,0,QMAX)
         ENDIF
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      IF(N.GE.28) IT=28
      NS=N
C...Generate spacelike parton shower (if required)
      IPU1=0
changed
      IF(IRES(1).EQ.1) IPU1=21
      IPU2=23

      XGTEST = (P(21,5)**2 + P(IT,5)**2+DBLE(PYPAR(22)))/PARL(22)

      IF(IPY(14).GE.1.AND.XGTEST.LT.0.999) THEN
         CALL PYSSPA(IPU1,IPU2)
         IF(LST(21).NE.0) THEN
            IF(LST(21).EQ.55) THEN
            ELSE
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF(.NOT.(IPY(14).GE.1.AND.XGTEST.LT.0.999).OR.LST(21).EQ.55) THEN
         IF(LST(21).EQ.55) LST(21)=0
         DO 120 I=NS+1,NS+4
            DO 120 J=1,5
               K(I,J)=0
               P(I,J)=0.D0
  120    V(I,J)=0.
Cadded
         K(NS+1,1)=11
         K(NS+1,3)=21
         K(NS+1,2)=KFL(2,1)
         DO 130 J=1,5
  130    P(NS+1,J)=P(21,J)
         K(NS+2,1)=-1
         K(NS+2,3)=NS+1
         K(NS+3,1)=13
         K(NS+3,2)=KFL(2,2)
         K(NS+3,3)=23
         K(NS+3,4)=23
         K(NS+3,5)=23
         P(NS+3,1)= 0.D0
         P(NS+3,2)= 0.D0
         P(NS+3,5)= 0.D0
         Q12 = 0.D0
         IF(Q2LM.GT.2.D0) THEN
            P(NS+3,3)=(P(IT,5)**2+Q2LM+Q12)*(P(21,4)-P(21,3))/(2.D0*
     +      Q2LM)
            P(NS+3,4)=-P(NS+3,3)
         ELSE
            P(NS+3,3)=-(P(IT,5)**2+Q2LM+Q12)/(P(21,4)+P(21,3))/2.D0
            P(NS+3,4)=-P(NS+3,3)
         ENDIF
         IF(NCALL.LT.NPRIN) THEN
            write(6,*) ' LMEPS NS+3,23,it ',P(NS+3,3),P(23,3),P(IT,5)
            write(6,*) ' LMEPS Q2 ',Q2,P(21,5)**2
         ENDIF

         K(NS+4,1)=-1
         K(NS+4,3)=NS+3
         P(NS+4,3)=23
         P(NS+4,4)=23
         P(24,1)=NS+3
         P(24,2)=NS+3
         K(23,4)=K(23,4)+(NS+3)*MSTU(5)
         K(23,5)=K(23,5)+(NS+3)*MSTU(5)
         IPU1=0
         IPU2=NS+3
         N=N+4
      ENDIF
C...Rotate and boost outgoing parton shower
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' before rotate and boost outgoing PS '
         CALL DULIST(2)
      ENDIF
changed      IF(ILEPTO.EQ.1.AND.IPRO.NE.12.AND.N.GT.31) THEN
      IF(ILEPTO.EQ.1.AND.IPRO.NE.12.AND.N.GT.31.and.IRES(1).EQ.0) THEN
         K(N+1,1)=0
         DO 140 J=1,4
  140    P(N+1,J)=P(NS+1,J)+P(NS+3,J)
         ROBO(1)=DLANGL(P(IT,3),SQRT(P(IT,1)**2+P(IT,2)**2))
         ROBO(2)=DLANGL(P(IT,1),P(IT,2))
         CALL DUDBRB(25,NS,0.D0,-ROBO(2),0.D0,0.D0,0.D0)
         CALL DUDBRB(25,NS,-ROBO(1),0.D0,0.D0,0.D0,0.D0)
         DELTAP(1) = P(N+1,1)
         DELTAP(2) = P(N+1,2)
         DELTAP(3) = P(N+1,3) - P(IT,3)
         DELTAP(4) = DSQRT(DELTAP(1)**2+DELTAP(2)**2+DELTAP(3)**2)
         DPLONG = -(P(IT,3)*DELTAP(3))/DELTAP(4)
         DBTOT = -(DPLONG*P(IT,4)- P(N+1,4)*DSQRT(DMAX1(0.D0,P(N+1,4)**
     +   2- P(IT,4)**2+DPLONG**2)))/(DPLONG**2+P(N+1,4)**2)
         DGAMMA = 1.D0/DSQRT(1.D0-DBTOT**2)
         DO 150 I = 1,3
            DROBO(I+2)=DELTAP(I)/(DGAMMA/(DGAMMA+1.D0)* (P(N+1,4)-
     +      DGAMMA*P(IT,4))+DGAMMA*P(IT,4))
  150    CONTINUE
         CALL DUDBRB(25,NS,0.D0,0.D0,DROBO(3),DROBO(4),DROBO(5))
C...End phi-correction
      ENDIF

      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' after rotate and boost outgoing PS '
         call dulist(1)
      ENDIF
C...Hadron remnant and primordial kt
      IPY(47)=N
      IPU1=0
      PYVAR(31)=XEL
      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' LMEPS before PYREMN: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
      CALL PYREMN(IPU1,IPU2)
      IF(NCALL.LE.NPRIN) CALL DULIST(2)
      IF(IPY(48).GE.1) THEN
         LST(21)=47+IPY(48)
         RETURN
      ENDIF
C...Rearrange partons along strings
      MSTU(24)=0
      MSTU(28)=0
      IF(IFPS.EQ.10) THEN
         MSTJO = MSTJ(105)
         MSTJ(105) = -1
      ENDIF
      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' LMEPS before DUPREP: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF

      CALL DUPREP(0)

cGB....
      IF(IFPS.EQ.10) THEN
         MSTJ(105) = MSTJO
      ENDIF
      IF(MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN
         if (iLMEPSPrint .le. MxLMEPSPrint) then
            iLMEPSPrint = iLMEPSPrint + 1
            WRITE(6,*) ' LMEPS: LUPREP error MSTU(24)= ',MSTU(24)
            write(6,*) ' color configuration ICOLORA,IRESPRO ',
     +           ICOLORA,IRESPRO
            write(6,*) ' res. photon ?  ',IRES(1)
            call DULIST(2)
         endif
         LST(21)=50
         RETURN
      ENDIF

C...Clean up event record -> order:
C...1=inc. lepton; 2=inc. nucleon; 3=exch boson; 4=scat. lepton;
C... (4+1) = rad. gamma if any
C...5+1=inc. parton before initial shower; 6+1=inc. parton at hard scattering
C...after shower; 7+1,8+1=first,second parton from hard scattering
C...before final shower
      IF(NCALL.LT.NPRIN) CALL DULIST(1)
      DO 160 J=1,5
         K(N+2,J)=K(3,J)
         P(N+2,J)=P(3,J)
         K(N+1,J)=K(4,J)
         P(N+1,J)=P(4,J)
         K(10,J)=0
         K(11,J)=0
         K(N+3,J)=KS(32,J)
         P(N+3,J)=PS(32,J)
         K(12,J)=KS(32,J)
         P(12,J)=PS(32,J)


         K(N+4,J)=KS(NIPH,J)
         P(N+4,J)=PS(NIPH,J)
  160 CONTINUE
      DO 170 I=1,20
         DO 170 J=1,5
            PS(I,J) = P(I,J)
  170 KS(I,J) = K(I,J)
      IRFL = 0
      DO 180 J=1,5
         IPOM=0
         IF(ILEPTO.EQ.1) THEN
            K(3,J)=K(N+4,J)
            P(3,J)=P(N+4,J)
            NIA1T = 3
            NIEL=4
            K(NIEL,J)=KS(9,J)
            P(NIEL,J)=PS(9,J)
            IF(KS(13,1).NE.0) THEN
               K(NIEL+1,J)=KS(13,J)
               P(NIEL+1,J)=PS(13,J)
               NIEL = 5
            ENDIF
            IF(KS(11,1).NE.0) THEN
               K(NIEL+1+IRFL,J)=KS(31,J)
               P(NIEL+1+IRFL,J)=PS(31,J)
               IPOM = 1
            ENDIF
            IF(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13) K(4,1)=1

         ENDIF
         K(NIEL+1+IRFL+IPOM,J)=K(N+1,J)
         P(NIEL+1+IRFL+IPOM,J)=P(N+1,J)
         K(NIEL+2+IRFL+IPOM,J)=K(NS+3,J)
         K(NIEL+2+IRFL+IPOM,3) = NIEL+1+IRFL+IPOM
         P(NIEL+2+IRFL+IPOM,J)=P(NS+3,J)
         NIA2 = NIEL + 2 + IRFL + IPOM
         K(NIEL+3+IRFL+IPOM,J)=K(25,J)
         P(NIEL+3+IRFL+IPOM,J)=P(25,J)
         NF1 = NIEL + 3 + IRFL + IPOM
         IF(IPRO.NE.12) THEN
            K(NIEL+4+IRFL+IPOM,J)=K(26,J)
            P(NIEL+4+IRFL+IPOM,J)=P(26,J)
            K(NIEL+5+IRFL+IPOM,J)=K(27,J)
            P(NIEL+5+IRFL+IPOM,J)=P(27,J)
            NF2 = NIEL + 5 + IRFL + IPOM
         ELSEIF(IPRO.EQ.12) THEN
            K(NIEL+4+IRFL+IPOM,J)=0
            NF2=NF1
         ENDIF
  180 CONTINUE
      K(3,3)=1
      K(4,3)=1
      K(NIEL+2+IRFL+IPOM,1)=21
      K(NIEL+2+IRFL+IPOM,3)=5+IRFL+IPOM
      K(NIEL+2+IRFL+IPOM,4)=0
      K(NIEL+2+IRFL+IPOM,5)=0
      K(NIEL+3+IRFL+IPOM,1)=21
      K(NIEL+3+IRFL+IPOM,3)=NIA1T
      K(NIEL+3+IRFL+IPOM,4)=0
      K(NIEL+3+IRFL+IPOM,5)=0
      IF(IPRO.NE.12) THEN
         K(NIEL+4+IRFL+IPOM,1)=21
         K(NIEL+4+IRFL+IPOM,3)=NIA1T
         K(NIEL+4+IRFL+IPOM,4)=0
         K(NIEL+4+IRFL+IPOM,5)=0
         K(NIEL+5+IRFL+IPOM,1)=21
         K(NIEL+5+IRFL+IPOM,3)=NIA1T
         K(NIEL+5+IRFL+IPOM,4)=0
         K(NIEL+5+IRFL+IPOM,5)=0
      ENDIF
      K(25,3) = NIEL + 3 + IRFL + IPOM
      K(26,3) = NIEL + 4 + IRFL + IPOM
      K(27,3) = NIEL + 5 + IRFL + IPOM
      IF(IRES(1).EQ.1) THEN
         P(NIEL+1,5)=-DSQRT(dabs(P(NIEL+1,1)**2 + P(NIEL+1,2)**2
     +   + P(NIEL+1,3)**2 - P(NIEL+1,4)**2))
         P(NIEL+2,5)=-DSQRT(dabs(P(NIEL+2,1)**2 + P(NIEL+2,2)**2
     +   + P(NIEL+2,3)**2 - P(NIEL+2,4)**2))
      ENDIF
      NIA1 = NIA1T
C...Deactivate obsolete lines 9, 10, 21, NS+1 (extra lines with boson)
      K(21,1)=0

      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
      K(13,1)=0
      K(14,1)=0
      K(15,1)=0
C...Zero irrelevant lines with K(I,1)<0
      DO 200 I=1,N
         IF(K(I,1).LT.0) THEN
            DO 190 J=1,5
               K(I,J)=0
  190       P(I,J)=0.D0
         ENDIF
  200 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' LMEPS: before LUEDIT '
         CALL DULIST(1)
      ENDIF
C...Delete internal parton lines, i.e. with K(I,1)=13,14
c      IF(MOD(1/10,10).EQ.0) THEN
      CALL DUEDIT(14)
c      ENDIF
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' LMEPS:after LUEDIT(14) '
         CALL DULIST(1)
      ENDIF
C...Delete empty lines
      CALL DUEDIT(12)
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' LMEPS:after LUEDIT(12) '
         CALL DULIST(1)
      ENDIF
      IF(IPRO.EQ.12.AND.PT2GEN.NE.0.D0) THEN
         phielt=DLANGL(P(1,1),P(1,2))
         CALL DUDBRB(0,0,0.D0,-phielt,0.D0,0.D0,0.D0)
         PZT = DSQRT(P(NF1,4)**2 - P(NF1,5)**2 - PT2GEN)
         STHARD = DLANGL(PZT,DSQRT(PT2GEN))
         CALL DUDBRB(NF1,NF1,STHARD,0.D0,0.D0,0.D0,0.D0)
         CALL DUDBRB(NF1+2,N,STHARD,0.D0,0.D0,0.D0,0.D0)
         phirot=PHIGEN
         CALL DUDBRB(NF1,NF1,0.D0,PHIROT,0.D0,0.D0,0.D0)
         CALL DUDBRB(NF1+2,N,0.D0,PHIROT,0.D0,0.D0,0.D0)
         CALL DUDBRB(0,0,0.D0,phielt,0.D0,0.D0,0.D0)
c       nre=5
c       te11 = P(nre,2)*p(1,3)-p(nre,3)*p(1,2)
c       te12 = P(nre,3)*p(1,1)-p(nre,1)*p(1,3)
c       te13 = P(nre,1)*p(1,2)-p(nre,2)*p(1,1)
c       te21 = P(nre,2)*p(nf1,3)-p(nre,3)*p(nf1,2)
c       te22 = P(nre,3)*p(nf1,1)-p(nre,1)*p(nf1,3)
c       te23 = P(nre,1)*p(nf1,2)-p(nre,2)*p(nf1,1)
c       prod = te11*te21 + te12*te22 + te13*te23
c       write(6,*) 'sc vec_prod ',prod
c       prod=prod/dsqrt(te11**2+te12**2+te13**2)
c       prod=prod/dsqrt(te21**2+te22**2+te23**2)
c       write(6,*) ' LMEPS NF1,phielt,phirot ',NF1,phielt,phirot
c       write(6,*) ' LMEPS checking pt2,phi,phi_calc',pt2gen,phigen
c     &    ,dacos(prod)
c       write(6,*) ' LMEPS checking cos(phi,phi_calc)',cos(phigen)
c     &    ,prod
c         call dulist(1)
      ENDIF

c check for enegy momentum conservation
      DO 220 I=1,N
         IF(P(I,5).LT.0.0D0.OR.K(I,1).GT.2) GOTO 220
         DO 210 J=1,4
            IF(ABS(P(I,J)).LE.1.E-6) THEN
               P(I,J)=0.0D0
            ENDIF
  210    CONTINUE
  220 CONTINUE
      CHEC = 5.D-3
      P1 =ABS(P(1,1)+P(2,1))
      P2 =ABS(P(1,2)+P(2,2))
      P3 =ABS(P(1,3)+P(2,3))
      P4 =ABS(P(1,4)+P(2,4))
      PC1 = (ABS(DPLU(0,1)) - P1)/P4
      PC2 = (ABS(DPLU(0,2)) - P2)/P4
      PC3 = (ABS(DPLU(0,3)) - P3)/P4
      PC4 = (ABS(DPLU(0,4)) - P4)/P4
      IF(DABS(PC1).GT.CHEC.
     +   OR.DABS(PC2).GT.CHEC.
     +   OR.DABS(PC3).GT.CHEC.
     +   OR.DABS(PC4).GT.CHEC) THEN
         write(6,*) ' LMEPS: energy of final particles not correct'
     +   ,chec
         write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
         write(6,*) ' Q2 = ',Q2,' Q2_calc = ',DOT1(3,3)
         write(6,*) ' yy = ',yy,' y_calc = ',DOT1(2,3)/DOT1(2,1)
         write(6,*) ' x_bj = ',Q2/PARL(22),' W = ',SQRT(-Q2+PARL(22))
         write(6,*) ' IPRO = ',IPRO
         call DULIST(1)
         LST(21) = 100
      ENDIF
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' in gamma pom system '
         CALL DULIST(1)
      ENDIF
      IF(ILEPTO.EQ.1) THEN
         CALL DUDBRB(0,0,STHETA(2),SPHI(2),0.D0,0.D0,0.D0)
         CALL DUDBRB(0,0,0.D0,0.D0,DBETA(2,1),DBETA(2,2),DBETA(2,3))
      ENDIF
c... change z coordinate from LEPTO61 standard to rapgap standard
      do 230 i=1,n
  230 p(i,3) = - p(i,3)
      DO 240 I=1,4
         IF(ABS(P(NIEL,I)).LE.1.E-6) THEN
            P(NIEL,I)=0.0D0
         ENDIF
         P(1,I)=PBEAM(1,I)
         P(2,I)=PBEAM(2,I)
  240 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' after boost dbeta in lmeps'
         CALL DULIST(1)
      ENDIF
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' LMEPS final '
         CALL DULIST(1)
      ENDIF
      RETURN
      END
