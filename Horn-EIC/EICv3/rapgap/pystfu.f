*CMZ :  2.08/00 06/06/99  16.42.15  by  Hannes Jung
*CMZ :  2.07/03 23/05/99  18.33.17  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
C*********************************************************************

      SUBROUTINE PYSTFU(KF,XI,Q2I,XPQ)
	Implicit None
	Integer KF
      REAL XI,Q2I,XPQ,SLO,SHI,EPS,RCCB,XPQHF
      DIMENSION XPQ(-6:6)
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

*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
      DOUBLE PRECISION QG2
      COMMON/SEMIH/ QG2
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
C...MODIFIED FOR USE IN HERACLES (BY H.SPIESBERGER, 8.2.93)
      Real PYSTOP,PYSLAM
	Integer NPYMOD,NPYMAX,NPYMIN
      COMMON /PYSTFUC/ PYSTOP,PYSLAM,NPYMOD,NPYMAX,NPYMIN
C Input:
C   PYSTOP = top mass
C   NPYMAX = maximal flavour
C   NPYMIN = minimal flavour
C   NPYMOD = Choice of parametrization
      Integer LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND,LPAR,LPARIN,IPART
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12),IPART
*KEEP,RGRAHER.
      REAL XPQDIF,XPQPI
	Integer IHERPYS
	Integer NBQ2,NBX
      PARAMETER (NBQ2=20)
      PARAMETER (NBX=20)
      COMMON /RAHER/ IHERPYS,XPQDIF(-6:6,NBX,NBQ2),XPQPI(-6:6,NBX,NBQ2)
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGLQ2.
      DOUBLE PRECISION Q2SUPP
      COMMON/LOWQ2S/Q2SUPP
*KEND.
      Integer NQPY,NXPY
      PARAMETER (NQPY=60,NXPY=60)

      REAL XPC(NXPY,NQPY),XXI(NXPY),Q2II(NQPY),XPRT,Q2T
      REAL CUT,XLP,YLP,W2LP,Q2LP,ULP
      DOUBLE PRECISION PARL
	Integer LLst
      COMMON /RAPTOU/ CUT(14),LLST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      LOGICAL FIRST
      REAL X,Q2
	Integer KPHF
      COMMON /HEAVYF/X,Q2,KPHF
	Real GEV2NB
      EXTERNAL XPQHF
	Double Precision dq,dx
	Real dezot,qg20,qgx,alamz,xz0,fp,phi1,xd,qd,x1p,x2p
	Integer I,J,ix,iq
      DATA GEV2NB/.3893857E+6/
      DATA FIRST/.TRUE./
      SLO = 0.0
      SHI = 1.0
      EPS = 1.E-4
      Q2 = Q2I
      X = XI
      IF(KF.EQ.22.AND.Q2LP.NE.0.0) THEN
         CALL RYSTGA(X,Q2,Q2LP,XPQ)
         RETURN
      ENDIF
      IF(ISEMIH.EQ.1.AND.IDIR.EQ.1) THEN
         IF(IPRO.EQ.10.OR.IPRO.EQ.11.OR.IPRO.EQ.13.OR.IPRO.EQ.14) THEN
c unintegrated gluon for semihard approach
            DO 10  I=-6,6
   10       XPQ(I)=0.0
            QG20= 2.0
            ALAMZ = 0.056
            XZ0 = 0.333333
            IF(X.LT.XZ0) THEN
c            write(6,*) 'XZ0,X,QG2',XZ0,x,QG2
               DEZOT = EXP(3.56*SQRT(ALOG(XZ0/X)))
c               write(6,*) 'ALAMZ,DEZOT,ALAMZ**2*DEZOT',
c     &                    ALAMZ,DEZOT,ALAMZ**2*DEZOT
               QGX = QG20 + ALAMZ**2*DEZOT
            ELSE
               QGX = QG20
            ENDIF
cccc               QGX = QG20
c               write(6,*) ' QGX ALAMZ DEZOT QG20',QGX,ALAMZ,DEZOT,QG20
            IF(SNGL(QG2).GT.QGX) THEN
               FP = (QGX/SNGL(QG2))**2
c               write(6,*) ' QGX>QG2 ,FP',QGX,QG2,FP
            ELSE
c               write(6,*) ' QGX<QG2 ,FP',QGX,QG2,FP
               FP = 1.0
            ENDIF
ccc g1
            PHI1 = 0.97E6 * 0.05/(X+0.05)*(1.0 -X)**3
ccc g2
ccc            PHI1 = 0.65E6 * (1.0 - X)**5
c divide by GEV2NB because we want the gluon dimensionless
            XPQ(0) = PHI1*FP/GEV2NB
c
c            write(6,*) PHI1,FP,XPQ(0)
            GOTO 20
         ENDIF
      ENDIF
      CALL RYSTFU(KF,XI,Q2I,XPQ)
c      write(6,*) ' after rystfu ',xpq(4)
   20 CONTINUE
      IF(IHF.EQ.1) THEN
         KPHF = 4
         IF(IHFLA.GE.4) KPHF=IHFLA
         Q2 = Q2I
         X = XI
c         write(6,*) ' pystfu x,q2,xpq(4) old',x,q2,xpq(4)
         IF(FIRST) THEN
            DQ = (LOG(SSS) - LOG(QMI))/DFLOAT(NQPY-1)
            DX = (LOG10(1.D0) - LOG10(QMI/SSS))/DFLOAT(NXPY-1)
            DO 30  I=1,NXPY
               DO 30  J=1,NQPY
                  XXI(I) = SNGL(10**(LOG10(QMI/SSS) + DX*DFLOAT(I-1)))
                  IF(XXI(I).GE.1.D0) XXI(I) = 0.999
                  X = XXI(I)
                  Q2II(J) = SNGL(EXP(LOG(QMI) + DQ*DFLOAT(J-1)))
                  Q2 = Q2II(J)
                  CALL INTGA(SLO,SHI,XPQHF,EPS,RCCB)
                  XPC(I,J)=RCCB/2.
   30       CONTINUE
            FIRST=.FALSE.
         ENDIF
         XPRT = XI
         Q2T= Q2I
         IF(Q2T.LT.Q2II(1)) Q2T=Q2II(1)
         IF(XPRT.LT.XXI(1).OR.XPRT.GT.XXI(NXPY) .OR.Q2T.LT.Q2II(1)
     +    .OR.Q2T.GT.Q2II(NQPY)) THEN
c            WRITE(6,*) 'PYSTFU: X or Q2 values outside grid '
c            WRITE(6,*) ' X_min ',XXI(1),' X_max ',XXI(20), ' actual '
c     +      //'X ', XPRT
c            WRITE(6,*) ' Q2_min ',Q2II(1),' Q2_max ',Q2II(15), ' '
c     +      //'actual Q2 ',Q2T
c             write(6,*) ' this was for IPRO =',ipro,idir,ng,npom
            IF(XPRT.LT.XXI(1)) XPRT=XXI(1)
            IF(XPRT.GT.XXI(NXPY)) XPRT=XXI(NXPY)
            IF(Q2T.LT.Q2II(1)) Q2T = Q2II(1)
            IF(Q2T.GT.Q2II(NQPY)) Q2T = Q2II(NQPY)
         ENDIF
         IX = 0
   40    IX = IX + 1
         IF(XPRT.GT.XXI(IX+1)) GOTO 40
         IQ = 0
   50    IQ = IQ + 1
         IF(Q2T.GT.Q2II(IQ+1)) GOTO 50

         XD = (XPRT - XXI(IX))/(XXI(IX+1)-XXI(IX))
         QD = (Q2T - Q2II(IQ))/(Q2II(IQ+1) - Q2II(IQ))
         X1P=(XPC(IX+1,IQ)-XPC(IX,IQ))*XD +XPC(IX,IQ)
         X2P=(XPC(IX+1,IQ+1)-XPC(IX,IQ+1))*XD + XPC(IX,IQ+1)
         XPQ(KPHF) = (X2P-X1P)*QD + X1P
         XPQ(-KPHF) = XPQ(KPHF)
         X=XI
         Q2=Q2I
c         CALL INTGA(SLO,SHI,XPQHF,EPS,RCCB)

c         write(6,*) ' xpq(4),charm new',xpq(4)
c         XPQ(4)=RCCB/2.
c         XPQ(-4)=XPQ(4)
c         write(6,*) ' xpq(4),charm new',xpq(4),rccb/2.
         DO 60 I=1,6
            IF(I.EQ.IABS(KPHF)) GOTO 60
            XPQ(I) = 0.0
            XPQ(-I) = 0.0
   60    CONTINUE
      ENDIF

      RETURN
      END
