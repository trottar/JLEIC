*CMZ :  2.08/02 20/07/99  11.39.37  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.21.38  by  Hannes Jung
*CMZ :  2.07/03 16/05/99  09.53.09  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE PARTI(IP1,YX,YVAL,WEIGHT,IFL,IST)
      IMPLICIT NONE
	Integer IP1,IFL,IST
	Double Precision YX,YVAL,WEIGHT
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

*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
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


*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.
      Integer NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
	Double Precision PHEP,VHKK
      PARAMETER (NMXHEP=2000)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     +                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
     +                PHEP(5,NMXHEP),VHKK(4,NMXHEP)
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEND.
C     from HERACLES;
	Integer ICHNN,LST,IRES
      COMMON /HSCHNN/ ICHNN
      COMMON/EPPARA/LST(30),IRES(2)
	Double Precision PHI,DBGAM
      COMMON/DIFFA/ PHI
      DOUBLE PRECISION ME
      DIMENSION DBGAM(4)
      REAL SNGL
	Double Precision XG1,ALPH_EM,FGAM,SQRTS,EC,X11,COST2,SPHE,CPHE
	Double Precision PEP,PEZ,PDLE,CSPHI,XP1,BOCHCK,THE,PHIE
	Double Precision SPHI,STHETA,PT
	Integer I,NIPH,IE,NFI,IG
C IST = 0, ONLY PARTICLE MOMENTA, BUT SCALE Q2Q IN STRUCTURE FUNCTION
C NOT YET DEFINED.
C IST = 1, CALCULATE ONLY STRUCTURE FUNCTION SINCE Q2Q NOW DEFINED
C
C IF PARTICLE IP1 = 11 (ELECTRON) USE EQUIVALENT PHOTON APPROXIMATION
C
c      write(6,*) 'parti  Q2 = ',Q2,' yx = ',yx
      XG1 = DBLE(XEL)
      YVAL =  999999.D0
      WEIGHT = -999999.D0
      ME = DBLE(ULMASS(IP1))
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))

C IF 1ST INCOMING PARTICLE IS A ELECTRON AND LOOKING FOR DIRECT PHOTON
      IF((IABS(IP1).EQ.11.OR.IABS(IP1).EQ.13)) THEN
C
         IF(IABS(IP1).EQ.11.OR.IABS(IP1).EQ.13) THEN
C...          WEIZAECKER WILLIAMS APPROXIMATION FOR BREMSTRAHLUNG (EPA)
C...          WITH Q2 DEPENDENCE
            FGAM = (1.D0 - YX +YX*YX/2.D0)/Q2/YX - ME*ME*YX/Q2/Q2

            YVAL = ALPH_EM * FGAM / PI
         ENDIF
         IF(QMI.EQ.0.D0) THEN
            WEIGHT = YX*Q2*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
         ELSE
c new try weighting with 1/q**4
            WEIGHT = YX*DLOG(YMAX/YMIN)*Q2**2*(Q2MAX-Q2MIN)/Q2MIN/
     +      Q2MAX
         ENDIF
C...

C...      CALCULATE COS OF ELECTRON VIA FORMULA GIVEN IN CAMPBELL ET AL
C...      CAN.J.PHYS. 39 1981 (1742)
C...      VIA X Q2 Y
C...      X11= Q2/(2P.Q)
         SQRTS=DSQRT(SSS)
         EC= 0.5D0*SQRTS
         IF(IHERAC.EQ.0) THEN
            X11= Q2/YX/SSS
            COST2= ((1.D0 - YX) - X11*YX)/((1.D0 - YX) + X11*YX)
            IF(COST2.GT.1.) WRITE(6,*) 'YX Q2 COST2 ',YX,Q2,COST2
C...      COST2 IS SCATTERING ANGLE OF ELECTRON
            THE= DACOS(COST2)
            PHIE = PHI
            IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
            SPHE=DSIN(PHIE)
            CPHE=DCOS(PHIE)
C...      PEL1 = 4 VECTOR OF SCATTERED ELECTRON IN EP CMS SYSTEM
            N=N+1
            K(N,1)=1
            IF(INTER.LT.2) THEN
               K(N,2)=K(1,2)
            ELSEIF(INTER.EQ.2) THEN
               K(N,2)=ISIGN(1,K(1,2))*12
            ENDIF
            K(N,3)=IFL
c. 2*P.l_e
c new gymnastics for numerical stability
            PEP = DSQRT(P(2,3)**2 + P(2,5)**2)
            PEZ = P(2,3)
            PDLE = 2.D0*DOT1(1,2)
            X11 = Q2/PDLE/YX
            P(N,5)=DBLE(ULMASS(K(N,2)))
            P(N,4) = P(1,4) + (PEZ*(Q2+P(N,5)**2)/ (2.D0*P(1,4)) -
     +      YX*PDLE/2.D0)/(PEP+PEZ)
            P(N,3) = P(1,4) - (PEP*(Q2+P(N,5)**2)/ (2.D0*P(1,4)) + YX*
     +      PDLE/2.0D0)/(PEP+PEZ)
            P(N,3)=-P(N,3)
            PT = Q2 - (Q2 + P(N,5)**2)/P(1,4)/(PEP+PEZ)* (YX*PDLE/2.D0 +
     +      P(2,5)**2*(Q2+P(N,5)**2)/4.D0/ P(1,4)/(PEP+PEZ))
            P(N,1)=DSQRT(DMAX1(0.D0,PT))*CPHE
            P(N,2)=DSQRT(DMAX1(0.D0,PT))*SPHE

C...      PPH = 4 VECTOR OF (VIRTUAL) PHOTON IN EP CMS SYSTEM
            N=N+1
            K(N,1)=21
            K(N,2)=KEB
            K(N,3)=IFL
            P(N,1)= P(1,1) - P(3,1)
            P(N,2)= P(1,2) - P(3,2)
            P(N,3)= P(1,3) - P(3,3)
            P(N,4)= P(1,4) - P(3,4)
            P(N,5)= -SQRT(ABS(DOT1(N,N)))
c            write(6,*) 'parti: ',P(N,5)**2,q2
            NIA1 = N
            NIPH = N
         ELSEIF(IHERAC.EQ.1) THEN
c            CALL PEPEVT
            IE=N+1
            NFI = N + 2
            IF(PHEP(4,3).NE.0.0) THEN
               IG=N+2
               NFI = N + 3
c radiative gamma
               K(IG,1)=1
               K(IG,2)=KPH
c radiated gamma has origin in e or e'
c               write(6,*) 'parti: ICHNN=',ICHNN
               IF(ICHNN.EQ.6.OR.ICHNN.EQ.12) THEN
                  K(IG,3)=1
               ELSE
                  K(IG,3)=4
               ENDIF
            ENDIF
c scattered electron
            K(IE,1)=1
            K(IE,2)=K(1,2)
            K(IE,3)=IFL
c exchanged  boson
            IF(IABS(IDHEP(1)).EQ.11) THEN
               KEB=KPH
            ELSEIF(IABS(IDHEP(1)).EQ.12) THEN
               KEB=-24*ISIGN(1,K(1,2))
            ELSE
               write(6,*) ' PARTI: scattered lepton not known ',IDHEP(1)
            ENDIF
            K(NFI,1)=21
            K(NFI,2)=KEB
            K(NFI,3)=IFL
            DO 10  I=1,5
               P(IE,I) = PHEP(I,1)
               P(NFI,I) = PHEP(I,5) - PHEP(I,1) - PHEP(I,3)
               IF(PHEP(4,3).NE.0.0) THEN
                  P(IG,I) = PHEP(I,3)
               ENDIF
   10       CONTINUE
c calculate phi
            CSPHI = P(IE,1)/SQRT(P(IE,1)**2+P(IE,2)**2)
            PHI = DACOS(CSPHI)
            N = NFI
            NIA1 = NFI
            NIPH = NIA1
            DO 20  I = 3,N
   20       P(I,3) = -P(I,3)
            P(NFI,5)= -SQRT(ABS(DOT1(NFI,NFI)))
c         write(6,*) ' parti: q2hs =',q2hs,p(nfi,5)**2
C BOOST TO EP CMS
            CALL DUDBRB(3,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/
     +      CM(4))
         ELSE
            write(6,*) ' not valid for IHERAC ',IHERAC
         ENDIF
         IF(IRES(1).EQ.0) THEN
            NIR1 = -99999
         ELSEIF(IRES(1).EQ.1) THEN
c          write(6,*) 'parti yx = ',yx,' xel ',xg1
C         XG1 = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XP1 = E_GLUON/E_PHOTON
            XP1=XG1/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
            N = N + 1
            K(N+1,1)=21
            K(N+1,2)=21
            K(N+1,3)=IFL
            P(N+1,1)= P(N-1,1)*XP1
            P(N+1,2)= P(N-1,2)*XP1
            P(N+1,3)= P(N-1,3)*XP1
cccc            P(N+1,4)= DABS(P(N-1,3))*XP1
            P(N+1,4)= DABS(P(N-1,4))*XP1
            P(N+1,5)= DBLE(ULMASS(21))
            NIA1 = N+1
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
            K(N,1)=1
            K(N,2)=21
            K(N,3)=IFL
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,5)= DBLE(ULMASS(21))
            P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
            P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
            NIR1 = N

c            write(6,*) ' parti yx,xg1,xp1',yx,xg1,xp1
c            call dulist(1)
C NOW BOOST TO GAMMA PROTON FRAME
            DBGAM(1) = P(NIPH,1) + P(2,1)
            DBGAM(2) = P(NIPH,2) + P(2,2)
            DBGAM(3) = P(NIPH,3) + P(2,3)
            DBGAM(4) = P(NIPH,4) + P(2,4)
            BOCHCK = (DBGAM(1)/DBGAM(4))**2 + (DBGAM(2)/DBGAM(4))**2 +
     +      (DBGAM(3)/DBGAM(4))**2
            BOCHCK = DSQRT(BOCHCK)
            IF(BOCHCK.GT.0.99999999D0) goto 30
            CALL DUDBRB(0,N,0.D0,0.D0,-DBGAM(1)/DBGAM(4),-DBGAM(2)/
     +      DBGAM(4), -DBGAM(3)/DBGAM(4))
            SPHI = DLANGL(P(NIPH,1),P(NIPH,2))
            call DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
            STHETA = DLANGL(P(NIPH,3),P(NIPH,1))
            call DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
            P(N+1,1) = 0.D0
            P(N+1,2) = 0.D0
            P(N+1,3) = XP1*P(NIPH,3)
            P(N+1,4) = XP1*(P(NIPH,4)*P(2,4)-P(NIPH,3)*P(2,3))/
     +                 (P(2,4)-P(2,3))
            P(N+1,3) = P(N+1,4)
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,4)= P(1,4) - P(N+1,4)- P(3,4)
            P(N,5)= DSQRT(DABS(P(N,1)**2+P(N,2)**2+P(N,3)**2-P(N,4)**2))
            N = N + 1
c            call dulist(1)
            call DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
            call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
            CALL DUDBRB(0,N,0.D0,0.D0,DBGAM(1)/DBGAM(4),DBGAM(2)/
     +      DBGAM(4), DBGAM(3)/DBGAM(4))
c            call dulist(1)
c            pause
         ENDIF
C IF 1ST INCOMING PARTICLE IS A PHOTON
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(1).EQ.0) THEN
C...      4 VECTOR OF (REAL) PHOTON IN CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KPH
         K(N,3)=IFL
         P(N,1)= P(1,1)
         P(N,2)= P(1,2)
         P(N,3)= P(1,3)
         P(N,4)= P(1,4)
         P(N,5)=DBLE(ULMASS(KPH))
         NIA1 = N
         NIR1 = -99999
         YVAL = 1.D0
         WEIGHT = 1.D0
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(1).EQ.1) THEN
         WEIGHT = 1.0D0
         YVAL = 1.0D0

C         XG1 = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XP1 = E_GLUON/E_PHOTON
         XP1=XG1/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
         N = N + 1
         K(N+1,1)=21
         K(N+1,2)=21
         K(N+1,3)=IFL
         P(N+1,1)= P(N-1,1)*XP1
         P(N+1,2)= P(N-1,2)*XP1
         P(N+1,3)= P(N-1,3)*XP1
cccc            P(N+1,4)= DABS(P(N-1,3))*XP1
         P(N+1,4)= DABS(P(N-1,4))*XP1
         P(N+1,5)= DBLE(ULMASS(21))
         NIA1 = N+1
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
         K(N,1)=1
         K(N,2)=21
         K(N,3)=IFL
         P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
         P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
         P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
         P(N,5)= DBLE(ULMASS(21))
         P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
         P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
         NIR1 = N
         N = N + 1
c            write(6,*) ' parti yx,xg1,xp1',yx,xg1,xp1
c            call dulist(1)
      ENDIF
      RETURN
   30 write(6,*) ' PARTI boost error '
      RETURN
      END
